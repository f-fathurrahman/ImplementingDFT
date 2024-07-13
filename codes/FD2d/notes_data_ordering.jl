# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
using Printf

# %%
using Statistics: mean

# %%
using SparseArrays: sparse

# %%
using LinearAlgebra: I, kron

# %%
speye(N::Int64) = sparse( Matrix(1.0I, N, N) );

# %%
const ⊗ = kron;

# %%
include("../FD1d/init_FD1d_grid.jl")
include("../FD1d/build_D2_matrix_11pt.jl");

# %% [markdown]
# A simpler definition of `FD2dGrid`:

# %%
struct FD2dGrid
    Npoints::Int64
    Lx::Float64
    Ly::Float64
    Nx::Int64
    Ny::Int64
    hx::Float64
    hy::Float64
    dVol::Float64
    x::Array{Float64,1}
    y::Array{Float64,1}
    r::Array{Float64,2}
    #
    idx_ip2xy::Array{Int64,2}
    idx_xy2ip::Array{Int64,2}
end

# %%
# The first version, loop over Ny, then Nx
function create_FD2dGrid_v1(
    x_domain::Tuple{Float64,Float64}, Nx::Int64,
    y_domain::Tuple{Float64,Float64}, Ny::Int64;
)
    x, hx = init_FD1d_grid(x_domain, Nx)
    y, hy = init_FD1d_grid(y_domain, Ny)
    Lx = x_domain[2] - x_domain[1]
    Ly = y_domain[2] - y_domain[1]
    dVol = hx*hy
    Npoints = Nx*Ny
    r = zeros(2,Npoints)
    ip = 0
    idx_ip2xy = zeros(Int64,2,Npoints)
    idx_xy2ip = zeros(Int64,Nx,Ny)
    # We loop over Ny then Nx
    for j in 1:Ny, i in 1:Nx
        ip = ip + 1
        r[1,ip] = x[i]
        r[2,ip] = y[j]
        idx_ip2xy[1,ip] = i
        idx_ip2xy[2,ip] = j
        idx_xy2ip[i,j] = ip
    end
    return FD2dGrid(Npoints, Lx, Ly, Nx, Ny, hx, hy, dVol, x, y, r, idx_ip2xy, idx_xy2ip)
end;

# %%
# The first version, loop over Nx, then Ny
function create_FD2dGrid_v2(
    x_domain::Tuple{Float64,Float64}, Nx::Int64,
    y_domain::Tuple{Float64,Float64}, Ny::Int64;
)
    x, hx = init_FD1d_grid(x_domain, Nx)
    y, hy = init_FD1d_grid(y_domain, Ny)
    Lx = x_domain[2] - x_domain[1]
    Ly = y_domain[2] - y_domain[1]
    dVol = hx*hy
    Npoints = Nx*Ny
    r = zeros(2,Npoints)
    ip = 0
    idx_ip2xy = zeros(Int64,2,Npoints)
    idx_xy2ip = zeros(Int64,Nx,Ny)
    # We loop over Nx then Ny
    for i in 1:Nx, j in 1:Ny
        ip = ip + 1
        r[1,ip] = x[i]
        r[2,ip] = y[j]
        idx_ip2xy[1,ip] = i
        idx_ip2xy[2,ip] = j
        idx_xy2ip[i,j] = ip
    end
    return FD2dGrid(Npoints, Lx, Ly, Nx, Ny, hx, hy, dVol, x, y, r, idx_ip2xy, idx_xy2ip)
end;

# %% [markdown]
# This is how we generate our 2d grid:

# %%
x_domain = (-5.0, 4.0); Nx = 3;
y_domain = (-3.0, 2.0); Ny = 4;

# %%
grid_v1 = create_FD2dGrid_v1(x_domain, Nx, y_domain, Ny);
grid_v2 = create_FD2dGrid_v2(x_domain, Nx, y_domain, Ny);

# %%
grid_v1.x

# %%
grid_v2.x

# %%
grid_v1.y

# %%
grid_v2.y

# %%
println("Grid v1")
for ip in 1:Npoints
    @printf("%4d %18.10f %18.10f\n", ip, grid_v1.r[1,ip], grid_v1.r[2,ip])
end;

# %%
println("Grid v2")
for ip in 1:Npoints
    @printf("%4d %18.10f %18.10f\n", ip, grid_v2.r[1,ip], grid_v2.r[2,ip])
end;

# %% [markdown]
# If we reshape `r[1,:]` to a 2d array:

# %%
reshape(grid_v1.r[1,:], (Nx,Ny))

# %% [markdown]
# But for grid_v2, the shape is not (Nx,Ny):

# %%
reshape(grid_v2.r[1,:], (Nx,Ny))

# %% [markdown]
# The proper way to reshape `grid_v2.r[1,:]` is by using shape `(Ny,Nx)` instead of `(Nx,Ny)`:

# %%
reshape(grid_v2.r[1,:], (Ny,Nx))

# %% [markdown]
# # Laplacian operator

# %% [markdown]
# These functions do not utilize `Nx` and `Ny` ordering, only the sequence of `kron`'ed matrices are different. This means that they can accept both `grid_v1` and `grid_v2`, only the sequence of outer products between D2x, D2y, Ix, and Iy are different. 

# %%
function build_nabla2_matrix_v1( grid )
    Nx = grid.Nx
    Ny = grid.Ny
    hx = grid.hx
    hy = grid.hy
    D2x = build_D2_matrix_11pt(Nx, hx)
    Ix = speye(Nx)
    D2y = build_D2_matrix_11pt(Ny, hy)
    Iy = speye(Ny)
    #∇2 = kron(D2y, Ix) + kron(Iy, D2x) # y at the most left
    ∇2 = D2y ⊗ Ix + Iy ⊗ D2x
    return ∇2
end;

# %%
function build_nabla2_matrix_v2( grid )
    Nx = grid.Nx
    Ny = grid.Ny
    hx = grid.hx
    hy = grid.hy
    D2x = build_D2_matrix_11pt(Nx, hx)
    Ix = speye(Nx)
    D2y = build_D2_matrix_11pt(Ny, hy)
    Iy = speye(Ny)
    #∇2 = kron(D2x, Iy) + kron(Ix, D2y) # at the most left
    ∇2 = D2x ⊗ Iy + Ix ⊗ D2y
    return ∇2
end;

# %% [markdown]
# Here is a function to create a Gaussian function. It loops over linear index `ip` instead of individual indices for x and y directions, so it should be applicable for both grid_v1 and grid_v2.

# %%
function my_gaussian( grid::FD2dGrid; α=1.0, center=(0.0, 0.0) )
    Npoints = grid.Npoints
    f = zeros(Npoints)
    for ip in 1:Npoints
        x = grid.r[1,ip] - center[1]
        y = grid.r[2,ip] - center[2]
        r2 = x^2 + y^2
        f[ip] = exp(-α*r2)
    end
    return f
end;

# %%
function my_d2_gaussian( grid::FD2dGrid; α=1.0, center=(0.0, 0.0) )
    Npoints = grid.Npoints
    f = zeros(Npoints)
    for ip in 1:Npoints
        x = grid.r[1,ip] - center[1]
        y = grid.r[2,ip] - center[2]
        r2 = x^2 + y^2
        f[ip] = 4*α*(r2 - 1)*exp(-α*r2)
    end
    return f
end;

# %%
x_domain = (-10.0, 10.5); Nx = 120;
y_domain = (-10.0, 10.2); Ny = 120;
center = (0.0, 0.1);

grid_v1 = create_FD2dGrid_v1(x_domain, Nx, y_domain, Ny)
grid_v2 = create_FD2dGrid_v2(x_domain, Nx, y_domain, Ny)

f_v1 = my_gaussian(grid_v1, center=center)
f_v2 = my_gaussian(grid_v2, center=center);

# %%
d2f_analytic_v1 = my_d2_gaussian(grid_v1, center=center)
d2f_v1 = build_nabla2_matrix_v1(grid_v1) * f_v1;
diff_v1 = abs.(d2f_v1 - d2f_analytic_v1);
mean(diff_v1)

# %%
d2f_analytic_v2 = my_d2_gaussian(grid_v2, center=center)
d2f_v2 = build_nabla2_matrix_v2(grid_v2) * f_v2;
diff_v2 = abs.(d2f_v2 - d2f_analytic_v2);
mean(diff_v2)

# %% [markdown]
# It seems that for the v1 the correct order for building Laplacian matrix is:
# $$
# \nabla^2 = \mathbb{I}_{y} \otimes D^{(2)}_{x} +  D^{(2)}_{y} \otimes \mathbb{I}_{x}
# $$
# i.e. x direction is at the most right, or y is the most left.

# %% [markdown]
# Anf for the v2 the correct order for building Laplacian matrix is:
# $$
# \nabla^2 = D^{(2)}_{x} \otimes \mathbb{I}_{y} + \mathbb{I}_{x} \otimes D^{(2)}_{y}
# $$
# i.e. x direction is at the most left, or y is the most right.

# %%

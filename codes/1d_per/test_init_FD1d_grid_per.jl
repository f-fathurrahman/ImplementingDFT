include("init_FD1d_grid_per.jl")

function main()
    A = -5.0
    B =  5.0
    Npoints = 8
    x, h = init_FD1d_grid_per( A, B, Npoints )
    println(x)
    println(h)
end

main()
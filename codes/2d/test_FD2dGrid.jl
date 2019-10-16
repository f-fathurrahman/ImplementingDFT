using Printf

include("FD2dGrid.jl")

function main()
    Nx = 50
    Ny = 50
    fdgrid = FD2dGrid( (-5.0,5.0), Nx, (-5.0,5.0), Ny )
    println(fdgrid)
end

main()
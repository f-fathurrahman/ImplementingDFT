include("build_D2_matrix_3pt_per.jl")
include("build_D2_matrix_5pt_per.jl")
#include("build_D2_matrix_7pt.jl")
#include("build_D2_matrix_9pt.jl")

function test_D2_3pt()
    D2 = build_D2_matrix_3pt_per(5, 1.0)
    display(D2)
    println()
end


function test_D2_5pt()
    D2 = build_D2_matrix_5pt_per(8, sqrt(1.0/12.0))
    display(D2)
    println()
end

#=
function test_D2_7pt()
    D2 = build_D2_matrix_7pt(9, sqrt(1.0/180.0))
    display(D2)
    println()
end

function test_D2_9pt()
    D2 = build_D2_matrix_9pt(9, sqrt(1.0/5040.0))
    display(D2)
    println()
end
=#

#test_D2_3pt()
test_D2_5pt()
#test_D2_7pt()
#test_D2_9pt()
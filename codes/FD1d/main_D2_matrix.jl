include("build_D2_matrix_3pt.jl")
include("build_D2_matrix_5pt.jl")
include("build_D2_matrix_7pt.jl")
include("build_D2_matrix_9pt.jl")
include("build_D2_matrix_11pt.jl")

function test_D2_3pt()
    D2 = build_D2_matrix_3pt(5, 1.0)
    display(D2)
    println()
end

function test_D2_5pt()
    D2 = build_D2_matrix_5pt(9, sqrt(1.0/12.0))
    display(D2)
    println()
end

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

function test_D2_11pt()
    D2 = build_D2_matrix_11pt(11, sqrt(1.0/25200.0))
    display(D2)
    println()
end

#test_D2_3pt()
#test_D2_5pt()
#test_D2_7pt()
test_D2_9pt()
test_D2_11pt()
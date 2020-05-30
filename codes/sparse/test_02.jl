using SparseArrays

include("../FD1d/build_D2_matrix_3pt.jl")

function test_D2_3pt()

    D2 = build_D2_matrix_3pt(5, 1.0)
    display(D2); println()

    println("size = ", Base.summarysize(D2))

    D2_sp = sparse(D2)
    println("size = ", Base.summarysize(D2_sp))

    println("nzval  = ", D2_sp.nzval)
    println("rowval = ", D2_sp.rowval)
    println("colptr = ", D2_sp.colptr)
end

test_D2_3pt()
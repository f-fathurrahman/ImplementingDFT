using SparseArrays
using Printf

function test_01()
    A = zeros(5,4)
    A[1,:] = [0.0, 0.0, 1.0, 0.0]
    A[2,:] = [5.0, 8.0, 0.0, 0.0]
    A[3,:] = [0.0, 0.0, 3.0, 0.0]
    A[4,:] = [0.0, 6.0, 0.0, 0.0]

    display(A); println()

    sp_A = sparse(A)
    println("nzval  = ", sp_A.nzval)
    println("rowval = ", sp_A.rowval)
    println("colptr = ", sp_A.colptr)
end
#test_01()

function test_02()
    A = zeros(6,4)
    A[1,:] = [0.0, 0.0, 1.0, 5.0]
    A[2,:] = [5.0, 0.0, 0.0, 0.0]
    A[3,:] = [0.0, 0.0, 3.0, 0.0]
    A[4,:] = [0.0, 0.0, 0.0, 0.0]
    A[5,:] = [0.0, 0.0, 1.0, 0.0]
    A[6,:] = [5.0, 0.0, 0.0, 6.0]

    sp_A = sparse(A)
    println("nzval  = ", sp_A.nzval)
    println("rowval = ", sp_A.rowval)
    println("colptr = ", sp_A.colptr)

end
test_02()

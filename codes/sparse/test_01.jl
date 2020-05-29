using SparseArrays
using Printf

function main()
    A = zeros(4,4)
    A[1,:] = [0.0, 0.0, 0.0, 0.0]
    A[2,:] = [5.0, 8.0, 0.0, 0.0]
    A[3,:] = [0.0, 0.0, 3.0, 0.0]
    A[4,:] = [0.0, 6.0, 0.0, 0.0]

    sp_A = sparse(A)

end
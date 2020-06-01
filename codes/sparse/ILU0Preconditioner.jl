# Poor man's preconditioner
# Adapted from SPARSKIT2
# Originally written by Prof. Yousef Saad

struct ILU0Preconditioner
    alu::Vector{Float64}
    jlu::Vector{Int64}
    ju::Vector{Int64}
end

function ILU0Preconditioner( M )

    n  = M.n
    a  = M.nzval
    ja = M.rowval
    ia = M.colptr

    Nnzval = length(a)
    alu = zeros( Float64, Nnzval+1 )
    jlu = zeros( Int64, Nnzval+1 )
    ju  = zeros( Int64, n )
    
    # Work array
    iw = zeros( Int64, n )

    SMALL = 1e-10

    ju0 = n + 2
    jlu[1] = ju0

    for ii in 1:n
        
        js = ju0
        
        # generating row number ii of L and U

        for j in ia[ii]:ia[ii+1]-1
            
            # copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
            jcol = ja[j]
            
            if jcol == ii
               alu[ii] = a[j]
               iw[jcol] = ii
               ju[ii]  = ju0
            else
               alu[ju0] = a[j]
               jlu[ju0] = ja[j]
               iw[jcol] = ju0
               ju0 = ju0 + 1
            end
        end

        jlu[ii+1] = ju0
        jf = ju0 - 1
        jm = ju[ii] - 1
        
        # exit if diagonal element is reached.
        for j in js:jm
            jrow = jlu[j]
            tl = alu[j]*alu[jrow]
            alu[j] = tl
            # perform  linear combination
            for jj in ju[jrow]:jlu[jrow+1]-1
                jw = iw[ jlu[jj] ]
                if jw != 0
                    alu[jw] = alu[jw] - tl*alu[jj]
                end
            end
        end 

        # invert  and store diagonal element.
        if abs(alu[ii]) <= SMALL
            println("ii = ", ii)
            println("alu[ii] = ", alu[ii])
            # zero pivot :
            #600 ierr = ii
            error("Small pivot element")
        end
        
        alu[ii] = 1.0/alu[ii]

        # reset pointer iw to zero
        iw[ii] = 0
        for i in js:jf
            iw[jlu[i]] = 0
        end

    end
      
    return ILU0Preconditioner( alu, jlu, ju )

end

import LinearAlgebra: ldiv!
function ldiv!( prec::ILU0Preconditioner, v )

    jlu = prec.jlu
    ju = prec.ju
    alu = prec.alu

    n = length(prec.ju)

    for i in 1:n
        for k in jlu[i]:ju[i]-1
            v[i] = v[i] - alu[k] * v[jlu[k]]
        end
    end

    # backward solve.
    for i = n:-1:1
        for k in ju[i]:jlu[i+1]-1
            v[i] = v[i] - alu[k] * v[jlu[k]]
        end
        v[i] = alu[i]*v[i]
    end
    return

end

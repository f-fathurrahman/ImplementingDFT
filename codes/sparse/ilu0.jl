function init_ilu0!( n, a, ja, ia, alu, jlu, ju, iw )

# Input:
# - n: matrix dimension
# - a: nonzero values
# - rowval:
# - colptr:

# Output:
#  ALLOCATE( alu_ilu0(Npoints*(Nx+Ny+Nz-2)) ) <= real
#  ALLOCATE( jlu_ilu0(Npoints*(Nx+Ny+Nz-2)) ) <= integer
#  ALLOCATE( ju_ilu0(Npoints) )
#  ALLOCATE( iw_ilu0(Npoints) )

#  CALL ilu0( Npoints, -0.5d0*nzval, rowval, colptr, alu_ilu0, jlu_ilu0, ju_ilu0, iw_ilu0, ierr )

    ju0 = n + 2
    jlu[1] = ju0
    
    # initialize work vector to zero's
    for i in 1:n
      iw[i] = 0
    end

    # main loop
    for ii in 1:n
        
        js = ju0
        
        # generating row number ii of L and U.
        for j in ia[ii]:ia[ii+1]-1
            
            #copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
            jcol = ja[j]
            
            if jcol == ii
               alu[ii] = a[j]
               iw[jcol] = ii
               ju[ii]  = ju0
            else
               alu[ju0] = a[j]
               jlu[ju0] = ja[j]
               iw[jcol] = ju0
               ju0 = ju0+1
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
                if jw .ne. 0
                    alu[jw] = alu[jw] - tl*alu[jj]
                end
            end
        end 

        # invert  and store diagonal element.
        if alu[ii] <= SMALL
            # zero pivot :
            #600 ierr = ii
            error("Small pivot element")
        end
        
        alu[ii] = 1.0/alu[ii]

        # reset pointer iw to zero
        iw[ii] = 0
        do i in js:jf
            iw[jlu[i]] = 0
        end

    end
      
    return

end

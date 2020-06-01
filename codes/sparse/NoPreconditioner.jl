struct NoPreconditioner end

function NoPreconditioner(M)
    return NoPreconditioner()
end

import LinearAlgebra: ldiv!
function ldiv!( prec::NoPreconditioner, v )
    return
end
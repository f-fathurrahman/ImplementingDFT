# Only works for orthogonal lattice vectors
function calc_dr_periodic!( LL, r, r0, dr )

  xx1    = abs( r[1] - r0[1] )
  xx2    = abs( r[1] - r0[1] + LL[1] )
  xx3    = abs( r[1] - r0[1] - LL[1] )
  dr[1]  = min( xx1, xx2, xx3 )
    
  yy1    = abs( r[2] - r0[2] )
  yy2    = abs( r[2] - r0[2] + LL[2] )
  yy3    = abs( r[2] - r0[2] - LL[2] )
  dr[2]  = min( yy1, yy2, yy3 )
    
  zz1    = abs( r[3] - r0[3] )
  zz2    = abs( r[3] - r0[3] + LL[3] )
  zz3    = abs( r[3] - r0[3] - LL[3] )
  dr[3]  = min( zz1, zz2, zz3 )
    
  return
end
import qm3.mol
import qm3.engines._colvar_v
try:
    import cStringIO as io
except:
    import io
import qm3.maths.matrix

f = io.StringIO( """6

O       0.12109      0.06944     -0.22458
H      -0.52694      0.16499     -0.94583
H       0.70159     -0.63565     -0.54677
O      -0.45114      1.12675      2.21102
H      -0.29157      0.59483      1.39876
H       0.05804      1.92714      2.01036
""" )
m = qm3.mol.molecule()
m.xyz_read( f )
f.close()
m.chrg = [ -0.834, 0.417, 0.417, -0.834, 0.417, 0.417 ]
m.func = 0.0
m.grad = [ 0.0 for i in range( 3 * m.natm ) ]
o = qm3.engines._colvar_v.coulomb( 0.01, 102., [ 0 ], [ 3, 4, 5 ], [ 1 ], [] )
print( o.get_grad( m ) )
print( m.func )
qm3.maths.matrix.mprint( m.grad, m.natm, 3 )

m.func = 0.0
m.grad = [ 0.0 for i in range( 3 * m.natm ) ]
o = qm3.engines._colvar_v.fswitch( 0.01, 102., 12., 14., [ 0 ], [ 3, 4, 5 ], [ 1 ], [] )
print( o.get_grad( m ) )
print( m.func )
qm3.maths.matrix.mprint( m.grad, m.natm, 3 )

import qm3.mol
import qm3.engines.mmint
import qm3.maths.matrix
try:
    import cStringIO as io
except:
    import io

f = io.StringIO( """6

O       0.12109      0.06944     -0.22458
H      -0.52694      0.16499     -0.94583
H       0.70159     -0.63565     -0.54677
O      -0.45114      1.12675      2.21102
H      -0.29157      0.59483      1.39876
H       0.05804      1.92714      2.01036
""" )
mol = qm3.mol.molecule()
mol.xyz_read( f )
f.close()
mol.chrg = [ -0.834, 0.417, 0.417, -0.834, 0.417, 0.417 ]
mol.anum = [ 8, 1, 1, 8, 1, 1 ]

#mol.epsi = [ 0.9374, 0.0000, 0.0000, 0.9374, 0.000, 0.000 ] # Sqrt( eps * 4.184 ) >> Sqrt(kJ/mol)
#mol.rmin = [ 1.6610, 0.0000, 0.0000, 1.6610, 0.000, 0.000 ] # rmin/2 >> Ang
mol.type = [ "OW", "HW", "HW", "OW", "HW", "HW" ]
f = io.StringIO( """
OW    0.21002    1.661
HW    0.00000    0.000
""" )
out = qm3.engines.mmint.non_bonded( mol, f )
print( list( zip( mol.epsi, mol.rmin ) ) )

obj = qm3.engines.mmint.QMLJ_MMEL( mol, [ 0, 1, 2 ], [ 3, 4, 5 ], [] )
mol.grad = [ 0.0 for i in range( 3 * mol.natm ) ]
mol.hess = [ 0.0 for i in range( 9 * 9 ) ]
obj.get_hess( mol )
qm3.maths.matrix.mprint( mol.grad, mol.natm, 3 )
qm3.maths.matrix.mprint( mol.hess, 9, 9 )

obj = qm3.engines.mmint.QMLJ_MMEL( mol, [ 0, 1, 2 ], [ 3, 4, 5 ], [ [0,3,0.5] ] )
mol.grad = [ 0.0 for i in range( 3 * mol.natm ) ]
mol.hess = [ 0.0 for i in range( 9 * 9 ) ]
obj.get_hess( mol )
qm3.maths.matrix.mprint( mol.grad, mol.natm, 3 )
qm3.maths.matrix.mprint( mol.hess, 9, 9 )

obj = qm3.engines.mmint.QMLJ( mol, [ 0, 1, 2 ], [ 3, 4, 5 ], [] )
mol.grad = [ 0.0 for i in range( 3 * mol.natm ) ]
mol.hess = [ 0.0 for i in range( 9 * 9 ) ]
obj.get_hess( mol )
qm3.maths.matrix.mprint( mol.hess, 9, 9 )
mol.hess = []
dsp = 0.001
for i in range( 9 ):
    crd = mol.coor[i]
    mol.coor[i] = crd + dsp
    mol.grad = [ 0.0 for j in range( 3 * mol.natm ) ]
    obj.get_grad( mol )
    gpd = mol.grad[0:9]
    mol.coor[i] = crd - dsp
    mol.grad = [ 0.0 for j in range( 3 * mol.natm ) ]
    obj.get_grad( mol )
    gmd = mol.grad[0:9]
    mol.coor[i] = crd
    mol.hess += [ ( gpd[j] - gmd[j] ) / ( 2.0 * dsp ) for j in range( 9 ) ]
qm3.maths.matrix.mprint( mol.hess, 9, 9 )

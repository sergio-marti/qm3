import qm3.mol
import qm3.engines.mmres
try:
    import cStringIO as io
except:
    import io
import qm3.maths.matrix

f = io.StringIO( """4

C           2.0462659688       -0.7390557351       -1.2415870721
C           3.5769772508       -0.7978398171       -1.3291287529
C           4.0651517497       -1.6344633400       -2.5498538145
C           5.5719572779       -1.6411989523       -2.5768153041
""" )
m = qm3.mol.molecule()
m.xyz_read( f )
f.close()
r = qm3.engines.mmres.distance( 222.5 * 4.184 * 2.0, 1.53, [ 0, 1 ] )
m.func = 0.0
m.grad = [ 0.0 for i in range( 3 * m.natm ) ]
m.hess = [ 0.0 for i in range( 9 * m.natm * m.natm ) ]
# -- analytical hessian (set hessian indexes of the restrained atoms)
r.hind = [ 0, 1 ]
r.get_hess( m )
qm3.maths.matrix.mprint( m.grad, m.natm, 3 )
qm3.maths.matrix.mprint( m.hess, 3 * m.natm, 3 * m.natm )
# -- central difference
nhes = []
disp = 1.e-4
for i in range( 3 * m.natm ):
    bak = m.coor[i]
    m.grad = [ 0.0 for jj in range( 3 * m.natm ) ]
    m.coor[i] = bak + disp
    r.get_grad( m )
    gf = m.grad[:]
    m.grad = [ 0.0 for jj in range( 3 * m.natm ) ]
    m.coor[i] = bak - disp
    r.get_grad( m )
    gb = m.grad[:]
    m.coor[i] = bak
    nhes += [ (ii-jj) / ( 2.0 * disp ) for ii,jj in zip( gf, gb ) ]
qm3.maths.matrix.mprint( nhes, 3 * m.natm, 3 * m.natm )


r = qm3.engines.mmres.angle( 58.35 * 4.184 * 2.0, 113.6, [ 0, 1, 2 ] )
m.func = 0.0
m.grad = [ 0.0 for i in range( 3 * m.natm ) ]
m.hess = [ 0.0 for i in range( 9 * m.natm * m.natm ) ]
r.hind = [ 0, 1, 2 ]
r.get_hess( m )
qm3.maths.matrix.mprint( m.hess, 3 * m.natm, 3 * m.natm )
# -- central difference
nhes = []
disp = 1.e-4
for i in range( 3 * m.natm ):
    bak = m.coor[i]
    m.grad = [ 0.0 for jj in range( 3 * m.natm ) ]
    m.coor[i] = bak + disp
    r.get_grad( m )
    gf = m.grad[:]
    m.grad = [ 0.0 for jj in range( 3 * m.natm ) ]
    m.coor[i] = bak - disp
    r.get_grad( m )
    gb = m.grad[:]
    m.coor[i] = bak
    nhes += [ (ii-jj) / ( 2.0 * disp ) for ii,jj in zip( gf, gb ) ]
qm3.maths.matrix.mprint( nhes, 3 * m.natm, 3 * m.natm )


r = qm3.engines.mmres.dihedral( { 3: [ 0.195 * 4.184, 0.0 ] }, [ 0, 1, 2, 3 ] )
m.func = 0.0
m.grad = [ 0.0 for i in range( 3 * m.natm ) ]
m.hess = [ 0.0 for i in range( 9 * m.natm * m.natm ) ]
r.hind = [ 0, 1, 2, 3 ]
r.get_hess( m )
qm3.maths.matrix.mprint( m.hess, 3 * m.natm, 3 * m.natm )
# -- central difference
nhes = []
disp = 1.e-4
for i in range( 3 * m.natm ):
    bak = m.coor[i]
    m.grad = [ 0.0 for jj in range( 3 * m.natm ) ]
    m.coor[i] = bak + disp
    r.get_grad( m )
    gf = m.grad[:]
    m.grad = [ 0.0 for jj in range( 3 * m.natm ) ]
    m.coor[i] = bak - disp
    r.get_grad( m )
    gb = m.grad[:]
    m.coor[i] = bak
    nhes += [ (ii-jj) / ( 2.0 * disp ) for ii,jj in zip( gf, gb ) ]
qm3.maths.matrix.mprint( nhes, 3 * m.natm, 3 * m.natm )

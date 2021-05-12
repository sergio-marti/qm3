#!/usr/bin/env python3
import  pickle
import  qm3.mol
import  qm3.engines

f = open( "molec.pk", "rb" )
m = pickle.load( f )
f.close()
m.coor = [ .0 for i in range( 3 * m.natm ) ]
m.xyz_read( "xyz", replace = True )
x  = []
x += list( m.indx["A"][244].values() )
x += list( m.indx["A"][433].values() )
x += list( m.indx["A"][434].values() )
x += list( m.indx["A"][435].values() )
x += list( m.indx["A"][436].values() )
x += list( m.indx["A"][437].values() )
x += list( m.indx["A"][438].values() )
x += list( m.indx["A"][439].values() )
x += list( m.indx["A"][440].values() )
x += list( m.indx["A"][441].values() )
x.sort()

f = open( "sele_QM.pk", "wb" )
pickle.dump( x, f )
f.close()
f = open( "bonds.pk", "rb" )
b = pickle.load( f )
f.close()
qm3.engines.exclusions( x, m, b )
f = open( "sele_LA.pk", "rb" )
l = pickle.load( f )
f.close()

box = [ m.coor[0:3][:], m.coor[0:3][:] ]
cen = [ .0, .0, .0 ]
for i in range( m.natm ):
    i3 = i * 3
    for j in [0, 1, 2 ]:
        box[0][j] = min( box[0][j], m.coor[i3+j] )
        box[1][j] = max( box[1][j], m.coor[i3+j] )
        cen[j] += m.coor[i3+j]
cen = [ i / m.natm for i in cen ]
print( cen )
m.boxl = [ i-j for i,j in zip( box[1], box[0] ) ]
print( m.boxl )
for i in range( m.natm ):
    i3 = i * 3
    for j in [0, 1, 2]:
        m.coor[i3+j] -= cen[j]

y = m.sph_sel( x, 20. )
#y = sorted( set( y ).difference( set( x + [ j for i,j in l ] ) ) )
y = sorted( set( y ).difference( set( x ) ) )
f = open( "sele_MM.pk", "wb" )
pickle.dump( y, f )
f.close()

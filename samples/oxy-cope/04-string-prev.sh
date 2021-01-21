py_exe=python3

$py_exe << EOD
import  qm3.mol
import  qm3.utils
import  qm3.maths.interpolation
import  matplotlib.pyplot as plt
import  math
import  glob
import  os

bnds = [ ( "A", 1, "C2",  "A", 1, "C3" ),
         ( "A", 1, "C5",  "A", 1, "C6" ) ]

D = len( bnds )
W = 64
g = open( "04.string.def", "wt" )
g.write( "%d %d\n"%( D, W ) )
f = True
T = []
C = sorted( glob.glob( "03.neb.??.pdb" ) )
for fn in C:
    print( fn )
    m = qm3.mol.molecule( fn )
    T.append( [] )
    for i1,j1,k1, i2,j2,k2 in bnds:
        a1 = m.indx[i1][j1][k1]
        a2 = m.indx[i2][j2][k2]
        if( f ):
            g.write( "dist%10d%10d%10.1lf\n"%( a1, a2, 6000. ) )
        T[-1].append( qm3.utils.distance( m.coor[3*a1:3*a1+3], m.coor[3*a2:3*a2+3] ) )
    f = False

S = []
N = len( T )
x = list( range( N ) )
X = [ ( N - 1 ) / ( W - 1 ) * i for i in range( W ) ]
plt.clf()
plt.grid( True )
for i in range( D ):
    y = [ T[j][i] for j in range( N ) ]
    plt.plot( x, y, 'o' )
    o = qm3.maths.interpolation.hermite_spline( x, y )
    Y = [ o.calc( i )[0] for i in X ]
    plt.plot( X, Y, '.-' )
    S.append( Y[:] )
plt.savefig( "04.strdef.pdf" )

for i in range( W ):
    for j in range( D ):
        g.write( "%12.3lf"%( S[j][i] ) )
    g.write( "\n" )
g.close()

for i in range( W ):
    tmp = []
    for j in range( N ):
        tmp.append( ( sum( [ math.fabs( S[k][i] - T[j][k] ) for k in range( D ) ] ), C[j] ) )
    tmp.sort()
    os.symlink( tmp[0][1], "04.str_seed.%02d.pdb"%( i ) )
EOD

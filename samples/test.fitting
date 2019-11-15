import math
import qm3.actions.fitting

def polyn( x, c ):
    return( sum( [ c[i]*math.pow(x,i) for i in range( len( c ) ) ] ) )

print( "Generic Least Squares Fitting:" )
o = qm3.actions.fitting.problem( [1.0,2,3,4,5] , [2.0,5,7,9,8], polyn, [1.0, 1.0, 1.0, 1.0] )
o.fit()
o.table()

print( "\n\nPolynomial Regression:" )
r, c, s = qm3.actions.fitting.poly_fit( [1.0,2,3,4,5] , [2.0,5,7,9,8], 3 )
print( "R^2: ", r )
print( "[a0, a1, ..., ak] : " )
print( c )
print( "[s0, s1, ..., sk] : " )
print( s )

y = [ 7, 7, 5, 4, 2 ]
x = [ [ 7, 4, 10, 16, 13 ], [ 7, 3, 5, 7, 3 ] ]
print( "\n\nMultiple Linear Regression:" )
r, c, s = qm3.actions.fitting.MLR( x, y )
print( "R^2: ", r )
print( "[a0, a1, ..., ak] : " )
print( c )
print( "[s0, s1, ..., sk] : " )
print( s )

print( "\n\nPartial Least Squares:" )
x = [ [ 7, 4, 10, 16, 13 ], [ 7, 3, 5, 7, 3 ], [ 13, 14, 12, 11, 10 ], [ 7, 7, 5, 3, 3 ] ]
r, c, mx, my = qm3.actions.fitting.PLS( x, y )
print( "R^2: ", r )
print( "[a1, ..., ak] : " )
print( c )
print( "[<x1>, ..., <xk>] : " )
print( mx )
print( "<y>: ", my )

import random
random.seed()
n = 64
x = [ 0 + 0.1 * i for i in range( n ) ]
y = [ math.sin( i ) for i in x ]
p = [ i + random.random() - 0.5 for i in y ]
f, g = qm3.actions.fitting.savitzky_golay( x, p )
try:
    import matplotlib.pyplot as plt
    plt.clf()
    plt.grid( True )
    plt.plot( x, y, '.' )
    plt.plot( x, p, 'o' )
    plt.plot( x, f, '-' )
    plt.show()
except:
    for i in range( n ):
        print( "%20.10lf%20.10lf%20.10lf"%( x[i], y[i], f[i] ) )

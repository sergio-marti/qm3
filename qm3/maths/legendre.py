# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math



# legendre( 0, x ) == 1.
# legendre( 1, x ) == x
def legendre_f( n, x ):
    f = [ 1.0, x ]
    for i in range( 2, n + 1 ):
        y = float( 2.0 * ( i - 1.0 ) + 1.0 ) / float( i ) * x * f[1] - float( i - 1.0 ) / float( i ) * f[0]
        del f[0]
        f.append( y )
    return( f )



def legendre_g( n, x ):
    f = legendre_f( n, x )
    return( float( n ) / ( 1.0 - x * x ) * ( f[0] - x * f[1] ) )



def legendre_h( n, x ):
    f = legendre_f( n, x )
    g = float( n ) / ( 1.0 - x * x ) * ( f[0] - x * f[1] )
    return( ( 2.0 * x * g - float( n * ( n + 1.0 ) ) * f[1] ) / ( 1.0 - x * x ) )



# Gauss-Legendre Quadrature:
# int_a^b{ f(x) dx } = {b - a} over {2} sum_{ -1 }^{ {}+1 }{ w_i · f left( {{b - a} over 2 } · x_i + {{b +a} over 2} right) }
#
def GL_coefficients( n, eps = 1.0e-8 ):
    xi = []
    wi = []
    for i in range( 1, n + 1 ):
        x = ( 1.0 - float( i - 1.0 ) / float( 8.0 * n * n * n ) ) * math.cos( math.pi * float( 4.0 * i - 1.0 ) / float( 4.0 * n + 2.0 ) )
        f = legendre_f( n, x )
        g = float( n ) / ( 1.0 - x * x ) * ( f[0] - x * f[1] )
        h = ( 2.0 * x * g - float( n * ( n + 1.0 ) ) * f[1] ) / ( 1.0 - x * x )
        j = 1
        while( j < 1000 and math.fabs( f[1] ) > eps ):
            x -= 2.0 * f[1] * g / ( 2.0 * g * g - f[1] * h )
            f = legendre_f( n, x )
            g = float( n ) / ( 1.0 - x * x ) * ( f[0] - x * f[1] )
            h = ( 2.0 * x * g - float( n * ( n + 1.0 ) ) * f[1] ) / ( 1.0 - x * x )
            j += 1
        xi.append( x )
        wi.append( 2.0 * ( 1.0 - x * x ) / ( float( n ) * float( n ) * f[0] * f[0] ) )
    return( xi, wi )
        



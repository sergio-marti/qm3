# -*- coding: iso-8859-1 -*-

# https://en.wikipedia.org/wiki/Euler_method
# https://en.wikipedia.org/wiki/Numerov%27s_method
# https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
# https://en.wikipedia.org/wiki/Midpoint_method

import math

try:
    import qm3.maths._ode
except:
    pass



# --------------------------------------
#
#    dx1 / dt = f1( t, x1, ..., xn )
#    ...
#    dxn / dt = fn( t, x1, ..., xn )
#
#    function returns a list with the value of the f1, ..., fn
#

def Euler_init( function, dt, t, x ):
    return( [ dt * i + j for i,j in zip( function( t, x ), x ) ] )


def Euler_step( function, dt, t, x, xo ):
    return( [ 2.0 * dt * i + j for i,j in zip( function( t, x ), xo ) ] )



def RK4_step( function, dt, t, x ):
    k1 = [ dt * i for i in function( t, x ) ]
    dx = [ i + 0.5 * j for i,j in zip( x, k1 ) ]
    k2 = [ dt * i for i in function( t + dt * 0.5, dx ) ]
    dx = [ i + 0.5 * j for i,j in zip( x, k2 ) ]
    k3 = [ dt * i for i in function( t + dt * 0.5, dx ) ]
    dx = [ i + j for i,j in zip( x, k3 ) ]
    k4 = [ dt * i for i in function( t + dt, dx ) ]
    nx = [ (i+(j+2.0*k+2.0*l+m)/6.0) for i,j,k,l,m in zip( x, k1, k2, k3, k4 ) ]
    return( nx )


# ----------------------------------------------------------------------
#
#    d^2x1 / dt^2 = f1( t, x1, ..., xn, v1, ..., vn )    vi = dxi / dt
#    ...
#    d^2xn / dt^2 = fn( t, x1, ..., xn, v1, ..., vn )
#
#    function returns a list with the value of the f1, ..., fn
#

def RK42_step( function, dt, t, x, v ):
    k1 = [ dt * i for i in v ]
    l1 = [ dt * i for i in function( t, x, v ) ]
    dx = [ i + 0.5 * j for i,j in zip( x, k1 ) ]
    dv = [ i + 0.5 * j for i,j in zip( v, l1 ) ]
    k2 = [ dt * i for i in dv ]
    l2 = [ dt * i for i in function( t + dt * 0.5, dx, dv ) ]
    dx = [ i + 0.5 * j for i,j in zip( x, k2 ) ]
    dv = [ i + 0.5 * j for i,j in zip( v, l2 ) ]
    k3 = [ dt * i for i in dv ]
    l3 = [ dt * i for i in function( t + dt * 0.5, dx, dv ) ]
    dx = [ i + j for i,j in zip( x, k3 ) ]
    dv = [ i + j for i,j in zip( v, l3 ) ]
    k4 = [ dt * i for i in dv ]
    l4 = [ dt * i for i in function( t + dt, dx, dv ) ]
    nx = [ (i+(j+2.0*k+2.0*l+m)/6.0) for i,j,k,l,m in zip( x, k1, k2, k3, k4 ) ]
    nv = [ (i+(j+2.0*k+2.0*l+m)/6.0) for i,j,k,l,m in zip( v, l1, l2, l3, l4 ) ]
    return( nx, nv )



def least_squares_finite_elements_2d( nx, dx, ny, dy, grd, step_size = 10.0, deltaQ = -1 ):
    try:
        return( qm3.maths._ode.least_squares_finite_elements_2d( nx, dx, ny, dy, grd, step_size, deltaQ ) )
    except:
        raise Exception( "ode.least_squares_finite_elements_2d: qm3.maths._ode.so not available..." )

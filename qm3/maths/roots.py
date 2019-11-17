# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.maths.matrix



# -- SINGLE FUNCTION of SINGLE VAR
def bisect( function, x0, xf, max_iter = 1000, eps = 1.0e-10 ):
    l0 = x0
    lf = xf
    f0 = function( l0 )
    ff = function( lf )
    ni = 0
    while( ni < max_iter and math.fabs( l0 - lf ) > eps ):
        lm = ( l0 + lf ) * 0.5
        fm = function( lm )
        if( f0 * fm <= 0.0 ):
            lf = lm
            ff = fm
        else:
            l0 = lm
            f0 = fm
        ni += 1
    if( ni >= max_iter ):
        return( None )
    return( lm )



# -- SINGLE FUNCTION of SINGLE VAR
def ridders( function, x0, xf, max_iter = 1000, eps = 1.0e-10 ):
    f0 = function( x0 )
    if( math.fabs( f0 ) <= eps ):
        return( x0 )
    ff = function( xf )
    if( math.fabs( ff ) <= eps ):
        return( xf )
    ni = 0
    fn = eps * 2.0
    while( ni < max_iter and math.fabs( fn ) > eps ):
        xm = 0.5 * ( x0 + xf )
        fm = function( xm )
        s = math.sqrt( fm * fm - f0 * ff )
        if( s == 0.0 ):
            return( None )
        d = ( xm - x0 ) * fm / s
        if( f0 - ff < 0.0 ):
            d = - d
        xn = xm + d
        fn = function( xn )
        if( math.fabs( fn ) <= eps ):
            return( xn )
        if( math.copysign( 1.0, fm ) == math.copysign( 1.0, fn ) ):
            if( math.copysign( 1.0, f0 ) == math.copysign( 1.0, fn ) ):
                xf = xn
                ff = fn
            else:
                x0 = xn
                f0 = fn
        else:
            x0 = xm
            f0 = fm
            xf = xn
            ff = fn
        ni += 1
    if( ni >= max_iter ):
        return( None )
    return( xn )



# -- SINGLE FUNCTION of SINGLE VAR (NUMERICAL GRADIENT)
def newton_raphson( function, x0, max_iter = 1000, max_stp = 1.0, eps = 1.0e-10, dsp = 1.0e-4 ):
    xc = x0
    fx = function( xc )
    ni = 0
    while( ni < max_iter and math.fabs( fx ) > eps ):
        gx = ( function( xc + dsp ) - fx ) / dsp
        dx = min( max_stp, math.fabs( fx / gx ) ) * (fx/gx) * math.fabs(gx/fx)
        xc -= dx
        fx = function( xc )
        ni += 1
    if( ni >= max_iter ):
        return( None )
    return( xc )
    


# -- SINGLE FUNCTION of SINGLE VAR (NUMERICAL HESSIAN)
def halley( function, gradient, x0, max_iter = 1000, max_stp = 1.0, eps = 1.0e-10, dsp = 1.0e-4 ):
    xc = x0
    fx = function( xc )
    ni = 0
    while( ni < max_iter and math.fabs( fx ) > eps ):
        gx = gradient( xc )
        hx = ( gradient( xc + dsp ) - gx ) / dsp
        dx = 2.0 * fx * gx / ( 2.0 * gx * gx - fx * hx )
#        dx = fx / ( gx - fx * hx / gx * 0.5 )
        if( dx == 0.0 ):
            dx = max_stp
        dx = min( max_stp, math.fabs( dx ) ) * dx / math.fabs( dx )
        xc -= dx
        fx = function( xc )
        ni += 1
    if( ni >= max_iter ):
        return( None )
    return( xc )
    


# --  n FUNCTIONS of n VARS (NUMERICAL GRADIENT)
def multi_newton_raphson( function, x0, max_iter = 1000, max_stp = 1.0, eps = 1.0e-10, dsp = 1.0e-4 ):
    nv = len( x0 )
    xc = x0[:]
    fx = [ function[i]( xc ) for i in range( nv ) ]
    mx = max( [ math.fabs( i ) for i in fx ] )
    ni = 0
    while( ni < max_iter and mx > eps ):
        gx = []
        for i in range( nv ):
            for j in range( nv ):
                xc[j] += dsp
                gx.append( ( function[i]( xc ) - fx[i] ) / dsp )
                xc[j] -= dsp
        dx = qm3.maths.matrix.mult( qm3.maths.matrix.inverse( gx, nv, nv ), nv, nv, fx, nv, 1 )
        mx = math.sqrt( sum( i*i for i in dx ) )
        if( mx > max_stp ):
            dx = [ i / mx * max_stp for i in dx ]
        xc = [ xc[i] - dx[i] for i in range( nv ) ]
        fx = [ function[i]( xc ) for i in range( nv ) ]
        mx = max( [ math.fabs( i ) for i in fx ] )
        ni += 1
    if( ni >= max_iter ):
        return( None )
    return( xc )
    




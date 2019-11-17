# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.maths.matrix
import qm3.constants
import qm3.elements
import os
import stat
import struct
try:
    import cPickle as pickle
except:
    import pickle

try:
    import qm3.utils._conn
    conn_so = True
except:
    conn_so = False


# ----------------------------------------------------------------------------------
# Logging default
#
def log_event( txt ):
    sys.stdout.write( txt + "\n" )
    sys.stdout.flush()



# ----------------------------------------------------------------------------------
# Basic Geometry
#
def distanceSQ( ci, cj, box = None ):
    if( box == None ):
        return( sum( [ (i-j)*(i-j) for i,j in zip( ci, cj ) ] ) )
    else:
        dr = [ ( ci[i] - cj[i] ) - box[i] * round( ( ci[i] - cj[i] ) / box[i], 0 ) for i in [0, 1, 2] ]
        return( dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2] )


def distance( ci, cj, box = None ):
    return( math.sqrt( distanceSQ( ci, cj, box ) ) )



def angleRAD( ci, cj, ck ):
    vij = [ (i-j) for i,j in zip( ci, cj ) ]
    mij = math.sqrt( sum( [ i*i for i in vij ] ) )
    vkj = [ (i-j) for i,j in zip( ck, cj ) ]
    mkj = math.sqrt( sum( [ i*i for i in vkj ] ) )
    try:
        return( math.acos( sum( [ i*j for i,j in zip( vij, vkj ) ] ) / ( mij * mkj ) )  )
    except:
        return( 0.0 )


def angle( ci, cj, ck ):
    return( angleRAD( ci, cj, ck ) * qm3.constants.R2D )



def dihedralRAD( ci, cj, ck, cl ):
    vij = [ (i-j) for i,j in zip( ci, cj ) ]
    vkj = [ (i-j) for i,j in zip( ck, cj ) ]
    vlj = [ (i-j) for i,j in zip( cl, cj ) ]
    pik = qm3.maths.matrix.cross_product( vij, vkj )
    plk = qm3.maths.matrix.cross_product( vlj, vkj )
    m1  = sum( [ (i*j) for i,j in zip( pik, plk ) ] )
    m2  = math.sqrt( sum( [ (i*j) for i,j in zip( pik, pik ) ] ) )
    m3  = math.sqrt( sum( [ (i*j) for i,j in zip( plk, plk ) ] ) )
    o   = 0.0
    if( m2 != 0.0 and m3 != 0.0 ):
        o = m1 / ( m2 * m3 )
        if( math.fabs( o ) > 1. ):
            o = math.fabs( o ) / o
        o = math.acos( o ) 
        if( sum( [ (i*j) for i,j in zip( vij, plk ) ] ) < 0.0 ):
            o = -o
    return( o )


def dihedral( ci, cj, ck, cl ):
    return( dihedralRAD( ci, cj, ck , cl ) * qm3.constants.R2D )



# ----------------------------------------------------------------------------------
# Geometric utilities
#
def center( mass, coor ):
    size = len( coor ) // 3
    mt = sum( mass )
    mc = [ 0.0, 0.0, 0.0 ]
    for i in range( size ):
        for j in [ 0, 1, 2 ]:
            mc[j] += coor[i*3+j] * mass[i]
    mc[0] /= mt; mc[1] /= mt; mc[2] /= mt
    for i in range( size ):
        i3 = i * 3
        for j in [0, 1, 2]:
            coor[i3+j] -= mc[j]
    return( mc )



def moments_of_inertia( mass, coor ):
    size = len( coor ) // 3
    xx = 0.0; xy = 0.0; xz = 0.0; yy = 0.0; yz = 0.0; zz = 0.0
    mc = center( mass, coor )
    for i in range( size ):
        i3 = i * 3
        xx += mass[i] * coor[i3]   * coor[i3]
        xy += mass[i] * coor[i3]   * coor[i3+1]
        xz += mass[i] * coor[i3]   * coor[i3+2]
        yy += mass[i] * coor[i3+1] * coor[i3+1]
        yz += mass[i] * coor[i3+1] * coor[i3+2]
        zz += mass[i] * coor[i3+2] * coor[i3+2]
    val, vec = qm3.maths.matrix.diag( qm3.maths.matrix.from_upper_diagonal_rows( [ yy + zz, -xy, -xz, xx + zz, -yz, xx + yy ], 3 ), 3 )
    if( qm3.maths.matrix.det( vec, 3 ) < 0.0 ):
        for i in [0, 3, 6]:
            vec[i] = - vec[i]
    for i in range( size ):
        i3 = i * 3
        t = qm3.maths.matrix.mult( coor[i3:i3+3], 1, 3, vec, 3, 3 )
        coor[i3:i3+3] = t[:]
    return( mc, vec )



def superimpose_quaternion( mass, coor_a, coor_b ):
    size = len( coor_a ) // 3
    t = coor_a[:]
    mc = center( mass, t )
    mx = center( mass, coor_b )
    m = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    for i in range( size ):
        i3 = i * 3
        c = [   2. * ( t[i3+1] * coor_b[i3+2] - t[i3+2] * coor_b[i3+1] ),
                2. * ( t[i3+2] * coor_b[i3  ] - t[i3  ] * coor_b[i3+2] ),
                2. * ( t[i3  ] * coor_b[i3+1] - t[i3+1] * coor_b[i3  ] ) ]
        f = sum( [ (j+k)*(j+k) for j,k in zip( t[i3:i3+3], coor_b[i3:i3+3] ) ] )
        m[0] += mass[i] * sum( [ (j-k)*(j-k) for j,k in zip( t[i3:i3+3], coor_b[i3:i3+3] ) ] )
        m[1] += mass[i] * c[0]
        m[2] += mass[i] * c[1]
        m[3] += mass[i] * c[2]
        m[4] += mass[i] * ( f - 4. * t[i3] * coor_b[i3] )
        m[5] -= mass[i] * 2. * ( t[i3] * coor_b[i3+1] + t[i3+1] * coor_b[i3] )
        m[6] -= mass[i] * 2. * ( t[i3] * coor_b[i3+2] + t[i3+2] * coor_b[i3] )
        m[7] += mass[i] * ( f - 4. * t[i3+1] * coor_b[i3+1] )
        m[8] -= mass[i] * 2. * ( t[i3+1] * coor_b[i3+2] + t[i3+2] * coor_b[i3+1] )
        m[9] += mass[i] * ( f - 4. * t[i3+2] * coor_b[i3+2] )
    val, vec = qm3.maths.matrix.diag( qm3.maths.matrix.from_upper_diagonal_rows( m, 4 ), 4 )
    q0 = vec[0]; q1 = vec[4]; q2 = vec[8]; q3 = vec[12]
    u = [ q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2. * (  q0 * q3 + q1 * q2 ), 2. * ( - q0 * q2 + q1 * q3 ),
            2. * ( - q0 * q3 + q1 * q2 ), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2. * (   q0 * q1 + q2 * q3 ),
            2. * (   q0 * q2 + q1 * q3 ), 2. * ( - q0 * q1 + q2 * q3 ), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3 ]
    if( math.fabs( qm3.maths.matrix.det( u, 3 ) - 1. ) > 1.0e-6 ):
        print( " - superimpose_quaternion: invalid rotation qm3.maths.matrix! " )
        return
    for i in range( size ):
        i3 = i * 3
        t = coor_b[i3:i3+3]
        coor_b[i3  ] = u[0] * t[0] + u[3] * t[1] + u[6] * t[2] + mc[0]
        coor_b[i3+1] = u[1] * t[0] + u[4] * t[1] + u[7] * t[2] + mc[1]
        coor_b[i3+2] = u[2] * t[0] + u[5] * t[1] + u[8] * t[2] + mc[2]
    return( mx, mc, u )



def superimpose_kabsch( mass, coor_a, coor_b ):
    size = len( coor_a ) // 3
    t = coor_a[:]
    mc = center( mass, t )
    mx = center( mass, coor_b )
    m = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    for i in [0, 1, 2]:
        for j in [0, 1, 2]:
            T = 0.0
            for k in range( size ):
                T += mass[k] * t[3*k+j] * coor_b[3*k+i]
            m[3*i+j] = T
    val, vec = qm3.maths.matrix.diag( qm3.maths.matrix.mult( qm3.maths.matrix.T( m, 3, 3 ), 3, 3, m , 3, 3 ), 3 )
    a = [ vec[2], vec[1], vec[5] * vec[7] - vec[4] * vec[8],
            vec[5], vec[4], vec[8] * vec[1] - vec[2] * vec[7],
            vec[8], vec[7], vec[2] * vec[4] - vec[1] * vec[5] ]
    ra = qm3.maths.matrix.mult( m, 3, 3, a, 3, 3 )
    b = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    t = math.sqrt( sum( [ ra[i] * ra[i] for i in [0, 3, 6] ] ) )
    b[0] = ra[0] / t; b[3] = ra[3] / t; b[6] = ra[6] / t
    t = math.sqrt( sum( [ ra[i] * ra[i] for i in [1, 4, 7] ] ) )
    b[1] = ra[1] / t; b[4] = ra[4] / t; b[7] = ra[7] / t
    b[2] = b[3] * b[7] - b[6] * b[4]
    b[5] = b[6] * b[1] - b[0] * b[7]
    b[8] = b[0] * b[4] - b[3] * b[1]
    t = sum( [ b[i]*ra[i] for i in [2, 5, 8] ] )
    if( t < 0.0 ):
        b[2] = - b[2]; b[5] = - b[5]; b[8] = - b[8]
    i = 0
    while( i < 2 ):
        u = qm3.maths.matrix.mult( b, 3, 3, qm3.maths.matrix.T( a, 3, 3 ), 3, 3 )
        if( math.fabs( qm3.maths.matrix.det( u, 3 ) - 1. ) > 1.0e-6 ):
            b[2] = - b[2]; b[5] = - b[5]; b[8] = - b[8]
            if( i == 1 ):
                print( " - superimpose_kabsch: invalid rotation qm3.maths.matrix! " )
                return
        else:
            i += 1
        i += 1
    for i in range( size ):
        i3 = i * 3
        t = coor_b[i3:i3+3]
        coor_b[i3  ] = u[0] * t[0] + u[3] * t[1] + u[6] * t[2] + mc[0]
        coor_b[i3+1] = u[1] * t[0] + u[4] * t[1] + u[7] * t[2] + mc[1]
        coor_b[i3+2] = u[2] * t[0] + u[5] * t[1] + u[8] * t[2] + mc[2]
    return( mx, mc, u )



# -----------------------------------------------------------------------------
# a = {1, 0, 0};
# b = {0, 1, 0};
# c = Dot[a, b]
# v = Cross[a, b]
# m = {{0, -v[[3]], v[[2]]}, {v[[3]], 0, -v[[1]]}, {-v[[2]], v[[1]], 0}};
# m // MatrixForm
# r = IdentityMatrix[3] + m + m.m/(1 + c);
# r // MatrixForm
# r.a
# -----------------------------------------------------------------------------
# -- coor_b should be initially centered...
def superimpose_vector( vec_a, vec_b, coor_b = None ):
    size = 0
    if( coor_b ):
        size = len( coor_b ) // 3
    nrm_a = qm3.maths.matrix.norm( vec_a )
    nrm_b = qm3.maths.matrix.norm( vec_b )
    fac_c = qm3.maths.matrix.dot_product( nrm_a, nrm_b )
    if( fac_c < 0.0 ):
        for i in range( size ):
            i3 = i * 3
            tmp = qm3.maths.matrix.dot_product( nrm_b, coor_b[i3:i3+3] )
            for j in [0, 1, 2]:
                coor_b[i3+j] = coor_b[i3+j] - 2. * nrm_b[j] * tmp
        if( fac_c == -1.0 ):
            return
        nrm_b = [ -i for i in nrm_b ]
        fac_c = qm3.maths.matrix.dot_product( nrm_a, nrm_b )
    fac_c = 1. / ( 1. + fac_c )
    vec_p = qm3.maths.matrix.cross_product( nrm_b, nrm_a )
    mat   = [ 1.-(vec_p[1]*vec_p[1]+vec_p[2]*vec_p[2])*fac_c, vec_p[0]*vec_p[1]*fac_c-vec_p[2], vec_p[0]*vec_p[2]*fac_c+vec_p[1], vec_p[0]*vec_p[1]*fac_c+vec_p[2], 1.-(vec_p[0]*vec_p[0]+vec_p[2]*vec_p[2])*fac_c, vec_p[1]*vec_p[2]*fac_c-vec_p[0], vec_p[0]*vec_p[2]*fac_c-vec_p[1], vec_p[1]*vec_p[2]*fac_c+vec_p[0], 1.-(vec_p[0]*vec_p[0]+vec_p[1]*vec_p[1])*fac_c ]
    for i in range( size ):
        i3 = i * 3
        t  = coor_b[i3:i3+3]
        coor_b[i3  ] = mat[0] * t[0] + mat[1] * t[1] + mat[2] * t[2]
        coor_b[i3+1] = mat[3] * t[0] + mat[4] * t[1] + mat[5] * t[2]
        coor_b[i3+2] = mat[6] * t[0] + mat[7] * t[1] + mat[8] * t[2]
#        coor_b[i3:i3+3] = qm3.maths.matrix.mult( mat, 3, 3, coor_b[i3:i3+3], 3, 1 )
    return( mat )



def rotate( coor, center, axis, theta ):
    size = len( coor ) // 3
    # clockwise rotation...
    cos = math.cos( - theta / qm3.constants.R2D )
    sin = math.sin( - theta / qm3.constants.R2D )
    m = math.sqrt( sum( [ i*i for i in axis ] ) )
    vr = [ i/m for i in axis ]
    if( vr[2] != 0.0 ):
        vp = [ 1., 1., -( vr[0] + vr[1] ) / vr[2] ]
    else:
        if( vr[1] != 0.0 ):
            vp = [ 1., -vr[0] / vr[1], 0.0 ]
        else:
            vp = [ 0.0, 1.0, 0.0 ]
    m = math.sqrt( sum( [ i*i for i in vp ] ) )
    vp = [ i/m for i in vp ]
    vu = [ vr[1] * vp[2] - vr[2] * vp [1], vr[2] * vp[0] - vr[0] * vp[2], vr[0] * vp[1] - vr[1] * vp[0] ]
    m = math.sqrt( sum( [ i*i for i in vu ] ) )
    vu = [ i/m for i in vu ]
    mcb = [ vp[0], vp[1], vp[2], vu[0], vu[1], vu[2], vr[0], vr[1], vr[2] ]
    mcbi = qm3.maths.matrix.inverse( mcb, 3, 3 )
#    mcbi = qm3.maths.matrix.T( mcb, 3, 3 )
    for i in range( size ):
        i3 = i * 3
        t = [ coor[i3+j] - center[j] for j in [0, 1, 2] ]
        t = qm3.maths.matrix.mult( mcb, 3, 3, t, 3, 1 )
        t = [ cos * t[0] + sin * t[1], cos * t[1] - sin * t[0], t[2] ]
        t = qm3.maths.matrix.mult( mcbi, 3, 3, t, 3, 1 )
        coor[i3:i3+3] = [ t[j] + center[j] for j in [0, 1, 2] ]



# ----------------------------------------------------------------------------------
# Normal Mode Utilities
#
def get_RT_modes( mass, coor ):
    size = len( coor )
    rt = []
    for i in range( 6 ):
        rt.append( [] )
    # use center of geometry instead
    if( mass == None ):
        mass = []
        for i in range( size // 3 ):
            mass.append( 1. )
        mt = float( size ) / 3.
    # use center of mass
    else:
        mt = sum( mass )
    mc = [ 0.0, 0.0, 0.0 ]
    for i in range( size // 3 ):
        for j in [ 0, 1, 2 ]:
            mc[j] += coor[i*3+j] * mass[i]
    mc[0] /= mt; mc[1] /= mt; mc[2] /= mt
    for i in range( size // 3 ):
        j = 3 * i
        mt = math.sqrt( mass[i] )
        rt[0] += [ mt, 0.0, 0.0 ]
        rt[1] += [ 0.0, mt, 0.0 ]
        rt[2] += [ 0.0, 0.0, mt ]
        rt[3] += [ 0.0, - ( coor[j+2] - mc[2] ) * mt, ( coor[j+1] - mc[1] ) * mt ]
        rt[4] += [ ( coor[j+2] - mc[2] ) * mt, 0.0, - ( coor[j  ] - mc[0] ) * mt ]
        rt[5] += [ - ( coor[j+1] - mc[1] ) * mt, ( coor[j  ] - mc[0] ) * mt, 0.0 ]
    # orthogonalize modes
    for i in range( 6 ):
        for j in range( i ):
            t = qm3.maths.matrix.dot_product( rt[i], rt[j] )
            for k in range( size ):
                rt[i][k] -= t * rt[j][k]
        t = math.sqrt( qm3.maths.matrix.dot_product( rt[i], rt[i] ) )
        if( t > 0.0 ):
            for k in range( size ):
                rt[i][k] /= t
    return( rt )



def project_RT_modes( mass, coor, grad, hess ): 
    size = len( coor )
    rt = get_RT_modes( mass, coor )
    # G' = G - Tx * G * Tx - ... - Rx * G * Rx - ...
    if( grad ):
        for i in range( 6 ):
            t = sum( [ ii * jj  for ii,jj in zip( grad, rt[i] ) ] )
            for j in range( size ):
                grad[j] -= t * rt[i][j]
    if( hess ):
        # P = I - Tx * Tx - ... - Rx * Rx - ... 
        ix = [ 0.0 for ii in range( size * size ) ]
        for i in range( size ):
            ix[size*i+i] += 1.
            for j in range( size ):
                for l in range( 6 ):
                    ix[size*i+j] -= rt[l][i] * rt[l][j]
        # H' = P * H * P
        t = qm3.maths.matrix.mult( ix, size, size, qm3.maths.matrix.mult( hess, size, size, ix, size, size ), size, size )
        for i in range( size * size ):
            hess[i] = t[i]



def raise_hessian_RT( mass, coor, hess, large = 5.0e4 ):
    size = len( coor )
    rt = get_RT_modes( mass, coor )
    # H' = H + large * Tx · Tx + ... + large * Rx · Rx + ...
    for i in range( size ):
        for j in range( size ):
            for k in range( 6 ):
                hess[i*size+j] += large * rt[k][i] * rt[k][j]
            


def hessian_frequencies( mass, coor, hess, project_RT = True ):
    size = len( coor )
    w = [ 1.0 / math.sqrt( mass[i] ) for i in range( size // 3 ) ]
    h = []
    for i in range( size ):
        for j in range( size ):
            h.append( hess[i*size+j] * w[i//3] * w[j//3] )
    if( project_RT ):
        project_RT_modes( mass, coor, None, h )
    freq, mods = qm3.maths.matrix.diag( h, size )
    wns = 1.0e11 / ( 2. * math.pi * qm3.constants.C )
    for i in range( size ):
        if( freq[i] < 0.0 ):
            freq[i] = - math.sqrt( math.fabs( freq[i] ) ) * wns
        else:
            freq[i] =   math.sqrt( math.fabs( freq[i] ) ) * wns
        # Mass-Weight the modes... (for vibrations >> normal_mode_view)
        for j in range( size ):
            mods[i*size+j] *= w[j//3]
    # freq: cm^-1, mods: 1/sqrt[g/mol]
    return( freq, mods )



def force_constants( mass, freq, mods ):
    size = len( freq )
    wns  = 1.0e11 / ( 2. * math.pi * qm3.constants.C )
    # eigenvalues: cm^-1 >> kJ/g.A^2
    val  = [ math.pow( freq[i] / wns, 2.0 ) for i in range( size ) ]
    # eigenvectors: dimensionless
    wei  = [ math.sqrt( mass[i] ) for i in range( size // 3 ) ]
    vec  = []
    for i in range( size ):
        for j in range( size ):
            vec.append( mods[i*size+j] * wei[j//3] )
    # reduced masses: g/mol
    rmas = [ 1.0 / sum( [ vec[i*size+j] * vec[i*size+j] / mass[i//3] for i in range( size ) ] ) for j in range( size ) ]
    # force qm3.constants: mDyne/A
    cte  = 1.e21 / qm3.constants.NA
    frce = [ cte * rmas[i] * val[i] for i in range( size ) ]
    return( rmas, frce )



def normal_mode_view( coor, freq, mods, symb, who, temp = 298.15, afac = 1. ):
    size = len( coor )
    ome = math.fabs( freq[who] )
    if( ome < 0.1 ): 
        return
    wns = 1.0e11 / ( 2. * math.pi * qm3.constants.C )
    amp = math.sqrt( 2. * 1.0e-3 * qm3.constants.KB * qm3.constants.NA * temp ) * ( wns / ome )
    fd = open( "nmode.%d"%( who ), "wt" )
    for i in range( 10 ):
        for j in range( 10 ):
            fd.write( "%d\n%10.3lf cm^-1\n"%( size / 3, freq[who] ) )
            fac = afac * amp * math.sin( 2. * math.pi * float(j) / 10. )
            for k in range( size // 3 ):
                fd.write( "%2s%8s%20.12lf%20.12lf%20.12lf\n"%( symb[k], "",
                      coor[3*k] + fac *   mods[who+size*(3*k)],
                    coor[3*k+1] + fac * mods[who+size*(3*k+1)],
                    coor[3*k+2] + fac * mods[who+size*(3*k+2)] ) )
    fd.close()



def gibbs_rrho( mass, coor, freq, temp = 298.15, press = 1.0, symm = 1.0, fcut = 10. ):
    # Translation (divided by N_Avogadro)
    mt = sum( mass ) 
    qt = math.pow( 2.0 * math.pi * mt / qm3.constants.NA * 1.0e-3 * qm3.constants.KB * temp / ( qm3.constants.H * qm3.constants.H ), 1.5 ) * qm3.constants.KB * temp / ( press * 1.013250E+5 )
    qt = math.log( qt )
    # Rotations
    mc = [ 0.0, 0.0, 0.0 ]
    for i in range( len( mass ) ):
        i3 = i * 3
        for j in [0, 1, 2]:
            mc[j] += mass[i] * coor[i3+j]
    mc[0] /= mt; mc[1] /= mt; mc[2] /= mt
    xx = 0.0; xy = 0.0; xz = 0.0; yy = 0.0; yz = 0.0; zz = 0.0
    for i in range( len( mass ) ):
        i3 = i * 3
        xx += mass[i] * ( coor[i3]   - mc[0] ) * ( coor[i3]   - mc[0] )
        xy += mass[i] * ( coor[i3]   - mc[0] ) * ( coor[i3+1] - mc[1] )
        xz += mass[i] * ( coor[i3]   - mc[0] ) * ( coor[i3+2] - mc[2] )
        yy += mass[i] * ( coor[i3+1] - mc[1] ) * ( coor[i3+1] - mc[1] )
        yz += mass[i] * ( coor[i3+1] - mc[1] ) * ( coor[i3+2] - mc[2] )
        zz += mass[i] * ( coor[i3+2] - mc[2] ) * ( coor[i3+2] - mc[2] )
    val, vec = qm3.maths.matrix.diag( qm3.maths.matrix.from_upper_diagonal_rows( [ yy + zz, -xy, -xz, xx + zz, -yz, xx + yy ], 3 ), 3 )
    t = ( 8.0 * math.pi * math.pi * qm3.constants.KB * temp ) / ( qm3.constants.H * qm3.constants.H * qm3.constants.NA ) * 1.0e-23
    qr = math.sqrt( math.pi * t * t * t * val[0] * val[1] * val[2] ) / symm
    qr = math.log( qr )
    # Vibrations
    t = 100.0 * qm3.constants.C * qm3.constants.H / ( qm3.constants.KB * temp )
    qv = 1.0
    nf = 0
    for f in freq:
        if( f >= fcut ):
            qv /= ( 1.0 - math.exp( - f * t ) )
        else:
            nf += 1
    qv = math.log( qv )
    # ZPE && Gibbs (kJ/mol)
    # G = F + PV = - RT Ln Q + nRT   <<  (nRT cancels out when calculating relative terms...)
    zz = sum( freq[nf:] ) * 0.5 * 100.0 * qm3.constants.C * qm3.constants.H * qm3.constants.NA * 1.0e-3
    gg = - qm3.constants.R * temp * ( qt + qr + qv ) * 1.0e-3
    return( ( zz, gg ) )



def sub_hessian( hess, atm_I, atm_J, size ):
    II  = 3 * atm_I
    JJ  = 3 * atm_J
    out = []
    for i in [ 0, 1, 2 ]:
        for j in [ 0, 1, 2 ]:
            out.append( hess[(II+i)*size+II+j] )
        for j in [ 0, 1, 2 ]:
            out.append( hess[(II+i)*size+JJ+j] )
    for i in [ 0, 1, 2 ]:
        for j in [ 0, 1, 2 ]:
            out.append( hess[(JJ+i)*size+II+j] )
        for j in [ 0, 1, 2 ]:
            out.append( hess[(JJ+i)*size+JJ+j] )
    return( out )



# ----------------------------------------------------------------------------------
# Hessian Updaters
#
def update_bfgs( dx, dg, hess ):
    size = len( dx )
    vec  = []
    tys  = 0.0
    tsvs = 0.0
    for i in range( size ):
        t0 = 0.0
        for j in range( size ):
            t0 += dx[j] * hess[i*size+j]
        vec.append( t0 )
        tsvs += dx[i] * t0
        tys += dg[i] * dx[i]
    for i in range( size ):
        for j in range( size ):
            hess[i*size+j] += dg[i] * dg[j] / tys - vec[i] * vec[j] / tsvs
    print( " * Hessian updated (BFGS)" )



def update_sr1( dx, dg, hess ):
    size = len( dx )
    vec = []
    tvs = 0.0
    for i in range( size ):
        t0 = 0.0
        for j in range( size ):
            t0 += dx[j] * hess[i*size+j]
        vec.append( t0 - dg[i] )
        tvs += dx[i] * vec[i]
    for i in range( size ):
        for j in range( size ):
            hess[i*size+j] += - vec[i] * vec[j] / tvs
    print( " * Hessian updated (SR1)" )



def update_psb( dx, dg, hess ):
    size = len( dx )
    vec = []
    tss = 0.0
    tvs = 0.0
    for i in range( size ):
        t0 = 0.0
        for j in range( size ):
            t0 += dx[j] * hess[i*size+j]
        vec.append( t0 - dg[i] )
        tss += dx[i] * dx[i]
        tvs += vec[i] * dx[i]
    t0 = tss * tss
    for i in range( size ):
        for j in range( size ):
            hess[i*size+j] += dx[i] * dx[j] * tvs / t0 - ( vec[i] * dx[j] + vec[j] * dx[i] ) / tss
    print( " * Hessian updated (PsB)" )



def update_bofill( dx, dg, hess ):
    size = len( dx )
    vec = []
    tvv = 0.0
    tss = 0.0
    tvs = 0.0
    for i in range( size ):
        t0 = 0.0
        for j in range( size ):
            t0 += dx[j] * hess[i*size+j]
        vec.append( t0 - dg[i] )
        tvv += vec[i] * vec[i]
        tss += dx[i] * dx[i]
        tvs += vec[i] * dx[i]
    ph = tvs * tvs / ( tvv * tss )
    t0 = tss * tss
    for i in range( size ):
        for j in range( size ):
            hess[i*size+j] += - ph * ( vec[i] * vec[j] / tvs ) + \
                    ( 1. - ph ) * ( dx[i] * dx[j] * tvs / t0 - ( vec[i] * dx[j] + vec[j] * dx[i] ) / tss )
    print( " * Hessian updated (combined Bofill, Phi = %.6lf)"%( ph ) )



def manage_hessian( coor, grad, hess, should_update = False, update_func = update_bofill, dump_name = "update.dump" ):
    size = len( coor )
    if( should_update and os.access( dump_name, os.R_OK ) and os.stat( dump_name )[stat.ST_SIZE] == size * ( size * 8 + 16 ) ):
        f = open( dump_name, "rb" )
        dx = list( struct.unpack( "%dd"%( size ), f.read( size * 8 ) )[:] )
        dg = list( struct.unpack( "%dd"%( size ), f.read( size * 8 ) )[:] )
        hh = list( struct.unpack( "%dd"%( size * size ), f.read( size * size * 8 ) )[:] )
        f.close()
        t = 0.0
        for i in range( size ):
            dx[i] = coor[i] - dx[i]
            dg[i] = grad[i] - dg[i]
            t += dx[i] * dx[i]
        if( math.sqrt( t ) > 1.0e-4 ):
            update_func( dx, dg, hh )
        else:
            print( " * Hessian read from previous calculation (update.dump)" )
        for i in range( size * size ):
            hess[i] = hh[i]
    f = open( dump_name, "wb" )
    for i in range( size ):
        f.write( struct.pack( "d", coor[i] ) )
    for i in range( size ):
        f.write( struct.pack( "d", grad[i] ) )
    for i in range( size * size ):
        f.write( struct.pack( "d", hess[i] ) )
    f.close()
    


# ----------------------------------------------------------------------------------
# Connectivity
#
def connectivity( molec ):
    # searches only by chain on current residue and with the next one ( m * n * (n-1) / 2 )
    if( len( molec.res_lim ) > 2 ):
        print( ">> connectivity: QUICK python search..." )
        bond = []
        x = len( molec.res_lim )
        for l in range( len( molec.seg_lim ) - 1 ):
            for i in range( molec.seg_lim[l], molec.seg_lim[l+1] ):
                for j in range( molec.res_lim[i], molec.res_lim[i+1] - 1 ):
                    j3 = j * 3
                    rj = qm3.elements.r_cov[molec.anum[j]] + 0.05
                    # current residue
                    for k in range( j+1, molec.res_lim[i+1] ):
                        if( molec.anum[j] == 1 and molec.anum[k] == 1 ):
                            continue
                        rk = qm3.elements.r_cov[molec.anum[k]] + 0.05
                        t  = ( rj + rk ) * ( rj + rk )
                        if( distanceSQ( molec.coor[j3:j3+3], molec.coor[k*3:k*3+3] ) <= t ):
                            bond.append( [ j, k ] )
                    # i+1 one
                    if( i+2 < x ):
                        for k in range( molec.res_lim[i+1], molec.res_lim[i+2] ):
                            if( molec.anum[j] == 1 and molec.anum[k] == 1 ):
                                continue
                            rk = qm3.elements.r_cov[molec.anum[k]] + 0.05
                            t  = ( rj + rk ) * ( rj + rk )
                            if( distanceSQ( molec.coor[j3:j3+3], molec.coor[k*3:k*3+3] ) <= t ):
                                bond.append( [ j, k ] )
    # all atoms: N * (N-1) / 2
    else:
        if( conn_so ):
            print( ">> connectivity: FULL binary threaded search..." )
            bond = qm3.utils._conn.connectivity( os.sysconf( 'SC_NPROCESSORS_ONLN' ), molec )
        else:
            print( ">> connectivity: FULL python sequential search..." )
            bond = []
            for i in range( molec.natm - 1 ):
                i3 = i * 3
                ri = qm3.elements.r_cov[molec.anum[i]] + 0.05
                for j in range( i + 1, molec.natm ):
                    if( molec.anum[i] == 1 and molec.anum[j] == 1 ):
                        continue
                    rj = qm3.elements.r_cov[molec.anum[j]] + 0.05
                    t  = ( ri + rj ) * ( ri + rj )
                    if( distanceSQ( molec.coor[i3:i3+3], molec.coor[j*3:j*3+3] ) <= t ):
                        bond.append( [ i, j ] )
    return( bond )

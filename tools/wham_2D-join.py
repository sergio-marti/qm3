#!/usr/bin/env python3

import numpy
import glob
import sys


data  = glob.glob( "dat.??.??" )
npt_x = 23 * 2
npt_y = 26 * 2
temp  = 300.0
skip  = 9000


##############################################################
kB    = 8.31447086363e-3  # in kJ/mol
inf   = 1.e30
skip += 1

class _index_expression_class:
    maxint = float( 'inf' )
    def __getitem__( self, item ):
        if( type( item ) != type( () ) ):
            return ( item, )
        else:
            return item

    def __len__( self ):
        return self.maxint

    def __getslice__( self, start, stop ):
        if( stop == self.maxint ):
            stop = None
        return self[start:stop:None]

index_expression = _index_expression_class()

def wham( data, bins, t, conv_limit = 0.001 ):
    beta = 1.0 / ( kB * t )
    dim = len( bins )
    nwin = len( data )
    axes = list( map( binAxis, bins ) )
    hist = []
    pot = []
    for window in data:
        rc0, fc, rc = window
        hist.append( histogram( rc, bins ) )
        pot.append( biasPotential( rc0, fc, axes ) )
    hist = numpy.array( hist )
    pot = numpy.array( pot )
    npoints = integrate( hist )
    index = index_expression[::] + dim*index_expression[numpy.newaxis]
    f = numpy.zeros( nwin, numpy.float )
    lastf = f
    num = 1.0 * numpy.add.reduce( hist, 0 )
    cnt = 0
    while( True ):
        expnt = beta * ( f[index] - pot )
        OK    = numpy.greater( expnt, -200.0 )
        den   = numpy.add.reduce( npoints[index] * numpy.exp ( numpy.where ( OK, expnt, -200.0 ) ) )
        good  = numpy.not_equal( den, 0.0 )
        rho   = num/den
        rho   = numpy.where( good, rho, 0.0 )
        expnt = - beta * pot
        OK    = numpy.greater( expnt, -200.0 )
        f     = - numpy.log( integrate( rho * numpy.exp( numpy.where( OK, expnt, -200.0 ) ) ) ) / beta
        conv  = numpy.maximum.reduce( numpy.fabs( f - lastf ) )
        cnt  += 1
        if( conv < conv_limit ): break
        lastf = f
    print( cnt, conv ) 
    mask = numpy.equal( rho, 0.0 )
    pmf = -numpy.log( rho + mask ) / beta + inf * mask
    pmf_min = numpy.minimum.reduce( numpy.ravel( pmf ) )
    pmf = pmf - pmf_min
    return axes, rho, pmf

def binAxis( bin_spec ):
    min, max, n = bin_spec
    w = ( max - min ) / ( n - 1 )
    return min + w*numpy.arange( 0.0, n )

def histogram( data, bins ):
    dim = data.shape[1]
    if( dim != len( bins ) ):
        raise ValueError( 'Inconsistent data' )
    min = numpy.array( list( map( lambda b: b[0], bins ) ) )
    max = numpy.array( list( map( lambda b: b[1], bins ) ) )
    n = numpy.array( list( map( lambda b: b[2], bins ) ) )
    w = ( max - min ) / ( n - 1 )
    indices = ( ( data - min + w / 2 ) / w ).astype( numpy.int32 )
    h = numpy.zeros( *( tuple( n ), numpy.int32 ) )
    for i in range( data.shape[0] ):
        ind = indices[i]
        if( numpy.logical_and.reduce( numpy.logical_and( numpy.greater_equal( ind, 0 ), numpy.less( ind, n ) ) ) ):
            ind = tuple( ind )
            h[ind] = h[ind] + 1
    return h

def biasPotential( rc0, fc, axes ):
    dim = rc0.shape[0]
    p = numpy.array( 0.0 )
    for i in range( dim ):
        p = numpy.add.outer( p, fc[i] * ( axes[i] - rc0[i] ) ** 2 )
    return p

def integrate( a ):
    while( len( a.shape ) > 1 ):
        a = numpy.add.reduce( a, -1 )
    return a


data.sort()

windows = []
for fn in data:
    print( fn )
    data = open( fn, "rt" ).readlines()
    t = data[0].split()
    rc0 = [ float( t[1] ), float( t[3] )  ]
    fc = [ 0.5 * float( t[0] ), 0.5 * float( t[2] ) ]
    rc = [ [], [] ]
    for i in range( skip, len( data ) ):
        t = data[i].split()
        for j in [0, 1]:
            rc[j].append( float( t[j] ) )
    windows.append( ( numpy.array( rc0 ), numpy.array( fc ), numpy.transpose( numpy.array( rc ) ) ) )


def combine( window_list ):
    new = []
    for w in window_list:
        done = 0
        for i in range( len( new ) ):
            if( numpy.logical_and.reduce( numpy.equal( new[i][0], w[0] ) ) and
                numpy.logical_and.reduce( numpy.equal( new[i][1], w[1] ) ) ):
                new[i] = ( new[i][0], new[i][1], numpy.concatenate( ( new[i][2], w[2] ) ) )
                done = 1
                break
        if( not done ):
            new.append( w )
    return new

min1 =  inf
min2 =  inf
max1 = -inf
max2 = -inf
for rc, fc, data in windows:
    for point in data:
        min1 = min( min1, point[0] )
        max1 = max( max1, point[0] )
        min2 = min( min2, point[1] )
        max2 = max( max2, point[1] )

print( windows[1] )
windows = combine( windows )
print( len( windows ) )

sys.stdout.flush()

c1_bins = ( min1, max1, npt_x )
c2_bins = ( min2, max2, npt_y )
axes, rho, pmf = wham( windows, ( c1_bins, c2_bins ), temp, 0.01 )

big = numpy.equal( pmf, inf )
pmf_min = numpy.minimum.reduce( numpy.ravel( pmf ) )
pmf_out = numpy.choose( big, ( pmf - pmf_min, 0 ) )
pmf_max = numpy.maximum.reduce( numpy.ravel( pmf ) )

f = open( "pmf.dat", "wt" )
for i in range( npt_x ):
    for j in range( npt_y ):
        f.write( "%20.10lf%20.10lf%20.10lf\n"%( axes[0][i], axes[1][j], pmf_out[i][j] ) )
    f.write( "\n" )
f.close()

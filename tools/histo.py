#!/usr/bin/env python3
import  sys
import  qm3.mol
import  qm3.fio.dcd
import  qm3.utils
import  math
import  time

m = qm3.mol.molecule( sys.argv[1] )
d = qm3.fio.dcd.dcd()
d.open_read( sys.argv[2] )

bnd = []
dst = []
MAX = []
MIN = []
for i in range( 3, len( sys.argv ), 2 ):
    try:
        _a = sys.argv[i].split( "/" )
        _b = sys.argv[i+1].split( "/" )
        bnd.append( ( 3 * m.indx[_a[0]][int(_a[1])][_a[2]], 3 * m.indx[_b[0]][int(_b[1])][_b[2]] ) )
    except:
        bnd.append( ( 3 * ( int( sys.argv[i] ) - 1 ), 3 * ( int( sys.argv[i+1] ) - 1 ) ) )
    dst.append( [] )
    MAX.append( 0 )
    MIN.append( 1.e99 )

t0 = time.time()
n = len( bnd )
N = 0.0
while( d.next( m ) ):
    for i in range( n ):
        t = qm3.utils.distance( m.coor[bnd[i][0]:bnd[i][0]+3], m.coor[bnd[i][1]:bnd[i][1]+3] )
        dst[i].append( t )
        MAX[i] = max( MAX[i], t )
        MIN[i] = min( MIN[i], t )
    N += 1.0
d.close()
print( time.time() - t0 )

RNG = []
for i in range( n ):
    _m = math.floor( MIN[i] * 10 ) // 5 * 0.5
    _M = math.ceil( MAX[i] * 10 ) / 10
    RNG.append( [] )
    j = 0
    while( _m + j < _M ):
        RNG[-1].append( _m + j )
        j += 0.5
    RNG[-1].append( _m + j )

f = open( "histo.log", "wt" )
for i in range( n ):
    lin = "\n%6d %-6s %6d %-6s: %12.6lf %12.6lf %12.6lf"%(
        m.resi[bnd[i][0]//3], m.labl[bnd[i][0]//3],
        m.resi[bnd[i][1]//3], m.labl[bnd[i][1]//3],
        MIN[i], sum( dst[i] ) / N, MAX[i] )
    print( lin )
    f.write( lin + "\n" )
    tmp = [ 0.0 for j in range( len( RNG[i] ) ) ]
    med = tmp[:]
    siz = RNG[i][-1] - RNG[i][0]
    for j in range( int( N ) ):
        who = int( ( ( dst[i][j] - RNG[i][0] ) * 10 ) // 5 )
        tmp[who] += 1.0
        med[who] += dst[i][j]
    for j in range( len( RNG[i] ) - 1 ):
        if( tmp[j] > 0 ):
            if( tmp[j] / N > 0.1 ):
                lin = "\t%2.1lf - %2.1lf: %6.4lf   %4.1lf%%"%( RNG[i][j], RNG[i][j+1], med[j] / tmp[j], tmp[j] / N * 100 )
                print( lin )
                f.write( lin + "\n" )
f.close()

f = open( "histo.dat", "wt" )
for i in range( int( N ) ):
    for j in range( n ):
        f.write( "%14.6lf"%( dst[j][i] ) )
    f.write( "\n" )
f.close()

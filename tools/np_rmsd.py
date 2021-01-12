#!/usr/bin/env python3
import  sys
import  numpy
import  qm3.mol
import  qm3.elements
import  qm3.fio.dcd

mol = qm3.mol.molecule( sys.argv[1] )

sel = [ i * 3 for i in range( mol.natm ) if mol.labl[i] in [ "C", "N", "CA" ] ]

ref = numpy.array( sum( [ mol.coor[i:i+3] for i in sel ], [] ) ).reshape( len( sel ), 3 )
ref -= numpy.average( ref, axis=0 )

dcd = qm3.fio.dcd.dcd()
dcd.open_read( sys.argv[2] )

while( dcd.next( mol ) ):
    cur = numpy.array( sum( [ mol.coor[i:i+3] for i in sel ], [] ) ).reshape( len( sel ), 3 )
    cur -= numpy.average( cur, axis=0 )

    cov = numpy.dot( cur.T, ref )
    r1, s, r2 = numpy.linalg.svd( cov )
    if( numpy.linalg.det( cov ) < 0 ):
        r2[2,:] *= -1.0
    mod = numpy.dot( cur, numpy.dot( r1, r2 ) )
    print( numpy.sqrt( numpy.average( numpy.sum( ( ref - mod ) ** 2, axis=1 ) ) ) )

dcd.close()

# -*- coding: iso-8859-1 -*-
from __future__ import print_function, division
import    sys
if( sys.version_info[0] == 2 ):
    range = xrange
import    random
import    time
import    math
import    os
import    qm3.utils.msi


# ---------------------------------------
# python samples/test.msi &
# sleep 1
# for i in {0..9}; do
#     python samples/test.msi $i &
# done
# ---------------------------------------


random.seed()

npc = 10
pth = os.path.join( os.getcwd(), "qm3_msi" )

if( len( sys.argv ) == 1 ):
    try:
        os.unlink( pth )
    except:
        pass
    qm3.utils.msi.server( npc, unix = pth )
else:
    who = int( sys.argv[1] )
    con = qm3.utils.msi.client( node = who, unix = pth )
    time.sleep( 10 * random.random() )
    con.barrier()
    time.sleep( 10 * random.random() )
    con.barrier()
    time.sleep( 10 * random.random() )
    con.barrier()
    nel = 2300
    if( con.node == 0 ):
        for i in range( 1, npc ):
            tmp = con.recv( i, nel )
            sys.stderr.write( "0 << %03d: read: %.3lf\n"%( i, sum( tmp ) ) )
            sys.stderr.flush()
            con.send( i, tmp )
    else:
        tmp = [ random.random() for i in range( nel ) ]
        con.send( 0, tmp )
        sys.stderr.write( "%03d >> 0: writ: %.3lf\n"%( con.node, sum( tmp ) ) )
        sys.stderr.flush()
        pmt = con.recv( 0, nel )
        sys.stderr.write( "%03d     : diff: %.3lf\n"%( con.node, sum( [ math.fabs( i - j ) for i,j in zip( tmp, pmt ) ] ) ) )
        sys.stderr.flush()
    time.sleep( 10 * random.random() )
    con.barrier()
    sys.stderr.write( "%03d: done!\n"%( con.node ) )
    sys.stderr.flush()

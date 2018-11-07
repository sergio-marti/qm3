import random
import qm3.maths.matrix
import qm3.utils._mpi
import sys

random.seed()
n = 10
m = [ random.random() for i in range( n*n ) ]

# 2 racks: mpirun -n 5 python script
grp = { 0: [ 1, None ], 1: [ None, 0 ], 2: [ 3, None ], 3: [ 4, 2 ], 4: [ None, 3 ] }

node, nwin = qm3.utils._mpi.init()

# sleep
if( grp[node][1] != None ):
    t = qm3.utils._mpi.recv_r8( grp[node][1], 1 )
    sys.stderr.write( "node: " + str( node ) + ", prev result: " + str( t[0] ) + "\n" )

# calculate
t = qm3.maths.matrix.det( m, n )
sys.stderr.write( "node: " + str( node ) + ", curr result: " + str( t ) + "\n" )

# relieve
if( grp[node][0] != None ):
    qm3.utils._mpi.send_r8( grp[node][0], [ t ] )

# sync
qm3.utils._mpi.barrier()

qm3.utils._mpi.stop()

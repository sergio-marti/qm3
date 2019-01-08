import qm3.utils._shm
import struct

n = 10
p = qm3.utils._shm.alloc( n * 8 )

qm3.utils._shm.write( p, b"".join( [ struct.pack( "d", i ) for i in range( n ) ] ) )
print( qm3.utils._shm.read_r8( p, n ) )

qm3.utils._shm.write_r8( p, [ float( i ) for i in range( n ) ] )
print( list( struct.unpack( "%dd"%( n ), qm3.utils._shm.read( p, n * 8 ) ) ) )

qm3.utils._shm.clean( p )

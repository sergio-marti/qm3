#!/usr/bin/env python3
import  sys
import  struct
import  pickle
import  numpy

f = open( sys.argv[1], "rb" )
inp = pickle.load( f ).flatten()
out = pickle.load( f ).flatten()
out -= numpy.min( out )
f.close()
siz = len( out )
dim = len( inp ) // siz
print( siz, dim )
f = open( "raw_data", "wb" )
f.write( struct.pack( "i", dim ) )
f.write( struct.pack( "i", siz ) )
for e in inp:
    f.write( struct.pack( "d", e ) )
for e in out:
    f.write( struct.pack( "d", e ) )
f.close()

#!/usr/bin/env python3
import  sys
import  struct
import  os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3" 
import  tensorflow as tf
import  pickle
import  numpy

f = open( "raw_model", "wb" )
x = tf.keras.models.load_model( sys.argv[1] ).trainable_weights
print( len( x ) )
f.write( struct.pack( "i", len( x ) ) )
for m in x:
    t = m.numpy()
    if( len( t.shape ) == 2 ):
        print( t.shape )
        f.write( struct.pack( "i", t.shape[0] ) )
        f.write( struct.pack( "i", t.shape[1] ) )
    else:
        print( (1,t.shape[0]) )
        f.write( struct.pack( "i", 1 ) )
        f.write( struct.pack( "i", t.shape[0] ) )
    for e in t.flatten():
        f.write( struct.pack( "d", e ) )
f.close()

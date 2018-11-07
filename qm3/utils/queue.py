# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	socket
import	codecs
try:
	import	cPickle as pickle
except:
	import	pickle


class Queue( object ):

	def __init__( self, cpus ):
		self.unix = "/tmp/xqueue.%d.%d"%( os.getuid(), os.getpid() )
		self.sckt = socket.socket( socket.AF_UNIX, socket.SOCK_STREAM )
		self.cpus = cpus
		self.slen = 1024 * 1024
		self.sckt.bind( self.unix )
		self.sckt.listen( self.cpus + 1 )
		self.data = []


	def serve( self ):
		for i in range( self.cpus ):
			sck, acc = self.sckt.accept()
			buf = sck.recv( self.slen ).decode( "ascii" )
			while( len( buf ) < 20 ):
				buf += sck.recv( self.slen ).decode( "ascii" )
			siz = int( buf[0:20] )
			buf = buf[20:]
			while( len( buf ) < siz ):
				buf += sck.recv( self.slen ).decode( "ascii" )
			self.data += pickle.loads( codecs.decode( buf.encode( "ascii" ), "base64" ) )
		self.sckt.close()
		os.unlink( self.unix )


	def client( self, data ):
		sckt = socket.socket( socket.AF_UNIX, socket.SOCK_STREAM )
		sckt.connect( self.unix )
		buf = codecs.encode( pickle.dumps( data ), "base64" ).decode( "ascii" )
		buf = "%020d"%( len( buf ) ) + buf
		sckt.send( buf.encode( "ascii" ) )
		sckt.close()




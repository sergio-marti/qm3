# -*- coding: iso-8859-1 -*-

# Message Socket interface

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	socket
import	struct


class server:
	def __init__( self, ncpu, pnum = 6969, srvr = None ):
		self.ncpu = ncpu
		self.node = 0
		self.slen = 1024
		if( srvr ):
			hstn = srvr
		else:
			hstn = socket.gethostbyname( socket.gethostname() )
		sys.stderr.write( "[server] listening at: %s\n"%( hstn ) )
		sys.stderr.flush()
		sckt = socket.socket( socket.AF_INET, socket.SOCK_STREAM )
		sckt.bind( ( hstn, pnum ) )
		self.chld = {}
		while( len( self.chld ) < ncpu - 1 ):
			sckt.listen( ncpu )
			chld, addr = sckt.accept()
			node = len( self.chld ) + 1
			self.chld[node] = chld
			self.__send( chld, struct.pack( "i", node ) )
			sys.stderr.write( "[server] added: %d/%d\n"%( node, ncpu - 1 ) )
			sys.stderr.flush()
		sckt.close()
		sys.stderr.write( "[server] done!\n" )
		sys.stderr.flush()


	def stop( self ):
		for i in sorted( self.chld ):
			self.chld[i].close()


	def __send( self, sck, msg ):
		s = len( msg )
		j = 0
		for i in range( s // self.slen ):
			sck.send( msg[j:j+self.slen] )
			j += self.slen
		s %= self.slen
		if( s != 0 ):
			sck.send( msg[j:] + b"0" * ( self.slen - s ) )


	def __recv( self, sck, siz ):
		msg = sck.recv( self.slen )
		while( len( msg ) < siz ):
			msg += sck.recv( self.slen )
		return( msg[0:siz] )


	def barrier( self ):
		out = []
		for i in sorted( self.chld ):
			out.append( struct.unpack( "i", self.__recv( self.chld[i], 4 ) )[0] )
		for i in sorted( self.chld ):
			self.__send( self.chld[i], struct.pack( "i", 1 ) )
		return( sum( out ) == len( self.chld ) )


	def send_r8( self, who, lst ):
		msg = b""
		for r in lst:
			msg += struct.pack( "d", r )
		self.__send( self.chld[who], msg )


	def recv_r8( self, who, siz ):
		return( list( struct.unpack( "%dd"%( siz ), self.__recv( self.chld[who], siz * 8 ) ) ) )


	def send_i4( self, who, lst ):
		msg = b""
		for r in lst:
			msg += struct.pack( "i", r )
		self.__send( self.chld[who], msg )


	def recv_i4( self, who, siz ):
		return( list( struct.unpack( "%di"%( siz ), self.__recv( self.chld[who], siz * 4 ) ) ) )





class client:
	def __init__( self, hstn, pnum = 6969 ):
		self.slen = 1024
		self.sckt = socket.socket( socket.AF_INET, socket.SOCK_STREAM )
		self.sckt.connect( ( hstn, pnum ) )
		self.node = struct.unpack( "i", self.__recv( 4 ) )[0]
		sys.stderr.write( "[client:%d] connected!\n"%( self.node ) )
		sys.stderr.flush()


	def stop( self ):
		self.sckt.close()


	def __send( self, msg ):
		s = len( msg )
		j = 0
		for i in range( s // self.slen ):
			self.sckt.send( msg[j:j+self.slen] )
			j += self.slen
		s %= self.slen
		if( s != 0 ):
			self.sckt.send( msg[j:] + b"0" * ( self.slen - s ) )


	def __recv( self, siz ):
		msg = self.sckt.recv( self.slen )
		while( len( msg ) < siz ):
			msg += self.sckt.recv( self.slen )
		return( msg[0:siz] )


	def barrier( self ):
		self.__send( struct.pack( "i", 1 ) )
		return( struct.unpack( "i", self.__recv( 4 ) )[0] == 1 )


	def send_r8( self, lst ):
		msg = b""
		for r in lst:
			msg += struct.pack( "d", r )
		self.__send( msg )


	def recv_r8( self, siz ):
		return( list( struct.unpack( "%dd"%( siz ), self.__recv( siz * 8 ) ) ) )


	def send_i4( self, lst ):
		msg = b""
		for r in lst:
			msg += struct.pack( "i", r )
		self.__send( msg )


	def recv_i4( self, siz ):
		return( list( struct.unpack( "%di"%( siz ), self.__recv( siz * 4 ) ) ) )

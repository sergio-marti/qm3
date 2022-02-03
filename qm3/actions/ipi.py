# -*- coding: iso-8859-1 -*-
import socket
import struct
import qm3.constants
import qm3.problem



class ipi_problem( qm3.problem.template ):
    def __init__( self ):
        qm3.problem.template.__init__( self )
        self.sckt = None
        self.bead = None
        self.stat = b"NEEDINIT    "
#        self.viri = b"".join( [ struct.pack( "d", 0.0 ) for i in range( 9 ) ] )
        self.viri = b"\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00"
        self._ce = 1. / qm3.constants.H2J
        self._cg = self._ce * qm3.constants.A0


    def __recv( self, siz ):
        out = self.sckt.recv( siz )
        rem = siz - len( out )
        while( rem > 0 ):
            out += self.sckt.recv( rem )
            rem = siz - len( out )
        return( out )


    def map_coor( self, coor ):
        self.coor = coor[:]


    def map_grad( self ):
        return( [ - i * self._cg for i in self.grad ] )


    def connect( self, info ):
        # path of the unix socket
        if( type( info ) == str ):
            self.sckt = socket.socket( socket.AF_UNIX, socket.SOCK_STREAM )
            self.sckt.connect( info )
        # ( hostname, port_number = 31415 )
        elif( ( type( info ) == tuple or type( info ) == list ) and len( info ) == 2 ):
            self.sckt = socket.socket( socket.AF_INET, socket.SOCK_STREAM )
            self.sckt.connect( ( socket.gethostbyaddr( info[0] )[0], int( info[1] ) ) )
        try:
            while( True ):
                msg = self.__recv( 12 )
                if( msg == b"EXIT        " ):
                    return
                elif( msg == b"STATUS      " ):
                    self.sckt.sendall( self.stat )
                elif( msg == b"POSDATA     " ):
                    rbox = self.__recv( 72 )
                    ibox = self.__recv( 72 )
                    size = struct.unpack( "i", self.__recv( 4 ) )[0] * 3
                    coor = [ i * qm3.constants.A0 for i in struct.unpack( "%dd"%( size ), self.__recv( size * 8 ) ) ]
                    self.map_coor( coor )
                    self.get_grad()
                    self.stat = b"HAVEDATA    "
                elif( msg == b"GETFORCE    " ):
                    self.sckt.sendall( b"FORCEREADY  " )
                    self.sckt.sendall( struct.pack( "d", self.func * self._ce ) )
                    grad = self.map_grad()
                    size = len( grad )
                    self.sckt.sendall( struct.pack( "i", size // 3 ) )
                    self.sckt.sendall( b"".join( [ struct.pack( "d", grad[i] ) for i in range( size ) ] ) )
                    self.sckt.sendall( self.viri )
                    self.sckt.sendall( struct.pack( "i", 1 ) )
                    self.sckt.sendall( b"\x00" )
                    self.stat = b"READY       "
                elif( msg == b"INIT        " ):
                    self.bead = struct.unpack( "i", self.__recv( 4 ) )[0]
                    print( ">> assgined bead: ", self.bead )
                    self.__recv( struct.unpack( "i", self.__recv( 4 ) )[0] )
                    self.stat = b"READY       "
        finally:
            self.sckt.close()

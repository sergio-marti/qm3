# -*- coding: iso-8859-1 -*-

# Message Socket interface

import socket
import struct
import time
import threading


class server:
    def __recv( self, sckt ):
        msg = sckt.recv( self.slen )
        try:
            cmd = msg[0:1]
            who = int( msg[1:4] )
            dst = int( msg[4:7] )
            siz = int( msg[7:17] )
            dat = msg[17:]
            while( len( dat ) < siz ):
                dat += sckt.recv( self.slen )
            return( cmd, who, dst, dat[0:siz] )
        except:
            return( b"#", -1, -1, b"" )


    def __send( self, sckt, msg ):
        s = len( msg )
        j = 0
        for i in range( s // self.slen ):
            sckt.sendall( msg[j:j+self.slen] )
            j += self.slen
        s %= self.slen
        if( s != 0 ):
            sckt.sendall( msg[j:] + b"0" * ( self.slen - s ) )


    def __serve( self, chld ):
        cmd, who, dst, dat = self.__recv( chld )
        if( cmd == b"#" or who < 0 or who >= self.ncpu ):
            chld.close()
            sys.stderr.write( "[server] rejected invalid node(%d) or malformed message...\n"%( who ) )
            sys.stderr.flush()
            return()
        self.lock.acquire()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if( cmd == b"B" ):
            self.barc[who] += 1
            if( not( self.barc[who] in self.barb and self.barc[who] in self.bard ) ):
                self.barb[self.barc[who]] = [ 0 for i in range( self.ncpu ) ]
                self.bard[self.barc[who]] = True
                sys.stderr.write( "[server] barrier %d +++ (node: %03d)\n"%( self.barc[who], who ) )
                sys.stderr.flush()
            self.barb[self.barc[who]][who] = 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        elif( cmd == b"P" ):
            if( self.bard[self.barc[who]] ):
                self.bard[self.barc[who]] = sum( self.barb[self.barc[who]] ) < self.ncpu
                if( self.bard[self.barc[who]] ):
                    self.__send( chld, b"10" )
                else:
                    self.barb[self.barc[who]][who] = 0
                    self.__send( chld, b"00" )
            else:
                self.barb[self.barc[who]][who] = 0
                if( sum( self.barb[self.barc[who]] ) == 0 ):
                    sys.stderr.write( "[server] barrier %d --- (node: %03d)\n"%( self.barc[who], who ) )
                    sys.stderr.flush()
                    del self.barb[self.barc[who]]
                    del self.bard[self.barc[who]]
                self.__send( chld, b"00" )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        elif( cmd == b"W" ):
            if( self.data[dst][who] == b"" ):
                self.data[dst][who] = dat
                self.__send( chld, b"10" )
            else:
                self.__send( chld, b"00" )
        elif( cmd == b"R" ):
            if( self.data[who][dst] != b"" ):
                self.__send( chld, b"1" + self.data[who][dst] )
                self.data[who][dst] = b""
            else:
                self.__send( chld, b"00" )
        self.lock.release()
        chld.close()


    def __init__( self, ncpu, inet = ( "", 6969 ), unix = None ):
        self.ncpu = ncpu
        self.slen = 1024
        self.data = [ [ b"" for j in range( self.ncpu ) ] for i in range( self.ncpu ) ]
        self.barb = {}
        self.bard = {}
        self.barc = [ -1 for i in range( self.ncpu ) ]
        self.lock = threading.Lock()
        if( unix ):
            sys.stderr.write( "[server] listening at: %s\n"%( unix ) )
            sys.stderr.flush()
            self.sckt = socket.socket( socket.AF_UNIX, socket.SOCK_STREAM )
            self.sckt.bind( unix )
        else:
            if( inet[0] == "" ):
                inet = ( socket.gethostbyname( socket.gethostname() ), inet[1] )
            sys.stderr.write( "[server] listening at: %s\n"%( str( inet ) ) )
            sys.stderr.flush()
            self.sckt = socket.socket( socket.AF_INET, socket.SOCK_STREAM )
            self.sckt.bind( inet )
        while( True ):
            self.sckt.listen( ncpu * 10 )
            chld, addr = self.sckt.accept()
            threading.Thread( target = self.__serve, args = ( chld, ) ).start()
#        self.sckt.close()
#        sys.stderr.write( "[server] done!\n" )
#        sys.stderr.flush()



class client:
    def __send( self, sckt, msg ):
        s = len( msg )
        j = 0
        for i in range( s // self.slen ):
            sckt.sendall( msg[j:j+self.slen] )
            j += self.slen
        s %= self.slen
        if( s != 0 ):
            sckt.sendall( msg[j:] + b"0" * ( self.slen - s ) )


    def __recv( self, sckt, siz ):
        msg = sckt.recv( self.slen )
        if( msg[0:1] == b"1" ):
            while( len( msg ) <= siz ):
                msg += sckt.recv( self.slen )
            return( True, msg[1:siz+1] )
        else:
            return( False, b"" )


    def __init__( self, node, inet = ( "", 6969 ), unix = None ):
        self.node = node
        self.slen = 1024
        if( unix ):
            self.kind = socket.AF_UNIX
            self.addr = unix
        else:
            if( inet[0] == "" ):
                inet = ( socket.gethostbyname( socket.gethostname() ), inet[1] )
            self.kind = socket.AF_INET
            self.addr = inet


    def barrier( self ):
        sckt = socket.socket( self.kind, socket.SOCK_STREAM )
        sckt.connect( self.addr )
        self.__send( sckt, b"B%03d00000"%( self.node ) )
        sckt.close()
        flg = True
        while( flg ):
            time.sleep( 0.1 )
            sckt = socket.socket( self.kind, socket.SOCK_STREAM )
            sckt.connect( self.addr )
            self.__send( sckt, b"P%03d00000"%( self.node ) )
            flg, dat = self.__recv( sckt, 1 )
            sckt.close()


    def send( self, dst, lst ):
        msg = b"W%03d%03d%010d"%( self.node, dst, len( lst ) * 8 )
        for r in lst:
            msg += struct.pack( "d", float( r ) )
        sckt = socket.socket( self.kind, socket.SOCK_STREAM )
        sckt.connect( self.addr )
        self.__send( sckt, msg )
        flg, tmp = self.__recv( sckt, 1 )
        sckt.close()
        while( not flg ):
            time.sleep( 0.1 )
            sckt = socket.socket( self.kind, socket.SOCK_STREAM )
            sckt.connect( self.addr )
            self.__send( sckt, msg )
            flg, tmp = self.__recv( sckt, 1 )
            sckt.close()


    def recv( self, src, siz ):
        sckt = socket.socket( self.kind, socket.SOCK_STREAM )
        sckt.connect( self.addr )
        self.__send( sckt, b"R%03d%03d00"%( self.node, src ) )
        flg, msg = self.__recv( sckt, siz * 8 )
        sckt.close()
        while( not flg ):
            time.sleep( 0.1 )
            sckt = socket.socket( self.kind, socket.SOCK_STREAM )
            sckt.connect( self.addr )
            self.__send( sckt, b"R%03d%03d00"%( self.node, src ) )
            flg, msg = self.__recv( sckt, siz * 8 )
            sckt.close()
        return( list( struct.unpack( "%dd"%( siz ), msg ) ) )


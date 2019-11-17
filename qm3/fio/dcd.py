# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange

import struct



class dcd( object ):

    def __init__( self, fname = None, qprint = True ):
        self.N  = 0
        self.X  = []
        self.Y  = []
        self.Z  = []
        self._S = []
        self._F = None
        self._E = None
        self._Q = None
        self._n = 0
        self._C = 0
        self._W = None
        if( fname ):
            self.read( fname, qprint )


    def read( self, fname, qprint = True ):
        self.__init__()
        self._F = open( fname, "rb" )    
        if( qprint ):
            print( "* [%s] "%( fname ) + (55-len(fname))*"-" )
        if( struct.unpack( "<i", self._F.read( 4 ) )[0] == 84 ):
            self._E = "<"
            if( qprint ):
                print( "+ \tLittle Endian" )
        else:
            self._E = ">"
            if( qprint ):
                print( "+ \tBig Endian" )
        self._F.read( 4 )
        t = struct.unpack( self._E + "i", self._F.read( 4 ) )[0]
        if( qprint ):
            print( "+ \tNFrames: %d"%( t ) )
        self._F.read( 28 )
        f = struct.unpack( self._E + "i", self._F.read( 4 ) )[0]
        self._F.read( 4 )
        self._Q = struct.unpack( self._E + "i", self._F.read( 4 ) )[0]
        self._F.read( 40 )
        self._F.read( struct.unpack( self._E + "i", self._F.read( 4 ) )[0] + 8 )
        self.N = struct.unpack( self._E + "i", self._F.read( 4 ) )[0]
        if( qprint ):
            print( "+ \tAtoms: %d"%( self.N ) )
        self._F.read( 4 )
        self._n = self.N - f
        self._S = []
        if( f > 0 ):
            if( qprint ):
                print( "+ \tFixed Atoms: %d"%( f ) )
                print( "+ \tFree  Atoms: %d"%( self._n ) )
            self._F.read( 4 )
            for i in struct.unpack( self._E + "%di"%( self._n ), self._F.read( 4 * self._n ) ):
                self._S.append( i - 1 )
            self._F.read( 4 )
        self._C = 0
        for i in range( self.N ):
            self.X.append( 9999.0 )
            self.Y.append( 9999.0 )
            self.Z.append( 9999.0 )
        if( qprint ):
            print( 60*"-" )
        return( t )


    def next( self ):
        try:
            self._F.read( 4 )
        except AttributeError:
            print( "* No DCD open or End-Of-File reached...\n" )
            return( False )
        else:
            if( self._Q == 1 ):
                if( len( self._F.read( 56 ) ) != 56 ):
                    return( False )
            if( self.N != self._n and self._C > 0 ):
                ii = 4 * self._n
                jj = 20 + 3 * ii
                bb = self._F.read( jj )
                if( len( bb ) != jj ):
                    return( False )
                kk = 0
                cx = struct.unpack( self._E + "%df"%( self._n ), bb[kk:kk+ii] )
                kk += ( 8 + ii )
                cy = struct.unpack( self._E + "%df"%( self._n ), bb[kk:kk+ii] )
                kk += ( 8 + ii )
                cz = struct.unpack( self._E + "%df"%( self._n ), bb[kk:kk+ii] )
                for i in range( self._n ):
                    self.X[self._S[i]] = cx[i]
                    self.Y[self._S[i]] = cy[i]
                    self.Z[self._S[i]] = cz[i]
            else:
                ii = 4 * self.N
                jj = 20 + 3 * ii
                bb = self._F.read( jj )
                if( len( bb ) != jj ):
                    return( False )
                kk = 0
                cx = struct.unpack( self._E + "%df"%( self.N ), bb[kk:kk+ii] )
                kk += ( 8 + ii )
                cy = struct.unpack( self._E + "%df"%( self.N ), bb[kk:kk+ii] )
                kk += ( 8 + ii )
                cz = struct.unpack( self._E + "%df"%( self.N ), bb[kk:kk+ii] )
                for i in range( self.N ):
                    self.X[i] = cx[i]
                    self.Y[i] = cy[i]
                    self.Z[i] = cz[i]
            del cx, cy, cz
            self._C += 1
            return( True )


    def close( self ):
        self._F.close()
        # fix number of frames added...
        if( self._W != None ):
            self._F = open( self._W, "rb+" )
            self._F.seek( 8 )
            self._F.write( struct.pack( "i", self._C ) )
            self._F.seek( 20 )
            self._F.write( struct.pack( "i", self._C ) )
            self._F.close()


    def write( self, fname, natoms, sele = None ):
        self.__init__()
        self._W = fname
        self._C = 0
        self._F = open( fname, "wb" )
        self._F.write( struct.pack( "i", 84 ) + b"CORD" )
        self._F.write( struct.pack( "i", -1 ) )
        self._F.write( struct.pack( "i", 1 ) )
        self._F.write( struct.pack( "i", 1 ) )
        self._F.write( struct.pack( "i", -1 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self.N = natoms
        if( sele != None and len( sele ) < self.N ):
            self._S = sele[:]
            self._n = len( sele )
            f = self.N - self._n
        else:
            self._S = []
            self._n = self.N
            f = 0
        self._F.write( struct.pack( "i", 3 * self._n ) )
        self._F.write( struct.pack( "i", f ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 0 ) )
        self._F.write( struct.pack( "i", 84 ) )
        self._F.write( struct.pack( "i", 4 + 80 ) )
        self._F.write( struct.pack( "i", 1 ) )
        self._F.write( b"* Created with QMCube (qm3)                                                     " )
        self._F.write( struct.pack( "i", 4 + 80 ) )
        self._F.write( struct.pack( "i", 4 ) + struct.pack( "i", self.N ) + struct.pack( "i", 4 ) )
        if( self._S != [] ):
            b = struct.pack( "i", 4 * len( sele ) )
            self._F.write( b )
            for i in self._S:
                self._F.write( struct.pack( "i", i+1 ) )
            self._F.write( b )
        for i in range( self.N ):
            self.X.append( 9999. )
            self.Y.append( 9999. )
            self.Z.append( 9999. )


    def append( self ):
        try:
            self._F.tell()
        except AttributeError:
            print( "* No DCD open..." )
            return( False )
        else:
            if( self.N != self._n and self._C > 0 ):
                t = struct.pack( "i", self._n * 4 )
                bx = b""
                by = b""
                bz = b""
                for i in self._S:
                    bx += struct.pack( "f", self.X[i] )
                    by += struct.pack( "f", self.Y[i] )
                    bz += struct.pack( "f", self.Z[i] )
            else:
                t = struct.pack( "i", self.N * 4 )
                bx = b""
                by = b""
                bz = b""
                for i in range( self.N ):
                    bx += struct.pack( "f", self.X[i] )
                    by += struct.pack( "f", self.Y[i] )
                    bz += struct.pack( "f", self.Z[i] )
            self._F.write( t + bx + t + t + by + t + t + bz + t )
            self._F.flush()
            self._C += 1
            return( True )


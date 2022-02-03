# -*- coding: iso-8859-1 -*-
import struct
import os
import stat

try:

    from qm3.fio._dcd import dcd

except:

    class dcd( object ):
    
        def __init__( self ):
            self.__clean()
                
    
        def __clean( self ):
            self.natm  = 0
            self.sele = []
            self.fdes = None
            self.fsiz = 0
            self.endi = ""
            self.crys = 0
            self.free = 0
            self.curr = 0
            self.writ = ""
            self.head = 0
    
    
        def open_read( self, fname, qprint = True ):
            self.__clean()
            self.fdes = open( fname, "rb" )    
            self.fsiz = os.stat( fname )[stat.ST_SIZE]
            if( qprint ):
                print( "* [%s] "%( fname ) + (55-len(fname))*"-" )
            if( struct.unpack( "<i", self.fdes.read( 4 ) )[0] == 84 ):
                self.endi = "<"
                if( qprint ):
                    print( "+ \tLittle Endian" )
            else:
                self.endi = ">"
                if( qprint ):
                    print( "+ \tBig Endian" )
            self.fdes.read( 4 )
            t = struct.unpack( self.endi + "i", self.fdes.read( 4 ) )[0]
            if( qprint ):
                print( "+ \tNFrames: %d"%( t ) )
            self.fdes.read( 28 )
            f = struct.unpack( self.endi + "i", self.fdes.read( 4 ) )[0]
            self.fdes.read( 4 )
            self.crys = struct.unpack( self.endi + "i", self.fdes.read( 4 ) )[0]
            self.fdes.read( 40 )
            self.fdes.read( struct.unpack( self.endi + "i", self.fdes.read( 4 ) )[0] + 8 )
            self.natm = struct.unpack( self.endi + "i", self.fdes.read( 4 ) )[0]
            if( qprint ):
                print( "+ \tAtoms: %d"%( self.natm ) )
            self.fdes.read( 4 )
            self.free = self.natm - f
            self.sele = []
            if( f > 0 ):
                if( qprint ):
                    print( "+ \tFixed Atoms: %d"%( f ) )
                    print( "+ \tFree  Atoms: %d"%( self.free ) )
                self.fdes.read( 4 )
                for i in struct.unpack( self.endi + "%di"%( self.free ), self.fdes.read( 4 * self.free ) ):
                    self.sele.append( i - 1 )
                self.fdes.read( 4 )
            self.curr = 0
            self.head = self.fdes.tell()
            if( qprint ):
                print( 60*"-" )
            return( t )
    
    
        def next( self, molec ):
            try:
                self.fdes.read( 4 )
            except AttributeError:
                print( "* No DCD open or End-Of-File reached...\n" )
                return( False )
            else:
                if( self.crys == 1 ):
                    if( len( self.fdes.read( 56 ) ) != 56 ):
                        return( False )
                if( self.natm != self.free and self.curr > 0 ):
                    ii = 4 * self.free
                    if( self.fsiz - self.fdes.tell() - 20 - 3 * ii < 0 ):
                        return( False )
                    cx = struct.unpack( self.endi + "%df"%( self.natm ), self.fdes.read( ii ) )
                    self.fdes.read( 8 )
                    cy = struct.unpack( self.endi + "%df"%( self.natm ), self.fdes.read( ii ) )
                    self.fdes.read( 8 )
                    cz = struct.unpack( self.endi + "%df"%( self.natm ), self.fdes.read( ii ) )
                    self.fdes.read( 4 )
                    for i in range( self.free ):
                        i3 = self.sele[i] * 3
                        molec.coor[i3  ] = cx[i]
                        molec.coor[i3+1] = cy[i]
                        molec.coor[i3+2] = cz[i]
                    del cx, cy, cz
                else:
                    ii = 4 * self.natm
                    if( self.fsiz - self.fdes.tell() - 20 - 3 * ii < 0 ):
                        return( False )
                    cx = struct.unpack( self.endi + "%df"%( self.natm ), self.fdes.read( ii ) )
                    self.fdes.read( 8 )
                    cy = struct.unpack( self.endi + "%df"%( self.natm ), self.fdes.read( ii ) )
                    self.fdes.read( 8 )
                    cz = struct.unpack( self.endi + "%df"%( self.natm ), self.fdes.read( ii ) )
                    self.fdes.read( 4 )
                    for i in range( self.natm ):
                        i3 = i * 3
                        molec.coor[i3  ] = cx[i]
                        molec.coor[i3+1] = cy[i]
                        molec.coor[i3+2] = cz[i]
                    del cx, cy, cz
                self.curr += 1
                return( True )
    
    
        def goto( self, num ):
            if( num >= 0 ):
                if( self.natm != self.free and num > 1 ):
                    dsp = 4 + self.crys * 56 + 3 * 4 * self.free + 20
                else:
                    dsp = 4 + self.crys * 56 + 3 * 4 * self.natm + 20
                self.curr = num
                self.fdes.seek( self.head + dsp * num )
    
    
        def close( self ):
            self.fdes.close()
            # fix number of frames added...
            if( self.writ != "" ):
                self.fdes = open( self.writ, "rb+" )
                self.fdes.seek( 8 )
                self.fdes.write( struct.pack( "i", self.curr ) )
                self.fdes.seek( 20 )
                self.fdes.write( struct.pack( "i", self.curr ) )
                self.fdes.close()
    
    
        def open_write( self, fname, natoms, sele = None ):
            self.__init__()
            self.writ = fname
            self.curr = 0
            self.fdes = open( fname, "wb" )
            self.fdes.write( struct.pack( "i", 84 ) + b"CORD" )
            self.fdes.write( struct.pack( "i", -1 ) )
            self.fdes.write( struct.pack( "i", 1 ) )
            self.fdes.write( struct.pack( "i", 1 ) )
            self.fdes.write( struct.pack( "i", -1 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.natm = natoms
            if( sele != None and len( sele ) < self.natm ):
                self.sele = sele[:]
                self.free = len( sele )
                f = self.natm - self.free
            else:
                self.sele = []
                self.free = self.natm
                f = 0
            self.fdes.write( struct.pack( "i", 3 * self.free ) )
            self.fdes.write( struct.pack( "i", f ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 0 ) )
            self.fdes.write( struct.pack( "i", 84 ) )
            self.fdes.write( struct.pack( "i", 4 + 80 ) )
            self.fdes.write( struct.pack( "i", 1 ) )
            self.fdes.write( b"* Created with QMCube (qm3)                                                     " )
            self.fdes.write( struct.pack( "i", 4 + 80 ) )
            self.fdes.write( struct.pack( "i", 4 ) + struct.pack( "i", self.natm ) + struct.pack( "i", 4 ) )
            if( self.sele != [] ):
                b = struct.pack( "i", 4 * len( sele ) )
                self.fdes.write( b )
                for i in self.sele:
                    self.fdes.write( struct.pack( "i", i+1 ) )
                self.fdes.write( b )
    
    
        def append( self, molec ):
            try:
                self.fdes.tell()
            except AttributeError:
                print( "* No DCD open..." )
                return( False )
            else:
                if( self.natm != self.free and self.curr > 0 ):
                    t = struct.pack( "i", self.free * 4 )
                    bx = b""
                    by = b""
                    bz = b""
                    for i in self.sele:
                        i3 = 3 * i
                        bx += struct.pack( "f", molec.coor[i3]   )
                        by += struct.pack( "f", molec.coor[i3+1] )
                        bz += struct.pack( "f", molec.coor[i3+2] )
                else:
                    t = struct.pack( "i", self.natm * 4 )
                    bx = b""
                    by = b""
                    bz = b""
                    for i in range( self.natm ):
                        i3 = 3 * i
                        bx += struct.pack( "f", molec.coor[i3]   )
                        by += struct.pack( "f", molec.coor[i3+1] )
                        bz += struct.pack( "f", molec.coor[i3+2] )
                self.fdes.write( t + bx + t + t + by + t + t + bz + t )
                self.fdes.flush()
                self.curr += 1
                return( True )
    

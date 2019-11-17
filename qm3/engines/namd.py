# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import re
import struct
import math
import qm3.constants
import qm3.fio
import os, stat
import multiprocessing
import time


__tail = ""


def coordinates_read( mol, fname ):
    f = open( fname, "rb" )
    b = f.read( 4 )
    if( mol.natm == struct.unpack( "<i", b )[0] ):
        e = "<"
    elif( mol.natm == struct.unpack( ">i", b )[0] ):
        e = ">"
    else:
        print( "- Wrong atom number... (or broken!)" )
        f.close()
        return
#    mol.coor = list( struct.unpack( e+"%dd"%( mol.natm * 3 ), f.read( 8 * mol.natm * 3 ) ) )
    for i in range( mol.natm ):
        i3 = i * 3
        for j in [0, 1, 2]:
            mol.coor[i3+j] = struct.unpack( e+"d", f.read( 8 ) )[0]
    f.close()
    

def coordinates_write( mol, fname ):
    f = open( fname, "wb" )
    f.write( struct.pack( "i", mol.natm ) )
    for i in range( mol.natm ):
        i3 = i * 3
        for j in [0, 1, 2]:
            f.write( struct.pack( "d", mol.coor[i3+j] ) )
    f.close()


def pdb_write( mol, fname = None, fixed = [] ):
    f = qm3.fio.open_w( fname )
    s = [ 0 for i in range( mol.natm ) ]
    for i in fixed:
        s[i] = 1
    for i in range( mol.natm ):
        i3 = i * 3
        f.write( "ATOM  %5d %-5s%-5s%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf      %-4s\n"%( ( i % 99999 ) + 1, 
            " " * ( len( mol.labl[i] ) < 4 ) + mol.labl[i],
            mol.resn[i], mol.resi[i], mol.coor[i3], mol.coor[i3+1], mol.coor[i3+2], 
            s[i], 0.0, mol.segn[i] ) )
    f.write( "END\n" )
    qm3.fio.close( f, fname )


def topology_read( mol, fname = None ):
    global    __tail
    __tail = ""
    mol.type = []
    mol.chrg = []
    mol.mass = []
    fd = qm3.fio.open_r( fname )
    if( fd.readline().split()[0] == "PSF" ):
        fd.readline()
        for i in range( int( fd.readline().split()[0] ) + 1 ):
            fd.readline()
        if( mol.natm == int( fd.readline().split()[0] ) ):
            for i in range( mol.natm ):
                t = fd.readline().split()
                if( mol.segn[i] == t[1] and mol.resi[i] == int( t[2] ) and mol.resn[i] == t[3] and mol.labl[i] == t[4]  ):
                    mol.type.append( t[5] )
                    mol.chrg.append( float( t[6] ) )
                    mol.mass.append( float( t[7] ) )
                else:
                    print( "- Wrong data (%d): %s/%s %d/%s %s/%s %s/%s"%( i+1, mol.segn[i], t[1], mol.resi[i], t[2], mol.resn[i], t[3], mol.labl[i], t[4] ) )
            __tail = fd.read()
        else:
            print( "- Invalid number of atoms in PSF!" )
    qm3.fio.close( fd, fname )


def topology_write( mol, fname = None ):
    global    __tail
    fd = qm3.fio.open_w( fname )
    fd.write( "PSF\n\n       1 !NTITLE\n REMARKS generated structure x-plor psf file\n\n%8d !NATOM\n"%( mol.natm ) )
    for i in range( mol.natm ):
        fd.write( "%8d %-5s%-5d%-5s%-5s%-5s%10.6lf%14.4lf%12d\n"%( i + 1, mol.segn[i], mol.resi[i], mol.resn[i], mol.labl[i], mol.type[i], mol.chrg[i], mol.mass[i], 0 ) )
    fd.write( __tail )
    qm3.fio.close( fd, fname )




class namd_pipe( object ):

    def __init__( self ):
        self.pfd = open( "namd.pipe", "wt" )


    def stop( self ):
        self.pfd.write( "exit\n" )
        self.pfd.flush()
        self.pfd.close()


    def update_chrg( self, mol ):
        f = open( "namd.chrg", "wt" )
        f.write( str.join( "\n", [ "%12.6lf"%( i ) for i in mol.chrg ] ) )
        f.close()
        self.pfd.write( "charges\n" )
        self.pfd.flush()


    def update_coor( self, mol ):
# patched (ScriptTcl.C) --------------------------
        f = open( "namd.coor", "wb" )
        f.write( struct.pack( "i", mol.natm ) )
        for i in range( mol.natm ):
            i3 = i * 3
            for j in [0, 1, 2]:
                f.write( struct.pack( "f", mol.coor[i3+j] ) )
        f.close()
# ------------------------------------------------
        self.pfd.write( "coordinates\n" )
        self.pfd.flush()


    def get_func( self, mol ):
        self.update_coor( mol )
        self.pfd.write( "energy\n" )
        self.pfd.flush()
# patched (Controller.C) -------------------------
        while( not os.access( "namd.ener", os.R_OK ) ):
            time.sleep( 0.01 )
        while( os.stat( "namd.ener" )[stat.ST_SIZE] < 8 ):
            time.sleep( 0.01 )
        f = open( "namd.ener", "rb" )
        mol.func += struct.unpack( "d", f.read( 8 ) )[0] * qm3.constants.K2J
        f.close()
        os.unlink( "namd.ener" )
# ------------------------------------------------


    def get_grad( self, mol ):
        self.update_coor( mol )
        self.pfd.write( "gradient\n" )
        self.pfd.flush()
# patched (Controller.C) -------------------------
        while( not os.access( "namd.ener", os.R_OK ) ):
            time.sleep( 0.01 )
        while( os.stat( "namd.ener" )[stat.ST_SIZE] < 8 ):
            time.sleep( 0.01 )
        f = open( "namd.ener", "rb" )
        mol.func += struct.unpack( "d", f.read( 8 ) )[0] * qm3.constants.K2J
        f.close()
        os.unlink( "namd.ener" )
# ------------------------------------------------
        while( not os.access( "namd.force", os.R_OK ) ):
            time.sleep( 0.01 )
        while( os.stat( "namd.force" )[stat.ST_SIZE] < mol.natm * 24 + 4 ):
            time.sleep( 0.01 )
        f = open( "namd.force", "rb" )
        if( struct.unpack( "i", f.read( 4 ) )[0] == mol.natm ):
            g = [ -i * qm3.constants.K2J for i in struct.unpack( "%dd"%( mol.natm * 3 ), f.read( mol.natm * 24 ) ) ]
            for i in range( mol.natm ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    mol.grad[i3+j] += g[i3+j]
        f.close()
        os.unlink( "namd.force" )



try:
    import qm3.utils._shm
    class namd_shm( object ):

        def __init__( self ):
            f = open( "namd.shmid", "rt" )
            self.shm = int( f.read() )
            f.close()
            self.pfd = open( "namd.pipe", "wt" )
            self.clr = struct.pack( "d", 0.0 )
            self.eok = struct.pack( "d", 1.0 )
            self.gok = struct.pack( "d", 2.0 )


        def stop( self ):
            self.pfd.write( "exit\n" )
            self.pfd.flush()
            self.pfd.close()


        def update_chrg( self, mol ):
            qm3.utils._shm.write_r8( self.shm, [ 0 ] + mol.chrg )
            self.pfd.write( "charges\n" )
            self.pfd.flush()


        def update_coor( self, mol ):
            qm3.utils._shm.write_r8( self.shm, [ 0 ] + mol.coor )
            self.pfd.write( "coordinates\n" )
            self.pfd.flush()


        def get_func( self, mol ):
            self.update_coor( mol )
            qm3.utils._shm.write( self.shm, self.clr )
            self.pfd.write( "energy\n" )
            self.pfd.flush()
            while( qm3.utils._shm.read( self.shm, 8 ) != self.eok ):
                time.sleep( 0.01 )
            mol.func += qm3.utils._shm.read_r8( self.shm, 2 )[1] * qm3.constants.K2J


        def get_grad( self, mol ):
            self.update_coor( mol )
            qm3.utils._shm.write( self.shm, self.clr )
            self.pfd.write( "gradient\n" )
            self.pfd.flush()
            while( qm3.utils._shm.read( self.shm, 8 ) != self.gok ):
                time.sleep( 0.01 )
            tmp = qm3.utils._shm.read_r8( self.shm, 2 + 3 * mol.natm )
            mol.func += tmp[1] * qm3.constants.K2J
            for i in range( mol.natm ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    mol.grad[i3+j] -= tmp[2+i3+j] * qm3.constants.K2J
except:
    pass



class namd( object ):

    def __init__( self ):
        self.exe = "bash r.namd"


    def update_coor( self, mol ):
        coordinates_write( mol, "namd.coor" )


    def update_chrg( self, mol, fname ):
        topology_write( mol, fname )


    def get_func( self, mol ):
        self.update_coor( mol )
        os.system( self.exe )
        f = open( "namd.out", "rt" )
        mol.func += float( re.compile( "ENERGY:       0.*" ).findall( f.read() )[0].split()[11] ) * qm3.constants.K2J
        f.close()


    def get_grad( self, mol ):
        self.update_coor( mol )
        self.get_func( mol )
        f = open( "namd.force", "rb" )
        if( struct.unpack( "i", f.read( 4 ) )[0] == mol.natm ):
            g = [ -i * qm3.constants.K2J for i in struct.unpack( "%dd"%( mol.natm * 3 ), f.read( mol.natm * 24 ) ) ]
            for i in range( mol.natm ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    mol.grad[i3+j] += g[i3+j]
        f.close()



# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import os
import sys
import qm3.fio
import qm3.mol
import struct
import stat
import time
import qm3.elements
import math



def coordinates_read( fname = None ):
    """
Subsystem     1  A
   452 ! # of residues.
!===============================================================================
Residue     1  SER
    11 ! # of atoms.
     1   N             7       -5.0950000000     27.6770000000    -14.8700000000
     2   H             1       -4.2860000000     28.3550000000    -14.8320000000
     3   CA            6       -6.1260000000     27.5020000000    -13.9140000000
    """
    mol = qm3.mol.molecule()
    f = qm3.fio.open_r( fname )
    t = f.readline().strip().split()
    while( t != [] ):
        if( t[0].lower() == "subsystem" ):
            segn = t[2]
        if( t[0].lower() == "orthorhombic" ):
            mol.boxl = [ float( t[1] ), float( t[2] ), float( t[3] ) ]
        if( t[0].lower() == "cubic" ):
            mol.boxl = [ float( t[1] ), float( t[1] ), float( t[1] ) ]
        if( t[0].lower() == "residue" ):
            resi = int( t[1] )
            resn = t[2]
            t = f.readline().split()
            while( t[0][0] == "!" ):
                t = f.readline().split()
            for i in range( int( t[0] ) ):
                t = f.readline().split()
                while( t[0][0] == "!" ):
                    t = f.readline().split()
                mol.segn.append( segn )
                mol.resn.append( resn )
                mol.resi.append( resi )
                mol.labl.append( t[1] )
                mol.anum.append( int( t[2] ) )
                mol.coor.append( float( t[3] ) )
                mol.coor.append( float( t[4] ) )
                mol.coor.append( float( t[5] ) )
                mol.natm += 1
        t = f.readline().split()
    qm3.fio.close( f, fname )
    mol.settle()
    return( mol )


def coordinates_write( mol, fname = None ):
    f = qm3.fio.open_w( fname )
    f.write( "%d %d %d\n"%( mol.natm, len( mol.res_lim ) - 1, len( mol.seg_lim ) - 1 ) )
    if( mol.boxl != [ qm3.mol.MXLAT, qm3.mol.MXLAT, qm3.mol.MXLAT ]  ):
        f.write( "Symmetry  1\n" )
        if( mol.boxl[0] == mol.boxl[1] and mol.boxl[0] == mol.boxl[2] ):
            f.write( "CUBIC  %lf\n"%( mol.boxl[0] ) )
        else:
            f.write( "ORTHORHOMBIC  %lf %lf %lf\n"%( mol.boxl[0], mol.boxl[1], mol.boxl[2] ) )
    l = 1
    for i in range( len( mol.seg_lim ) - 1 ):
        f.write ( "Subsystem%6d  %s\n"%( i + 1, mol.segn[mol.res_lim[mol.seg_lim[i]]] ) )
        f.write ( "%6d\n"%( mol.seg_lim[i+1] - mol.seg_lim[i] ) )
        for j in range( mol.seg_lim[i], mol.seg_lim[i+1] ):
            f.write ( "Residue%6d  %s\n"%( mol.resi[mol.res_lim[j]], mol.resn[mol.res_lim[j]] ) )
            f.write ( "%6d\n"%( mol.res_lim[j+1] - mol.res_lim[j] ) )
            for k in range( mol.res_lim[j], mol.res_lim[j+1] ):
                f.write ( "%6d   %-10s%5d%20.10lf%20.10lf%20.10lf\n"%( l, mol.labl[k], int( mol.anum[k] ),
                    mol.coor[3*k], mol.coor[3*k+1], mol.coor[3*k+2] ) )
                l += 1
    qm3.fio.close( f, fname )


def topology_read( mol, fname ):
    mol.mass = []
    mol.chrg = []
    mol.epsi = []
    mol.rmin = []
    f = open( fname, "rb" )
    for i in range( 6 ):
        f.read( struct.unpack( "i", f.read( 4 ) )[0] + 4 )
    f.read( 4 )
    if( mol.natm == struct.unpack( "i", f.read( 4 ) )[0] ):
        f.read( 8 )
        mol.mass = list( struct.unpack( "%dd"%( mol.natm ), f.read( 8 * mol.natm ) ) )
        f.read( 4 )
        for i in range( 3 ):
            f.read( struct.unpack( "i", f.read( 4 ) )[0] + 4 )
        f.read( 4 )
        mol.chrg = list( struct.unpack( "%dd"%( mol.natm ), f.read( 8 * mol.natm ) ) )
        f.read( 4 )
        f.read( struct.unpack( "i", f.read( 4 ) )[0] + 4 )
        f.read( 4 )
        mol.epsi = [ i*0.5 for i in struct.unpack( "%dd"%( mol.natm ), f.read( 8 * mol.natm ) ) ]
        f.read( 4 )
        f.read( struct.unpack( "i", f.read( 4 ) )[0] + 4 )
        f.read( 4 )
        c = 0.5 * math.pow( 2, 1.0 / 6.0 )
        mol.rmin = [ i*i*c for i in struct.unpack( "%dd"%( mol.natm ), f.read( 8 * mol.natm ) ) ]
    f.close()


def sequence( mol, fname = None ):
    f = qm3.fio.open_w( fname )
    f.write( "Sequence\n%d\n\n"%( len( mol.seg_lim ) - 1 ) )
    for i in range( len( mol.seg_lim ) - 1 ):
        f.write ( "Subsystem  %s\n"%( mol.segn[mol.res_lim[mol.seg_lim[i]]] ) )
        f.write ( "%6d\n"%( mol.seg_lim[i+1] - mol.seg_lim[i] ) )
        k = 0
        for j in range( mol.seg_lim[i], mol.seg_lim[i+1] ):
            f.write ( " %s ;"%( mol.resn[mol.res_lim[j]] ) )
            k += 1
            if( k%12 == 0 ):
                f.write( "\n" )
        if( k%12 != 0 ):
            f.write( "\n" )
        f.write( "! Variant TYPE  RESN RESI\n" )
        f.write( "End\n\n" )
    f.write( "\n! Link TYPE  SEGN_A RESN_A RESI_A   SEGN_B RESN_B RESI_B\n" )
    f.write( "End" )
    qm3.fio.close( f, fname )


def selection( mol, sele, fname = None ):
    f = qm3.fio.open_w( fname )
    t = [ False for i in range( mol.natm ) ]
    for i in sele:
        t[i] = True
    f.write( """subroutine my_sele( selection )
    use atoms,             only : natoms
    logical, dimension(1:natoms), intent( inout ) :: selection
    selection = .false.
""" )
    for i in range( len( mol.res_lim ) - 1 ):
        if( sum( t[mol.res_lim[i]:mol.res_lim[i+1]] ) == mol.res_lim[i+1] - mol.res_lim[i] ):
            f.write( "\t! %s/%d\n"%( mol.segn[mol.res_lim[i]], mol.resi[mol.res_lim[i]] ) )
            f.write( "\tselection(%d:%d) = .true.\n"%( mol.res_lim[i]+1, mol.res_lim[i+1] ) )
        else:
            for j in range( mol.res_lim[i], mol.res_lim[i+1] ):
                if( t[j] ):
                    f.write( "\t! %s/%d/%s\n"%( mol.segn[j], mol.resi[j], mol.labl[j] ) )
                    f.write( "\tselection(%d) = .true.\n"%( j+1 ) )
    f.write( "end subroutine" )
    qm3.fio.close( f, fname )



try:
    import qm3.engines._dynamo
    class py_dynamo( qm3.engines._dynamo.dynamo ):
        def __init__( self, path = "./dynamo.so" ):
            qm3.engines._dynamo.dynamo.__init__( self, path )
except:
    pass



class dynamo_pipe( object ):

    def __init__( self ):
        self.pfd = open( "dynamo.pipe", "wt" )


    def stop( self ):
        self.pfd.write( "exit\n" )
        self.pfd.flush()
        self.pfd.close()
        try:
            os.unlink( self.inp )
        except:
            pass


    def update_chrg( self, mol ):
        f = open( "dynamo.chrg", "wb" )
        f.write( struct.pack( "i", 8 * mol.natm ) )
        for i in range( mol.natm ):
            f.write( struct.pack( "d", mol.chrg[i] ) )
        f.write( struct.pack( "i", 8 * mol.natm ) )
        f.close()
        self.pfd.write( "charges\n" )
        self.pfd.flush()


    def update_coor( self, mol ):
        f = open( "dynamo.crd", "wb" )
        f.write( struct.pack( "i", 24 * mol.natm ) )
        for i in range( 3 * mol.natm ):
            f.write( struct.pack( "d", mol.coor[i] ) )
        f.write( struct.pack( "i", 24 * mol.natm ) )
        f.close()
        self.pfd.write( "coordinates\n" )
        self.pfd.flush()


    def get_func( self, mol ):
        self.update_coor( mol )
        try:
            os.unlink( "dynamo.dat" )
        except:
            pass
        self.pfd.write( "energy\n" )
        self.pfd.flush()
        while( not os.path.isfile( "dynamo.dat" ) ):
            time.sleep( 0.01 )
        while( os.stat( "dynamo.dat" )[stat.ST_SIZE] < 16 ):
            time.sleep( 0.01 )
        f = open( "dynamo.dat", "rb" )
        f.read( 4 )
        mol.func += struct.unpack( "d", f.read( 8 ) )[0]
        f.close()


    def get_grad( self, mol ):
        self.update_coor( mol )
        try:
            os.unlink( "dynamo.dat" )
        except:
            pass
        self.pfd.write( "gradient\n" )
        self.pfd.flush()
        while( not os.path.isfile( "dynamo.dat" ) ):
            time.sleep( 0.01 )
        t = 24 + 24 * mol.natm
        while( os.stat( "dynamo.dat" )[stat.ST_SIZE] < t ):
            time.sleep( 0.01 )
        f = open( "dynamo.dat", "rb" )
        f.read( 4 )
        mol.func += struct.unpack( "d", f.read( 8 ) )[0]
        f.read( 8 )
        t = list( struct.unpack( "%dd"%( 3 * mol.natm ), f.read( 24 * mol.natm ) ) )
        f.close()
        for i in range( 3 * mol.natm ):
            mol.grad[i] += t[i]



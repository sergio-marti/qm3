# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import qm3.constants
import qm3.fio
import os
import re
import struct


__topo = ""
__frmt = re.compile( "[aAiIeEdD]([0-9]+)" )



def coordinates_read( mol, fname = None ):
    f = qm3.fio.open_r( fname )
    f.readline()
    n = int( f.readline().strip() )
    if( mol.natm != n ):
        print( "- Wrong atom number... (or broken!)" )
        f.close()
        return
    n *= 3
    mol.coor = []
    while( len( mol.coor ) < n ):
        mol.coor += [ float( j ) for j in f.readline().split() ]
    mol.boxl = [ float( j ) for j in f.readline().split()[0:3] ]
    qm3.fio.close( f, fname )


def coordinates_write( mol, fname = None ):
    f = qm3.fio.open_w( fname )
    f.write( "comment\n%d\n"%( mol.natm ) )
    n = 3 * mol.natm
    for i in range( n ):
        f.write( "%12.7lf"%( mol.coor[i] ) )
        if( (i+1)%6 == 0 ):
            f.write( "\n" )
    if( n%6 != 0 ):
        f.write( "\n" )
    f.write( "%12.7lf%12.7lf%12.7lf%12.7lf%12.7lf%12.7lf\n"%( mol.boxl[0], mol.boxl[1], mol.boxl[2], 90.0, 90.0, 90.0 ) )
    qm3.fio.close( f, fname )


def topology_read( mol, fname = None ):
    global    __topo, __frmt
    __topo = ""
    f = qm3.fio.open_r( fname )
    if( mol.natm > 0 ):
        l = f.readline()
        while( l != "" ):
            if( l[0:12].upper() == "%FLAG CHARGE" ):
                dsp = int( __frmt.findall( f.readline() )[0] )
                mol.chrg = []
                while( len( mol.chrg ) < mol.natm ):
                    l = f.readline()
                    mol.chrg += [ float( l[i:i+dsp] ) / 18.2223 for i in range( 0, len( l ) - 1, dsp ) ]
            else:
                __topo += l
            l = f.readline()
    else:
        mol.initialize()
        nres = 0
        lres = []
        l = f.readline()
        while( l != "" ):
            if( l[0:12].upper() == "%FLAG CHARGE" ):
                dsp = int( __frmt.findall( f.readline() )[0] )
                while( len( mol.chrg ) < mol.natm ):
                    l = f.readline()
                    mol.chrg += [ float( l[i:i+dsp] ) / 18.2223 for i in range( 0, len( l ) - 1, dsp ) ]
#                print( "chrg", mol.chrg[0], mol.chrg[-1], sum( mol.chrg ) )
            elif( l[0:14].upper() == "%FLAG POINTERS" ):
                __topo += l
                l = f.readline()
                __topo += l
                dsp = int( __frmt.findall( l )[0] )
                l = f.readline()
                __topo += l
                mol.natm = int( l[0:dsp] )
                l = f.readline()
                __topo += l
                nres = int( l[dsp:2*dsp] )
#                print( "natm", mol.natm, nres )
            elif( l[0:15].upper() == "%FLAG ATOM_NAME" ):
                __topo += l
                l = f.readline()
                __topo += l
                dsp = int( __frmt.findall( l )[0] )
                while( len( mol.labl ) < mol.natm ):
                    l = f.readline()
                    __topo += l
                    mol.labl += [ l[i:i+dsp].strip() for i in range( 0, len( l ) - 1, dsp ) ]
#                print( "labl", mol.labl[0], mol.labl[-1] )
            elif( l[0:21].upper() == "%FLAG AMBER_ATOM_TYPE" ):
                __topo += l
                l = f.readline()
                __topo += l
                dsp = int( __frmt.findall( l )[0] )
                while( len( mol.type ) < mol.natm ):
                    l = f.readline()
                    __topo += l
                    mol.type += [ l[i:i+dsp].strip() for i in range( 0, len( l ) - 1, dsp ) ]
#                print( "type", mol.type[0], mol.type[-1] )
            elif( l[0:19].upper() == "%FLAG ATOMIC_NUMBER" ):
                __topo += l
                l = f.readline()
                __topo += l
                dsp = int( __frmt.findall( l )[0] )
                while( len( mol.anum ) < mol.natm ):
                    l = f.readline()
                    __topo += l
                    mol.anum += [ int( l[i:i+dsp] ) for i in range( 0, len( l ) - 1, dsp ) ]
#                print( "anum", mol.anum[0], mol.anum[-1] )
            elif( l[0:10].upper() == "%FLAG MASS" ):
                __topo += l
                l = f.readline()
                __topo += l
                dsp = int( __frmt.findall( l )[0] )
                while( len( mol.mass ) < mol.natm ):
                    l = f.readline()
                    __topo += l
                    mol.mass += [ float( l[i:i+dsp] ) for i in range( 0, len( l ) - 1, dsp ) ]
#                print( "mass", mol.mass[0], mol.mass[-1], sum( mol.mass ) )
            elif( l[0:19].upper() == "%FLAG RESIDUE_LABEL" ):
                __topo += l
                l = f.readline()
                __topo += l
                dsp = int( __frmt.findall( l )[0] )
                while( len( lres ) < nres ):
                    l = f.readline()
                    __topo += l
                    lres += [ l[i:i+dsp].strip() for i in range( 0, len( l ) - 1, dsp ) ]
#                print( "resn", lres[0], lres[-1] )
            elif( l[0:21].upper() == "%FLAG RESIDUE_POINTER" ):
                __topo += l
                l = f.readline()
                __topo += l
                dsp = int( __frmt.findall( l )[0] )
                while( len( mol.res_lim ) < nres ):
                    l = f.readline()
                    __topo += l
                    mol.res_lim += [ int( l[i:i+dsp] ) - 1 for i in range( 0, len( l ) - 1, dsp ) ]
                mol.res_lim.append( mol.natm )
#                print( "rlim", len( mol.res_lim ) )
            elif( l[0:20].upper() == "%FLAG BOX_DIMENSIONS" ):
                __topo += l
                l = f.readline()
                __topo += l
                dsp = int( __frmt.findall( l )[0] )
                l = f.readline()
                __topo += l
                mol.boxl = [ float( l[i:i+dsp] ) for i in range( 0, len( l ) - 1, dsp ) ][1:4]
#                print( "boxl", mol.boxl )
            else:
                __topo += l
            l = f.readline()
        mol.segn = [ "A" for i in range( mol.natm ) ]
        mol.seg_lim = [ 0, mol.natm ]
        for i in range( nres ):
            for j in range( mol.res_lim[i], mol.res_lim[i+1] ):
                mol.resi.append( i + 1 )
                mol.resn.append( lres[i] )
    qm3.fio.close( f, fname )


def topology_write( mol, fname = None ):
    global    __topo
    if( __topo != "" and mol.chrg != [] ):
        f = qm3.fio.open_w( fname )
        f.write( __topo )
        f.write( "%FLAG CHARGE                                                                    \n" )
        f.write( "%FORMAT(5E16.8)                                                                 \n" )
        for i in range( mol.natm ):
            f.write( "%16.8lf"%( mol.chrg[i] * 18.2223 ) )
            if( (i+1)%5 == 0 ):
                f.write( "\n" )
        qm3.fio.close( f, fname )



try:
    import qm3.engines._sander
    class py_sander( qm3.engines._sander.sander ):
        def __init__( self, mol, prmtop, cutoff = 10.0, PBC = True, qmsel = None, method = "AM1", charge = 0 ):
            qm3.engines._sander.sander.__init__( self, mol, prmtop, cutoff, PBC, qmsel, method, charge )
except:
    pass



class sander( object ):

    def __init__( self ):
        self.exe = "bash r.sander"


    def update_coor( self, mol ):
        coordinates_write( mol, "inpcrd" )


    def update_chrg( self, mol ):
        topology_write( mol, "prmtop" )


    def get_func( self, mol ):
        self.update_coor( mol )
        os.system( self.exe )
        f = open( "mden", "rt" )
        l = f.readline()
        k = True
        while( l and k ):
            if( l[0:2] == "L6" ):
                try:
                    mol.func += qm3.constants.K2J * float( l.strip().split()[2] )
                    k = False
                except:
                    pass
            l = f.readline()
        f.close()


    def get_grad( self, mol ):
        self.update_coor( mol )
        self.get_func( mol )
        # dirty NETCDF wrapper...
        f = open( "mdfrc", "rb" )
        b = f.read()
        f.close()
        o_x = ord( "x" ); o_y = ord( "y" ); o_z = ord( "z" )
        n = len( b )
        w = 0
        f = True
        while( w < n-2 and f ):
            # bytes (python3) / str (python2) compliant...
            if( ( b[w] == o_x and b[w+1] == o_y and b[w+2] == o_z ) or ( b[w] == "x" and b[w+1] == "y" and b[w+2] == "z" ) ):
                w += 7
                f  = False
            w += 1
        g = struct.unpack( ">%df"%( mol.natm * 3 ), b[w:] )
        for i in range( mol.natm ):
            i3 = i * 3
            for j in [0, 1, 2]:
                mol.grad[i3+j] -= g[i3+j] * qm3.constants.K2J




# ------------------------------------------------------------------------------------
# - ORCA_FAKE - ("orca" in current folder, for UNMODIFIED sander versions...)
# ------------------------------------------------------------------------------------
##!/usr/bin/env python
#
#f = open( "inpfile.xyz", "rt" )
#nQM = int( f.readline().strip() )
#f.close()
#f = open( "orc_job.engrad", "wt" )
#f.write( """#
## Number of atoms
##
# %d
##
## The current total energy in Eh
##
#      0.000000000000
##
## The current gradient in Eh/bohr
##
#"""%( nQM ) )
#for i in range( nQM ):
#    f.write( "       0.000000000000\n       0.000000000000\n       0.000000000000\n" )
#f.close()
#
#f = open( "ptchrg.xyz", "rt" )
#nMM = int( f.readline().strip() )
#f.close()
#f = open( "orc_job.pcgrad", "wt" )
#f.write( "%d\n"%( nMM ) )
#for i in range( nMM ):
#    f.write( "   0.000000000000   0.000000000000   0.000000000000\n" )
#f.close()
#
#f = open( "orc_job.dat", "wt" )
#f.write( "FINAL SINGLE POINT ENERGY         0.000000000000\n" )
#f.write( "Total Dipole Moment    :      0.00000       0.00000       0.00000\n" )
#f.write( "                        -----------------------------------------\n" )
#f.write( "Magnitude (a.u.)       :      0.00000\n" )
#f.close()
# ------------------------------------------------------------------------------------

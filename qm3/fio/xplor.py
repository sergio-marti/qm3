# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.mol
import qm3.fio


def psf_read( mol, fname ):
    mol.type = []
    mol.chrg = []
    mol.mass = []
    f = qm3.fio.open_r( fname )
    if( f.readline().split()[0] == "PSF" ):
        f.readline()
        for i in range( int( f.readline().split()[0] ) + 1 ):
            f.readline()
        n = int( f.readline().split()[0] )
        if( mol.natm == 0 ):
            mol.natm = n
            mol.segn = []
            mol.resi = []
            mol.resn = []
            mol.labl = []
            for i in range( n ):
                t = f.readline().split()
                mol.segn.append( t[1] )
                mol.resi.append( int( t[2] ) )
                mol.resn.append( t[3] )
                mol.labl.append( t[4] )
                mol.type.append( t[5] )
                mol.chrg.append( float( t[6] ) )
                mol.mass.append( float( t[7] ) )
        elif( mol.natm == n ):
            for i in range( mol.natm ):
                t = f.readline().split()
                if( mol.segn[i] == t[1] and mol.resi[i] == int( t[2] ) and mol.resn[i] == t[3] and mol.labl[i] == t[4]  ):
                    mol.type.append( t[5] )
                    mol.chrg.append( float( t[6] ) )
                    mol.mass.append( float( t[7] ) )
                else:
                    print( "- Wrong PSF data (%d): %s/%s %d/%s %s/%s %s/%s"%( i+1, mol.segn[i], t[1], mol.resi[i], t[2], mol.resn[i], t[3], mol.labl[i], t[4] ) )
                    mol.type = []
                    mol.chrg = []
                    mol.mass = []
                    return( [] )
        else:
            print( "- Invalid number of atoms in PSF!" )
            return( [] )
    bonds = []
    f.readline()
    n = int( f.readline().strip().split()[0] )
    while( len( bonds ) < n ):
        t = [ int( i ) - 1 for i in f.readline().strip().split() ]
        for i in range( len( t ) // 2 ):
            bonds.append( [ t[2*i], t[2*i+1] ] )
    qm3.fio.close( f, fname )
    return( bonds )

# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.mol
import qm3.fio


__topo = ""


def psf_read( mol, fname ):
    global  __topo
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
            mol.settle()
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
    bonds  = []
    __topo = f.readline()
    l = f.readline()
    n = int( l.strip().split()[0] )
    __topo += l
    while( len( bonds ) < n ):
        l = f.readline()
        t = [ int( i ) - 1 for i in l.strip().split() ]
        for i in range( len( t ) // 2 ):
            bonds.append( [ t[2*i], t[2*i+1] ] )
        __topo += l
    __topo += f.read()
    qm3.fio.close( f, fname )
    return( bonds )



def psf_write( mol, fname = None ):
    global    __topo
    fd = qm3.fio.open_w( fname )
    if( sum( [ mol.seg_lim[i] - mol.seg_lim[i-1] >= 10000 for i in range( 1, len( mol.seg_lim ) ) ] ) > 0 ):
        # expanded format EXT:
        # II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I) == 0
        # (I10,1X, [A8,1X], [A8,1X], [A8,1X], [A8,1X], [A4,1X], 2G14.6,I8) XPLOR
        fd.write( "PSF EXT\n\n       1 !NTITLE\n REMARKS generated structure x-plor extended psf file\n\n%8d !NATOM\n"%( mol.natm ) )
        for i in range( mol.natm ):
            fd.write( "%8d %-9s%-9d%-9s%-9s%-5s%14.6lf%14.6lf%8d\n"%( i + 1, mol.segn[i], mol.resi[i], mol.resn[i], mol.labl[i],
                mol.type[i], mol.chrg[i], mol.mass[i], 0 ) )
    else:
        # standard format:
        # II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I) == 0
        # (I8,1X, [A4,1X], [A4,1X], [A4,1X], [A4,1X], [A4,1X], 2G14.6,I8)  XPLOR
        fd.write( "PSF\n\n       1 !NTITLE\n REMARKS generated structure x-plor psf file\n\n%8d !NATOM\n"%( mol.natm ) )
        for i in range( mol.natm ):
            fd.write( "%8d %-5s%-5d%-5s%-5s%-5s%14.6lf%14.6lf%8d\n"%( i + 1, mol.segn[i], mol.resi[i], mol.resn[i], mol.labl[i],
                mol.type[i], mol.chrg[i], mol.mass[i], 0 ) )
    fd.write( __topo )
    qm3.fio.close( fd, fname )


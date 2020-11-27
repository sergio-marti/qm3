# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import os
import re
import math
import struct
import qm3.constants
import qm3.elements
import qm3.engines


def lio_input( obj, mol ):
    f = open( "lio.xyz", "wt" )
    j = 0
    for i in obj.sel:
        i3 = i * 3
        f.write( "%4d%20.10lf%20.10lf%20.10lf\n"%( qm3.elements.rsymbol[obj.smb[j]],
            mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
            mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
            mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
        j += 1
    if( len( obj.lnk ) > 0 ):
        obj.vla = []
        k = len( obj.sel )
        for i,j in obj.lnk:
            c, v = qm3.engines.LA_coordinates( i, j, mol )
            # To allow the interaction of the Link-Atom with the environment change atomic number to "1"
            f.write( "%4d%20.10lf%20.10lf%20.10lf\n"%( 1, c[0], c[1], c[2] ) )
            obj.vla.append( ( obj.sel.index( i ), k, v[:] ) )
            k += 1
    if( len( obj.nbn ) > 0 ):
        for i in obj.nbn:
            i3 = i * 3
            f.write( "%12.4lf%20.10lf%20.10lf%20.10lf\n"%( mol.chrg[i],
                mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
    f.close()
    f = open( "lio.inp", "wt" )
    buf = obj.inp.replace( "qm3_natm", str( len( obj.sel ) + len( obj.lnk ) ) )
    buf = buf.replace( "qm3_nchg", str( len( obj.nbn ) ) )
    f.write( buf )
    f.close()



class lio( qm3.engines.qmbase ):

    def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
        qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
        self.exe = "bash r.lio"


    def mk_input( self, mol, run ):
        lio_input( self, mol )


    def parse_log( self, mol, run ):
        f = open( "lio.log", "rt" )
        mol.func += self._ce * float( re.compile( "Total energy =[\ ]*([0-9\.\-]+)" ).findall( f.read() )[0] )
        f.close()
        if( run == "grad" ):
            f = open( "lio.grad", "rt" )
            g = []
            for i in range( len( self.sel ) + len( self.lnk ) ):
                g += [ float( j ) * self._cg for j in f.readline().strip().split()[1:4] ]
            qm3.engines.LA_gradient( self.vla, g )
            for i in range( len( self.sel ) ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    mol.grad[3*self.sel[i]+j] += g[i3+j]
            if( len( self.nbn ) > 0 ):
                for i in self.nbn:
                    i3 = i * 3
                    g = [ float( j ) * self._cg for j in f.readline().strip().split()[1:4] ]
                    for j in [0, 1, 2]:
                        mol.grad[i3+j] += g[j]
            f.close()



try:
    import ctypes
    class dl_lio( qm3.engines.qmbase ):
    
        def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
            qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )

            self.nQM = len( self.sel ) + len( self.lnk )
            self.siz = 1 + 3 * ( self.nQM + len( self.nbn ) ) + self.nQM
            self.vec = ( ctypes.c_double * self.siz )()
            self.lib = ctypes.CDLL( os.getenv( "QM3_LIBLIO" ) )
            self.lib.qm3_lio_calc_.argtypes = [ ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_double ) ]
            self.lib.qm3_lio_calc_.restype = None

            lio_input( self, mol )
            self.lib.qm3_lio_init_()


        def update_coor( self, mol ):
            for i in range( len( self.sel ) ):
                i3 = self.sel[i] * 3
                j3 = i * 3
                for j in [0, 1, 2]:
                    self.vec[j3+j] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) 
            self.vla = []
            k = len( self.sel )
            for i in range( len( self.lnk ) ):
                j3 = k * 3
                c, v = qm3.engines.LA_coordinates( self.lnk[i][0], self.lnk[i][1], mol )
                for j in [0, 1, 2]:
                    self.vec[j3+j] = c[j]
                self.vla.append( ( self.sel.index( self.lnk[i][0] ), k, v[:] ) )
                k += 1
            for i in range( len( self.nbn ) ):
                i3 = self.nbn[i] * 3
                j3 = ( self.nQM + i ) * 3
                for j in [0, 1, 2]:
                    self.vec[j3+j] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) 


        def get_func( self, mol ):
            self.update_coor( mol )
            self.lib.qm3_lio_calc_( ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0] * self._ce
            for i in range( len( self.sel ) ):
                mol.chrg[self.sel[i]] = self.vec[i+1]


        def get_grad( self, mol ):
            self.update_coor( mol )
            self.lib.qm3_lio_calc_( ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0] * self._ce
            for i in range( len( self.sel ) ):
                mol.chrg[self.sel[i]] = self.vec[i+1]
            g = [ j * self._cg for j in self.vec[self.nQM+1:] ]
            qm3.engines.LA_gradient( self.vla, g )
            for i in range( len( self.sel ) ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    mol.grad[3*self.sel[i]+j] += g[i3+j]
            for i in range( len( self.nbn ) ):
                i3 = ( self.nQM + i ) * 3
                for j in [0, 1, 2]:
                    mol.grad[3*self.nbn[i]+j] += g[i3+j]

except:
    pass

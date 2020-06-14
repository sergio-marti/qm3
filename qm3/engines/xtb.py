# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import os
import re
import qm3.engines
import qm3.constants
import qm3.elements



class xtb( qm3.engines.qmbase ):

    def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
        qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
        self.exe = "bash r.xtb"
        f = open( "xtb.inp", "wt" )
        f.write( self.inp )
        f.close()


    def mk_input( self, mol ):
        f = open( "xtb_qm.xyz", "wt" )
        f.write( "%d\n\n"%( len( self.sel ) + len( self.lnk ) ) )
        j = 0
        for i in self.sel:
            i3 = i * 3
            f.write( "%2s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j], 
                mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
            j += 1
        if( self.lnk ):
            self.vla = []
            k = len( self.sel )
            for i,j in self.lnk:
                c, v = qm3.engines.LA_coordinates( i, j, mol )
                f.write( "%-2s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] ) )
                self.vla.append( ( self.sel.index( i ), k, v[:] ) )
                k += 1
        f.close()
        if( self.nbn ):
            f = open( "xtb_mm", "wt" )
            f.write( "%d\n"%( len( self.nbn ) ) )
            for i in self.nbn:
                i3 = i * 3
                f.write( "%12.4lf%20.10lf%20.10lf%20.10lf%12d\n"%( mol.chrg[i],
                    mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3] / mol.boxl[0], 0 ), 
                    mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                    mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), 99 ) )
            f.close()


    def parse_log( self, mol, run ):
        f = open( "energy", "rt" )
        mol.func += float( re.compile( "[0-9\.\-]+" ).findall( f.read() )[1] ) * self._ce
        f.close()
        if( run == "grad" ):
            f = open( "gradient", "rt" )
            g = [ float( i ) * self._cg for i in re.compile( "[0-9\.\-]+E[0-9\.\-]+" ).findall( f.read() ) ]
            f.close()
            qm3.engines.LA_gradient( self.vla, g )
            for i in range( len( self.sel ) ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    mol.grad[3*self.sel[i]+j] += g[i3+j]
            if( self.nbn and os.access( "pcgrad", os.R_OK ) ):
                f = open( "pcgrad", "rt" )
                g = [ float( i ) * self._cg for i in re.compile( "[0-9\.\-]+" ).findall( f.read() ) ]
                f.close()
                for i in range( len( self.nbn ) ):
                    i3 = i * 3
                    for j in [0, 1, 2]:
                        mol.grad[3*self.nbn[i]+j] += g[i3+j]


    def get_func( self, mol ):
        self.mk_input( mol )
        os.system( self.exe )
        self.parse_log( mol, "ener" )


    def get_grad( self, mol ):
        self.mk_input( mol )
        os.system( self.exe )
        self.parse_log( mol, "grad" )




try:
    import ctypes

    class dl_xtb:
    
        def __init__( self, mol, chrg, sele, nbnd = [], link = [] ):
            self._cx = qm3.constants.A0
            self._ce = qm3.constants.H2J
            self._cg = self._ce / qm3.constants.A0
            self._ch = self._cg / qm3.constants.A0

            self.sel = sorted( sele )
            self.lnk = link[:]
            self.nbn = sorted( set( nbnd ).difference( set( sele + sum( link, [] ) ) ) )
            self.vla = []

            if( not mol.chrg ):
                mol.chrg = [ 0.0 for i in range( mol.natm ) ]

            self.nQM = len( self.sel ) + len( self.lnk )
            self.nMM = len( self.nbn )
            # 2 + nQM [QM_chg] + 3 * nQM [QM_crd/grd] + nQM [QM_mul] + nMM [MM_chg] + 3 * nMM [MM_crd/grd]
            self.siz = 2 + 5 * self.nQM + 4 * self.nMM
            self.vec = ( ctypes.c_double * self.siz )()

            self.vec[1] = chrg
            j = 2
            for i in self.sel:
                self.vec[j] = mol.anum[i]
                j += 1
            for i in range( len( self.lnk ) ):
                self.vec[j] = 1
                j += 1
            j  = 2 + 5 * self.nQM
            for i in range( self.nMM ):
                self.vec[j+i] = mol.chrg[self.nbn[i]]
            
            self.lib = ctypes.CDLL( os.getenv( "QM3_LIBXTB" ) )
            self.lib.xtb_calc_.argtypes = [ 
                ctypes.POINTER( ctypes.c_int ),
                ctypes.POINTER( ctypes.c_int ),
                ctypes.POINTER( ctypes.c_int ),
                ctypes.POINTER( ctypes.c_double ) ]
            self.lib.xtb_calc_.restype = None

    
        def update_coor( self, mol ):
            j3 = 2 + self.nQM
            for i in range( len( self.sel ) ):
                i3 = self.sel[i] * 3
                for j in [0, 1, 2]:
                    self.vec[j3] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 )
                    j3 += 1
            self.vla = []
            k = len( self.sel )
            for i in range( len( self.lnk ) ):
                c, v = qm3.engines.LA_coordinates( self.lnk[i][0], self.lnk[i][1], mol )
                for j in [0, 1, 2]:
                    self.vec[j3] = c[j]
                    j3 += 1
                self.vla.append( ( self.sel.index( self.lnk[i][0] ), k, v[:] ) )
                k += 1
            j3 = 2 + 5 * self.nQM + self.nMM
            for i in range( self.nMM ):
                i3 = self.nbn[i] * 3
                for j in [0, 1, 2]:
                    self.vec[j3] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 )
                    j3 += 1


        def get_func( self, mol ):
            self.update_coor( mol )
            self.lib.xtb_calc_( ctypes.c_int( self.nQM ), ctypes.c_int( self.nMM ), ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0] * self._ce
            k = 2 + 4 * self.nQM
            for i in range( len( self.sel ) ):
                mol.chrg[self.sel[i]] = self.vec[k+i]
    
    
        def get_grad( self, mol ):
            self.update_coor( mol )
            self.lib.xtb_calc_( ctypes.c_int( self.nQM ), ctypes.c_int( self.nMM ), ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0] * self._ce
            k = 2 + 4 * self.nQM
            for i in range( len( self.sel ) ):
                mol.chrg[self.sel[i]] = self.vec[k+i]
            k = 2 + self.nQM
            g = [ self.vec[k+j] * self._cg for j in range( 3 * self.nQM ) ]
            qm3.engines.LA_gradient( self.vla, g )
            j3 = 0
            for i in range( len( self.sel ) ):
                i3 = 3 * self.sel[i]
                for j in [0, 1, 2]:
                    mol.grad[i3+j] += g[j3]
                    j3 += 1
            j3 = 2 + 5 * self.nQM + self.nMM
            for i in range( self.nMM ):
                i3 = 3 * self.nbn[i]
                for j in [0, 1, 2]:
                    mol.grad[i3+j] += self.vec[j3] * self._cg
                    j3 += 1

except:
    pass

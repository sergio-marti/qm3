# -*- coding: iso-8859-1 -*-
import os
import re
import qm3.engines
import qm3.elements

try:
    import ctypes

    class dl_mopac:
    
        def __init__( self, mol, meth, chrg, mult, sele, nbnd = [], link = [], con = -1, cof = -1 ):
            self.sel = sorted( sele )
            self.lnk = link[:]
            self.nbn = sorted( set( nbnd ).difference( set( sele + sum( link, [] ) ) ) )
            self.vla = []

            hami = { "MNDO": 0, "AM1": 1, "RM1": 2, "PM3": 3, "PDDG": 4 }

            if( not mol.chrg ):
                mol.chrg = [ 0.0 for i in range( mol.natm ) ]

            self.nQM = len( self.sel ) + len( self.lnk )
            self.nMM = len( self.nbn )
            # 1 + 3 * nQM [QM_crd/grd] + nQM [QM_mul] + nMM [MM_chg] + 3 * nMM [MM_crd/grd]
            self.siz = 1 + 4 * self.nQM + 4 * self.nMM
            self.vec = ( ctypes.c_double * self.siz )()

            self.lib = ctypes.CDLL( os.getenv( "QM3_LIBMOPAC" ) )

            self.lib.qm3_mopac_setup_.argtypes = [ ctypes.POINTER( ctypes.c_int ),
                ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_int ),
                ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_int ),
                ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_double ),
                ctypes.POINTER( ctypes.c_double ), ctypes.POINTER( ctypes.c_double ) ]
            self.lib.qm3_mopac_setup_.restype = None

            self.lib.qm3_mopac_calc_.argtypes = [ ctypes.POINTER( ctypes.c_int ),
                ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_double ) ]
            self.lib.qm3_mopac_calc_.restype = None

            j = 1
            for i in self.sel:
                self.vec[j] = mol.anum[i]
                j += 1
            for i in range( len( self.lnk ) ):
                self.vec[j] = 1
                j += 1
            self.lib.qm3_mopac_setup_( ctypes.c_int( self.nQM ), ctypes.c_int( self.nMM ),
                    ctypes.c_int( hami[meth] ),
                    ctypes.c_int( chrg ),
                    ctypes.c_int( mult ),
                    ctypes.c_int( self.siz ), self.vec,
                    ctypes.c_double( con ), ctypes.c_double( cof ) )


        def stop( self ):
            self.lib.qm3_mopac_clean_()

    
        def update_coor( self, mol ):
            j3 = 1
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
            j3 = 1 + 4 * self.nQM + self.nMM
            for i in range( self.nMM ):
                i3 = self.nbn[i] * 3
                for j in [0, 1, 2]:
                    self.vec[j3] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 )
                    j3 += 1
            j  = 1 + 4 * self.nQM
            for i in range( self.nMM ):
                self.vec[j+i] = mol.chrg[self.nbn[i]]


        def get_func( self, mol, maxit = 200 ):
            self.update_coor( mol )
            self.lib.qm3_mopac_calc_( ctypes.c_int( maxit ), ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0]
            k = 1 + 3 * self.nQM
            for i in range( len( self.sel ) ):
                mol.chrg[self.sel[i]] = self.vec[k+i]


        def get_grad( self, mol, maxit = 200 ):
            self.update_coor( mol )
            self.lib.qm3_mopac_calc_( ctypes.c_int( maxit ), ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0]
            k = 1 + 3 * self.nQM
            for i in range( len( self.sel ) ):
                mol.chrg[self.sel[i]] = self.vec[k+i]
            g = [ self.vec[j] for j in range( 1, 3 * self.nQM + 1 ) ]
            qm3.engines.LA_gradient( self.vla, g )
            j3 = 0
            for i in range( len( self.sel ) ):
                i3 = 3 * self.sel[i]
                for j in [0, 1, 2]:
                    mol.grad[i3+j] += g[j3]
                    j3 += 1
            j3 = 1 + 4 * self.nQM + self.nMM
            for i in range( self.nMM ):
                i3 = 3 * self.nbn[i]
                for j in [0, 1, 2]:
                    mol.grad[i3+j] += self.vec[j3]
                    j3 += 1

except:
    pass

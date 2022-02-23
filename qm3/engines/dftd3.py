# -*- coding: iso-8859-1 -*-
import os
import re
import qm3.constants



class run_single( object ):

    def __init__( self, mol, sele ):
        self.exe = "bash r.dftd3"
        self._ce = qm3.constants.H2J
        self._cg = self._ce / qm3.constants.A0
        self.sel = sorted( sele )
        self.smb = mol.guess_symbols( self.sel )
        self.pat = re.compile( "Edisp /kcal,au:[\ ]+[0-9\.\-]+[\ ]+([0-9\.\-]+)" )


    def mk_input( self, mol ):
        f = open( "dftd3.xyz", "wt" )
        f.write( "%d\n\n"%( len( self.sel ) ) )
        j = 0
        for i in self.sel:
            i3 = i * 3
            f.write( "%2s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j],
                mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
            j += 1
        f.close()


    def parse_log( self, mol, run ):
        f = open( "dftd3.log", "rt" )
        mol.func += float( self.pat.findall( f.read() )[0] ) * self._ce
        f.close()
        if( run == "grad" and os.path.isfile( "dftd3_gradient" ) ):
            f = open( "dftd3_gradient", "rt" )
            g = [ float( i ) * self._cg for i in f.read().split() ]
            for i in range( len( self.sel ) ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    mol.grad[3*self.sel[i]+j] += g[i3+j]
            f.close()


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

    class run_dynlib( object ):
    
        def __init__( self, mol, parm, sele ):
            self._ce = qm3.constants.H2J
            self._cg = self._ce / qm3.constants.A0
            self.sel = sorted( sele )
            self.anu = [ mol.anum[i] for i in self.sel ]
            self.nat = len( self.sel )
            self.siz = self.nat * 4
            self.vec = ( ctypes.c_double * self.siz )()
            self.lib = ctypes.CDLL( os.getenv( "QM3_LIBDFTD3" ) )
            self.lib.qm3_dftd3_init_.argtypes = [ ctypes.POINTER( ctypes.c_double ) ]
            self.lib.qm3_dftd3_init_.restype = None
            self.lib.qm3_dftd3_calc_.argtypes = [ ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_double ) ]
            self.lib.qm3_dftd3_calc_.restype = None
            tmp = ( ctypes.c_double * 6 )()
            # -------------------------------------------------
            # see subroutine setfuncpar at core.f90 (dftd3-lib)
            tmp[0] = parm["version"]
            tmp[1] = parm["s6"]
            tmp[2] = parm["rs6"]
            tmp[3] = parm["s18"]
            tmp[4] = parm["rs18"]
            tmp[5] = parm["alp"]
            # -------------------------------------------------
            self.lib.qm3_dftd3_init_( tmp )
    
    
        def update_coor( self, mol ):
            for i in range( self.nat ):
                self.vec[i] = self.anu[i]
            for i in range( self.nat ):
                i3 = self.sel[i] * 3
                j3 = self.nat + i * 3
                for j in [0, 1, 2]:
                    self.vec[j3+j] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) 


        def get_func( self, mol ):
            self.update_coor( mol )
            self.lib.qm3_dftd3_calc_( ctypes.c_int( self.nat ), ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0] * self._ce
    
    
        def get_grad( self, mol ):
            self.update_coor( mol )
            self.lib.qm3_dftd3_calc_( ctypes.c_int( self.nat ), ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0] * self._ce
            for i in range( self.nat ):
                i3 = 1 + 3 * i
                j3 = 3 * self.sel[i]
                for j in [0, 1, 2]:
                    mol.grad[j3+j] += self.vec[i3+j] * self._cg

except:
    pass

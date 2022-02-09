# -*- coding: iso-8859-1 -*-
import os
import re
import qm3.constants
import qm3.engines



class run_single():

    def __init__( self, mol, sele, link = [] ):
        self.exe = "bash r.dftd4"
        self._ce = qm3.constants.H2J
        self._cg = self._ce / qm3.constants.A0
        self.sel = sorted( sele )
        self.smb = mol.guess_symbols( self.sel )
        self.pat = re.compile( "Edisp[\ ]+/kcal,au:[\ ]+[0-9\.\-]+[\ ]+([0-9\.\-]+)" )
        self.lnk = link[:]
        self.vla = []
        self.nat = len( self.sel )
        self.all = self.nat + len( self.lnk )


    def mk_input( self, mol ):
        f = open( "dftd4.xyz", "wt" )
        f.write( "%d\n\n"%( self.all ) )
        j = 0
        for i in self.sel:
            i3 = i * 3
            f.write( "%2s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j],
                mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
            j += 1
        if( len( self.lnk ) > 0 ):
            self.vla = []
            k = self.nat
            for i,j in self.lnk:
                c, v = qm3.engines.LA_coordinates( i, j, mol )
                f.write( "%-2s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] ) )
                self.vla.append( ( self.sel.index( i ), k, v[:] ) )
                k += 1
        f.close()


    def parse_log( self, mol, run ):
        f = open( "dftd4.log", "rt" )
        mol.func += float( self.pat.findall( f.read() )[0] ) * self._ce
        f.close()
        if( run == "grad" and os.path.isfile( "gradient" ) ):
            g = []
            f = open( "gradient", "rt" )
            for i in range( 2 + self.all ):
                f.readline()
            for i in range( self.all ):
                g += [ float( j.replace( "D", "E" ) ) * self._cg for j in f.readline().strip().split() ]
            f.close()
            qm3.engines.LA_gradient( self.vla, g )
            for i in range( self.nat ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    mol.grad[3*self.sel[i]+j] += g[i3+j]


    def get_func( self, mol ):
        self.mk_input( mol )
        os.system( self.exe )
        self.parse_log( mol, "ener" )


    def get_grad( self, mol ):
        self.mk_input( mol )
        try:
            os.unlink( "gradient" )
        except:
            pass
        os.system( self.exe )
        self.parse_log( mol, "grad" )




try:
    import ctypes

    class run_dynlib():
    
        def __init__( self, mol, parm, sele, link = [] ):
            self._ce = qm3.constants.H2J
            self._cg = self._ce / qm3.constants.A0
            self.sel = sorted( sele )
            self.lnk = link[:]
            self.vla = []
            self.nat = len( self.sel )
            self.all = self.nat + len( self.lnk )
            self.siz = max( 6, self.all * 3 + 1 )
            self.vec = ( ctypes.c_double * self.siz )()
            self.lib = ctypes.CDLL( os.getenv( "QM3_LIBDFTD4" ) )
            self.lib.qm3_dftd4_init_.argtypes = [ ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_double ) ]
            self.lib.qm3_dftd4_init_.restype = None
            self.lib.qm3_dftd4_calc_.argtypes = [ ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_double ) ]
            self.lib.qm3_dftd4_calc_.restype = None
            # -------------------------------------------------
            self.vec[0] = parm["chrg"]
            k = 1
            for i in range( self.nat ):
                self.vec[k] = mol.anum[self.sel[i]]
                k += 1
            if( len( self.lnk ) > 0 ):
                for i,j in self.lnk:
                    self.vec[k] = 1
                    k += 1
            self.vec[self.all+1] = parm["s6"]
            self.vec[self.all+2] = parm["s8"]
            self.vec[self.all+3] = parm["a1"]
            self.vec[self.all+4] = parm["a2"]
            self.lib.qm3_dftd4_init_( ctypes.c_int( self.all ), ctypes.c_int( self.siz ), self.vec )
    
    
        def update_coor( self, mol ):
            for i in range( self.nat ):
                i3 = self.sel[i] * 3
                j3 = i * 3
                for j in [0, 1, 2]:
                    self.vec[j3+j] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) 
            if( len( self.lnk ) > 0 ):
                self.vla = []
                k = self.nat 
                for i,j in self.lnk:
                    k3 = 3 * k
                    c, v = qm3.engines.LA_coordinates( i, j, mol )
                    for m in [0,1,2]:
                        self.vec[k3+m] = c[m]
                    self.vla.append( ( self.sel.index( i ), k, v[:] ) )
                    k += 1


        def get_func( self, mol ):
            self.update_coor( mol )
            self.lib.qm3_dftd4_calc_( ctypes.c_int( self.all ), ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0] * self._ce
    
    
        def get_grad( self, mol ):
            self.update_coor( mol )
            self.lib.qm3_dftd4_calc_( ctypes.c_int( self.all ), ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0] * self._ce
            g = [ i * self._cg for i in self.vec[1:] ]
            qm3.engines.LA_gradient( self.vla, g )
            for i in range( self.nat ):
                i3 = 3 * i
                j3 = 3 * self.sel[i]
                for j in [0, 1, 2]:
                    mol.grad[j3+j] += g[i3+j]

except:
    pass

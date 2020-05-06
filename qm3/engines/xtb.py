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
# -- xtb.h -----------------------------
#typedef struct SCC_options {
#   int prlevel;
#   int parallel;
#   double acc;
#   double etemp;
#   bool grad;
#   bool restart;
#   bool ccm;
#   int maxiter;
#   char solvent[20];
#} SCC_options;
# --------------------------------------
    class xtb__sccopt( ctypes.Structure ):
        _fields_ = [ ( 'print_level', ctypes.c_int ),
            ( 'parallel', ctypes.c_int ),
            ( 'accuracy', ctypes.c_double ),
            ( 'electronic_temperature', ctypes.c_double ),
            ( 'gradient', ctypes.c_bool ),
            ( 'restart', ctypes.c_bool ),
            ( 'ccm', ctypes.c_bool ),
            ( 'max_iterations', ctypes.c_int ),
            ( 'solvent', ctypes.c_char*20 ) ]

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
            self.smb = mol.guess_symbols( self.sel )
            self.chg = chrg

            if( not mol.chrg ):
                mol.chrg = [ 0.0 for i in range( mol.natm ) ]

            self.nQM = len( self.sel ) + len( self.lnk )
            self.nMM = len( self.nbn )

            self.qm_atn = ( ctypes.c_int * self.nQM )()
            j = 0
            for i in range( len( self.sel ) ):
                self.qm_atn[j] = qm3.elements.rsymbol[self.smb[i]]
                j += 1
            for i in range( len( self.lnk ) ):
                self.qm_atn[j] = 1
                j += 1

            self.qm_crd = ( ctypes.c_double * ( 3 * self.nQM ) )()
            self.qm_grd = ( ctypes.c_double * ( 3 * self.nQM ) )()
            self.mm_chg = ( ctypes.c_double * self.nMM )()
            self.mm_atn = ( ctypes.c_int * self.nMM )()
            self.mm_gam = ( ctypes.c_double * self.nMM )()
            for i in range( len( self.nbn ) ):
                self.mm_chg[i] = mol.chrg[self.nbn[i]]
                self.mm_atn[i] = 0
                self.mm_gam[i] = 99.0
            self.mm_crd = ( ctypes.c_double * ( 3 * self.nMM ) )()
            self.mm_grd = ( ctypes.c_double * ( 3 * self.nMM ) )()

            self.lib = ctypes.CDLL( os.getenv( "QM3_LIBXTB" ) )
# -- xtb.h -----------------------------
#extern int
#GFN2_QMMM_calculation(
#      const int* natoms,
#      const int* attyp,
#      const double* charge,
#      const int* uhf,
#      const double* coord,
#      const SCC_options* opt,
#      const char* output,
#      const int* npc,
#      const double* pc_q,
#      const int* pc_at,
#      const double* pc_gam,
#      const double* pc_coord,
#      double* energy,
#      double* grad,
#      double* pc_grad);
# --------------------------------------
            self.lib.GFN2_QMMM_calculation.argtypes = [ 
                ctypes.POINTER( ctypes.c_int ),     # number of atoms
                ctypes.POINTER( ctypes.c_int ),     # atomic numbers
                ctypes.POINTER( ctypes.c_double ),  # molecular charge
                ctypes.POINTER( ctypes.c_int ),     # number of unpaired electrons
                ctypes.POINTER( ctypes.c_double ),  # coordinates
                ctypes.POINTER( xtb__sccopt ),      # SCC Options
                ctypes.c_char_p,                    # output file name
                ctypes.POINTER( ctypes.c_int ),     # number of MM point charges
                ctypes.POINTER( ctypes.c_double ),  # MM charges
                ctypes.POINTER( ctypes.c_int ),     # NULL
                ctypes.POINTER( ctypes.c_double ),  # NULL
                ctypes.POINTER( ctypes.c_double ),  # MM coordinates
                ctypes.POINTER( ctypes.c_double ),  # energy
                ctypes.POINTER( ctypes.c_double ),  # gradient
                ctypes.POINTER( ctypes.c_double ) ] # MM gradient
            self.lib.GFN2_QMMM_calculation.restype = ctypes.c_int

    
        def update_coor( self, mol ):
            for i in range( len( self.sel ) ):
                i3 = self.sel[i] * 3
                j3 = i * 3
                for j in [0, 1, 2]:
                    self.qm_crd[j3+j] = ( mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) ) / self._cx
            self.vla = []
            k = len( self.sel )
            for i in range( len( self.lnk ) ):
                j3 = k * 3
                c, v = qm3.engines.LA_coordinates( self.lnk[i][0], self.lnk[i][1], mol )
                for j in [0, 1, 2]:
                    self.qm_crd[j3+j] = c[j] / self._cx
                self.vla.append( ( self.sel.index( self.lnk[i][0] ), k, v[:] ) )
                k += 1
            for i in range( len( self.nbn ) ):
                i3 = self.nbn[i] * 3
                j3 = i * 3
                for j in [0, 1, 2]:
                    self.mm_crd[j3+j] = ( mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) ) / self._cx


        def get_func( self, mol ):
            self.update_coor( mol )
            chg = ctypes.c_double( self.chg )
            mul = ctypes.c_int( 0 )
            opt = xtb__sccopt( 0, 0, 1.0, 300.0, False, False, False, 1000, "none".encode( "utf-8" ) )
            ene = ctypes.c_double( 0.0 )
            ret = self.lib.GFN2_QMMM_calculation(
                ctypes.c_int( self.nQM ),
                self.qm_atn,
                chg,
                mul,
                self.qm_crd,
                opt,
                "-".encode( "utf-8" ),
                ctypes.c_int( self.nMM ),
                self.mm_chg,
                self.mm_atn,
                self.mm_gam,
                self.mm_crd,
                ene,
                self.qm_grd,
                self.mm_grd )
            if( ret == 0 ):
                mol.func += ene.value * self._ce
            else:
                sys.exit( -1 )
    
    
        def get_grad( self, mol ):
            self.update_coor( mol )
            chg = ctypes.c_double( self.chg )
            mul = ctypes.c_int( 0 )
            opt = xtb__sccopt( 0, 0, 1.0, 300.0, True, False, False, 1000, "none".encode( "utf-8" ) )
            ene = ctypes.c_double( 0.0 )
            ret = self.lib.GFN2_QMMM_calculation(
                ctypes.c_int( self.nQM ),
                self.qm_atn,
                chg,
                mul,
                self.qm_crd,
                opt,
                "-".encode( "utf-8" ),
                ctypes.c_int( self.nMM ),
                self.mm_chg,
                self.mm_atn,
                self.mm_gam,
                self.mm_crd,
                ene,
                self.qm_grd,
                self.mm_grd )
            if( ret == 0 ):
                mol.func += ene.value * self._ce
                g = [ j * self._cg for j in self.qm_grd ]
                qm3.engines.LA_gradient( self.vla, g )
                for i in range( len( self.sel ) ):
                    i3 = self.sel[i] * 3
                    j3 = i * 3
                    for j in [0, 1, 2]:
                        mol.grad[i3+j] += g[j3+j]
                for i in range( len( self.nbn ) ):
                    i3 = self.nbn[i] * 3
                    j3 = i * 3
                    for j in [0, 1, 2]:
                        mol.grad[i3+j] += self.mm_grd[j3+j] * self._cg
            else:
                sys.exit( -1 )

except:
    pass

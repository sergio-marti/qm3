# -*- coding: iso-8859-1 -*-
import os
import qm3.engines



def dftb_input( obj, mol, run ):
    s_qm = "  %d C\n  %s\n"%( len( obj.sel ) + len( obj.lnk ), str.join( " ", obj.tbl ) )
    j = 0
    for i in obj.sel:
        i3 = i * 3
        s_qm += "  %4d%4d%20.10lf%20.10lf%20.10lf\n"%( j + 1, obj.tbl.index( obj.smb[j] ) + 1,
            mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
            mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
            mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) )
        j += 1
    if( len( obj.lnk ) > 0 ):
        obj.vla = []
        k = len( obj.sel )
        w = obj.tbl.index( "H" ) + 1
        for i,j in obj.lnk:
            c, v = qm3.engines.LA_coordinates( i, j, mol )
            s_qm += "  %4d%4d%20.10lf%20.10lf%20.10lf\n"%( k + 1, w, c[0], c[1], c[2] )
            obj.vla.append( ( obj.sel.index( i ), k, [-v[0], -v[1], -v[2]] ) )
            k += 1
    s_wf = ""
    if( os.access( "charges.bin", os.R_OK ) ):
        s_wf = "  ReadInitialCharges = Yes"
    s_rn = "  CalculateForces = No"
    if( run == "grad" ):
        s_rn = "  CalculateForces = Yes"
    s_nq = ""
    if( len( obj.nbn ) > 0 ):
        s_nq = str( len( obj.nbn ) )
        g = open( "charges.dat", "wt" )
        for i in obj.nbn:
            i3 = i * 3
            g.write( "%20.10lf%20.10lf%20.10lf%12.6lf\n"%(
                mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), mol.chrg[i] ) )
        g.close()
    f = open( "dftb_in.hsd", "wt" )
    buf = obj.inp.replace( "qm3_atoms", s_qm[:-1] )
    buf = buf.replace( "qm3_guess", s_wf )
    buf = buf.replace( "qm3_job", s_rn )
    buf = buf.replace( "qm3_nchg", s_nq )
    f.write( buf )
    f.close()



class run_single( qm3.engines.qmbase ):

    def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
        qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
        self.exe = "bash r.dftb"
        self.tbl = { i:None for i in self.smb }
        if( len( self.lnk ) > 0 ):
            self.tbl["H"] = None
        self.tbl = list( self.tbl )


    def mk_input( self, mol, run ):
        dftb_input( self, mol, run )


    def parse_log( self, mol, run ):
        f = open( "detailed.out", "rt" )
        l = f.readline()
        while( l != "" ):
            if( l.strip() == "Net atomic charges (e)" or l.strip() == "Atomic gross charges (e)" ):
                f.readline()
                for i in range( len( self.sel ) ):
                    mol.chrg[self.sel[i]] = float( f.readline().split()[1] )
            if( l[0:20].strip() == "Total energy:" ):
                mol.func += float( l.split()[2] ) * self._ce
            if( l.strip() == "Total Forces" and run == "grad" ):
                g = []
                for i in range( len( self.sel ) + len( self.vla ) ):
                    g += [ - float( j ) * self._cg for j in f.readline().split()[1:] ]
                qm3.engines.LA_gradient( self.vla, g )
                for i in range( len( self.sel ) ):
                    i3 = i * 3
                    for j in [0, 1, 2]:
                        mol.grad[3*self.sel[i]+j] += g[i3+j]
            if( l.strip() == "Forces on external charges" and run == "grad" ):
                for i in range( len( self.nbn ) ):
                    g = [ - float( j ) * self._cg for j in f.readline().split() ]
                    for j in [0, 1, 2]:
                        mol.grad[3*self.nbn[i]+j] += g[j]
            l = f.readline()
        f.close()



try:
    import ctypes

    class run_dynlib( qm3.engines.qmbase ):
    
        def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
            qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
            self.exe = "bash r.dftb"
            self.tbl = { i:None for i in self.smb }
            if( len( self.lnk ) > 0 ):
                self.tbl["H"] = None
            self.tbl = list( self.tbl )

            self.nQM = len( self.sel ) + len( self.lnk )
            self.nMM = len( self.nbn )
            self.siz = 1 + 3 * ( self.nQM + self.nMM ) + self.nMM + self.nQM
            self.vec = ( ctypes.c_double * self.siz )()
            self.lib = ctypes.CDLL( os.getenv( "QM3_LIBDFTB" ) )
            self.lib.qm3_dftb_calc_.argtypes = [ ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_double ) ]
            self.lib.qm3_dftb_calc_.restype = None

            dftb_input( self, mol, "grad" )
            self.lib.qm3_dftb_init_()
    
    
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
            k = 3 * ( self.nQM + self.nMM )
            for i in range( self.nMM ):
                i3 = self.nbn[i] * 3
                j3 = ( self.nQM + i ) * 3
                for j in [0, 1, 2]:
                    self.vec[j3+j] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) 
#            for i in range( self.nMM ):
                self.vec[k+i] = mol.chrg[self.nbn[i]]


        def get_func( self, mol ):
            self.update_coor( mol )
            self.lib.qm3_dftb_calc_( ctypes.c_int( self.nQM ), ctypes.c_int( self.nMM ), ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0] * self._ce
            for i in range( len( self.sel ) ):
                mol.chrg[self.sel[i]] = self.vec[i+1]
    
    
        def get_grad( self, mol ):
            self.update_coor( mol )
            self.lib.qm3_dftb_calc_( ctypes.c_int( self.nQM ), ctypes.c_int( self.nMM ), ctypes.c_int( self.siz ), self.vec )
            mol.func += self.vec[0] * self._ce
            for i in range( len( self.sel ) ):
                mol.chrg[self.sel[i]] = self.vec[i+1]
            g = [ j * self._cg for j in self.vec[self.nQM+1:] ]
            qm3.engines.LA_gradient( self.vla, g )
            for i in range( len( self.sel ) ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    mol.grad[3*self.sel[i]+j] += g[i3+j]
            for i in range( self.nMM ):
                i3 = ( self.nQM + i ) * 3
                for j in [0, 1, 2]:
                    mol.grad[3*self.nbn[i]+j] += g[i3+j]

except:
    pass

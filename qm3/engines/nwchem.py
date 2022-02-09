# -*- coding: iso-8859-1 -*-
import os
import math
import qm3.engines



class run_single( qm3.engines.qmbase ):

    def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
        qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
        self.exe = "bash r.nwchem"


    def mk_input( self, mol, run ):
        s_qm = ""
        j = 0
        for i in self.sel:
            i3 = i * 3
            s_qm += "%4s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j],
                mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) )
            j += 1
        if( len( self.lnk ) > 0 ):
            self.vla = []
            k = len( self.sel )
            for i,j in self.lnk:
                c, v = qm3.engines.LA_coordinates( i, j, mol )
                s_qm += "%-4s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] )
                self.vla.append( ( self.sel.index( i ), k, v[:] ) )
                k += 1
        s_wf = ""
        if( os.access( "nwchem.movecs", os.R_OK ) ):
            s_wf = "vectors input nwchem.movecs"
        s_mm = ""
        if( len( self.nbn ) > 0 ):
            s_mm = "set bq:max_nbq %d\n"%( len( self.nbn ) + 1 )
            s_mm += "bq units angstroms\n  force nwchem.mmgrad\n  load nwchem.mmchrg units angstroms format 1 2 3 4\nend"
            g = open( "nwchem.mmchrg", "wt" )
            for i in self.nbn:
                i3 = i * 3
                g.write( "%20.10lf%20.10lf%20.10lf%12.4lf\n"%(
                    mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                    mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                    mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), mol.chrg[i] ) )
            g.close()
        if( self.inp.lower().find( "dft" ) > -1 ):
            s_rn = "task dft"
        else:
            s_rn = "task scf"
        if( run == "grad" ):
            s_rn += " gradient"
        f = open( "nwchem.nw", "wt" )
        buf = self.inp.replace( "qm3_atoms", s_qm[:-1] )
        buf = buf.replace( "qm3_job", s_rn )
        buf = buf.replace( "qm3_guess", s_wf )
        buf = buf.replace( "qm3_charges", s_mm )
        f.write( buf )
        f.close()


    def parse_log( self, mol, run ):
        f = open( "nwchem.log", "rt" )
        l = f.readline()
        while( l != "" ):
            L = l.strip()
            # read energy
            if( L.find( "Total " ) > -1 and L.find( " energy " ) > -1 ):
                t = L.split()
                if( len( t ) == 5 ):
                    mol.func += float( t[4] ) * self._ce
            # read gradient and LAs projection
            if( run == "grad" and L.find( "ENERGY GRADIENTS" ) > -1 ):
                f.readline(); f.readline(); f.readline()
                g = []
                for i in range( len( self.sel ) + len( self.lnk ) ):
                    g += [ float( j ) * self._cg for j in f.readline().strip().split()[-3:] ]
                qm3.engines.LA_gradient( self.vla, g )
                for i in range( len( self.sel ) ):
                    i3 = i * 3
                    for j in [0, 1, 2]:
                        mol.grad[3*self.sel[i]+j] += g[i3+j]
            l = f.readline()
        f.close()
        if( self.nbn and os.access( "nwchem.mmgrad", os.R_OK ) ):
            f = open( "nwchem.mmgrad", "rt" )
            f.readline()
            for i in range( len( self.nbn ) ):
                t = [ float( j ) for j in f.readline().strip().split() ]
                mol.grad[3*self.nbn[i]]   += t[0] * self._cg
                mol.grad[3*self.nbn[i]+1] += t[1] * self._cg
                mol.grad[3*self.nbn[i]+2] += t[2] * self._cg
            f.close()


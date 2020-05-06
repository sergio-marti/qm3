# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import os
import re
import qm3.engines



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
                f.write( "%12.4lf%20.10lf%20.10lf%20.10lf\n"%( mol.chrg[i],
                    mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3] / mol.boxl[0], 0 ), 
                    mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                    mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
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


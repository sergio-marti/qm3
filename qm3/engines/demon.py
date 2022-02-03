# -*- coding: iso-8859-1 -*-
import os
import qm3.engines



class demon( qm3.engines.qmbase ):

    def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
        qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
        self.exe = "bash r.demon"


    def mk_input( self, mol, run ):
        s_wf = ""
        if( os.access( "deMon.rst", os.R_OK ) ):
            s_wf = "guess restart"
        s_qm = ""
        j = 0
        for i in self.sel:
            i3 = i * 3
            s_qm += "%4s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j],
                mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) )
            j += 1
        if( len( self.lnk ) > 0):
            self.vla = []
            k = len( self.sel )
            for i,j in self.lnk:
                c, v = qm3.engines.LA_coordinates( i, j, mol )
                s_qm += "%-4s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] )
                self.vla.append( ( self.sel.index( i ), k, v[:] ) )
                k += 1
        f = open( "deMon.inp", "wt" )
        buf = self.inp.replace( "qm3_atoms", s_qm )
        buf = buf.replace( "qm3_guess", s_wf )
        f.write( buf )
        f.close()
        if( len( self.nbn ) > 0 ):
            f = open( "deMon.cub", "wt" )
            for i in self.nbn:
                i3 = i * 3
                f.write( "%20.10lf%20.10lf%20.10lf%12.4lf\n"%(
                    mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                    mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                    mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), mol.chrg[i] ) )
            f.close()


    def parse_log( self, mol, run ):
        f = open( "deMon.qmm", "rt" )
        mol.func += float( f.readline().split()[-1] ) * self._ce
        if( run == "grad" ):
            l = f.readline()
            while( l != "" ):
                # read gradient and LAs projection
                if( l.strip() == "QMFORCES" ):
                    g = []
                    for i in range( len( self.sel ) + len( self.lnk ) ):
                        g += [ float( j ) * self._cg for j in f.readline().strip().split()[-3:] ]
                    qm3.engines.LA_gradient( self.vla, g )
                    for i in range( len( self.sel ) ):
                        i3 = i * 3
                        for j in [0, 1, 2]:
                            mol.grad[3*self.sel[i]+j] += g[i3+j]
                # read MM gradient
                if( self.nbn and l.strip() == "EMBEDFORCES" ):
                    for i in range( len( self.nbn ) ):
                        t = [ float( j ) * self._cg for j in f.readline().strip().split() ]
                        mol.grad[3*self.nbn[i]]   += t[0]
                        mol.grad[3*self.nbn[i]+1] += t[1]
                        mol.grad[3*self.nbn[i]+2] += t[2]
                l = f.readline()
        f.close()

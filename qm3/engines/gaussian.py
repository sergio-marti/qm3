# -*- coding: iso-8859-1 -*-
import os
import qm3.maths.matrix
import qm3.engines



#
# detects whether return the electrostatic charges gradient (_mmint.QMLJ) searching for "prop=(field,read)" 
# or it will be calcualted after classically (_mmint.QMLJ_MMEL)
#
class gaussian( qm3.engines.qmbase ):

    def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
        qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
        self.exe = "bash r.gauss"
        self.gmm = ( self.inp.lower().find( "prop=(field,read)" ) > -1 )


    def mk_input( self, mol, run ):
        s_qm = ""
        j = 0
        for i in self.sel:
            i3 = i * 3
            s_qm += "%2s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j],
                    mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                    mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                    mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) )
            j += 1
        if( len( self.lnk ) > 0 ):
            self.vla = []
            k = len( self.sel )
            for i,j in self.lnk:
                c, v = qm3.engines.LA_coordinates( i, j, mol )
                s_qm += "%-2s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] )
                self.vla.append( ( self.sel.index( i ), k, v[:] ) )
                k += 1
        s_mm = ""
        s_nq = ""
        if( len( self.nbn ) > 0 ):
            for i in self.nbn:
                i3 = i * 3
                s_mm += "%20.10lf%20.10lf%20.10lf%12.4lf\n"%( 
                    mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                    mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                    mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), 
                    mol.chrg[i] )
                s_nq += "%20.10lf%20.10lf%20.10lf\n"%( 
                    mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                    mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                    mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) )
        s_rn = ""
        if( run == "grad" ):
            s_rn = "force"
        elif( run == "hess" ):
            s_rn = "freq=noraman cphf(maxinv=10000)"
        s_wf = ""
        if( os.access( "gauss.chk", os.R_OK ) ):
            s_wf = "guess=(read)"
        f = open( "gauss.com", "wt" )
        buf = self.inp.replace( "qm3_atoms", s_qm[:-1] )
        buf = buf.replace( "qm3_job", s_rn )
        buf = buf.replace( "qm3_guess", s_wf )
        buf = buf.replace( "qm3_charges", s_mm[:-1] )
        buf = buf.replace( "qm3_field", s_nq[:-1] )
        f.write( buf )
        f.write( "\n\n\n\n\n" )
        f.close()


    def parse_log( self, mol, run ):
        fd = open( "Test.FChk", "rt" )
        l = fd.readline()
        while( l != "" ):
            # read energy
            if( l[0:12] == "Total Energy" ):
                mol.func += float( l.split()[3] ) * self._ce
            # read gadient
            if( run in [ "grad", "hess" ] and l[0:18] == "Cartesian Gradient" ):
                i = int( l.split()[-1] )
                j = int( i // 5 ) + ( i%5 != 0 )
                i = 0
                g = []
                while( i < j ):
                    l = fd.readline()
                    for itm in l.split():
                        g.append( float( itm ) * self._cg )
                    i += 1
                # remove LAs from gradient
                qm3.engines.LA_gradient( self.vla, g )
                # copy array
                for i in range( len( self.sel ) ):
                    i3 = i * 3
                    for j in [0, 1, 2]:
                        mol.grad[3*self.sel[i]+j] += g[i3+j]
                # read hessian (columns)
                if( run == "hess" ):
                    l = fd.readline()
                    i = int( l.split()[-1] )
                    j = int( i // 5 ) + ( i % 5 != 0 )
                    i = 0
                    h = []
                    while( i < j ):
                        l = fd.readline()
                        for itm in l.split():
                            h.append( float( itm ) * self._ch )
                        i += 1
                    # truncate LAs and swap hessian (cols>>rows)
                    i = 3 * len( self.sel )
                    j = i * ( i + 1 ) // 2
##                    mol.hess = qm3.maths.matrix.from_upper_diagonal_columns( h[0:j], i )
                    t = qm3.maths.matrix.from_upper_diagonal_columns( h[0:j], i )
                    for j in range( i * i ):
                        mol.hess[j] += t[j]
            # read charges
            if( l[0:11] == "ESP Charges" ):
                i = int( l.split()[-1] )
                j = int( i // 5 ) + ( i % 5 != 0 )
                i = 0
                k = 0
                while( i < j ):
                    l = fd.readline()
                    for itm in l.split():
                        if( k < len( self.sel ) ):
                            mol.chrg[self.sel[k]] = float( itm )
                            k += 1
                    i += 1
            l = fd.readline()
        fd.close()
        # remove autoenergy of the charges, and calculate MM gradient
        if( len( self.nbn ) > 0 ):
            fd = open( "gauss.log", "rt" )
            l = fd.readline()
            fl = True
            while( l != "" and fl ):
                if( l[0:29] == " Self energy of the charges =" ):
                    fl = False
                    mol.func -= float( l.split()[-2] ) * self._ce
                l = fd.readline()
            if( run in [ "grad", "hess" ] and self.gmm ):
                fl = True
                while( l != "" and fl ):
                    if( l.strip() == "Potential          X             Y             Z" ):
                        fl = False
                        for i in range( 1 + len( self.sel ) + len( self.lnk ) ):
                            fd.readline()
                        for i in self.nbn:
                            i3 = i * 3
                            t = fd.readline().split()[2:]
                            for j in [0, 1, 2]:
                                mol.grad[i3+j] += - self._cg * mol.chrg[i] * float( t[j] )
                    l = fd.readline()
            fd.close()
        # return
        os.unlink( "Test.FChk" )



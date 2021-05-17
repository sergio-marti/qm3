# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.constants
import qm3.utils
import qm3.fio
import qm3.maths.matrix


# 10.1103/PhysRevLett.100.020603


class distance( object ):
    def __init__( self, indx, vpot, sigm, nstp, wall, name = None ):
        """
        vpot [kJ/mol]: 2 k_B T [ 2.5 kJ/mol @ 300 K ]
        sigm [A]: rms of the variable (md run in the minimia)
        nstp [md steps]: 50 - 100
        wall [kumb, llim, ulim]: kumb [kJ/mol], llim + sigm/2 [A], ulim - sigm/2 [A]
        """
        self.vpot = vpot
        self.sigm = sigm
        self.indx = indx[:]
        self.wall = wall[:]
        self.nstp = nstp
        self.accu = 0.0
        self.step = 0
        self.valu = []
        self.fd   = None
        self.name = name
        if( name != None ):
            self.fd = qm3.fio.open_w( name )
            self.fd.write( "%20.10lf%20.10lf\n"%( vpot, sigm ) )
        else:
            self.fd = None


    def get_grad( self, molec ):
        # calculate current coordinate [cv] and its derivative [dr]
        ai = 3 * self.indx[0]
        aj = 3 * self.indx[1]
        dr = [ i-j for i,j in zip( molec.coor[ai:ai+3], molec.coor[aj:aj+3] ) ]
        cv = math.sqrt( sum( [ i * i for i in dr ] ) )
        dr = [ i / cv for i in dr ]
        # increment the step number and accumulate
        if( cv >= self.wall[1] and cv <= self.wall[2] ):
            self.step += 1
            self.accu += cv
            if( self.step % self.nstp == 0 ):
                nv = self.accu / self.nstp
                self.valu.append( nv )
                if( self.fd != None ):
                    self.fd.write( "%20.10lf\n"%( nv ) )
                    self.fd.flush()
                self.accu = 0.0
                self.step = 0
        # add gaussians
        df = [ .0, .0, .0 ]
        for rv in self.valu:
            ds = ( cv - rv ) / self.sigm
            ss = self.vpot * math.exp( - 0.5 * ds * ds )
            molec.func += ss
            for k in [0, 1, 2]:
                df[k] -= ss * ds * dr[k] / self.sigm
        # add walls
        if( cv < self.wall[1] ):
            ds = self.wall[0] * ( cv - self.wall[1] - self.sigm * 0.5 )
            molec.func += 0.5 * ds * ( cv - self.wall[1] - self.sigm * 0.5 )
            for k in [0, 1, 2]:
                df[k] += ds * dr[k]
        if( cv > self.wall[2] ):
            ds = self.wall[0] * ( cv - self.wall[2] + self.sigm * 0.5 )
            molec.func += 0.5 * ds * ( cv - self.wall[2] + self.sigm * 0.5 )
            for k in [0, 1, 2]:
                df[k] += ds * dr[k]
        # update gradients
        for k in [0, 1, 2]:
            molec.grad[ai+k] += df[k]
            molec.grad[aj+k] -= df[k]


    def close( self ):
        qm3.fio.close( self.fd, self.name )



class multiple_distance( object ):
    def __init__( self, indx, weig, vpot, sigm, nstp, wall, name = None ):
        """
        vpot [kJ/mol]: 2 k_B T [ 2.5 kJ/mol @ 300 K ]
        sigm [A]: rms of the variable (md run in the minimia)
        nstp [md steps]: 50 - 100
        wall [kumb, llim, ulim]: kumb [kJ/mol], llim + sigm/2 [A], ulim - sigm/2 [A]
        """
        if( len( weig ) * 2 != len( indx ) ):
            print( "- metadyn.multiple_distance: Number of ATOMS should be TWICE the number of WEIGHTS!" )
        self.vpot = vpot
        self.sigm = sigm
        self.indx = indx[:]
        self.weig = weig[:]
        self.size = len( weig )
        self.wall = wall[:]
        self.nstp = nstp
        self.accu = 0.0
        self.step = 0
        self.valu = []
        self.fd   = None
        self.name = name
        if( name != None ):
            self.fd = qm3.fio.open_w( name )
            self.fd.write( "%20.10lf%20.10lf\n"%( vpot, sigm ) )
        else:
            self.fd = None


    def get_grad( self, molec ):
        # calculate current coordinate [cv] and its derivative [dr]
        dr = []
        cv = 0.0
        for i in range( self.size ):
            i3 = i * 3
            ii = 3 * self.indx[2*i]
            jj = 3 * self.indx[2*i+1]
            dr += [ j-k for j,k in zip( molec.coor[ii:ii+3], molec.coor[jj:jj+3] ) ]
            rr = math.sqrt( sum( [ j * j for j in dr[i3:i3+3] ] ) )
            cv += rr * self.weig[i]
            for j in range( i3, i3 + 3 ):
                dr[j] /= rr
        # increment the step number and accumulate
        if( cv >= self.wall[1] and cv <= self.wall[2] ):
            self.step += 1
            self.accu += cv
            if( self.step % self.nstp == 0 ):
                nv = self.accu / self.nstp
                self.valu.append( nv )
                if( self.fd != None ):
                    self.fd.write( "%20.10lf\n"%( nv ) )
                    self.fd.flush()
                self.accu = 0.0
                self.step = 0
        # add gaussians
        df = [ .0 for i in range( 3 * self.size ) ]
        for rv in self.valu:
            ds = ( cv - rv ) / self.sigm
            ss = self.vpot * math.exp( - 0.5 * ds * ds )
            molec.func += ss
            for i in range( self.size ):
                i3 = i * 3
                for k in [0, 1, 2]:
                    df[i3+k] -= ss * ds * dr[i3+k] / self.sigm * self.weig[i]
        # add walls
        if( cv < self.wall[1] ):
            ds = self.wall[0] * ( cv - self.wall[1] - self.sigm * 0.5 )
            molec.func += 0.5 * ds * ( cv - self.wall[1] - self.sigm * 0.5 )
            for i in range( self.size ):
                i3 = i * 3
                for k in [0, 1, 2]:
                    df[i3+k] += ds * dr[i3+k] * self.weig[i]
        if( cv > self.wall[2] ):
            ds = self.wall[0] * ( cv - self.wall[2] + self.sigm * 0.5 )
            molec.func += 0.5 * ds * ( cv - self.wall[2] + self.sigm * 0.5 )
            for i in range( self.size ):
                i3 = i * 3
                for k in [0, 1, 2]:
                    df[i3+k] += ds * dr[i3+k] * self.weig[i]
        # update gradients
        for i in range( self.size ):
            i3 = i * 3
            ii = 3 * self.indx[2*i]
            jj = 3 * self.indx[2*i+1]
            for k in [0, 1, 2]:
                molec.grad[ii+k] += df[i3+k]
                molec.grad[jj+k] -= df[i3+k]


    def close( self ):
        qm3.fio.close( self.fd, self.name )



class colvar_s( object ):
    def __init__( self, conf, str_crd, str_met, vpot, sigm, nstp, wall, name = None ):
        """
        vpot [kJ/mol]: 2 k_B T [ 2.5 kJ/mol @ 300 K ]
        sigm [A]: rms of the variable (md run in the minimia)
        nstp [md steps]: 50 - 100
        wall [kumb, llim, ulim]: kumb [kJ/mol], llim + sigm/2 [A], ulim - sigm/2 [A]
------------------------------------------------------------------------
ncrd      nwin
dist      atom_i    atom_j
...
dist      atom_i    atom_j
------------------------------------------------------------------------
"""
        f = qm3.fio.open_r( conf )
        t = f.readline().strip().split()
        self.ncrd = int( t[0] )
        self.nwin = int( t[1] )
        self.jidx = {}
        self.func = []
        self.atom = []
        for i in range( self.ncrd ):
            t = f.readline().strip().split()
            if( t[0][0:4] == "dist" and len( t ) == 3 ):
                self.func.append( self.distance )
                a_i = int( t[1] )
                a_j = int( t[2] )
                self.atom.append( ( a_i, a_j ) )
                self.jidx[a_i] = True
                self.jidx[a_j] = True
        qm3.fio.close( f, conf )
        self.jidx = { jj: ii for ii,jj in zip( range( len( self.jidx ) ), sorted( self.jidx ) ) }
        self.idxj = { self.jidx[ii]: ii for ii in iter( self.jidx ) }
        self.jcol = 3 * len( self.jidx )
        # load (previous) equi-destributed string
        f = qm3.fio.open_r( str_crd )
        self.rcrd = [ float( i ) for i in f.read().split() ]
        qm3.fio.close( f, str_crd )
        # load (previous) string metrics
        nc2 = self.ncrd * self.ncrd
        if( str_met != None ):
            f = qm3.fio.open_r( str_met )
            self.rmet = [ float( i ) for i in f.read().split() ]
            qm3.fio.close( f, str_met )
        # use identity for the metrics
        else:
            tmp = [ 0.0 for i in range( nc2 ) ]
            for i in range( self.ncrd ):
                tmp[i*self.ncrd+i] = 1.0
            self.rmet = []
            for i in range( self.nwin ):
                self.rmet += tmp[:]
        # get the arc length of the current string...
        self.arcl = []
        for i in range( 1, self.nwin ):
            tmp = [ self.rcrd[i*self.ncrd+j] - self.rcrd[(i-1)*self.ncrd+j] for j in range( self.ncrd ) ]
            mat = qm3.maths.matrix.inverse( [ 0.5 * ( self.rmet[i*nc2+j] + self.rmet[(i-1)*nc2+j] ) for j in range( nc2 ) ], self.ncrd, self.ncrd )
            mat = qm3.maths.matrix.mult( mat, self.ncrd, self.ncrd, tmp, self.ncrd, 1 )
            self.arcl.append( math.sqrt( sum( [ tmp[j] * mat[j] for j in range( self.ncrd ) ] ) ) )
        self.delz = sum( self.arcl ) / float( self.nwin - 1.0 )
        print( "Colective variable s range: [%.3lf - %.3lf: %.6lf] _Ang"%( 0.0, sum( self.arcl ), self.delz ) )
        # store inverse metrics from references...
        tmp = []
        for i in range( self.nwin ):
            tmp += qm3.maths.matrix.inverse( self.rmet[i*nc2:(i+1)*nc2], self.ncrd, self.ncrd )
        self.rmet = tmp[:]
        # metadynamics stuff
        self.vpot = vpot
        self.sigm = sigm
        self.wall = wall[:]
        self.nstp = nstp
        self.accu = 0.0
        self.step = 0
        self.valu = []
        self.fd   = None
        self.name = name
        if( name != None ):
            self.fd = qm3.fio.open_w( name )
            self.fd.write( "%20.10lf%20.10lf\n"%( vpot, sigm ) )
        else:
            self.fd = None


    def get_grad( self, molec ):
        ccrd = []
        jaco = [ 0.0 for i in range( self.ncrd * self.jcol ) ]
        for i in range( self.ncrd ):
            ccrd.append( self.func[i]( i, molec, jaco ) )
        nc2  = self.ncrd * self.ncrd
        cdst = []
        jder = []
        for i in range( self.nwin ):
            vec = [ ccrd[j] - self.rcrd[i*self.ncrd+j] for j in range( self.ncrd ) ]
            mat = qm3.maths.matrix.mult( self.rmet[i*nc2:(i+1)*nc2], self.ncrd, self.ncrd, vec, self.ncrd, 1 )
            cdst.append( math.sqrt( sum( [ vec[j] * mat[j] for j in range( self.ncrd ) ] ) ) )
            sd1 = qm3.maths.matrix.mult( mat, 1, self.ncrd, jaco, self.ncrd, self.jcol )
            sd2 = qm3.maths.matrix.mult( self.rmet[i*nc2:(i+1)*nc2], self.ncrd, self.ncrd, jaco, self.ncrd, self.jcol )
            sd3 = qm3.maths.matrix.mult( vec, 1, self.ncrd, sd2, self.ncrd, self.jcol )
            jder.append( [ 0.5 * ( sd1[j] + sd3[j] ) / cdst[-1] for j in range( self.jcol ) ] )
        cexp = [ math.exp( - cdst[i] / self.delz ) for i in range( self.nwin ) ]
        sumn = sum( [ i * self.delz * cexp[i] for i in range( self.nwin ) ] )
        sumd = sum( cexp )
        cval = sumn / sumd
        # increment the step number and accumulate
        if( cval >= self.wall[1] and cval <= self.wall[2] ):
            self.step += 1
            self.accu += cval
            if( self.step % self.nstp == 0 ):
                nv = self.accu / self.nstp
                self.valu.append( nv )
                if( self.fd != None ):
                    self.fd.write( "%20.10lf\n"%( nv ) )
                    self.fd.flush()
                self.accu = 0.0
                self.step = 0
        # add gaussians
        sder = [ 0.0 for i in range( self.jcol ) ]
        for rv in self.valu:
            ds = ( cval - rv ) / self.sigm
            ss = self.vpot * math.exp( - 0.5 * ds * ds )
            molec.func += ss
            for i in range( self.jcol ):
                for j in range( self.nwin ):
                    sder[i] -= ss * ds / self.sigm * jder[j][i] * ( cval / self.delz - j ) * cexp[j] / sumd
        # add walls
        if( cval < self.wall[1] ):
            ds = self.wall[0] * ( cval - self.wall[1] - self.sigm * 0.5 )
            molec.func += 0.5 * ds * ( cval - self.wall[1] - self.sigm * 0.5 )
            for i in range( self.jcol ):
                for j in range( self.nwin ):
                    sder[i] += ds * jder[j][i] * ( cval / self.delz - j ) * cexp[j] / sumd
        if( cval > self.wall[2] ):
            ds = self.wall[0] * ( cval - self.wall[2] + self.sigm * 0.5 )
            molec.func += 0.5 * ds * ( cval - self.wall[2] + self.sigm * 0.5 )
            for i in range( self.jcol ):
                for j in range( self.nwin ):
                    sder[i] += ds * jder[j][i] * ( cval / self.delz - j ) * cexp[j] / sumd
        # update gradients
        for i in range( len( self.jidx ) ):
            i3 = i * 3
            j3 = self.idxj[i] * 3
            for j in [0, 1, 2]:
                molec.grad[j3+j] += sder[i3+j]
        return( cval, ccrd )


    def distance( self, icrd, molec, jacob ):
        ai = self.atom[icrd][0]
        aj = self.atom[icrd][1]
        dd = [ (jj-ii) for ii,jj in zip( molec.coor[3*ai:3*ai+3], molec.coor[3*aj:3*aj+3] ) ]
        vv = math.sqrt( sum( [ ii*ii for ii in dd ] ) )
        for k in [0, 1, 2]:
            jacob[icrd*self.jcol+3*self.jidx[ai]+k] -= dd[k] / vv
            jacob[icrd*self.jcol+3*self.jidx[aj]+k] += dd[k] / vv
        return( vv )


    def close( self ):
        qm3.fio.close( self.fd, self.name )

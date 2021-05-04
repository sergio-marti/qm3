# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.constants
import qm3.utils
import qm3.maths.matrix


#10.1103/PhysRevLett.100.020603


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
        if( name != None ):
            self.fd = open( name, "wt" )
            self.fd.write( "%20.10lf%20.10lf\n"%( vpot, sigm ) )


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
        if( self.fd != None ):
            self.fd.close()



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
        if( name != None ):
            self.fd = open( name, "wt" )
            self.fd.write( "%20.10lf%20.10lf\n"%( vpot, sigm ) )


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
        if( self.fd != None ):
            self.fd.close()

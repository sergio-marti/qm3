# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.maths.rand
import qm3.utils
import qm3.constants
import time



def default_log( txt ):
    sys.stdout.write( txt + "\n" )
    sys.stdout.flush()



class model:
    def __init__( self, molec ):
        """
        define as much engines as needed based on the molec
        """
        pass

    def get_grad( self, molec ):
        """
        sequencially accumulate all the engines.get_grad
        """
        pass



class md_template:
    def __init__( self ):
        """
        setup a molecule (with masses):     self.mole
        setup a working model:              self.engn
        setup the active atoms:             self.sele
        setup the problem size:             self.size = 3 * len( self.sele )
        map the active coordinates on:      self.coor
        map the active masses on:           self.mass
        """
        self.mole = None
        self.engn = None
        self.sele = None
        self.size = 0.0
        self.coor = []
        self.mass = []
        self.func = 0.0
        self.grad = []
        self.velo = []

        self.get_grad = self.get_grad_aver



    def current_step( self, istep ):
        pass



    def setup( self, pi_atoms, num_beads = 8, temperature = 300.0 ):
        self.temp = temperature
        self.__kt = 1000.0 * qm3.constants.NA * qm3.constants.KB * self.temp
        self.__cc = 10.0 / ( qm3.constants.KB * qm3.constants.NA )
        # -- translational modes ----------------------
        self.prjT = [ [], [], [] ]
        for i in range( len( self.sele ) ):
            sm = math.sqrt( self.mass[i] )
            self.prjT[0] += [ sm, 0.0, 0.0 ]
            self.prjT[1] += [ 0.0, sm, 0.0 ]
            self.prjT[2] += [ 0.0, 0.0, sm ]
        for i in [0, 1, 2]:
            t = 0.0
            for k in range( self.size ):
                t += self.prjT[i][k] * self.prjT[i][k]
            for k in range( self.size ):
                self.prjT[i][k] /= t
        # -- velocities -------------------------------
        self.velo = []
        for i in range( self.size // 3 ):
            SD = math.sqrt( self.__kt / self.mass[i] )
            self.velo += [ qm3.maths.rand.gauss( 0.0, SD ) * 1.0e-2 for j in [0, 1, 2] ]
        for i in [0, 1, 2]:
            t = 0.0
            for j in range( self.size ):
                t += self.prjT[i][j] * self.velo[j]
            for j in range( self.size ):
                self.velo[j] -= t * self.prjT[i][j]
        scl = math.sqrt( self.temp / self.current_temperature()[0] )
        for i in range( self.size ):
            self.velo[i] *= scl
        # -- ring polymer -----------------------------
        self.rp_bead = num_beads
        self.rp_dime = len( pi_atoms )
        self.rp_scal = 1. / ( self.rp_bead * self.rp_dime )
        self.rp_atom = pi_atoms[:]
        self.rp_mota = [ self.sele.index( i ) for i in self.rp_atom ]
        cons = 2.0 * self.rp_bead * math.pow( self.temp * qm3.constants.KB * math.pi / qm3.constants.H, 2.0 ) * 1.0e-26
        self.rp_coor = {}
        self.rp_grad = {}
        self.rp_velo = {}
        self.rp_kumb = {}
        for i in range( self.rp_dime ):
            self.rp_kumb[self.rp_atom[i]] = cons * self.mass[self.rp_mota[i]]
            self.rp_coor[self.rp_atom[i]] = []
            self.rp_grad[self.rp_atom[i]] = []
            self.rp_velo[self.rp_atom[i]] = []
            i3  = 3 * self.rp_mota[i]
            cvm = math.sqrt( sum( [ j * j for j in self.velo[i3:i3+3] ] ) ) 
            for j in range( self.rp_bead ):
                self.rp_coor[self.rp_atom[i]] += self.coor[i3:i3+3]
                self.rp_grad[self.rp_atom[i]] += [ .0, .0, .0 ]
                tmp = [ -1.0 + qm3.maths.rand.random() * 2.0,
                        -1.0 + qm3.maths.rand.random() * 2.0,
                        -1.0 + qm3.maths.rand.random() * 2.0 ]
                pmt = math.sqrt( tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2] )
# -- remove translation on the initial beads velocities?
                for k in [0, 1, 2]:
                    self.rp_velo[self.rp_atom[i]].append( tmp[k] * cvm / pmt )



    def current_temperature( self ):
        Kin = 0.0
        for i in range( self.size // 3 ):
            i3  = i * 3
            Kin += self.mass[i] * ( self.velo[i3] * self.velo[i3] + self.velo[i3+1] * self.velo[i3+1] + self.velo[i3+2] * self.velo[i3+2] )
        return( Kin * self.__cc / float( self.size - 3.0 ), Kin * 0.005 )



    def get_grad_aver( self ):
        # ---------------------------------------------
        for i in range( len( self.sele ) ):
            i3 = 3 * i
            j3 = 3 * self.sele[i]
            for j in [0, 1, 2]:
                self.mole.coor[i3+j] = self.coor[j3+j]
        # ---------------------------------------------
        zero = [ 0.0 for i in range( 3 * self.mole.natm ) ]
        self.func = 0.0
        self.grad = [ 0.0 for i in range( self.size ) ]
        for i in self.rp_atom:
            i3 = 3 * i
            for j in range( self.rp_bead ):
                j3 = 3 * j
                self.mole.func = 0.0
                self.mole.grad = zero[:]
                self.mole.coor[i3:i3+3] = self.rp_coor[i][j3:j3+3]
                self.engn.get_grad( self.mole )
                self.func += self.mole.func * self.rp_scal
                for k in range( len( self.sele ) ):
                    k3 = 3 * k
                    m3 = 3 * self.sele[k]
                    for m in [0, 1, 2]:
                        self.grad[k3+m] += self.mole.grad[m3+m] * self.rp_scal
                # ---------------------------------------------
                if( j == 0 ):
                    mm = 3 * ( self.rp_bead - 1 )
                    pp = 3
                elif( j == self.rp_bead - 1 ):
                    mm = 3 * ( j - 1 )
                    pp = 0
                else:
                    mm = 3 * ( j - 1 )
                    pp = 3 * ( j + 1 )
                for k in [0, 1, 2]:
                    self.func += self.rp_kumb[i] * ( self.rp_coor[i][j3+k] - self.rp_coor[i][mm+k] )
                    self.rp_grad[i][j3+k] = self.mole.grad[i3+k] / self.rp_bead + 2.0 * self.rp_kumb[i] * ( 
                        2.0 * self.rp_coor[i][j3+k] - self.rp_coor[i][mm+k] - self.rp_coor[i][pp+k] )



    def get_grad_splt( self ):
        # ---------------------------------------------
        for i in range( len( self.sele ) ):
            i3 = 3 * i
            j3 = 3 * self.sele[i]
            for j in [0, 1, 2]:
                self.mole.coor[i3+j] = self.coor[j3+j]
        # ---------------------------------------------
        zero = [ 0.0 for i in range( 3 * self.mole.natm ) ]
        self.mole.func = 0.0
        self.mole.grad = zero[:]
        self.engn.get_grad( self.mole )
        self.func = self.mole.func
        self.grad = []
        for i in self.sele:
            i3 = 3 * i
            self.grad += self.mole.grad[i3:i3+3]
        # ---------------------------------------------
        for i in self.rp_atom:
            i3 = 3 * i
            for j in range( self.rp_bead ):
                j3 = 3 * j
                self.mole.grad = zero[:]
                self.mole.coor[i3:i3+3] = self.rp_coor[i][j3:j3+3]
                self.engn.get_grad( self.mole )
                if( j == 0 ):
                    mm = 3 * ( self.rp_bead - 1 )
                    pp = 3
                elif( j == self.rp_bead - 1 ):
                    mm = 3 * ( j - 1 )
                    pp = 0
                else:
                    mm = 3 * ( j - 1 )
                    pp = 3 * ( j + 1 )
                for k in [0, 1, 2]:
                    self.func += self.rp_kumb[i] * ( self.rp_coor[i][j3+k] - self.rp_coor[i][mm+k] )
                    self.rp_grad[i][j3+k] = self.mole.grad[i3+k] / self.rp_bead + 2.0 * self.rp_kumb[i] * (
                        2.0 * self.rp_coor[i][j3+k] - self.rp_coor[i][mm+k] - self.rp_coor[i][pp+k] )



    def integrate( self, step_size = 0.001, step_number = 1000, gamma_factor = 50.0, print_frequency = 100, log_function = default_log ):
        log_function( "----------------- Ring Polymer Molecular Dynamics: Langevin-Verlet (RPMD-NVT)\n" )
        log_function( "Step Size:          %20.10lg (ps)"%( step_size ) )
        log_function( "Temperature:        %20.10lg (K)"%( self.temp ) )
        log_function( "Gamma Factor:       %20.10lg (ps^-1)"%( gamma_factor ) )
        log_function( "Step Number:        %20d"%( step_number ) )
        log_function( "Print Frequency:    %20d"%( print_frequency ) )
        __ff  = step_size * gamma_factor
        if( __ff < 0.01 ):
            __ff = 0.01
            log_function( "\n>> Gamma factor:    %20.10lg (ps^-1)"%( 0.01 / step_size ) )
        log_function( "\n%20s%20s%20s%20s%20s"%( "Time (ps)", "Potential (kJ/mol)", "Kinetic (kJ/mol)", "Total (kJ/mol)", "Temperature (K)" ) )
        log_function( 100 * "-" )

        __c0  = math.exp( - __ff )
        __c1  = ( 1.0 - __c0 ) / __ff
        __c2  = ( 1.0 - __c1 ) / __ff
        __sr  = step_size * math.sqrt( ( 2.0 - ( 3.0 - 4.0 * __c0 + __c0 * __c0 ) / __ff ) / __ff )
        __sv  = math.sqrt( 1.0 - __c0 * __c0 )
        __cv1 = step_size * ( 1.0 - __c0 ) * ( 1.0 - __c0 ) / ( __ff * __sr * __sv )
        __cv2 = math.sqrt( 1.0 - __cv1 * __cv1 )
        __ff  = [ 1.0e-2 * math.sqrt( self.__kt / m ) for m in self.mass ]
        __fr1 = step_size * __c1
        __fv1 = step_size * ( __c1 - __c2 )
        __fv2 = step_size * __c2
        __fr2 = step_size * __fv2

        xavr = [ 0.0, 0.0, 0.0, 0.0 ]
        xrms = [ 0.0, 0.0, 0.0, 0.0 ]
        time = 0.0

        o_accl = []
        c_accl = []
        self.get_grad()
        for i in range( self.size // 3 ):
            i3 = 3 * i
            for j in [0, 1, 2]:
                c_accl.append( - self.grad[i3+j] / self.mass[i] * 100.0 )
                o_accl.append( 0.0 )
        bak = [ .0, .0, .0 ]
        for i in [0, 1, 2]:
            t = 0.0
            for j in range( self.size ):
                t += c_accl[j] * self.prjT[i][j]
            for j in range( self.size ):
                c_accl[j] -= t * self.prjT[i][j]
            bak[i] = t

        kab = bak[:]
        crp_accl = {}
        orp_accl = {}
        for i in range( self.rp_dime ):
            crp_accl[self.rp_atom[i]] = []
            orp_accl[self.rp_atom[i]] = []
            i3 = 3 * self.rp_mota[i]
            bak = kab[:]
            for j in [0, 1, 2]:
                bak[j] -= self.prjT[j][i3+j] * c_accl[i3+j]
            for j in range( self.rp_bead ):
                j3 = 3 * j
                for k in [0, 1, 2]:
                    crp_accl[self.rp_atom[i]].append( - self.rp_grad[self.rp_atom[i]][j3+k] / self.mass[self.rp_mota[i]] * 100.0 )
                    orp_accl[self.rp_atom[i]].append( 0.0 )
                for ii in [0, 1, 2]:
                    t = bak[ii]
                    for jj in [0, 1, 2]:
                        t += self.prjT[ii][i3+jj] * crp_accl[self.rp_atom[i]][j3+jj]
                    for jj in [0, 1, 2]:
                        crp_accl[self.rp_atom[i]][j3+jj] -= t * self.prjT[ii][i3+jj]

        T, Kin = self.current_temperature()
        log_function( "%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf"%( 0.0, self.func, Kin, self.func + Kin, T ) )
        xavr[0] += self.func
        xavr[1] += Kin
        xavr[2] += self.func + Kin
        xavr[3] += T
        xrms[0] += self.func * self.func
        xrms[1] += Kin * Kin
        xrms[2] += ( self.func + Kin ) * ( self.func + Kin )
        xrms[3] += T * T

        for istep in range( step_number ):
            time += step_size

            for i in range( self.size ):
                self.coor[i] += __fr1 * self.velo[i] + __fr2 * c_accl[i]
                r1 = qm3.maths.rand.gauss( 0.0, 1.0 )
                r2 = qm3.maths.rand.gauss( 0.0, 1.0 )
                self.coor[i] += __ff[i//3] * __sr * r1
                o_accl[i] = __c0 * self.velo[i] + __fv1 * c_accl[i] + __ff[i//3] * __sv * ( __cv1 * r1 + __cv2 * r2 )

            for i in range( self.rp_dime ):
                for j in range( 3 * self.rp_bead ):
                    self.rp_coor[self.rp_atom[i]][j] += __fr1 * self.rp_velo[self.rp_atom[i]][j] + __fr2 * crp_accl[self.rp_atom[i]][j]
                    r1 = qm3.maths.rand.gauss( 0.0, 1.0 )
                    r2 = qm3.maths.rand.gauss( 0.0, 1.0 )
                    self.rp_coor[self.rp_atom[i]][j] += __ff[self.rp_mota[i]] * __sr * r1
                    orp_accl[self.rp_atom[i]][j] = __c0 * self.rp_velo[self.rp_atom[i]][j] + __fv1 * crp_accl[self.rp_atom[i]][j] + __ff[self.rp_mota[i]] * __sv * ( __cv1 * r1 + __cv2 * r2 )

            self.get_grad()

            for i in range( self.size ):
                c_accl[i] = - self.grad[i] / self.mass[i//3] * 100.0

            for i in [0, 1, 2]:
                t = 0.0
                for j in range( self.size ):
                    t += c_accl[j] * self.prjT[i][j]
                for j in range( self.size ):
                    c_accl[j] -= t * self.prjT[i][j]
                bak[i] = t

            kab = bak[:]
            for i in range( self.rp_dime ):
                i3 = 3 * self.rp_mota[i]
                bak = kab[:]
                for j in [0, 1, 2]:
                    bak[j] -= self.prjT[j][i3+j] * c_accl[i3+j]
                for j in range( self.rp_bead ):
                    j3 = 3 * j
                    for k in [0, 1, 2]:
                        crp_accl[self.rp_atom[i]][j3+k] = - self.rp_grad[self.rp_atom[i]][j3+k] / self.mass[self.rp_mota[i]] * 100.0
                    for ii in [0, 1, 2]:
                        t = bak[ii]
                        for jj in [0, 1, 2]:
                            t += self.prjT[ii][i3+jj] * crp_accl[self.rp_atom[i]][j3+jj]
                        for jj in [0, 1, 2]:
                            crp_accl[self.rp_atom[i]][j3+jj] -= t * self.prjT[ii][i3+jj]

            for i in range( self.size ):
                self.velo[i] = o_accl[i] + __fv2 * c_accl[i]
            for i in [0, 1, 2]:
                t = 0.0
                for j in range( self.size ):
                    t += self.velo[j] * self.prjT[i][j]
                for j in range( self.size ):
                    self.velo[j] -= t * self.prjT[i][j]
                bak[i] = t

            kab = bak[:]
            for i in range( self.rp_dime ):
                i3 = 3 * self.rp_mota[i]
                bak = kab[:]
                for j in [0, 1, 2]:
                    bak[j] -= self.prjT[j][i3+j] * self.velo[i3+j]
                for j in range( self.rp_bead ):
                    j3 = 3 * j
                    for k in [0, 1, 2]:
                        self.rp_velo[self.rp_atom[i]][j3+k] = orp_accl[self.rp_atom[i]][j3+k] + __fv2 * crp_accl[self.rp_atom[i]][j3+k]
                for ii in [0, 1, 2]:
                    t = bak[ii]
                    for jj in [0, 1, 2]:
                        t += self.prjT[ii][i3+jj] * self.rp_velo[self.rp_atom[i]][j3+jj]
                    for jj in [0, 1, 2]:
                        self.rp_velo[self.rp_atom[i]][j3+jj] -= t * self.prjT[ii][i3+jj]

            t = 1. / self.rp_bead
            for i in range( self.rp_dime ):
                i3 = 3 * self.rp_mota[i]
                self.coor[i3:i3+3] = [ .0, .0, .0 ]
                self.velo[i3:i3+3] = [ .0, .0, .0 ]
                o_accl[i3:i3+3] = [ .0, .0, .0 ]
                c_accl[i3:i3+3] = [ .0, .0, .0 ]
                for j in range( self.rp_bead ):
                    j3 = 3 * j
                    for k in [0, 1, 2]:
                        self.coor[i3+k] += t * self.rp_coor[self.rp_atom[i]][j3+k]
                        self.velo[i3+k] += t * self.rp_velo[self.rp_atom[i]][j3+k]
                        o_accl[i3+k]    += t * orp_accl[self.rp_atom[i]][j3+k]
                        c_accl[i3+k]    += t * crp_accl[self.rp_atom[i]][j3+k]

            T, Kin = self.current_temperature()
            xavr[0] += self.func
            xavr[1] += Kin
            xavr[2] += self.func + Kin
            xavr[3] += T
            xrms[0] += self.func * self.func
            xrms[1] += Kin * Kin
            xrms[2] += ( self.func + Kin ) * ( self.func + Kin )
            xrms[3] += T * T
            if( istep % print_frequency == 0 or True ):
                log_function( "%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf"%( time, self.func, Kin, self.func + Kin, T ) )
            self.current_step( istep )

        log_function( 100 * "-" )
        savr = "%-20s"%( "Averages:" )
        srms = "%-20s"%( "RMS Deviations:" )
        for i in range( len( xavr ) ):
            tavr  = xavr[i] / float( step_number + 1 )
            savr += "%20.5lf"%( tavr )
            trms  = xrms[i] / float( step_number + 1 )
            srms += "%20.5lf"%( math.sqrt( math.fabs( trms - tavr * tavr ) ) )
        log_function( savr )
        log_function( srms )
        log_function( 100 * "-" + "\n" )

# -*- coding: iso-8859-1 -*-
import  sys
import  math
import  qm3.maths.rand
import  qm3.maths.matrix
import  qm3.utils
import  qm3.constants
import  time
import  os



def default_log( txt ):
    sys.stdout.write( txt + "\n" )
    sys.stdout.flush()



class model( object ):
    def __init__( self, molec ):
        """
        define as much engines as needed based on the molec
        """
        self.sele = []
        self.size = 3 * len( self.sele )


    def get_grad( self, molec ):
        """
        sequencially accumulate all the engines.get_grad
        """
        molec.func = 0
        molec.grad = [ 0.0 for i in range( 3 * molec.natm ) ]


    def get_hess( self, molec, nder = 1.0e-4 ):
        """
        defaults to numerical hessian
        """
        hh = []
        for j in self.sele:
            j3 = j * 3
            for k in [0, 1, 2]:
                bak = molec.coor[j3+k]
                molec.coor[j3+k] = bak + nder
                self.get_grad( molec )
                gp = []
                for l in self.sele:
                    l3 = l * 3
                    gp += molec.grad[l3:l3+3][:]
                molec.coor[j3+k] = bak - nder
                self.get_grad( molec )
                gm = []
                for l in self.sele:
                    l3 = l * 3
                    gm += molec.grad[l3:l3+3][:]
                molec.coor[j3+k] = bak
                hh.append( [ ( gp[l] - gm[l] ) / ( 2.0 * nder ) for l in range( self.size ) ] )
        for j in range( self.size ):
            for k in range( self.size ):
                if( j != k ):
                    t = 0.5 * ( hh[j][k] + hh[k][j] )
                    hh[j][k] = t
                    hh[k][j] = t
        molec.hess = []
        for j in range( self.size ):
            molec.hess += hh[j]




class dynamics( object ):
    def __init__( self, mole, sele, engn ):
        self.mole = mole
        self.sele = sele[:]
        self.engn = engn
        self.size = 3 * len( sele )
        self.mass = []
        self.coor = []
        for i in self.sele:
            i3 = i * 3
            self.mass.append( self.mole.mass[i] )
            self.coor += self.mole.coor[i3:i3+3]



    def current_step( self, istep ):
#        f = open( "output", "at" )
#        f.write( "%d\n\n"%( self.mole.natm + ( self.rp_bead - 1 ) * self.rp_dime ) )
#        for i in range( self.mole.natm ):
#            if( not i in self.rp_atom ):
#                f.write( "%-2s%20.10lf%20.10lf%20.10lf\n"%( self.mole.labl[i],
#                    self.coor[3*i], self.coor[3*i+1], self.coor[3*i+2] ) )
#        for i in self.rp_atom:
#            for j in range( self.rp_bead ):
#                f.write( "%-2s%20.10lf%20.10lf%20.10lf\n"%( self.mole.labl[i],
#                    self.rp_coor[i][3*j], self.rp_coor[i][3*j+1], self.rp_coor[i][3*j+2] ) )
#        f.close()
        pass



    def setup( self, pi_atoms, num_beads = 8, temperature = 300.0 ):
        self.temp = temperature
        self.__kt = 1000.0 * qm3.constants.NA * qm3.constants.KB * self.temp
        self.__cc = 10.0 / ( qm3.constants.KB * qm3.constants.NA )
        # -- translational modes ----------------------
        self.prjT = [ [], [], [] ]
        mt = math.sqrt( sum( self.mass ) )
        for i in range( len( self.sele ) ):
            sm = math.sqrt( self.mass[i] )
            self.prjT[0] += [ sm / mt, 0.0, 0.0 ]
            self.prjT[1] += [ 0.0, sm / mt, 0.0 ]
            self.prjT[2] += [ 0.0, 0.0, sm / mt ]
        # -- velocities -------------------------------
        self.velo = []
        for i in range( self.size // 3 ):
            SD = math.sqrt( self.__kt / self.mass[i] )
            self.velo += [ qm3.maths.rand.gauss( 0.0, SD ) * 1.0e-2 for j in [0, 1, 2] ]
        for i in [0, 1, 2]:
            t = qm3.maths.matrix.dot_product( self.prjT[i], self.velo )
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



    def get_grad( self ):
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
                    self.func += self.rp_kumb[i] * math.pow( self.rp_coor[i][j3+k] - self.rp_coor[i][mm+k], 2.0 )
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
            t = qm3.maths.matrix.dot_product( self.prjT[i], c_accl )
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
                t = qm3.maths.matrix.dot_product( self.prjT[i], c_accl )
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
                t = qm3.maths.matrix.dot_product( self.prjT[i], self.velo )
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




#
# J. Phys. Chem. Lett. v7, p4374 (2016) [10.1021/acs.jpclett.6b02115]
# J. Chem. Phys. v134, p184107 (2011) [10.1063/1.3587240]
# J. Chem. Phys. v148, p102334 (2018) [10.1063/1.5007180]
#
class instanton( object ):
    def __init__( self, mole, sele, engn, num_beads = 64, temperature = 300.0 ):
        """
U = 1 over P sum from{ i=1 } to { P }{ V_i left(x_{i,1}, dotslow ,x_{i,3N} right)}+
sum from{ i=1 } to { P }{ sum from{ j=1 } to { 3N }{
k color gray { left lbrace  {2 P %pi^2 k_B^2 T^2 } over {h^2} 10^-26 right rbrace }
m_{j/3} left( x_{i,j} - x_{i-1,j} right)^2  } }
newline
{partial U }over{partial x_{i,k} } = 1 over P {partial V_i }over{partial x_{i,k} }left(x_{i,1}, dotslow
,x_{i,3N} right)+2 k m_{ k/3 }left( 2x_{i,k} - x_{i-1,k} - x_{i+1,k} right)
newline
{partial^2 U} over {partial x_{i,k} partial x_{j,l}} =
%delta_{j=i,l,k} over P 
{partial^2 V_i} over {partial x_{i,k} partial x_{j,l}}
left(x_{i,1}, dotslow ,x_{i,3N} right) +4 k m_{ k/3 } %delta_{j=i,l=k} -2 k m_{ k/3 } %delta_{j=i+1,l=k}
-2 k m_{ k/3 } %delta_{j=i-1,l=k}
        """
        self.mole = mole
        self.sele = sele[:]
        self.engn = engn
        self.mass = [ self.mole.mass[i] for i in self.sele ]
        self.temp = temperature
        self.bead = num_beads
        self.half = self.bead // 2
        self.kumb = 2.0 * self.bead * math.pow( self.temp * qm3.constants.KB * math.pi / qm3.constants.H, 2.0 ) * 1.0e-26
        self.disp = 3 * len( self.sele )
        self.size = self.disp * ( self.half + 1 )
        print( "[RP] bead:", self.bead, self.half + 1 )
        print( "[RP] temp: %.2lf _K"%( self.temp ) )
        print( "[RP] kumb: %.2lf _kJ/(mol A^2)"%( self.kumb ) )



    def __rotations( self, coor, symm = 1.0 ):
        kk = ( 8.0 * math.pi * math.pi * qm3.constants.KB * self.temp ) / ( qm3.constants.H * qm3.constants.H * qm3.constants.NA ) * 1.0e-23
        mt = sum( self.mass )
        mc = [ 0.0, 0.0, 0.0 ]
        for j in self.sele:
            j3 = j * 3
            for k in [0, 1, 2]:
                mc[k] += self.mole.mass[j] * coor[j3+k]
        mc[0] /= mt; mc[1] /= mt; mc[2] /= mt
        xx = 0.0; xy = 0.0; xz = 0.0; yy = 0.0; yz = 0.0; zz = 0.0
        for j in self.sele:
            j3 = j * 3
            xx += self.mole.mass[j] * ( coor[j3]   - mc[0] ) * ( coor[j3]   - mc[0] )
            xy += self.mole.mass[j] * ( coor[j3]   - mc[0] ) * ( coor[j3+1] - mc[1] )
            xz += self.mole.mass[j] * ( coor[j3]   - mc[0] ) * ( coor[j3+2] - mc[2] )
            yy += self.mole.mass[j] * ( coor[j3+1] - mc[1] ) * ( coor[j3+1] - mc[1] )
            yz += self.mole.mass[j] * ( coor[j3+1] - mc[1] ) * ( coor[j3+2] - mc[2] )
            zz += self.mole.mass[j] * ( coor[j3+2] - mc[2] ) * ( coor[j3+2] - mc[2] )
        val, vec = qm3.maths.matrix.diag( qm3.maths.matrix.from_upper_diagonal_rows( [ yy + zz, -xy, -xz, xx + zz, -yz, xx + yy ], 3 ), 3 )
        return( math.log( math.sqrt( math.pi * kk * kk * kk * val[0] * val[1] * val[2] ) / symm ) )



    def calc_tst( self, r_coor, r_func, r_hess, t_coor, t_func, t_hess, t_symm = 1.0 ):
        """
Q_rot = 1 over %sigma left( %pi left( {8 %pi^2 k_B T 10^-23} over {h^2 N_A} right)^3 det I right)^{ 1 over 2 }
~~~~~~~
size 10 { I =  sum from{j=1} to{N}{ m_j left[ left( {vec{r}}_{j} cdot {vec{r}}_{j} right) I_3 - {vec{r}}_{j} times {vec{r}}_{j} right] } }
~~~~~~~
size 10 { {vec{r}}_{j} = left( x_{j},y_{j},z_{j} right) - left( x_{CM},y_{CM},z_{CM} right) }
newline
Q_vib =  prod from{k=1} to{3N-6/7} { 1 over {2 sinh left( 1 over 2 {h  %ípsilon_k 100 c} over {k_B T} right) } }
~~~~~~~
k_TST = {k_B T} over h { Q_rot^{%Ux2021 } · Q_vib^{%Ux2021 } } over { Q_rot^{R } · Q_vib^{R } } e^{ - {{V^{ %Ux2021 } - V^R} over {k_B T}} 10^3 }
        """
        # activation potential energy
        efunc = - ( t_func - r_func ) * 1000.0 / ( self.temp * qm3.constants.KB * qm3.constants.NA )
        print( "[TS]dfunc: %20.10le (%.2lf _kJ/mol)"%( efunc, t_func - r_func ) )
        # rotational partition function
        rQR = self.__rotations( r_coor )
        print( "[TS]l_rQR: %20.10le"%( rQR ) )
        tQR = self.__rotations( t_coor, t_symm )
        print( "[TS]l_tQR: %20.10le"%( tQR ) )
        # vibrational partition function
        kk = 100.0 * qm3.constants.C * qm3.constants.H / ( qm3.constants.KB * self.temp )
        cc = []
        for j in self.sele:
            j3 = j * 3
            cc += r_coor[j3:j3+3][:]
        frq = qm3.utils.hessian_frequencies( self.mass, cc, r_hess )[0]
        rQV = 0.0
        for f in frq[6:]:
            rQV -= math.log( 2.0 * math.sinh( f * kk * 0.5 ) )
        print( "[TS]rfreq: " + ", ".join( [ "%.1lf"%( math.fabs( i ) ) for i in frq[0:7] ] ) + " _cm^-1" )
        print( "[TS]l_rQV: %20.10le"%( rQV ) )
        cc = []
        for j in self.sele:
            j3 = j * 3
            cc += t_coor[j3:j3+3][:]
        frq, vec = qm3.utils.hessian_frequencies( self.mass, cc, t_hess )
        self.tst_mode = [ vec[i*self.disp] for i in range( self.disp ) ]
        t = math.sqrt( sum( [ i * i for i in self.tst_mode ] ) )
        self.tst_mode = [ i / t for i in self.tst_mode ]
        tQV = 0.0
        for f in frq[7:]:
            tQV -= math.log( 2.0 * math.sinh( f * kk * 0.5 ) )
        print( "[TS]tfreq: " + ", ".join( [ "%.1lf"%( math.fabs( i ) ) for i in frq[0:8] ] ) + " _cm^-1" )
        print( "[RP] T_co: %.2lf _K"%( math.fabs( frq[0] ) * 100.0 * qm3.constants.C * qm3.constants.H / ( 2.0 * math.pi * qm3.constants.KB ) ) )
        print( "[TS]l_tQV: %20.10le"%( tQV ) )
        # kinetic constant
        self.T_cons = qm3.constants.KB * self.temp / qm3.constants.H * math.exp( tQV + tQR - rQV - rQR + efunc )
        print( "[TS] kcin: %20.10le _1/s"%( self.T_cons  ) )
        # initially map TS coordinates
        self.mole.coor = t_coor



    def setup( self, step_size = 0.3 ):
        # build the ring polymers supermolecule (= num_beads / 2 + 1 molecules)
        self.coor = []
        dsp = 2 * step_size / self.half
        for i in range( self.half + 1 ):
            for j in range( len( self.sele ) ):
                j3 = j * 3
                J3 = self.sele[j] * 3
                for k in [0, 1, 2]:
                    self.coor.append( self.mole.coor[J3+k] + ( i * dsp - step_size ) * self.tst_mode[j3+k] )



    def current_step( self, istep ):
        pass
#        f = open( "output", "at" )
#        f.write( "%d\n\n"%( self.size // 3 ) )
#        for i in range( self.half + 1 ):
#            i_cc = i * self.disp
#            for j in range( len( self.sele ) ):
#                j3 = i_cc + j * 3
#                f.write( "%-2s%20.10lf%20.10lf%20.10lf\n"%( self.mole.labl[self.sele[j]][0],
#                    self.coor[j3], self.coor[j3+1], self.coor[j3+2] ) )
#        f.close()



    def get_grad( self ):
        self.ener = []
        self.func = 0.0
        self.grad = [ 0.0 for i in range( self.size ) ]
        for i in range( self.half + 1 ):
            i_cc = i * self.disp
            if( i == 0 ):
                scal = 1.0
                i_mm = ( i + 1 ) * self.disp
                i_pp = ( i + 1 ) * self.disp
            elif( i == self.half ):
                scal = 1.0
                i_mm = ( i - 1 ) * self.disp
                i_pp = ( i - 1 ) * self.disp
            else:
                scal = 2.0
                i_mm = ( i - 1 ) * self.disp
                i_pp = ( i + 1 ) * self.disp
            # map current polymer into molecule coordinates
            for j in range( len( self.sele ) ):
                j3 = i_cc + j * 3
                J3 = self.sele[j] * 3
                self.mole.coor[J3:J3+3] = self.coor[j3:j3+3]
            # get B.O. potential for current polymer
            self.engn.get_grad( self.mole )
            self.ener.append( self.mole.func )
            self.func += scal * self.mole.func / self.bead
            for j in range( len( self.sele ) ):
                j3 = i_cc + j * 3
                J3 = self.sele[j] * 3
                for k in [0, 1, 2]:
                    self.grad[j3+k] += self.mole.grad[J3+k] / self.bead
            # get R.P. potential among polymers
            for j in range( len( self.sele ) ):
                j3 = j * 3
                kk = self.mass[j] * self.kumb
                for k in [0, 1, 2]:
                    self.func += scal * kk * math.pow( self.coor[i_cc+j3+k] - self.coor[i_mm+j3+k], 2.0 )
                    self.grad[i_cc+j3+k] += 2.0 * kk * ( 2.0 * self.coor[i_cc+j3+k] - self.coor[i_mm+j3+k] - self.coor[i_pp+j3+k] )



    def get_hess( self ):
        if( not os.path.isfile( "update.dump" ) ):
            z    = [ 0.0 for i in range( self.size ) ]
            hess = [ z[:] for i in range( self.size ) ]
            for i in range( self.half + 1 ):
                i_cc = i * self.disp
                for j in range( len( self.sele ) ):
                    j3 = i_cc + j * 3
                    J3 = self.sele[j] * 3
                    self.mole.coor[J3:J3+3] = self.coor[j3:j3+3]
                # ---------------------------------
                self.engn.get_hess( self.mole )
                hh = []
                l  = 0
                for j in range( self.disp ):
                    hh.append( [] )
                    for k in range( self.disp ):
                        hh[-1].append( self.mole.hess[l] / self.bead )
                        l += 1
                # ---------------------------------
                if( i == 0 ):
                    i_mm = self.size - self.disp
                    i_pp = self.disp
                elif( i == self.half ):
                    i_mm = i_cc - self.disp
                    i_pp = 0
                else:
                    i_mm = i_cc - self.disp
                    i_pp = i_cc + self.disp
                for j in range( self.disp ):
                    for k in range( self.disp ):
                        hess[i_cc+j][i_cc+k] = hh[j][k]
                        if( j == k ):
                            t = self.kumb * self.mass[j//3]
                            hess[i_cc+j][i_cc+k] += 4.0 * t
                            hess[i_cc+j][i_mm+k] = -2.0 * t
                            hess[i_cc+j][i_pp+k] = -2.0 * t
            self.get_grad()
            self.hess = []
            for j in range( self.size ):
                for k in range( self.size ):
                    self.hess.append( hess[j][k] )
            qm3.utils.manage_hessian( self.coor, self.grad, self.hess, should_update = False )
        else:
            self.get_grad()
            self.hess = [ .0 for i in range( self.size * self.size ) ]
            qm3.utils.manage_hessian( self.coor, self.grad, self.hess, should_update = True )
        qm3.utils.raise_hessian_RT( self.mass * ( self.half + 1 ), self.coor, self.hess )



    def calc_rpt( self, r_coor, r_func, r_hess ):
        """
Q_rot = left[ 1 over %sigma right]_R left( %pi left( {8 %pi^2 P k_B T 10^-23} over{h^2 N_A} right)^3 det I right)^{ 1 over 2 }
~~~~~
size 10 { 
I = sum from{ i=1 } to{ P }{ sum from{j=1} to{N}{ m_j left[ left( {vec{r}}_{i,j} cdot
{vec{r}}_{i,j} right) I_3 - {vec{r}}_{i,j} times {vec{r}}_{i,j} right] } } }
~~~~~
size 9{ {vec{r}}_{i,j} = left( x_{i,j},y_{i,j},z_{i,j} right) - left( x_{CM},y_{CM},z_{CM} right) }
newline
Q_vib = prod from{k=1} to{3N cdot P-6} { 1 over {2 sinh left( 1 over 2 {h lline %ípsilon_k rline 100 c } over {P k_B T} right) } }
newline
k_{RP} = {k_B T P} over h
left( {2 B %pi k_B T P 10^-23} over{h^2 N_A} right)^{1 over 2} ~~ 
{ Q_rot^{%Ux2021 } · Q_vib^{%Ux2021 } } over { Q_rot^{R } · Q_vib^{R } } e^{ - {{U^{ %Ux2021 } - V^R} over {k_B T}} 10^3 }
~~~~~~~~
size 10 { B= sum from{ i=1 } to { P }{ sum from{j=1} to{3N} { m_{j/3} left( x_{i,j} - x_{i-1,j} right)^2 } } }
        """
        # activation potential energy
        efunc = - ( self.func - r_func ) * 1000.0 / ( self.temp * qm3.constants.KB * qm3.constants.NA )
        print( "[RP]dfunc: %20.10le (%.2lf _kJ/mol)"%( efunc, self.func - r_func ) )
        # rotational partition function
        rQR = self.__rotations( r_coor ) + 3.0 * math.log( self.bead )
        print( "[RP]l_rQR: %20.10le"%( rQR ) )
        mt = self.bead * sum( self.mass )
        mc = [ 0.0, 0.0, 0.0 ]
        for i in range( self.half + 1 ):
            i_cc = i * self.disp
            if( i == 0 ):
                scal = 1.0
            elif( i == self.half ):
                scal = 1.0
            else:
                scal = 2.0
            for j in range( len( self.sele ) ):
                j3 = i_cc + j * 3
                for k in [0, 1, 2]:
                    mc[k] += scal * self.mass[j] * self.coor[j3+k]
        mc[0] /= mt; mc[1] /= mt; mc[2] /= mt
        xx = 0.0; xy = 0.0; xz = 0.0; yy = 0.0; yz = 0.0; zz = 0.0
        for i in range( self.half + 1 ):
            i_cc = i * self.disp
            if( i == 0 ):
                scal = 1.0
            elif( i == self.half ):
                scal = 1.0
            else:
                scal = 2.0
            for j in range( len( self.sele ) ):
                j3 = i_cc + j * 3
                xx += scal * self.mass[j] * ( self.coor[j3]   - mc[0] ) * ( self.coor[j3]   - mc[0] )
                xy += scal * self.mass[j] * ( self.coor[j3]   - mc[0] ) * ( self.coor[j3+1] - mc[1] )
                xz += scal * self.mass[j] * ( self.coor[j3]   - mc[0] ) * ( self.coor[j3+2] - mc[2] )
                yy += scal * self.mass[j] * ( self.coor[j3+1] - mc[1] ) * ( self.coor[j3+1] - mc[1] )
                yz += scal * self.mass[j] * ( self.coor[j3+1] - mc[1] ) * ( self.coor[j3+2] - mc[2] )
                zz += scal * self.mass[j] * ( self.coor[j3+2] - mc[2] ) * ( self.coor[j3+2] - mc[2] )
        val, vec = qm3.maths.matrix.diag( qm3.maths.matrix.from_upper_diagonal_rows( [ yy + zz, -xy, -xz, xx + zz, -yz, xx + yy ], 3 ), 3 )
        kk = ( 8.0 * math.pi * math.pi * qm3.constants.KB * self.temp * self.bead ) / ( qm3.constants.H * qm3.constants.H * qm3.constants.NA ) * 1.0e-23
        tQR = math.log( math.sqrt( math.pi * kk * kk * kk * val[0] * val[1] * val[2] ) )
        print( "[RP]l_tQR: %20.10le"%( tQR ) )
        # vibrational partition function
        size = self.bead * self.disp
        z    = [ 0.0 for i in range( size ) ]
        kk   = 100.0 * qm3.constants.C * qm3.constants.H / ( self.bead * qm3.constants.KB * self.temp )
        # ---- reactants
        hess = [ z[:] for i in range( size ) ]
        hh   = []
        l    = 0
        for j in range( self.disp ):
            hh.append( [] )
            for k in range( self.disp ):
                hh[-1].append( r_hess[l] / self.bead )
                l += 1
        for i in range( self.bead ):
            i_cc = i * self.disp
            if( i == 0 ):
                i_mm = size - self.disp
                i_pp = i_cc + self.disp
            elif( i == self.bead - 1 ):
                i_mm = i_cc - self.disp
                i_pp = 0
            else:
                i_mm = i_cc - self.disp
                i_pp = i_cc + self.disp
            for j in range( self.disp ):
                for k in range( self.disp ):
                    hess[i_cc+j][i_cc+k] = hh[j][k]
                    if( j == k ):
                        t = self.kumb * self.mass[j//3]
                        hess[i_cc+j][i_cc+k] += 4.0 * t
                        hess[i_cc+j][i_mm+k] = -2.0 * t
                        hess[i_cc+j][i_pp+k] = -2.0 * t
        HESS = []
        for j in range( size ):
            for k in range( size ):
                HESS.append( hess[j][k] )
        coor = []
        for j in self.sele:
            j3 = j * 3
            coor += r_coor[j3:j3+3]
        frq = qm3.utils.hessian_frequencies( self.mass * self.bead, coor * self.bead, HESS, True )[0]
        rQV = 0.0
        for f in frq[6:]:
            rQV -= math.log( 2.0 * math.sinh( f * kk * 0.5 ) )
        print( "[RP]rfreq: " + ", ".join( [ "%.1lf"%( math.fabs( i ) ) for i in frq[0:7] ] ) + " _cm^-1" )
        print( "[RP]l_rQV: %20.10le"%( rQV ) )
        # ---- instanton
        hess = [ z[:] for i in range( size ) ]
        t0 = time.time()
        for i in range( self.half + 1 ):
            i_cc = i * self.disp
            for j in range( len( self.sele ) ):
                j3 = i_cc + j * 3
                J3 = self.sele[j] * 3
                self.mole.coor[J3:J3+3] = self.coor[j3:j3+3]
            # ---------------------------------
            self.engn.get_hess( self.mole )
            hh = []
            l  = 0
            for j in range( self.disp ):
                hh.append( [] )
                for k in range( self.disp ):
                    hh[-1].append( self.mole.hess[l] / self.bead )
                    l += 1
            # ---------------------------------
            if( i == 0 ):
                i_mm = size - self.disp
                i_pp = i_cc + self.disp
            else:
                i_mm = i_cc - self.disp
                i_pp = i_cc + self.disp
            for j in range( self.disp ):
                for k in range( self.disp ):
                    hess[i_cc+j][i_cc+k] = hh[j][k]
                    if( j == k ):
                        t = self.kumb * self.mass[j//3]
                        hess[i_cc+j][i_cc+k] += 4.0 * t
                        hess[i_cc+j][i_mm+k] = -2.0 * t
                        hess[i_cc+j][i_pp+k] = -2.0 * t
            if( i != 0 and i != self.half ):
                i_cc = ( self.bead - i ) * self.disp
                if( i == 1 ):
                    i_mm = i_cc - self.disp
                    i_pp = 0
                else:
                    i_mm = i_cc - self.disp
                    i_pp = i_cc + self.disp
                for j in range( self.disp ):
                    for k in range( self.disp ):
                        hess[i_cc+j][i_cc+k] = hh[j][k]
                        if( j == k ):
                            t = self.kumb * self.mass[j//3]
                            hess[i_cc+j][i_cc+k] += 4.0 * t
                            hess[i_cc+j][i_mm+k] = -2.0 * t
                            hess[i_cc+j][i_pp+k] = -2.0 * t
        print( "(RP_hess: %.1lf _seg)"%( time.time() - t0 ) )
        HESS = []
        for j in range( size ):
            for k in range( size ):
                HESS.append( hess[j][k] )
        coor = self.coor[:]
        for i in range( self.half - 1, 0, -1 ):
            i_cc = i * self.disp
            coor += self.coor[i_cc:i_cc+self.disp]
        t0 = time.time()
        frq = qm3.utils.hessian_frequencies( self.mass * self.bead, coor, HESS, True )[0]
        print( "(RP_diaq: %.1lf _seg)"%( time.time() - t0 ) )
        tQV = 0.0
        for f in frq[8:]:
            tQV -= math.log( 2.0 * math.sinh( f * kk * 0.5 ) )
        tQV -= math.log( 2.0 * math.sinh( math.fabs( frq[0] ) * kk * 0.5 ) )
        print( "[RP]tfreq: " + ", ".join( [ "%.1lf"%( math.fabs( i ) ) for i in frq[0:9] ] ) + " _cm^-1" )
        print( "[RP]l_tQV: %20.10le"%( tQV ) )
        # beads rotational exchange
        tQP = 0.0
        for i in range( self.half + 1 ):
            i_cc = i * self.disp
            if( i == 0 ):
                scal = 1.0
                i_mm = ( i + 1 ) * self.disp
            elif( i == self.half ):
                scal = 1.0
                i_mm = ( i - 1 ) * self.disp
            else:
                scal = 2.0
                i_mm = ( i - 1 ) * self.disp
            for j in range( len( self.sele ) ):
                j3 = j * 3
                for k in [0, 1, 2]:
                    tQP += scal * self.mass[j] * math.pow( self.coor[i_cc+j3+k] - self.coor[i_mm+j3+k], 2.0 )
        kk = 2.0 * math.pi * self.bead * self.temp * qm3.constants.KB * 1.0e-23 / ( qm3.constants.NA * qm3.constants.H * qm3.constants.H )
        tQP = 0.5 * math.log( kk * tQP )
        print( "[RP]l_tQP: %20.10le"%( tQP ) )
        # kinetic constant
        self.R_cons = qm3.constants.KB * self.temp * self.bead / qm3.constants.H * math.exp( tQP + tQV + tQR - rQV - rQR + efunc )
        print( "[RP] kcin: %20.10le _1/s"%( self.R_cons  ) )
        print( "    kappa: %20.10le"%( self.R_cons / self.T_cons ) )
            


# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.maths.rand
import qm3.maths.matrix
import qm3.utils
import qm3.constants



_Tfac = 10.0 / ( qm3.constants.KB * qm3.constants.NA )
_Gfac = 1000.0 * qm3.constants.NA  # g/mol > kg



def default_log( txt ):
    sys.stdout.write( txt + "\n" )
    sys.stdout.flush()



def current_temperature( obj, project = True ):
    KEf = 0.0
    for i in range( obj.size // 3 ):
        i3  = i * 3
        KEf += obj.mass[i] * ( obj.velo[i3] * obj.velo[i3] + obj.velo[i3+1] * obj.velo[i3+1] + obj.velo[i3+2] * obj.velo[i3+2] )
#    KEf = sum( [ obj.mass[ii] * sum( [ jj*jj for jj in obj.velo[3*ii:3*ii+3] ] ) for ii in range( obj.size // 3 ) ] )
    if( project ):
#        T = KEf * _Tfac / float( obj.size - 6 )
        T = KEf * _Tfac / float( obj.size - 3 )
    else:
        T = KEf * _Tfac / float( obj.size )
    Kin = KEf * 0.005
    return( T, Kin )



def assign_velocities( obj, temperature = 300.0, project = True ):
    obj.velo = []
    KT = qm3.constants.KB * temperature * _Gfac
    for i in range( obj.size // 3 ):
        SD = math.sqrt( KT / obj.mass[i] )
        obj.velo += [ qm3.maths.rand.gauss( 0.0, SD ) * 1.0e-2 for ii in [0, 1, 2] ]     # ang/ps
    if( project ):
#        qm3.utils.project_RT_modes( obj.mass, obj.coor, obj.velo, [] )
        obj.prjT = [ [], [], [] ]
        mt = math.sqrt( sum( obj.mass ) )
        for i in range( obj.size // 3 ):
            sm = math.sqrt( obj.mass[i] )
            obj.prjT[0] += [ sm / mt, 0.0, 0.0 ]
            obj.prjT[1] += [ 0.0, sm / mt, 0.0 ]
            obj.prjT[2] += [ 0.0, 0.0, sm / mt ]
        for i in [0, 1, 2]:
            t = qm3.maths.matrix.dot_product( obj.prjT[i], obj.velo )
            for j in range( obj.size ):
                obj.velo[j] -= t * obj.prjT[i][j]
    T, Kin = current_temperature( obj, project )
    scf = math.sqrt( temperature / T )
    for i in range( obj.size ):
        obj.velo[i] *= scf



class velocity_verlet( object ):
    def __init__( self, obj, step_size = 0.001,
                            temperature = 300.0,
                            scale_frequency = 100, 
                            print_frequency = 100,
                            project = True,
                            step_number = -1,
                            temperature_coupling = 0.1,
                            log_function = default_log ):
        self.obj = obj
        self.temperature = temperature
        self.scale_frequency = scale_frequency
        self.print_frequency = print_frequency
        self.project = project
        self.log_function = log_function

        if( scale_frequency > 0 ):
            self.log_function( "---------------------------------------- Dynamics: Velocity-Verlet (NVT)\n" )
            self.log_function( "Temperature:        %20.10lg (K)"%( temperature ) )
            self.log_function( "Scale Frequency:    %20d"%( scale_frequency ) )
            self.integrate = self.integrate_scale
        elif( temperature_coupling > 0.0 ):
            self.log_function( "---------------------------------------- Dynamics: Velocity-Verlet (NVT, Berendsen)\n" )
            self.log_function( "Temperature:        %20.10lg (K)"%( temperature ) )
            self.log_function( "Temp. coupling:     %20.10lg (ps)"%( temperature_coupling ) )
            self.tc = step_size / temperature_coupling
            self.integrate = self.integrate_berendsen
        else:
            self.log_function( "---------------------------------------- Dynamics: Velocity-Verlet (NVE)\n" )
            self.integrate = self.integrate_nve
        self.log_function( "Step Size:          %20.10lg (ps)"%( step_size ) )
        if( step_number > 0 ):
            self.log_function( "Step Number:        %20d"%( step_number ) )
        self.log_function( "Print Frequency:    %20d"%( print_frequency ) )
        self.log_function( "\n%20s%20s%20s%20s%20s"%( "Time (ps)", "Potential (kJ/mol)", "Kinetic (kJ/mol)", "Total (kJ/mol)", "Temperature (K)" ) )
        self.log_function( 100 * "-" )
        self.fc = step_size
        self.fv = self.fc * 0.5
        self.fa = self.fc * self.fv
        self.xavr = [ 0.0, 0.0, 0.0, 0.0 ]
        self.xrms = [ 0.0, 0.0, 0.0, 0.0 ]
        self.time = 0.0
        self.istp = 0
        self.cacc = []
        self.oacc = []
        self.obj.get_grad()
        for i in range( self.obj.size // 3 ):
            i3 = i * 3
            for j in [0, 1, 2]:
                self.cacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
                self.oacc.append( 0.0 )
        if( self.project ):
#            qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.cacc, [] )
            for i in [0, 1, 2]:
                t = qm3.maths.matrix.dot_product( self.obj.prjT[i], self.cacc )
                for j in range( self.obj.size ):
                    self.cacc[j] -= t * self.obj.prjT[i][j]
        self.T, self.Kin = current_temperature( self.obj, self.project )
        log_function( "%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf"%( 0.0, self.obj.func, self.Kin, self.obj.func + self.Kin, self.T ) )
        self.xavr[0] += self.obj.func
        self.xavr[1] += self.Kin
        self.xavr[2] += self.obj.func + self.Kin
        self.xavr[3] += self.T
        self.xrms[0] += self.obj.func * self.obj.func
        self.xrms[1] += self.Kin * self.Kin
        self.xrms[2] += ( self.obj.func + self.Kin ) * ( self.obj.func + self.Kin )
        self.xrms[3] += self.T * self.T
        if( step_number > 0 ):
            for i in range( step_number ):
                self.integrate()
            self.stats()


    def integrate_nve( self ):
        self.time += self.fc
        self.istp += 1
        for i in range( self.obj.size ):
            self.obj.coor[i] += self.fc * self.obj.velo[i] + self.fa * self.cacc[i]
        self.obj.get_grad()
        for i in range( self.obj.size // 3 ):
            i3 = i * 3
            for j in [0, 1, 2]:
                self.oacc[i3+j] = self.cacc[i3+j]
                self.cacc[i3+j] = - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0
        if( self.project ):
#            qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.cacc, [] )
            for i in [0, 1, 2]:
                t = qm3.maths.matrix.dot_product( self.obj.prjT[i], self.cacc )
                for j in range( self.obj.size ):
                    self.cacc[j] -= t * self.obj.prjT[i][j]
        for i in range( self.obj.size ):
            self.obj.velo[i] += self.fv * ( self.oacc[i] + self.cacc[i] )
        self.T, self.Kin = current_temperature( self.obj, self.project )
        self.xavr[0] += self.obj.func
        self.xavr[1] += self.Kin
        self.xavr[2] += self.obj.func + self.Kin
        self.xavr[3] += self.T
        self.xrms[0] += self.obj.func * self.obj.func
        self.xrms[1] += self.Kin * self.Kin
        self.xrms[2] += ( self.obj.func + self.Kin ) * ( self.obj.func + self.Kin )
        self.xrms[3] += self.T * self.T
        if( self.istp % self.print_frequency == 0 ):
            self.log_function( "%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf"%( self.time, self.obj.func, self.Kin, self.obj.func + self.Kin, self.T ) )
        self.obj.current_step( self.istp )


    def integrate_scale( self ):
        self.time += self.fc
        self.istp += 1
        for i in range( self.obj.size ):
            self.obj.coor[i] += self.fc * self.obj.velo[i] + self.fa * self.cacc[i]
        self.obj.get_grad()
        for i in range( self.obj.size // 3 ):
            i3 = i * 3
            for j in [0, 1, 2]:
                self.oacc[i3+j] = self.cacc[i3+j]
                self.cacc[i3+j] = - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0
        if( self.project ):
#            qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.cacc, [] )
            for i in [0, 1, 2]:
                t = qm3.maths.matrix.dot_product( self.obj.prjT[i], self.cacc )
                for j in range( self.obj.size ):
                    self.cacc[j] -= t * self.obj.prjT[i][j]
        for i in range( self.obj.size ):
            self.obj.velo[i] += self.fv * ( self.oacc[i] + self.cacc[i] )
        self.T, self.Kin = current_temperature( self.obj, self.project )
        if( self.istp%self.scale_frequency == 0 ): 
            scf = math.sqrt( self.temperature / self.T )
            for i in range( self.obj.size ):
                self.obj.velo[i] *= scf
            self.T, self.Kin = current_temperature( self.obj, self.project )
        self.xavr[0] += self.obj.func
        self.xavr[1] += self.Kin
        self.xavr[2] += self.obj.func + self.Kin
        self.xavr[3] += self.T
        self.xrms[0] += self.obj.func * self.obj.func
        self.xrms[1] += self.Kin * self.Kin
        self.xrms[2] += ( self.obj.func + self.Kin ) * ( self.obj.func + self.Kin )
        self.xrms[3] += self.T * self.T
        if( self.istp % self.print_frequency == 0 ):
            self.log_function( "%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf"%( self.time, self.obj.func, self.Kin, self.obj.func + self.Kin, self.T ) )
        self.obj.current_step( self.istp )


    def integrate_berendsen( self ):
        self.time += self.fc
        self.istp += 1
        for i in range( self.obj.size ):
            self.obj.coor[i] += self.fc * self.obj.velo[i] + self.fa * self.cacc[i]
        self.obj.get_grad()
        for i in range( self.obj.size // 3 ):
            i3 = i * 3
            for j in [0, 1, 2]:
                self.oacc[i3+j] = self.cacc[i3+j]
                self.cacc[i3+j] = - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0
        if( self.project ):
#            qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.cacc, [] )
            for i in [0, 1, 2]:
                t = qm3.maths.matrix.dot_product( self.obj.prjT[i], self.cacc )
                for j in range( self.obj.size ):
                    self.cacc[j] -= t * self.obj.prjT[i][j]
        for i in range( self.obj.size ):
            self.obj.velo[i] += self.fv * ( self.oacc[i] + self.cacc[i] )
        self.T, self.Kin = current_temperature( self.obj, self.project )
        scv = math.sqrt( 1.0 + self.tc * ( self.temperature / self.T - 1.0 ) )
        scv = min( max( scv, 0.9 ), 1.1 )
        for i in range( self.obj.size ):
            self.obj.velo[i] *= scv
        self.T, self.Kin = current_temperature( self.obj, self.project )
        self.xavr[0] += self.obj.func
        self.xavr[1] += self.Kin
        self.xavr[2] += self.obj.func + self.Kin
        self.xavr[3] += self.T
        self.xrms[0] += self.obj.func * self.obj.func
        self.xrms[1] += self.Kin * self.Kin
        self.xrms[2] += ( self.obj.func + self.Kin ) * ( self.obj.func + self.Kin )
        self.xrms[3] += self.T * self.T
        if( self.istp % self.print_frequency == 0 ):
            self.log_function( "%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf"%( self.time, self.obj.func, self.Kin, self.obj.func + self.Kin, self.T ) )
        self.obj.current_step( self.istp )


    def stats( self ):
        self.log_function( 100 * "-" )
        savr = "%-20s"%( "Averages:" )
        srms = "%-20s"%( "RMS Deviations:" )
        for i in range( len( self.xavr ) ):
            tavr  = self.xavr[i] / float( self.istp + 1 )
            savr += "%20.5lf"%( tavr )
            trms  = self.xrms[i] / float( self.istp + 1 )
            srms += "%20.5lf"%( math.sqrt( math.fabs( trms - tavr * tavr ) ) )
        self.log_function( savr )
        self.log_function( srms )
        self.log_function( 100 * "-" + "\n" )



class langevin_verlet( object ):
    def __init__( self, obj, step_size = 0.001,
                            temperature = 300.0,
                            gamma_factor = 50.0, 
                            print_frequency = 100,
                            project = True,
                            step_number = -1,
                            log_function = default_log ):
        self.obj = obj
        self.step_size = step_size
        self.temperature = temperature
        self.print_frequency = print_frequency
        self.project = project
        self.log_function = log_function

        self.log_function( "---------------------------------------- Dynamics: Langevin-Verlet (NVT)\n" )
        self.log_function( "Step Size:          %20.10lg (ps)"%( step_size ) )
        self.log_function( "Temperature:        %20.10lg (K)"%( temperature ) )
        self.log_function( "Gamma Factor:       %20.10lg (ps^-1)"%( gamma_factor ) )
        if( step_number > 0 ):
            self.log_function( "Step Number:        %20d"%( step_number ) )
        self.log_function( "Print Frequency:    %20d"%( print_frequency ) )
        ff  = step_size * gamma_factor
        if( ff < 0.01 ):
            ff = 0.01
            self.log_function( "\n>> Gamma factor:    %20.10lg (ps^-1)"%( 0.01 / step_size ) )
        self.log_function( "\n%20s%20s%20s%20s%20s"%( "Time (ps)", "Potential (kJ/mol)", "Kinetic (kJ/mol)", "Total (kJ/mol)", "Temperature (K)" ) )
        self.log_function( 100 * "-" )
        self.c0  = math.exp( - ff )
        self.c1  = ( 1.0 - self.c0 ) / ff
        self.c2  = ( 1.0 - self.c1 ) / ff
        self.sr  = self.step_size * math.sqrt( ( 2.0 - ( 3.0 - 4.0 * self.c0 + self.c0 * self.c0 ) / ff ) / ff )
        self.sv  = math.sqrt( 1.0 - self.c0 * self.c0 )
        self.cv1 = self.step_size * ( 1.0 - self.c0 ) * ( 1.0 - self.c0 ) / ( ff * self.sr * self.sv )
        self.cv2 = math.sqrt( 1.0 - self.cv1 * self.cv1 )
        self.ff  = [ 1.0e-2 * math.sqrt( qm3.constants.KB * self.temperature * _Gfac / ii ) for ii in self.obj.mass ]
        self.fr1 = self.step_size * self.c1
        self.fv1 = self.step_size * ( self.c1 - self.c2 )
        self.fv2 = self.step_size * self.c2
        self.fr2 = self.step_size * self.fv2
        self.xavr = [ 0.0, 0.0, 0.0, 0.0 ]
        self.xrms = [ 0.0, 0.0, 0.0, 0.0 ]
        self.time = 0.0
        self.istp = 0
        self.oacc = []
        self.cacc = []
        self.obj.get_grad()
        for i in range( self.obj.size // 3 ):
            i3 = i * 3
            for j in [0, 1, 2]:
                self.cacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
                self.oacc.append( 0.0 )
        if( self.project ):
#            qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.cacc, [] )
            for i in [0, 1, 2]:
                t = qm3.maths.matrix.dot_product( self.obj.prjT[i], self.cacc )
                for j in range( self.obj.size ):
                    self.cacc[j] -= t * self.obj.prjT[i][j]
        self.T, self.Kin = current_temperature( self.obj, self.project )
        log_function( "%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf"%( 0.0, self.obj.func, self.Kin, self.obj.func + self.Kin, self.T ) )
        self.xavr[0] += self.obj.func
        self.xavr[1] += self.Kin
        self.xavr[2] += self.obj.func + self.Kin
        self.xavr[3] += self.T
        self.xrms[0] += self.obj.func * self.obj.func
        self.xrms[1] += self.Kin * self.Kin
        self.xrms[2] += ( self.obj.func + self.Kin ) * ( self.obj.func + self.Kin )
        self.xrms[3] += self.T * self.T
        if( step_number > 0 ):
            for i in range( step_number ):
                self.integrate()
            self.stats()


    def integrate( self ):
        self.time += self.step_size
        self.istp += 1
        for i in range( self.obj.size ):
            self.obj.coor[i] += self.fr1 * self.obj.velo[i] + self.fr2 * self.cacc[i]
            # random forces
            r1 = qm3.maths.rand.gauss( 0.0, 1.0 )
            r2 = qm3.maths.rand.gauss( 0.0, 1.0 )
            self.obj.coor[i] += self.ff[i//3] * self.sr * r1
            self.oacc[i] = self.c0 * self.obj.velo[i] + self.fv1 * self.cacc[i] + self.ff[i//3] * self.sv * ( self.cv1 * r1 + self.cv2 * r2 )
        self.obj.get_grad()
        for i in range( self.obj.size // 3 ):
            i3 = i * 3
            for j in [0, 1, 2]:
                self.cacc[i3+j] = - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0
        if( self.project ):
#            RT = qm3.utils.get_RT_modes( self.obj.mass, self.obj.coor )
#            for i in range( 6 ):
#                t = sum( [ ii * jj  for ii,jj in zip( self.cacc, RT[i] ) ] )
#                for j in range( self.obj.size ):
#                    self.cacc[j] -= t * RT[i][j]
            for i in [0, 1, 2]:
                t = qm3.maths.matrix.dot_product( self.obj.prjT[i], self.cacc )
                for j in range( self.obj.size ):
                    self.cacc[j] -= t * self.obj.prjT[i][j]
        for i in range( self.obj.size ):
            self.obj.velo[i] = self.oacc[i] + self.fv2 * self.cacc[i]
        if( self.project ):
#            for i in range( 6 ):
#                t = sum( [ ii * jj  for ii,jj in zip( self.obj.velo, RT[i] ) ] )
#                for j in range( self.obj.size ):
#                    self.obj.velo[j] -= t * RT[i][j]
            for i in [0, 1, 2]:
                t = qm3.maths.matrix.dot_product( self.obj.prjT[i], self.obj.velo )
                for j in range( self.obj.size ):
                    self.obj.velo[j] -= t * self.obj.prjT[i][j]
        self.T, self.Kin = current_temperature( self.obj, self.project )
        self.xavr[0] += self.obj.func
        self.xavr[1] += self.Kin
        self.xavr[2] += self.obj.func + self.Kin
        self.xavr[3] += self.T
        self.xrms[0] += self.obj.func * self.obj.func
        self.xrms[1] += self.Kin * self.Kin
        self.xrms[2] += ( self.obj.func + self.Kin ) * ( self.obj.func + self.Kin )
        self.xrms[3] += self.T * self.T
        if( self.istp % self.print_frequency == 0 ):
            self.log_function( "%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf"%( self.time, self.obj.func, self.Kin, self.obj.func + self.Kin, self.T ) )
        self.obj.current_step( self.istp )


    def stats( self ):
        self.log_function( 100 * "-" )
        savr = "%-20s"%( "Averages:" )
        srms = "%-20s"%( "RMS Deviations:" )
        for i in range( len( self.xavr ) ):
            tavr  = self.xavr[i] / float( self.istp + 1 )
            savr += "%20.5lf"%( tavr )
            trms  = self.xrms[i] / float( self.istp + 1 )
            srms += "%20.5lf"%( math.sqrt( math.fabs( trms - tavr * tavr ) ) )
        self.log_function( savr )
        self.log_function( srms )
        self.log_function( 100 * "-" + "\n" )



# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.utils
import qm3.actions.dynamics



def default_log( txt ):
    sys.stdout.write( txt + "\n" )
    sys.stdout.flush()



class grote_hynes( object ):
    def __init__( self, obj, coordinate,
                            step_size = 0.0005,
                            print_frequency = 100,
                            project_RT = True,
                            step_number = -1,
                            log_function = default_log ):
        self.obj = obj
        self.print_frequency = print_frequency
        self.project_RT = project_RT
        self.log_function = log_function
        self.coordinate = coordinate
        self.reference = self.coordinate.value( obj.coor )
        self.log_function( "---------------------------------------- Dynamics: Grote-Hynes (NVE)\n" )
        self.log_function( "Step Size:          %20.10lg (ps)"%( step_size ) )
        if( step_number > 0 ):
            self.log_function( "Step Number:        %20d"%( step_number ) )
        self.log_function( "Print Frequency:    %20d"%( print_frequency ) )
        self.log_function( "Coordinate ref.:    %20.10lg"%( self.reference ) )
        self.log_function( "\n%20s%20s%20s%20s%20s"%( "Time (ps)", "Potential (kJ/mol)", "Kinetic (kJ/mol)", "Total (kJ/mol)", "Temperature (K)" ) )
        self.log_function( 100 * "-" )
        self.fc = step_size
        self.fv = self.fc * 0.5
        self.fa = self.fc * self.fv
        self.xavr = [ 0.0, 0.0, 0.0, 0.0 ]
        self.xrms = [ 0.0, 0.0, 0.0, 0.0 ]
        self.time = 0.0
        self.istp = 0
        self.vacc = []
        self.obj.get_grad()
        if( self.project_RT ):
            qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.obj.grad, None )
        for i in range( self.obj.size // 3 ):
            i3 = i * 3
            for j in [0, 1, 2]:
                self.vacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
        self.coordinate.force( self.obj.mass, self.obj.coor, self.obj.grad )
        self.T, self.Kin = qm3.actions.dynamics.current_temperature( self.obj, self.project_RT )
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
            self.coordinate.stop()


    def integrate( self ):
        self.time += self.fc
        self.istp += 1
        # constraint positions
        self.coordinate.constraint( self.obj.mass, self.obj.coor, self.fc, self.obj.velo,
            self.fa, self.vacc, 100.0 * self.fa, self.reference )
        for i in range( self.obj.size ):
            self.obj.coor[i] += self.fc * self.obj.velo[i] + self.fa * self.vacc[i]
        self.obj.get_grad()
        if( self.project_RT ):
            qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.obj.grad, None )
        self.coordinate.force( self.obj.mass, self.obj.coor, self.obj.grad )
        vtmp = self.vacc[:]
        self.vacc = []
        for i in range( self.obj.size // 3 ):
            i3 = i * 3
            for j in [0, 1, 2]:
                self.vacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
        self.coordinate.constraint( self.obj.mass, self.obj.coor, self.fc, self.obj.velo,
            self.fa, self.vacc, 100.0 * self.fc * self.fc, self.reference )
        for i in range( self.obj.size ):
            self.obj.velo[i] += self.fv * ( vtmp[i] + self.vacc[i] )
        self.T, self.Kin = qm3.actions.dynamics.current_temperature( self.obj, self.project_RT )
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



#
# Coordinates
#
class coordinate_antisymmetric( object ):
    """
    [ DONOR, TRANSFERRED, ACCEPTOR ]
    """

    def __init__( self, atoms ):
        self.atoms = atoms[:]
        self.data  = []
        self.__fd  = open( "GH_antisymmetric.log", "wt" )


    def stop( self ):
        self.__fd.close()


    def value( self, coor ):
        return( qm3.utils.distance( coor[3*self.atoms[0]:3*self.atoms[0]+3], coor[3*self.atoms[1]:3*self.atoms[1]+3] ) -
            qm3.utils.distance( coor[3*self.atoms[1]:3*self.atoms[1]+3], coor[3*self.atoms[2]:3*self.atoms[2]+3] ) )


    def force( self, mass, coor, grad ):
        d1 = qm3.utils.distance( coor[3*self.atoms[0]:3*self.atoms[0]+3], coor[3*self.atoms[1]:3*self.atoms[1]+3] )
        d2 = qm3.utils.distance( coor[3*self.atoms[1]:3*self.atoms[1]+3], coor[3*self.atoms[2]:3*self.atoms[2]+3] )
        B  = []
        for i in [0, 1, 2]:
            B.append( ( coor[3*self.atoms[0]+i] - coor[3*self.atoms[1]+i] ) / d1 )
        for i in [0, 1, 2]:
            B.append( - ( coor[3*self.atoms[0]+i] - coor[3*self.atoms[1]+i] ) / d1 
                        - ( coor[3*self.atoms[2]+i] - coor[3*self.atoms[1]+i] ) / d2 )
        for i in [0, 1, 2]:
            B.append( ( coor[3*self.atoms[2]+i] - coor[3*self.atoms[1]+i] ) / d2 )
        rms = 1.0 / ( B[0] * B[0] / mass[self.atoms[0]] + B[1] * B[1] / mass[self.atoms[0]] + B[2] * B[2] / mass[self.atoms[0]] +
                B[3] * B[3] / mass[self.atoms[1]] + B[4] * B[4] / mass[self.atoms[1]] + B[5] * B[5] / mass[self.atoms[1]] +
                B[6] * B[6] / mass[self.atoms[2]] + B[7] * B[7] / mass[self.atoms[2]] + B[8] * B[8] / mass[self.atoms[2]]  )
        frz = - rms * ( B[0] * grad[3*self.atoms[0]] / mass[self.atoms[0]] + 
            B[1] * grad[3*self.atoms[0]+1] / mass[self.atoms[0]] +
            B[2] * grad[3*self.atoms[0]+2] / mass[self.atoms[0]] +
            B[3] * grad[3*self.atoms[1]]   / mass[self.atoms[1]] + 
            B[4] * grad[3*self.atoms[1]+1] / mass[self.atoms[1]] +
            B[5] * grad[3*self.atoms[1]+2] / mass[self.atoms[1]] +
            B[6] * grad[3*self.atoms[2]]   / mass[self.atoms[2]] + 
            B[7] * grad[3*self.atoms[2]+1] / mass[self.atoms[2]] +
            B[8] * grad[3*self.atoms[2]+2] / mass[self.atoms[2]] )
        self.__fd.write( "%20.10lf%20.10lf%20.10lf\n"%( d1 - d2, frz, rms ) )
        self.__fd.flush()
        self.data.append( frz )


    def constraint( self, mass, coor, fvel, velo, facc, acce, flmb, refc ):
        m = [ mass[self.atoms[0]], mass[self.atoms[1]], mass[self.atoms[2]] ]
        r = []
        v = []
        a = []
        for i in self.atoms:
            i3 = 3 * i
            for j in [0, 1, 2]:
                r.append( coor[i3+j] )
                v.append( velo[i3+j] )
                a.append( acce[i3+j] )
        qq = True
        it = 0
        ll = 1.0
        ff = self.__lambda_func( m, r, fvel, v, facc, a, flmb, ll, refc )
        while( qq ):
# ------------------------------------------------------------------------------------------------------
            # numerical
#            gg = ( self.__lambda_func( m, r, fvel, v, facc, a, flmb, ll + 1.e-6, refc ) -
#                self.__lambda_func( m, r, fvel, v, facc, a, flmb, ll - 1.e-6, refc ) ) / 2.e-6
            # analytical
            gg = self.__lambda_grad( m, r, fvel, v, facc, a, flmb, ll, refc )
# ------------------------------------------------------------------------------------------------------
            dd = ff / gg
            if( math.fabs( dd ) > 100.0 ):
                dd = math.copysign( 100.0, dd )
            ll -= dd
            ff = self.__lambda_func( m, r, fvel, v, facc, a, flmb, ll, refc )
            it += 1
            qq = it < 10000 and math.fabs( ff ) > 1.e-14
        r1 = [ coor[3*self.atoms[0]+i] - coor[3*self.atoms[1]+i] for i in [0, 1, 2] ]
        m1 = math.sqrt( sum( [ i*i for i in r1 ] ) )
        r2 = [ coor[3*self.atoms[2]+i] - coor[3*self.atoms[1]+i] for i in [0, 1, 2] ]
        m2 = math.sqrt( sum( [ i*i for i in r1 ] ) )
        for i in [0, 1, 2]:
            acce[3*self.atoms[0]+i] -= 100.0 * ll / mass[self.atoms[0]] * r1[i] / m1
            acce[3*self.atoms[2]+i] += 100.0 * ll / mass[self.atoms[2]] * r2[i] / m2
            acce[3*self.atoms[1]+i] += 100.0 * ll / mass[self.atoms[1]] * ( r1[i] / m1 - r2[i] / m2 )


    def __lambda_func( self, m, r, fv, v, fa, a, fl, lmb, c ):
        return( -c + math.sqrt(math.pow(a[0]*fa - a[3]*fa + r[0] - r[3] - (fl*lmb*(r[0] - r[3])) /
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (fl*lmb*(-((r[0] - r[3]) / math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (-r[3] + r[6]) / math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] + fv*v[0] - fv*v[3],2) + math.pow(a[1]*fa - a[4]*fa + r[1] - r[4] - 
            (fl*lmb*(r[1] - r[4])) / (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (fl*lmb*(-((r[1] - r[4]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + (-r[4] + r[7]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] + fv*v[1] - fv*v[4],2) + math.pow(a[2]*fa - a[5]*fa + r[2] - 
            (fl*lmb*(r[2] - r[5])) / (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - r[5] + (fl*lmb*(-((r[2] - r[5]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + (-r[5] + r[8]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] + fv*v[2] - fv*v[5],2))
            - math.sqrt(math.pow(-(a[3]*fa) + a[6]*fa - r[3] + r[6] + (fl*lmb*(-r[3] + r[6])) /
            (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2))) + 
            (fl*lmb*(-((r[0] - r[3]) / math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (-r[3] + r[6]) / math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] - fv*v[3] + fv*v[6],2) + math.pow(-(a[4]*fa) + a[7]*fa - r[4] + r[7] + 
            (fl*lmb*(-r[4] + r[7])) / (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))) + (fl*lmb*(-((r[1] - r[4]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + (-r[4] + r[7]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] - fv*v[4] + fv*v[7],2) + math.pow(-(a[5]*fa) + a[8]*fa - r[5] + r[8] + 
            (fl*lmb*(-r[5] + r[8])) / (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))) + (fl*lmb*(-((r[2] - r[5]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[5] + r[8]) / math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] - fv*v[5] + fv*v[8],2)) )


    def __lambda_grad( self, m, r, fv, v, fa, a, fl, lmb, c ):
        return( (2*(-((fl*(r[0] - r[3])) / (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2)))) + (fl*(-((r[0] - r[3]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[3] + r[6]) / math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1])* (a[0]*fa - a[3]*fa + r[0] - r[3] - 
            (fl*lmb*(r[0] - r[3])) / (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (fl*lmb*(-((r[0] - r[3]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[3] + r[6]) / math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] + fv*v[0] - fv*v[3]) + 2*(-((fl*(r[1] - r[4])) /
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2)))) + (fl*(-((r[1] - r[4]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (-r[4] + r[7]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1])* (a[1]*fa - a[4]*fa + r[1] - r[4] - 
            (fl*lmb*(r[1] - r[4])) / (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (fl*lmb*(-((r[1] - r[4]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (-r[4] + r[7]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] + fv*v[1] - fv*v[4]) + 2*(-((fl*(r[2] - r[5])) /
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2)))) + (fl*(-((r[2] - r[5]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (-r[5] + r[8]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1])* (a[2]*fa - a[5]*fa + r[2] - 
            (fl*lmb*(r[2] - r[5])) / (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - r[5] + (fl*lmb*(-((r[2] - r[5]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[5] + r[8]) / math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] + fv*v[2] - fv*v[5])) /
            (2.*math.sqrt(math.pow(a[0]*fa - a[3]*fa + r[0] - r[3] - 
            (fl*lmb*(r[0] - r[3])) / (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (fl*lmb*(-((r[0] - r[3]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[3] + r[6]) / math.sqrt(math.pow(-r[3] + r[6],2) + 
            math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2)))) / m[1] + fv*v[0] - fv*v[3],2) + 
            math.pow(a[1]*fa - a[4]*fa + r[1] - r[4] - (fl*lmb*(r[1] - r[4])) /
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (fl*lmb*(-((r[1] - r[4]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[4] + r[7]) / math.sqrt(math.pow(-r[3] + r[6],2) + 
            math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2)))) / m[1] + fv*v[1] - fv*v[4],2) + 
            math.pow(a[2]*fa - a[5]*fa + r[2] - (fl*lmb*(r[2] - r[5])) /
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - r[5] + (fl*lmb*(-((r[2] - r[5]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[5] + r[8]) / math.sqrt(math.pow(-r[3] + r[6],2) + 
            math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2)))) / m[1] + fv*v[2] - fv*v[5],2))) - 
            (2*((fl*(-r[3] + r[6])) / (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))) + (fl*(-((r[0] - r[3]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[3] + r[6]) / math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1])* (-(a[3]*fa) + a[6]*fa - r[3] + r[6] + 
            (fl*lmb*(-r[3] + r[6])) / (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))) + (fl*lmb*(-((r[0] - r[3]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[3] + r[6]) / math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] - fv*v[3] + fv*v[6]) + 2*((fl*(-r[4] + r[7])) /
            (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))) + (fl*(-((r[1] - r[4]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[4] + r[7]) / math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1])* (-(a[4]*fa) + a[7]*fa - r[4] + r[7] + 
            (fl*lmb*(-r[4] + r[7])) / (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))) + (fl*lmb*(-((r[1] - r[4]) /
            math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + 
            (-r[4] + r[7]) / math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] - fv*v[4] + fv*v[7]) + 2*((fl*(-r[5] + r[8])) /
            (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2))) + 
            (fl*(-((r[2] - r[5]) / math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (-r[5] + r[8]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2))))/m[1])*
            (-(a[5]*fa) + a[8]*fa - r[5] + r[8] + (fl*lmb*(-r[5] + r[8])) /
            (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2))) + 
            (fl*lmb*(-((r[2] - r[5]) / math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + (-r[5] + r[8]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))))/m[1] - fv*v[5] + fv*v[8])) /
            (2.*math.sqrt(math.pow(-(a[3]*fa) + a[6]*fa - r[3] + r[6] + (fl*lmb*(-r[3] + r[6])) /
            (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2))) + 
            (fl*lmb*(-((r[0] - r[3]) / math.sqrt(math.pow(r[0] - r[3],2) + 
            math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + (-r[3] + r[6]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2)))) /
            m[1] - fv*v[3] + fv*v[6],2) + math.pow(-(a[4]*fa) + a[7]*fa - r[4] + r[7] + 
            (fl*lmb*(-r[4] + r[7])) / (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))) + (fl*lmb*(-((r[1] - r[4]) / math.sqrt(math.pow(r[0] - r[3],2) + 
            math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + (-r[4] + r[7]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2)))) /
            m[1] - fv*v[4] + fv*v[7],2) + math.pow(-(a[5]*fa) + a[8]*fa - r[5] + r[8] + 
            (fl*lmb*(-r[5] + r[8])) / (m[2]*math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + 
            math.pow(-r[5] + r[8],2))) + (fl*lmb*(-((r[2] - r[5]) / math.sqrt(math.pow(r[0] - r[3],2) + 
            math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) + (-r[5] + r[8]) /
            math.sqrt(math.pow(-r[3] + r[6],2) + math.pow(-r[4] + r[7],2) + math.pow(-r[5] + r[8],2)))) /
            m[1] - fv*v[5] + fv*v[8],2))) )




class coordinate_bond( object ):
    """
    f[r_, d_] := Sqrt[(r[[1]] - r[[4]])^2 + (r[[2]] - r[[5]])^2 + (r[[3]] - r[[6]])^2] - d
    rdt[mi_, r_, fv_, v_, fa_, a_, fl_, lmb_] := r + fv*v + fa*a - fl/mi*lmb*{ D[f[r, d], r[[1]]], D[f[r, d], r[[2]]], D[f[r, d], r[[3]]], D[f[r, d], r[[4]]], D[f[r, d], r[[5]]], D[f[r, d], r[[6]]] }
    CForm[f[rdt[{m1, m1, m1, m2, m2, m2}, {r1, r2, r3, r4, r5, r6}, fv, {v1, v2, v3, v4, v5, v6}, fa, {a1, a2, a3, a4, a5, a6}, fl, lmb], c]]
    CForm[D[f[ rdt[{m1, m1, m1, m2, m2, m2}, {r1, r2, r3, r4, r5, r6}, fv, {v1, v2, v3, v4, v5, v6}, fa, {a1, a2, a3, a4, a5, a6}, fl, lmb], c], lmb]]
    """

    def __init__( self, atoms ):
        self.atoms = atoms[:]
        self.data  = []
        self.__fd  = open( "GH_bond.log", "wt" )


    def stop( self ):
        self.__fd.close()


    def value( self, coor ):
        return( qm3.utils.distance( coor[3*self.atoms[0]:3*self.atoms[0]+3], coor[3*self.atoms[1]:3*self.atoms[1]+3] ) )


    def force( self, mass, coor, grad ):
        d = [ coor[3*self.atoms[0]+i] - coor[3*self.atoms[1]+i] for i in [0, 1, 2] ]
        m = math.sqrt( sum( [ i*i for i in d ] ) )
        B = [ d[0] / m, d[1] / m, d[2] / m, - d[0] / m, - d[1] / m, - d[2] / m ]
        rms = 1.0 / ( B[0] * B[0] / mass[self.atoms[0]] + B[1] * B[1] / mass[self.atoms[0]] + B[2] * B[2] / mass[self.atoms[0]] +
                B[3] * B[3] / mass[self.atoms[1]] + B[4] * B[4] / mass[self.atoms[1]] + B[5] * B[5] / mass[self.atoms[1]] )
        frz = - rms * ( B[0] * grad[3*self.atoms[0]] / mass[self.atoms[0]] + 
            B[1] * grad[3*self.atoms[0]+1] / mass[self.atoms[0]] +
            B[2] * grad[3*self.atoms[0]+2] / mass[self.atoms[0]] +
            B[3] * grad[3*self.atoms[1]]   / mass[self.atoms[1]] + 
            B[4] * grad[3*self.atoms[1]+1] / mass[self.atoms[1]] +
            B[5] * grad[3*self.atoms[1]+2] / mass[self.atoms[1]] )
        self.__fd.write( "%20.10lf%20.10lf%20.10lf\n"%( m, frz, rms ) )
        self.__fd.flush()
        self.data.append( frz )


    def constraint( self, mass, coor, fvel, velo, facc, acce, flmb, refc ):
        m = [ mass[self.atoms[0]], mass[self.atoms[1]] ]
        r = []
        v = []
        a = []
        for i in self.atoms:
            i3 = 3 * i
            for j in [0, 1, 2]:
                r.append( coor[i3+j] )
                v.append( velo[i3+j] )
                a.append( acce[i3+j] )
        qq = True
        it = 0
        ll = 1.0
        ff = self.__lambda_func( m, r, fvel, v, facc, a, flmb, ll, refc )
        while( qq ):
# ------------------------------------------------------------------------------------------------------
            # numerical
#            gg = ( self.__lambda_func( m, r, fvel, v, facc, a, flmb, ll + 1.e-6, refc ) -
#                self.__lambda_func( m, r, fvel, v, facc, a, flmb, ll - 1.e-6, refc ) ) / 2.e-6
            # analytical
            gg = self.__lambda_grad( m, r, fvel, v, facc, a, flmb, ll, refc )
# ------------------------------------------------------------------------------------------------------
            dd = ff / gg
            if( math.fabs( dd ) > 100.0 ):
                dd = math.copysign( 100.0, dd )
            ll -= dd
            ff = self.__lambda_func( m, r, fvel, v, facc, a, flmb, ll, refc )
            it += 1
            qq = it < 10000 and math.fabs( ff ) > 1.e-14
        print( it, ff, gg, ll )
        rr = [ coor[3*self.atoms[0]+i] - coor[3*self.atoms[1]+i] for i in [0, 1, 2] ]
        mm = math.sqrt( sum( [ i*i for i in rr ] ) )
        for i in [0, 1, 2]:
            acce[3*self.atoms[0]+i] -= 100.0 * ll / mass[self.atoms[0]] * rr[i] / mm
            acce[3*self.atoms[1]+i] += 100.0 * ll / mass[self.atoms[1]] * rr[i] / mm


    def __lambda_func( self, m, r, fv, v, fa, a, fl, lmb, c ):
        return( -c + math.sqrt(math.pow(a[0]*fa - a[3]*fa + r[0] - r[3] - (fl*lmb*(r[0] - r[3]))/
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))) - 
            (fl*lmb*(r[0] - r[3]))/ (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + fv*v[0] - fv*v[3],2) + math.pow(a[1]*fa - a[4]*fa + r[1] - r[4] - 
            (fl*lmb*(r[1] - r[4]))/ (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - (fl*lmb*(r[1] - r[4]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + fv*v[1] - fv*v[4],2) + math.pow(a[2]*fa - a[5]*fa + r[2] - 
            (fl*lmb*(r[2] - r[5]))/ (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - (fl*lmb*(r[2] - r[5]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - r[5] + fv*v[2] - fv*v[5],2)) )


    def __lambda_grad( self, m, r, fv, v, fa, a, fl, lmb, c ):
        return( (2*(-((fl*(r[0] - r[3]))/ (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2)))) - (fl*(r[0] - r[3]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + math.pow(r[2] - r[5],2))))*
            (a[0]*fa - a[3]*fa + r[0] - r[3] - (fl*lmb*(r[0] - r[3]))/
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - (fl*lmb*(r[0] - r[3]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + fv*v[0] - fv*v[3]) + 2*(-((fl*(r[1] - r[4]))/
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2)))) - (fl*(r[1] - r[4]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))))* (a[1]*fa - a[4]*fa + r[1] - r[4] - 
            (fl*lmb*(r[1] - r[4]))/ (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - (fl*lmb*(r[1] - r[4]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + fv*v[1] - fv*v[4]) + 2*(-((fl*(r[2] - r[5]))/
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2)))) - (fl*(r[2] - r[5]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))))* (a[2]*fa - a[5]*fa + r[2] - (fl*lmb*(r[2] - r[5]))/
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - (fl*lmb*(r[2] - r[5]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - r[5] + fv*v[2] - fv*v[5]))/
            (2.*math.sqrt(math.pow(a[0]*fa - a[3]*fa + r[0] - r[3] - (fl*lmb*(r[0] - r[3]))/
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - (fl*lmb*(r[0] - r[3]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + fv*v[0] - fv*v[3],2) + 
            math.pow(a[1]*fa - a[4]*fa + r[1] - r[4] - (fl*lmb*(r[1] - r[4]))/
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - (fl*lmb*(r[1] - r[4]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) + fv*v[1] - fv*v[4],2) + 
            math.pow(a[2]*fa - a[5]*fa + r[2] - (fl*lmb*(r[2] - r[5]))/
            (m[0]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - (fl*lmb*(r[2] - r[5]))/
            (m[1]*math.sqrt(math.pow(r[0] - r[3],2) + math.pow(r[1] - r[4],2) + 
            math.pow(r[2] - r[5],2))) - r[5] + fv*v[2] - fv*v[5],2))) )




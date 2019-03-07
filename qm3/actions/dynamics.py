# -*- coding: iso-8859-1 -*-


from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.maths.rand
import	qm3.utils
import	qm3.constants



_Tfac = 10.0 / ( qm3.constants.KB * qm3.constants.NA )
_Gfac = 1000.0 * qm3.constants.NA					# g/mol > kg



def default_log( txt ):
	sys.stdout.write( txt + "\n" )
	sys.stdout.flush()



def current_temperature( obj, project_RT = True ):
	KEf = 0.0
	for i in range( obj.size // 3 ):
		i3  = i * 3
		KEf += obj.mass[i] * ( obj.velo[i3] * obj.velo[i3] + obj.velo[i3+1] * obj.velo[i3+1] + obj.velo[i3+2] * obj.velo[i3+2] )
#	KEf = sum( [ obj.mass[ii] * sum( [ jj*jj for jj in obj.velo[3*ii:3*ii+3] ] ) for ii in range( obj.size // 3 ) ] )
	if( project_RT ):
		T = KEf * _Tfac / float( obj.size - 6.0 )
	else:
		T = KEf * _Tfac / float( obj.size )
	Kin = KEf * 0.005
	return( T, Kin )



def assign_velocities( obj, temperature = 300.0, project_RT = True ):
	obj.velo = []
	KT = qm3.constants.KB * temperature * _Gfac
	for i in range( obj.size // 3 ):
		SD = math.sqrt( KT / obj.mass[i] )
		obj.velo += [ qm3.maths.rand.gauss( 0.0, SD ) * 1.0e-2 for ii in [0, 1, 2] ]	 # ang/ps
	if( project_RT ):
		qm3.utils.project_RT_modes( obj.mass, obj.coor, obj.velo, None )
	T, Kin = current_temperature( obj, project_RT )
	scf = math.sqrt( temperature / T )
	for i in range( obj.size ):
		obj.velo[i] *= scf



class velocity_verlet( object ):
	def __init__( self, obj, step_size = 0.001,
							temperature = 300.0,
							scale_frequency = 100, 
							print_frequency = 100,
							project_RT = True,
							step_number = -1,
							temperature_coupling = 0.1,
							log_function = default_log ):
		self.obj = obj
		self.temperature = temperature
		self.scale_frequency = scale_frequency
		self.print_frequency = print_frequency
		self.project_RT = project_RT
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
		self.vacc = []
		self.obj.get_grad()
		for i in range( self.obj.size // 3 ):
			i3 = i * 3
			for j in [0, 1, 2]:
				self.vacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
		if( self.project_RT ):
			qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.vacc, None )
		self.T, self.Kin = current_temperature( self.obj, self.project_RT )
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
		vtmp = []
		for i in range( self.obj.size ):
			self.obj.coor[i] += self.fc * self.obj.velo[i] + self.fa * self.vacc[i]
		self.obj.get_grad()
		vtmp = self.vacc[:]
		self.vacc = []
		for i in range( self.obj.size // 3 ):
			i3 = i * 3
			for j in [0, 1, 2]:
				self.vacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
		if( self.project_RT ):
			qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.vacc, None )
		for i in range( self.obj.size ):
			self.obj.velo[i] += self.fv * ( vtmp[i] + self.vacc[i] )
		self.T, self.Kin = current_temperature( self.obj, self.project_RT )
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
		vtmp = []
		for i in range( self.obj.size ):
			self.obj.coor[i] += self.fc * self.obj.velo[i] + self.fa * self.vacc[i]
		self.obj.get_grad()
		vtmp = self.vacc[:]
		self.vacc = []
		for i in range( self.obj.size // 3 ):
			i3 = i * 3
			for j in [0, 1, 2]:
				self.vacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
		if( self.project_RT ):
			qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.vacc, None )
		for i in range( self.obj.size ):
			self.obj.velo[i] += self.fv * ( vtmp[i] + self.vacc[i] )
		self.T, self.Kin = current_temperature( self.obj, self.project_RT )
		if( self.istp%self.scale_frequency == 0 ): 
			scf = math.sqrt( self.temperature / self.T )
			for i in range( self.obj.size ):
				self.obj.velo[i] *= scf
			self.T, self.Kin = current_temperature( self.obj, self.project_RT )
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
		vtmp = []
		for i in range( self.obj.size ):
			self.obj.coor[i] += self.fc * self.obj.velo[i] + self.fa * self.vacc[i]
		self.obj.get_grad()
		vtmp = self.vacc[:]
		self.vacc = []
		for i in range( self.obj.size // 3 ):
			i3 = i * 3
			for j in [0, 1, 2]:
				self.vacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
		if( self.project_RT ):
			qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.vacc, None )
		for i in range( self.obj.size ):
			self.obj.velo[i] += self.fv * ( vtmp[i] + self.vacc[i] )
		self.T, self.Kin = current_temperature( self.obj, self.project_RT )
		scv = math.sqrt( 1.0 + self.tc * ( self.temperature / self.T - 1.0 ) )
		scv = min( max( scv, 0.9 ), 1.1 )
		for i in range( self.obj.size ):
			self.obj.velo[i] *= scv
		self.T, self.Kin = current_temperature( self.obj, self.project_RT )
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
							project_RT = True,
							step_number = -1,
							log_function = default_log ):
		self.obj = obj
		self.step_size = step_size
		self.temperature = temperature
		self.print_frequency = print_frequency
		self.project_RT = project_RT
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
		self.vacc = []
		self.obj.get_grad()
		for i in range( self.obj.size // 3 ):
			i3 = i * 3
			for j in [0, 1, 2]:
				self.vacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
		if( self.project_RT ):
			qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.vacc, None )
		self.T, self.Kin = current_temperature( self.obj, self.project_RT )
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
		vtmp = []
		for i in range( self.obj.size ):
			self.obj.coor[i] += self.fr1 * self.obj.velo[i] + self.fr2 * self.vacc[i]
			# random forces
			r1 = qm3.maths.rand.gauss( 0.0, 1.0 )
			r2 = qm3.maths.rand.gauss( 0.0, 1.0 )
			self.obj.coor[i] += self.ff[i//3] * self.sr * r1
			vtmp.append( self.c0 * self.obj.velo[i] + self.fv1 * self.vacc[i] + self.ff[i//3] * self.sv * ( self.cv1 * r1 + self.cv2 * r2 ) )
		self.obj.get_grad()
		self.vacc = []
		for i in range( self.obj.size // 3 ):
			i3 = i * 3
			for j in [0, 1, 2]:
				self.vacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
		if( self.project_RT ):
			RT = qm3.utils.get_RT_modes( self.obj.mass, self.obj.coor )
			for i in range( 6 ):
				t = sum( [ ii * jj  for ii,jj in zip( self.vacc, RT[i] ) ] )
				for j in range( self.obj.size ):
					self.vacc[j] -= t * RT[i][j]
		for i in range( self.obj.size ):
			self.obj.velo[i] = vtmp[i] + self.fv2 * self.vacc[i]
		if( self.project_RT ):
			for i in range( 6 ):
				t = sum( [ ii * jj  for ii,jj in zip( self.obj.velo, RT[i] ) ] )
				for j in range( self.obj.size ):
					self.obj.velo[j] -= t * RT[i][j]
		self.T, self.Kin = current_temperature( self.obj, self.project_RT )
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



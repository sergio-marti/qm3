# -*- coding: iso-8859-1 -*-


from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.actions.dynamics
import	qm3.utils
import	qm3.constants



class leapfrog_verlet( object ):
	def __init__( self, obj, step_size = 0.001,
							temperature = 300.0,
							temperature_coupling = 0.1,
							pressure = 1.0,
							pressure_coupling = 1.0,
							print_frequency = 100,
							project_RT = True,
							step_number = -1,
							log_function = qm3.actions.dynamics.default_log ):
		self.obj = obj
		self.step_size = step_size
		self.temperature = temperature
		self.pressure = pressure
		self.print_frequency = print_frequency
		self.project_RT = project_RT
		self.log_function = log_function

		self.log_function( "---------------------------------------- Dynamics: Leapfrog-Verlet (NpT, Berendsen)\n" )
		self.log_function( "Step Size:          %20.10lg (ps)"%( step_size ) )
		self.log_function( "Temperature:        %20.10lg (K)"%( temperature ) )
		self.log_function( "Temp. coupling:     %20.10lg (ps)"%( temperature_coupling ) )
		self.log_function( "Pressure:           %20.10lg (atm)"%( pressure ) )
		self.log_function( "Press. Coupling:    %20.10lg"%( pressure_coupling ) )
		if( step_number > 0 ):
			self.log_function( "Step Number:        %20d"%( step_number ) )
		self.log_function( "Print Frequency:    %20d"%( print_frequency ) )
		self.log_function( "\n%20s%20s%20s%20s%20s%20s%20s"%( "Time (ps)", "Potential (kJ/mol)", "Kinetic (kJ/mol)", "Total (kJ/mol)", "Temperature (K)", "Pressure (atm)", "Volume (A^3)" ) )
		self.log_function( 140 * "-" )
		self.tc   = step_size / temperature_coupling
		self.pc   = step_size / pressure_coupling
		self.pcnv = 1.0e33 / ( qm3.constants.NA * 1.01325e5 ) # kJ/mol.A^3 >> atm
		self.xavr = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
		self.xrms = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
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
		func = self.obj.func
		self.T, self.Kin = qm3.actions.dynamics.current_temperature( self.obj, self.project_RT )
		self.volu = self.obj.boxl[0] * self.obj.boxl[1] * self.obj.boxl[2]
		self.pres = ( 2.0 * self.Kin / self.volu - self.__potential_vs_volume() ) / self.obj.size * self.pcnv 
		self.log_function( "%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf"%( 0.0, func, self.Kin, self.obj.func + self.Kin, self.T, self.pres, self.volu ) )
		self.xavr[0] += func
		self.xavr[1] += self.Kin
		self.xavr[2] += func + self.Kin
		self.xavr[3] += self.T
		self.xavr[4] += self.pres
		self.xavr[5] += self.volu
		self.xrms[0] += func * func
		self.xrms[1] += self.Kin * self.Kin
		self.xrms[2] += ( func + self.Kin ) * ( func + self.Kin )
		self.xrms[3] += self.T * self.T
		self.xrms[4] += self.pres * self.pres
		self.xrms[5] += self.volu * self.volu
		if( step_number > 0 ):
			for i in range( step_number ):
				self.integrate()
			self.stats()


	# [ U(Vol+d) - U(Vol-d) ] / 2d 
	def __potential_vs_volume( self, displacement = 3.5183e-8 ):
		bcrd = self.obj.coor[:]
		# forward
		t = 1.0 + displacement
		for i in range( self.obj.size ):
			self.obj.coor[i] = bcrd[i] * t
		self.obj.get_func()
		ff = self.obj.func
		# backward
		t = 1.0 - displacement
		for i in range( self.obj.size ):
			self.obj.coor[i] = bcrd[i] * t
		self.obj.get_func()
		fb = self.obj.func
		# volume
		dvol = ( math.exp( math.log( 1 + displacement ) * 3.0 ) - 1.0 ) * self.obj.boxl[0] * self.obj.boxl[1] * self.obj.boxl[2]
		for i in range( self.obj.size ):
			self.obj.coor[i] = bcrd[i]
		return( 0.5 * ( ff - fb ) / dvol )


	def integrate( self ):
		self.time += self.step_size
		self.istp += 1
		# update velocities
		for i in range( self.obj.size ):
			self.obj.velo[i] += self.step_size * self.vacc[i]
		# scale temperature
		self.T, self.Kin = qm3.actions.dynamics.current_temperature( self.obj, self.project_RT )
		scv = math.sqrt( 1.0 + self.tc * ( self.temperature / self.T - 1.0 ) )
		scv = min( max( scv, 0.9 ), 1.1 )
		for i in range( self.obj.size ):
			self.obj.velo[i] *= scv
		self.T, self.Kin = qm3.actions.dynamics.current_temperature( self.obj, self.project_RT )
		# update coordinates
		for i in range( self.obj.size ):
			self.obj.coor[i] += self.step_size * self.obj.velo[i]
		# scale pressure
		self.pres = ( 2.0 * self.Kin / self.volu - self.__potential_vs_volume() ) / self.obj.size * self.pcnv 
		scp = math.exp( math.log( 1.0 - self.pc * ( self.pressure - self.pres ) / math.fabs( self.pressure - self.pres ) ) / 3.0 )
		for i in [0, 1, 2]:
			self.obj.boxl[i] *= scp
		self.volu = self.obj.boxl[0] * self.obj.boxl[1] * self.obj.boxl[2]
		for i in range( self.obj.size ):
			self.obj.coor[i] *= scp
		# accelerations
		self.vacc = []
		self.obj.get_grad()
		for i in range( self.obj.size // 3 ):
			i3 = i * 3
			for j in [0, 1, 2]:
				self.vacc.append( - self.obj.grad[i3+j] / self.obj.mass[i] * 100.0 )
		if( self.project_RT ):
			qm3.utils.project_RT_modes( self.obj.mass, self.obj.coor, self.vacc, None )
		# stats
		self.xavr[0] += self.obj.func
		self.xavr[1] += self.Kin
		self.xavr[2] += self.obj.func + self.Kin
		self.xavr[3] += self.T
		self.xavr[4] += self.pres
		self.xavr[5] += self.volu
		self.xrms[0] += self.obj.func * self.obj.func
		self.xrms[1] += self.Kin * self.Kin
		self.xrms[2] += ( self.obj.func + self.Kin ) * ( self.obj.func + self.Kin )
		self.xrms[3] += self.T * self.T
		self.xrms[4] += self.pres * self.pres
		self.xrms[5] += self.volu * self.volu
		if( self.istp % self.print_frequency == 0 ):
			self.log_function( "%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf%20.5lf"%( self.time, self.obj.func, self.Kin, self.obj.func + self.Kin, self.T, self.pres, self.volu ) )
		self.obj.current_step( self.istp )


	def stats( self ):
		self.log_function( 140 * "-" )
		savr = "%-20s"%( "Averages:" )
		srms = "%-20s"%( "RMS Deviations:" )
		for i in range( len( self.xavr ) ):
			tavr  = self.xavr[i] / float( self.istp + 1 )
			savr += "%20.5lf"%( tavr )
			trms  = self.xrms[i] / float( self.istp + 1 )
			srms += "%20.5lf"%( math.sqrt( math.fabs( trms - tavr * tavr ) ) )
		self.log_function( savr )
		self.log_function( srms )
		self.log_function( 140 * "-" + "\n" )





if( __name__ == "__main__" ):
	import qm3.actions.minimize
	import qm3.elements
	import qm3.problem
	import qm3.maths.rand
	import qm3.io.dcd
	
	
	class argon_box( qm3.problem.template ):
		def __init__( self, number_of_particles = 100, liquid = True ):
			self.epsi = 1.0451 # kJ/mol
			self.sigm = 3.3450 # A
			# densities in g/cm^3
			if( liquid ):
				density = 1.40
				who = "liq"
			else:
				density = 1.784e-3
				who = "gas"
			# -----------------------------------------
			self.natm = number_of_particles
			lattice   = round( math.exp( math.log( self.natm * qm3.elements.mass[18] * 1.0e24 / ( qm3.constants.NA * density ) ) / 3.0 ), 2 )
			print( ">> Lattice (%d/%s): "%( self.natm, who ) , lattice, "_A" )
			self.boxl = [ lattice, lattice, lattice ]
			self.mass = [ qm3.elements.mass[18] ] * self.natm
			self.size = 3 * self.natm
			self.coor = []
			for i in range( self.size ):
				self.coor.append( lattice * qm3.maths.rand.random() )
			self.func = 0.0
			self.grad = []
	
			self.dcd = qm3.io.dcd.dcd()
			self.dcd.write( "dcd", self.natm )
			self.flg = False
	
			self.s6 = self.sigm * self.sigm * self.sigm
			self.s6 *= self.s6
			self.ef = 4.0 * self.epsi * self.s6
			self.gf = 6.0 * self.ef
	
	
		def get_func( self ):
			self.func = 0.0
			for i in range( self.natm - 1 ):
				i3 = i * 3
				for j in range( i+1, self.natm ):
					j3 = j * 3
					dr = [ ( self.coor[i3+k] - self.coor[j3+k] ) - self.boxl[k] * round( ( self.coor[i3+k] - self.coor[j3+k] ) / self.boxl[k], 0 ) for k in [0, 1, 2] ]
					dd = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]
					ss = 1.0 / ( dd * dd * dd )
					self.func += ss * ( self.s6 * ss - 1.0 )
			self.func *= self.ef
	
		
		def get_grad( self ):
			virial    = 0.0
			self.func = 0.0
			self.grad = [ 0.0 for i in range( self.size ) ]
			for i in range( self.natm - 1 ):
				i3 = i * 3
				for j in range( i+1, self.natm ):
					j3 = j * 3
					dr = [ ( self.coor[i3+k] - self.coor[j3+k] ) - self.boxl[k] * round( ( self.coor[i3+k] - self.coor[j3+k] ) / self.boxl[k], 0 ) for k in [0, 1, 2] ]
					dd = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]
					ss = 1.0 / ( dd * dd * dd )
					self.func += ss * ( self.s6 * ss - 1.0 )
					tt = self.gf * ss / dd * ( 2.0 * self.s6 * ss - 1.0 )
					for k in [0, 1, 2]:
						self.grad[i3+k] -= tt * dr[k]
						self.grad[j3+k] += tt * dr[k]
						virial += tt * dr[k] * dr[k]
			self.func *= self.ef
			return( virial )
	
	
		def xyz_write( self, fname ):
			f = open( fname, "wt" )
			f.write( "%d\n\n"%( self.natm ) )
			for i in range( self.natm ):
				i3 = i * 3
				f.write( "Ar%12.6lf%12.6lf%12.6lf\n"%( self.coor[i3], self.coor[i3+1], self.coor[i3+2] ) )
			f.close()
	
	
		def xyz_read( self, fname ):
			f = open( fname, "rt" )
			f.readline(); f.readline()
			for i in range( self.natm ):
				i3 = i * 3
				self.coor[i3:i3+3] = [ float( j ) for j in f.readline().strip().split()[1:] ]
			f.close()
	
	
		def current_step( self, istep ):
			if( self.flg and istep%1 == 0 ):
				for i in range( self.natm ):
					i3 = i * 3
					self.dcd.X[i] = self.coor[i3]
					self.dcd.Y[i] = self.coor[i3+1]
					self.dcd.Z[i] = self.coor[i3+2]
				self.dcd.append()
	
	
	
	
	obj = argon_box( 200 )
	
#	qm3.actions.minimize.steepest_descent( obj, step_number = 100, step_size = 1., gradient_tolerance = 10. )
#	qm3.actions.minimize.fire( obj, step_number = 1000, step_size = 0.1, gradient_tolerance = 1. )
#	obj.xyz_write( "xyz" )
	
	obj.xyz_read( "xyz" )
	qm3.actions.dynamics.assign_velocities( obj, 80.0, project_RT = True )
	obj.flg = True
	
	leapfrog_verlet( obj, step_size = 0.001, temperature = 80.0, pressure_coupling = 1.0,
		print_frequency = 1, project_RT = True, step_number = 10000 )
	obj.dcd.close()
	print( obj.boxl, "_A" )

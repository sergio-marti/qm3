# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import  sys
if( sys.version_info[0] == 2 ):
        range = xrange
import  math
import	random

random.seed()


def default_log( txt ):
	sys.stdout.write( txt + "\n" )
	sys.stdout.flush()


def differential_evolution( obj,
							boundaries,
							step_number = 1000,
							print_frequency = 10,
							step_tolerance = 1.0e-4,
							population_size = 20,
							mutation_factor = 0.8,
							crossover_probability = 0.7,
							maximum_generations = 100,
							log_function = default_log ):
	# -------------------------------------------------------------------------
	population_size = max( population_size, 10 )
	mutation_factor = min( max( mutation_factor, 0.5 ), 2.0 )
	crossover_probability = min( max( crossover_probability, 0.0 ), 1.0 )
	# -------------------------------------------------------------------------
	log_function( "---------------------------------------- Genetic Minimization (DE: rand/1+bin)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Step Tolerance:     %20.10lg\n"%( step_tolerance ) )
	log_function( "Population Size:    %20d"%( population_size ) )
	log_function( "Mutation Factor:    %20.10lg"%( mutation_factor ) )
	log_function( "Crossover Prob.:    %20.10lg"%( crossover_probability ) )
	log_function( "Max. Generations:   %20d\n"%( maximum_generations ) )
	log_function( "%10s%30s"%( "Step", "Function" ) )
	log_function( "-" * 40 )
	minc = []
	disp = []
	for i in range( obj.size ):
		minc.append( min( boundaries[i] ) )
		disp.append( math.fabs( boundaries[i][1] - boundaries[i][0] ) )
	coor = []
	func = []
	for i in range( population_size ):
		coor.append( [] )
		for j in range( obj.size ):
			coor[-1].append( random.random() )
		obj.coor = [ minc[j] + disp[j] * coor[-1][j] for j in range( obj.size ) ]
		obj.get_func()
		func.append( obj.func )
	ok_fun = min( func )
	ok_crd = coor[func.index( ok_fun )][:]
	log_function( "%10s%30.10lf"%( "", ok_fun ) )
	it = 0
	df = step_tolerance * 2.0
	while( it < step_number and df > step_tolerance ):
		for i in range( population_size ):
			# -------------------------------------------------------------------------
			# rand/1 + binomial
			a, b, c = random.sample( [ j for j in range( population_size ) if j != i ], 3 )
			trial = []
			for j in range( obj.size ):
				if( random.random() < crossover_probability ):
					trial.append( min( max( coor[a][j] + mutation_factor * ( coor[b][j] - coor[c][j] ), 0.0 ), 1.0 ) )
				else:
					trial.append( coor[i][j] )
			# -------------------------------------------------------------------------
			obj.coor = [ minc[j] + disp[j] * trial[j] for j in range( obj.size ) ]
			obj.get_func()
			if( obj.func < func[i] ):
				func[i] = obj.func
				coor[i] = trial[:]
				if( obj.func < ok_fun ):
					df = math.sqrt( sum( [ math.pow( ok_crd[j] - trial[j], 2.0 ) for j in range( obj.size ) ] ) / float( obj.size ) )
					ok_fun = obj.func
					ok_crd = trial[:]
		it += 1
		if( it % print_frequency == 0 ):
			log_function( "%10d%30.10lf"%( it, ok_fun ) )
		obj.current_step( it )
	if( it % print_frequency != 0 ):
		log_function( "%10d%30.10lf"%( it, ok_fun ) )
	log_function( "-" * 40 )
	obj.coor = [ minc[j] + disp[j] * ok_crd[j] for j in range( obj.size ) ]




try:
	import	qm3.utils._mpi
	def mpi_differential_evolution( obj,
							boundaries,
							node, ncpu,
							step_number = 1000,
							print_frequency = 10,
							step_tolerance = 1.0e-4,
							population_size = 10,
							mutation_factor = 0.8,
							crossover_probability = 0.7,
							log_function = default_log ):
		# -------------------------------------------------------------------------
		mutation_factor = min( max( mutation_factor, 0.5 ), 2.0 )
		crossover_probability = min( max( crossover_probability, 0.0 ), 1.0 )
		# -------------------------------------------------------------------------
		if( node == 0 ):
			log_function( "---------------------------------------- Genetic Minimization (DE: rand-to-best/1+bin)\n" )
			log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
			log_function( "Step Number:        %20d"%( step_number ) )
			log_function( "Print Frequency:    %20d"%( print_frequency ) )
			log_function( "Step Tolerance:     %20.10lg\n"%( step_tolerance ) )
			log_function( "Population:         %20d"%( population_size * ncpu ) )
			log_function( "Mutation Factor:    %20.10lg"%( mutation_factor ) )
			log_function( "Crossover Prob.:    %20.10lg\n"%( crossover_probability ) )
			log_function( "%10s%30s"%( "Step", "Function" ) )
			log_function( "-" * 40 )
		minc = []
		disp = []
		for i in range( obj.size ):
			minc.append( min( boundaries[i] ) )
			disp.append( math.fabs( boundaries[i][1] - boundaries[i][0] ) )
		coor = []
		func = []
		for i in range( population_size ):
			coor.append( [] )
			for j in range( obj.size ):
				coor[-1].append( random.random() )
			obj.coor = [ minc[j] + disp[j] * coor[-1][j] for j in range( obj.size ) ]
			obj.get_func()
			func.append( obj.func )
		ok_fun = min( func )
		ok_crd = coor[func.index( ok_fun )][:]
		# sync ----------------------------------------------------------------
		qm3.utils._mpi.barrier()
		if( node == 0 ):
			FUNC = [ ok_fun ]
			COOR = [ ok_crd[:] ]
			for i in range( 1, ncpu ):
				FUNC += qm3.utils._mpi.recv_r8( i, 1 )
				COOR.append( qm3.utils._mpi.recv_r8( i, obj.size ) )
			ok_fun = min( FUNC )
			ok_crd = COOR[FUNC.index( ok_fun )][:]
			for i in range( 1, ncpu ):
				qm3.utils._mpi.send_r8( i, [ ok_fun ] )
				qm3.utils._mpi.send_r8( i, ok_crd )
		else:
			qm3.utils._mpi.send_r8( 0, [ ok_fun ] )
			qm3.utils._mpi.send_r8( 0, ok_crd )
			ok_fun = qm3.utils._mpi.recv_r8( 0, 1 )[0]
			ok_crd = qm3.utils._mpi.recv_r8( 0, obj.size )
		# ---------------------------------------------------------------- sync
		if( node == 0 ):
			log_function( "%10s%30.10lf"%( "", ok_fun ) )
		it = 0
		df = step_tolerance * 2.0
		while( it < step_number and df > step_tolerance ):
			for i in range( population_size ):
				# -------------------------------------------------------------------------
				# rand-to-best/1 + binomial
				a, b, c = random.sample( [ j for j in range( population_size ) if j != i ], 3 )
				trial = []
				for j in range( obj.size ):
					if( random.random() < crossover_probability ):
						trial.append( min( max( coor[a][j] + mutation_factor * ( coor[b][j] - coor[c][j] ) + mutation_factor * ( ok_crd[j] - coor[a][j] ) , 0.0 ), 1.0 ) )
					else:
						trial.append( coor[i][j] )
				# -------------------------------------------------------------------------
				obj.coor = [ minc[j] + disp[j] * trial[j] for j in range( obj.size ) ]
				obj.get_func()
				if( obj.func < func[i] ):
					func[i] = obj.func
					coor[i] = trial[:]
					if( obj.func < ok_fun ):
						ok_fun = obj.func
						ok_crd = trial[:]
			# sync ----------------------------------------------------------------
			qm3.utils._mpi.barrier()
			if( node == 0 ):
				FUNC = [ ok_fun ]
				COOR = [ ok_crd[:] ]
				for i in range( 1, ncpu ):
					FUNC += qm3.utils._mpi.recv_r8( i, 1 )
					COOR.append( qm3.utils._mpi.recv_r8( i, obj.size ) )
				ok_fun = min( FUNC )
				ok_crd = COOR[FUNC.index( ok_fun )][:]
				for i in range( 1, ncpu ):
					qm3.utils._mpi.send_r8( i, [ ok_fun ] )
					qm3.utils._mpi.send_r8( i, ok_crd )
				if( ok_fun < FUNC[0] ):
					df = math.sqrt( sum( [ math.pow( ok_crd[j] - COOR[0][j], 2.0 ) for j in range( obj.size ) ] ) / float( obj.size ) )
				for i in range( 1, ncpu ):
					qm3.utils._mpi.send_r8( i, [ df ] )
			else:
				qm3.utils._mpi.send_r8( 0, [ ok_fun ] )
				qm3.utils._mpi.send_r8( 0, ok_crd )
				ok_fun = qm3.utils._mpi.recv_r8( 0, 1 )[0]
				ok_crd = qm3.utils._mpi.recv_r8( 0, obj.size )
				df = qm3.utils._mpi.recv_r8( 0, 1 )[0]
			# ---------------------------------------------------------------- sync
			it += 1
			if( node == 0 ):
				if( it % print_frequency == 0 ):
					log_function( "%10d%30.10lf"%( it, ok_fun ) )
				obj.current_step( it )
		if( it % print_frequency != 0 and node == 0 ):
			log_function( "%10d%30.10lf"%( it, ok_fun ) )
		if( node == 0 ):
			log_function( "-" * 40 )
			obj.coor = [ minc[j] + disp[j] * ok_crd[j] for j in range( obj.size ) ]
except:
	pass

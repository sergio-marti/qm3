# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import  sys
if( sys.version_info[0] == 2 ):
        range = xrange
import  math
import	qm3.maths.rand
import	qm3.utils.queue
import	multiprocessing
try:
	import	cPickle as pickle
except:
	import	pickle


def default_log( txt ):
	sys.stdout.write( txt + "\n" )
	sys.stdout.flush()


"""
the larger the population, the better the performance (but the larger the calculation...)
	a value of: population = 5 * size   would be ideal whether it is affordable...

diffevo [differential_evolution] is the standard (serialized) implementation, and works
	(pretty) well for large population_sizes. During the 'intercourses', the improved
	new generations are incorporated to the 'active' population... (allowing mixing).

mpi_diffevo is a MPI version, in which the global population is splitted into mpi_ncpu
	groups, and at each generation, the best of each 'village' are compared to choose
	the global best one... ('intercourses' are 'local' with no interaction among 'villages').

smp_diffevo is a SMP implementation, which allows parallel evaluation of the new genetarions...
	Two approaches:
		1) build all the trials, and evaluate them using a multiprocessing.Pool( ncpu )
		2) build slots of ncpu trials, and evaluate them using multiprocessing.Process/qm3.utils.Queque

	While (1) is faster, it does not allow to the improved trials to 'interact' with the
	remaining parents... In (2) best trials and parents are swapped each size//ncpu chunks.

"""


def differential_evolution( obj,
						boundaries,
						step_number = 200,
						print_frequency = 10,
						step_tolerance = 1.0e-6,
						population_size = 10,
						mutation_factor = 0.8,
						crossover_probability = 0.75,
						checkpointing = True,
						log_function = default_log ):
	# -------------------------------------------------------------------------
	# the larger the population, the better the performance (but the larger the calculation...)
	population_size = max( population_size, obj.size * 2 )
	mutation_factor = min( max( mutation_factor, 0.1 ), 1.0 )
	crossover_probability = min( max( crossover_probability, 0.1 ), 1.0 )
	# -------------------------------------------------------------------------
	log_function( "---------------------------------------- Genetic Minimization (DE: rand/1+bin)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Step Tolerance:     %20.10lg\n"%( step_tolerance ) )
	log_function( "Population:         %20d"%( population_size ) )
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
			coor[-1].append( qm3.maths.rand.random() )
		obj.coor = [ minc[j] + disp[j] * coor[-1][j] for j in range( obj.size ) ]
		obj.get_func()
		func.append( obj.func )
	ok_fun = min( func )
	ok_crd = coor[func.index( ok_fun )][:]
	ok_stp = 2.0 * step_tolerance
	log_function( "%10s%30.10lf"%( "", ok_fun ) )
	it = 0
	ff = ok_fun
	while( it < step_number and ok_stp > step_tolerance ):
		for i in range( population_size ):
			# -------------------------------------------------------------------------
			# rand/1 + binomial
			a, b, c = qm3.maths.rand.sample( [ j for j in range( population_size ) if j != i ], 3 )
			trial = []
			for j in range( obj.size ):
				if( qm3.maths.rand.random() < crossover_probability ):
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
					ok_stp = math.sqrt( sum( [ math.pow( ok_crd[j] - trial[j], 2.0 ) for j in range( obj.size ) ] ) / float( obj.size ) )
					ok_fun = obj.func
					ok_crd = trial[:]
		if( ff > ok_fun ):
			ff = ok_fun
			it += 1
			if( it % print_frequency == 0 ):
				log_function( "%10d%30.10lf"%( it, ok_fun ) + " (%.1le)"%( ok_stp ) )
				if( checkpointing ):
					fd = open( "diffevo.chk", "wb" )
					pickle.dump( [ minc[j] + disp[j] * ok_crd[j] for j in range( obj.size ) ], fd )
					fd.close()
			obj.current_step( it )
	if( it % print_frequency != 0 ):
		log_function( "%10d%30.10lf"%( it, ok_fun ) + " (%.1le)"%( ok_stp ) )
	log_function( "-" * 40 )
	obj.coor = [ minc[j] + disp[j] * ok_crd[j] for j in range( obj.size ) ]
	obj.func = ok_fun


diffevo = differential_evolution



try:
	import	qm3.utils._mpi
	def mpi_diffevo( obj,
				boundaries,
				mpi_node, mpi_ncpu,
				step_number = 200,
				print_frequency = 10,
				step_tolerance = 1.0e-6,
				population_size = 10,
				mutation_factor = 0.8,
				crossover_probability = 0.75,
				checkpointing = True,
				log_function = default_log ):
		# -------------------------------------------------------------------------
		mutation_factor = min( max( mutation_factor, 0.1 ), 1.0 )
		crossover_probability = min( max( crossover_probability, 0.1 ), 1.0 )
		# -------------------------------------------------------------------------
		if( mpi_node == 0 ):
			log_function( "---------------------------------------- Genetic Minimization (MPI/DE: rand-to-best/1+bin)\n" )
			log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
			log_function( "Step Number:        %20d"%( step_number ) )
			log_function( "Print Frequency:    %20d"%( print_frequency ) )
			log_function( "Step Tolerance:     %20.10lg\n"%( step_tolerance ) )
			log_function( "MPI Population:     %20s"%( "%d x %d"%( population_size, mpi_ncpu ) ) )
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
				coor[-1].append( qm3.maths.rand.random() )
			obj.coor = [ minc[j] + disp[j] * coor[-1][j] for j in range( obj.size ) ]
			obj.get_func()
			func.append( obj.func )
		ok_fun = min( func )
		ok_crd = coor[func.index( ok_fun )][:]
		ok_stp = 2.0 * step_tolerance
		# sync ----------------------------------------------------------------
		if( mpi_ncpu > 1 ):
			qm3.utils._mpi.barrier()
			if( mpi_node == 0 ):
				FUNC = [ ok_fun ]
				COOR = [ ok_crd[:] ]
				for i in range( 1, mpi_ncpu ):
					FUNC += qm3.utils._mpi.recv_r8( i, 1 )
					COOR.append( qm3.utils._mpi.recv_r8( i, obj.size ) )
				ok_fun = min( FUNC )
				ok_crd = COOR[FUNC.index( ok_fun )][:]
				for i in range( 1, mpi_ncpu ):
					qm3.utils._mpi.send_r8( i, [ ok_fun ] )
					qm3.utils._mpi.send_r8( i, ok_crd )
			else:
				qm3.utils._mpi.send_r8( 0, [ ok_fun ] )
				qm3.utils._mpi.send_r8( 0, ok_crd )
				ok_fun = qm3.utils._mpi.recv_r8( 0, 1 )[0]
				ok_crd = qm3.utils._mpi.recv_r8( 0, obj.size )
		# ---------------------------------------------------------------- sync
		if( mpi_node == 0 ):
			log_function( "%10s%30.10lf"%( "", ok_fun ) )
		it = 0
		qq = 1
		ff = ok_fun
		while( qq == 1 ):
			for i in range( population_size ):
				# -------------------------------------------------------------------------
#				# rand/1 + binomial
#				a, b, c = qm3.maths.rand.sample( [ j for j in range( population_size ) if j != i ], 3 )
#				trial = []
#				for j in range( obj.size ):
#					if( qm3.maths.rand.random() < crossover_probability ):
#						trial.append( min( max( coor[a][j] + mutation_factor * ( coor[b][j] - coor[c][j] ), 0.0 ), 1.0 ) )
#					else:
#						trial.append( coor[i][j] )
				# -------------------------------------------------------------------------
				# rand-to-best/1 + binomial
				a, b = qm3.maths.rand.sample( [ j for j in range( population_size ) if j != i ], 2 )
				trial = []
				for j in range( obj.size ):
					if( qm3.maths.rand.random() < crossover_probability ):
						trial.append( min( max( coor[i][j] + mutation_factor * ( ok_crd[j] - coor[i][j] + coor[a][j] - coor[b][j] ), 0.0 ), 1.0 ) )
					else:
						trial.append( coor[i][j] )
				# -------------------------------------------------------------------------
				obj.coor = [ minc[j] + disp[j] * trial[j] for j in range( obj.size ) ]
				obj.get_func()
				if( obj.func < func[i] ):
					func[i] = obj.func
					coor[i] = trial[:]
					if( obj.func < ok_fun ):
						ok_stp = math.sqrt( sum( [ math.pow( ok_crd[j] - trial[j], 2.0 ) for j in range( obj.size ) ] ) / float( obj.size ) )
						ok_fun = obj.func
						ok_crd = trial[:]
			# sync ----------------------------------------------------------------
			if( mpi_ncpu > 1 ):
				qm3.utils._mpi.barrier()
				if( mpi_node == 0 ):
					STEP = [ ok_stp ]
					FUNC = [ ok_fun ]
					COOR = [ ok_crd[:] ]
					for i in range( 1, mpi_ncpu ):
						STEP += qm3.utils._mpi.recv_r8( i, 1 )
						FUNC += qm3.utils._mpi.recv_r8( i, 1 )
						COOR.append( qm3.utils._mpi.recv_r8( i, obj.size ) )
					ok_fun = min( FUNC )
					ok_crd = COOR[FUNC.index( ok_fun )][:]
					ok_stp = STEP[FUNC.index( ok_fun )]
					for i in range( 1, mpi_ncpu ):
						qm3.utils._mpi.send_r8( i, [ ok_stp ] )
						qm3.utils._mpi.send_r8( i, [ ok_fun ] )
						qm3.utils._mpi.send_r8( i, ok_crd )
				else:
					qm3.utils._mpi.send_r8( 0, [ ok_stp ] )
					qm3.utils._mpi.send_r8( 0, [ ok_fun ] )
					qm3.utils._mpi.send_r8( 0, ok_crd )
					ok_stp = qm3.utils._mpi.recv_r8( 0, 1 )[0]
					ok_fun = qm3.utils._mpi.recv_r8( 0, 1 )[0]
					ok_crd = qm3.utils._mpi.recv_r8( 0, obj.size )
			# ---------------------------------------------------------------- sync
			if( mpi_node == 0 ):
				if( ff > ok_fun ):
					ff = ok_fun
					it += 1
					qq = int( it < step_number and ok_stp > step_tolerance )
					if( it % print_frequency == 0 ):
						log_function( "%10d%30.10lf"%( it, ok_fun ) + " (%.1le)"%( ok_stp ) )
					obj.current_step( it )
			if( mpi_ncpu > 1 ):
				qm3.utils._mpi.barrier()
				if( mpi_node == 0 ):
					for i in range( 1, mpi_ncpu ):
						qm3.utils._mpi.send_i4( i, [ qq ] )
				else:
					qq = qm3.utils._mpi.recv_i4( 0, 1 )[0]
		if( mpi_node == 0 ):
			if( it % print_frequency != 0 ):
				log_function( "%10d%30.10lf"%( it, ok_fun ) + " (%.1le)"%( ok_stp ) )
				if( checkpointing ):
					fd = open( "diffevo.chk", "wb" )
					pickle.dump( [ minc[j] + disp[j] * ok_crd[j] for j in range( obj.size ) ], fd )
					fd.close()
			log_function( "-" * 40 )
			obj.coor = [ minc[j] + disp[j] * ok_crd[j] for j in range( obj.size ) ]
			obj.func = ok_fun
		else:
			log_function( "mpi_node = %d: done!"%( mpi_node ) )
except:
	pass



def _smp_diffevo_fitness_pool( args ):
	who, obj, vec = args
	obj.coor = vec[:]
	obj.get_func()
	return( ( who, obj.func ) )


def __smp_diffevo_fitness_process( que, who, obj, vec ):
	obj.coor = vec[:]
	obj.get_func()
	que.client( [ ( who, obj.func ) ] )


def smp_diffevo( objs,
			boundaries,
			step_number = 200,
			print_frequency = 10,
			step_tolerance = 1.0e-6,
			population_size = 10,
			mutation_factor = 0.8,
			crossover_probability = 0.75,
			checkpointing = True,
			log_function = default_log ):
	# -------------------------------------------------------------------------
	ncpu = len( objs )
	size = objs[0].size
	population_size = max( population_size, size * 2 )
	population_size = ( population_size // ncpu + ( population_size % ncpu != 0 ) ) * ncpu
	chnk = population_size // ncpu
	mutation_factor = min( max( mutation_factor, 0.1 ), 1.0 )
	crossover_probability = min( max( crossover_probability, 0.1 ), 1.0 )
	# -------------------------------------------------------------------------
	log_function( "---------------------------------------- Genetic Minimization (SMP/DE: rand-to-best/1+bin)\n" )
	log_function( "Degrees of Freedom: %20ld"%( size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Step Tolerance:     %20.10lg\n"%( step_tolerance ) )
	log_function( "CPUs Number:        %20d"%( ncpu ) )
	log_function( "Population:         %20d"%( population_size ) )
	log_function( "Pop. Chunks:        %20d"%( chnk ) )
	log_function( "Mutation Factor:    %20.10lg"%( mutation_factor ) )
	log_function( "Crossover Prob.:    %20.10lg\n"%( crossover_probability ) )
	log_function( "%10s%30s"%( "Step", "Function" ) )
	log_function( "-" * 40 )

	minc = []
	disp = []
	for i in range( size ):
		minc.append( min( boundaries[i] ) )
		disp.append( math.fabs( boundaries[i][1] - boundaries[i][0] ) )

	coor = []
	func = []
	args = []
	for i in range( chnk ):
		for j in range( ncpu ):
			coor.append( [] )
			for k in range( size ):
				coor[-1].append( qm3.maths.rand.random() )
			args.append( [ i * ncpu + j, objs[j], [ minc[k] + disp[k] * coor[-1][k] for k in range( size ) ] ] )
	work = multiprocessing.Pool( ncpu )
	func = [ j for i,j in sorted( work.map( _smp_diffevo_fitness_pool, args ) ) ]
	work.close()
	work.join()

	ok_fun = min( func )
	ok_crd = coor[func.index( ok_fun )][:]
	ok_stp = 2.0 * step_tolerance
	log_function( "%10s%30.10lf"%( "", ok_fun ) )
	it = 0
	ff = ok_fun
	while( it < step_number and ok_stp > step_tolerance ):
# =============================================================================================
# Model: 1
		args = []
		T    = []
		for i in range( chnk ):
			for j in range( ncpu ):
				T.append( [] )
				w = i * ncpu + j
				# -------------------------------------------------------------------------
#				# rand/1 + binomial
#				a, b, c = qm3.maths.rand.sample( [ k for k in range( population_size ) if k != w ], 3 )
#				for k in range( size ):
#					if( qm3.maths.rand.random() < crossover_probability ):
#						T[-1].append( min( max( coor[a][k] + mutation_factor * ( coor[b][k] - coor[c][k] ), 0.0 ), 1.0 ) )
#					else:
#						T[-1].append( coor[w][k] )
				# -------------------------------------------------------------------------
				# rand-to-best/1 + binomial
				a, b = qm3.maths.rand.sample( [ k for k in range( population_size ) if k != i ], 2 )
				for k in range( size ):
					if( qm3.maths.rand.random() < crossover_probability ):
						T[-1].append( min( max( coor[i][k] + mutation_factor * ( ok_crd[k] - coor[i][k] + coor[a][k] - coor[b][k] ), 0.0 ), 1.0 ) )
					else:
						T[-1].append( coor[i][k] )
				args.append( [ w, objs[j], [ minc[k] + disp[k] * T[-1][k] for k in range( size ) ] ] )
		work = multiprocessing.Pool( ncpu )
		F = [ j for i,j in sorted( work.map( _smp_diffevo_fitness_pool, args ) ) ]
		work.close()
		work.join()
		for i in range( population_size ):
			if( F[i] < func[i] ):
				func[i] = F[i]
				coor[i] = T[i][:]
				if( F[i] < ok_fun ):
					ok_stp = math.sqrt( sum( [ math.pow( ok_crd[k] - T[i][k], 2.0 ) for k in range( size ) ] ) / float( size ) )
					ok_fun = F[i]
					ok_crd = T[i][:]
# =============================================================================================
# Model: 2
#		for i in range( chnk ):
#			Q = qm3.utils.queue.Queue( ncpu )
#			T = []
#			for j in range( ncpu ):
#				T.append( [] )
#				w = i * ncpu + j
#				# -------------------------------------------------------------------------
#				# rand/1 + binomial
#				a, b, c = qm3.maths.rand.sample( [ k for k in range( population_size ) if k != w ], 3 )
#				for k in range( size ):
#					if( qm3.maths.rand.random() < crossover_probability ):
#						T[-1].append( min( max( coor[a][k] + mutation_factor * ( coor[b][k] - coor[c][k] ), 0.0 ), 1.0 ) )
#					else:
#						T[-1].append( coor[w][k] )
#				multiprocessing.Process( target = __smp_diffevo_fitness_process,
#					args = ( Q, w, objs[j], [ minc[k] + disp[k] * T[-1][k] for k in range( size ) ] ) ).start()
#			Q.serve()
#			F = sorted( Q.data )
#			for j in range( ncpu ):
#				if( F[j][1] < func[F[j][0]] ):
#					func[F[j][0]] = F[j][1]
#					coor[F[j][0]] = T[j][:]
#					if( F[j][1] < ok_fun ):
#						ok_stp = math.sqrt( sum( [ math.pow( ok_crd[k] - T[j][k], 2.0 ) for k in range( size ) ] ) / float( size ) )
#						ok_fun = F[j][1]
#						ok_crd = T[j][:]
		if( ff > ok_fun ):
			ff = ok_fun
			it += 1
			if( it % print_frequency == 0 ):
				log_function( "%10d%30.10lf"%( it, ok_fun ) + " (%.1le)"%( ok_stp ) )
				if( checkpointing ):
					fd = open( "diffevo.chk", "wb" )
					pickle.dump( [ minc[j] + disp[j] * ok_crd[j] for j in range( size ) ], fd )
					fd.close()
	if( it % print_frequency != 0 ):
		log_function( "%10d%30.10lf"%( it, ok_fun ) + " (%.1le)"%( ok_stp ) )
	log_function( "-" * 40 )
	objs[0].coor = [ minc[j] + disp[j] * ok_crd[j] for j in range( size ) ]
	objs[0].func = ok_fun



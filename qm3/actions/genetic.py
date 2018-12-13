# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import  sys
if( sys.version_info[0] == 2 ):
        range = xrange
import  math
import	qm3.maths.rand


def default_log( txt ):
	sys.stdout.write( txt + "\n" )
	sys.stdout.flush()


try:
	import	qm3.utils._mpi
	__has_mpi = True
except:
	__has_mpi = False


def differential_evolution( obj,
						boundaries,
						step_number = 200,
						print_frequency = 10,
						step_tolerance = 1.0e-6,
						population_size = 10,
						mutation_factor = 0.8,
						crossover_probability = 0.75,
						mpi_node = 0, mpi_ncpu = 1,
						log_function = default_log ):
	# -------------------------------------------------------------------------
	population_size = max( population_size, 10 )
	mutation_factor = min( max( mutation_factor, 0.1 ), 1.0 )
	crossover_probability = min( max( crossover_probability, 0.1 ), 1.0 )
	# -------------------------------------------------------------------------
	if( mpi_node == 0 ):
		log_function( "---------------------------------------- Genetic Minimization (DE: rand/1+bin)\n" )
		log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
		log_function( "Step Number:        %20d"%( step_number ) )
		log_function( "Print Frequency:    %20d"%( print_frequency ) )
		log_function( "Step Tolerance:     %20.10lg\n"%( step_tolerance ) )
		if( mpi_ncpu > 1 ):
			log_function( "MPI Population:     %20s"%( "%d x %d"%( population_size, mpi_ncpu ) ) )
		else:
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
	# sync ----------------------------------------------------------------
	if( __has_mpi and mpi_ncpu > 1 ):
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
			# rand-to-best/1 + binomial
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
		# sync ----------------------------------------------------------------
		if( __has_mpi and mpi_ncpu > 1 ):
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
		if( __has_mpi and mpi_ncpu > 1 ):
			qm3.utils._mpi.barrier()
			if( mpi_node == 0 ):
				for i in range( 1, mpi_ncpu ):
					qm3.utils._mpi.send_i4( i, [ qq ] )
			else:
				qq = qm3.utils._mpi.recv_i4( 0, 1 )[0]
	if( mpi_node == 0 ):
		if( it % print_frequency != 0 ):
			log_function( "%10d%30.10lf"%( it, ok_fun ) + " (%.1le)"%( ok_stp ) )
		log_function( "-" * 40 )
		obj.coor = [ minc[j] + disp[j] * ok_crd[j] for j in range( obj.size ) ]
	else:
		log_function( "mpi_node = %d: done!"%( mpi_node ) )



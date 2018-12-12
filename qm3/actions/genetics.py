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
							step_number = 3000,
							print_frequency = 10,
							step_tolerance = 1.0e-4,
							population_size = 20,
							mutation_factor = 0.8,
							crossover_probability = 0.7,
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
	if( it % print_frequency != 0 ):
		log_function( "%10d%30.10lf"%( it, ok_fun ) )
	log_function( "-" * 40 )
	obj.coor = [ minc[j] + disp[j] * ok_crd[j] for j in range( obj.size ) ]

# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.maths.matrix



def default_log( txt ):
	sys.stdout.write( txt + "\n" )
	sys.stdout.flush()



def __project_RT_modes( w, x, g, h ):
	s  = len( x )
	mt = 0.0
	mc = [ 0.0, 0.0, 0.0 ]
	for i in range( s // 3 ):
		i3 = i * 3
		mt += w[i3] * w[i3]
		for j in [0, 1, 2]:
			mc[j] += x[i3+j] * w[i3]
	mc = [ mc[i] / mt for i in [0, 1, 2] ]
	rt = [ 0.0 for i in range( 6 * s ) ]
	for i in range( s // 3 ):
		j = 3 * i
		rt[j]	    = w[j]
		rt[s+j+1]   = w[j]
		rt[2*s+j+2] = w[j]
		rt[3*s+j+1] = - ( x[j+2] - mc[2] * w[j] )
		rt[3*s+j+2] =   ( x[j+1] - mc[1] * w[j] )
		rt[4*s+j]   =   ( x[j+2] - mc[2] * w[j] )
		rt[4*s+j+2] = - ( x[j  ] - mc[0] * w[j] )
		rt[5*s+j]   = - ( x[j+1] - mc[1] * w[j] )
		rt[5*s+j+1] =   ( x[j  ] - mc[0] * w[j] )
	for i in range( 6 ):
		for j in range( i ):
			tmp = sum( [ rt[i*s+k] * rt[j*s+k] for k in range( s ) ] )
			for k in range( s ):
				rt[i*s+k] -= tmp * rt[j*s+k]
		tmp = math.sqrt( sum( [ rt[i*s+k] * rt[i*s+k] for k in range( s ) ] ) )
		for k in range( s ):
			rt[i*s+k] /= tmp
	# gradient
	for i in range( 6 ):
		tmp = sum( [ g[k] * rt[i*s+k] for k in range( s ) ] )
		for k in range( s ):
			g[k] -= tmp * rt[i*s+k]
	# hessian
	ix = [ 0.0 for i in range( s * s ) ]
	for i in range( s ):
		ix[s*i+i] += 1.
		for j in range( s ):
			for k in range( 6 ):
				ix[s*i+j] -= rt[k*s+i] * rt[k*s+j]
	t = qm3.maths.matrix.mult( ix, s, s, qm3.maths.matrix.mult( h, s, s, ix, s, s ), s, s )
	for i in range( s * s ):
		h[i] = t[i]



def steepest_descent( obj, 
			step_number = 100,
			step_size = 0.0028,			# use positive/forward or negative/reverse
			gradient_tolerance = 0.1,
			print_frequency = 10,
			project_RT = True,
			from_saddle = True,
			avoid_recrossing = True,
			log_function = default_log ):
	log_function( "\n---------------------------------------- Minimum Path (SD)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Step Size:          %20.10lg"%( step_size ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Gradient Tolerance: %20.10lg"%( gradient_tolerance ) )
	log_function( "Project RT modes:   %20s"%( project_RT ) )
	log_function( "From Saddle:        %20s"%( from_saddle ) )
	log_function( "Avoid Recrossing:   %20s\n"%( avoid_recrossing ) )
	log_function( "%10s%20s%20s%10s"%( "Step", "Function", "Gradient", "Nskip" ) )
	log_function( "-" * 60 )
	# -- TODO --
	log_function( "-" * 60 + "\n" )



def baker( obj, 
			step_number = 100,
			step_size = 0.0028,			# use positive/forward or negative/reverse
			gradient_tolerance = 0.1,
			print_frequency = 10,
			project_RT = True,
			from_saddle = True,
			avoid_recrossing = True,
			log_function = default_log ):
	log_function( "\n---------------------------------------- Minimum Path (Baker)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Step Size:          %20.10lg"%( step_size ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Gradient Tolerance: %20.10lg"%( gradient_tolerance ) )
	log_function( "Project RT modes:   %20s"%( project_RT ) )
	log_function( "From Saddle:        %20s"%( from_saddle ) )
	log_function( "Avoid Recrossing:   %20s\n"%( avoid_recrossing ) )
	log_function( "%10s%20s%20s%10s"%( "Step", "Function", "Gradient", "Nskip" ) )
	log_function( "-" * 60 )
	vcut = 0.00035481432270250985
	lrge = 1.0e+6
	step = 50.0
	tol2 = 1.0e-8
	mxit = 999
	s    = min( len( obj.mass ), obj.size )
	k    = obj.size // s
	w    = [ 0.0 for i in range( obj.size ) ]
	for i in range( s ):
		w[i*k] = math.sqrt( obj.mass[i] )
		for j in range( 1, k ):
			w[i*k+j] = w[i*k]
	dx = [ 0.0 for i in range( obj.size ) ]
	gx = [ 0.0 for i in range( obj.size ) ]
	x  = [ obj.coor[i] * w[i] for i in range( obj.size ) ]
	if( from_saddle ):
		obj.get_hess()
		h = []
		k = 0
		for i in range( obj.size ):
			for j in range( obj.size ):
				h.append( obj.hess[k] / ( w[i] * w[j] ) )
				k += 1
		g = [ obj.grad[i] / w[i] for i in range( obj.size ) ]
		if( project_RT ):
			__project_RT_modes( w, x, g, h )
		val, vec = qm3.maths.matrix.diag( h, obj.size )
		nskp = sum( [ 1 for i in range( obj.size ) if val[i] < vcut ] )
		tmp  = step_size / math.sqrt( sum( [ vec[i*obj.size] * vec[i*obj.size] for i in range( obj.size ) ] ) )
		dx   = [ vec[i*obj.size] * tmp for i in range( obj.size ) ]
		if( avoid_recrossing ):
			ox   = dx[:]
	else:
		nskp = 7
	step_size = math.fabs( step_size )
	mskp      = 6 * project_RT
	grms      = gradient_tolerance * 2.0
	it1       = 0
	flg       = True
	while( it1 < step_number and ( grms > gradient_tolerance or nskp > mskp ) and flg ):
		for i in range( obj.size ):
			x[i] += dx[i]
			obj.coor[i] = x[i] / w[i]
		obj.current_step( it1 )
		obj.get_hess()
		h = []
		k = 0
		for i in range( obj.size ):
			for j in range( obj.size ):
				h.append( obj.hess[k] / ( w[i] * w[j] ) )
				k += 1
		g = [ obj.grad[i] / w[i] for i in range( obj.size ) ]
		if( project_RT ):
			__project_RT_modes( w, x, g, h )
		val, vec = qm3.maths.matrix.diag( h, obj.size )
		nskp = sum( [ 1 for i in range( obj.size ) if val[i] < vcut ] )
		grms = math.sqrt( sum( [ g[i] * g[i] for i in range( obj.size ) ] ) )
		# transform gradient vector to the local hessian modes
		for i in range( obj.size ):
			dx[i]  = 0.0
			val[i] /= grms
			gx[i]  = sum( [ g[j] * vec[j*obj.size+i] for j in range( obj.size ) ] )
		# minimize along the selected modes (skip first)
		lmbd = 0.0
		if( val[0] < 0.0 ):
			lmbd = val[0] - step
			l1   = val[0]
			l2   = - lrge
		ovr = sum( [ gx[j] * gx[j] / ( lmbd - val[j] ) for j in range( obj.size ) ] )
		i   = 0
		while( i < mxit and math.fabs( lmbd - ovr ) >= tol2 ):
			if( val[0] > 0.0 ):
				lmbd = ovr;
			else:
				if( ovr < lmbd ):
					l1 = lmbd;
				if( ovr > lmbd ):
					l2 = lmbd;
				if( l2 > - lrge ):
					lmbd = 0.5 * ( l1 + l2 )
				elif( l2 == - lrge ):
					lmbd -= step;
			ovr = sum( [ gx[j] * gx[j] / ( lmbd - val[j] ) for j in range( obj.size ) ] )
			i += 1
		if( i > mxit ):
			log_function( "\n -- Too much lambda iterations..." )
			flg = False
		# check final step (too small or large...)
		for i in range( obj.size ):
			dx[i] += sum( [ vec[i*obj.size+j] * gx[j] / ( lmbd - val[j] ) for j in range( obj.size ) ] )
		ovr = math.sqrt( sum( [ dx[i] * dx[i] for i in range( obj.size ) ] ) )
		if( ovr < tol2 ):
			log_function( "\n -- The step size is *very* small..." )
			flg = False
		# scale long steps...
		if( ovr > step_size ):
			for i in range( obj.size ):
				dx[i] *= step_size / ovr
		# avoid recrossing
		if( avoid_recrossing and nskp > mskp ):
			tmp = sum( [ ox[i] * dx[i] for i in range( obj.size ) ] )
			if( tmp < 0.0 ):
				dx = [ -dx[i] for i in range( obj.size ) ]
		grms /= math.sqrt( float( obj.size ) )
		it1 += 1
		if( it1%print_frequency == 0 ):
			log_function( "%10ld%20.5lf%20.10lf%10ld"%( it1, obj.func, grms, nskp ) )
	if( it1%print_frequency != 0 ):
		log_function( "%10ld%20.5lf%20.10lf%10ld"%( it1, obj.func, grms, nskp ) )
	log_function( "-" * 60 + "\n" )



def page_mciver( obj, 
			step_number = 100,
			step_size = 0.0028,			# use positive/forward or negative/reverse
			gradient_tolerance = 0.1,
			print_frequency = 10,
			project_RT = True,
			from_saddle = True,
			avoid_recrossing = True,
			log_function = default_log ):
	log_function( "\n---------------------------------------- Minimum Path (Page-McIver)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Step Size:          %20.10lg"%( step_size ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Gradient Tolerance: %20.10lg"%( gradient_tolerance ) )
	log_function( "Project RT modes:   %20s"%( project_RT ) )
	log_function( "From Saddle:        %20s"%( from_saddle ) )
	log_function( "Avoid Recrossing:   %20s\n"%( avoid_recrossing ) )
	log_function( "%10s%20s%20s%10s"%( "Step", "Function", "Gradient", "Nskip" ) )
	log_function( "-" * 60 )
	vcut = 0.00035481432270250985
	it2m = 1000
	it3m = 100000
	s    = min( len( obj.mass ), obj.size )
	k    = obj.size // s
	w    = [ 0.0 for i in range( obj.size ) ]
	for i in range( s ):
		w[i*k] = math.sqrt( obj.mass[i] )
		for j in range( 1, k ):
			w[i*k+j] = w[i*k]
	dx = [ 0.0 for i in range( obj.size ) ]
	v  = [ 0.0 for i in range( obj.size ) ]
	x  = [ obj.coor[i] * w[i] for i in range( obj.size ) ]
	if( from_saddle ):
		obj.get_hess()
		h = []
		k = 0
		for i in range( obj.size ):
			for j in range( obj.size ):
				h.append( obj.hess[k] / ( w[i] * w[j] ) )
				k += 1
		g = [ obj.grad[i] / w[i] for i in range( obj.size ) ]
		if( project_RT ):
			__project_RT_modes( w, x, g, h )
		val, vec = qm3.maths.matrix.diag( h, obj.size )
		nskp = sum( [ 1 for i in range( obj.size ) if val[i] < vcut ] )
		tmp  = step_size / math.sqrt( sum( [ vec[i*obj.size] * vec[i*obj.size] for i in range( obj.size ) ] ) )
		dx   = [ vec[i*obj.size] * tmp for i in range( obj.size ) ]
		if( avoid_recrossing ):
			ox   = dx[:]
	else:
		nskp = 7
	step_size = math.fabs( step_size )
	mskp      = 6 * project_RT
	grms      = gradient_tolerance * 2.0
	it1       = 0
	flg       = True
	while( it1 < step_number and ( grms > gradient_tolerance or nskp > mskp ) and flg ):
		for i in range( obj.size ):
			x[i] += dx[i]
			obj.coor[i] = x[i] / w[i]
		obj.current_step( it1 )
		obj.get_hess()
		h = []
		k = 0
		for i in range( obj.size ):
			for j in range( obj.size ):
				h.append( obj.hess[k] / ( w[i] * w[j] ) )
				k += 1
		g = [ obj.grad[i] / w[i] for i in range( obj.size ) ]
		if( project_RT ):
			__project_RT_modes( w, x, g, h )
		val, vec = qm3.maths.matrix.diag( h, obj.size )
		nskp = sum( [ 1 for i in range( obj.size ) if val[i] < vcut ] )
		grms = math.sqrt( sum( [ g[i] * g[i] for i in range( obj.size ) ] ) )
		# transform gradient vector to the local hessian modes
		for i in range( obj.size ):
			g[i]   /= grms
			val[i] /= grms
		for i in range( obj.size ):
			v[i]  = sum( [ g[j] * vec[j*obj.size+i] for j in range( obj.size ) ] )
		# search for the PM step
		pm_dt = 0.2 * step_size
		pm_t  = 1.e10
		pm_ot = 0.0
		it2   = 0
		it3   = 0
		while( math.fabs( 1.0 - pm_ot / pm_t ) > 1.0e-6 and it2 < it2m and it3 < it3m ):
			it2  += 1
			pm_ot = pm_t
			pm_dt = 0.5 * pm_dt
			pm_t  = 0.0
			pm_ft = math.sqrt( sum( [ math.pow( v[i] * math.exp( - val[i] * pm_t ), 2.0 ) for i in range( obj.size ) ] ) )
			pm_s  = 0.0;
			it3   = 0;
			while( pm_s < step_size and it3 < it3m ):
				it3  += 1
				pm_os = pm_s
				pm_of = pm_ft
				pm_t += pm_dt
				pm_ft = math.sqrt( sum( [ math.pow( v[i] * math.exp( - val[i] * pm_t ), 2.0 ) for i in range( obj.size ) ] ) )
				pm_s += 0.5 * pm_dt * ( pm_ft + pm_of )
			# does not converge if it cannot reach the minimum with the current step-size (we have ended...)
			if( pm_os != pm_s ):
				pm_t -= ( step_size - pm_s ) * pm_dt / ( pm_os - pm_s )
			else:
				flg = False
		if( math.fabs( 1.0 - pm_t / pm_ot ) <= 1.0e-6 and flg ):
			for i in range( obj.size ):
				tmp = val[i] * pm_t
				if( math.fabs( tmp ) < 1.e-8 ):
					v[i] *= - pm_t * ( 1.0 - tmp  * ( 0.5 - tmp / 6.0 ) )
				else:
					v[i] *= ( math.exp( - tmp ) - 1.0 ) / val[i]
			for i in range( obj.size ):
				dx[i] = sum( [ v[j] * vec[i*obj.size+j] for j in range( obj.size ) ] )
			# avoid recrossing
			if( avoid_recrossing and nskp > mskp ):
				tmp = sum( [ ox[i] * dx[i] for i in range( obj.size ) ] )
				if( tmp < 0.0 ):
					dx = [ -dx[i] for i in range( obj.size ) ]
		grms /= math.sqrt( float( obj.size ) )
		it1 += 1
		if( it1%print_frequency == 0 ):
			log_function( "%10ld%20.5lf%20.10lf%10ld"%( it1, obj.func, grms, nskp ) )
	if( it1%print_frequency != 0 ):
		log_function( "%10ld%20.5lf%20.10lf%10ld"%( it1, obj.func, grms, nskp ) )
	log_function( "-" * 60 + "\n" )




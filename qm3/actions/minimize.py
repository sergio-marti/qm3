# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.maths.matrix


try:
	import qm3.actions._minimize
	has_minimize_so = True
except:
	has_minimize_so = False



def __grms( vec ):
#	o = 0.0
#	for i in vec:
#		o += i * i
#	o = math.sqrt( o )
	o = math.sqrt( sum( [ i*i for i in vec ] ) )
	return( o, o / math.sqrt( len( vec ) ) )



def default_log( txt ):
	sys.stdout.write( txt + "\n" )
	sys.stdout.flush()



def downhill( obj, 
				step_number = 1000,
				step_size = 0.1,
				print_frequency = 10,
				step_tolerance = 1.0e-6,
				log_function = default_log ):
	log_function( "---------------------------------------- Minimization (DH)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Step Size:          %20.10lg"%( step_size ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Step Tolerance:     %20.10lg\n"%( step_tolerance ) )
	log_function( "%10s%20s%20s"%( "Step", "Function", "Displacement" ) )
	log_function( "-" * 50 )
	x = [ obj.coor[:] ]
	for i in range( obj.size ):
		x.append( obj.coor[:] )
		x[i+1][i] += step_size
	f = []
	for i in range( obj.size + 1 ):
		obj.coor = x[i][:]
		obj.get_func()
		f.append( obj.func )
	n  = float( obj.size )
	k  = 0
	df = step_tolerance * 2.0
	while( k < step_number and df > step_tolerance ):
		iLo = f.index( min( f ) )
		iHi = f.index( max( f ) )
		d = []
		for i in range( obj.size ):
			d.append( sum( [ x[j][i] for j in range( obj.size + 1 ) ] ) )
		for i in range( obj.size ):
			d[i] = ( d[i] - (n + 1.0) * x[iHi][i] ) / n
		df = math.sqrt( sum( [ d[i] * d[i] for i in range( obj.size ) ] ) / n )
		if( df > step_tolerance ):
			# Try Reflection
			xNew = [ x[iHi][i] + 2.0 * d[i] for i in range( obj.size ) ]
			obj.coor = xNew[:]
			obj.get_func()
			fNew = obj.func
			if( fNew <= f[iLo] ): # (accept reflection)
				x[iHi] = xNew[:]
				f[iHi] = fNew
				# Try Expanding the Reflection
				xNew = [ x[iHi][i] + d[i] for i in range( obj.size ) ]
				obj.coor = xNew[:]
				obj.get_func()
				fNew = obj.func
				if( fNew <= f[iLo] ): # (accept expansion)
					x[iHi] = xNew[:]
					f[iHi] = fNew
			else:
				# Try Reflection
				if( fNew <= f[iHi] ): # (accept reflection)
					x[iHi] = xNew[:]
					f[iHi] = fNew
				else:
					# Try Contraction
					xNew = [ x[iHi][i] + 0.5 * d[i] for i in range( obj.size ) ]
					obj.coor = xNew[:]
					obj.get_func()
					fNew = obj.func
					if( fNew <= f[iHi] ): # (accept contraction)
						x[iHi] = xNew[:]
						f[iHi] = fNew
					else:
						# Use Shrinkage
						for i in range( obj.size + 1 ):
							if( i != iLo ):
								x[i] = [ x[i][j] - 0.5 * x[iLo][j] for j in range( obj.size ) ]
								obj.coor = x[i][:]
								obj.get_func()
								f[i] = obj.func
		k += 1
		if( k%print_frequency == 0 ):
			log_function( "%10d%20.5lf%20.10lf"%( k, f[iLo], df ) )
		obj.current_step( k )
	if( k%print_frequency != 0 ):
		log_function( "%10d%20.5lf%20.10lf"%( k, f[iLo], df ) )
	log_function( "-" * 50 )
	obj.coor = x[iLo][:]
	obj.get_func()



def steepest_descent( obj, 
						step_number = 100,
						step_size = 0.1,
						print_frequency = 10,
						gradient_tolerance = 15.,
						log_function = default_log ):
	log_function( "---------------------------------------- Minimization (SD)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Step Size:          %20.10lg"%( step_size ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Gradient Tolerance: %20.10lg\n"%( gradient_tolerance ) )
	obj.get_grad()
	norm, grms = __grms( obj.grad )
	if( norm > step_size ):
		ssiz = step_size
	elif( norm > gradient_tolerance ):
		ssiz = norm
	else:
		ssiz = gradient_tolerance
	log_function( "%10s%20s%20s%20s"%( "Step", "Function", "Gradient", "Displacement" ) )
	log_function( "-" * 70 )
	log_function( "%10s%20.5lf%20.8lf%20.10lf"%( "", obj.func, grms, ssiz ) )
	i = 0
	while( i < step_number and grms > gradient_tolerance ):
		# -- perform step
		for j in range( obj.size ):
			obj.coor[j] -= obj.grad[j] / norm * ssiz
		# -- check new point
		obj.get_grad()
		norm, grms = __grms( obj.grad )
		fcur = round( obj.func, 0 )
		if( norm > step_size ):
			ssiz = step_size
		elif( norm > gradient_tolerance ):
			ssiz = norm
		else:
			ssiz = gradient_tolerance
		i = i + 1
		if( i%print_frequency == 0 ):
			log_function( "%10d%20.5lf%20.10lf%20.10lf"%( i, obj.func, grms, ssiz ) )
		obj.current_step( i )
	if( i%print_frequency != 0 ):
		log_function( "%10d%20.5lf%20.10lf%20.10lf"%( i + 1, obj.func, grms, ssiz ) )
	log_function( "-" * 70 + "\n" )



def adam( obj, 
				step_number = 100,
				step_size = 0.1,
				print_frequency = 10,
				gradient_tolerance = 15.,
				log_function = default_log ):
	log_function( "---------------------------------------- Minimization (ADAM)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Step Size:          %20.10lg"%( step_size ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Gradient Tolerance: %20.10lg\n"%( gradient_tolerance ) )
	obj.get_grad()
	norm, grms = __grms( obj.grad )
	beta = 0.9
	gamm = 0.999
	epsi = 1.e-8
	log_function( "%10s%20s%20s"%( "Step", "Function", "Gradient" ) )
	log_function( "-" * 50 )
	log_function( "%10s%20.5lf%20.8lf"%( "", obj.func, grms ) )
	v = [ 0.0 for i in range( obj.size ) ]
	s = [ 0.0 for i in range( obj.size ) ]
	i = 0
	while( i < step_number and grms > gradient_tolerance ):
		# -- perform step
		pbet = 1.0 / ( 1.0 - math.pow( beta, i + 1 ) )
		pgam = 1.0 / ( 1.0 - math.pow( gamm, i + 1 ) )
		for j in range( obj.size ):
			v[j] = beta * v[j] + ( 1.0 - beta ) * obj.grad[j]
			s[j] = gamm * s[j] + ( 1.0 - gamm ) * obj.grad[j] * obj.grad[j]
			obj.coor[j] -= step_size * v[j] * pbet / ( math.sqrt( s[j] * pgam ) + epsi )
		# -- check new point
		obj.get_grad()
		norm, grms = __grms( obj.grad )
		i = i + 1
		if( i%print_frequency == 0 ):
			log_function( "%10d%20.5lf%20.10lf"%( i, obj.func, grms ) )
		obj.current_step( i )
	if( i%print_frequency != 0 ):
		log_function( "%10d%20.5lf%20.10lf"%( i + 1, obj.func, grms ) )
	log_function( "-" * 50 )



def fire( obj,
			step_number = 100,
			step_size = 0.1,
			print_frequency = 10,
			gradient_tolerance = 1.5,
			log_function = default_log ):
	log_function( "---------------------------------------- Minimization (FIRE)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Step Size:          %20.10lg"%( step_size ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Gradient Tolerance: %20.10lg\n"%( gradient_tolerance ) )
	log_function( "%10s%20s%20s%20s"%( "Step", "Function", "Gradient", "Displacement" ) )
	log_function( "-" * 70 )
	nstp = 0
	ssiz = step_size
	alph = 0.1
	velo = [ 0.0 for i in range( obj.size ) ]
	step = [ 0.0 for i in range( obj.size ) ]
	obj.get_grad()
	norm, grms = __grms( obj.grad )
	log_function( "%10s%20.5lf%20.8lf%20.10lf"%( "", obj.func, grms, ssiz ) )
	i = 0
	while( i < step_number and grms > gradient_tolerance ):
		vsiz = math.sqrt( sum( [ velo[j] * velo[j] for j in range( obj.size ) ] ) )
		vfac = sum( [ - velo[j] * obj.grad[j] for j in range( obj.size ) ] )
		if( vfac > 0.0 ):
			velo = [ ( 1.0 - alph ) * velo[j] - alph * obj.grad[j] / norm * vsiz for j in range( obj.size ) ]
			if( nstp > 5 ):
				ssiz = min( ssiz * 1.1, step_size )
				alph *= 0.99
			nstp += 1
		else:
			velo = [ 0.0 for j in range( obj.size ) ]
			alph = 0.1
			ssiz *= 0.5
			nstp = 0

		for j in range( obj.size ):
			velo[j] -= ssiz * obj.grad[j]
			step[j] = ssiz * velo[j]
		tmp = math.sqrt( sum( [ step[j] * step[j] for j in range( obj.size ) ] ) )
		if( tmp > ssiz ):
			for j in range( obj.size ):
				step[j] *= ssiz / tmp
		for j in range( obj.size ):
			obj.coor[j] += step[j]

		obj.get_grad()
		norm, grms = __grms( obj.grad )

		i = i + 1
		if( i%print_frequency == 0 ):
			log_function( "%10d%20.5lf%20.10lf%20.10lf"%( i, obj.func, grms, ssiz ) )
		obj.current_step( i )
	if( i%print_frequency != 0 ):
		log_function( "%10d%20.5lf%20.10lf%20.10lf"%( i + 1, obj.func, grms, ssiz ) )
	log_function( "-" * 70 + "\n" )



#def conjugate_gradient( obj, 
#			step_number = 100,
#			step_size = 0.1,
#			print_frequency = 10,
#			gradient_tolerance = 1.5,
#			log_function = default_log ):
#	log_function( "---------------------------------------- Minimization (CG)\n" )
#	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
#	log_function( "Step Number:        %20d"%( step_number ) )
#	log_function( "Step Size:          %20.10lg"%( step_size ) )
#	log_function( "Print Frequency:    %20d"%( print_frequency ) )
#	log_function( "Gradient Tolerance: %20.10lg\n"%( gradient_tolerance ) )
#	log_function( "%10s%20s%20s"%( "Step", "Function", "Gradient" ) )
#	log_function( "-" * 50 )
# --------------------------------------
#	log_function( "-" * 50 + "\n" )



def l_bfgs( obj, 
			step_number = 100,
			step_size = 0.1,
			print_frequency = 10,
			gradient_tolerance = 1.5,
			history = 9,
			log_function = default_log ):
	log_function( "---------------------------------------- Minimization (L-BFGS)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Step Size:          %20.10lg"%( step_size ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Gradient Tolerance: %20.10lg"%( gradient_tolerance ) )
	log_function( "Number of Updates:  %20d\n"%( history ) )
	log_function( "%10s%20s%20s"%( "Step", "Function", "Gradient" ) )
	log_function( "-" * 50 )
	aux  = [ 0.0 for i in range( history ) ]
	rho  = [ 0.0 for i in range( history ) ]
	og   = [ 0.0 for i in range( obj.size ) ]
	ox   = [ 0.0 for i in range( obj.size ) ]
	step = [ 0.0 for i in range( obj.size ) ]
	dg   = []
	dx   = []
	for j in range( history ):
		dx.append( [ 0.0 for ii in range( obj.size ) ] )
		dg.append( [ 0.0 for ii in range( obj.size ) ] )
	obj.get_grad()
	norm, grms = __grms( obj.grad )
	log_function( "%10s%20.5lf%20.8lf"%( "", obj.func, grms ) )
	i = 0
	while( i < step_number and grms > gradient_tolerance ):
		if( i > history ):
			tmp =  dx.pop( 0 );  dx.append( tmp[:] )
			tmp =  dg.pop( 0 );  dg.append( tmp[:] )
			tmp = rho.pop( 0 ); rho.append( tmp    )
		if( i > 0 ):
			j   = min( i, history ) - 1
			hgx = 0.0
			hgg = 0.0
			for k in range( obj.size ):
				dx[j][k] = obj.coor[k] - ox[k]
				dg[j][k] = obj.grad[k] - og[k]
				hgx += dg[j][k] * dx[j][k]
				hgg += dg[j][k] * dg[j][k]
			rho[j] = 1.0 / hgx
			hscal  = hgx / hgg
		ox = obj.coor[:]
		og = obj.grad[:]
		step = [ -ii for ii in obj.grad ]
		if( i == 0 ):
			step = [ ii/norm for ii in step ]
		else:
			for j in range( min( i, history ) - 1, -1, -1 ):
				aux[j] = rho[j] * sum( [ ii*jj for ii,jj in zip( step, dx[j] ) ] )
				for k in range( obj.size ):
					step[k] -= aux[j] * dg[j][k]
			for k in range( obj.size ):
				step[k] *= hscal
			for j in range( min( i, history )  ):
				aux[j] -= rho[j] * sum( [ ii*jj for ii,jj in zip( step, dg[j] ) ] )
				for k in range( obj.size ):
					step[k] += aux[j] * dx[j][k]
		tmp = math.sqrt( sum( [ ii*ii for ii in step ] ) )
		if( tmp > step_size ):
			for j in range( obj.size ):
				step[j] *= step_size / tmp
		for j in range( obj.size ):
			obj.coor[j] += step[j]

		obj.get_grad()
		norm, grms = __grms( obj.grad )
		i = i + 1
		if( i%print_frequency == 0 ):
			log_function( "%10d%20.5lf%20.10lf"%( i, obj.func, grms ) )
		obj.current_step( i )
	if( i%print_frequency != 0 ):
		log_function( "%10d%20.5lf%20.10lf"%( i + 1, obj.func, grms ) )
	log_function( "-" * 50 + "\n" )



def conjugate_gradient_plus( obj, 
			step_number = 100,
			print_frequency = 10,
			gradient_tolerance = 1.5,
			method = "Polak-Ribiere", 
			restart = True,
			log_function = default_log ):
	log_function( "------------------------------------------ Minimization (CG+)\n" )
	log_function( "Degrees of Freedom:   %20ld"%( obj.size ) )
	log_function( "Step Number:          %20d"%( step_number ) )
	log_function( "Print Frequency:      %20d"%( print_frequency ) )
	log_function( "Gradient Tolerance:   %20.10lg"%( gradient_tolerance ) )
	log_function( "Method:             %22s\n"%( method ) )
	log_function( "%10s%20s%20s"%( "Step", "Function", "Gradient" ) )
	log_function( "-" * 50 )
	irest = int( restart )
	dmeth = { "Fletcher-Reeves" : 1, "Polak-Ribiere" : 2, "Positive Polak-Ribiere": 3 }
	if( method in dmeth ):
		imeth = dmeth[method]
	else:
		imeth = 3
	if( has_minimize_so ):
		qm3.actions._minimize.cgp( obj, step_number, gradient_tolerance, print_frequency, irest, imeth, log_function )
	else:
		raise Exception( "minimize.conjugate_gradient_plus: qm3.actions._minimize.so not available..." )
	log_function( "-" * 50 + "\n" )



def baker( obj, 
			step_number = 100,
			step_size = 0.1,
			print_frequency = 10,
			gradient_tolerance = 1.5,
			follow_mode = -1,		# -1 : minimum / 0... : mode following (TS)
			allow_overlap = False,
			log_function = default_log ):
	if( follow_mode >= obj.size or follow_mode < -1 ):
		follow_mode = -1
	log_function( "---------------------------------------- Minimization (Baker)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Following Mode:     %20d"%( follow_mode ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Step Size:          %20.10lg"%( step_size ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Gradient Tolerance: %20.10lg"%( gradient_tolerance ) )
	log_function( "Allow Overlap:      %20s\n"%( allow_overlap ) )
	if( follow_mode > -1 ):
		log_function( "%10s%20s%20s%20s"%( "Step", "Function", "Gradient", "Nneg,Fmode,Eval" ) )
		log_function( "-" * 70 )
	else:
		log_function( "%10s%20s%20s%5s"%( "Step", "Function", "Gradient", "Nneg" ) )
		log_function( "-" * 55 )
	mstp = 1.0e-1
	lrge = 1.0e+6
	step = 50.0
	tol1 = 1.0e-4
	tol2 = 1.0e-8
	emax = 1.0e5
	emin = 1.0e-3
	mxit = 999
	dx   = [ 0.0 for i in range( obj.size ) ]
	grms = gradient_tolerance * 2.0
	k    = 0
	flg  = True
	while( k < step_number and grms > gradient_tolerance and flg ):
		# update coordinates
		for i in range( obj.size ):
			obj.coor[i] += dx[i]

		# get into for the new point
		obj.get_hess()
		ei, ev = qm3.maths.matrix.diag( obj.hess, obj.size )

		# scale eigenvalues and take the number of negative modes...
		for i in range( obj.size ):
			if( math.fabs( ei[i] ) < emin  ):
				if( ei[i] < 0.0 ):
					ei[i] = - emin
				else:
					ei[i] = emin
			if( math.fabs( ei[i] ) > emax  ):
				if( ei[i] < 0.0 ):
					ei[i] = - emax
				else:
					ei[i] = emax
		nneg = sum( [ 1 for i in range( obj.size ) if ei[i] < 0.0 ] )

		# transform gradient vector to the local hessian modes, and startup dx
		dx   = []
		gx   = []
		grms = 0.0
		for i in range( obj.size ):
			grms += obj.grad[i] * obj.grad[i]
			dx.append( 0.0 )
			gx.append( 0.0 )
			for j in range( obj.size ):
				gx[i] += obj.grad[j] * ev[j*obj.size+i]
		grms = math.sqrt( grms )

		# check whether we are searching a specific mode or not
		if( follow_mode > -1 ):

			# check for changes in current mode via overlapping
			lowr = follow_mode;
			if( k > 0 and allow_overlap == 1 ):
				ovr = 0.0;
				for i in range( obj.size ):
					ovr += ev[i*obj.size+follow_mode] * mvec[i]
#				ovr = math.fabs( ovr )
				for j in range( obj.size ):
					if( j != follow_mode ):
						tmp = 0.0
						for i in range( obj.size ):
							tmp += ev[i*obj.size+j] * mvec[i]
#						tmp = math.fabs( tmp )
						if( tmp > ovr ):
							ovr = tmp
							follow_mode = j
				if( lowr != follow_mode ):
					log_function( "[Allow_Overlap] Selected following mode: %ld (%.6lf)"%( follow_mode, ovr ) )
			eiv  = ei[follow_mode];
			mvec = []
			for i in range( obj.size ):
				mvec.append( ev[i*obj.size+follow_mode] )

			# Calculate the step for the maximization
			if( math.fabs( gx[follow_mode] ) > tol1 ):
				tmp = 0.5 * ( ei[follow_mode] + math.sqrt( ei[follow_mode] * ei[follow_mode] + 4.0 * gx[follow_mode] * gx[follow_mode] ) ) 
				lmbd = gx[follow_mode] / ( tmp - ei[follow_mode] )
			else:
				if( nneg == 1 ):
					lmbd = - gx[follow_mode] / ei[follow_mode]
				else:
					lmbd = mstp
			for i in range( obj.size ):
				dx[i] = lmbd * ev[i*obj.size+follow_mode];

		# minimize along the selected modes (skip first if followed or not)
		if( follow_mode == 0 ):
			lowr = 1
		else:
			lowr = 0 
		lmbd = 0.0
		if( ei[lowr] < 0.0 ):
			lmbd = ei[lowr] - step
			l1   = ei[lowr]
			l2   = - lrge
		ovr  = 0.0;
		for j in range( obj.size ):
			if( j != follow_mode ):
				ovr += ( gx[j] * gx[j] ) / ( lmbd - ei[j] )
		i = 0
		while( i < mxit and math.fabs( lmbd - ovr ) >= tol2 ):
			if( ei[lowr] > 0.0 ):
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
			ovr  = 0.0;
			for j in range( obj.size ):
				if( j != follow_mode ):
					ovr += ( gx[j] * gx[j] ) / ( lmbd - ei[j] )
			i += 1
		if( i > mxit ):
			log_function( "\n -- Too much lambda iterations..." )
			flg = False

		# modify follow_mode eigenvalues and vectors...
		if( follow_mode > -1 ):
			ei[follow_mode] = lmbd - 1.0
			for i in range( obj.size ):
				ev[i*obj.size+follow_mode] = 0.0

		# check final step (may be too small or too large...)
		for i in range( obj.size ):
			dx[i] += sum( [ ev[i*obj.size+j] * gx[j] / ( lmbd - ei[j] ) for j in range( obj.size ) ] )
		ovr = math.sqrt( sum( [ dx[i] * dx[i] for i in range( obj.size ) ] ) )

		# checking for a small step
		if( ovr < tol2 ):
			log_function( "\n -- The step size is *very* small..." )
			flg = False

		# scale long steps...
		if( ovr > step_size ):
			for i in range( obj.size ):
				dx[i] *= step_size / ovr

		# next step...
		k    += 1
		grms /= math.sqrt( float( obj.size ) )

		# print something...
		if( k%print_frequency == 0 ):
			if( follow_mode < 0 ):
				log_function( "%10ld%20.5lf%20.10lf%5ld"%( k, obj.func, grms, nneg ) )
			else:
				log_function( "%10ld%20.5lf%20.10lf%5ld%5ld%10.2lf"%( k, obj.func, grms, nneg, follow_mode, eiv ) )
		obj.current_step( k )

	if( k%print_frequency != 0 ):
		if( follow_mode < 0 ):
			log_function( "%10ld%20.5lf%20.10lf%5ld"%( k, obj.func, grms, nneg ) )
		else:
			log_function( "%10ld%20.5lf%20.10lf%5ld%5ld%10.2lf"%( k, obj.func, grms, nneg, follow_mode, ei[follow_mode] ) )

	if( follow_mode > -1 ):
		log_function( "-" * 70 + "\n" )
	else:
		log_function( "-" * 55 )



def rfo( obj, 
			step_number = 100,
			step_size = 0.1,
			print_frequency = 10,
			gradient_tolerance = 1.5,
			follow_mode = -1,		# -1 : minimum / 0... : mode following (TS)
			log_function = default_log ):
	if( follow_mode >= obj.size or follow_mode < -1 ):
		follow_mode = -1
	log_function( "---------------------------------------- Minimization (RFO)\n" )
	log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
	log_function( "Following Mode:     %20d"%( follow_mode ) )
	log_function( "Step Number:        %20d"%( step_number ) )
	log_function( "Step Size:          %20.10lg"%( step_size ) )
	log_function( "Print Frequency:    %20d"%( print_frequency ) )
	log_function( "Gradient Tolerance: %20.10lg\n"%( gradient_tolerance ) )
	log_function( "%10s%20s%20s"%( "Step", "Function", "Gradient" ) )
	log_function( "-" * 50 )
	tol2 = 1.0e-8
	dx   = [ 0.0 for i in range( obj.size ) ]
	dd   = [ 0.5 for i in range( obj.size ) ]
	if( follow_mode > -1 ):
		dd[follow_mode] *= -1.0
	grms = gradient_tolerance * 2.0
	k    = 0
	flg  = True
	while( k < step_number and grms > gradient_tolerance and flg ):
		# update coordinates
		for i in range( obj.size ):
			obj.coor[i] -= dx[i]
	
		# get into for the new point
		obj.get_hess()
		ei, ev = qm3.maths.matrix.diag( obj.hess, obj.size )

		# calcualte new step
		vg = qm3.maths.matrix.mult( qm3.maths.matrix.T( ev, obj.size, obj.size ), obj.size, obj.size, obj.grad, obj.size, 1 )
		lg = [ dd[i] * ( math.fabs( ei[i] ) + math.sqrt( ei[i] * ei[i] + 4.0 * vg[i] * vg[i] ) ) for i in range( obj.size ) ]
		# just using lambda, and not the shifted version: vg[i] / ( ei[i] - lg[i] )
		# pure newton-raphson step would use only the corresponding eigenvalue (ei[i])
		dx = qm3.maths.matrix.mult( ev, obj.size, obj.size, [ vg[i] / lg[i] for i in range( obj.size ) ], obj.size, 1 )
		tt = math.sqrt( sum( [ dx[i] * dx[i] for i in range( obj.size ) ] ) )
	
		# checking for a small step
		if( tt < tol2 ):
			log_function( "\n -- The step size is *very* small..." )
			flg = False
	
		# scale long steps...
		if( tt > step_size ):
			for i in range( obj.size ):
				dx[i] *= step_size / tt
	
		# next step...
		k   += 1
		grms = math.sqrt( sum( [ obj.grad[i] * obj.grad[i] for i in range( obj.size ) ] ) / float( obj.size ) )
	
		# print something...
		if( k%print_frequency == 0 ):
			log_function( "%10ld%20.5lf%20.10lf"%( k, obj.func, grms ) )
		obj.current_step( k )
	
	if( k%print_frequency != 0 ):
		log_function( "%10ld%20.5lf%20.10lf"%( k, obj.func, grms ) )

	log_function( "-" * 50 )



###############################################################################
# Iterable version of the minimizers
#
class stepped_fire( object ):

	def __init__( self, obj, step_size = 0.1, print_frequency = 10, log_function = default_log ):

		self.obj  = obj
		self.nstp = 0
		self.ssiz = step_size
		self.alph = 0.1
		self.velo = [ 0.0 for i in range( obj.size ) ]
		self.step = [ 0.0 for i in range( obj.size ) ]
		self.log_function = log_function
		self.print_frequency = print_frequency
		self.step_size = step_size

		self.log_function( "---------------------------------------- Minimization (FIRE)\n" )
		self.log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
		self.log_function( "Step Size:          %20.10lg"%( step_size ) )
		self.log_function( "Print Frequency:    %20d\n"%( print_frequency ) )
		self.log_function( "%10s%20s%20s"%( "Step", "Function", "Gradient" ) )
		self.log_function( "-" * 50 )

		self.obj.get_grad()
		self.norm = math.sqrt( sum( [ i*i for i in self.obj.grad ] ) )
		self.grms = self.norm / math.sqrt( self.obj.size )
		self.log_function( "%10s%20.5lf%20.8lf"%( "", self.obj.func, self.grms  ) )
		self.i = 0


	def iterate( self ):
		vsiz = math.sqrt( sum( [ self.velo[j] * self.velo[j] for j in range( self.obj.size ) ] ) )
		vfac = sum( [ - self.velo[j] * self.obj.grad[j] for j in range( self.obj.size ) ] )
		if( vfac > 0.0 ):
			self.velo = [ ( 1.0 - self.alph ) * self.velo[j] - self.alph * self.obj.grad[j] / self.norm * vsiz for j in range( self.obj.size ) ]
			if( self.nstp > 5 ):
				self.ssiz = min( self.ssiz * 1.1, self.step_size )
				self.alph *= 0.99
			self.nstp += 1
		else:
			self.velo = [ 0.0 for j in range( self.obj.size ) ]
			self.alph = 0.1
			self.ssiz *= 0.5
			self.nstp = 0

		for j in range( self.obj.size ):
			self.velo[j] -= self.ssiz * self.obj.grad[j]
			self.step[j]  = self.ssiz * self.velo[j]
		tmp = math.sqrt( sum( [ self.step[j] * self.step[j] for j in range( self.obj.size ) ] ) )
		if( tmp > self.ssiz ):
			for j in range( self.obj.size ):
				self.step[j] *= self.ssiz / tmp
		for j in range( self.obj.size ):
			self.obj.coor[j] += self.step[j]

		self.obj.get_grad()
		self.norm = math.sqrt( sum( [ i*i for i in self.obj.grad ] ) )
		self.grms = self.norm / math.sqrt( self.obj.size )
	
		self.i += 1
		if( self.i%self.print_frequency == 0 ):
			self.log_function( "%10d%20.5lf%20.10lf"%( self.i, self.obj.func, self.grms ) )
		self.obj.current_step( self.i )

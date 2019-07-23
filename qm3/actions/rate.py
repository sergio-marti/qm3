# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.constants
import	qm3.maths.matrix
import	qm3.maths.roots
import	qm3.maths.integration
import	qm3.maths.interpolation



# mass weighted:  xyz * sqrt(m)  ;  grd / sqrt(m)  ;  hes / sqrt(mi * mj)
def __project( mas, crd, grd, hes, prjgrad = False ):
	siz = len( crd )
	mtt = 0.0
	cen = [ 0.0, 0.0, 0.0 ]
	for i in range( siz // 3 ):
		k = i * 3
		mtt += mas[k] * mas[k]
		for j in [0, 1, 2]:
			cen[j] += crd[k+j] * mas[k]
	cen = [ cen[i] / mtt for i in [0, 1, 2] ]
	rtm = [ 0.0 for i in range( 6 * siz ) ]
	for i in range( siz // 3 ):
		k              = i * 3
		rtm[k]	       = mas[k]
		rtm[siz+k+1]   = mas[k]
		rtm[2*siz+k+2] = mas[k]
		rtm[3*siz+k+1] = - ( crd[k+2] - cen[2] * mas[k] )
		rtm[3*siz+k+2] =   ( crd[k+1] - cen[1] * mas[k] )
		rtm[4*siz+k]   =   ( crd[k+2] - cen[2] * mas[k] )
		rtm[4*siz+k+2] = - ( crd[k  ] - cen[0] * mas[k] )
		rtm[5*siz+k]   = - ( crd[k+1] - cen[1] * mas[k] )
		rtm[5*siz+k+1] =   ( crd[k  ] - cen[0] * mas[k] )
	for i in range( 6 ):
		for j in range( i ):
			tmp = sum( [ rtm[i*siz+k] * rtm[j*siz+k] for k in range( siz ) ] )
			for k in range( siz ):
				rtm[i*siz+k] -= tmp * rtm[j*siz+k]
		tmp = math.sqrt( sum( [ rtm[i*siz+k] * rtm[i*siz+k] for k in range( siz ) ] ) )
		for k in range( siz ):
			rtm[i*siz+k] /= tmp
	# gradient
	for i in range( 6 ):
		tmp = sum( [ grd[k] * rtm[i*siz+k] for k in range( siz ) ] )
		for k in range( siz ):
			grd[k] -= tmp * rtm[i*siz+k]
	# hessian
	ixx = [ 0.0 for i in range( siz * siz ) ]
	for i in range( siz ):
		ixx[siz*i+i] += 1.
		for j in range( siz ):
			for k in range( 6 ):
				ixx[siz*i+j] -= rtm[k*siz+i] * rtm[k*siz+j]
	tmp = qm3.maths.matrix.mult( ixx, siz, siz, qm3.maths.matrix.mult( hes, siz, siz, ixx, siz, siz ), siz, siz )
	for i in range( siz * siz ):
		hes[i] = tmp[i]
	val, vec = qm3.maths.matrix.diag( tmp, siz )
	# -- remove the reaction coordinate from the frequencies...
	c = 1.0e13 / ( 2.0 * math.pi )
	t = []
	for i in range( siz ):
		if( val[i] < 0.0 ):
			t.append( - math.sqrt( math.fabs( val[i] ) ) * c )
		else:
			t.append(   math.sqrt( math.fabs( val[i] ) ) * c )
	skp  = sum( [ 1 for i in val if i < 1 ] )
	if( prjgrad ):
		tmp = sum( [ i * i for i in grd ] )
		ixx = [ 0.0 for i in range( siz * siz ) ]
		for i in range( siz ):
			ixx[siz*i+i] += 1.
			for j in range( siz ):
				ixx[siz*i+j] -= grd[i] * grd[j] / tmp
		tmp = qm3.maths.matrix.mult( ixx, siz, siz, qm3.maths.matrix.mult( hes, siz, siz, ixx, siz, siz ), siz, siz )
		val, tmp = qm3.maths.matrix.diag( tmp, siz )
	return( skp, val, vec )



def path_write( fd, obj ):
	# -- func
	fd.write( "%16.8le\n"%( obj.func ) )
	# -- coor
	for i in range( obj.size // 3 ):
		i3 = i * 3
		for j in [0, 1, 2]:
			fd.write( "%16.8le"%( obj.coor[i3+j] ) )
		fd.write( "\n" )
	# -- grad
	for i in range( obj.size // 3 ):
		i3 = i * 3
		for j in [0, 1, 2]:
			fd.write( "%16.8le"%( obj.grad[i3+j] ) )
		fd.write( "\n" )
	# -- hess (upper diagonal by rows)
	k = 1
	for i in range( obj.size ):
		ix = i * obj.size
		for j in range( i, obj.size ):
			fd.write( "%16.8le"%( obj.hess[ix+j] ) )
			if( k%6 == 0 ):
				fd.write( "\n" )
			k += 1
	if( k%6 != 0 ):
		fd.write( "\n" )
	fd.flush()



def path_read( fd, obj ):
	nh = obj.size * ( obj.size + 1 ) // 2
	try:
		func = float( fd.readline().strip() )
		coor = []
		for i in range( obj.size // 3 ):
			coor += [ float( j ) for j in fd.readline().strip().split() ]
		grad = []
		for i in range( obj.size // 3 ):
			grad += [ float( j ) for j in fd.readline().strip().split() ]
		temp = []
		while( len( temp ) < nh ):
			temp += [ float( j ) for j in fd.readline().strip().split() ]
		hess = qm3.maths.matrix.from_upper_diagonal_rows( temp, obj.size )
		# ---------------------------------------------------------
		obj.func = func
		k = 0
		for i in range( obj.size ):
			obj.coor[i] = coor[i]
			obj.grad[i] = grad[i]
			for j in range( obj.size ):
				obj.hess[k] = hess[k]
				k += 1
		return( True )
	except:
		return( False )



def curvature( obj, prjgrad = False ):
	"""
		frq			cm^-1
		zpe			kJ/mol		
		eta			1 / [ ang * sqrt( g/mol ) ]
		tau			ang * sqrt( g/mol )
	"""
	w = []
	for i in range( obj.size // 3 ):
		t = math.sqrt( obj.mass[i] )
		w += [ t, t, t ]
	x = [ obj.coor[i] * w[i] for i in range( obj.size ) ]
	g = [ obj.grad[i] / w[i] for i in range( obj.size ) ]
	h = []
	k = 0
	for i in range( obj.size ):
		for j in range( obj.size ):
			h.append( obj.hess[k] / ( w[i] * w[j] ) )
			k += 1
	neg, val, vec = __project( w, x, g, h, prjgrad )
	c = 1.0e13 / ( 2.0 * math.pi )
	for i in range( obj.size ):
		if( val[i] < 0.0 ):
			val[i] = - math.sqrt( math.fabs( val[i] ) ) * c
		else:
			val[i] =   math.sqrt( math.fabs( val[i] ) ) * c
	# -----------------------------------------------------------------------
	wcn  = 100.0 * qm3.constants.C
	skp  = sum( [ 1 for i in val if i < wcn ] )
	zpe  = 0.5 * qm3.constants.H * 1.0e-3 * qm3.constants.NA * sum( val[skp:] )
	try:
		G    = math.sqrt( sum( [ i * i for i in g ] ) )
		dgds = [ i / G for i in qm3.maths.matrix.mult( h, obj.size, obj.size, [ -j / G for j in g ], obj.size, 1 ) ]
		Bm   = []
		for i in range( skp, obj.size ):
			Bm.append( sum( [ vec[i+obj.size*j] * dgds[j] for j in range( obj.size ) ] ) )
		eta  = math.sqrt( sum( [ i * i for i in Bm ] ) )
		Bmw2 = math.pow( sum( [ i*i * j*j for i,j in zip( Bm, val[skp:] ) ] ), 0.25 )
		tau  = math.sqrt( eta * qm3.constants.H / ( 2.0 * math.pi ) * 1.0e23 * qm3.constants.NA ) / Bmw2
	except:
		eta  = None
		tau  = None
	return( neg, [ i / wcn for i in val ], zpe, eta, tau )



def effective_reduced_masses( s_coor, eta, tau ):
	i_tau = qm3.maths.interpolation.hermite_spline( s_coor, tau )
	erm = []
	for i in range( len( s_coor ) ):
		tmp = i_tau.calc( s_coor[i] )[1]
		tmp = - 2.0 * eta[i] * tau[i] - eta[i] * tau[i] * eta[i] * tau[i] + tmp * tmp
		if( tmp < 0.0 ):
			erm.append( min( 1.0, math.exp( tmp ) ) )
		else:
			erm.append( 1.0 )
	print()
	print( "%20s%20s"%( "s ", "mu_cd-sc " ) )
	print( 40 * "-" )
	for i in range( len( s_coor ) ):
		print( "%20.10lf%20.10lf"%( s_coor[i], erm[i] ) )
	return( erm )



def transmission_coefficient( s_coor, v_adia, r_mass = 1.0, temp = 298.15 ):
	"""
		s_coor		ang * sqrt( g/mol )
		v_adia		adiabatic energy (V+ZPE) relative to reactants potential (V) in kJ/mol
		r_mass		relative reduced masses (float: ZCT / list: SCT)
	"""
	# ------------------------------------------------------------------------------------
	def __turning_points( s_reac, s_prod, s_max, vad, ecr ):
		if( vad.calc( s_reac )[0] > ecr ):
			s_lt = s_reac
		else:
			s_lt = qm3.maths.roots.bisect( lambda x: vad.calc( x )[0] - ecr, s_reac, s_max, max_iter = 10000 )
		if( vad.calc( s_prod )[0] > ecr ):
			s_gt = s_prod
		else:
			s_gt = qm3.maths.roots.bisect( lambda x: vad.calc( x )[0] - ecr, s_max, s_prod, max_iter = 10000 )
		sys.stdout.flush()
		return( s_lt, s_gt )
	# ------------------------------------------------------------------------------------
	# -- maximum of the adiabatic potential ----------------------------------------------
	who = s_coor.index( 0.0 )
	vad  = qm3.maths.interpolation.hermite_spline( s_coor, v_adia )
	d    = 1.e-4
	smax = max( -0.1, s_coor[0] )
	f, g = vad.calc( smax )
	h    = - math.fabs( ( vad.calc( smax + d )[1] - vad.calc( smax - d )[1] ) / ( 2.0 * d ) )
	i    = 0
	while( i < 10000 and math.fabs( g ) > 1.0e-6 ):
		ds   = g / h
		smax -= ds / math.fabs( ds ) * min( math.fabs( ds ), 0.001 )
		f, g = vad.calc( smax )
		h    = - math.fabs( ( vad.calc( smax + d )[1] - vad.calc( smax - d )[1] ) / ( 2.0 * d ) )
		i += 1
	print()
	if( i < 10000 ):
		vmax = vad.calc( smax )[0]
		print( "V_max(s_max: %.4lf): %.4lf kJ/mol [spline]"%( smax, vmax ) )
	else:
		# try quadratic interpolation
		tmp = qm3.maths.matrix.mult( 
				qm3.maths.matrix.inverse( [ 1.0, s_coor[who-1], s_coor[who-1]*s_coor[who-1], 1.0, s_coor[who], s_coor[who]*s_coor[who], 1.0, s_coor[who+1], s_coor[who+1]*s_coor[who+1] ], 3, 3 ), 3, 3, 
				[ v_adia[who-1], v_adia[who], v_adia[who+1] ], 3, 1 )
		smax = - tmp[1] / ( 2.0 * tmp[2] )
		vmax = vad.calc( smax )[0]
		print( "V_max(s_max: %.4lf): %.4lf kJ/mol [quadint]"%( smax, vmax ) )
	# -- transmission probabilities ------------------------------------------------------
	c   = 2.0e-10 * math.pi / ( qm3.constants.H * qm3.constants.NA )
	emx = max( v_adia[0], v_adia[-1] ) + 1.0
#	npt = 100
#	eds = ( vmax - emx ) / ( npt - 1 )
	npt = len( qm3.maths.integration.gauss_legendre_xi )
	eds = vmax - emx
	print()
	print( "%16s%16s%10s%10s"%( "V_adi ", "Probability ", "s< ", "s> " ) )
	print( 56 * "-" )
	if( type( r_mass ) == float ):
		mef = qm3.maths.interpolation.hermite_spline( s_coor, [ r_mass for i in range( len( s_coor ) ) ] )
	else:
		mef = qm3.maths.interpolation.hermite_spline( s_coor, r_mass )
	ener = []
	prob = []
	for i in range( npt ):
#		ecr = emx + i * eds
		ecr = emx + 0.5 * ( 1.0 + qm3.maths.integration.gauss_legendre_xi[i] ) * eds
		ener.append( ecr )
		s_lt, s_gt = __turning_points( s_coor[0], s_coor[-1], smax, vad, ecr )
		tmp = c * qm3.maths.integration.Gauss( lambda x: math.sqrt( 2.0 * math.fabs( mef.calc( x )[0] ) * math.fabs( vad.calc( x )[0] - ecr ) ), s_lt, s_gt )
		prob.append( 1.0 / ( 1.0 + math.exp( 2.0 * tmp ) ) )
		print( "%16.3lf%16.4le%10.4lf%10.4lf"%( ecr, prob[-1], s_lt, s_gt ) )
	# -- transmission coefficient ------------------------------------------------------
	r = 1.0e3 / ( qm3.constants.R * temp )
	e = qm3.maths.interpolation.hermite_spline( ener, [ i * math.exp( - j * r ) for i,j in zip( prob, ener ) ] )
	K = 1.0 + r * qm3.maths.integration.Gauss( lambda x: e.calc( x )[0], ener[0], ener[-1] ) * math.exp( ener[-1] * r )
	print()
	if( type( r_mass ) == float ):
		print( "Kappa_ZCT: %.4lf (%.2lf)"%( K, temp ) )
	else:
		print( "Kappa_SCT: %.4lf (%.2lf)"%( K, temp ) )
	return( K, ener, prob )


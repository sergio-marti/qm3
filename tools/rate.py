# -*- coding: iso-8859-1 -*-
from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.elements
import	qm3.constants
import	qm3.maths.matrix
import	qm3.maths.roots
import	qm3.maths.integration
import	qm3.maths.interpolation



################################################################################################################
def __projections( size, mass, coor, grad, hess ): 
	mc = [ 0.0, 0.0, 0.0 ]
	t = 0.0
	for i in range( size // 3 ):
		t += ( mass[3*i] * mass[3*i] )
		for j in [0, 1, 2]:
			mc[j] += coor[3*i+j] * mass[3*i]
	mc = [ i / t for i in mc ]
	rt = []
	for i in range( 6 ):
		rt.append( [ 0.0 for i in range( size ) ] )
	for i in range( size // 3 ):
		j = 3 * i
		rt[0][j:j+3] = [ mass[j], 0.0, 0.0 ]
		rt[1][j:j+3] = [ 0.0, mass[j], 0.0 ]
		rt[2][j:j+3] = [ 0.0, 0.0, mass[j] ]
		rt[3][j:j+3] = [ 0.0, - ( coor[j+2] - mc[2] * mass[j] ), ( coor[j+1] - mc[1] * mass[j] ) ]
		rt[4][j:j+3] = [ ( coor[j+2] - mc[2] * mass[j] ), 0.0, - ( coor[j  ] - mc[0] * mass[j] ) ]
		rt[5][j:j+3] = [ - ( coor[j+1] - mc[1] * mass[j] ), ( coor[j  ] - mc[0] * mass[j] ), 0.0 ]
	for i in range( 6 ):
		for j in range( i ):
			t = sum( [ ii * jj for ii,jj in zip( rt[i], rt[j] ) ] )
			for k in range( size ):
				rt[i][k] -= t * rt[j][k]
		t = math.sqrt( sum( [ ii * ii for ii in rt[i] ] ) )
		for k in range( size ):
			rt[i][k] /= t
	# G' = G - Tx * G * Tx - ... - Rx * G * Rx - ...
	if( grad ):
		for i in range( 6 ):
			t = sum( [ ii * jj  for ii,jj in zip( grad, rt[i] ) ] )
			for j in range( size ):
				grad[j] -= t * rt[i][j]
	# P = I - Tx * Tx - ... - Rx * Rx - ... - G * G
	ix = [ 0.0 for ii in range( size * size ) ]
	for i in range( size ):
		for j in range( size ):
			for l in range( 6 ):
				ix[size*i+j] -= rt[l][i] * rt[l][j]
	if( grad ):
		t = sum( [ grad[i] * grad[i] for i in range( size ) ] )
		for i in range( size ):
			for j in range( size ):
				ix[size*i+j] -= grad[i] * grad[j] / t
	for i in range( size ):
		ix[size*i+i] += 1.0
	# H' = P * H * P
	hx = qm3.maths.matrix.mult( ix, size, size, qm3.maths.matrix.mult( hess, size, size, ix, size, size ), size, size )
	for i in range( size * size ):
		hess[i] = hx[i]
	v, V = qm3.maths.matrix.diag( hx, size )
	t = 1.0e13 / ( 2.0 * math.pi )
	for i in range( size ):
		if( v[i] < 0.0 ):
			v[i] = - math.sqrt( math.fabs( v[i] ) ) * t
		else:
			v[i] =   math.sqrt( math.fabs( v[i] ) ) * t
	return( v, V )



def __mep_analysis( f, d_, e0_, x_, w_ ):
	"""
	B_m ( s ) = 1 over lline g(s) rline L^T_m(s) cdot { left[ H(s) cdot { {-g(s)} over lline g(s) rline } right]  } ~~~~ %eta (s) = left(  sum_{m}^{3N-7}{ B_m^2 (s) } right)^{ 1 over 2 }
	newline
	t(s) = %eta^{ 1 over 2 }(s) left( %Ux210F over %my right)^{ 1 over 2 } left[ sum_m^{3N-7}{ B_m^2(s) %omega_m^2(s)} right]^{ - {1 over 4} } ~~~~ a( s ) = t( s ) %eta( s )
	"""
	wcm = 100.0 * qm3.constants.C
	n3_ = len( w_ )
	n_  = n3_ // 3
	nh_ = n3_ * ( n3_ + 1 ) // 2
	o = []
	x = x_[:]
	a = 0.0
	e = f.readline()
	while( e and e[0] != "*" ):
		y = x[:]
		e = float( e )
		x = []
		for i in range( n_ ):
			x += [ float( j ) * w_[3*i] for j in f.readline().split() ]
		a += d_ * math.sqrt( sum( [ (ii-jj)*( ii-jj) for ii,jj in zip( x, y ) ] ) )
		g = []
		for i in range( n_ ):
			g += [ float( j ) / w_[3*i] for j in f.readline().split() ]
		h = []
		while( len( h ) < nh_ ):
			h += [ float( j ) for j in f.readline().split() ]

		# >> cambiar esto por algo mas "QM3"
		h = qm3.maths.matrix.from_upper_diagonal_columns( h, n3_ )

		for i in range( n3_ ):
			for j in range( n3_ ):
				h[i*n3_+j] /= ( w_[i] * w_[j] )
		v, V = __projections( n3_, w_, x, g, h )
		nn = sum( [ 1 for i in v if i < wcm ] )
#		if( nn == 7 ):
#		if( nn > 6 ):
		if( True ):
			z = 0.5 * qm3.constants.H * 1.0e-3 * qm3.constants.NA * sum( v[nn:] )
			# curvature analysis  ----------------------------------------------------------------------------------------------
			G = math.sqrt( sum( [ i * i for i in g ] ) )
			dgds = [ i / G for i in qm3.maths.matrix.mult( h, n3_, n3_, [ -i / G for i in g ], n3_, 1 ) ]
			Bm = []
			for i in range( nn, n3_ ):
				Bm.append( sum( [ V[i+n3_*j] * dgds[j] for j in range( n3_ ) ] ) )
			eta  = math.sqrt( sum( [ i * i for i in Bm ] ) )
			BmW2 = math.pow( sum( [ i*i * j*j for i,j in zip( Bm, v[nn:] ) ] ), 0.25 )
			t_s  = math.sqrt( eta * qm3.constants.H / ( 2.0 * math.pi ) ) / BmW2 * 1.0e10 * math.sqrt( 1.0e3 * qm3.constants.NA )
			# ------------------------------------------------------------------------------------------------------------------
			o.append( ( a, e - e0_, z, t_s, eta ) )
			print( "%10.4lf%16.3lf%16.3lf%16.3lf  %8.1lf%8.1lf%8.1lf%8.1lf%8.1lf%8.1lf%8.1lf%8.1lf"%( 
				a, ( e - e0_), z, ( e - e0_ + z ),
				v[0] / wcm, v[1] / wcm, v[2] / wcm, v[3] / wcm, v[4] / wcm, v[5] / wcm, v[6] / wcm, v[7] / wcm ) )
		e = f.readline()
	return( o )



def parse_MEP( fname, temperature = 298.15, isotopes = [] ):
	"""
	MEP structure:

	Z		atomic number
	*reac
	vpot	kJ/mol
	coor	A
	hess 	kJ/(mol*A^2)
	*sadd
	vpot
	coor
	hess
	*>>reac
	vpot
	coor
	grad    kJ/(mol*A)
	[hess]
	...
	*>>prod
	vpot
	coor
	grad
	[hess]
	...

	returns:
		ZPE(r), s_coordinate, V(i)-V(reac), ZPE(i), relative centrifugal-dominant masses [t(s),eta(s)]
	"""
	print( 80 * "=" )
	print( " s in (g/mol)^0.5 * A\n dE, ZPE, dVa in kJ/mol\n mu_eff is adimensional (mu_eff/mu)\n Frequencies in cm^-1" )
	print( 80 * "-", "\n" )
	f = open( fname, "rt" )
	# masses and system size
	n_ = 0
	w_ = []
	for i in f.readline().split():
		t = math.sqrt( qm3.elements.mass[int(i)] )
		w_.append( t ); w_.append( t ); w_.append( t )
		n_ += 1
	for i,j in isotopes:
		t = math.sqrt( j )
		w_[3*i]   = t; w_[3*i+1] = t; w_[3*i+2] = t
	s_ = []
	# =========================================================================
	# Reactant
	f.readline()
	e0_ = float( f.readline() )
	x = []
	for i in range( n_ ):
		x += [ float( j ) * w_[3*i] for j in f.readline().split() ]
	h = []
	n3_ = n_ * 3
	nh_ = n3_ * ( n3_ + 1 ) // 2
	while( len( h ) < nh_ ):
		h += [ float( j ) for j in f.readline().split() ]

	# >> cambiar esto por algo mas "QM3"
	h = qm3.maths.matrix.from_upper_diagonal_columns( h, n3_ )

	for i in range( n3_ ):
		for j in range( n3_ ):
			h[i*n3_+j] /= ( w_[i] * w_[j] )
	v = __projections( n3_, w_, x, None, h )[0]
	print( "Reactant frequencies:" )
	for i in range( n3_ ):
		print( "%12.2lf"%( v[i] / ( 100.0 * qm3.constants.C ) ), end = "" )
		if( (i+1)%6 == 0 ):
			print()
	z0_ = 0.5 * qm3.constants.H * 1.0e-3 * qm3.constants.NA * sum( [ i for i in v[6:] if i > 1.0 ] )
	print()
	print( "Reactant ZPE: %16.4lf"%( z0_ ) )
	print()
	# =========================================================================
	f.readline()
	# TS
	e = float( f.readline() )
	x_ = []
	for i in range( n_ ):
		x_ += [ float( j ) * w_[3*i] for j in f.readline().split() ]
	h = []
	while( len( h ) < nh_ ):
		h += [ float( j ) for j in f.readline().split() ]

	# >> cambiar esto por algo mas "QM3"
	h = qm3.maths.matrix.from_upper_diagonal_columns( h, n3_ )

	for i in range( n3_ ):
		for j in range( n3_ ):
			h[i*n3_+j] /= ( w_[i] * w_[j] )
	v = __projections( n3_, w_, x_, None, h )[0]
	print( "TS frequencies:" )
	for i in range( n3_ ):
		print( "%12.2lf"%( v[i] / ( 100.0 * qm3.constants.C ) ), end = "" )
		if( (i+1)%6 == 0 ):
			print()
	z = 0.5 * qm3.constants.H * 1.0e-3 * qm3.constants.NA * sum( [ i for i in v[7:] if i > 1.0 ] )
	print()
	print( "TS ZPE: %16.4lf"%( z ) )
	print( "TS dE : %16.4lf"%( e - e0_ ) )
	print( "TS dVa: %16.4lf"%( e - e0_ + z ) )
	print( "WIGNER: %16.4lf"%( 1.0 + 1.0 / 24.0 * math.pow( v[0] * qm3.constants.H / ( qm3.constants.KB * temperature ), 2.0 ) ) )
	print()
	# =========================================================================
	# IRC
	print( "%10s%16s%16s%16s  %64s"%( "s ", "dE ", "ZPE ", "dVa ", "Frequencies " ) )
	print( 124 * "-" )
	d_ = 1.0
	if( f.readline().lower().strip() == "*>>reac" ):
		d_ = -1.0
	s_ += __mep_analysis( f,  d_, e0_, x_, w_ )
	s_ += __mep_analysis( f, -d_, e0_, x_, w_ )
	f.close()
	# =========================================================================
	s_.sort()
	rs = []; re = []; rz = []; rt = []; rk = []
	for i,j,k,l,m in s_:
		rs.append( i )
		re.append( j )
		rz.append( k )
		rt.append( l )
		rk.append( m )
	i = 0
	while( rs[i] < 0.0 ):
		i += 1
	rs.insert( i, 0.0 )
	re.insert( i, e - e0_ )
	rz.insert( i, z )
	# fit of t(s) and k(s) at the transition state
	rt.insert( i, 0.5 * ( rt[i] + rt[i-1] ) )
	rk.insert( i, 0.5 * ( rk[i] + rk[i-1] ) )
	print()
	print( "%20s%20s%20s%20s%20s"%( "s ", "dE ", "ZPE ", "t ", "k " ) )
	print( 100 * "-" )
	for i in range( len( rs ) ):
		print( "%20.10lf%20.10lf%20.10lf%20.10lf%20.10lf"%( rs[i], re[i], rz[i], rt[i], rk[i] ) )
	print()
	return( z0_, rs, re, rz, rt, rk )



################################################################################################################
def effective_reduced_masses( s_coor, m_t, m_k ):
	"""
	%my_eff(s) = %my cdot min left lbrace stack{ e^{- 2 a(s) - left[ a(s) right]^2 + left( dt over ds right)^2 } # 1} right none ~~~~ %my_{ a.u. } = 1 over {  N_A m_e 10^3 }
	"""
	dtds = qm3.maths.interpolation.hermite_spline( s_coor, m_t, "akima" )
	rm = []
	for i in range( len( s_coor ) ):
		t = dtds.calc( s_coor[i] )[1]
		t = - 2.0 * m_t[i] * m_k[i] - m_t[i] * m_t[i] * m_k[i] * m_k[i] + t * t
		if( t < 0.0 ):
			rm.append( min( 1.0, math.exp( t ) ) )
		else:
			rm.append( 1.0 )
	print()
	print( "%20s%20s"%( "s ", "mu_cd-sc " ) )
	print( 40 * "-" )
	for i in range( len( s_coor ) ):
		print( "%20.10lf%20.10lf"%( s_coor[i], rm[i] ) )
	print()
	return( rm )



################################################################################################################
def transmission_probabilities( s_coor, v_adia, s_maxi, v_maxi, r_zpe, r_mass ):
	"""
	%theta left( E right) = %Ux210F ^{ {}-1 } int_{ s<{} }^{ s>{} }{ left(2 %my_eff left[ V_a left( s right)- E right] right)^{1 over 2} ds }
	newline
	P left( E right) = left(  1 + e^{2 cdot %theta left( E right)} right)^{ {}-1 } ~~~~ E <= V_a^{s=0}

		s_coor			in amu^0.5*A
		v_adia			adiabatic energy (V+ZPE) relative to reactants potential (V) in kJ/mol
		v_maxi			maximum of the adiabatic energy (V+ZPE) relative to reactants potential (V) in kJ/mol
		r_zpe			ZPE of the reactants in kJ/mol
		r_mass			scaling factor (float: ZCT) or relative reduced masses (list: SCT)
	"""
	def __get_turning_points( s_coor, v_adia, s_maxi, e_cr, vag ):
# -------------------------------------------------------------------------------------------
#		t = [ i - e_cr for i in v_adia ]
#		T = t[:]; T.sort()
#		if( t[0] > 0.0 ):
#			s_lt = s_coor[0]
#		else:
#			if( T[-1] >= 0.0 ):
#				i = 0
#				while( t[i] < 0.0 ):
#					i += 1
#				s_lt = qm3.maths.roots.bisect( lambda x: aki_vag.calc( x )[0] - e_cr, s_coor[i-1], s_coor[i] )
#			else:
#				s_lt = qm3.maths.roots.bisect( lambda x: aki_vag.calc( x )[0] - e_cr, s_coor[0], s_maxi,  max_iter = 10000 )
#		if( t[-1] > 0.0 ):
#			s_gt = s_coor[-1]
#		else:
#			if( T[-1] >= 0.0 ):
#				i = -1
#				while( t[i] < 0.0 ):
#					i -= 1
#				s_gt = qm3.maths.roots.bisect( lambda x: aki_vag.calc( x )[0] - e_cr, s_coor[i], s_coor[i+1] )
#			else:
#				s_gt = qm3.maths.roots.bisect( lambda x: aki_vag.calc( x )[0] - e_cr, s_maxi, s_coor[-1], max_iter = 10000 )
# -------------------------------------------------------------------------------------------
		if( vag.calc( s_coor[0] )[0] > e_cr ):
			s_lt = s_coor[0]
		else:
			s_lt = qm3.maths.roots.bisect( lambda x: vag.calc( x )[0] - e_cr, s_coor[0], s_maxi,  max_iter = 10000 )
		if( vag.calc( s_coor[-1] )[0] > e_cr ):
			s_gt = s_coor[-1]
		else:
			s_gt = qm3.maths.roots.bisect( lambda x: vag.calc( x )[0] - e_cr, s_maxi, s_coor[-1], max_iter = 10000 )
# -------------------------------------------------------------------------------------------
		return( s_lt, s_gt )
	npts = len( qm3.maths.integration.gauss_legendre_xi )
	ener = []
	prob = []
	# -- setup masses
	if( type( r_mass ) == list ):
		erm_ = [ i for i in r_mass ]
	else:
		erm_ = [ 1.0 for i in range( len( s_coor ) ) ]
	# -- integrate probabilities
	c = 2.0e-10 * math.pi / ( qm3.constants.H * qm3.constants.NA )
# -------------------------------------- #
#	e_ds = ( v_maxi - r_zpe )
# -------------------------------------- #
	e_mx = max( v_adia[0], v_adia[-1] )
	e_ds = v_maxi - e_mx
# -------------------------------------- #
	print()
	print( "%16s%20s%10s%10s"%( "dVa", "Trans. Probability", "s<", "s>" ) )
	print( 56 * "-" )
	vag = qm3.maths.interpolation.hermite_spline( s_coor, v_adia, "akima" )
	erm = qm3.maths.interpolation.hermite_spline( s_coor,   erm_, "akima" )
	for i in range( npts ):
# -------------------------------------- #
#		e_cr = r_zpe + 0.5 * ( 1.0 + qm3.maths.integration.gauss_legendre_xi[i] ) * e_ds
# -------------------------------------- #
		e_cr = e_mx + 0.5 * ( 1.0 + qm3.maths.integration.gauss_legendre_xi[i] ) * e_ds
# -------------------------------------- #
		ener.append( e_cr )
		s_lt, s_gt = __get_turning_points( s_coor, v_adia, s_maxi, e_cr, vag )
		o = c * qm3.maths.integration.Gauss( lambda x: math.sqrt( 2.0 * math.fabs( erm.calc( x )[0] ) * math.fabs( vag.calc( x )[0] - e_cr ) ), s_lt, s_gt )
		prob.append( 1.0 / ( 1.0 + math.exp( 2.0 * o ) ) )
		print( "%16.3lf%20.4le%10.4lf%10.4lf"%( e_cr - r_zpe, prob[-1], s_lt, s_gt ) )
	return( ener, prob )
	


################################################################################################################
def calc_kappa( energy, probability, temperature = 298.15 ):
	"""
	%cappa = e^{{V^max_a over {kT} }} 1 over {kT} int_{0}^{%infinito}{P left( V_a right) e^{-{V_a over {kT} }} dV_a}

		energy			adiabatic energy (V+ZPE) relative to reactants potential (V) in kJ/mol
		probability		in adimensional
		temperature		in kelvin
	"""
	r = 1.0e3 / ( qm3.constants.R * temperature )
# -- don't know why... (subroutine boltz / polyag.f)
#	n = len( qm3.maths.integration.gauss_legendre_xi )
#	d = energy[-1] - energy[0]
#	o = 0.0
#	for i in range( n ):
#		t = math.exp( 0.5 * d * r * ( 1.0 - qm3.maths.integration.gauss_legendre_xi[i] ) )
#		t -= 1.0 / t
#		o += qm3.maths.integration.gauss_legendre_wi[i] * t * probability[i]
#	return( 1.0 + 0.5 * d * r * o )
# --------------------------------------------------------
# The Journal of Physical Chemistry, Vol. 84, No. 13, 1980
#
	ene = qm3.maths.interpolation.hermite_spline( energy, [ i * math.exp( - j * r ) for i,j in zip( probability, energy ) ], "akima" )
	o = qm3.maths.integration.Gauss( lambda x: ene.calc( x )[0], energy[0], energy[-1] )
	return( 1.0 + r * o * math.exp( energy[-1] * r ) )



################################################################################################################
def calc_sigma( s_coor, v_adia, temperature = 298.15 ):
	"""
	%GAMMA = e^{- { lline {%DELTA V^max_a - %DELTA V^{s=0}_a  } rline over {k T} }}

		s_coor			in amu^0.5*A
		v_adia			adiabatic energy (V+ZPE) relative to reactants potential (V) in kJ/mol
		temperature		in kelvin
	"""
	def __find_max( s_, v_ ):
		o = qm3.maths.interpolation.hermite_spline( s_, v_, "akima" )
		d = 0.1
		x = 0.0
		f, g = o.calc( x )
		h = ( o.calc( x + d )[1] - o.calc( x - d )[1] ) / ( 2.0 * d )
		i = 0
		while( i < 1000 and math.fabs( g ) > 1.0e-6 ):
			X = g / h
# hessian should be ALWAYS negative...
#			print( "%5d%20.10lf%20.10lf%20.10lf%20.10lf%20.10lf%20.10lf"%( i, x, f, g, h, X, X / math.fabs( X ) * min( math.fabs( X ), d )  ) )
			x -= X / math.fabs( X ) * min( math.fabs( X ), d )
			f, g = o.calc( x )
			h = ( o.calc( x + d )[1] - o.calc( x - d )[1] ) / ( 2.0 * d )
			i += 1
		return( i < 1000, x, f )
	v_ts = v_adia[s_coor.index( 0.0 )]
	flg, s_mx, v_mx = __find_max( s_coor, v_adia )
	print( "Akima-fitted maximum: %14.6lf%14.6lf"%( s_mx, v_mx ) )
	# -- quadratic interpolation...
	t = v_adia[:]
	t.sort()
	w = v_adia.index( t[-1] )
	v = qm3.maths.matrix.mult( 
			qm3.maths.matrix.inverse( [ 1.0, s_coor[w-1], s_coor[w-1]*s_coor[w-1], 1.0, s_coor[w], s_coor[w]*s_coor[w], 1.0, s_coor[w+1], s_coor[w+1]*s_coor[w+1] ], 3, 3 ), 3, 3, 
			[ v_adia[w-1], v_adia[w], v_adia[w+1] ], 3, 1 )
	s_MX = - v[1] / ( 2.0 * v[2] )
	v_MX = v[0] + v[1] * s_MX + v[2] * s_MX * s_MX
	print( "Quad.Interp. maximum: %14.6lf%14.6lf"%( s_MX, v_MX ) )
	# ----------------------------------------------
	if( flg ):
		return( s_mx, v_mx, math.exp( - math.fabs( v_mx - v_ts ) * 1.0e3 / ( qm3.constants.R * temperature ) ) )
	else:
		return( s_MX, v_MX, math.exp( - math.fabs( v_MX - v_ts ) * 1.0e3 / ( qm3.constants.R * temperature ) ) )





################################################################################################################
if( __name__ == "__main__" ):

	T_ = 298.0

	z_, s_, t_e, t_z, t_t, t_k = parse_MEP( "rate.data", T_ )
	v_ = [ i + j for i,j in zip( t_e, t_z ) ]
	sx_, vx_, S_ = calc_sigma( s_, v_, T_ )
	print()
	print( "Maximum of the adiabatic (s = %.3lf): %.4lf"%( sx_, vx_ ) )
	print( "Sigma (T = %.2lf) = %.3lf"%( T_, S_ ) )

	e_, p_ = transmission_probabilities( s_, v_, sx_, vx_, z_, 1.0 )
	K0_ = calc_kappa( e_, p_, T_ )
	print()
	print( "Kappa_ZCT (T = %.2lf): %.3lf"%( T_, K0_ ) )

	t_m = effective_reduced_masses( s_, t_t, t_k )
	# -------
	i = s_.index( 0.0 )
	j = i - 2
	while( j > 1 and t_m[j-1] > t_m[j] ):
		j -= 1
	t_n = j
	j = i + 2
	while( j < len( s_ ) - 1 and t_m[j+1] > t_m[j] ):
		j += 1
	t_p = j
	m_ = []
	for i in range( t_n ):
		m_.append( t_m[t_n] )
	for i in range( t_n, t_p ):
		m_.append( t_m[i] )
	for i in range( t_p, len( s_ ) ):
		m_.append( t_m[t_p] )
	# -------
	e_, p_ = transmission_probabilities( s_, v_, sx_, vx_, z_, m_ )
	K1_ = calc_kappa( e_, p_, T_ )
	print()
	print( "Kappa_SCT (T = %.2lf): %.3lf"%( T_, K1_ ) )



# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.maths.integration
import	qm3.maths.roots
import	random
import	qm3.maths.matrix



def stats( v ):
	if( type( v ) == float ):
		m = v
		s = 0.0
	elif( len( v ) == 1 ):
		m = v[0]
		s = 0.0
	else:
		n = float( len( v ) )
		m = sum( v ) / n
		s = math.sqrt( sum( [ (i-m)*(i-m) for i in v ] ) / float( n - 1.0 ) )
	return( m, s )



def autocorrelation( v, k = 1 ):
	n = len( v )
	if( k < n ):
		m = sum( v ) / float( n )
		t = sum( [ (v[i]-m)*(v[i]-m) for i in range( n ) ] )
		o = sum( [ (v[i]-m)*(v[i+k]-m) for i in range( n - k ) ] )
		return( o / t )
	else:
		return( 0.0 )



# The value of the sampling ratio that arises from any given data sequence is the factor 
# by which the number of configurations sampled must be increased in order to obtain the
# same precision that would result from randomly distributed data points.
def sampling_ratio( v ):
	n = len( v )
	m = sum( v ) / float( n )
	o = 0.0
	for i in range( 1, n ):
		o += ( v[i] - m ) * ( v[i-1] - m )
	o /= sum( [ (i-m)*(i-m) for i in v ] )
	return( ( 1.0 + o ) / ( 1.0 - o ) )



#	http://home.deib.polimi.it/matteucc/Clustering/tutorial_html/kmeans.html
#	http://datasciencelab.wordpress.com/2014/01/15/improved-seeding-for-clustering-with-k-means/
class k_means( object ):

	def __init__( self, data ):
		self._x = data[:]
		self._k = None
		self._mu = None
		self._omu = None
		self.clusters = None
		self.indexes  = None
 

	def _cluster_points( self ):
		self.clusters = {}
		self.indexes  = {}
		for j in range( len( self._x ) ):
			key = min( [ ( i[0], math.fabs( self._x[j] - self._mu[i[0]] ) ) for i in enumerate( self._mu ) ], key = lambda t: t[1] )[0]
			try:
				self.clusters[key].append( self._x[j] )
				self.indexes[key].append( j )
			except KeyError:
				self.clusters[key] = [ self._x[j] ]
				self.indexes[key] = [ j ]
 

	def _reevaluate_centers( self ):
		self._mu = []
		for k in sorted( self.clusters ):
			self._mu.append( sum( self.clusters[k] ) / float( len( self.clusters[k] ) ) )
 

	def _has_converged( self ):
		return( set( self._mu ) == set( self._omu ) and len( self._mu ) == len( self._omu ) )


	def _choose_next_center( self ):
		d2 = [ min( [ (x-c)*(x-c) for c in self._mu ] ) for x in self._x ]
		s = sum( d2 )
		p = [ i / s for i in d2 ]
		c = [ sum( p[0:i+1] ) for i in range( len( p ) ) ]
		r = random.random()
		i = c.index( [ j for j in c if j >= r ][0] )
		return( self._x[i] )
 

	def _init_centers(self):
# ---------------------------------------------------
#		self._mu = random.sample( self._x, self._k )
# ---------------------------------------------------
		self._mu = random.sample( self._x, 1 )
		while( len( self._mu ) < self._k ):
			self._mu.append( self._choose_next_center() )
# ---------------------------------------------------
 

	def find_centers( self, k ):
		self._k = k
		self._omu = random.sample( self._x, self._k )
		self._init_centers()
		while( not self._has_converged() ):
			self._omu = self._mu[:]
			self._cluster_points()
			self._reevaluate_centers()
		return( self.clusters )



# Numerical Recipes in C, python math module documentation
def gammaln( n ):
	n = float( n )
#	c = [ 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 ]
#	y = n
#	t = n + 5.5
#	t -= ( n + 0.5 ) * math.log( t )
#	o = 1.000000000190015
#	for i in c:
#		y += 1.
#		o += i / y
#	return( - t + math.log( 2.5066282746310005 * o / n ) )
	return( math.lgamma( n ) )



def beta( a, b ):
	return( math.exp( gammaln( a ) + gammaln( b ) - gammaln( a + b ) ) )



def __betacf( x, a, b ):
	x = float( x )
	a = float( a )
	b = float( b )
	MAXIT = 100
	EPS   = 3.0e-7
	FPMIN = 1.0e-30
	qab = a + b
	qap = a + 1.0
	qam = a - 1.0
	c = 1.0
	d = 1.0 - qab * x / qap
	if( math.fabs( d ) < FPMIN ):
		d = FPMIN
	d = 1.0 / d
	h = d
	for m in range( 1, MAXIT + 1 ):
		m2 = 2.0 * float( m )
		aa = float( m ) * ( b - float( m ) ) * x / ( ( qam + m2 ) * ( a + m2 ) )
		d = 1.0 + aa * d
		if( math.fabs( d ) < FPMIN ):
			d = FPMIN
		c = 1.0 + aa / c
		if( math.fabs( c ) < FPMIN ):
			c = FPMIN
		d = 1.0 / d
		h *= d * c
		aa = - ( a + float( m ) ) * ( qab + float( m ) ) * x / ( ( a + m2 ) * ( qap + m2 ) )
		d = 1.0 + aa * d
		if( math.fabs( d ) < FPMIN ):
			d = FPMIN
		c = 1.0 + aa / c
		if( math.fabs( c ) < FPMIN ):
			c = FPMIN
		d = 1.0 / d
		de = d * c
		h *= de
		if( math.fabs( de - 1.0 ) < EPS):
			break
	if( m == MAXIT ):
		return( None )
	else:
		return( h )



def betai( x, a, b ):
	x = float( x )
	a = float( a )
	b = float( b )
	if( x < 0.0 or x > 1.0 ):
		return( None )
	if( x == 0.0 or x == 1.0 ):
		 bt = 0.0;
	else:
		bt = math.exp( gammaln( a + b ) - gammaln( a ) - gammaln( b ) + a * math.log( x ) + b * math.log( 1.0 - x ) )
	if( x < ( a + 1.0 ) / ( a + b + 2.0 ) ):
		return( bt * __betacf( x, a, b ) / a )
	else:
		return( 1.0 - bt * __betacf( 1.0 - x, b, a ) / b )



# http://en.wikipedia.org/wiki/Student's_t-distribution
def student( x, v ):
	x = float( x )
	v = float( v )
	return( math.pow( 1.0 + x * x / v, - 0.5 * ( v + 1.0 ) ) / ( math.sqrt( v ) * beta( 0.5, v * 0.5 ) ) )



# Mathematica:  Series[Hypergeometric2F1[a, b, c, x], {x, 0, 3}], Hypergeometric2F1[0.5, (v + 1)/2, 1.5, -x^2/v]
def student_cdf( x, v ):
	x = float( x )
	v = float( v )
	return( None )



# http://en.wikipedia.org/wiki/F-distribution
def fisher( x, v1, v2 ):
	x  = float( x )
	v1 = float( v1 )
	v2 = float( v2 )
	if( x == 0.0 ):
		return( 0.0 )
	return( math.pow( v1 / v2, 0.5 * v1 ) * math.pow( x, 0.5 * v1 - 1.0 ) * math.pow( 1.0 + v1 / v2 * x, - 0.5 * ( v1 + v2 ) ) / beta( v1 * 0.5, v2 * 0.5 ) )



def fisher_cdf( x, v1, v2 ):
	x  = float( x )
	v1 = float( v1 )
	v2 = float( v2 )
	if( x == 0.0 ):
		return( 0.0 )
	return( betai( v1 * x / ( v1 * x + v2 ), 0.5 * v1, 0.5 * v2 ) )



# http://en.wikipedia.org/wiki/Normal_distribution
def normal( x, m = 0.0, s = 1.0 ):
	x = float( x )
	m = float( m )
	s = float( s )
	return( 1.0 / ( s * math.sqrt( 2.0 * math.pi ) ) * math.exp( -0.5 * math.pow( ( x - m ) / s, 2.0 ) ) )



def normal_cdf( x, m = 0.0, s = 1.0 ):
	x = float( x )
	m = float( m )
	s = float( s )
	return( 0.5 * ( 1.0 + math.erf( ( x - m ) / ( math.sqrt( 2.0 ) * s ) ) ) )



# http://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm
def inverse_student( alfa, v ):
	return( round( qm3.maths.roots.newton_raphson( lambda t: 0.5 * ( 1.0 - 2.0 * alfa ) - 
		qm3.maths.integration.Gauss( lambda x: student( x, v ), 0.0, t ), 0.5, max_iter = 1000, eps = 1.0e-6 ), 3 ) )
#		qm3.maths.integration.Simpson_f( lambda x: student( x, v ), 0.0, t, n = 1000 )[0], 0.5, max_iter = 1000, eps = 1.e-6 ), 3 ) )



# http://www.itl.nist.gov/div898/handbook/eda/section3/eda3673.htm
def inverse_fisher( alfa, v1, v2 ):
	x = 0.0
	while( fisher( x, v1, v2 ) < 0.1 ):
		x += 0.1
	return( round( qm3.maths.roots.newton_raphson( lambda t: ( 1.0 - alfa ) - fisher_cdf( t, v1, v2 ), x, max_iter = 1000, eps = 1.0e-6 ), 3 ) )



# http://www.itl.nist.gov/div898/handbook/eda/section3/eda3671.htm
def inverse_normal( alfa, m = 0.0, s = 1.0 ):
	return( round( qm3.maths.roots.newton_raphson( lambda t: ( 1.0 - 0.5 * alfa ) - normal_cdf( t, m, s ), m, max_iter = 1000, eps = 1.0e-6 ), 3 ) )
		


#
# check for BOTH SUBSETS being RELATED
#
def i_test( na, ma, sa, nb, mb, sb, a = 0.05 ):		# definition based on t-Student, so kinda compatible results
	na = float( na )
	ma = float( ma )
	sa = float( sa )
	nb = float( nb )
	mb = float( mb )
	sb = float( sb )
	if( na >= 30 ):
		ca = inverse_normal( a, ma, sa )
	elif( na > 1 ):
		ca = inverse_student( a, na - 1 )
	else:
		ca = 0.0
	ea = ca * sa / math.sqrt( na )
	if( nb >= 30 ):
		cb = inverse_normal( a, mb, sb )
	elif( nb > 1 ):
		cb = inverse_student( a, nb - 1 )
	else:
		cb = 0.0
	eb = cb * sb / math.sqrt( nb )
	if( ma < mb ):
		o = ( ma + ea >= mb - eb )
	else:
		o = ( mb + eb >= ma - ea )
	print( "(%.6lf,%.6lf,%.6lf) overlaps (%.6lf,%.6lf,%.6lf): %s"%( ma-ea, ma, ma+ea, mb-eb, mb, mb+eb, o ) )
	return( o )



# Check for both samples having the same mean (unpaired data without common values between subsets)
def t_test( na, ma, sa, nb, mb, sb, a = 0.05 ):
	na = float( na )
	ma = float( ma )
	sa = float( sa )
	nb = float( nb )
	mb = float( mb )
	sb = float( sb )
	tc = math.fabs( ma - mb ) / math.sqrt( sa * sa / na + sb * sb / nb )
	df = math.pow( sa * sa / na + sb * sb / nb, 2.0 )  / ( sa * sa * sa * sa / ( na * na * ( na - 1.0 ) ) + sb * sb * sb * sb / ( nb * nb * ( nb - 1.0 ) ) )
	df = max( 1, int( round( df, 0 ) ) )
	tt = inverse_student( a, df )
	o = ( tc <= tt )
	print( "t_calc(%.6lf) <= t_tab(%d,%.4lf = %.6lf): %s"%( tc, df, a, tt, o ) )
	return( o )



# Check for both samples having the same variance (unpaired data without common values between subsets)
def f_test( na, sa, nb, sb, a = 0.05 ):
	na = float( na )
	sa = float( sa )
	nb = float( nb )
	sb = float( sb )
	if( sa > sb ):
		tc = sa * sa / ( sb * sb )
		tt = inverse_fisher( a, na - 1, nb - 1 )
		o = tc <= tt
		print( "f_calc(%.6lf) <= f_tab(%d,%d,%.4lf = %.6lf): %s"%( tc, na - 1, nb - 1, a, tt, o ) )
	else:
		tc = sb * sb / ( sa * sa )
		tt = inverse_fisher( a, nb - 1, na - 1 )
		o = tc <= tt
		print( "f_calc(%.6lf) <= f_tab(%d,%d,%.4lf = %.6lf): %s"%( tc, nb - 1, na - 1, a, tt, o ) )
	return( o )



# - Principal Components Analysis
# X: [ [x1], ..., [xk] ] (k,N)
#
class PCA( object ):

	def __init__( self, x ):
		self.col = len( x )
		self.row = len( x[0] )
		if( sum( [ len( x[i] ) for i in range( self.col ) ] ) != self.col * self.row ):
			raise Exception( "PCA: Invalid dimensions: X_{N,k}" )
		nm = float( self.row )
		self.med = [ sum( x[i] ) / nm for i in range( self.col ) ]
		self.dat = [ [ x[i][j] - self.med[i] for j in range( self.row ) ] for i in range( self.col ) ]
		cov = []
		for i in range( self.col ):
			for j in range( i, self.col ):
				cov.append( sum( [ self.dat[i][l]*self.dat[j][l] for l in range( self.row ) ] ) / nm )
		print( "[Cov]" )
		cov = qm3.maths.matrix.from_upper_diagonal_rows( cov, self.col )
		qm3.maths.matrix.mprint( cov, self.col, self.col )
		self.val, self.vec = qm3.maths.matrix.diag( cov, self.col )
#		self.val, self.vec = qm3.maths._matrix.jacobi( cov, self.col )
		print( "Eig-val: ", self.val )
		print( "Eig-vec:" )
		qm3.maths.matrix.mprint( self.vec, self.col, self.col )


	# Selection indexes: [ 0-k ]
	def select( self, sel ):
		if( len( sel ) > self.col ):
			raise Exception( "PCA: Invalid dimensions: sel_{max(k),1}" )
		self.red = len( sel )
		self.QT = []
		for i in sel:
			for j in range( self.col ):
				self.QT.append( self.vec[j*self.col+i] )
		self.Q = qm3.maths.matrix.T( self.QT, self.red, self.col )
		print( "[Q]" )
		qm3.maths.matrix.mprint( self.Q, self.col, self.red )
		print( "[reduced]" )
		self.out = qm3.maths.matrix.mult( qm3.maths.matrix.T( sum( self.dat, [] ), self.col, self.row ), self.row, self.col, self.Q, self.col, self.red )
		#normalized (-mean)
		#qm3.maths.matrix.mprint( self.out, self.row, self.red )
		for i in range( self.row ):
			for j in range( self.red ):
				self.out[i*self.red+j] += self.med[sel[j]]
		#original axis
		qm3.maths.matrix.mprint( self.out, self.row, self.red )
		

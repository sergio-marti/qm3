# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.maths.rand
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


# REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE
# REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE
# http://home.deib.polimi.it/matteucc/Clustering/tutorial_html/kmeans.html
# http://datasciencelab.wordpress.com/2014/01/15/improved-seeding-for-clustering-with-k-means/
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
		r = qm3.maths.rand.random()
		i = c.index( [ j for j in c if j >= r ][0] )
		return( self._x[i] )
 

	def _init_centers(self):
# ---------------------------------------------------
#		self._mu = qm3.maths.rand.random.sample( self._x, self._k )
# ---------------------------------------------------
		self._mu = qm3.maths.rand.sample( self._x, 1 )
		while( len( self._mu ) < self._k ):
			self._mu.append( self._choose_next_center() )
# ---------------------------------------------------
 

	def find_centers( self, k ):
		self._k = k
		self._omu = qm3.maths.rand.sample( self._x, self._k )
		self._init_centers()
		while( not self._has_converged() ):
			self._omu = self._mu[:]
			self._cluster_points()
			self._reevaluate_centers()
		return( self.clusters )
# REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE
# REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE - REWRITE



# - Principal Components Analysis
# X: [ [x1_N], ..., [xk_N] ] (vars:k, data:N)
#
class PCA( object ):

	def __init__( self, x ):
		self.var = len( x )
		self.dim = len( x[0] )
		if( sum( [ len( x[i] ) for i in range( self.var ) ] ) != self.var * self.dim ):
			raise Exception( "PCA: All variables (rows) have not the same dimensions" )
		self.med = [ sum( x[i] ) / float( self.dim ) for i in range( self.var ) ]
		self.dat = [ [ x[i][j] - self.med[i] for j in range( self.dim ) ] for i in range( self.var ) ]
		cov = []
		for i in range( self.var ):
			for j in range( i, self.var ):
				cov.append( sum( [ self.dat[i][l] * self.dat[j][l] for l in range( self.dim ) ] ) / float( self.dim ) )
		cov = qm3.maths.matrix.from_upper_diagonal_rows( cov, self.var )
		try:
			self.val, self.vec, conv = qm3.maths._matrix.jacobi( cov, self.var )
			self.flg = True
		except:
			# diag SORTS eigenvalues: reduced systems can not be fully recovered from mean data...
			self.val, self.vec = qm3.maths.matrix.diag( cov, self.var )
			self.flg = False


	# Selection indexes: [ 0-k ]
	def select( self, sel, reduced = True ):
		ind = sorted( list( set( [ i for i in sel if i >= 0 and i < self.var ] ) ) )
		red = len( ind )
		if( red == 0 ):
			raise Exception( "PCA: Invalid selection" )
		mat = []
		for i in ind:
			for j in range( self.var ):
				mat.append( self.vec[j*self.var+i] )
		if( reduced ):
			out = qm3.maths.matrix.mult( mat, red, self.var, sum( self.dat, [] ), self.var, self.dim  )
			if( self.flg ):
				for i in range( red ):
					for j in range( self.dim ):
						out[i*self.dim+j] += self.med[ind[i]]
		else:
			mat = qm3.maths.matrix.mult( qm3.maths.matrix.T( mat, red, self.var ), self.var, red, mat, red, self.var )
			out = qm3.maths.matrix.mult( mat, self.var, self.var, sum( self.dat, [] ), self.var, self.dim  )
			for i in range( self.var ):
				for j in range( self.dim ):
					out[i*self.dim+j] += self.med[i]
		return( out )



try:
	import numpy
	def npk_means( data, K ):
		M = [ data[numpy.random.randint(data.shape[0])] ]
		while( len( M ) < K ):
			d2 = numpy.array( [ min( [ numpy.power( numpy.linalg.norm( x - c ), 2.0 ) for c in M ] ) for x in data ] )
			cp = ( d2 / d2.sum() ).cumsum()
			r  = numpy.random.random()
			M.append( data[numpy.where( cp >= r )[0][0]] )
		M = numpy.array( M )
		C = None
		I = None
		o = data[numpy.random.choice( range( data.shape[0] ), K, replace = False )]
		while( C == None or numpy.setdiff1d( numpy.unique( o ), numpy.unique( M ) ).size != 0 ):
			o =M 
			C = {}
			I = {}
			for j in range( data.shape[0] ):
				w = min( [ ( numpy.linalg.norm( data[j] - M[i] ), i ) for i in range( M.shape[0] ) ] )[1]
				try:
					C[w].append( data[j] )
					I[w].append( j )
				except:
					C[w] = [ data[j] ]
					I[w] = [ j ]
			M = numpy.array( [ numpy.mean( C[k], axis = 0 ) for k in iter( C ) ] )
		return( C, I )
	

	class npPCA( object ):
		def __init__( self, data ):
			self.var = data.shape[0]
			self.dim = data.shape[1]
			self.med = data.mean( axis = 1 )
			self.dat = numpy.array( [ data[i,:] - self.med[i] for i in range( self.var ) ] )
			cov = numpy.dot( self.dat, self.dat.T ) / self.dim
			self.val, self.vec = numpy.linalg.eigh( cov )
	
		def select( self, sel, reduced = True ):
			ind = sorted( list( set( [ i for i in sel if i >= 0 and i < self.var ] ) ) )
			if( reduced ):
				red = self.vec[:,ind].T
				out = red.dot( self.dat )
			else:
				red = self.vec[:,ind].dot( self.vec[:,ind].T )
				out = red.dot( self.dat )
				for i in range( self.var ):
					out[i,:] += self.med[i]
			return( out )
except:
	pass

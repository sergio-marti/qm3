# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.actions.minimize
import	qm3.problem
import	qm3.maths.matrix




"""
                                Samples for PES fitting

			[http://www.originlab.com/doc/Origin-Help/Curve-Fitting-Function]


def function( x, c ):
	return( c[0] * ( 1 + c[1] * math.exp( - x / c[2] ) ) / ( 1 + c[3] * math.exp( - x / c[4] ) ) )


c[2] / c[5] can be fixed at minima

def function( x, c ):
	return( c[0] * ( 1 - math.exp( - c[1] * ( x - c[2] ) ) ** 2 + c[3] * math.exp( - c[4] * ( x - c[5] ) ** 2 ) ) )


"""


class problem( qm3.problem.template ):

	def __init__( self, x, y, ext_func, ext_parm ):
		self.efun = ext_func
		self.coor = ext_parm[:]
		self.size = len( self.coor )
		self.func = 0
		self.grad = [ 0.0 for i in range( self.size ) ]
		self.hess = None
		self.__x = x[:]
		self.__y = y[:]
		self.__n = len( self.__x )


	def get_func( self ):
		self.func = 0.0
		for i in range( self.__n ):
			self.func += math.pow( self.efun( self.__x[i], self.coor ) - self.__y[i], 2.0 )


	def get_grad( self ):
		self.get_func()
		be = self.func
		ic = 1.0e-6
		for i in range( self.size ):
			bc = self.coor[i]
			self.coor[i] += ic
			self.get_func()
			self.grad[i] = ( self.func - be ) / ic
			self.coor[i] = bc
		self.func = be
	

	def correlation( self ):
		y_m = sum( self.__y ) / float( self.__n )
		y_d = 0.0
		y_s = 0.0
		for i in range( self.__n ):
			p_v = self.__y[i] - self.efun( self.__x[i], self.coor )
			y_d += p_v * p_v
			y_s += math.pow( self.__y[i] - y_m, 2.0 )
		self.Rsq = 1.0 - y_d / y_s
		return( self.Rsq )

	
	def table( self, fname = None ):
		if( fname ):
			fd = open( fname, "wt" )
		else:
			fd = sys.stdout
		fd.write( "# %.10lf\n"%( self.Rsq ) )
		fd.write( "# " + str.join( " ", [ "%.10lf"%( self.coor[i] ) for i in range( self.size ) ] ) + "\n" )
		fd.write( "#%19s%20s%20s\n"%( "X", "Y", "Fitted" ) )
		for i in range( self.__n ):
			fd.write( "%20.10lf%20.10lf%20.10lf\n"%( self.__x[i], self.__y[i], self.efun( self.__x[i], self.coor ) ) )
		if( fname ):
			fd.close()


	def fit( self, minimize_func = None ):
		if( minimize_func ):
			minimize_func( self )
		else:
			qm3.actions.minimize.steepest_descent( self, step_number = 1000, print_frequency = 1000, gradient_tolerance = 100.0, step_size = 0.1 )
			qm3.actions.minimize.l_bfgs( self, step_number = 1000, print_frequency = 100, gradient_tolerance = 0.01, step_size = 0.1 )
		self.correlation()
		return( self.Rsq, self.coor )





# - Multiple Linear Regression:
# Y: [] (N)
# X: [ [x1], ..., [xk] ] (k,N)
#
# correlation (R^2)
# coefficients: [ao, a1, ..., ak]
# standard deviations for coefficients: [s(ao), s(a1), ..., s(ak)]
def MLR( x, y, normalize = False ):
	n = len( y )
	k = len( x )
	if( sum( [ len( x[i] ) for i in range( k ) ] ) != k * n ):
		raise Exception( "MLR: Invalid dimensions: X_{N,k} vs Y_{N,1}" )
	if( normalize ):
		# --------------------------------------------------------------------------------
		# Normalized X & Y: [a1, ..., ak]
		#
		my = sum( y ) / float( n )
		ny = [ float( y[i] ) - my for i in range( n ) ]
		mx = []
		for i in range( k ):
			mx.append( sum( x[i] ) / float( n ) )
		nx = []
		for j in range( k ):
			t = []
			for i in range( n ):
				t.append( float( x[j][i] ) - mx[j] )
			nx.append( t[:] )
		v = []
		for j in range( k ):
			v.append( sum( [ ny[i]*nx[j][i] for i in range( n ) ] ) )
		m = []
		for i in range( k ):
			for j in range( i, k ):
				m.append( sum( [ nx[i][l]*nx[j][l] for l in range( n ) ] ) )
		b = qm3.maths.matrix.inverse( qm3.maths.matrix.from_upper_diagonal_rows( m, k ), k, k )
		c = qm3.maths.matrix.mult( b, k, k, v, k, 1 )
		r = sum( [ ny[i]*ny[i] for i in range( n ) ] )
		t = 0.0
		for i in range( n ):
			tt = 0.0
			for j in range( k ):
				tt += nx[j][i] * c[j]
			t += ( tt - ny[i] ) * ( tt - ny[i] )
		r = 1 - t / r
		return( r, c, [ b[i*k+i] for i in range( k ) ], mx, my )
	else:
	# --------------------------------------------------------------------------------
	# UN-Normalized X & Y: [a0,a1, ..., ak]
	#
		v = [ float( sum( y ) ) ]
		for i in range( k ):
			v.append( float( sum( [ ii*jj for ii,jj in zip( y, x[i] ) ] ) ) )
		m = [ float( n ) ]
		for i in range( k ):
			m.append( float( sum( x[i] ) ) )
		for i in range( k ):
			for j in range( i, k ):
				m.append( float( sum( [ ii*jj for ii,jj in zip( x[i], x[j] ) ] ) ) )
		b = qm3.maths.matrix.inverse( qm3.maths.matrix.from_upper_diagonal_rows( m, k+1 ), k+1, k+1 )
		c = qm3.maths.matrix.mult( b, k+1, k+1, v, k+1, 1 )
		r = sum( [ (float(i)-v[0]/float(n))*(float(i)-v[0]/float(n)) for i in y ] )
		t = 0.0
		for i in range( n ):
			tt = c[0]
			for j in range( k ):
				tt += x[j][i] * c[j+1]
			t += ( tt - y[i] ) * ( tt - y[i] )
		r = 1 - t / r
		if( n > k + 1 ):
			return( r, c, [ t * b[i*(k+1)+i] / ( n - k - 1.0 ) for i in range( k+1 ) ] )
		else:
			return( r, c, None )



# - Partial Least Squares
# Y: [] (N)
# X: [ [x1], ..., [xk] ] (k,N)
#
# correlation (R^2)
# coefficients: [a1, ..., ak]
# mean-X: [ <x1>, ..., <xk> ]
# mean-Y: <y>
def PLS( x, y ):
	n = len( y )
	k = len( x )
	if( sum( [ len( i ) for i in x ] ) != k * n or len( y ) != n ):
		raise Exception( "PLS: Invalid dimensions: X_{N,k} vs Y_{N,1}" )
	mx = [ float( sum( x[i] ) ) / float( n ) for i in range( k ) ]
	nx = []
	for j in range( n ):
		for i in range( k ):
			nx.append( float( x[i][j] ) - mx[i] )
	ox = nx[:]
	my = float( sum( y ) ) / float( n )
	ny = [ float( y[i] - my ) for i in range( n ) ]
	u = ny[:]
	D = []
	C = []
	P = []
	ex = 1e9
	while( ex > 0.0001 ):
		eu = 1e9
		while( eu > 0.0001 ):
			w = qm3.maths.matrix.norm( qm3.maths.matrix.mult( qm3.maths.matrix.T( nx, n, k ), k, n, u, n, 1 ) )
			t = qm3.maths.matrix.norm( qm3.maths.matrix.mult( nx, n, k, w, k, 1 ) )
			c = sum( [ ny[i]*t[i] for i in range( n ) ] )
			c /= math.fabs( c )
			u = [ ny[i] * c for i in range( n ) ]
			eu = math.sqrt( sum( [ (ny[i]-u[i])*(ny[i]-u[i]) for i in range( n ) ] ) )
		C.append( c )
		b = sum( [ t[i]*u[i] for i in range( n )  ] )
		D.append( b )
#		_p = sum( [ t[i]*t[i] for i in range( n ) ] )
#		p = [ ii/_p for ii in qm3.maths.matrix.mult( qm3.maths.matrix.T( nx, n, k ), k, n, t, n, 1 ) ]
		p = qm3.maths.matrix.mult( qm3.maths.matrix.T( nx, n, k ), k, n, t, n, 1 )
		P += p[:]
#		q = sum( [ ny[i]*u[i] for i in range( n ) ] ) / sum( [ u[i]*u[i] for i in range( n ) ] )
		nx = [ ii-jj for ii,jj in zip( nx, qm3.maths.matrix.mult( t, n, 1, p, 1, k ) ) ]
		ny = [ ny[i]-b*t[i]*c for i in range( n ) ]
		ex = math.sqrt( sum( [ nx[i]*nx[i] for i in range( n * k ) ] ) )
	m = len( C )
	B = qm3.maths.matrix.from_diagonal( D, m )
	coef = qm3.maths.matrix.mult( qm3.maths.matrix.inverse( P, m, k ), k, m, qm3.maths.matrix.mult( B, m, m, C, m, 1 ), m, 1 )
	sr = 0.0
	sy = 0.0
	for i in range( n ):
		sy += ( y[i] - my ) * ( y[i] - my )
		yc = 0.0
		for j in range( k ):
			yc += ox[i*k+j] * coef[j]
		sr += ( yc + my - y[i] ) * ( yc + my - y[i] )
	return( 1 - sr / sy, coef, mx, my )



def poly_val( v, x ):
	return( sum( [ v[i] * math.pow( x, i ) for i in range( len( v ) ) ] ) )



#  [ a0, a1, ..., an ]  ==  a0 + a1 * x + ... + an * x^n
def poly_fit( vec_x, vec_y, order ):
	def __gen_sum( v, k ):
		return( sum( [ math.pow( v[i], k ) for i in range( len( v ) ) ] ) )

	def __gen_yxsum( vy, vx, k ):
		return( sum( [ vy[i]*math.pow( vx[i], k ) for i in range( len( vy ) ) ] ) )

	siz = len( vec_x )
	ror = order + 1
	if( siz < ror or siz != len( vec_y ) ):
		return( None )
	# Independent terms...
	ind = [ __gen_sum( vec_y, 1 ) ]
	for i in range( 1, ror ):
		ind.append( __gen_yxsum( vec_y, vec_x, i ) )
	# Matrix terms...
	t = [ siz ]
	for i in range( 1, 2 * ror - 1 ):
		t.append( __gen_sum( vec_x, i ) )
	mat = []
	for i in range( ror ):
		for j in range( ror ):
			mat.append( t[i+j] )
	sdc = qm3.maths.matrix.inverse( mat, ror, ror )
	cof = qm3.maths.matrix.mult( sdc, ror, ror, ind, ror, 1 )
	# Calculate correlation of the fit
	y_m = __gen_sum( vec_y, 1 ) / float( siz )
	y_d = 0.0
	y_s = 0.0
	for i in range( siz ):
		t = vec_y[i] - poly_val( cof, vec_x[i] )
		y_d += t * t
		y_s += ( vec_y[i] - y_m ) * ( vec_y[i] - y_m )
	return( 1.0 - y_d / y_s, cof, [ sdc[i*ror+i] for i in range( ror ) ] )



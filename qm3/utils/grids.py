# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	qm3.io
import	re
import	math
import	qm3.maths.interpolation
import	os


sys.path.append( os.getenv( "QM3_MPLOT3D" ) )
try:
	import	matplotlib.pyplot
	from	mplot3d import axes3d
	from	mplot3d import proj3d
	import	numpy
	has_mplot3d = True
except:
	has_mplot3d = False


try:
	import	qm3.utils._grids
except:
	pass



def get_ranges( fname ):
	t_x = {}
	t_y = {}
	f = qm3.io.open_r( fname )
	for l in f:
		t = l.split()
		if( len( t ) == 3 ):
			t_x[float( t[0] )] = None
			t_y[float( t[1] )] = None
	qm3.io.close( f, fname )
	t_x = sorted( t_x )
	d_x = 0.0
	for i in range( 1, len( t_x ) ):
		d_x += t_x[i] - t_x[i-1]
	d_x = d_x / float( len( t_x ) - 1 )
	print( t_x[0], t_x[-1], d_x, int( round( ( t_x[-1] - t_x[0] ) / d_x, 0 ) ) + 1 )
	t_y = sorted( t_y )
	d_y = 0.0
	for i in range( 1, len( t_y ) ):
		d_y += t_y[i] - t_y[i-1]
	d_y = d_y / float( len( t_y ) - 1 )
	print( t_y[0], t_y[-1], d_y, int( round(  ( t_y[-1] - t_y[0] ) / d_y, 0 ) ) + 1 )
	return( int( round( ( t_x[-1] - t_x[0] ) / d_x, 0 ) ) + 1, int( round(  ( t_y[-1] - t_y[0] ) / d_y, 0 ) ) + 1 )




class grid( object ):

	__number = re.compile( "^[0-9\.\-eE]+$" )

	def __init__( self, fname = None, interpolant = qm3.maths.interpolation.cubic_spline ):
		self.x = []
		self.y = []
		self.z = []
		self.__intp = interpolant
		self.__spln = None
		if( fname ):
			self.parse( fname )


	def calc( self, rx, ry ):
		return( None )


	# Whatever the ordering ALWAYS returns: fixed_X, changing_Y
	def parse( self, fname = None ):
		f = qm3.io.open_r( fname )
		d = []
		for l in f:
			t = l.strip().split()
			if( len( t ) == 3 and self.__number.match( t[0] ) and self.__number.match( t[1] ) and self.__number.match( t[2] ) ):
				d.append( [ float( t[0] ), float( t[1] ), float( t[2] ) ] )
		qm3.io.close( f, fname )
		d.sort()
		ny = 0
		t  = d[ny][0]
		while( d[ny][0] == t ):
			ny += 1
		nz = len( d )
		nx = nz // ny
		self.x = [ d[i*ny][0] for i in range( nx ) ]
		self.y = [ d[i][1] for i in range( ny ) ]
		self.z = [ d[i][2] for i in range( nz ) ]
		if( self.__intp ):
			self.__spln = qm3.maths.interpolation.interpolate_2d( self.x, self.y, self.z, self.__intp )
			self.calc = self.__spln.calc


	# Transform Z into: changing_X, fixed_Y
	def rotate( self ):
		t = []
		k = 0
		for i in self.x:
			for j in self.y:
				t.append( [ j, i, self.z[k] ] )
				k += 1
		t.sort()
		return( [ t[i][2] for i in range( k ) ] )


	def regular( self, fname = None, points = ( 10, 10 ), gauss = ( 0.1, 0.1 ) ):
		def __pythag( dx, dy ):
			x = math.fabs( dx )
			y = math.fabs( dy )
			if( x > y ):
				return( x * math.sqrt( 1.0 + y * y / ( x * x ) ) )
			if( y == 0.0 ):
				return( 0.0 )
			return( y * math.sqrt( 1.0 + x * x / ( y * y ) ) )
		dat = []
		min_x = None
		min_y = None
		max_x = None
		max_y = None
		f = qm3.io.open_r( fname )
		for l in f:
			t = l.split()
			if( len( t ) == 3 ):
				if( self.__number.match( t[0] ) and self.__number.match( t[1] ) and self.__number.match( t[2] ) ):
					rx = float( t[0] )
					ry = float( t[1] )
					if( min_x and min_y and max_x and max_y ):
						min_x = min( min_x, rx )
						min_y = min( min_y, ry )
						max_x = max( max_x, rx )
						max_y = max( max_y, ry )
					else:
						min_x = rx
						min_y = ry
						max_x = rx
						max_y = ry
					dat.append( [ rx, ry, float( t[2] ) ] )
		qm3.io.close( f, fname )
		dx = ( max_x - min_x ) / float( points[0] - 1.0 )
		print( "[X] delta: %.4lf  points: %3d  range: %8.2lf / %8.2lf"%( dx, points[0], min_x, max_x ) )
		dy = ( max_y - min_y ) / float( points[1] - 1.0 )
		print( "[Y] delta: %.4lf  points: %3d  range: %8.2lf / %8.2lf"%( dy, points[1], min_y, max_y ) )
		self.x = []
		for i in range( points[0] ):
			self.x.append( min_x + dx * i )
		self.y = []
		for i in range( points[1] ):
			self.y.append( min_y + dy * i )
		try:
			self.z = qm3.utils._grids.regular( self.x, self.y, dat, gauss )
		except:
			self.z = []
			for i in self.x:
				for j in self.y:
					rz = 0.0
					rw = 0.0
					for a,b,c in dat:
						dst = __pythag( ( a - i ) / gauss[0], ( b - j ) / gauss[1] )
						w = math.exp( - dst * dst )
						rz += c * w
						rw += w
					self.z.append( rz / rw )
		if( self.__intp ):
			self.__spln = qm3.maths.interpolation.interpolate_2d( self.x, self.y, self.z, self.__intp )
			self.calc = self.__spln.calc


	def save( self, fname = None ):
		k = 0
		f = qm3.io.open_w( fname )
		for i in self.x:
			for j in self.y:
				f.write( "%18.6lf%18.6lf%18.6lf\n"%( i, j, self.z[k] ) )
				k += 1
			f.write( "\n" )
		qm3.io.close( f, fname )


	def plot( self ):
		if( has_mplot3d ):
			def __orthogonal_proj(zfront, zback):
				a = (zfront+zback)/(zfront-zback)
				b = -2*(zfront*zback)/(zfront-zback)
				return numpy.array([[1,0,0,0], [0,1,0,0], [0,0,a,b], [0,0,-0.0001,zback]])
			proj3d.persp_transformation = __orthogonal_proj
			fig = matplotlib.pyplot.figure()
			axs = fig.gca( projection = "3d" )
#			axs = axes3d.Axes3D( matplotlib.pyplot.figure() )
			rz  = self.rotate()
			nx  = len( self.x )
			ny  = len( self.y )
			lx  = []
			ly  = []
			lz  = []
			for i in range( ny ):
				lx.append( self.x[:] )
				ly.append( nx * [ self.y[i] ] )
				lz.append( rz[i*nx:(i+1)*nx][:] )
			z_min = min( self.z )
			z_max = max( self.z )
			z_lvl = [ z_min + ( z_max - z_min ) / 30. * float( i ) for i in range( 31 ) ]
			lx = numpy.array( lx, dtype=float )
			ly = numpy.array( ly, dtype=float )
			lz = numpy.array( lz, dtype=float )
			axs.plot_surface( lx, ly, lz, rstride = 1, cstride = 1, cmap = "coolwarm", linewidths = 0.1 )
#			axs.contour( lx, ly, lz, zdir = "z", levels = z_lvl, linewidths = 2, cmap = "coolwarm" )
			axs.contour( lx, ly, lz, zdir = "z", offset = z_min, levels = z_lvl, linewidths = 2, cmap = "coolwarm" )
			axs.view_init( 90, -89 )
			matplotlib.pyplot.show()
		else:
			return


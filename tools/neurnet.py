#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

import	math
import	qm3.maths.rand
import	qm3.maths.matrix




class Perceptron( object ):
	def __init__( self, x, y ):
		self.k = len( x[0] )
		self.n = len( y )
		# -- X_n,k
		self.x = sum( x, [] )
		# -- Y_n,1
		self.y = y[:]
		# -- W_k,1 && B_1,1
		self.w = [ qm3.maths.rand.random() for i in range( self.k ) ]
		self.b = qm3.maths.rand.random()


	@staticmethod
	def fun( x ):
		return( 1.0 / ( 1.0 + math.exp( - x ) ) )


	@staticmethod
	def dfun( s ):
		return( s * ( 1.0 - s ) )


	def train( self, steps = 100000, learning_ratio = 0.01, tolerance = 1.e-8 ):
		l = 1e99
		for istep in range( steps ):
			r = [ self.fun( i + self.b ) for i in qm3.maths.matrix.mult( self.x, self.n, self.k, self.w, self.k, 1 ) ]
			e = [ r[i] - self.y[i] for i in range( self.n ) ]
			t = [ e[i] * self.dfun( r[i] ) for i in range( self.n ) ]
			d = qm3.maths.matrix.mult( t, 1, self.n, self.x, self.n, self.k )
			self.w = [ self.w[i] - learning_ratio * d[i] for i in range( self.k ) ]
			self.b -= learning_ratio * sum( t )
			o = l
			l = sum( [ i * i for i in e ] ) / float( self.n )
			if( math.fabs( l - o ) < tolerance ):
				break
		return( l )


	def calc( self, x ):
		return( self.fun( sum( [ x[i] * self.w[i] for i in range( self.k ) ] ) + self.b ) )




class Perceptron_1L( object ):
	def __init__( self, x, y, nodes, output_dim = None ):
		self.n = len( y )
		self.k = len( x[0] )
		if( output_dim ):
			self.o = output_dim
		else:
			self.o = len( set( y ) )
		self.h = max( nodes, self.o + 1 )
		# -- X_n,k
		self.x = sum( x, [] )
		# -- Y_n,o
		self.y = [ 0 for i in range( self.n * self.o ) ]
		for i in range( self.n ):
			self.y[i*self.o+y[i]] = 1
		# -- WH_k,h  / BH_h,1 / WO_h,o / BO_o,1
		self.wh = [ qm3.maths.rand.random() for i in range( self.k * self.h ) ]
		self.bh = [ qm3.maths.rand.random() for i in range( self.h ) ]
		self.wo = [ qm3.maths.rand.random() for i in range( self.h * self.o ) ]
		self.bo = [ qm3.maths.rand.random() for i in range( self.o ) ]


	@staticmethod
	def fun( x ):
		return( 1 / ( 1 + math.exp( - x ) ) )


	@staticmethod
	def dfun( x ):
		def f( x ):
			return( 1 / ( 1 + math.exp( - x ) ) )
		return( f( x ) * ( 1 - f( x ) ) )


	def train( self, steps = 100000, learning_ratio = 0.01, tolerance = 1.e-8 ):
		l = 1e99
		for istep in range( steps ):
			# xh_n,h
			xh = [ i + j for i,j in zip( qm3.maths.matrix.mult( self.x, self.n, self.k, self.wh, self.k, self.h ), self.bh * self.n ) ]
			# fh_n,h
			fh = [ self.fun( i ) for i in xh ]
			# xo_n,o
			xo = [ i + j for i,j in zip( qm3.maths.matrix.mult( fh, self.n, self.h, self.wo, self.h, self.o ), self.bo * self.n ) ]
			# fo_n,o
			fo = [ self.fun( i ) for i in xo ]
			# ee_n,o
			ee = [ fo[i] - self.y[i] for i in range( self.n * self.o ) ]

			# dh_n,h
			dh = [ self.dfun( i ) for i in xh ]
			# do_n,o
			do = [ self.dfun( i ) for i in xo ]

			# tt_n,o
			tt  = [ i * j for i,j in zip( ee, do ) ]

			self.bo = [ i - learning_ratio * j for i,j in zip( self.bo, qm3.maths.matrix.mult( [ 1.0 ] * self.n, 1, self.n, tt, self.n, self.o ) ) ]
			self.wo = [ i - learning_ratio * j for i,j in zip( self.wo, qm3.maths.matrix.mult( qm3.maths.matrix.T( fh, self.n, self.h ), self.h, self.n, tt, self.n, self.o ) ) ]
			# tt_n,h
			tt  = [ i * j for i,j in zip ( qm3.maths.matrix.mult( tt, self.n, self.o, qm3.maths.matrix.T( self.wo, self.h, self.o ), self.o, self.h ), dh ) ]

			self.bh = [ i - learning_ratio * j for i,j in zip( self.bh, qm3.maths.matrix.mult( [ 1.0 ] * self.n, 1, self.n, tt, self.n, self.h ) ) ]
			self.wh = [ i - learning_ratio * j for i,j in zip( self.wh, qm3.maths.matrix.mult( qm3.maths.matrix.T( self.x, self.n, self.k ), self.k, self.n, tt, self.n, self.h ) ) ]

			o = l
			l = sum( [ i * i for i in ee ] ) / float( self.n )
			if( math.fabs( l - o ) < tolerance ):
				break
		return( l )


	def calc( self, x ):
		fh = [ self.fun( i + j ) for i,j in zip( qm3.maths.matrix.mult( x, 1, self.k, self.wh, self.k, self.h ), self.bh ) ]
		return( [ self.fun( i + j ) for i,j in zip( qm3.maths.matrix.mult( fh, 1, self.h, self.wo, self.h, self.o ), self.bo ) ] )




if( __name__ == "__main__" ):
	import	matplotlib.pyplot as plt
	# -----------------------------------------------------------------------------------------
	# binary classification:
	if( False ):
		x = []
		y = []
		f = 0.40
		for i in range( 50 ):
			x.append( [ qm3.maths.rand.random() - f, qm3.maths.rand.random() - f ] )
			y.append( 0 )
		for i in range( 50 ):
			x.append( [ qm3.maths.rand.random() + f, qm3.maths.rand.random() + f ] )
			y.append( 1 )
		
		_x = [ x[i][0] for i in range( 100 ) ]
		_y = [ x[i][1] for i in range( 100 ) ]
		plt.clf()
		plt.grid( True )
		plt.xlim( -0.5, 1.5 )
		plt.ylim( -0.5, 1.5 )
		plt.scatter( _x, _y, c = y, marker = "o", cmap = plt.cm.winter, edgecolors='black' )
		plt.savefig( "2.pdf" )
		
		NN = Perceptron( x, y )
		l = NN.train()
		print( "loss:", l )
		
		_X = []
		_Y = []
		_Z = []
		mx = min( _x )
		dx = ( max( _x ) - mx ) / 500
		my = min( _y )
		dy = ( max( _y ) - my ) / 500
		for i in range( 500 ):
			for j in range( 500 ):
				_X.append( mx + dx * i )
				_Y.append( my + dy * j )
				_Z.append( NN.calc( [ _X[-1], _Y[-1] ] ) )
		plt.clf()
		plt.grid( True )
		plt.xlim( -0.5, 1.5 )
		plt.ylim( -0.5, 1.5 )
		plt.scatter( _X, _Y, c = _Z, marker = ",", cmap = plt.cm.winter, edgecolors = 'face' )
		plt.scatter( _x, _y, c = y, marker = "o", cmap = plt.cm.winter, edgecolors='black' )
		plt.savefig( "2_%.6lf.pdf"%( l ) )
	
	# -----------------------------------------------------------------------------------------
	# ternary classification:
	if( False ):
		x = []
		y = []
		for i in range( 50 ):
			x.append( [ qm3.maths.rand.random() - 0.5, qm3.maths.rand.random() - 1.0 ] )
			y.append( 0 )
		for i in range( 50 ):
			x.append( [ qm3.maths.rand.random() - 1.0, qm3.maths.rand.random() ] )
			y.append( 1 )
		for i in range( 50 ):
			x.append( [ qm3.maths.rand.random(), qm3.maths.rand.random() ] )
			y.append( 2 )
		
		_x = [ x[i][0] for i in range( 150 ) ]
		_y = [ x[i][1] for i in range( 150 ) ]
		plt.clf()
		plt.grid( True )
		plt.xlim( -1.1, 1.1 )
		plt.ylim( -1.1, 1.1 )
		plt.scatter( _x, _y, c = y, marker = "o", cmap = plt.cm.hot, edgecolors='black' )
		plt.savefig( "3.pdf" )
		
		NN = Perceptron_1L( x, y, 4 )
		l = NN.train( learning_ratio = 0.1 )
		print( "loss:", l )
	
		_X = []
		_Y = []
		_Z = []
		mx = min( _x )
		dx = ( max( _x ) - mx ) / 500
		my = min( _y )
		dy = ( max( _y ) - my ) / 500
		for i in range( 500 ):
			for j in range( 500 ):
				_X.append( mx + dx * i )
				_Y.append( my + dy * j )
				r = NN.calc( [ _X[-1], _Y[-1] ] )
				t = max( r )
				_Z.append( r.index( t ) - 1 + t )
		plt.clf()
		plt.grid( True )
		plt.xlim( -1.1, 1.1 )
		plt.ylim( -1.1, 1.1 )
		plt.scatter( _X, _Y, c = _Z, marker = ",", cmap = plt.cm.hot, edgecolors = 'face' )
		plt.scatter( _x, _y, c = y, marker = "o", cmap = plt.cm.hot, edgecolors='black' )
		plt.savefig( "3_%.6lf.pdf"%( l ) )



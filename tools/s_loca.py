#!/usr/bin/env python

import	sys
import	qm3.actions.minimize
import	qm3.maths.interpolation
import	qm3.utils.gridss

class pmf2d:
	def __init__( self, data ):
		self.size = 2
		self.func = .0
		self.coor = [ 0, 0 ]
		self.grad = [ 0, 0 ]
		self.hess = [ 0, 0, 0 ]
		x, y, z = qm3.utils.gridss.parse( data )
		self.spln = qm3.maths.interpolation.interpolate_2d( x, y, z, qm3.maths.interpolation.cubic_spline )
#		self.spln = qm3.maths.interpolation.interpolate_2d( x, y, z, lambda x,y: qm3.maths.interpolation.hermite_spline( x, y, "akima" ) )

	def current_step( self, istep ):
		pass

	def check_coor( self ):
		self.coor[0] = min( max( self.coor[0], self.spln.x[0] ), self.spln.x[-1] )
		self.coor[1] = min( max( self.coor[1], self.spln.y[0] ), self.spln.y[-1] )

	def get_func( self ):
		self.check_coor()
		o = self.spln.calc( self.coor[0], self.coor[1] )
		self.func = o[0]

	def get_grad( self ):
		self.check_coor()
		o = self.spln.calc( self.coor[0], self.coor[1] )
		self.func = o[0]
		self.grad = [ o[1], o[2] ]

	def get_hess( self ):
		self.check_coor()
		o = self.spln.calc( self.coor[0], self.coor[1] )
		self.func = o[0]
		self.grad = [ o[1], o[2] ]
		dh = 1.e-6
		of = self.spln.calc( self.coor[0] + dh, self.coor[1] )
		ob = self.spln.calc( self.coor[0] - dh, self.coor[1] )
		r1 = [ ( of[1] - ob[1] ) / ( 2. * dh ), ( of[2] - ob[2] ) / ( 2. * dh ) ]
		of = self.spln.calc( self.coor[0], self.coor[1] + dh )
		ob = self.spln.calc( self.coor[0], self.coor[1] - dh )
		r2 = [ ( of[1] - ob[1] ) / ( 2. * dh ), ( of[2] - ob[2] ) / ( 2. * dh ) ]
		self.hess = [ r1[0], ( r1[1] + r2[0] ) / 2., ( r1[1] + r2[0] ) / 2., r2[1] ]



x = pmf2d( sys.argv[1] )

x.coor = [ float( sys.argv[2] ), float( sys.argv[3] ) ]
x.get_func()
print x.coor, x.func

#qm3.actions.minimize.baker( x, 
#	step_number = 1000, 
#	step_size = 0.01,
#	follow_mode = int( sys.argv[4] ),
#	print_frequency = 100, 
#	allow_overlap = True,
#	gradient_tolerance = 0.01 )

qm3.actions.minimize.rfo( x, 
	step_number = 1000, 
	step_size = 0.01,
	follow_mode = int( sys.argv[4] ),
	print_frequency = 100, 
	gradient_tolerance = 0.01 )

x.check_coor()
print x.coor, x.func

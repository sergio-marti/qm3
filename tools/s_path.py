#!/usr/bin/env python3

import	sys
import	qm3.actions.paths
import	qm3.maths.interpolation
import	qm3.maths.grids

class pmf2d:
	def __init__( self, data ):
		self.size = 2
		self.mass = [ 1., 1. ]
		self.func = .0
		self.coor = [ 0, 0 ]
		self.grad = [ 0, 0 ]
		self.hess = [ 0, 0, 0 ]
		self.grid = qm3.maths.grids.grid( data, qm3.maths.interpolation.cubic_spline )
#		self.grid = qm3.maths.grids.grid( data, lambda x,y: qm3.maths.interpolation.hermite_spline( x, y, "akima" ) )
		self.fd = None

	def current_step( self, istep ):
		self.check_coor()
		self.fd.write( "%20.10lf%20.10lf%20.10lf\n"%( self.coor[0], self.coor[1], self.func ) )

	def check_coor( self ):
		self.coor[0] = min( max( self.coor[0], self.grid.x[0] ), self.grid.x[-1] )
		self.coor[1] = min( max( self.coor[1], self.grid.y[0] ), self.grid.y[-1] )

	def get_func( self ):
		self.check_coor()
		o = self.grid.calc( self.coor[0], self.coor[1] )
		self.func = o[0]

	def get_grad( self ):
		self.check_coor()
		o = self.grid.calc( self.coor[0], self.coor[1] )
		self.func = o[0]
		self.grad = [ o[1], o[2] ]

	def get_hess( self ):
		self.check_coor()
		o = self.grid.calc( self.coor[0], self.coor[1] )
		self.func = o[0]
		self.grad = [ o[1], o[2] ]
		dh = 1.e-6
		of = self.grid.calc( self.coor[0] + dh, self.coor[1] )
		ob = self.grid.calc( self.coor[0] - dh, self.coor[1] )
		r1 = [ ( of[1] - ob[1] ) / ( 2. * dh ), ( of[2] - ob[2] ) / ( 2. * dh ) ]
		of = self.grid.calc( self.coor[0], self.coor[1] + dh )
		ob = self.grid.calc( self.coor[0], self.coor[1] - dh )
		r2 = [ ( of[1] - ob[1] ) / ( 2. * dh ), ( of[2] - ob[2] ) / ( 2. * dh ) ]
		self.hess = [ r1[0], ( r1[1] + r2[0] ) / 2., ( r1[1] + r2[0] ) / 2., r2[1] ]



x = pmf2d( sys.argv[1] )

x.fd = file( "path.log", "wt" )
x.fd.write( "# Forward\n" )
x.coor = [ float( sys.argv[2] ), float( sys.argv[3] ) ]
x.get_func()
x.current_step( 0 )
qm3.actions.paths.page_mciver( x, step_number = 1000, step_size = 0.01, print_frequency = 1, project_RT = False, avoid_recrossing = False )
#qm3.actions.paths.baker( x, step_number = 1000, step_size = 0.01, print_frequency = 1, project_RT = False, avoid_recrossing = False )
x.fd.write( "# Reverse\n" )
x.coor = [ float( sys.argv[2] ), float( sys.argv[3] ) ]
x.get_func()
x.current_step( 0 )
qm3.actions.paths.page_mciver( x, step_number = 1000, step_size = -0.01, print_frequency = 1, project_RT = False, avoid_recrossing = False )
#qm3.actions.paths.baker( x, step_number = 1000, step_size = -0.01, print_frequency = 1, project_RT = False, avoid_recrossing = False )
x.fd.close()

#!/usr/bin/env python

import	pickle
import	qm3.mol
import	qm3.io.dcd
import	qm3.problem
import	qm3.engines.mol_mech
import	qm3.actions.minimize
import	qm3.actions.dynamics


class my_problem( qm3.problem.template ):
	def __init__( self ):
		qm3.problem.template.__init__( self )
		self.mol = qm3.mol.molecule()
		self.mol.xyz_read( "last.xyz" )
		self.mol.guess_atomic_numbers()
		self.mol.fill_masses()
		self.mol.type = [ self.mol.labl[i] + str( i+1 ) for i in range( self.mol.natm ) ]
		f = open( "last.chg", "rb" )
		self.mol.chrg = pickle.load( f ) 
		f.close()
		self.mol.epsi = [0.5411838874172069, 0.5411838874172069, 0.5254940532489402, 0.5254940532489402, 0.4709055106919009,
			0.8678248671246982, 0.31688483712541377, 0.6033307550589477, 0.8283574107835338, 0.31688483712541377,
			0.31688483712541377, 0.31688483712541377, 0.0, 0.9373579892442374, 0.9373579892442374, 0.0, 0.0, 0.0,
			0.843374175559105, 0.31688483712541377, 0.0]
		self.mol.rmin = [1.992, 1.992, 1.956, 1.956, 1.675, 1.839, 1.381, 1.824, 1.741, 1.381, 1.381, 1.381, 0.0, 1.661,
			1.661, 0.0, 0.0, 0.0, 1.824, 1.381, 0.0]

		self.eng = qm3.engines.mol_mech.simple_force_field()
		self.eng.cut_on   = -1
		self.eng.cut_off  = -1
		self.eng.cut_list = -1
		self.eng.guess_bonds( self.mol )
		self.eng.calc_connectivity()
		self.eng.guess_angles()
		self.eng.guess_dihedrals()
		self.eng.impr = [ [ 7, 13, 14, 3, 52.5, 1.9 ] ]
		self.eng.load_parameters( self.mol, "parm" )

		self.size = 3 * self.mol.natm
		self.coor = self.mol.coor
		self.mass = self.mol.mass
		self.func = 0.0
		self.grad = []

		self.flg = False
		self.dcd = qm3.io.dcd.dcd()
		self.dcd.write( "test.dcd", self.mol.natm )


	def current_step( self, istep ):
		if( self.flg and istep%100 == 0 ):
			self.mol.dcd_write( self.dcd )


	def get_func( self ):
		self.mol.func = 0.0
		self.eng.get_func( self.mol )
		self.func = self.mol.func


	def get_grad( self ):
		self.mol.func = 0.0
		self.mol.grad = [ 0.0 for i in range( self.size )  ]
		self.eng.get_grad( self.mol )
		self.func = self.mol.func
		self.grad = self.mol.grad





o = my_problem()

qm3.actions.minimize.fire( o, step_number = 1000, gradient_tolerance = 0.1 )
o.mol.xyz_write( "test.xyz" )
import qm3.maths.matrix
print( o.func )
qm3.maths.matrix.mprint( o.grad, o.mol.natm, 3 )

o.flg = True
qm3.actions.dynamics.assign_velocities( o, 300.0 )
qm3.actions.dynamics.langevin_verlet( o, temperature = 300.0, step_number = 100000 )
o.dcd.close()

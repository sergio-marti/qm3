#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.mol
import	qm3.elements
import	qm3.problem
import	qm3.constants
import	qm3.actions.minimize
import	qm3.utils
import	qm3.engines.dftb
try:
	import	cPickle as pickle
except:
	import	pickle



class qwerty( qm3.problem.template ):
	def __init__( self ):
		qm3.problem.template.__init__( self )
		self.mol = qm3.mol.molecule()
		self.mol.xyz_read( "xyz" )
		self.mol.chrg = [ 0.0 for i in range( self.mol.natm ) ]
		self.mol.anum = [ qm3.elements.rsymbol[self.mol.labl[i].title()] for i in range( self.mol.natm ) ]

		self.eng = qm3.engines.dftb.dftb( self.mol, list( range( self.mol.natm ) ) )
		self.eng.prm = "/storage/projects/vlc78/xexo/prueba/3ob-3-1/"
		self.eng.exe = "/storage/projects/vlc78/xexo/prueba/dftb+_git > dftb_in.log"
		self.eng.chg = 0

		self.size = 3 * self.mol.natm
		self.coor = self.mol.coor
		self.func = 0.0
		self.grad = []
		self.hess = []


	def get_func( self ):
		self.mol.func = 0.0
		self.eng.get_func( self.mol )
		self.func = self.mol.func


	def get_grad( self ):
		self.mol.func = 0.0
		self.mol.grad = [ 0.0 for i in range( self.size ) ]
		self.eng.get_grad( self.mol )
		self.func = self.mol.func
		self.grad = self.mol.grad[:]


	def get_hess( self ):
		self.num_hess()
		


obj = qwerty()
qm3.actions.minimize.l_bfgs( obj, step_size = 0.1, gradient_tolerance = 0.1, print_frequency = 100, step_number = 1000 )
obj.mol.xyz_write( "last.xyz" )

f = open( "last.chg", "wb" )
pickle.dump( obj.mol.chrg, f )
f.close()

obj.get_hess()
obj.mol.fill_masses()

f = open( "last.hes", "wb" )
pickle.dump( obj.hess, f )
f.close()

val, vec = qm3.utils.hessian_frequencies( obj.mol.mass, obj.mol.coor, obj.hess )
print( val )

#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

import math
import qm3.utils
import qm3.problem



#
# J. Chem. Phys. v113, p9978 (2000) [10.1063/1.1323224]
#

class serial_neb( qm3.problem.template ):
	def __init__( self, problem, sele, kumb, nodes, crd_ini, crd_end ):
		"""
	problem should be a qm3.problem.template working on the 'sele' atoms (see core/envi optimization)
	and should provide the following method:

			def neb_data( self, node )

	ideally it should store current coordinates of the node:
			
				f = open( "node.%02d"%( node ) )
				f.write( "REMARK %12.6lf%12.6lf%20.3lf\n"%( distance_1, distance_2, self.func ) )
				self.mol.pdb_write( f )
				f.close()

	as well as information about the current node as a REMARK within the PDB  (geometry, energy, ...):

	the value of the 'kumb' should be approx the same of the potential energy barrier
	when optimizing the whole band, set the 'gradient_tolerance' equal to 0.5 * nodes (_kJ/mol.A)
		"""
		qm3.problem.template.__init__( self )

		self.node = nodes
		self.kumb = kumb
		self.sele = sele[:]
		self.crd0 = crd_ini[:]
		self.crdf = crd_end[:]
		self.prob = problem

		self.dime = len( self.crd0 )
		self.size = self.dime * self.node
		delt      = [ ( jj - ii ) / float( self.node + 1 ) for ii,jj in zip( self.crd0, self.crdf ) ]
		self.coor = []
		for k in range( 1, self.node + 1 ):
			self.coor += [ ii + jj * k for ii,jj in zip( self.crd0, delt ) ]
		self.mass = [ self.prob.mol.mass[i] for i in self.sele ]
		self.func = 0
		self.grad = []

		# end-points potential
		self.prob.coor = self.crd0[:]
		self.prob.get_func()
		self.pot0 = self.prob.func
		self.prob.coor = self.crdf[:]
		self.prob.get_func()
		self.potf = self.prob.func
		print( self.pot0, self.potf )


	def get_grad_OWN( self ):
		self.func = 0.0
		self.grad = []
		for who in range( self.node ):
			if( who == 0 ):
				# first node
				ref = [ ( ii + jj ) / 2.0 for ii,jj in zip( self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime], self.crd0 ) ]
				gum = [ self.kumb * ( ii - jj ) for ii,jj in zip( self.coor[self.dime*who:self.dime*who+self.dime], ref ) ]
				tau = [ ii-jj for ii,jj in zip( self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime], self.crd0 ) ]
			elif( who == self.node - 1 ):
				# last node
				ref = [ ( ii + jj ) / 2.0 for ii,jj in zip( self.crdf, self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime] ) ]
				gum = [ self.kumb * ( ii - jj ) for ii,jj in zip( self.coor[self.dime*who:self.dime*who+self.dime], ref ) ]
				tau = [ ii-jj for ii,jj in zip( self.crdf, self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime] ) ]
			else:
				# intermediates
				ref = [ ( ii + jj ) / 2.0 for ii,jj in zip( self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime],
						self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime] ) ]
				gum = [ self.kumb * ( ii - jj ) for ii,jj in zip( self.coor[self.dime*who:self.dime*who+self.dime], ref ) ]
				if( who >= 2 and who < self.node - 2 ):
					tau = []
					for j in range( self.dime ):
						# interpolated: a0 + a1 x + a2 x^2 + a3 x^3 + a4 x^4
#						tau.append( ( -8.0 * self.coor[self.dime*(who-2)+j] + self.coor[self.dime*(who-1)+j]
#									+ 8.0 * self.coor[self.dime*(who+1)+j] - self.coor[self.dime*(who+2)+j] ) )
						# LS-fitted: a0 + a1 x + a2 x^2
						tau.append( ( -2.0 * self.coor[self.dime*(who-2)+j] - self.coor[self.dime*(who-1)+j]
									+ self.coor[self.dime*(who+1)+j] + 2.0 *  self.coor[self.dime*(who+2)+j] ) )
				else:
					tau = [ ii-jj for ii,jj in zip( self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime],
							self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime] ) ]
			# common to all nodes
			tmp = math.sqrt( sum( [ ii * ii for ii in tau ] ) )
			tau = [ ii / tmp for ii in tau ]
			self.prob.coor = self.coor[self.dime*who:self.dime*who+self.dime][:]
			self.prob.get_grad()
			self.func += self.prob.func
			try:
				self.prob.neb_data( who )
			except:
				pass
			usp = sum( [ ii * jj for ii,jj in zip( tau, gum ) ] )
			gsp = sum( [ ii * jj for ii,jj in zip( tau, self.prob.grad ) ] )
			grd = [ ii - gsp * jj + usp * jj for ii,jj in zip( self.prob.grad, tau ) ]
			qm3.utils.project_RT_modes( self.mass, self.prob.coor, grd, None )
			self.grad += grd[:]


	def get_grad( self ):
		# ----------------------------------------------------------------------
		def __calc_tau( potm, poti, potp, crdm, crdi, crdp ):
			dcM = [ ii-jj for ii,jj in zip( crdp, crdi ) ]
			dcm = [ ii-jj for ii,jj in zip( crdi, crdm ) ]
			dpM = max( math.fabs( potp - poti ), math.fabs( potm - poti ) )
			dpm = min( math.fabs( potp - poti ), math.fabs( potm - poti ) )
			if( potp > poti and poti > potm ):
				tau = dcM[:]
			elif( potp < poti and poti < potm ):
				tau = dcm[:]
			else:
				if( potp > potm ):
					tau = [ dpM * ii + dpm * jj for ii,jj in zip( dcM, dcm ) ]
				else:
					tau = [ dpm * ii + dpM * jj for ii,jj in zip( dcM, dcm ) ]
			tmp = math.sqrt( sum( [ ii * ii for ii in tau ] ) )
			tau = [ ii / tmp for ii in tau ]
			mcM = math.sqrt( sum( [ ii * ii for ii in dcM ] ) )
			mcm = math.sqrt( sum( [ ii * ii for ii in dcm ] ) )
			gum = [ - self.kumb * ( mcM - mcm ) * ii for ii in tau ]
			return( tau, gum )
		# ----------------------------------------------------------------------
		vpot = []
		gpot = []
		# get potential energy and gradients from the chain elements
		for who in range( self.node ):
			self.prob.coor = self.coor[self.dime*who:self.dime*who+self.dime][:]
			self.prob.get_grad()
			vpot.append( self.prob.func )
			gpot += self.prob.grad[:]
			try:
				self.prob.neb_data( who )
			except:
				pass
		self.func = sum( vpot )
		self.grad = []
		# calculate neb components	
		for who in range( self.node ):
			if( who == 0 ):
				# first node
				tau, gum = __calc_tau( self.pot0, vpot[who], vpot[who+1],
							self.crd0,
							self.coor[self.dime*who:self.dime*who+self.dime],
							self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime] )
			elif( who == self.node - 1 ):
				# last node
				tau, gum = __calc_tau( vpot[who-1], vpot[who], self.potf,
							self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime],
							self.coor[self.dime*who:self.dime*who+self.dime],
							self.crdf )
			else:
				# intermediates
				tau, gum = __calc_tau( vpot[who-1], vpot[who], vpot[who+1],
							self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime],
							self.coor[self.dime*who:self.dime*who+self.dime],
							self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime] )
			# common to all nodes
			gsp = sum( [ ii * jj for ii,jj in zip( tau, gpot[self.dime*who:self.dime*who+self.dime] ) ] )
			grd = [ ii - gsp * jj + kk for ii,jj,kk in zip( gpot[self.dime*who:self.dime*who+self.dime], tau, gum ) ]
			qm3.utils.project_RT_modes( self.mass, self.coor[self.dime*who:self.dime*who+self.dime], grd, None )
			self.grad += grd[:]




try:
	import	qm3.utils._mpi

	class parall_neb( qm3.problem.template ):
		def __init__( self, problem, sele, kumb, node, ncpu, chnk, crd_ini, crd_end ):
			"""
	problem should be a qm3.problem.template working on the 'sele' atoms (see core/envi optimization)
	and should provide the following method:

			def neb_data( self, node )

	ideally it should store current coordinates of the node:
			
				f = open( "node.%02d"%( node ) )
				f.write( "REMARK %12.6lf%12.6lf%20.3lf\n"%( distance_1, distance_2, self.func ) )
				self.mol.pdb_write( f )
				f.close()

	as well as information about the current node as a REMARK within the PDB  (geometry, energy, ...):

	the value of the 'kumb' should be approx the same of the potential energy barrier
	when optimizing the whole band, set the 'gradient_tolerance' equal to 0.5 * nodes (_kJ/mol.A)
			"""
			qm3.problem.template.__init__( self )
	
			self.node = ncpu * chnk
			self.wami = node
			self.ncpu = ncpu
			self.chnk = chnk

			self.kumb = kumb
			self.sele = sele[:]
			self.crd0 = crd_ini[:]
			self.crdf = crd_end[:]
			self.prob = problem
	
			self.dime = len( self.crd0 )
			self.size = self.dime * self.node
			delt      = [ ( jj - ii ) / float( self.node + 1 ) for ii,jj in zip( self.crd0, self.crdf ) ]
			self.coor = []
			for k in range( 1, self.node + 1 ):
				self.coor += [ ii + jj * k for ii,jj in zip( self.crd0, delt ) ]
			self.mass = [ self.prob.mol.mass[i] for i in self.sele ]
			self.func = 0
			self.grad = []

			# end-points potential
			self.prob.coor = self.crd0[:]
			self.prob.get_func()
			self.pot0 = self.prob.func
			self.prob.coor = self.crdf[:]
			self.prob.get_func()
			self.potf = self.prob.func
			if( self.wami == 0 ):
				print( self.pot0, self.potf )
	
	
	
		def get_grad_OWN( self ):
			# calculate assigned nodes
			self.func = 0
			self.grad = []
			for who in range( self.wami * self.chnk, ( self.wami + 1 ) * self.chnk ):
				if( who == 0 ):
					# first node
					ref = [ ( ii + jj ) / 2.0 for ii,jj in zip( self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime], self.crd0 ) ]
					gum = [ self.kumb * ( ii - jj ) for ii,jj in zip( self.coor[self.dime*who:self.dime*who+self.dime], ref ) ]
					tau = [ ii-jj for ii,jj in zip( self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime], self.crd0 ) ]
				elif( who == self.node - 1 ):
					# last node
					ref = [ ( ii + jj ) / 2.0 for ii,jj in zip( self.crdf, self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime] ) ]
					gum = [ self.kumb * ( ii - jj ) for ii,jj in zip( self.coor[self.dime*who:self.dime*who+self.dime], ref ) ]
					tau = [ ii-jj for ii,jj in zip( self.crdf, self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime] ) ]
				else:
					# intermediates
					ref = [ ( ii + jj ) / 2.0 for ii,jj in zip( self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime],
							self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime] ) ]
					gum = [ self.kumb * ( ii - jj ) for ii,jj in zip( self.coor[self.dime*who:self.dime*who+self.dime], ref ) ]
					if( who >= 2 and who < self.node - 2 ):
						tau = []
						for j in range( self.dime ):
							tau.append( ( -2.0 * self.coor[self.dime*(who-2)+j] - self.coor[self.dime*(who-1)+j]
										+ self.coor[self.dime*(who+1)+j] + 2.0 *  self.coor[self.dime*(who+2)+j] ) )
					else:
						tau = [ ii-jj for ii,jj in zip( self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime],
								self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime] ) ]
				# common to all nodes
				tmp = math.sqrt( sum( [ ii * ii for ii in tau ] ) )
				tau = [ ii / tmp for ii in tau ]
				self.prob.coor = self.coor[self.dime*who:self.dime*who+self.dime][:]
				self.prob.get_grad()
				self.func += self.prob.func
				try:
					self.prob.neb_data( who )
				except:
					pass
				usp = sum( [ ii * jj for ii,jj in zip( tau, gum ) ] )
				gsp = sum( [ ii * jj for ii,jj in zip( tau, self.prob.grad ) ] )
				grd = [ ii - gsp * jj + usp * jj for ii,jj in zip( self.prob.grad, tau ) ]
				qm3.utils.project_RT_modes( self.mass, self.prob.coor, grd, None )
				self.grad += grd[:]
			# sync func from and to nodes
			qm3.utils._mpi.barrier()
			if( self.wami == 0 ):
				for i in range( 1, self.ncpu ):
					self.func += qm3.utils._mpi.recv_r8( i, 1 )[0]
				for i in range( 1, self.ncpu ):
					qm3.utils._mpi.send_r8( i, [ self.func ] )
			else:
				qm3.utils._mpi.send_r8( 0, [ self.func ] )
				self.func = qm3.utils._mpi.recv_r8( 0, 1 )[0]
			# sync grad from and to nodes
			qm3.utils._mpi.barrier()
			if( self.wami == 0 ):
				for i in range( 1, self.ncpu ):
					self.grad += qm3.utils._mpi.recv_r8( i, self.dime * self.chnk )
				for i in range( 1, self.ncpu ):
					qm3.utils._mpi.send_r8( i, self.grad )
			else:
				qm3.utils._mpi.send_r8( 0, self.grad )
				self.grad = qm3.utils._mpi.recv_r8( 0, self.size )


		def get_grad( self ):
			# ----------------------------------------------------------------------
			def __calc_tau( potm, poti, potp, crdm, crdi, crdp ):
				dcM = [ ii-jj for ii,jj in zip( crdp, crdi ) ]
				dcm = [ ii-jj for ii,jj in zip( crdi, crdm ) ]
				dpM = max( math.fabs( potp - poti ), math.fabs( potm - poti ) )
				dpm = min( math.fabs( potp - poti ), math.fabs( potm - poti ) )
				if( potp > poti and poti > potm ):
					tau = dcM[:]
				elif( potp < poti and poti < potm ):
					tau = dcm[:]
				else:
					if( potp > potm ):
						tau = [ dpM * ii + dpm * jj for ii,jj in zip( dcM, dcm ) ]
					else:
						tau = [ dpm * ii + dpM * jj for ii,jj in zip( dcM, dcm ) ]
				tmp = math.sqrt( sum( [ ii * ii for ii in tau ] ) )
				tau = [ ii / tmp for ii in tau ]
				mcM = math.sqrt( sum( [ ii * ii for ii in dcM ] ) )
				mcm = math.sqrt( sum( [ ii * ii for ii in dcm ] ) )
				gum = [ - self.kumb * ( mcM - mcm ) * ii for ii in tau ]
				return( tau, gum )
			# ----------------------------------------------------------------------
			vpot = []
			gpot = []
			# get potential energy and gradients from the chain elements
			for who in range( self.wami * self.chnk, ( self.wami + 1 ) * self.chnk ):
				self.prob.coor = self.coor[self.dime*who:self.dime*who+self.dime][:]
				self.prob.get_grad()
				vpot.append( self.prob.func )
				gpot += self.prob.grad[:]
				try:
					self.prob.neb_data( who )
				except:
					pass
			# sync vpot from and to nodes
			qm3.utils._mpi.barrier()
			if( self.wami == 0 ):
				for i in range( 1, self.ncpu ):
					vpot += qm3.utils._mpi.recv_r8( i, self.chnk )
				for i in range( 1, self.ncpu ):
					qm3.utils._mpi.send_r8( i, vpot )
			else:
				qm3.utils._mpi.send_r8( 0, vpot )
				vpot = qm3.utils._mpi.recv_r8( 0, self.node )
			# sync gpot from and to nodes
			qm3.utils._mpi.barrier()
			if( self.wami == 0 ):
				for i in range( 1, self.ncpu ):
					gpot += qm3.utils._mpi.recv_r8( i, self.dime * self.chnk )
				for i in range( 1, self.ncpu ):
					qm3.utils._mpi.send_r8( i, gpot )
			else:
				qm3.utils._mpi.send_r8( 0, gpot )
				gpot = qm3.utils._mpi.recv_r8( 0, self.size )
			self.func = sum( vpot )
			self.grad = []
			# calculate neb components	
			for who in range( self.node ):
				if( who == 0 ):
					# first node
					tau, gum = __calc_tau( self.pot0, vpot[who], vpot[who+1],
								self.crd0,
								self.coor[self.dime*who:self.dime*who+self.dime],
								self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime] )
				elif( who == self.node - 1 ):
					# last node
					tau, gum = __calc_tau( vpot[who-1], vpot[who], self.potf,
								self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime],
								self.coor[self.dime*who:self.dime*who+self.dime],
								self.crdf )
				else:
					# intermediates
					tau, gum = __calc_tau( vpot[who-1], vpot[who], vpot[who+1],
								self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime],
								self.coor[self.dime*who:self.dime*who+self.dime],
								self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime] )
				# common to all nodes
				gsp = sum( [ ii * jj for ii,jj in zip( tau, gpot[self.dime*who:self.dime*who+self.dime] ) ] )
				grd = [ ii - gsp * jj + kk for ii,jj,kk in zip( gpot[self.dime*who:self.dime*who+self.dime], tau, gum ) ]
				qm3.utils.project_RT_modes( self.mass, self.coor[self.dime*who:self.dime*who+self.dime], grd, None )
				self.grad += grd[:]
	
except:
	pass


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
import	qm3.actions.genetic
#import	qm3.utils._mpi
try:
	import	cPickle as pickle
except:
	import	pickle



def guess_bond( mol ):
	bond = []
	for i in range( mol.natm - 1 ):
		i3 = i * 3
		ri = qm3.elements.r_cov[mol.anum[i]] + 0.05
		for j in range( i + 1, mol.natm ):
			if( mol.anum[i] == 1 and mol.anum[j] == 1 ):
				continue
			rj = qm3.elements.r_cov[mol.anum[j]] + 0.05
			t  = ( ri + rj ) * ( ri + rj )
			if( qm3.utils.distanceSQ( mol.coor[i3:i3+3], mol.coor[j*3:j*3+3] ) <= t ):
				bond.append( [ i, j ] )
	return( bond )


def guess_angl( bond ):
	angl = []
	for i in range( len( bond ) - 1 ):
		for j in range( i + 1, len( bond ) ):
			if( bond[i][0] == bond[j][0] ):
				angl.append( [ bond[i][1], bond[i][0], bond[j][1] ] )
			elif( bond[i][0] == bond[j][1] ):
				angl.append( [ bond[i][1], bond[i][0], bond[j][0] ] )
			elif( bond[i][1] == bond[j][0] ):
				angl.append( [ bond[i][0], bond[i][1], bond[j][1] ] )
			elif( bond[i][1] == bond[j][1] ):
				angl.append( [ bond[i][0], bond[i][1], bond[j][0] ] )
	return( angl )


def guess_dihe( angl ):
	dihe = []
	for i in range( len( angl ) - 1 ):
		for j in range( i + 1, len( angl ) ):
			if( angl[i][1] == angl[j][0] and angl[i][2] == angl[j][1] ):
				dihe.append( [ angl[i][0], angl[i][1], angl[i][2], angl[j][2] ] )
			elif( angl[i][1] == angl[j][2] and angl[i][2] == angl[j][1] ):
				dihe.append( [ angl[i][0], angl[i][1], angl[i][2], angl[j][0] ] )
			elif( angl[i][1] == angl[j][0] and angl[i][0] == angl[j][1] ):
				dihe.append( [ angl[i][2], angl[i][1], angl[i][0], angl[j][2] ] )
			elif( angl[i][1] == angl[j][2] and angl[i][0] == angl[j][1] ):
				dihe.append( [ angl[i][2], angl[i][1], angl[i][0], angl[j][0] ] )
	return( dihe )


class auto_parametrize( qm3.problem.template ):
	"""
		bond:	fc * ( r - r0 ) ^ 2
				[ ( i, j ), ... ]

		angl:	fc * ( t - t0 ) ^ 2
				[ ( i, j, k ), ... ]

		dihe:	fc1 * ( 1 + cos( n=1 * p0 - d=0   ) ) + 
				fc2 * ( 1 + cos( n=2 * p0 - d=180 ) ) + 
				fc3 * ( 1 + cos( n=3 * p0 - d=0   ) )
				[ ( i, j, k, l ), ... ]

		impr:	fc * ( x - x0 ) ^ 2
				[ ( i_cent, j, k, l ), ... ]

		lenj:	epsi_ij * ( ( rmin_ij / r_ij ) ^ 12 - 2 * ( rmin_ij / r_ij ) ^ 6 )
				epsi_ij = Sqrt( epsi_i * epsi_j )
				rmin_ij = rmin_i/2 + rmin_j/2

		coul:	chrg_i * chrg_j / r_ij

		reference values for each type (r0, t0, p0, x0) are derived from the initial molecule...
		force constans in kJ/mol
	"""
	def __init__( self, mole, bond, angl, dihe, impr, hess ):
		qm3.problem.template.__init__( self )
		self.mole = mole
		self.bond = bond
		self.angl = angl
		self.dihe = dihe
		self.impr = impr
		self.targ = hess[:]

		self.n3 = 3 * self.mole.natm

		self.size  = len( self.bond ) + len( self.angl ) + len( self.dihe ) * 3 + len( self.impr )
		self.coor  = [ 1000.0 for i in range( len( self.bond ) ) ]
		self.coor += [  200.0 for i in range( len( self.angl ) ) ]
		self.coor += [    5.0 for i in range( len( self.dihe ) * 3 ) ]
		self.coor += [  100.0 for i in range( len( self.impr ) ) ]
		self.func  = 0.0
		self.grad  = []

		self.boun  = [ ( 500.0, 4000.0 ) for i in range( len( self.bond ) ) ]
		self.boun += [ (  50.0, 1000.0 ) for i in range( len( self.angl ) ) ]
		self.boun += [ (   0.0,   50.0 ) for i in range( len( self.dihe ) ) ]
		self.boun += [ (   0.0,   50.0 ) for i in range( len( self.dihe ) ) ]
		self.boun += [ (   0.0,   50.0 ) for i in range( len( self.dihe ) ) ]
		self.boun += [ (   0.0,  500.0 ) for i in range( len( self.impr ) ) ]

		self.nb14 = [ ( self.dihe[i][0], self.dihe[i][3] ) for i in range( len( self.dihe ) ) ]
		self.nbnd = []
		for i in range( self.mole.natm - 1 ):
			for j in range( i + 1, self.mole.natm ):
				xx = False
				for ii,jj in self.bond:
					xx |= ( ( i == ii ) and ( j == jj ) ) or ( ( i == jj ) and ( j == ii ) )
				for ii,jj,kk in self.angl:
					xx |= ( ( i == ii ) and ( j == kk ) ) or ( ( i == kk ) and ( j == ii ) )
				for ii,jj,kk,ll in self.dihe:
					xx |= ( ( i == ii ) and ( j == ll ) ) or ( ( i == ll ) and ( j == ii ) )
				if( not xx ):
					self.nbnd.append( ( i, j ) )

		self.r_bond = []
		for i,j in self.bond:
			self.r_bond.append( qm3.utils.distance( self.mole.coor[3*i:3*i+3], self.mole.coor[3*j:3*j+3] ) )
		self.r_angl = []
		for i,j,k in self.angl:
			self.r_angl.append( qm3.utils.angleRAD( self.mole.coor[3*i:3*i+3], self.mole.coor[3*j:3*j+3],
								self.mole.coor[3*k:3*k+3] ) )
		self.r_impr = []
		for i,j,k,l in self.impr:
			self.r_impr.append( qm3.utils.dihedral( self.mole.coor[3*i:3*i+3], self.mole.coor[3*j:3*j+3],
								self.mole.coor[3*k:3*k+3], self.mole.coor[3*l:3*l+3] ) )


	def mole_grad( self ):
		grad = [ 0.0 for i in range( self.n3 ) ]
		who  = 0
		# bond
		for w in range( len( self.bond ) ):
			i3  = 3 * self.bond[w][0]
			j3  = 3 * self.bond[w][1]
			ref = self.r_bond[w]
			dij = [ ii-jj for ii,jj in zip( self.mole.coor[i3:i3+3], self.mole.coor[j3:j3+3] ) ]
			val = math.sqrt( dij[0] * dij[0] + dij[1] * dij[1] + dij[2] * dij[2] )
			der = 2.0 * self.coor[who] * ( val - ref ) / val
			for m in [ 0, 1, 2 ]:
				grad[i3+m] += der * dij[m]
				grad[j3+m] -= der * dij[m]
			who += 1
		# angl
		for w in range( len( self.angl ) ):
			i3  = 3 * self.angl[w][0]
			j3  = 3 * self.angl[w][1]
			k3  = 3 * self.angl[w][2]
			ref = self.r_angl[w]
			dij = [ ii-jj for ii,jj in zip( self.mole.coor[i3:i3+3], self.mole.coor[j3:j3+3] ) ]
			rij = math.sqrt( sum( [ t*t for t in dij ] ) )
			dij = [ t / rij for t in dij ]
			dkj = [ ii-jj for ii,jj in zip( self.mole.coor[k3:k3+3], self.mole.coor[j3:j3+3] ) ]
			rkj = math.sqrt( sum( [ t*t for t in dkj ] ) )
			dkj = [ t / rkj for t in dkj ]
			fac = sum( [ ii*jj for ii,jj in zip( dij, dkj ) ] )
			fac = min( math.fabs( fac ), 1.0 - 1.0e-6 ) * fac / math.fabs( fac )
			val = math.acos( fac )
			dtx =  - 1.0 / math.sqrt( 1.0 - fac * fac )
			der = 2.0 * self.coor[who] * ( val - ref ) / dtx
			dti = [ ( ii - fac * jj ) / rij for ii,jj in zip( dkj, dij ) ]
			dtk = [ ( ii - fac * jj ) / rkj for ii,jj in zip( dij, dkj ) ]
			dtj = [ - ( ii + jj ) for ii,jj in zip( dti, dtk ) ]
			for m in [ 0, 1, 2 ]:
				grad[i3+m] += der * dti[m]
				grad[j3+m] += der * dtj[m]
				grad[k3+m] += der * dtk[m]
			who += 1
		# dihe
		for w in range( len( self.dihe ) ):
			i3  = 3 * self.dihe[w][0]
			j3  = 3 * self.dihe[w][1]
			k3  = 3 * self.dihe[w][2]
			l3  = 3 * self.dihe[w][3]
			dji = [ ii-jj for ii,jj in zip( self.mole.coor[j3:j3+3], self.mole.coor[i3:i3+3] ) ]
			dkj = [ ii-jj for ii,jj in zip( self.mole.coor[k3:k3+3], self.mole.coor[j3:j3+3] ) ]
			rkj = math.sqrt( sum( [ t*t for t in dkj ] ) )
			dlk = [ ii-jj for ii,jj in zip( self.mole.coor[l3:l3+3], self.mole.coor[k3:k3+3] ) ]
			vt  = [ dji[1] * dkj[2] - dkj[1] * dji[2], dji[2] * dkj[0] - dkj[2] * dji[0], dji[0] * dkj[1] - dkj[0] * dji[1] ]
			rt2 = sum( [ t*t for t in vt ] )
			vu  = [ dkj[1] * dlk[2] - dlk[1] * dkj[2], dkj[2] * dlk[0] - dlk[2] * dkj[0], dkj[0] * dlk[1] - dlk[0] * dkj[1] ]
			ru2 = sum( [ t*t for t in vu ] )
			vtu = [ vt[1] * vu[2] - vu[1] * vt[2], vt[2] * vu[0] - vu[2] * vt[0], vt[0] * vu[1] - vu[0] * vt[1] ]
			rtu = math.sqrt( rt2 * ru2 )
			if( rtu == 0.0 ):
				who += 3
				continue
			cs1 = sum( [ ii*jj for ii,jj in zip( vt, vu ) ] ) / rtu
			sn1 = sum( [ ii*jj for ii,jj in zip( dkj, vtu ) ] ) / ( rkj * rtu )
			cs2 = cs1 * cs1 - sn1 * sn1
			sn2 = 2.0 * cs1 * sn1
			sn3 = cs1 * sn2 + sn1 * cs2
			dph  = - self.coor[who]   *       sn1
			dph +=   self.coor[who+1] * 2.0 * sn2
			dph += - self.coor[who+2] * 3.0 * sn3
			dki = [ ii-jj for ii,jj in zip( self.mole.coor[k3:k3+3], self.mole.coor[i3:i3+3] ) ]
			dlj = [ ii-jj for ii,jj in zip( self.mole.coor[l3:l3+3], self.mole.coor[j3:j3+3] ) ]
			dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
					( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
					( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
			dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
					( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
					( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
			grad[i3]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph
			grad[i3+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph
			grad[i3+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph
			grad[j3]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph
			grad[j3+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph
			grad[j3+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph
			grad[k3]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph
			grad[k3+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph
			grad[k3+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph
			grad[l3]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph
			grad[l3+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph
			grad[l3+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph
			who += 3
		# impr
		for w in range( len( self.impr ) ):
			i3  = 3 * self.impr[w][0]
			j3  = 3 * self.impr[w][1]
			k3  = 3 * self.impr[w][2]
			l3  = 3 * self.impr[w][3]
			ref = self.r_impr[w]
			dji = [ ii-jj for ii,jj in zip( self.mole.coor[j3:j3+3], self.mole.coor[i3:i3+3] ) ]
			dkj = [ ii-jj for ii,jj in zip( self.mole.coor[k3:k3+3], self.mole.coor[j3:j3+3] ) ]
			rkj = math.sqrt( sum( [ t*t for t in dkj ] ) )
			dlk = [ ii-jj for ii,jj in zip( self.mole.coor[l3:l3+3], self.mole.coor[k3:k3+3] ) ]
			vt  = [ dji[1] * dkj[2] - dkj[1] * dji[2], dji[2] * dkj[0] - dkj[2] * dji[0], dji[0] * dkj[1] - dkj[0] * dji[1] ]
			rt2 = sum( [ t*t for t in vt ] )
			vu  = [ dkj[1] * dlk[2] - dlk[1] * dkj[2], dkj[2] * dlk[0] - dlk[2] * dkj[0], dkj[0] * dlk[1] - dlk[0] * dkj[1] ]
			ru2 = sum( [ t*t for t in vu ] )
			vtu = [ vt[1] * vu[2] - vu[1] * vt[2], vt[2] * vu[0] - vu[2] * vt[0], vt[0] * vu[1] - vu[0] * vt[1] ]
			rtu = math.sqrt( rt2 * ru2 )
			if( rtu == 0.0 ):
				who += 1
				continue
			cos = sum( [ ii*jj for ii,jj in zip( vt, vu ) ] ) / rtu
			sin = sum( [ ii*jj for ii,jj in zip( dkj, vtu ) ] ) / ( rkj * rtu )
			cos = min( 1.0, max( -1.0, cos ) )
			ang = qm3.constants.R2D * math.acos( cos )
			if( sin <= 0.0 ):
				ang = -ang
			if( math.fabs( ang + ref ) < math.fabs( ang - ref ) ):
				ref = -ref
			dt  = ang - ref
			while( dt >  180.0 ):
				dt -= 360.0
			while( dt < -180.0 ):
				dt += 360.0
			dt  /= qm3.constants.R2D
			dph = 2.0 * self.coor[who] * dt
			dki = [ ii-jj for ii,jj in zip( self.mole.coor[k3:k3+3], self.mole.coor[i3:i3+3] ) ]
			dlj = [ ii-jj for ii,jj in zip( self.mole.coor[l3:l3+3], self.mole.coor[j3:j3+3] ) ]
			dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
					( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
					( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
			dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
					( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
					( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
			grad[i3]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph
			grad[i3+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph
			grad[i3+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph
			grad[j3]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph
			grad[j3+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph
			grad[j3+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph
			grad[k3]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph
			grad[k3+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph
			grad[k3+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph
			grad[l3]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph
			grad[l3+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph
			grad[l3+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph
			who += 1
		# nbnd
		eps = 1389.35484620709144110151
		for i,j in self.nbnd:
			i3  = i * 3
			j3  = j * 3
			dr  = [ ii-jj for ii,jj in zip( self.mole.coor[i3:i3+3], self.mole.coor[j3:j3+3] ) ]
			r2  = sum( [ t*t for t in dr ] )
			s   = 1.0 / math.sqrt( r2 )
			s6  = math.pow( ( self.mole.rmin[i] + self.mole.rmin[j] ) * s, 6.0 )
			df  = 12.0 * self.mole.epsi[i] * self.mole.epsi[j] * s6 * ( 1.0 - s6 ) / r2
			df -= eps * self.mole.chrg[i] * self.mole.chrg[j] * s / r2
			for m in [ 0, 1, 2 ]:
				grad[i3+m] += df * dr[m]
				grad[j3+m] -= df * dr[m]
		# nb14
		for i,j in self.nb14:
			i3  = i * 3
			j3  = j * 3
			dr  = [ ii-jj for ii,jj in zip( self.mole.coor[i3:i3+3], self.mole.coor[j3:j3+3] ) ]
			r2  = sum( [ t*t for t in dr ] )
			s   = 1.0 / math.sqrt( r2 )
			s6  = math.pow( ( self.mole.rmin[i] + self.mole.rmin[j] ) * s, 6.0 )
			df  = 12.0 * self.mole.epsi[i] * self.mole.epsi[j] * s6 * ( 1.0 - s6 ) / r2
			df -= eps * self.mole.chrg[i] * self.mole.chrg[j] * s / r2
			for m in [ 0, 1, 2 ]:
				grad[i3+m] += df * dr[m] * 0.5
				grad[j3+m] -= df * dr[m] * 0.5
		return( grad )


	def mole_hess( self ):
		hh = []
		for i in range( self.mole.natm ):
			i3 = i * 3
			for j in [ 0, 1, 2 ]:
				coor = self.mole.coor[i3+j]
				self.mole.coor[i3+j] = coor + 1.0e-4
				g_for = self.mole_grad()
				self.mole.coor[i3+j] = coor - 1.0e-4
				g_bak = self.mole_grad()
				self.mole.coor[i3+j] = coor
				hh.append( [ ( g_for[k] - g_bak[k] ) / 2.0e-4 for k in range( self.n3 ) ] )
		hess = []
		for i in range( self.n3 ):
			for j in range( self.n3 ):
				hess.append( ( hh[i][j] + hh[j][i] ) * 0.5 )
		return( hess )


	def get_func( self ):
		hess = self.mole_hess()
		self.func = 0.0
		for i in range( self.n3 ):
			for j in range( i, self.n3 ):
				self.func += math.pow( hess[i*self.n3+j] - self.targ[i*self.n3+j], 2.0 )


	def get_grad( self ):
		self.num_grad()


	def minimize( self ):
		nn = 1.0 / math.sqrt( self.size )
		# ----------------------------------------------------------------------------------
		ss = 100.
		tt = ss * 10.
		it = 0
		self.get_grad()
		gg = math.sqrt( sum( [ t*t for t in self.grad ] ) )
		gr = gg * nn
		bb = [ round( self.func, 0 ) ]
		print( 80 * "-" )
		print( "%-20d%20.8lf%20.8lf%20.8lf"%( it, self.func, gr, ss ) )
		while( it < 1000 and ( gr > tt or self.func > tt  ) and ss > 1.0 / tt ):
			self.coor = [ self.coor[i] - self.grad[i] / gg * ss for i in range( self.size ) ]
			self.get_grad()
			gg = math.sqrt( sum( [ t*t for t in self.grad ] ) )
			gr = gg * nn
			ff = round( self.func, 0 )
			if( not ff in bb ):
				bb.append( ff )
			else:
				ss *= 0.75
			it += 1
#			if( it % 100 == 0 ):
			print( "%-20d%20.8lf%20.8lf%20.8lf"%( it, self.func, gr, ss ) )
#		if( it % 100 != 0 ):
#			print( "%-20d%20.8lf%20.8lf%20.8lf"%( it, self.func, gr, ss ) )
		# ----------------------------------------------------------------------------------


	def flush_parm( self ):
		# ---------------------------------------
		self.mole.type = []
		for i in range( self.mole.natm ):
			self.mole.type.append( self.mole.labl[i] + str( i + 1 ) )
		# ---------------------------------------
		f = open( "parm", "wt" )
		for i in range( self.mole.natm ):
			f.write( "%8s%12.3lf%12.3lf\n"%( self.mole.type[i], self.mole.epsi[i] * self.mole.epsi[i] / qm3.constants.K2J , self.mole.rmin[i] ) )	
		f.write( "\n" )
		c =0
		for i in range( len( self.bond ) ):
			f.write( "%8s%8s%12.1lf%12.3lf\n"%( self.mole.type[self.bond[i][0]],
				self.mole.type[self.bond[i][1]], self.coor[c] / qm3.constants.K2J, self.r_bond[i] ) )
			c += 1
		f.write( "\n" )
		for i in range( len( self.angl ) ):
			f.write( "%8s%8s%8s%12.1lf%12.1lf\n"%( self.mole.type[self.angl[i][0]],
				self.mole.type[self.angl[i][1]], self.mole.type[self.angl[i][2]],
				self.coor[c] / qm3.constants.K2J, self.r_angl[i] * qm3.constants.R2D ) )
			c += 1
		f.write( "\n" )
		for i in range( len( self.dihe ) ):
			f.write( "%8s%8s%8s%8s%12.3lf  1    0.0%12.3lf  2  180.0%12.3lf  3    0.0 \n"%(
				self.mole.type[self.dihe[i][0]], self.mole.type[self.dihe[i][1]],
				self.mole.type[self.dihe[i][2]], self.mole.type[self.dihe[i][3]],
				self.coor[c]   / qm3.constants.K2J,
				self.coor[c+1] / qm3.constants.K2J,
				self.coor[c+2] / qm3.constants.K2J ) )
			c += 3
		f.write( "\n" )
		for i in range( len( self.impr ) ):
			f.write( "%8s%8s%8s%8s%12.1lf%12.1lf\n"%(
				self.mole.type[self.impr[i][0]], self.mole.type[self.impr[i][1]],
				self.mole.type[self.impr[i][2]], self.mole.type[self.impr[i][3]],
				self.coor[c] / qm3.constants.K2J, self.r_impr[i] ) )
			c += 1
		f.close()
















node = 0
#node, ncpu = qm3.utils._mpi.init()

mol = qm3.mol.molecule()
mol.xyz_read( "last.xyz" )
mol.anum = [ qm3.elements.rsymbol[mol.labl[i].title()] for i in range( mol.natm ) ]
mol.fill_masses()

f = open( "last.chg", "rb" )
mol.chrg = pickle.load( f ) 
f.close()

mol.epsi = [0.5411838874172069, 0.5411838874172069, 0.5254940532489402, 0.5254940532489402, 0.4709055106919009,
			0.8678248671246982, 0.31688483712541377, 0.6033307550589477, 0.8283574107835338, 0.31688483712541377,
			0.31688483712541377, 0.31688483712541377, 0.0, 0.9373579892442374, 0.9373579892442374, 0.0, 0.0, 0.0,
			0.843374175559105, 0.31688483712541377, 0.0]
mol.rmin = [1.992, 1.992, 1.956, 1.956, 1.675, 1.839, 1.381, 1.824, 1.741, 1.381, 1.381, 1.381, 0.0, 1.661,
			1.661, 0.0, 0.0, 0.0, 1.824, 1.381, 0.0]

f = open( "last.hes", "rb" )
hess = pickle.load( f )
f.close()

bond = guess_bond( mol )
angl = guess_angl( bond )
dihe = guess_dihe( angl )
impr = [ [ 7, 13, 14, 3 ] ]

obj = auto_parametrize( mol, bond, angl, dihe, impr, hess )

#qm3.actions.genetic.diffevo( obj, obj.boun,
#	step_number = 2,
#	print_frequency = 1,
#	population_size = 2 * obj.size )

#qm3.actions.genetic.mpi_diffevo( obj, obj.boun,
#	mpi_node = node, mpi_ncpu = ncpu,
#	step_number = 2,
#	print_frequency = 1,
#	population_size = obj.size )

tmp = [ obj, 
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ),
	auto_parametrize( mol, bond, angl, dihe, impr, hess ) ]
qm3.actions.genetic.smp_diffevo( tmp, obj.boun,
	step_number = 200,
	print_frequency = 1,
	population_size = obj.size * 8 )

if( node == 0 ):
	val, vec = qm3.utils.hessian_frequencies( mol.mass, mol.coor, obj.mole_hess() )
	print( val )
	f = open( "resultados", "wb" )
	pickle.dump( obj.coor, f )
	f.close()
	obj.flush_parm()

#qm3.utils._mpi.barrier()
#qm3.utils._mpi.stop()

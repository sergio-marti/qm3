# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.constants



class distance( object ):

	def __init__( self, kumb, xref, indx, skip_LE = 0.0, skip_BE = 9.e99 ):
		self.kumb = kumb
		self.xref = xref
		self.indx = indx[:]
		self.skpL = skip_LE
		self.skpB = skip_BE


	def get_func( self, molec ):
		ii = 3 * self.indx[0]
		jj = 3 * self.indx[1]
		dr = [ i-j for i,j in zip( molec.coor[ii:ii+3], molec.coor[jj:jj+3] ) ]
		vv = math.sqrt( sum( [ i * i for i in dr ] ) )
		if( vv >= self.skpL and vv <= self.skpB ):
			molec.func += 0.5 * self.kumb * ( vv - self.xref ) * ( vv - self.xref )
		return( vv )


	def get_grad( self, molec ):
		ii = 3 * self.indx[0]
		jj = 3 * self.indx[1]
		dr = [ i-j for i,j in zip( molec.coor[ii:ii+3], molec.coor[jj:jj+3] ) ]
		vv = math.sqrt( sum( [ i * i for i in dr ] ) )
		df = self.kumb * ( vv - self.xref )
		if( vv >= self.skpL and vv <= self.skpB ):
			molec.func += 0.5 * df * ( vv - self.xref )
			df /= vv
			for i in [0, 1, 2]:
				molec.grad[ii+i] += df * dr[i]
				molec.grad[jj+i] -= df * dr[i]
		return( vv )



class angle( object ):

	def __init__( self, kumb, xref, indx ):
		self.kumb = kumb
		self.xref = xref
		self.indx = indx[:]


	def get_func( self, molec ):
		ii = 3 * self.indx[0]
		jj = 3 * self.indx[1]
		kk = 3 * self.indx[2]
		dij = [ i-j for i,j in zip( molec.coor[ii:ii+3], molec.coor[jj:jj+3] ) ]
		rij = math.sqrt( sum( [ i * i for i in dij ] ) )
		dij = [ i / rij for i in dij ]
		dkj = [ k-j for k,j in zip( molec.coor[kk:kk+3], molec.coor[jj:jj+3] ) ]
		rkj = math.sqrt( sum( [ i * i for i in dkj ] ) )
		dkj = [ i / rkj for i in dkj ]
		dot = sum( [ i * j for i,j in zip( dij, dkj ) ] )
		if( dot >= 0.0 ):
			dot = min( dot, 1.0 - 1.0e-6 )
		else:
			dot = - min( math.fabs( dot ), 1.0 - 1.0e-6 )
		vv = math.acos( dot ) * qm3.constants.R2D
		molec.func += 0.5 * self.kumb * ( vv - self.xref ) * ( vv - self.xref )
		return( vv )


	def get_grad( self, molec ):
		ii = 3 * self.indx[0]
		jj = 3 * self.indx[1]
		kk = 3 * self.indx[2]
		dij = [ i-j for i,j in zip( molec.coor[ii:ii+3], molec.coor[jj:jj+3] ) ]
		rij = math.sqrt( sum( [ i * i for i in dij ] ) )
		dij = [ i / rij for i in dij ]
		dkj = [ k-j for k,j in zip( molec.coor[kk:kk+3], molec.coor[jj:jj+3] ) ]
		rkj = math.sqrt( sum( [ i * i for i in dkj ] ) )
		dkj = [ i / rkj for i in dkj ]
		dot = sum( [ i * j for i,j in zip( dij, dkj ) ] )
		if( dot >= 0.0 ):
			dot = min( dot, 1.0 - 1.0e-6 )
		else:
			dot = - min( math.fabs( dot ), 1.0 - 1.0e-6 )
		vv = math.acos( dot ) * qm3.constants.R2D
		df = self.kumb * ( vv - self.xref )
		molec.func += 0.5 * df * ( vv - self.xref )
		df = - qm3.constants.R2D * df / math.sqrt( 1.0 - dot * dot )
		dti = [ ( dkj[i] - dot * dij[i] ) / rij for i in [ 0, 1, 2 ] ]
		dtk = [ ( dij[i] - dot * dkj[i] ) / rkj for i in [ 0, 1, 2 ] ]
		dtj = [ - ( dti[i] + dtk[i] ) for i in [ 0, 1, 2 ] ]
		for i in [0, 1, 2]:
			molec.grad[ii+i] += df * dti[i]
			molec.grad[jj+i] += df * dtj[i]
			molec.grad[kk+i] += df * dtk[i]
		return( vv )



class multiple_distance( object ):

	def __init__( self, kumb, xref, indx, weigh ):
		if( len( weigh ) * 2 != len( indx ) ):
			print( "- restraints.multiple_distance: Number of ATOMS should be TWICE the number of WEIGHTS!" )
			return( None )
		self.kumb = kumb
		self.xref = xref
		self.indx = indx[:]
		self.weig = weigh[:]
		self.size = len( weigh )


	def get_func( self, molec ):
		dr = []
		r  = []
		for i in range( self.size ):
			i3 = i * 3
			ii = 3 * self.indx[2*i]
			jj = 3 * self.indx[2*i+1]
			dr += [ j-k for j,k in zip( molec.coor[ii:ii+3], molec.coor[jj:jj+3] ) ]
			r.append( math.sqrt( sum( [ j * j for j in dr[i3:i3+3] ] ) ) )
		vv = sum( [ i * j for i,j in zip( r, self.weig ) ] )
		df = self.kumb * ( vv - self.xref )
		molec.func += 0.5 * df * ( vv - self.xref )
		return( vv )


	def get_grad( self, molec ):
		dr = []
		r  = []
		for i in range( self.size ):
			i3 = i * 3
			ii = 3 * self.indx[2*i]
			jj = 3 * self.indx[2*i+1]
			dr += [ j-k for j,k in zip( molec.coor[ii:ii+3], molec.coor[jj:jj+3] ) ]
			r.append( math.sqrt( sum( [ j * j for j in dr[i3:i3+3] ] ) ) )
		vv = sum( [ i * j for i,j in zip( r, self.weig ) ] )
		df = self.kumb * ( vv - self.xref )
		molec.func += 0.5 * df * ( vv - self.xref )
		for i in range( self.size ):
			i3 = i * 3
			ii = 3 * self.indx[2*i]
			jj = 3 * self.indx[2*i+1]
			tt = self.weig[i] * df / r[i]
			for j in [0, 1, 2]:
				molec.grad[ii+j] += tt * dr[i3+j]
				molec.grad[jj+j] -= tt * dr[i3+j]
		return( vv )



class tether( object ):

	def __init__( self, molec, kumb, indx ):
		self.kumb = kumb
		self.indx = {}
		for i in indx:
			self.indx[i] = molec.coor[3*i:3*i+3][:]


	def get_func( self, molec ):
		for w in self.indx:
			w3 = w * 3
			dr = [ i-j for i,j in zip( molec.coor[w3:w3+3], self.indx[w] ) ]
			molec.func += 0.5 * self.kumb * sum( [ i * i for i in dr ] )
		return( None )


	def get_grad( self, molec ):
		for w in self.indx:
			w3 = w * 3
			dr = [ i-j for i,j in zip( molec.coor[w3:w3+3], self.indx[w] ) ]
			molec.func += 0.5 * self.kumb * sum( [ i * i for i in dr ] ) 
			for i in [0, 1, 2]:
				molec.grad[w3+i] += self.kumb * dr[i]
		return( None )



class dihedral( object ):
	
	def __init__( self, molec, data, indx ):
	"""
	dihedral = force_constant * ( 1 + cos( periodicity * angle - displacement ) )

	data = {  periodicity: [ force_constant [kJ/mol], displacement [degrees] ], ... }

	X - C_sp3 - C_sp3 - X   =>  { 1: [ -0.3347, 0.0 ], 2: [ -0.5858, 84.8 ], 3: [ 1.3263, 0.0 ] }

	valid periodicity = [ 1 : 6 ]
	"""
	self.data = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	for i in range( 6 ):
		if( i+1 in data ):
			self.data[2*i]   = data[i+1][0]
			self.data[2*i+1] = data[i+1][1] / qm3.constants.R2D


	def get_func( self, molec ):
		ai  = 3 * self.indx[0]
		aj  = 3 * self.indx[1]
		ak  = 3 * self.indx[2]
		al  = 3 * self.indx[3]
		dji = [ ii-jj for ii,jj in zip( molec.coor[aj:aj+3], molec.coor[ai:ai+3] ) ]
		dkj = [ ii-jj for ii,jj in zip( molec.coor[ak:ak+3], molec.coor[aj:aj+3] ) ]
		rkj = math.sqrt( sum( [ ii*ii for ii in dkj ] ) )
		dlk = [ ii-jj for ii,jj in zip( molec.coor[al:al+3], molec.coor[ak:ak+3] ) ]
		vt  = [ dji[1] * dkj[2] - dkj[1] * dji[2], dji[2] * dkj[0] - dkj[2] * dji[0], dji[0] * dkj[1] - dkj[0] * dji[1] ]
		rt2 = sum( [ ii*ii for ii in vt ] )
		vu  = [ dkj[1] * dlk[2] - dlk[1] * dkj[2], dkj[2] * dlk[0] - dlk[2] * dkj[0], dkj[0] * dlk[1] - dlk[0] * dkj[1] ]
		ru2 = sum( [ ii*ii for ii in vu ] )
		vtu = [ vt[1] * vu[2] - vu[1] * vt[2], vt[2] * vu[0] - vu[2] * vt[0], vt[0] * vu[1] - vu[0] * vt[1] ]
		rtu = math.sqrt( rt2 * ru2 )
		if( rtu == 0.0 ):
			continue
		cs1 = sum( [ ii*jj for ii,jj in zip( vt, vu ) ] ) / rtu
		sn1 = sum( [ ii*jj for ii,jj in zip( dkj, vtu ) ] ) / ( rkj * rtu )
		cs2 = cs1 * cs1 - sn1 * sn1
		sn2 = 2.0 * cs1 * sn1
		cs3 = cs1 * cs2 - sn1 * sn2
		sn3 = cs1 * sn2 + sn1 * cs2
		cs4 = cs1 * cs3 - sn1 * sn3
		sn4 = cs1 * sn3 + sn1 * cs3
		cs5 = cs1 * cs4 - sn1 * sn4
		sn5 = cs1 * sn4 + sn1 * cs4
		cs6 = cs1 * cs5 - sn1 * sn5
		sn6 = cs1 * sn5 + sn1 * cs5
		if( self.data[0] != 0.0 ):
			cd = math.cos( self.data[1] )
			sd = math.sin( self.data[1] )
			molec.func += self.data[0] * ( 1.0 + cs1 * cd + sn1 * sd )
		if( self.data[2] != 0.0 ):
			cd = math.cos( self.data[3] )
			sd = math.sin( self.data[3] )
			molec.func += self.data[2] * ( 1.0 + cs2 * cd + sn2 * sd )
		if( self.data[4] != 0.0 ):
			cd = math.cos( self.data[5] )
			sd = math.sin( self.data[5] )
			molec.func += self.data[4] * ( 1.0 + cs3 * cd + sn3 * sd )
		if( self.data[6] != 0.0 ):
			cd = math.cos( self.data[7] )
			sd = math.sin( self.data[7] )
			molec.func += self.data[6] * ( 1.0 + cs4 * cd + sn4 * sd )
		if( self.data[8] != 0.0 ):
			cd = math.cos( self.data[9] )
			sd = math.sin( self.data[9] )
			molec.func += self.data[8] * ( 1.0 + cs5 * cd + sn5 * sd )
		if( self.data[10] != 0.0 ):
			cd = math.cos( self.data[11] )
			sd = math.sin( self.data[11] )
			molec.func += self.data[10] * ( 1.0 + cs6 * cd + sn6 * sd )


	def get_grad( self, molec ):
		ai  = 3 * self.indx[0]
		aj  = 3 * self.indx[1]
		ak  = 3 * self.indx[2]
		al  = 3 * self.indx[3]
		dji = [ ii-jj for ii,jj in zip( molec.coor[aj:aj+3], molec.coor[ai:ai+3] ) ]
		dkj = [ ii-jj for ii,jj in zip( molec.coor[ak:ak+3], molec.coor[aj:aj+3] ) ]
		rkj = math.sqrt( sum( [ ii*ii for ii in dkj ] ) )
		dlk = [ ii-jj for ii,jj in zip( molec.coor[al:al+3], molec.coor[ak:ak+3] ) ]
		vt  = [ dji[1] * dkj[2] - dkj[1] * dji[2], dji[2] * dkj[0] - dkj[2] * dji[0], dji[0] * dkj[1] - dkj[0] * dji[1] ]
		rt2 = sum( [ ii*ii for ii in vt ] )
		vu  = [ dkj[1] * dlk[2] - dlk[1] * dkj[2], dkj[2] * dlk[0] - dlk[2] * dkj[0], dkj[0] * dlk[1] - dlk[0] * dkj[1] ]
		ru2 = sum( [ ii*ii for ii in vu ] )
		vtu = [ vt[1] * vu[2] - vu[1] * vt[2], vt[2] * vu[0] - vu[2] * vt[0], vt[0] * vu[1] - vu[0] * vt[1] ]
		rtu = math.sqrt( rt2 * ru2 )
		if( rtu == 0.0 ):
			continue
		cs1 = sum( [ ii*jj for ii,jj in zip( vt, vu ) ] ) / rtu
		sn1 = sum( [ ii*jj for ii,jj in zip( dkj, vtu ) ] ) / ( rkj * rtu )
		cs2 = cs1 * cs1 - sn1 * sn1
		sn2 = 2.0 * cs1 * sn1
		cs3 = cs1 * cs2 - sn1 * sn2
		sn3 = cs1 * sn2 + sn1 * cs2
		cs4 = cs1 * cs3 - sn1 * sn3
		sn4 = cs1 * sn3 + sn1 * cs3
		cs5 = cs1 * cs4 - sn1 * sn4
		sn5 = cs1 * sn4 + sn1 * cs4
		cs6 = cs1 * cs5 - sn1 * sn5
		sn6 = cs1 * sn5 + sn1 * cs5
		dph = 0.0
		if( self.data[0] != 0.0 ):
			cd  = math.cos( self.data[1] )
			sd  = math.sin( self.data[1] )
			dph += self.data[0] * ( cs1 * sd - sn1 * cd )
			molec.func += self.data[0] * ( 1.0 + cs1 * cd + sn1 * sd )
		if( self.data[2] != 0.0 ):
			cd  = math.cos( self.data[3] )
			sd  = math.sin( self.data[3] )
			dph += self.data[2] * 2.0 * ( cs2 * sd - sn2 * cd )
			molec.func += self.data[2] * ( 1.0 + cs2 * cd + sn2 * sd )
		if( self.data[4] != 0.0 ):
			cd  = math.cos( self.data[5] )
			sd  = math.sin( self.data[5] )
			dph += self.data[4] * 3.0 * ( cs3 * sd - sn3 * cd )
			molec.func += self.data[4] * ( 1.0 + cs3 * cd + sn3 * sd )
		if( self.data[6] != 0.0 ):
			cd  = math.cos( self.data[7] )
			sd  = math.sin( self.data[7] )
			dph += self.data[6] * 4.0 * ( cs4 * sd - sn4 * cd )
			molec.func += self.data[6] * ( 1.0 + cs4 * cd + sn4 * sd )
		if( self.data[8] != 0.0 ):
			cd  = math.cos( self.data[9] )
			sd  = math.sin( self.data[9] )
			dph += self.data[8] * 5.0 * ( cs5 * sd - sn5 * cd )
			molec.func += self.data[8] * ( 1.0 + cs5 * cd + sn5 * sd )
		if( self.data[10] != 0.0 ):
			cd  = math.cos( self.data[11] )
			sd  = math.sin( self.data[11] )
			dph += self.data[10] * 6.0 * ( cs6 * sd - sn6 * cd )
			molec.func += self.data[10] * ( 1.0 + cs6 * cd + sn6 * sd )
		dki = [ ii-jj for ii,jj in zip( molec.coor[ak:ak+3], molec.coor[ai:ai+3] ) ]
		dlj = [ ii-jj for ii,jj in zip( molec.coor[al:al+3], molec.coor[aj:aj+3] ) ]
		dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
				( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
				( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
		dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
				( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
				( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
		molec.grad[ai]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph
		molec.grad[ai+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph
		molec.grad[ai+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph
		molec.grad[aj]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph
		molec.grad[aj+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph
		molec.grad[aj+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph
		molec.grad[ak]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph
		molec.grad[ak+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph
		molec.grad[ak+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph
		molec.grad[al]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph
		molec.grad[al+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph
		molec.grad[al+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph

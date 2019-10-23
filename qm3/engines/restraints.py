# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.constants


def mm_bond( molec, kumb, xref, a_i, a_j, skip_LE = 0.0, skip_BE = 9.e99,
			ffac = 1.0, grad = False, gfac = [ 1.0, 1.0 ], hess = False, hfac = None ):
	"""
	bond = force_constant * ( distance - reference )^2

	force_constant [kJ/mol.A^2]
	reference [A]
	"""
	ai = 3 * a_i
	aj = 3 * a_j
	dr = [ i-j for i,j in zip( molec.coor[ai:ai+3], molec.coor[aj:aj+3] ) ]
	r2 = sum( [ i * i for i in dr ] )
	vv = math.sqrt( r2 )
	df = kumb * ( vv - xref )
	if( vv >= skip_LE and vv <= skip_BE ):
		molec.func += df * ( vv - xref ) * ffac
		if( grad ):
			df *= 2.0 / vv
			for i in [0, 1, 2]:
				molec.grad[ai+i] += df * dr[i] * gfac[0]
				molec.grad[aj+i] -= df * dr[i] * gfac[1]
		if( hess ):
			tt  = ( 2.0 * kumb - df ) / r2
			hxx = ( tt * dr[0] * dr[0] + df )
			hxy =   tt * dr[0] * dr[1]
			hxz =   tt * dr[0] * dr[2]
			hyy = ( tt * dr[1] * dr[1] + df )
			hyz =   tt * dr[1] * dr[2]
			hzz = ( tt * dr[2] * dr[2] + df )
			# ii & jj -- hessian should have been previously initialized...
			n  = int( math.sqrt( len( molec.hess ) ) )
			n2 = n * 2
			for ii in hfac:
				if( ii > -1 ):
					jj = 3 * ( n * ii + ii )
					molec.hess[jj]      += hxx
					molec.hess[jj+1]    += hxy
					molec.hess[jj+2]    += hxz
					molec.hess[jj+n]    += hxy
					molec.hess[jj+n+1]  += hyy
					molec.hess[jj+n+2]  += hyz
					molec.hess[jj+n2]   += hxz
					molec.hess[jj+n2+1] += hyz
					molec.hess[jj+n2+2] += hzz
			# ij & ji (only if both atoms are involved...)
			if( hfac[0] > -1 and hfac[1] > -1 ):
				for ii,jj in [ ( hfac[0], hfac[1] ), ( hfac[1], hfac[0] ) ]:
					kk = 3 * ( n * ii + jj )
					molec.hess[kk]      -= hxx
					molec.hess[kk+1]    -= hxy
					molec.hess[kk+2]    -= hxz
					molec.hess[kk+n]    -= hxy
					molec.hess[kk+n+1]  -= hyy
					molec.hess[kk+n+2]  -= hyz
					molec.hess[kk+n2]   -= hxz
					molec.hess[kk+n2+1] -= hyz
					molec.hess[kk+n2+2] -= hzz
	return( vv )


def mm_angle( molec, kumb, xref, a_i, a_j, a_k,
			ffac = 1.0, grad = False, gfac = [ 1.0, 1.0, 1.0 ], hess = False, hfac = None ):
	"""
	angle = force_constant * ( angle - reference )^2

	force_constant [kJ/mol.rad^2]
	reference & return_value [rad]
	"""
	ai = 3 * a_i
	aj = 3 * a_j
	ak = 3 * a_k
	dij = [ i-j for i,j in zip( molec.coor[ai:ai+3], molec.coor[aj:aj+3] ) ]
	rij = math.sqrt( sum( [ i * i for i in dij ] ) )
	dij = [ i / rij for i in dij ]
	dkj = [ k-j for k,j in zip( molec.coor[ak:ak+3], molec.coor[aj:aj+3] ) ]
	rkj = math.sqrt( sum( [ i * i for i in dkj ] ) )
	dkj = [ i / rkj for i in dkj ]
	dot = sum( [ i * j for i,j in zip( dij, dkj ) ] )
	dot = min( 1.0, max( -1.0, dot ) )
	vv = math.acos( dot )
	df = kumb * ( vv - xref )
	molec.func += df * ( vv - xref ) * ffac
	if( grad ):
		df *= - 2.0 / math.sqrt( 1.0 - dot * dot )
		dti = [ ( dkj[i] - dot * dij[i] ) / rij for i in [ 0, 1, 2 ] ]
		dtk = [ ( dij[i] - dot * dkj[i] ) / rkj for i in [ 0, 1, 2 ] ]
		dtj = [ - ( dti[i] + dtk[i] ) for i in [ 0, 1, 2 ] ]
		for i in [0, 1, 2]:
			molec.grad[ai+i] += df * dti[i] * gfac[0]
			molec.grad[aj+i] += df * dtj[i] * gfac[1]
			molec.grad[ak+i] += df * dtk[i] * gfac[2]
	if( hess ):
		pass
	return( vv )


def mm_dihedral( molec, data, a_i, a_j, a_k, a_l,
				ffac = 1.0, grad = False, gfac = [ 1.0, 1.0, 1.0, 1.0 ], hess = False, hfac = None ):
	"""
	dihedral = force_constant * ( 1 + cos( periodicity * angle - displacement ) )

	force_constant [kJ/mol]
	displacement [rad]

	data = [ frc_per=1, dsp_per=1, frc_per=2, dsp_per=2, ..., frc_per=6, dsp_per=6 ]
	"""
	ai  = 3 * a_i
	aj  = 3 * a_j
	ak  = 3 * a_k
	al  = 3 * a_l
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
	cs1 = sum( [ ii*jj for ii,jj in zip( vt, vu ) ] ) / rtu
	cs1 = min( 1.0, max( -1.0, cs1 ) )
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
	if( data[0] != 0.0 ):
		cd  = math.cos( data[1] )
		sd  = math.sin( data[1] )
		dph += data[0] * ( cs1 * sd - sn1 * cd )
		molec.func += data[0] * ( 1.0 + cs1 * cd + sn1 * sd ) * ffac
	if( data[2] != 0.0 ):
		cd  = math.cos( data[3] )
		sd  = math.sin( data[3] )
		dph += data[2] * 2.0 * ( cs2 * sd - sn2 * cd )
		molec.func += data[2] * ( 1.0 + cs2 * cd + sn2 * sd ) * ffac
	if( data[4] != 0.0 ):
		cd  = math.cos( data[5] )
		sd  = math.sin( data[5] )
		dph += data[4] * 3.0 * ( cs3 * sd - sn3 * cd )
		molec.func += data[4] * ( 1.0 + cs3 * cd + sn3 * sd ) * ffac
	if( data[6] != 0.0 ):
		cd  = math.cos( data[7] )
		sd  = math.sin( data[7] )
		dph += data[6] * 4.0 * ( cs4 * sd - sn4 * cd )
		molec.func += data[6] * ( 1.0 + cs4 * cd + sn4 * sd ) * ffac
	if( data[8] != 0.0 ):
		cd  = math.cos( data[9] )
		sd  = math.sin( data[9] )
		dph += data[8] * 5.0 * ( cs5 * sd - sn5 * cd )
		molec.func += data[8] * ( 1.0 + cs5 * cd + sn5 * sd ) * ffac
	if( data[10] != 0.0 ):
		cd  = math.cos( data[11] )
		sd  = math.sin( data[11] )
		dph += data[10] * 6.0 * ( cs6 * sd - sn6 * cd )
		molec.func += data[10] * ( 1.0 + cs6 * cd + sn6 * sd ) * ffac
	if( grad ):
		dki = [ ii-jj for ii,jj in zip( molec.coor[ak:ak+3], molec.coor[ai:ai+3] ) ]
		dlj = [ ii-jj for ii,jj in zip( molec.coor[al:al+3], molec.coor[aj:aj+3] ) ]
		dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
				( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
				( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
		dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
				( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
				( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
		molec.grad[ai]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph * gfac[0]
		molec.grad[ai+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph * gfac[0]
		molec.grad[ai+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph * gfac[0]
		molec.grad[aj]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph * gfac[1]
		molec.grad[aj+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph * gfac[1]
		molec.grad[aj+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph * gfac[1]
		molec.grad[ak]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph * gfac[2]
		molec.grad[ak+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph * gfac[2]
		molec.grad[ak+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph * gfac[2]
		molec.grad[al]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph * gfac[3]
		molec.grad[al+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph * gfac[3]
		molec.grad[al+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph * gfac[3]
	if( hess ):
		pass
	ang = qm3.constants.R2D * math.acos( cs1 )
	if( sn1 <= 0.0 ):
		ang = -ang
	return( ang )


def mm_improper( molec, kumb, xref, a_i, a_j, a_k, a_l,
				ffac = 1.0, grad = False, gfac = [ 1.0, 1.0, 1.0, 1.0 ], hess = False, hfac = None ):
	"""
	improper = force_constant * ( angle - reference )^2

	force_constant [kJ/mol.rad^2]
	reference [deg]
	a_i should be central atom
	"""
	ai  = 3 * a_i
	aj  = 3 * a_j
	ak  = 3 * a_k
	al  = 3 * a_l
	dji = [ ii-jj for ii,jj in zip( molec.coor[aj:aj+3], molec.coor[ai:ai+3] ) ]
	dkj = [ ii-jj for ii,jj in zip( molec.coor[ak:ak+3], molec.coor[aj:aj+3] ) ]
	dlk = [ ii-jj for ii,jj in zip( molec.coor[al:al+3], molec.coor[ak:ak+3] ) ]
	vt  = [ dji[1] * dkj[2] - dkj[1] * dji[2], dji[2] * dkj[0] - dkj[2] * dji[0], dji[0] * dkj[1] - dkj[0] * dji[1] ]
	vu  = [ dkj[1] * dlk[2] - dlk[1] * dkj[2], dkj[2] * dlk[0] - dlk[2] * dkj[0], dkj[0] * dlk[1] - dlk[0] * dkj[1] ]
	vtu = [ vt[1] * vu[2] - vu[1] * vt[2], vt[2] * vu[0] - vu[2] * vt[0], vt[0] * vu[1] - vu[0] * vt[1] ]
	rt2 = sum( [ ii*ii for ii in vt ] )
	ru2 = sum( [ ii*ii for ii in vu ] )
	rtu = math.sqrt( rt2 * ru2 )
	rkj = math.sqrt( sum( [ ii*ii for ii in dkj ] ) )
	cos = sum( [ ii*jj for ii,jj in zip( vt, vu ) ] ) / rtu
	sin = sum( [ ii*jj for ii,jj in zip( dkj, vtu ) ] ) / ( rkj * rtu )
	cos = min( 1.0, max( -1.0, cos ) )
	ang = qm3.constants.R2D * math.acos( cos )
	if( sin <= 0.0 ):
		ang = -ang
	if( math.fabs( ang + xref ) < math.fabs( ang - xref ) ):
		xref = -xref
	dt  = ang - xref
	while( dt >  180.0 ):
		dt -= 360.0
	while( dt < -180.0 ):
		dt += 360.0
	dt /= qm3.constants.R2D
	molec.func += kumb * dt * dt * ffac
	if( grad ):
		dph = 2.0 * kumb * dt
		dki = [ ii-jj for ii,jj in zip( molec.coor[ak:ak+3], molec.coor[ai:ai+3] ) ]
		dlj = [ ii-jj for ii,jj in zip( molec.coor[al:al+3], molec.coor[aj:aj+3] ) ]
		dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
				( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
				( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
		dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
				( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
				( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
		molec.grad[ai]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph * gfac[0]
		molec.grad[ai+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph * gfac[0]
		molec.grad[ai+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph * gfac[0]
		molec.grad[aj]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph * gfac[1]
		molec.grad[aj+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph * gfac[1]
		molec.grad[aj+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph * gfac[1]
		molec.grad[ak]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph * gfac[2]
		molec.grad[ak+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph * gfac[2]
		molec.grad[ak+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph * gfac[2]
		molec.grad[al]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph * gfac[3]
		molec.grad[al+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph * gfac[3]
		molec.grad[al+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph * gfac[3]
	if( hess ):
		pass
	return( ang )




class distance( object ):
	def __init__( self, kumb, xref, indx, skip_LE = 0.0, skip_BE = 9.e99 ):
		self.kumb = kumb
		self.xref = xref
		self.indx = indx[:]
		self.skpL = skip_LE
		self.skpB = skip_BE
		self.ffac = 1.0
		self.gfac = [ 1.0, 1.0 ]
		self.hfac = [ -1, -1 ]


	def get_func( self, molec ):
		return( mm_bond( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.skpL, self.skpB, self.ffac ) )


	def get_grad( self, molec ):
		return( mm_bond( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.skpL, self.skpB,
				self.ffac, True, self.gfac ) )


	def get_hess( self, molec ):
		return( mm_bond( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.skpL, self.skpB,
				self.ffac, True, self.gfac, True, self.hfac ) )



class angle( object ):
	def __init__( self, kumb, xref, indx ):
		self.kumb = kumb
		self.xref = xref / qm3.constants.R2D
		self.indx = indx[:]
		self.ffac = 1.0
		self.gfac = [ 1.0, 1.0, 1.0 ]
		self.hfac = [ -1, -1, -1 ]


	def get_func( self, molec ):
		return( mm_angle( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2], self.ffac ) * qm3.constants.R2D )


	def get_grad( self, molec ):
		return( mm_angle( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2],
				self.ffac, True, self.gfac ) * qm3.constants.R2D )


	def get_hess( self, molec ):
		return( mm_angle( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2],
				self.ffac, True, self.gfac, True, self.hfac ) * qm3.constants.R2D )



class dihedral( object ):
	def __init__( self, data, indx ):
		"""
	data = {  periodicity: [ force_constant [kJ/mol], displacement [degrees] ], ... }

	X - C_sp3 - C_sp3 - X   =>  { 1: [ -0.3347, 0.0 ], 2: [ -0.5858, 84.8 ], 3: [ 1.3263, 0.0 ] }

	valid periodicities = [ 1 : 6 ]
		"""
		self.data = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
		for i in range( 6 ):
			if( i+1 in data ):
				self.data[2*i]   = data[i+1][0]
				self.data[2*i+1] = data[i+1][1] / qm3.constants.R2D
		self.ffac = 1.0
		self.gfac = [ 1.0, 1.0, 1.0, 1.0 ]
		self.hfac = [ -1, -1, -1, -1 ]


	def get_func( self, molec ):
		return( mm_dihedral( molec, data, self.indx[0], self.indx[1], self.indx[2], self.indx[3], self.ffac ) )


	def get_grad( self, molec ):
		return( mm_dihedral( molec, data, self.indx[0], self.indx[1], self.indx[2], self.indx[3],
				self.ffac, True, self.gfac ) )


	def get_hess( self, molec ):
		return( mm_dihedral( molec, data, self.indx[0], self.indx[1], self.indx[2], self.indx[3],
				self.ffac, True, self.gfac, True, self.hfac ) )



class improper( object ):
	def __init__( self, kumb, xref, indx ):
		self.kumb = kumb
		self.xref = xref
		self.indx = indx[:]
		self.ffac = 1.0
		self.gfac = [ 1.0, 1.0, 1.0, 1.0 ]
		self.hfac = [ -1, -1, -1, -1 ]


	def get_func( self, molec ):
		return( mm_improper( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2], self.indx[3], self.ffac ) )


	def get_grad( self, molec ):
		return( mm_improper( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2], self.indx[3],
				self.ffac, True, self.gfac ) )


	def get_hess( self, molec ):
		return( mm_improper( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2], self.indx[3],
				self.ffac, True, self.gfac, True, self.hfac ) )



class multiple_distance( object ):
	def __init__( self, kumb, xref, indx, weigth ):
		"""
	multiple_distance = force_constant * ( value - reference )^2

	value = SUM weigth_i * distance_i

	force_constant [kJ/mol.A^2]
	reference [A]
		"""
		if( len( weigh ) * 2 != len( indx ) ):
			print( "- restraints.multiple_distance: Number of ATOMS should be TWICE the number of WEIGHTS!" )
			return( None )
		self.kumb = kumb
		self.xref = xref
		self.indx = indx[:]
		self.weig = weigth[:]
		self.size = len( weigth )


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
		molec.func += df * ( vv - self.xref )
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
		molec.func += df * ( vv - self.xref )
		for i in range( self.size ):
			i3 = i * 3
			ii = 3 * self.indx[2*i]
			jj = 3 * self.indx[2*i+1]
			tt = 2.0 * self.weig[i] * df / r[i]
			for j in [0, 1, 2]:
				molec.grad[ii+j] += tt * dr[i3+j]
				molec.grad[jj+j] -= tt * dr[i3+j]
		return( vv )



class tether( object ):
	def __init__( self, molec, kumb, indx ):
		"""
	thether = force_constant * SUM ( cartesian - reference )^2

	force_constant [kJ/mol.A^2]
	reference [A]
		"""
		self.kumb = kumb
		self.indx = {}
		for i in indx:
			self.indx[i] = molec.coor[3*i:3*i+3][:]


	def get_func( self, molec ):
		for w in self.indx:
			w3 = w * 3
			dr = [ i-j for i,j in zip( molec.coor[w3:w3+3], self.indx[w] ) ]
			molec.func += self.kumb * sum( [ i * i for i in dr ] )
		return( None )


	def get_grad( self, molec ):
		for w in self.indx:
			w3 = w * 3
			dr = [ i-j for i,j in zip( molec.coor[w3:w3+3], self.indx[w] ) ]
			molec.func += self.kumb * sum( [ i * i for i in dr ] ) 
			for i in [0, 1, 2]:
				molec.grad[w3+i] += 2.0 * self.kumb * dr[i]
		return( None )


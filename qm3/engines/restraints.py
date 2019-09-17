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




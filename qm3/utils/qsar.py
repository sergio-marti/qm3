# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.elements
import	qm3.maths.matrix
import	qm3.utils



def remove_non_polar_H( mol, bnd = [] ):
	if( bnd == [] ):
		bnd = qm3.utils.connectivity( mol )
	con = [ [] for i in range( mol.natm ) ]
	for i,j in bnd:
		con[i].append( j )
		con[j].append( i )
	out = qm3.mol.molecule()
	rem = []
	chg = []
	for i in range( mol.natm ):
		if( mol.anum[i] != 1 or mol.anum[con[i][0]] != 6 ):
			out.labl.append( mol.labl[i] )
			out.resi.append( mol.resi[i] )
			out.resn.append( mol.resn[i] )
			out.segn.append( mol.segn[i] )
			out.coor += mol.coor[3*i:3*i+3]
			out.anum.append( mol.anum[i] )
			out.chrg.append( mol.chrg[i] )
			if( mol.type ):
				out.type.append( mol.type[i] )
			if( mol.mass ):
				out.mass.append( mol.mass[i] )
			if( mol.epsi ):
				out.epsi.append( mol.epsi[i] )
			if( mol.rmin ):
				out.rmin.append( mol.rmin[i] )
			out.natm += 1
			rem.append( i )
		else:
			chg.append( ( con[i][0], mol.chrg[i] ) )
	for w,q in chg:
		out.chrg[rem.index( w )] += q
	out.settle()
	return( out )



class topological_index( object ):
	def __init__( self, mol, bond = [] ):
		self.natm = mol.natm
		self.anum = mol.anum
		if( bond ):
			self.bond = bond[:]
		else:
			self.bond = qm3.utils.connectivity( mol )
		self.conn = [ [] for i in range( self.natm ) ]
		for i,j in self.bond:
			self.conn[i].append( j )
			self.conn[j].append( i )
		self.QNA_U = None
		self.QNA_V = None


	def Randic( self ):
		t = [ 1.0 / math.sqrt( float( len( self.conn[i] ) ) ) for i in range( self.natm ) ]
		return( sum( [ t[i] * t[j] for i,j in self.bond ] ) )


	@staticmethod
	def __hosoya( bnd, hal = []  ):
		if( hal == [] ):
			hal = bnd[:]
		lst = []
		k = 2 + len( hal[0] )
		for i in range( len( hal ) ):
			for j in range( len( bnd ) ):
				t = hal[i] + bnd[j]
				if( len( { m:None for m in t } ) == k ):
					t.sort()
					lst.append( t[:] )
		if( len( lst ) > 0 ):
			out = [ lst[0] ]
			for i in range( 1, len( lst ) ):
				if( not lst[i] in out ):
					out.append( lst[i][:] )
		else:
			out = []
		return( out )
	
	def Hosoya( self ):
		t = self.__hosoya( self.bond )
		s = 1 + len( self.bond ) + len( t )
		while( t != [] ):
			t = self.__hosoya( self.bond, t )
			s += len( t )
		return( s )


	@staticmethod
	def __shortest_path( con, ii, jj ):
		tmp = [ [ ii ] ]
		siz = len( tmp )
		cur = 0
		while( cur < siz ):
			tt = [ i for i in con[tmp[cur][-1]] if not i in tmp[cur] ]
			if( len( tt ) == 0 ):
				cur += 1
			elif( jj in tt ):
				for i in range( len( tt ) ):
					if( tt[i] != jj ):
						tmp.append( tmp[cur][:] + [ tt[i] ] )
						siz += 1
				tmp[cur].append( jj )
				cur += 1
			else:
				for i in range( 1, len( tt ) ):
					tmp.append( tmp[cur][:] + [ tt[i] ] )
					siz += 1
				tmp[cur].append( tt[0] )
		s = -1
		for i in range( len( tmp ) ):
			if( ( tmp[i][0] == ii and tmp[i][-1] == jj ) and ( s < 0 or len( tmp[i] ) < s ) ):
				s = len( tmp[i] )
		return( s - 1 )

	def Wiener( self ):
		o = []
		for i in range( self.natm - 1 ):
			for j in range( i+1, self.natm ):
				o.append( self.__shortest_path( self.conn, i, j ) )
		return( sum( o ) )


	# Chebyshev polynomials Tuv / Filimonov et al [10.1080/10629360903438370]
	def __calc_QNA( self ):
		con = []
		for i in range( self.natm ):
			for j in range( self.natm ):
				if( i == j ):
					con.append( 0.0 )
				elif( j in self.conn[i] ):
					con.append( 1.0 )
				else:
					con.append( 0.0 )
		val, vec = qm3.maths.matrix.diag( con, self.natm )
		val = qm3.maths.matrix.from_diagonal( [ math.exp( -0.5 * val[i] ) for i in range( self.natm ) ], self.natm )
		exp = qm3.maths.matrix.mult( qm3.maths.matrix.mult( vec, self.natm, self.natm, val, self.natm, self.natm ), self.natm, self.natm, qm3.maths.matrix.T( vec, self.natm, self.natm ), self.natm, self.natm )
		A = []
		B = []
		for i in range( self.natm ):
			A.append( 0.5 * ( qm3.elements.ionpot[self.anum[i]] + qm3.elements.eafin[self.anum[i]] ) )
			B.append( 1.0 / math.sqrt( qm3.elements.ionpot[self.anum[i]] - qm3.elements.eafin[self.anum[i]] ) )
		P = []
		Q = []
		for i in range( self.natm ):
			P.append( B[i] * sum( [ B[j] * exp[i*self.natm+j] for j in range( self.natm ) ] ) )
			Q.append( B[i] * sum( [ A[j] * B[j] * exp[i*self.natm+j] for j in range( self.natm ) ] ) )
		mP = sum( P ) / float( self.natm )
		sP = math.sqrt( sum( [ math.pow( P[i] - mP, 2.0 ) for i in range( self.natm ) ] ) / float( self.natm - 1 ) )
		mQ = sum( Q ) / float( self.natm )
		sQ = math.sqrt( sum( [ math.pow( Q[i] - mQ, 2.0 ) for i in range( self.natm ) ] ) / float( self.natm - 1 ) )
		PQ = sum( [ ( P[i] - mP ) * ( Q[i] - mQ ) for i in range( self.natm ) ] ) / ( float( self.natm - 1 ) * sP * sQ )
		P  = [ ( P[i] - mP ) / sP for i in range( self.natm ) ]
		Q  = [ ( Q[i] - mQ ) / sQ for i in range( self.natm ) ]
		tp = 1.0 / math.sqrt( 2.0 * ( 1.0 + PQ ) )
		tm = 1.0 / math.sqrt( 2.0 * ( 1.0 - PQ ) )
		self.QNA_U  = [ ( P[i] + Q[i] ) * tp for i in range( self.natm ) ]
		self.QNA_V  = [ ( P[i] - Q[i] ) * tm for i in range( self.natm ) ]

	def QNA( self, u, v ):
		if( not self.QNA_U or not self.QNA_V ):
			self.__calc_QNA()
		return( sum( [ math.cos( u * math.acos( math.tanh( self.QNA_U[i] ) ) ) * math.cos( v * math.acos( math.tanh( self.QNA_V[i] ) ) ) for i in range( self.natm ) ] ) / float( self.natm ) )

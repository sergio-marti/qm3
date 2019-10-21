# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	os
import	inspect
import	numpy


class sff( object ):
	def __init__( self ):
		self.__initialize()
		self.__K2J = 4.184
		self.__R2D = 180.0 / math.pi
		self.__RCV = { 0: 0.31, 1: 0.31, 2: 0.28, 3: 1.28, 4: 0.96, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57, 10: 0.58,
			11: 1.66, 12: 1.41, 13: 1.21, 14: 1.11, 15: 1.07, 16: 1.05, 17: 1.02, 18: 1.06, 19: 2.03, 20: 1.76,
			21: 1.70, 22: 1.60, 23: 1.53, 24: 1.39, 25: 1.39, 26: 1.32, 27: 1.26, 28: 1.24, 29: 1.32, 30: 1.22,
			31: 1.22, 32: 1.20, 33: 1.19, 34: 1.20, 35: 1.20, 36: 1.16, 37: 2.20, 38: 1.95, 39: 1.90, 40: 1.75,
			41: 1.64, 42: 1.54, 43: 1.47, 44: 1.46, 45: 1.42, 46: 1.39, 47: 1.45, 48: 1.44, 49: 1.42, 50: 1.39,
			51: 1.39, 52: 1.38, 53: 1.39, 54: 1.40, 55: 2.44, 56: 2.15, 57: 2.07, 58: 2.04, 59: 2.03, 60: 2.01,
			61: 1.99, 62: 1.98, 63: 1.98, 64: 1.96, 65: 1.94, 66: 1.92, 67: 1.92, 68: 1.89, 69: 1.90, 70: 1.87,
			71: 1.87, 72: 1.75, 73: 1.70, 74: 1.62, 75: 1.51, 76: 1.44, 77: 1.41, 78: 1.36, 79: 1.36, 80: 1.32,
			81: 1.45, 82: 1.46, 83: 1.48, 84: 1.40, 85: 1.50, 86: 1.50, 87: 2.60, 88: 2.21, 89: 2.15, 90: 2.06,
			91: 2.00, 92: 1.96, 93: 1.90, 94: 1.87, 95: 1.80, 96: 1.69, 97: 1.60, 98: 1.60, 99: 1.60, 100: 1.60,
			101: 1.60, 102: 1.60, 103: 1.60, 104: 1.60, 105: 1.60, 106: 1.60, 107: 1.60, 108: 1.60, 109: 1.60 }
		self.__SMB = { 0 : "X", 1 : "H", 2 : "He", 3 : "Li", 4 : "Be", 5 : "B", 6 : "C", 7 : "N", 8 : "O", 9 : "F", 10 : "Ne",
			11 : "Na", 12 : "Mg", 13 : "Al", 14 : "Si", 15 : "P", 16 : "S", 17 : "Cl", 18 : "Ar", 19 : "K", 20 : "Ca",
			21 : "Sc", 22 : "Ti", 23 : "V", 24 : "Cr", 25 : "Mn", 26 : "Fe", 27 : "Co", 28 : "Ni", 29 : "Cu", 30 : "Zn",
			31 : "Ga", 32 : "Ge", 33 : "As", 34 : "Se", 35 : "Br", 36 : "Kr", 37 : "Rb", 38 : "Sr", 39 : "Y", 40 : "Zr",
			41 : "Nb", 42 : "Mo", 43 : "Tc", 44 : "Ru", 45 : "Rh", 46 : "Pd", 47 : "Ag", 48 : "Cd", 49 : "In", 50 : "Sn",
			51 : "Sb", 52 : "Te", 53 : "I", 54 : "Xe", 55 : "Cs", 56 : "Ba", 57 : "La", 58 : "Ce", 59 : "Pr", 60 : "Nd",
			61 : "Pm", 62 : "Sm", 63 : "Eu", 64 : "Gg", 65 : "Tb", 66 : "Dy", 67 : "Ho", 68 : "Er", 69 : "Tm", 70 : "Yb",
			71 : "Lu", 72 : "Hf", 73 : "Ta", 74 : "W", 75 : "Re", 76 : "Os", 77 : "Ir", 78 : "Pt", 79 : "Au", 80 : "Hg",
			81 : "Tl", 82 : "Pb", 83 : "Bi", 84 : "Po", 85 : "At", 86 : "Rn", 87 : "Fr", 88 : "Ra", 89 : "Ac", 90 : "Th",
			91 : "Pa", 92 : "U", 93 : "Np", 94 : "Pu", 95 : "Am", 96 : "Cm", 97 : "Bk", 98 : "Cf", 99 : "Es", 100 : "Fm",
			101 : "Md", 102 : "No", 103 : "Lr", 104 : "Rf", 105 : "Db", 106 : "Sg", 107 : "Bh", 108 : "Hs", 109 : "Mt" }


	def __initialize( self ):
		# ----------------- Molecule
		self.natm = 0
		self.coor = []
		self.chrg = []
		self.type = []
		self.anum = []
		self.epsi = []
		self.rmin = []
		self.free = []
		# ----------------- MM
		self.func = 0.0
		self.grad = []
		self.bond = []
		self.conn = []
		self.angl = []
		self.dihe = []
		self.impr = []
		self.bond_data = []
		self.bond_indx = []
		self.angl_data = []
		self.angl_indx = []
		self.dihe_data = []
		self.dihe_indx = []
		self.nbnd = []


	def load_top( self, fname ):
		self.__initialize()
		f = open( fname, "rt" )
		for l in f:
			t = l.strip().split()
			self.type.append( t[0] )
			self.anum.append( int( t[1] ) )
			self.chrg.append( float( t[2] ) )
		f.close()
		self.natm = len( self.type )


	def load_xyz( self, fname ):
		f = open( fname, "rt" )
		n = int( f.readline().strip() )
		if( n != self.natm ):
			print( "- wrong number of atoms... %d (top) vs %d (xyz)"%( self.natm, n ) )
		else:
			self.coor = []
			f.readline()
			for i in range( self.natm ):
				t = f.readline().strip().split()
				self.coor += [ float( t[1] ), float( t[2] ), float( t[3] ) ]
		f.close()


	def write_xyz( self, fname, append = False ):
		if( append ):
			f = open( fname, "at" )
		else:
			f = open( fname, "wt" )
		f.write( "%d\n\n"%( self.natm ) )
		for i in range( self.natm ):
			f.write( "%-2s%20.10lf%20.10lf%20.10lf\n"%( self.__SMB[self.anum[i]], self.coor[3*i], self.coor[3*i+1], self.coor[3*i+2] ) )
		f.close()


	def parse_molecule( self, guess_bonds = True ):
		if( guess_bonds ):
			self.bond = []
			for i in range( self.natm - 1 ):
				i3 = i * 3
				ri = self.__RCV[self.anum[i]] + 0.05
				for j in range( i + 1, self.natm ):
					if( self.anum[i] == 1 and self.anum[j] == 1 ):
						continue
					rj = self.__RCV[self.anum[j]] + 0.05
					t  = ( ri + rj ) * ( ri + rj )
					if( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip( self.coor[i3:i3+3], self.coor[j*3:j*3+3] ) ] ) <= t ):
						self.bond.append( [ i, j ] )
		self.conn = [ [] for i in range( self.natm ) ]
		for i,j in self.bond:
			self.conn[i].append( j )
			self.conn[j].append( i )
		self.angl = []
		for i in range( len( self.bond ) - 1 ):
			for j in range( i + 1, len( self.bond ) ):
				if( self.bond[i][0] == self.bond[j][0] ):
					self.angl.append( [ self.bond[i][1], self.bond[i][0], self.bond[j][1] ] )
				elif( self.bond[i][0] == self.bond[j][1] ):
					self.angl.append( [ self.bond[i][1], self.bond[i][0], self.bond[j][0] ] )
				elif( self.bond[i][1] == self.bond[j][0] ):
					self.angl.append( [ self.bond[i][0], self.bond[i][1], self.bond[j][1] ] )
				elif( self.bond[i][1] == self.bond[j][1] ):
					self.angl.append( [ self.bond[i][0], self.bond[i][1], self.bond[j][0] ] )
		for i in range( len( self.angl )-1, -1, -1 ):
			if( self.angl[i][0] in self.conn[self.angl[i][2]] ):
				del self.angl[i]
		self.dihe = []
		for i in range( len( self.angl ) - 1 ):
			for j in range( i + 1, len( self.angl ) ):
				if( self.angl[i][1] == self.angl[j][0] and self.angl[i][2] == self.angl[j][1] ):
					self.dihe.append( [ self.angl[i][0], self.angl[i][1], self.angl[i][2], self.angl[j][2] ] )
				elif( self.angl[i][1] == self.angl[j][2] and self.angl[i][2] == self.angl[j][1] ):
					self.dihe.append( [ self.angl[i][0], self.angl[i][1], self.angl[i][2], self.angl[j][0] ] )
				elif( self.angl[i][1] == self.angl[j][0] and self.angl[i][0] == self.angl[j][1] ):
					self.dihe.append( [ self.angl[i][2], self.angl[i][1], self.angl[i][0], self.angl[j][2] ] )
				elif( self.angl[i][1] == self.angl[j][2] and self.angl[i][0] == self.angl[j][1] ):
					self.dihe.append( [ self.angl[i][2], self.angl[i][1], self.angl[i][0], self.angl[j][0] ] )
		for i in range( len( self.dihe )-1, -1, -1 ):
			if( self.dihe[i][0] in self.conn[self.dihe[i][3]] ):
				del self.dihe[i]


	def partial_charges( self, params = None ):
		if( params != None ):
			f = open( ffield, "rt" )
		else:
			f = open( os.path.abspath( os.path.dirname( inspect.getfile( self.__class__ ) ) ) + os.sep + "sff.eem", "rt" )
		kap = float( f.readline().strip() )
		prm = {}
		for l in f:
			t = l.strip().split()
			prm[t[0]] = [ float( t[1] ), float( t[2] ) ]
		f.close()
		crd =  numpy.array( self.coor ).reshape( ( self.natm, 3 ) )
		mat = []
		vec = []
		for i in range( self.natm ):
			for j in range( self.natm ):
				if( j == i ):
					mat.append( prm[self.type[i]][1] )
				else:
					mat.append( kap / numpy.linalg.norm( crd[i] - crd[j] ) )
			mat.append( -1.0 )
			vec.append( - prm[self.type[i]][0] )
		mat += [ 1.0 ] * self.natm + [ 0.0 ]
		mat = numpy.array( mat ).reshape( ( self.natm + 1, self.natm + 1 ) )
		vec.append( sum( self.chrg ) )
		vec = numpy.array( vec ).reshape( ( self.natm + 1, 1 ) )
		chg = numpy.linalg.solve( mat, vec )
		self.chrg = chg[0:self.natm].flatten().tolist()


	def load_parameters( self, ffield = None ):
		out = True
		self.bond_data = []
		self.bond_indx = []
		self.angl_data = []
		self.angl_indx = []
		self.dihe_data = []
		self.dihe_indx = []
		if( ffield != None ):
			f = open( ffield, "rt" )
		else:
			f = open( os.path.abspath( os.path.dirname( inspect.getfile( self.__class__ ) ) ) + os.sep + "sff.prm", "rt" )
		tmp_typ = {}
		cnt_bnd = 0
		tmp_bnd = {}
		cnt_ang = 0
		tmp_ang = {}
		cnt_dih = 0
		tmp_dih = {}
		for l in f:
			t = l.strip().split()
			if( len( t ) > 0 and t[0][0] != "#" ):
				if( len( t ) == 3 ):
					tmp_typ[t[0]] = [ math.sqrt( float( t[1] ) * self.__K2J ), float( t[2] ) ]
				elif( len( t ) == 4 ):
					self.bond_data.append( [ float( t[2] ) * self.__K2J, float( t[3] ) ] )
					tmp_bnd["%s:%s"%( t[0], t[1] )] = cnt_bnd
					cnt_bnd += 1
				elif( len( t ) == 5 ):
					self.angl_data.append( [ float( t[3] ) * self.__K2J, float( t[4] ) / self.__R2D ] )
					tmp_ang["%s:%s:%s"%( t[0], t[1], t[2] )] = cnt_ang
					cnt_ang += 1
				elif( len( t ) >= 7 ):
					tmp = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
					for i in range( 4, len( t ), 3 ):
						n = int( t[i+1] ) - 1
						if( n >= 0 and n < 6 ):
							tmp[2*n]   = float( t[i]   ) * self.__K2J
							tmp[2*n+1] = float( t[i+2] ) / self.__R2D
					self.dihe_data.append( tmp[:] )
					tmp_dih["%s:%s:%s:%s"%( t[0], t[1], t[2], t[3] )] = cnt_dih
					cnt_dih += 1
		f.close()
		self.epsi = []
		self.rmin = []
		for i in range( self.natm ):
			if( self.type[i] in tmp_typ ):
				self.epsi.append( tmp_typ[self.type[i]][0] )
				self.rmin.append( tmp_typ[self.type[i]][1] )
			else:
				self.epsi.append( None )
				self.rmin.append( None )
				print( "- missing atom type [%s]: %d"%( self.type[i], i+1 ) )
				out = False
		for i,j in self.bond:
			td = "%s:%s"%( self.type[i], self.type[j] )
			ti = "%s:%s"%( self.type[j], self.type[i] )
			if( td in tmp_bnd ):
				self.bond_indx.append( tmp_bnd[td] )
			elif( ti in tmp_bnd ):
				self.bond_indx.append( tmp_bnd[ti] )
			else:
				self.bond_indx.append( None )
				print( "- missing parameter [bond]: ", td )
				out = False
		for i,j,k in self.angl:
			td = "%s:%s:%s"%( self.type[i], self.type[j], self.type[k] )
			ti = "%s:%s:%s"%( self.type[k], self.type[j], self.type[i] )
			ts = "*:%s:*"%( self.type[j] )
			if( td in tmp_ang ):
				self.angl_indx.append( tmp_ang[td] )
			elif( ti in tmp_ang ):
				self.angl_indx.append( tmp_ang[ti] )
			elif( ts in tmp_ang ):
				self.angl_indx.append( tmp_ang[ts] )
			else:
				self.angl_indx.append( None )
				print( "- missing parameter [angl]: ", td )
				out = False
		for i,j,k,l in self.dihe:
			td = "%s:%s:%s:%s"%( self.type[i], self.type[j], self.type[k], self.type[l] )
			ti = "%s:%s:%s:%s"%( self.type[l], self.type[k], self.type[j], self.type[i] )
			ts = "*:%s:%s:*"%( self.type[j], self.type[k] )
			tz = "*:%s:%s:*"%( self.type[k], self.type[j] )
			if( td in tmp_dih ):
				self.dihe_indx.append( tmp_dih[td] )
			elif( ti in tmp_dih ):
				self.dihe_indx.append( tmp_dih[ti] )
			elif( ts in tmp_dih ):
				self.dihe_indx.append( tmp_dih[ts] )
			elif( tz in tmp_dih ):
				self.dihe_indx.append( tmp_dih[tz] )
			else:
				self.dihe_indx.append( None )
				print( "- missing parameter [dihe]: ", td )
				out = False
		return( out )


	def __ebond( self, gradient ):
		if( self.bond == [] ):
			return( 0.0 )
		out = 0.0
		for i in range( len( self.bond ) ):
			ai  = 3 * self.bond[i][0]
			aj  = 3 * self.bond[i][1]
			vec = [ ii-jj for ii,jj in zip( self.coor[ai:ai+3], self.coor[aj:aj+3] ) ]
			val = math.sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] )
			dif = val - self.bond_data[self.bond_indx[i]][1]
			tmp = dif * self.bond_data[self.bond_indx[i]][0]
			out += tmp * dif * ( self.free[self.bond[i][0]] or self.free[self.bond[i][1]] )
			if( gradient ):
				tmp *= 2.0 / val
				for j in [0, 1, 2]:
					self.grad[ai+j] += tmp * vec[j] * self.free[self.bond[i][0]]
					self.grad[aj+j] -= tmp * vec[j] * self.free[self.bond[i][1]]
		return( out )


	def __eangle( self, gradient ):
		if( self.angl == [] ):
			return( 0.0 )
		out = 0.0
		for i in range( len( self.angl ) ):
			ai  = 3 * self.angl[i][0]
			aj  = 3 * self.angl[i][1]
			ak  = 3 * self.angl[i][2]
			dij = [ ii-jj for ii,jj in zip( self.coor[ai:ai+3], self.coor[aj:aj+3] ) ]
			rij = math.sqrt( sum( [ j*j for j in dij ] ) )
			dij = [ j / rij for j in dij ]
			dkj = [ ii-jj for ii,jj in zip( self.coor[ak:ak+3], self.coor[aj:aj+3] ) ]
			rkj = math.sqrt( sum( [ j*j for j in dkj ] ) )
			dkj = [ j / rkj for j in dkj ]
			fac = sum( [ ii*jj for ii,jj in zip( dij, dkj ) ] )
			fac = min( math.fabs( fac ), 1.0 - 1.0e-6 ) * fac / math.fabs( fac )
			val = math.acos( fac )
			dif = val - self.angl_data[self.angl_indx[i]][1]
			tmp = dif * self.angl_data[self.angl_indx[i]][0]
			out += tmp * dif * ( self.free[self.angl[i][0]] or self.free[self.angl[i][1]] or self.free[self.angl[i][2]] )
			if( gradient ):
				dtx =  - 1.0 / math.sqrt( 1.0 - fac * fac )
				tmp *= 2.0 / dtx
				dti = [ ( ii - fac * jj ) / rij for ii,jj in zip( dkj, dij ) ]
				dtk = [ ( ii - fac * jj ) / rkj for ii,jj in zip( dij, dkj ) ]
				dtj = [ - ( ii + jj ) for ii,jj in zip( dti, dtk ) ]
				for j in [0, 1, 2]:
					self.grad[ai+j] += tmp * dti[j] * self.free[self.angl[i][0]]
					self.grad[aj+j] += tmp * dtj[j] * self.free[self.angl[i][1]]
					self.grad[ak+j] += tmp * dtk[j] * self.free[self.angl[i][2]]
		return( out )


	# Generic dihedral ("adapated" from Tinker)
	def __edihedral( self, gradient ):
		if( self.dihe == [] ):
			return( 0.0 )
		out = 0.0
		for i in range( len( self.dihe ) ):
			ai  = 3 * self.dihe[i][0]
			aj  = 3 * self.dihe[i][1]
			ak  = 3 * self.dihe[i][2]
			al  = 3 * self.dihe[i][3]
			dji = [ ii-jj for ii,jj in zip( self.coor[aj:aj+3], self.coor[ai:ai+3] ) ]
			dkj = [ ii-jj for ii,jj in zip( self.coor[ak:ak+3], self.coor[aj:aj+3] ) ]
			rkj = math.sqrt( sum( [ ii*ii for ii in dkj ] ) )
			dlk = [ ii-jj for ii,jj in zip( self.coor[al:al+3], self.coor[ak:ak+3] ) ]
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
			FF  = ( self.free[self.dihe[i][0]] or self.free[self.dihe[i][1]] or self.free[self.dihe[i][2]] or self.free[self.dihe[i][3]] )
			if( self.dihe_data[self.dihe_indx[i]][0] != 0.0 ):
				cd  = math.cos( self.dihe_data[self.dihe_indx[i]][1] )
				sd  = math.sin( self.dihe_data[self.dihe_indx[i]][1] )
				dph += self.dihe_data[self.dihe_indx[i]][0] * ( cs1 * sd - sn1 * cd )
				out += self.dihe_data[self.dihe_indx[i]][0] * ( 1.0 + cs1 * cd + sn1 * sd ) * FF
			if( self.dihe_data[self.dihe_indx[i]][2] != 0.0 ):
				cd  = math.cos( self.dihe_data[self.dihe_indx[i]][3] )
				sd  = math.sin( self.dihe_data[self.dihe_indx[i]][3] )
				dph += self.dihe_data[self.dihe_indx[i]][2] * 2.0 * ( cs2 * sd - sn2 * cd )
				out += self.dihe_data[self.dihe_indx[i]][2] * ( 1.0 + cs2 * cd + sn2 * sd ) * FF
			if( self.dihe_data[self.dihe_indx[i]][4] != 0.0 ):
				cd  = math.cos( self.dihe_data[self.dihe_indx[i]][5] )
				sd  = math.sin( self.dihe_data[self.dihe_indx[i]][5] )
				dph += self.dihe_data[self.dihe_indx[i]][4] * 3.0 * ( cs3 * sd - sn3 * cd )
				out += self.dihe_data[self.dihe_indx[i]][4] * ( 1.0 + cs3 * cd + sn3 * sd ) * FF
			if( self.dihe_data[self.dihe_indx[i]][6] != 0.0 ):
				cd  = math.cos( self.dihe_data[self.dihe_indx[i]][7] )
				sd  = math.sin( self.dihe_data[self.dihe_indx[i]][7] )
				dph += self.dihe_data[self.dihe_indx[i]][6] * 4.0 * ( cs4 * sd - sn4 * cd )
				out += self.dihe_data[self.dihe_indx[i]][6] * ( 1.0 + cs4 * cd + sn4 * sd ) * FF
			if( self.dihe_data[self.dihe_indx[i]][8] != 0.0 ):
				cd  = math.cos( self.dihe_data[self.dihe_indx[i]][9] )
				sd  = math.sin( self.dihe_data[self.dihe_indx[i]][9] )
				dph += self.dihe_data[self.dihe_indx[i]][8] * 5.0 * ( cs5 * sd - sn5 * cd )
				out += self.dihe_data[self.dihe_indx[i]][8] * ( 1.0 + cs5 * cd + sn5 * sd ) * FF
			if( self.dihe_data[self.dihe_indx[i]][10] != 0.0 ):
				cd  = math.cos( self.dihe_data[self.dihe_indx[i]][11] )
				sd  = math.sin( self.dihe_data[self.dihe_indx[i]][11] )
				dph += self.dihe_data[self.dihe_indx[i]][10] * 6.0 * ( cs6 * sd - sn6 * cd )
				out += self.dihe_data[self.dihe_indx[i]][10] * ( 1.0 + cs6 * cd + sn6 * sd ) * FF
			if( gradient ):
				dki = [ ii-jj for ii,jj in zip( self.coor[ak:ak+3], self.coor[ai:ai+3] ) ]
				dlj = [ ii-jj for ii,jj in zip( self.coor[al:al+3], self.coor[aj:aj+3] ) ]
				dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
						( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
						( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
				dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
						( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
						( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
				self.grad[ai]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph * self.free[self.dihe[i][0]]
				self.grad[ai+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph * self.free[self.dihe[i][0]]
				self.grad[ai+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph * self.free[self.dihe[i][0]]
				self.grad[aj]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph * self.free[self.dihe[i][1]]
				self.grad[aj+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph * self.free[self.dihe[i][1]]
				self.grad[aj+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph * self.free[self.dihe[i][1]]
				self.grad[ak]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph * self.free[self.dihe[i][2]]
				self.grad[ak+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph * self.free[self.dihe[i][2]]
				self.grad[ak+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph * self.free[self.dihe[i][2]]
				self.grad[al]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph * self.free[self.dihe[i][3]]
				self.grad[al+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph * self.free[self.dihe[i][3]]
				self.grad[al+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph * self.free[self.dihe[i][3]]
		return( out )


	# Impropers (also "adapated" from Tinker)
	def __eimproper( self, gradient ):
		if( self.impr == [] ):
			return( 0.0 )
		out = 0.0
		# self.impr = [ [ central_i, j, k, l, kmb (kcal/mol), ref (deg) ], ... ]
		for i in range( len( self.impr ) ):
			ai  = 3 * self.impr[i][0]
			aj  = 3 * self.impr[i][1]
			ak  = 3 * self.impr[i][2]
			al  = 3 * self.impr[i][3]
			dji = [ ii-jj for ii,jj in zip( self.coor[aj:aj+3], self.coor[ai:ai+3] ) ]
			dkj = [ ii-jj for ii,jj in zip( self.coor[ak:ak+3], self.coor[aj:aj+3] ) ]
			dlk = [ ii-jj for ii,jj in zip( self.coor[al:al+3], self.coor[ak:ak+3] ) ]
			vt  = [ dji[1] * dkj[2] - dkj[1] * dji[2], dji[2] * dkj[0] - dkj[2] * dji[0], dji[0] * dkj[1] - dkj[0] * dji[1] ]
			vu  = [ dkj[1] * dlk[2] - dlk[1] * dkj[2], dkj[2] * dlk[0] - dlk[2] * dkj[0], dkj[0] * dlk[1] - dlk[0] * dkj[1] ]
			vtu = [ vt[1] * vu[2] - vu[1] * vt[2], vt[2] * vu[0] - vu[2] * vt[0], vt[0] * vu[1] - vu[0] * vt[1] ]
			rt2 = sum( [ ii*ii for ii in vt ] )
			ru2 = sum( [ ii*ii for ii in vu ] )
			rtu = math.sqrt( rt2 * ru2 )
			if( rtu == 0.0 ):
				continue
			rkj = math.sqrt( sum( [ ii*ii for ii in dkj ] ) )
			cos = sum( [ ii*jj for ii,jj in zip( vt, vu ) ] ) / rtu
			sin = sum( [ ii*jj for ii,jj in zip( dkj, vtu ) ] ) / ( rkj * rtu )
			cos = min( 1.0, max( -1.0, cos ) )
			ang = self.__R2D * math.acos( cos )
			if( sin <= 0.0 ):
				ang = -ang
			kmb = self.impr[i][4] * self.__K2J
			ref = self.impr[i][5]
			if( math.fabs( ang + ref ) < math.fabs( ang - ref ) ):
				ref = -ref
			dt  = ang - ref
			while( dt >  180.0 ):
				dt -= 360.0
			while( dt < -180.0 ):
				dt += 360.0
			dt  /= self.__R2D
			FI  = self.free[self.dihe[i][0]]
			FJ  = self.free[self.dihe[i][1]]
			FK  = self.free[self.dihe[i][2]]
			FL  = self.free[self.dihe[i][3]]
			out += kmb * dt * dt * ( self.free[self.dihe[i][0]] or self.free[self.dihe[i][1]] or self.free[self.dihe[i][2]] or self.free[self.dihe[i][3]] )
			if( gradient ):
				dph = 2.0 * kmb * dt
				dki = [ ii-jj for ii,jj in zip( self.coor[ak:ak+3], self.coor[ai:ai+3] ) ]
				dlj = [ ii-jj for ii,jj in zip( self.coor[al:al+3], self.coor[aj:aj+3] ) ]
				dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
						( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
						( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
				dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
						( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
						( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
				self.grad[ai]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph * self.free[self.dihe[i][0]]
				self.grad[ai+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph * self.free[self.dihe[i][0]]
				self.grad[ai+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph * self.free[self.dihe[i][0]]
				self.grad[aj]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph * self.free[self.dihe[i][1]]
				self.grad[aj+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph * self.free[self.dihe[i][1]]
				self.grad[aj+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph * self.free[self.dihe[i][1]]
				self.grad[ak]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph * self.free[self.dihe[i][2]]
				self.grad[ak+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph * self.free[self.dihe[i][2]]
				self.grad[ak+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph * self.free[self.dihe[i][2]]
				self.grad[al]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph * self.free[self.dihe[i][3]]
				self.grad[al+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph * self.free[self.dihe[i][3]]
				self.grad[al+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph * self.free[self.dihe[i][3]]
		return( out )


	def non_bonding( self ):
		self.nbnd = []
		for i in range( self.natm - 1 ):
			for j in range( i+1, self.natm ):
				f = False
				n = len( self.bond )
				k = 0
				while( k < n and not f ):
					f |= ( ( i == self.bond[k][0] and j == self.bond[k][1] ) or ( i == self.bond[k][1] and j == self.bond[k][0] )  )
					k += 1
				n = len( self.angl )
				k = 0
				while( k < n and not f ):
					f |= ( ( i == self.angl[k][0] and j == self.angl[k][2] ) or ( i == self.angl[k][2] and j == self.angl[k][0] )  )
					k += 1
				n = len( self.dihe )
				k = 0
				while( k < n and not f ):
					f |= ( ( i == self.dihe[k][0] and j == self.dihe[k][3] ) or ( i == self.dihe[k][3] and j == self.dihe[k][0] )  )
					k += 1
				if( not f ):
					self.nbnd.append( [ i, j, 1.0 ] )
		for i,j,k,l in self.dihe:
			self.nbnd.append( [ i, l, 0.5 ] )


	def __enbonded( self, gradient, epsilon = 1.0 ):
		epsf = 1389.35484620709144110151 / epsilon
		if( self.nbnd == [] ):
			return( 0.0, 0.0 )
		oel = 0.0
		olj = 0.0
		for i,j,sc in self.nbnd:
			ai  = 3 * i
			aj  = 3 * j
			dr  = [ ii-jj for ii,jj in zip( self.coor[ai:ai+3], self.coor[aj:aj+3] ) ]
			r2  = sum( [ ii*ii for ii in dr ] )
			eij = self.epsi[i] * self.epsi[j]
			sij = self.rmin[i] + self.rmin[j]
			qij = self.chrg[i] * self.chrg[j] * epsf
			r   = 1.0 / math.sqrt( r2 )
			s6  = math.pow( sij * r, 6.0 )
			tmp = qij * r
			FF  = ( self.free[i] or self.free[j] ) * sc
			oel += tmp * FF
			olj += eij * s6 * ( s6 - 2.0 ) * FF
			if( gradient ):
				df = sc * ( 12.0 * eij * s6 * ( 1.0 - s6 ) - tmp ) / r2;
				for k in [0, 1, 2]:
					self.grad[ai+k] += df * dr[k] * self.free[i]
					self.grad[aj+k] -= df * dr[k] * self.free[j]
		return( oel, olj )


	def get_func( self, qprint = False ):
		e_bond = self.__ebond( False )
		e_angl = self.__eangle( False )
		e_dihe = self.__edihedral( False )
		e_impr = self.__eimproper( False )
		e_elec, e_vdwl = self.__enbonded( False )
		self.func = e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl
		if( qprint ):
			print( "ETot:", e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl, "_kJ/mol" )
			print( "   Bond:%18.4lf   Angl:%18.4lf   Dihe:%18.4lf"%( e_bond, e_angl, e_dihe ) )
			print( "   Impr:%18.4lf   Elec:%18.4lf   VdWl:%18.4lf"%( e_impr, e_elec, e_vdwl ) )


	def get_grad( self, qprint = False ):
		self.grad = [ 0.0 for i in range( 3 * self.natm ) ]
		e_bond = self.__ebond( True )
		e_angl = self.__eangle( True )
		e_dihe = self.__edihedral( True )
		e_impr = self.__eimproper( True )
		e_elec, e_vdwl = self.__enbonded( True )
		self.func = e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl
		if( qprint ):
			print( "ETot:", e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl, "_kJ/mol" )
			print( "   Bond:%18.4lf   Angl:%18.4lf   Dihe:%18.4lf"%( e_bond, e_angl, e_dihe ) )
			print( "   Impr:%18.4lf   Elec:%18.4lf   VdWl:%18.4lf"%( e_impr, e_elec, e_vdwl ) )


	def minimize( self, step_number = 1000, step_size = 0.1, print_frequency = 100, gradient_tolerance = 0.5 ):
		def __grms( vec, siz ):
			o = math.sqrt( sum( [ i*i for i in vec ] ) )
			return( o, o / math.sqrt( siz ) )
		print( "\n---------------------------------------- Minimization (FIRE)\n" )
		print( "Degrees of Freedom: %20ld"%( 3 * sum( self.free ) ) )
		print( "Step Number:        %20d"%( step_number ) )
		print( "Step Size:          %20.10lg"%( step_size ) )
		print( "Print Frequency:    %20d"%( print_frequency ) )
		print( "Gradient Tolerance: %20.10lg\n"%( gradient_tolerance ) )
		print( "%10s%20s%20s%20s"%( "Step", "Function", "Gradient", "Displacement" ) )
		print( "-" * 70 )
		nstp = 0
		ssiz = step_size
		alph = 0.1
		size = 3 * self.natm
		velo = [ 0.0 for i in range( size ) ]
		step = [ 0.0 for i in range( size ) ]
		self.get_grad()
		norm, grms = __grms( self.grad, 3.0 * sum( self.free ) )
		print( "%10s%20.5lf%20.8lf%20.10lf"%( "", self.func, grms, ssiz ) )
		i = 0
		while( i < step_number and grms > gradient_tolerance ):
			vsiz = math.sqrt( sum( [ velo[j] * velo[j] for j in range( size ) ] ) )
			vfac = sum( [ - velo[j] * self.grad[j] for j in range( size ) ] )
			if( vfac > 0.0 ):
				velo = [ ( 1.0 - alph ) * velo[j] - alph * self.grad[j] / norm * vsiz for j in range( size ) ]
				if( nstp > 5 ):
					ssiz = min( ssiz * 1.1, step_size )
					alph *= 0.99
				nstp += 1
			else:
				velo = [ 0.0 for j in range( size ) ]
				alph = 0.1
				ssiz *= 0.5
				nstp = 0
			for j in range( size ):
				velo[j] -= ssiz * self.grad[j]
				step[j] = ssiz * velo[j]
			tmp = math.sqrt( sum( [ step[j] * step[j] for j in range( size ) ] ) )
			if( tmp > ssiz ):
				for j in range( size ):
					step[j] *= ssiz / tmp
			for j in range( size ):
				self.coor[j] += step[j]
			self.get_grad()
			norm, grms = __grms( self.grad, 3.0 * sum( self.free ) )
			i = i + 1
			if( i%print_frequency == 0 ):
				print( "%10d%20.5lf%20.10lf%20.10lf"%( i, self.func, grms, ssiz ) )
		if( i%print_frequency != 0 ):
			print( "%10d%20.5lf%20.10lf%20.10lf"%( i + 1, self.func, grms, ssiz ) )
		print( "-" * 70 + "\n" )


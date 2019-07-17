# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	os
import	inspect
import	qm3.elements
import	qm3.maths.matrix
import	qm3.io
import	qm3.mol
import	qm3.utils
import	qm3.constants

try:
	import qm3.engines._mol_mech
	mol_mech_so = True
except:
	mol_mech_so = False



###################################################################################################
# SFF: Simple Force Field
#
class simple_force_field( object ):
	def __init__( self ):
		if( mol_mech_so ):
			self.ncpu = os.sysconf( 'SC_NPROCESSORS_ONLN' )
		else:
			self.ncpu = 1
		self.cut_on   = 10.0
		self.cut_off  = 12.0
		self.cut_list = 14.0
		self.initialize()


	def initialize( self ):
		self.natm = 0
		self.bond = []
		self.conn = []
		self.angl = []
		self.dihe = []
		self.impr = []
		self.nbnd = []
		self.nb14 = []
		self.bond_data = []
		self.bond_indx = []
		self.angl_data = []
		self.angl_indx = []
		self.dihe_data = []
		self.dihe_indx = []
		self.QNA_U = []
		self.QNA_V = []

	def guess_bonds( self, mol, quick = False ):
		self.initialize()
		self.natm = mol.natm
		self.bond = qm3.utils.connectivity( mol, quick )


	def calc_connectivity( self ):
		self.conn = [ [] for i in range( self.natm ) ]
		for i,j in self.bond:
			self.conn[i].append( j )
			self.conn[j].append( i )

	
	# SYBYL atom types (kinda)
	# http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
	# uses FORMAL CHARGES present in MOL.CHRG
	def guess_types( self, mol ):
		def __any( lst_a, lst_b ):
			return( len( set( lst_a ).intersection( set( lst_b ) ) ) > 0 )

		mol.type = []
		# default to atomic symbol...
		for i in range( mol.natm ):
			mol.type.append( qm3.elements.symbol[mol.anum[i]] )
			nb = len( self.conn[i] )
			if( mol.anum[i] == 1 ):
				if( mol.anum[self.conn[i][0]] == 6 ):
					mol.type[i] = "Hn"
			elif( mol.anum[i] == 14 ):
				mol.type[i] = "Si.3"
			elif( mol.anum[i] == 15 ):
				mol.type[i] = "P.3"
			elif( mol.anum[i] in [ 6, 7, 8, 16 ] ):
				mol.type[i] += "_%d"%( len( self.conn[i] ) )
		# 2nd pass
		for i in range( mol.natm ):
			if( mol.type[i] in [ "C_2", "C_1" ] ):
				mol.type[i] = "C.1"
			elif( mol.type[i] == "C_4" ):
				mol.type[i] = "C.3"

			elif( mol.type[i] == "C_3" ):
				if( __any( [ "C_3", "C.ar" ], [ mol.type[j] for j in self.conn[i] ] ) ):
					mol.type[i] = "C.ar"
				else:
					if( __any( [ "O_1", "O.co2", "O.2" ], [ mol.type[j] for j in self.conn[i] ] ) ):
						mol.type[i] = "C.co"
					else:
						mol.type[i] = "C.2"
			elif( mol.type[i] == "O_1" ):
				if( mol.type[self.conn[i][0]] in [ "C_3", "C.co" ] and mol.chrg[i] == -0.5 ):
					mol.type[i] = "O.co2"
				elif( mol.chrg[i] == -1.0 ):
					mol.type[i] = "O.x"
				else:
					mol.type[i] = "O.2"
			elif( mol.type[i] == "O_2" ):
				if( 1 in [ mol.anum[j] for j in self.conn[i] ] ):
					mol.type[i] = "O.h"
				else:
					mol.type[i] = "O.3"
			elif( mol.type[i] == "N_4" ):
				mol.type[i] = "N.4"
			elif( mol.type[i] == "N_3" ):
				if( mol.chrg[i] == 1.0 ):
					mol.type[i] = "N.pl"
				else:
					mol.type[i] = "N.3"
			elif( mol.type[i] == "N_2" ):
				mol.type[i] = "N.2"
			elif( mol.type[i] == "N_1" ):
				mol.type[i] = "N.1"
			elif( mol.type[i] == "S_2" ):
				if( 1 in [ mol.anum[j] for j in self.conn[i] ] ):
					mol.type[i] = "S.h"
				else:
					mol.type[i] = "S.3"
			elif( mol.type[i] == "S_1" ):
				if( mol.chrg[i] == -1.0 ):
					mol.type[i] = "S.x"
				else:
					mol.type[i] = "S.2"
			elif( mol.type[i] == "S_3" and _any( [ "O_1", "O.2" ], [ mol.type[j] for j in self.conn[i] ] ) ):
				mol.type[i] = "S.o"
			elif( mol.type[i] == "S_4" and _any( [ "O_1", "O.2" ], [ mol.type[j] for j in self.conn[i] ] ) ):
				mol.type[i] = "S.o2"


	def guess_angles( self ):
		if( mol_mech_so ):
			self.angl = qm3.engines._mol_mech.guess_angles( self )
		else:
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


	def guess_dihedrals( self ):
		if( mol_mech_so ):
			self.dihe = qm3.engines._mol_mech.guess_dihedrals( self )
		else:
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


	def parse_molecule( self, mol ):
		self.guess_bonds( mol )
		self.calc_connectivity()
		self.guess_angles()
		self.guess_dihedrals()
		self.guess_types( mol )


	# uses FORMAL CHARGES present in MOL.CHRG
	def guess_partial_charges( self, mol, method = "eem" ):
		if( method == "gasteiger" ):
			# Gasteiger partial charges ("adapted" from AmberTools/antechamber/charge.c)
			gas = {
				"H":     [  7.17,  6.24, -0.56,  20.02 ],	# polar hydrogen
				"Hn":    [  7.17,  6.24, -0.56,  20.02 ],	# non-polar hydrogen
				"C.1":   [ 10.39,  9.45,  0.73,  20.57 ],	# C sp1
				"C.2":   [  8.79,  9.32,  1.51,  19.62 ],	# C sp2 in single
				"C.ar":  [  8.79,  9.32,  1.51,  19.62 ],	# C sp2 in aromatic/conjugated
				"C.co":  [  8.79,  9.32,  1.51,  19.62 ],	# C sp2 in C=O
				"C.3":   [  7.98,  9.18,  1.88,  19.04 ],	# C sp3
				"N.1":   [ 15.68, 11.70, -0.27,  27.11 ],	# N sp1
				"N.2":   [ 12.87, 11.15,  0.85,  24.87 ],	# N sp2	in C=N
				"N.pl":  [ 12.32, 11.20,  1.34,  24.86 ],	# N sp2 (+)
				"N.3":   [ 11.54, 10.82,  1.36,  23.72 ],	# N sp3
				"N.4":   [  0.00, 11.86, 11.86,  23.72 ],	# N sp3 (+)
				"O.2":   [ 17.07, 13.79,  0.47,  31.33 ],	# O sp2 in C=O
				"O.3":   [ 14.18, 12.92,  1.39,  28.49 ],	# O sp3
				"O.h":   [ 17.07, 13.79,  0.47,  31.33 ],	# O sp3 in O-H
				"O.x":   [ 17.07, 13.79,  0.47,  31.33 ],	# O sp3 in O(-)
				"O.co2": [ 17.07, 13.79,  0.47,  31.33 ],	# O sp2/sp3 in CO2(-0.5 * 2)
				"S.2":   [ 10.88,  9.485, 1.325, 21.69 ],	# S sp2	
				"S.3":   [ 10.14,  9.13,  1.38,  20.65 ],	# S sp3
				"S.h":   [ 10.88,  9.485, 1.325, 21.69 ],	# S sp3 in S-H
				"S.x":   [ 10.88,  9.485, 1.325, 21.69 ],	# S sp3 in S(-)
				"S.o":   [ 10.14,  9.13,  1.38,  20.65 ],	# S sp2d in S=O
				"S.o2":  [ 12.00, 10.805, 1.195, 24.00 ],	# S spd2 in O=S=O
				"F":     [ 14.66, 13.85,  2.31,  30.82 ],
				"Cl":    [ 11.00,  9.69,  1.35,  22.04 ],
				"Br":    [ 10.08,  8.47,  1.16,  19.71 ],
				"I":     [  9.90,  7.96,  0.96,  18.82 ],
				"P.3":   [  8.90,  8.24,  0.96,  18.10 ]	# P (any)
			}
			ep = mol.chrg[:]
			ea = mol.chrg[:]
			ff = True
			it = 0
			df = 0.5
			while( it < 1000 and ff ):
				x = []
				for i in range( mol.natm ):
					x.append( gas[mol.type[i]][0] + ep[i] * ( gas[mol.type[i]][1] + ep[i] * gas[mol.type[i]][2] ) )
					if( x[i] == 0.0 ):
						x[i] = 1.0e-10
				for i,j in self.bond:
					if( x[i] <= x[j] ):
						q = ( x[j] - x[i] ) / gas[mol.type[i]][3] * df
						ea[i] += q
						ea[j] -= q
					else:
						q = ( x[i] - x[j] ) / gas[mol.type[j]][3] * df
						ea[i] -= q
						ea[j] += q
				s = 0.0
				for i in range( mol.natm ):
					s += ( ep[i] - ea[i] ) * ( ep[i] - ea[i] )
					ep[i] = ea[i]
				ff = math.sqrt( s / float( mol.natm ) ) > 1.0e-5
				df *= 0.5
				it += 1
			mol.chrg = ea[:]
		else:
			# Electronegativity Equalization Method (B3LYP_6-311G_NPA.par) [10.1186/s13321-015-0107-1]
			kap = 0.2509
			prm = {
				"H":     [ 2.3864, 0.6581 ],	# polar hydrogen
				"Hn":    [ 2.3864, 0.6581 ],	# non-polar hydrogen
				"C.1":   [ 2.4617, 0.3489 ],	# C sp1
				"C.2":   [ 2.5065, 0.3173 ],	# C sp2 in single
				"C.ar":  [ 2.5065, 0.3173 ],	# C sp2 in aromatic/conjugated
				"C.co":  [ 2.5065, 0.3173 ],	# C sp2 in C=O
				"C.3":   [ 2.4992, 0.3220 ],	# C sp3
				"N.1":   [ 2.5348, 0.4025 ],	# N sp1
				"N.2":   [ 2.5568, 0.2949 ],	# N sp2	in C=N
				"N.pl":  [ 2.5568, 0.2949 ],	# N sp2 (+)
				"N.3":   [ 2.5891, 0.4072 ],	# N sp3
				"N.4":   [ 2.5891, 0.4072 ],	# N sp3 (+)
				"O.2":   [ 2.6588, 0.4232 ],	# O sp2 in C=O
				"O.3":   [ 2.6342, 0.4041 ],	# O sp3
				"O.h":   [ 2.6342, 0.4041 ],	# O sp3 in O-H
				"O.x":   [ 2.6342, 0.4041 ],	# O sp3 in O(-)
				"O.co2": [ 2.6588, 0.4232 ],	# O sp2/sp3 in CO2(-0.5 * 2)
				"S.2":   [ 2.4884, 0.2043 ],	# S sp2	
				"S.3":   [ 2.4506, 0.2404 ],	# S sp3
				"S.h":   [ 2.4506, 0.2404 ],	# S sp3 in S-H
				"S.x":   [ 2.4506, 0.2404 ],	# S sp3 in S(-)
				"S.o":   [ 2.4884, 0.2043 ],	# S sp2d in S=O
				"S.o2":  [ 2.4884, 0.2043 ],	# S spd2 in O=S=O
				"F":     [ 3.0028, 1.2433 ],
				"Cl":    [ 2.5104, 0.8364 ],
				"Br":    [ 2.4244, 0.7511 ],
				"I":     [ 2.3272, 0.9303 ],
				"P.2":   [ 2.2098, 0.3281 ],	# P sp2
				"P.3":   [ 2.3898, 0.1902 ] 	# P sp3
			}
			mat = []
			vec = []
			for i in range( mol.natm ):
				for j in range( mol.natm ):
					if( j == i ):
						mat.append( prm[mol.type[i]][1] )
					else:
						mat.append( kap / qm3.utils.distance( mol.coor[3*i:3*i+3], mol.coor[3*j:3*j+3] ) )
				mat.append( -1 )
				vec.append( - prm[mol.type[i]][0] )
			mat += [ 1 ] * mol.natm + [ 0 ]
			vec.append( sum( mol.chrg ) )
			mol.chrg = qm3.maths.matrix.solve( mat, vec )[0:mol.natm]


	def load_parameters( self, mol, ffield = None ):
		out = True
		self.bond_data = []
		self.bond_indx = []
		self.angl_data = []
		self.angl_indx = []
		self.dihe_data = []
		self.dihe_indx = []
		if( ffield ):
			f = qm3.io.open_r( ffield )
		else:
			f = open( os.path.abspath( os.path.dirname( inspect.getfile( self.__class__ ) ) ) + os.sep + "mol_mech.prm", "rt" )
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
					tmp_typ[t[0]] = [ math.sqrt( float( t[1] ) * qm3.constants.K2J ), float( t[2] ) ]
				elif( len( t ) == 4 ):
					self.bond_data.append( [ float( t[2] ) * qm3.constants.K2J, float( t[3] ) ] )
					tmp_bnd["%s:%s"%( t[0], t[1] )] = cnt_bnd
					cnt_bnd += 1
				elif( len( t ) == 5 ):
					self.angl_data.append( [ float( t[3] ) * qm3.constants.K2J, float( t[4] ) / qm3.constants.R2D ] )
					tmp_ang["%s:%s:%s"%( t[0], t[1], t[2] )] = cnt_ang
					cnt_ang += 1
				elif( len( t ) >= 7 ):
					tmp = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
					for i in range( 4, len( t ), 3 ):
						n = int( t[i+1] ) - 1
						if( n >= 0 and n < 6 ):
							tmp[2*n]   = float( t[i]   ) * qm3.constants.K2J
							tmp[2*n+1] = float( t[i+2] ) / qm3.constants.R2D
					self.dihe_data.append( tmp[:] )
					tmp_dih["%s:%s:%s:%s"%( t[0], t[1], t[2], t[3] )] = cnt_dih
					cnt_dih += 1
		qm3.io.close( f, ffield )
		mol.epsi = []
		mol.rmin = []
		for i in range( mol.natm ):
			if( mol.type[i] in tmp_typ ):
				mol.epsi.append( tmp_typ[mol.type[i]][0] )
				mol.rmin.append( tmp_typ[mol.type[i]][1] )
			else:
				mol.epsi.append( None )
				mol.rmin.append( None )
				print( "- missing atom type [%s]: %d"%( mol.type[i], i+1 ) )
				out = False
		for i,j in self.bond:
			td = "%s:%s"%( mol.type[i], mol.type[j] )
			ti = "%s:%s"%( mol.type[j], mol.type[i] )
			if( td in tmp_bnd ):
				self.bond_indx.append( tmp_bnd[td] )
			elif( ti in tmp_bnd ):
				self.bond_indx.append( tmp_bnd[ti] )
			else:
				self.bond_indx.append( None )
				print( "- missing parameter [bond]: ", td )
				out = False
		for i,j,k in self.angl:
			td = "%s:%s:%s"%( mol.type[i], mol.type[j], mol.type[k] )
			ti = "%s:%s:%s"%( mol.type[k], mol.type[j], mol.type[i] )
			ts = "*:%s:*"%( mol.type[j] )
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
			td = "%s:%s:%s:%s"%( mol.type[i], mol.type[j], mol.type[k], mol.type[l] )
			ti = "%s:%s:%s:%s"%( mol.type[l], mol.type[k], mol.type[j], mol.type[i] )
			ts = "*:%s:%s:*"%( mol.type[j], mol.type[k] )
			tz = "*:%s:%s:*"%( mol.type[k], mol.type[j] )
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


	def energy_bond( self, mol, gradient = False ):
		if( self.bond == [] ):
			return( 0.0 )
		if( mol_mech_so ):
			out = qm3.engines._mol_mech.energy_bond( self, mol, gradient )
		else:
			out = 0.0
			for i in range( len( self.bond ) ):
				ai  = 3 * self.bond[i][0]
				aj  = 3 * self.bond[i][1]
				vec = [ ii-jj for ii,jj in zip( mol.coor[ai:ai+3], mol.coor[aj:aj+3] ) ]
				val = math.sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] )
				dif = val - self.bond_data[self.bond_indx[i]][1]
				tmp = dif * self.bond_data[self.bond_indx[i]][0]
				out += tmp * dif
				if( gradient ):
					tmp *= 2.0 / val
					for j in [0, 1, 2]:
						mol.grad[ai+j] += tmp * vec[j]
						mol.grad[aj+j] -= tmp * vec[j]
		return( out )


	def energy_angle( self, mol, gradient = False ):
		if( self.angl == [] ):
			return( 0.0 )
		if( mol_mech_so ):
			out = qm3.engines._mol_mech.energy_angle( self, mol, gradient )
		else:
			out = 0.0
			for i in range( len( self.angl ) ):
				ai  = 3 * self.angl[i][0]
				aj  = 3 * self.angl[i][1]
				ak  = 3 * self.angl[i][2]
				dij = [ ii-jj for ii,jj in zip( mol.coor[ai:ai+3], mol.coor[aj:aj+3] ) ]
				rij = math.sqrt( sum( [ j*j for j in dij ] ) )
				dij = [ j / rij for j in dij ]
				dkj = [ ii-jj for ii,jj in zip( mol.coor[ak:ak+3], mol.coor[aj:aj+3] ) ]
				rkj = math.sqrt( sum( [ j*j for j in dkj ] ) )
				dkj = [ j / rkj for j in dkj ]
				fac = sum( [ ii*jj for ii,jj in zip( dij, dkj ) ] )
				fac = min( math.fabs( fac ), 1.0 - 1.0e-6 ) * fac / math.fabs( fac )
				val = math.acos( fac )
				dif = val - self.angl_data[self.angl_indx[i]][1]
				tmp = dif * self.angl_data[self.angl_indx[i]][0]
				out += tmp * dif
				if( gradient ):
					dtx =  - 1.0 / math.sqrt( 1.0 - fac * fac )
					tmp *= 2.0 / dtx
					dti = [ ( ii - fac * jj ) / rij for ii,jj in zip( dkj, dij ) ]
					dtk = [ ( ii - fac * jj ) / rkj for ii,jj in zip( dij, dkj ) ]
					dtj = [ - ( ii + jj ) for ii,jj in zip( dti, dtk ) ]
					for j in [0, 1, 2]:
						mol.grad[ai+j] += tmp * dti[j]
						mol.grad[aj+j] += tmp * dtj[j]
						mol.grad[ak+j] += tmp * dtk[j]
		return( out )


	# Generic dihedral ("adapated" from Tinker)
	def energy_dihedral( self, mol, gradient = False ):
		if( self.dihe == [] ):
			return( 0.0 )
		if( mol_mech_so ):
			out = qm3.engines._mol_mech.energy_dihedral( self, mol, gradient )
		else:
			out = 0.0
			for i in range( len( self.dihe ) ):
				ai  = 3 * self.dihe[i][0]
				aj  = 3 * self.dihe[i][1]
				ak  = 3 * self.dihe[i][2]
				al  = 3 * self.dihe[i][3]
				dji = [ ii-jj for ii,jj in zip( mol.coor[aj:aj+3], mol.coor[ai:ai+3] ) ]
				dkj = [ ii-jj for ii,jj in zip( mol.coor[ak:ak+3], mol.coor[aj:aj+3] ) ]
				rkj = math.sqrt( sum( [ ii*ii for ii in dkj ] ) )
				dlk = [ ii-jj for ii,jj in zip( mol.coor[al:al+3], mol.coor[ak:ak+3] ) ]
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
				if( self.dihe_data[self.dihe_indx[i]][0] > 0.0 ):
					cd   = math.cos( self.dihe_data[self.dihe_indx[i]][1] )
					sd   = math.sin( self.dihe_data[self.dihe_indx[i]][1] )
					dph += self.dihe_data[self.dihe_indx[i]][0] * ( cs1 * sd - sn1 * cd )
					out += self.dihe_data[self.dihe_indx[i]][0] * ( 1.0 + cs1 * cd + sn1 * sd )
				if( self.dihe_data[self.dihe_indx[i]][2] > 0.0 ):
					cd   = math.cos( self.dihe_data[self.dihe_indx[i]][3] )
					sd   = math.sin( self.dihe_data[self.dihe_indx[i]][3] )
					dph += self.dihe_data[self.dihe_indx[i]][2] * 2.0 * ( cs2 * sd - sn2 * cd )
					out += self.dihe_data[self.dihe_indx[i]][2] * ( 1.0 + cs2 * cd + sn2 * sd )
				if( self.dihe_data[self.dihe_indx[i]][4] > 0.0 ):
					cd  = math.cos( self.dihe_data[self.dihe_indx[i]][5] )
					sd  = math.sin( self.dihe_data[self.dihe_indx[i]][5] )
					dph += self.dihe_data[self.dihe_indx[i]][4] * 3.0 * ( cs3 * sd - sn3 * cd )
					out += self.dihe_data[self.dihe_indx[i]][4] * ( 1.0 + cs3 * cd + sn3 * sd )
				if( self.dihe_data[self.dihe_indx[i]][6] > 0.0 ):
					cd  = math.cos( self.dihe_data[self.dihe_indx[i]][7] )
					sd  = math.sin( self.dihe_data[self.dihe_indx[i]][7] )
					dph = self.dihe_data[self.dihe_indx[i]][6] * 4.0 * ( cs4 * sd - sn4 * cd )
					out += self.dihe_data[self.dihe_indx[i]][6] * ( 1.0 + cs4 * cd + sn4 * sd )
				if( self.dihe_data[self.dihe_indx[i]][8] > 0.0 ):
					cd  = math.cos( self.dihe_data[self.dihe_indx[i]][9] )
					sd  = math.sin( self.dihe_data[self.dihe_indx[i]][9] )
					dph += self.dihe_data[self.dihe_indx[i]][8] * 5.0 * ( cs5 * sd - sn5 * cd )
					out += self.dihe_data[self.dihe_indx[i]][8] * ( 1.0 + cs5 * cd + sn5 * sd )
				if( self.dihe_data[self.dihe_indx[i]][10] > 0.0 ):
					cd  = math.cos( self.dihe_data[self.dihe_indx[i]][11] )
					sd  = math.sin( self.dihe_data[self.dihe_indx[i]][11] )
					dph += self.dihe_data[self.dihe_indx[i]][10] * 6.0 * ( cs6 * sd - sn6 * cd )
					out += self.dihe_data[self.dihe_indx[i]][10] * ( 1.0 + cs6 * cd + sn6 * sd )
				if( gradient ):
					dki = [ ii-jj for ii,jj in zip( mol.coor[ak:ak+3], mol.coor[ai:ai+3] ) ]
					dlj = [ ii-jj for ii,jj in zip( mol.coor[al:al+3], mol.coor[aj:aj+3] ) ]
					dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
							( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
							( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
					dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
							( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
							( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
					mol.grad[ai]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph
					mol.grad[ai+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph
					mol.grad[ai+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph
					mol.grad[aj]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph
					mol.grad[aj+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph
					mol.grad[aj+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph
					mol.grad[ak]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph
					mol.grad[ak+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph
					mol.grad[ak+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph
					mol.grad[al]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph
					mol.grad[al+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph
					mol.grad[al+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph
		return( out )


	# Impropers (also "adapated" from Tinker)
	def energy_improper( self, mol, gradient = False ):
		if( self.impr == [] ):
			return( 0.0 )
		out = 0.0
		# self.impr = [ [ central_i, j, k, l, kmb (kcal/mol), ref (deg) ], ... ]
		for i in range( len( self.impr ) ):
			ai  = 3 * self.impr[i][0]
			aj  = 3 * self.impr[i][1]
			ak  = 3 * self.impr[i][2]
			al  = 3 * self.impr[i][3]
			dji = [ ii-jj for ii,jj in zip( mol.coor[aj:aj+3], mol.coor[ai:ai+3] ) ]
			dkj = [ ii-jj for ii,jj in zip( mol.coor[ak:ak+3], mol.coor[aj:aj+3] ) ]
			dlk = [ ii-jj for ii,jj in zip( mol.coor[al:al+3], mol.coor[ak:ak+3] ) ]
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
			ang = qm3.constants.R2D * math.acos( cos )
			if( sin <= 0.0 ):
				ang = -ang
			kmb = self.impr[i][4] * qm3.constants.K2J
			ref = self.impr[i][5]
			if( math.fabs( ang + ref ) < math.fabs( ang - ref ) ):
				ref = -ref
			dt  = ang - ref
			while( dt >  180.0 ):
				dt -= 360.0
			while( dt < -180.0 ):
				dt += 360.0
			dt  /= qm3.constants.R2D
			out += kmb * dt * dt
			if( gradient ):
				dph = 2.0 * kmb * dt
				dki = [ ii-jj for ii,jj in zip( mol.coor[ak:ak+3], mol.coor[ai:ai+3] ) ]
				dlj = [ ii-jj for ii,jj in zip( mol.coor[al:al+3], mol.coor[aj:aj+3] ) ]
				dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
						( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
						( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
				dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
						( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
						( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
				mol.grad[ai]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph
				mol.grad[ai+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph
				mol.grad[ai+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph
				mol.grad[aj]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph
				mol.grad[aj+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph
				mol.grad[aj+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph
				mol.grad[ak]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph
				mol.grad[ak+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph
				mol.grad[ak+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph
				mol.grad[al]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph
				mol.grad[al+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph
				mol.grad[al+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph
		return( out )


	def update_non_bonded( self, mol ):
		if( mol_mech_so ):
			self.nbnd = qm3.engines._mol_mech.update_non_bonded( self, mol )
		else:
			self.nbnd = []
			if( self.cut_list > 0.0 ):
				c2l = self.cut_list * self.cut_list
			else:
				c2l = 1.0e99
			for i in range( mol.natm - 1 ):
				i3 = 3 * i
				i_crd = [ mol.coor[i3+k] - mol.boxl[k] * round( mol.coor[i3+k] / mol.boxl[k], 0 ) for k in [ 0, 1, 2 ] ]
				for j in range( i+1, mol.natm ):
					j_crd = [ mol.coor[j*3+k] - mol.boxl[k] * round( mol.coor[j*3+k] / mol.boxl[k], 0 ) for k in [ 0, 1, 2 ] ]
					if( qm3.utils.distanceSQ( i_crd, j_crd ) <= c2l ):
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
							self.nbnd.append( [ i, j ] )
		self.nb14 = []
		for i,j,k,l in self.dihe:
			self.nb14.append( [ i, l ] )
						

	def __non_bonded_interactions( self, mol, lst_nbnd, scale = 1.0, gradient = False, epsilon = 1.0 ):
		epsf = 1389.35484620709144110151 / epsilon
		oel  = 0.0
		olj  = 0.0
		if( self.cut_on > 0.0 and self.cut_off > self.cut_on ):
			c2on =  self.cut_on * self.cut_on
			c2of = self.cut_off * self.cut_off
			_g   = math.pow( c2of - c2on, 3.0 )
			_a   = c2of * c2of * ( c2of - 3.0 * c2on ) / _g
			_b   = 6.0 * c2of * c2on / _g
			_c   = - ( c2of + c2on ) / _g
			_d   = 0.4 / _g
			_el1 = 8.0 * ( c2of * c2on * ( self.cut_off - self.cut_on ) - 0.2 * ( self.cut_off * c2of * c2of - self.cut_on * c2on * c2on ) ) / _g
			_el2 = - _a / self.cut_off + _b * self.cut_off + _c * self.cut_off * c2of + _d * self.cut_off * c2of * c2of
			k6   = ( self.cut_off * c2of ) / ( self.cut_off * c2of - self.cut_on * c2on )
			k12  = math.pow( c2of, 3.0 ) / ( math.pow( c2of, 3.0 ) - math.pow( c2on, 3.0 ) )
			for i,j in lst_nbnd:
				ai = i * 3
				aj = j * 3
				dr = [ ii-jj for ii,jj in zip( mol.coor[ai:ai+3], mol.coor[aj:aj+3] ) ]
				dr = [ dr[k] - mol.boxl[k] * round( dr[k] / mol.boxl[k], 0 ) for k in [ 0, 1, 2] ]
				r2 = sum( [ ii*ii for ii in dr ] )
				if( r2 > c2of ):
					continue
				eij = mol.epsi[i] * mol.epsi[j]
				sij = mol.rmin[i] + mol.rmin[j]
				qij = mol.chrg[i] * mol.chrg[j] * epsf
				r   = math.sqrt( r2 )
				s   = 1.0 / r
				s3  = math.pow( sij * s, 3.0 )
				s6  = s3 * s3
				if( r2 <= c2on ):
					tmp = qij * s
					oel += tmp + qij * _el1
					s12  = s6 * s6
					_lj1 = math.pow( sij / self.cut_off * sij / self.cut_on, 3.0 )
					_lj2 = _lj1 * _lj1
					olj += eij * ( ( s12 - _lj2 ) - 2.0 * ( s6 - _lj1 ) )
					df   = ( 12.0 * eij * ( s6 - s12 ) - tmp ) / r2
				else:
					r3   = r * r2
					r5   = r3 * r2
					oel += qij * ( _a * s - _b * r - _c * r3 - _d * r5 + _el2 )
					_lj1 = math.pow( sij / self.cut_off, 3.0 )
					_lj2 = _lj1 * _lj1
					olj += eij * ( k12 * math.pow( s6 - _lj2, 2.0 ) - 2.0 * k6 * math.pow( s3 - _lj1, 2.0 ) )
					df   = - qij * ( _a / r3 + _b * s + 3.0 * _c * r + 5.0 * _d * r3 ) 
					df  -= 12.0 * eij * ( k12 * s6 * ( s6 - _lj2 ) - k6 * s3 * ( s3 - _lj1 ) ) / r2
				if( gradient ):
					for j in [0, 1, 2]:
						mol.grad[ai+j] += scale * df * dr[j]
						mol.grad[aj+j] -= scale * df * dr[j]
		else:
			for i,j in lst_nbnd:
				ai  = i * 3
				aj  = j * 3
				dr  = [ ii-jj for ii,jj in zip( mol.coor[ai:ai+3], mol.coor[aj:aj+3] ) ]
				dr = [ dr[k] - mol.boxl[k] * round( dr[k] / mol.boxl[k], 0 ) for k in [ 0, 1, 2] ]
				r2  = sum( [ ii*ii for ii in dr ] )
				eij = mol.epsi[i] * mol.epsi[j]
				sij = mol.rmin[i] + mol.rmin[j]
				qij = mol.chrg[i] * mol.chrg[j] * epsf
				r   = 1.0 / math.sqrt( r2 )
				s6  = math.pow( sij * r, 6.0 )
				tmp = qij * r
				oel += tmp
				olj += eij * s6 * ( s6 - 2.0 )
				if( gradient ):
					df = scale * ( 12.0 * eij * s6 * ( 1.0 - s6 ) - tmp ) / r2;
					for j in [0, 1, 2]:
						mol.grad[ai+j] += df * dr[j]
						mol.grad[aj+j] -= df * dr[j]
		oel *= scale
		olj *= scale
		return( oel, olj )

	def energy_non_bonded( self, mol, gradient = False, epsilon = 1.0 ):
		if( not self.nbnd ):
			self.update_non_bonded( mol )
		if( mol_mech_so ):
			oel, olj = qm3.engines._mol_mech.energy_non_bonded( self, mol, gradient )
		else:
			oel,   olj   = self.__non_bonded_interactions( mol, self.nbnd, 1.0, gradient, epsilon )
			oel14, olj14 = self.__non_bonded_interactions( mol, self.nb14, 0.5, gradient, epsilon )
			oel += oel14
			olj += olj14
		return( oel, olj )


	def get_func( self, mol, qprint = False ):
		e_bond = self.energy_bond( mol, gradient = False )
		e_angl = self.energy_angle( mol, gradient = False )
		e_dihe = self.energy_dihedral( mol, gradient = False )
		e_impr = self.energy_improper( mol, gradient = False )
		e_elec, e_vdwl = self.energy_non_bonded( mol, gradient = False )
		mol.func += e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl
		if( qprint ):
			print( "ETot:", e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl, "_kJ/mol" )
			print( "   Bond:%18.4lf   Angl:%18.4lf   Dihe:%18.4lf"%( e_bond, e_angl, e_dihe ) )
			print( "   Impr:%18.4lf   Elec:%18.4lf   VdWl:%18.4lf"%( e_impr, e_elec, e_vdwl ) )


	def get_grad( self, mol, qprint = False ):
		e_bond = self.energy_bond( mol, gradient = True )
		e_angl = self.energy_angle( mol, gradient = True )
		e_dihe = self.energy_dihedral( mol, gradient = True )
		e_impr = self.energy_improper( mol, gradient = True )
		e_elec, e_vdwl = self.energy_non_bonded( mol, gradient = True )
		mol.func += e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl
		if( qprint ):
			print( "ETot:", e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl, "_kJ/mol" )
			print( "   Bond:%18.4lf   Angl:%18.4lf   Dihe:%18.4lf"%( e_bond, e_angl, e_dihe ) )
			print( "   Impr:%18.4lf   Elec:%18.4lf   VdWl:%18.4lf"%( e_impr, e_elec, e_vdwl ) )



	###################################################################################################
	# QSAR topological stuff
	#
	def remove_non_polar_H( self, mol ):
		if( mol.anum and self.conn ):
			out = qm3.mol.molecule()
			rem = []
			chg = []
			for i in range( mol.natm ):
				if( mol.anum[i] != 1 or mol.anum[self.conn[i][0]] != 6 ):
					out.labl.append( mol.labl[i] )
					out.resi.append( mol.resi[i] )
					out.resn.append( mol.resn[i] )
					out.segn.append( mol.segn[i] )
					out.coor += mol.coor[3*i:3*i+3]
					out.anum.append( mol.anum[i] )
					if( mol.type ):
						out.type.append( mol.type[i] )
					if( mol.mass ):
						out.mass.append( mol.mass[i] )
					if( mol.chrg ):
						out.chrg.append( mol.chrg[i] )
					if( mol.epsi ):
						out.epsi.append( mol.epsi[i] )
					if( mol.rmin ):
						out.rmin.append( mol.rmin[i] )
					out.natm += 1
					rem.append( i )
				else:
					chg.append( ( self.conn[i][0], mol.chrg[i] ) )
			for w,q in chg:
				out.chrg[rem.index( w )] += q
			out.settle()
			self.initialize()
			self.guess_bonds( out )
			self.calc_connectivity()
			return( out )

	def topind_Randic( self ):
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
	
	def topind_Hosoya( self ):
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

	def topind_Wiener( self ):
		o = []
		for i in range( self.natm - 1 ):
			for j in range( i+1, self.natm ):
				o.append( self.__shortest_path( self.conn, i, j ) )
		return( sum( o ) )


	# Chebyshev polynomials Tuv / Filimonov et al [10.1080/10629360903438370]
	def __calc_QNA( self, mol ):
		con = []
		for i in range( mol.natm ):
			for j in range( mol.natm ):
				if( i == j ):
					con.append( 0.0 )
				elif( j in self.conn[i] ):
					con.append( 1.0 )
				else:
					con.append( 0.0 )
		val, vec = qm3.maths.matrix.diag( con, mol.natm )
		val = qm3.maths.matrix.from_diagonal( [ math.exp( -0.5 * val[i] ) for i in range( mol.natm ) ], mol.natm )
		exp = qm3.maths.matrix.mult( qm3.maths.matrix.mult( vec, mol.natm, mol.natm, val, mol.natm, mol.natm ), mol.natm, mol.natm, qm3.maths.matrix.T( vec, mol.natm, mol.natm ), mol.natm, mol.natm )
		A = []
		B = []
		for i in range( mol.natm ):
			A.append( 0.5 * ( qm3.elements.ionpot[mol.anum[i]] + qm3.elements.eafin[mol.anum[i]] ) )
			B.append( 1.0 / math.sqrt( qm3.elements.ionpot[mol.anum[i]] - qm3.elements.eafin[mol.anum[i]] ) )
		P = []
		Q = []
		for i in range( mol.natm ):
			P.append( B[i] * sum( [ B[j] * exp[i*mol.natm+j] for j in range( mol.natm ) ] ) )
			Q.append( B[i] * sum( [ A[j] * B[j] * exp[i*mol.natm+j] for j in range( mol.natm ) ] ) )
		mP = sum( P ) / float( mol.natm )
		sP = math.sqrt( sum( [ math.pow( P[i] - mP, 2.0 ) for i in range( mol.natm ) ] ) / float( mol.natm - 1 ) )
		mQ = sum( Q ) / float( mol.natm )
		sQ = math.sqrt( sum( [ math.pow( Q[i] - mQ, 2.0 ) for i in range( mol.natm ) ] ) / float( mol.natm - 1 ) )
		PQ = sum( [ ( P[i] - mP ) * ( Q[i] - mQ ) for i in range( mol.natm ) ] ) / ( float( mol.natm - 1 ) * sP * sQ )
		P  = [ ( P[i] - mP ) / sP for i in range( mol.natm ) ]
		Q  = [ ( Q[i] - mQ ) / sQ for i in range( mol.natm ) ]
		tp = 1.0 / math.sqrt( 2.0 * ( 1.0 + PQ ) )
		tm = 1.0 / math.sqrt( 2.0 * ( 1.0 - PQ ) )
		self.QNA_U  = [ ( P[i] + Q[i] ) * tp for i in range( mol.natm ) ]
		self.QNA_V  = [ ( P[i] - Q[i] ) * tm for i in range( mol.natm ) ]

	def topind_QNA( self, mol, u, v ):
		if( not self.QNA_U or not self.QNA_V ):
			self.__calc_QNA( mol )
		return( sum( [ math.cos( u * math.acos( math.tanh( self.QNA_U[i] ) ) ) * math.cos( v * math.acos( math.tanh( self.QNA_V[i] ) ) ) for i in range( mol.natm ) ] ) / float( mol.natm ) )

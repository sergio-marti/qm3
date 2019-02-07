#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

try:
	import cPickle as pickle
except:
	import pickle

try:
	import	tkinter
	import	tkinter.ttk as ttk
except:
	import	Tkinter as tkinter
	import	ttk

import	re
import	collections
import	qm3.mol


# ===============================================================================================


ENG =  collections.OrderedDict( [
# -- MM engines

#	[ "qm3.engines.dynamo.dynamo_pipe", [ [ "", "", "" ] ] ],
	
	[ "qm3.engines.dynamo.py_dynamo", [ [ "", "", "" ] ] ],

#	[ "qm3.engines.namd.namd", [
#		[ "psf", "PSF", "name of the PSF file used to define the model: str" ],
#		[ "cut", "10 12 14", "force-switching (floats): cut-on cut-off cut-list" ] ]
#	],

	[ "qm3.engines.namd.namd_pipe", [ 
		[ "psf", "PSF", "name of the PSF file used to define the model: str" ],
		[ "cut", "10 12 14", "force-switching (floats): cut-on cut-off cut-list" ] ]
	],

	[ "qm3.engines.namd.namd_shm", [ 
		[ "psf", "PSF", "name of the PSF file used to define the model: str" ],
		[ "cut", "10 12 14", "force-switching (floats): cut-on cut-off cut-list" ] ]
	],


	[ "qm3.engines.sander.sander", [
		[ "cut", "10", "cut-off: float" ],
		[ "qm_sel", ":1", "sander selection mask to set up the QM atoms: str" ],
		[ "qm_met", "AM1", "semi-empirical hamiltionan inplemented in sander: str" ],
		[ "qm_chg", "0", "charge of the QM atoms: int" ] ]
	],

	[ "qm3.engines.sander.py_sander", [
		[ "cut", "10", "cut-off: float" ],
		[ "qm_met", "AM1", "semi-empirical hamiltionan inplemented in sander: str" ],
		[ "qm_chg", "0", "charge of the QM atoms: int" ] ]
	],

	[ "qm3.engines.charmm.charmm_pipe", [ [ "", "", "" ] ] ],

	[ "qm3.engines.charmm.charmm_shm", [ [ "", "", "" ] ] ],

#	[ "qm3.engines._lammps.lammps", [ [ "", "", "" ] ] ],

	[ "qm3.engines._lammps.lammps_pipe", [ [ "", "", "" ] ] ],

	[ "qm3.engines._lammps.py_lammps", [ [ "", "", "" ] ] ],

# -- QM engines

	[ "qm3.engines.dftb.dftb", [
		[ "chg", "0", "charge of the QM atoms: int" ],
		[ "prm", "PRM", "folder of DFTB+ parameters set: str" ] ]
	],

	[ "qm3.engines.dftb.dl_dftb", [
		[ "chg", "0", "charge of the QM atoms: int" ],
		[ "prm", "PRM", "folder of DFTB+ parameters set: str" ] ]
	],

	[ "qm3.engines.sqm.sqm", [ [ "", "", "" ] ] ],

	[ "qm3.engines.gamess.gamess", [ [ "", "", "" ] ] ],

	[ "qm3.engines.nwchem.nwchem", [ [ "", "", "" ] ] ],

	[ "qm3.engines.demon.demon", [ [ "", "", "" ] ] ],

	[ "qm3.engines.orca.orca", [ [ "", "", "" ] ] ],
								 
#	[ "qm3.engines.smash.smash", [ [ "", "", "" ] ] ],
								 
	[ "qm3.engines.tchem.tchem", [ [ "", "", "" ] ] ],

	[ "qm3.engines.tchem.tchem_sckt", [ [ "", "", "" ] ] ],

	[ "qm3.engines.gaussian.gaussian", [ [ "", "", "" ] ] ],

	[ "qm3.engines.gaussian.gaussian_MMEL", [ [ "", "", "" ] ] ],

	[ "qm3.engines._psi4.py_psi4", [
		[ "opt", """{ "reference": "rks", "basis": "6-31g*", "d_convergence": 6, "scf_type": "direct", "guess": "read", "output": False, "charge": 0, "method": "b3lyp", "ncpus": 1, "memory": "1024 MB" }""", "python dictonary with Psi-4 options: { key: val }" ] ]
	],

# --- Internals / fixings...

	[ "qm3.engines._qmmm.Int_QMLJ", [
		[ "typ", "NBND", "name of the file with the Lennard-Jones types: str" ] ]
	],

	[ "qm3.engines._qmmm.Int_QMLJ_MMEL", [
		[ "typ", "NBND", "name of the file with the Lennard-Jones types: str" ] ]
	],

	[ "qm3.engines.restraints.distance", [ 
		[ "kmb", "", "force constant in kJ/mol.A^2: float" ],
		[ "ref", "", "reference value for distance in Angstroms: float" ],
		[ "idx", "", "C-indexes of the atoms defining the distance: list" ] ]
	],

	[ "qm3.engines.restraints.angle", [ 
		[ "kmb", "", "force constant in kJ/mol.deg^2: float" ],
		[ "ref", "", "reference value for angle in degrees: float" ],
		[ "idx", "", "C-indexes of the atoms defining the angle: list" ] ]
	],

	[ "qm3.engines.restraints.multiple_distance", [ 
		[ "kmb", "", "force constant in kJ/mol.A^2: float" ],
		[ "ref", "", "reference value for distance in Angstroms: float" ],
		[ "idx", "", "C-indexes of the atoms defining the distance: list" ],
		[ "wei", "", "weights for the linear combination of distances: list" ] ]
	] ] )


JOB = collections.OrderedDict( [
	[ "qm3.actions.minimize.steepest_descent", [
		[ "stp", "1000", "step number: int" ],
		[ "siz", "0.1", "step size: float" ],
		[ "prt", "100", "print frequency: int" ],
		[ "tol", "1.0", "gradient tolerance: float" ] ]
	],

	[ "qm3.actions.minimize.fire", [
		[ "stp", "1000", "step number: int" ],
		[ "siz", "0.1", "step size: float" ],
		[ "prt", "100", "print frequency: int" ],
		[ "tol", "1.0", "gradient tolerance: float" ] ]
	],

	[ "qm3.actions.minimize.l_bfgs", [
		[ "stp", "1000", "step number: int" ],
		[ "siz", "0.1", "step size: float" ],
		[ "prt", "100", "print frequency: int" ],
		[ "tol", "1.0", "gradient tolerance: float" ] ]
	],

	[ "qm3.actions.minimize.conjugate_gradient_plus", [
		[ "stp", "1000", "step number: int" ],
		[ "prt", "100", "print frequency: int" ],
		[ "tol", "1.0", "gradient tolerance: float" ] ]
	],

	[ "qm3.actions.dynamics.velocity_verlet", [
		[ "stp", "1000", "step number: int" ],
		[ "siz", "0.001", "step size in ps: float" ],
		[ "prt", "100", "print frequency: int" ],
		[ "scl", "100", """[NVT] temperature scaling frequency in steps: int (>0)""" ],
		[ "cpl", "0.1", """[NVT:Berendsen] temperature coupling frequency in ps: float (>0)""" ],
		[ "tmp", "300.0", "temperature in K: float" ] ]
	],

	[ "qm3.actions.dynamics.langevin_verlet", [
		[ "stp", "1000", "step number: int" ],
		[ "siz", "0.001", "step size in ps: float" ],
		[ "prt", "100", "print frequency: int" ],
		[ "gam", "50.0", "friction gamma factor in ps^-1: float" ],
		[ "tmp", "300.0", "temperature in K: float" ] ]
	] ] )


# ===============================================================================================


SP0 = re.compile( "^([0-9]+)$" )
SP1 = re.compile( "^([0-9]+)-([0-9]+)$" )
SP2 = re.compile( "^([A-Z0-9]+):([0-9]+)$" )
SP3 = re.compile( "^([A-Z0-9]+):([0-9]+)-([0-9]+)$" )
SP4 = re.compile( "^([A-Z0-9]+)/([0-9]+)/(.+)$" )
SP5 = re.compile( "^([A-Z0-9]+):([0-9]+)@([0-9\.]+)$" )


def parse_selection( mol, buf ):
	global	SP0, SP1, SP2, SP3, SP4, SP5
	sel  = [ False for i in range( mol.natm ) ]
	_not = False
	for itm in buf.split():
		# -- all
		if( itm == "*" ):
			for i in range( mol.natm ):
				sel[i] = True
		# -- negate ALL selection (at the end...)
		elif( itm == "not" ):
			_not = True
		# -- atom number (C-indexing)
		elif( SP0.match( itm ) ):
#			print( "SP0:", int(itm) )
			sel[int(itm)] = True
		# -- range of atom numbers
		elif( SP1.match( itm ) ):
			rng = [ int( i ) for i in SP1.findall( itm )[0] ]
#			print( "SP1:", rng )
			for i in range( rng[0], rng[1] + 1 ):
				sel[i] = True
		# -- residue number (by chain)
		elif( SP2.match( itm ) ):
			rng = SP2.findall( itm )[0]
#			print( "SP2:", rng )
			for i in list( mol.indx[rng[0]][int(rng[1])].values() ):
				sel[i] = True
		# -- range of residue numbers (by chain)
		elif( SP3.match( itm ) ):
			rng = SP3.findall( itm )[0]
#			print( "SP3:", rng )
			for i in range( int( rng[1] ), int( rng[2] ) + 1 ):
				for j in list( mol.indx[rng[0]][i].values() ):
					sel[j] = True
		# -- chain / residue_number / atom_label
		elif( SP4.match( itm ) ):
			rng = SP4.findall( itm )[0]
#			print( "SP4:", rng )
			sel[mol.indx[rng[0]][int(rng[1])][rng[2]]] = True
		# -- radial selection by residue around chain / residue_number
		elif( SP5.match( itm ) ):
			rng = SP5.findall( itm )[0]
#			print( "SP5:", rng )
			for i in mol.sph_sel( list( mol.indx[rng[0]][int(rng[1])].values() ), float( rng[2] ) ):
				sel[i] = True
		# -- backbone
		elif( itm == "backbone" ):
			for i in range( mol.natm ):
				if( mol.labl[i] in [ "C", "O", "CA", "N" ] ):
					sel[i] = True
	if( _not ):
		for i in range( mol.natm ):
			sel[i] = not sel[i]

	return( [ i for i in range( mol.natm ) if sel[i] ] )


# ===============================================================================================


class Hint( object ):
	def __init__( self, widget, text = "widget info", width = 800, background = "#FFFF66" ):
		self.bg = background
		self.widget = widget
		self.width = width
		self.text = text
		self.widget.bind( "<Enter>", self.showtip )
		self.widget.bind( "<Leave>", self.hidetip )
		self.tw = None

	def showtip( self, event = None ):
		x = y = 0
		x, y, cx, cy = self.widget.bbox( "insert" )
		x += self.widget.winfo_rootx() + 25
		y += self.widget.winfo_rooty() + 20
		self.hidetip()
		self.tw = tkinter.Toplevel( self.widget )
		self.tw.wm_overrideredirect( True )
		self.tw.wm_geometry( "+%d+%d"%( x, y ) )
		label = tkinter.Label( self.tw, text = self.text, justify = tkinter.LEFT,
					font = ( "Courier New", "14" ),
					background = self.bg, relief = tkinter.FLAT, borderwidth = 1,
					wraplength = self.width )
		label.pack( ipadx = 5 )

	def hidetip( self, event = None ):
		if( self.tw ):
			self.tw.destroy()
			self.tw = None


class gui_application( object ):
	def __update_canvas( self, event ):
		self.cnv.config( scrollregion = ( 0, 0, self.frm.winfo_width(), self.frm.winfo_height() ) )


	def __rem_eng( self ):
		key = self.__wid["eng"].get()
		if( key in self.__eng ):
			self.__eng[key][0].destroy()
			del self.__eng[key]


	def __add_eng( self ):
		global	ENG
		tbg = "#CCFFCC"
		key = self.__wid["eng"].get()
		if( not key in self.__eng and key in ENG ):
			frm = tkinter.LabelFrame( self.e_frm, text = key, labelanchor = tkinter.NW, borderwidth = 2,
				bg = self.__loc, relief = tkinter.GROOVE )
			frm.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 20, pady = 2 )
			self.__eng[key] = [ frm ]
			for p,d,h in ENG[key]:
				f = tkinter.Frame( frm, bg = self.__loc )
				f.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 4, pady = 4 )
				l = tkinter.Label( f, text = "%6s"%( p ), font = self.__fnt, bg = self.__loc, justify = tkinter.CENTER )
				l.pack( side = tkinter.LEFT, expand = 0, fill = None )
				e = tkinter.Entry( f, font = self.__fnt, justify = tkinter.LEFT )
				if( p != "" ):
					Hint( l, h, background = tbg )
					e.insert( 0, d )
					e.pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )
				self.__eng[key].append( ( p, e ) )


	def __rem_job( self ):
		key = self.__wid["job"].get()
		if( key in self.__job ):
			self.__job[key][0].destroy()
			del self.__job[key]


	def __add_job( self ):
		global	JOB
		tbg = "#CCFFCC"
		key = self.__wid["job"].get()
		if( not key in self.__job and key in JOB ):
			frm = tkinter.LabelFrame( self.j_frm, text = key, labelanchor = tkinter.NW, borderwidth = 2,
				bg = self.__loc, relief = tkinter.GROOVE )
			frm.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 20, pady = 2 )
			self.__job[key] = [ frm ]
			for p,d,h in JOB[key]:
				f = tkinter.Frame( frm, bg = self.__loc )
				f.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 4, pady = 4 )
				l = tkinter.Label( f, text = "%6s"%( p ), font = self.__fnt, bg = self.__loc, justify = tkinter.CENTER )
				l.pack( side = tkinter.LEFT, expand = 0, fill = None )
				e = tkinter.Entry( f, font = self.__fnt, justify = tkinter.LEFT )
				if( p != "" ):
					Hint( l, h, background = tbg )
					e.insert( 0, d )
					e.pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )
				self.__job[key].append( ( p, e ) )


	def __mkscript( self ):
		__box = None
		__sel = []
		__sqm = []
		__nbn = []
		__lnk = []
		# -- box size (if any...)
		if( self.__wid["box_x"].get() != "" ):
			__box = [ float( self.__wid["box_x"].get() ) ] * 3
			if( self.__wid["box_y"].get() != "" ):
				__box[1] = float( self.__wid["box_y"].get() )
			if( self.__wid["box_z"].get() != "" ):
				__box[2] = float( self.__wid["box_z"].get() )

		# -- selected atoms (or all of them/* defaulting to "[]")
		tmp = self.__wid["sel"].get( 1.0, tkinter.END ).strip().replace( "\n", " " )
		if( tmp != "" ):
			__sel = parse_selection( self.__mol, tmp )
			if( len( __sel ) == self.__mol.natm ):
				__sel = []
			else:
				f = open( "sele.pk", "wb" )
				pickle.dump( __sel, f )
				f.close()

		# -- QM atoms
		tmp = self.__wid["sqm"].get( 1.0, tkinter.END ).strip().replace( "\n", " " )
		if( tmp != "" ):
			__sqm = parse_selection( self.__mol, tmp )
			f = open( "sele_QM.pk", "wb" )
			pickle.dump( __sqm, f )
			f.close()

			# Link atoms
			tmp = self.__wid["lnk"].get( 1.0, tkinter.END ).strip().replace( "\n", " " )
			if( tmp != "" ):
				pmt = []
				for atm in tmp.split():
					i,j,k = atm.split( "/" )
					pmt.append( self.__mol.indx[i][int(j)][k] )
				__lnk = [ [ pmt[2*i], pmt[2*i+1] ] for i in range( len( pmt ) // 2 ) ]
			f = open( "sele_LA.pk", "wb" )
			pickle.dump( __lnk, f )
			f.close()

			# NBND interactions (remove always QM and LA)
			tmp = self.__wid["nbn"].get( 1.0, tkinter.END ).strip().replace( "\n", " " )
			if( tmp != "" ):
				__nbn = parse_selection( self.__mol, tmp )
			else:
				if( __sel != [] ):
					__nbn = __sel[:]
				else:
					__nbn = list( range( self.__mol.natm ) )
			__nbn = list( set( __nbn ).difference( set( __sqm + sum( __lnk, [] ) ) ) )
			f = open( "sele_MM.pk", "wb" )
			pickle.dump( __nbn, f )
			f.close()

			# calculate exclussion lists (calculate connectivity and assign scale factors...)
			if( __sqm != [] and __lnk != [] and __nbn != [] ):
				import	qm3.elements
				import	qm3.utils
				__exc = []
				# calculate (QM_residue + MM_residue) connectivity
				for ii,jj in __lnk:
					idx = sorted( list( set( self.__mol.indx[self.__mol.segn[ii]][self.__mol.resi[ii]].values() + self.__mol.indx[self.__mol.segn[jj]][self.__mol.resi[jj]].values() ) ) )
					nat = len( idx )
					anu = [ qm3.elements.rsymbol[self.__mol.labl[i][0]] for i in idx ]
					atm = [ [] for i in range( nat ) ]	
					smm = [ True for i in range( nat ) ]
					sqm = []
					for i in __sqm:
						try:
							smm[idx.index( i )] = False
							sqm.append( idx.index( i ) )
						except:
							pass
					for i in range( nat - 1 ):
						i3 = idx[i] * 3
						ri = qm3.elements.r_cov[anu[i]] + 0.05
						for j in range( i + 1, nat ):
							if( anu[i] == 1 and anu[j] == 1 ):
								continue
							j3 = idx[j] * 3
							rj = qm3.elements.r_cov[anu[j]] + 0.05 
							d2 = ( ri + rj ) * ( ri + rj )
							if( qm3.utils.distanceSQ( self.__mol.coor[i3:i3+3], self.__mol.coor[j3:j3+3] ) <= d2 ):
								atm[i].append( j )
								atm[j].append( i )
					# 1-2 (scale factor = 0.0)
					for i in sqm:
						for j in atm[i]:
							if( j != i and smm[j] ):
								__exc.append( [ idx[i], idx[j], 0.0 ] )
					# 1-3 (scale factor = 0.0)
					for i in sqm:
						for j in atm[i]:
							for k in atm[j]:
								if( k != i and smm[k] ):
									__exc.append( [ idx[i], idx[k], 0.0 ] )
					# 1-4 (scale factor = 0.5)
					for i in sqm:
						for j in atm[i]:
							for k in atm[j]:
								for l in atm[k]:
									if( k != i and l != j and l != i and smm[l] ):
										__exc.append( [ idx[i], idx[l], 0.5 ] )
				__exc.sort()
				f = open( "sele_EX.pk", "wb" )
				pickle.dump( __exc, f )
				f.close()

		__mod = { "qm3.mol": None, "qm3.io.dcd": None, "qm3.problem": None }
		for k in iter( self.__eng ):
			__mod[".".join( k.split( "." )[:-1])] = None
		for k in iter( self.__job ):
			__mod[".".join( k.split( "." )[:-1])] = None

		f = open( "run", "wt" )
		f.write( """from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
try:
	import cPickle as pickle
except:
	import pickle
import os
import time
""" )
		for m in iter( __mod ):
			f.write( "import %s\n"%( m ) )

		f.write( """

class my_problem( qm3.problem.template ):
	def __init__( self ):
		qm3.problem.template.__init__( self )

		self.mole = qm3.mol.molecule( "%s" )
"""%( self.__pdb ) )

		if( __box != None ):
			f.write( """
		self.mole.boxl = [ %.2lf, %.2lf, %.2lf ]
"""%( __box[0], __box[1], __box[2] ) )

		if( __sel == [] ):
			f.write( """
		self.size = 3 * self.mole.natm
		self.coor = self.mole.coor
		self.sele = []
""" )
		else:
			f.write( """
		f = open( "sele.pk", "rb" )
		self.sele = pickle.load( f )
		f.close()
		self.size = 3 * len( self.sele )
		self.coor = []
		self.mass = []
		for i in self.sele:
			i3 = i * 3
			self.coor += self.mole.coor[i3:i3+3][:]
			self.mass.append( self.mole.mass[i] )
""" )

		f.write( """
		self.func = 0.0
		self.grad = []
		self.hess = None
		self.deco = {}
""" )

		# == ENGINES ==
		qcl = ""
		who = 0
		for key in iter( self.__eng ):
			obj = { p: e.get() for p,e in self.__eng[key][1:] }
			knd = key.split( "." )[-1]

			# DYNAMO_PIPE
			if( knd == "dynamo_pipe" ):
				f.write( """
		try:
			os.unlink( "dynamo.pipe" )
		except Exception as e:
			print( e )
		os.mkfifo( "dynamo.pipe" )
		os.system( "bash r.dynamo &" )
		time.sleep( 10 )
		self.e%02d = %s()
"""%( who, key ) )
				qcl = """\t\tself.e%02d.update_chrg( self.mole )\n"""%( who )

			# PY_DYNAMO
			if( knd == "py_dynamo" ):
				f.write( """
		self.e%02d = %s( "./dynamo.so" )
"""%( who, key ) )
				qcl = """\t\tself.e%02d.update_chrg( self.mole )\n"""%( who )

			# NAMD
			if( knd == "namd" ):
				f.write( """
		self.e%02d = %s()
		self.e%02d.exe = "bash r.namd"
		self.mole.psf_read( "%s" )
"""%( who, key, who, obj["psf"] ) )
				qcl = """\t\tself.e%02d.update_chrg( self.mole, "%s" )\n"""%( who, obj["psf"] )
	
			# NAMD_PIPE / NAMD_SHM
			if( knd in [ "namd_pipe", "namd_shm" ] ):
				f.write( """
		try:
			os.unlink( "namd.pipe" )
		except Exception as e:
			print( e )
		os.mkfifo( "namd.pipe" )
		os.system( "bash r.namd &" )
		time.sleep( 10 )
		self.e%02d = %s()
		self.mole.psf_read( "%s" )
"""%( who, key, obj["psf"] ) )
				qcl = """\t\tself.e%02d.update_chrg( self.mole )\n"""%( who )

			# NAMD / NAMD_PIPE / NAMD_SHM
			if( knd in [ "namd", "namd_pipe", "namd_shm" ] ):
				import	qm3.engines.namd
				qm3.engines.namd.coordinates_write( self.__mol, "namd.coor" )
				if( __sel == [] ):
					tmp = [ False for i in range( self.__mol.natm ) ]
				else:
					tmp = [ True for i in range( self.__mol.natm ) ]
					for i in __sel:
						tmp[i] = False
				for i in __sqm:
					tmp[i] = True
				qm3.engines.namd.pdb_write( self.__mol, "namd.fix", [ i for i in range( self.__mol.natm ) if tmp[i] ] )
	
			# SANDER
			if( knd == "sander" ):
				f.write( """
		self.e%02d = %s()
		self.e%02d.exe = "bash r.sander"
		qm3.sander.topology_read( self.mole, "prmtop" )
"""%( who, key, who ) )
				qcl = """\t\tself.e%02d.update_chrg( self.mole )\n"""%( who )
	
			# PY_SANDER
			if( knd == "py_sander" ):
				if( __sqm == [] ):
					f.write( """
		self.e%02d = %s( self.mole, prmtop = "prmtop", cutoff = %.1lf, PBC = %s )
"""%( who, key, float( obj["cut"] ), __box != None ) )
				else:
					f.write( """
		f = open( "sele_QM.pk", "rb" )
		s_qm = pickle.load( f )
		f.close()
		self.e%02d = %s( self.mole, prmtop = "prmtop", cutoff = %.1lf, PBC = %s,
			qmsel = s_qm, method = "%s", charge = %d  )
		qm3.sander.topology_read( self.mole, "prmtop" )
"""%( who, key, float( obj["cut"] ), __box != None, obj["qm_met"], int( obj["qm_chg"] ) ) )

			# CHARMM_PIPE
			if( knd == "charmm_pipe" ):
				f.write( """
		try:
			os.unlink( "charmm.pipe" )
		except Exception as e:
			print( e )
		os.mkfifo( "charmm.pipe" )
		os.system( "bash r.charmm &" )
		self.e%02d = %s( "charmm.inp" )
		time.sleep( 10 )
"""%( who, key ) )
				qcl = """\t\tself.e%02d.update_chrg( self.mole )\n"""%( who )

			# CHARMM_SHM
			if( knd == "charmm_shm" ):
				f.write( """
		try:
			os.unlink( "charmm.pipe" )
		except Exception as e:
			print( e )
		os.mkfifo( "charmm.pipe" )
		os.system( "bash r.charmm &" )
		self.e%02d = %s( "charmm.inp" )
"""%( who, key ) )
				qcl = """\t\tself.e%02d.update_chrg( self.mole )\n"""%( who )

			# LAMMPS
			if( knd == "lammps" ):
				f.write( """
		self.e%02d = %s( "lammps.inp" )
		self.e%02d.exe = "bash r.lammps"
"""%( who, key, who ) )
				qcl = """\t\tself.e%02d.update_chrg( self.mole )\n"""%( who )
	
			# LAMMPS_PIPE
			if( knd == "lammps_pipe" ):
				f.write( """
		try:
			os.unlink( "lammps.pipe" )
		except Exception as e:
			print( e )
		os.mkfifo( "lammps.pipe" )
		os.system( "bash r.lammps &" )
		self.e%02d = %s( "lammps.inp" )
		time.sleep( 10 )
"""%( who, key ) )
				qcl = """\t\tself.e%02d.update_chrg( self.mole )\n"""%( who )
	
			# PY_LAMMPS
			if( knd == "py_lammps" ):
				f.write( """
		self.e%02d = %s( "lammps.inp" )
"""%( who, key ) )
				qcl = """\t\tself.e%02d.update_chrg( self.mole )\n"""%( who )

			# Common QM engines
			if( knd in [ "dftb", "dl_dftb", "sqm", "gamess", "nwchem", "demon", "orca", "py_psi4", 
					"gaussian", "gaussian_MMEL", "tchem_sckt", "tchem", "smash" ] ):
				f.write( """
		f = open( "sele_QM.pk", "rb" )
		s_qm = pickle.load( f )
		f.close()
		f = open( "sele_LA.pk", "rb" )
		s_la = pickle.load( f )
		f.close()
		f = open( "sele_MM.pk", "rb" )
		s_mm = pickle.load( f )
		f.close()
		if( self.mole.chrg == [] ):
			self.mole.chrg = [ 0.0 for i in range( self.mole.natm ) ]
		for i in s_qm:
			self.mole.chrg[i] = 0.0
""" )
				f.write( qcl )
	
			# DFTB
			if( knd == "dftb" ):
				f.write( """
		self.e%02d = %s( self.mole, s_qm, s_mm, s_la )
		self.e%02d.exe = "bash r.dftb"
		self.e%02d.chg = %d
		self.e%02d.prm = "%s/"
"""%( who, key, who, who, int( obj["chg"] ), who, obj["prm"] ) )
	
			# DL_DFTB
			if( knd == "dl_dftb" ):
				f.write( """
		self.e%02d = %s( self.mole, s_qm, s_mm, s_la )
		self.e%02d.chg = %d
		self.e%02d.prm = "%s/"
"""%( who, key, who, int( obj["chg"] ), who, obj["prm"] ) )
	
			# TCHEM_SCKT
			if( knd == "tchem_sckt" ):
				f.write( """
		self.e%02d = %s( self.mole, "sckt_", s_qm, s_mm, s_la )
"""%( who, key ) )
	
			# SQM / GAMESS / NWCHEM / DEMON / ORCA / TCHEM / SMASH
			if( knd in [ "sqm", "gamess", "nwchem", "demon", "orca", "tchem", "smash" ] ):
				f.write( """
		f = open( "i.%s", "rt" )
		i_qm = f.read()
		f.close()
		self.e%02d = %s( self.mole, i_qm, s_qm, s_mm, s_la )
		self.e%02d.exe = "bash r.%s"
"""%( knd, who, key, who, knd ) )

			# GAUSSIANs
			if( knd in [ "gaussian", "gaussian_MMEL" ] ):
				f.write( """
		f = open( "i.gauss.ini", "rt" )
		i_qm = f.read()
		f.close()
		f = open( "i.gauss.mid", "rt" )
		m_qm = f.read()
		f.close()
		f = open( "i.gauss.end", "rt" )
		e_qm = f.read()
		f.close()
		self.e%02d = %s( self.mole, i_qm, m_qm, e_qm, s_qm, s_mm, s_la )
		self.e%02d.exe = "bash r.g09"
"""%( who, key, who ) )
	
			# PSI4
			if( knd == "py_psi4" ):
				f.write( """
		o_qm = %s
		self.e%02d = %s( self.mole, s_qm, o_qm, s_mm, s_la )
"""%( obj["opt"], who, key ) )
	
			# Int_QMLJ / _MMEL
			if( knd in [ "Int_QMLJ", "Int_QMLJ_MMEL" ] ):
				f.write( """
		self.mole.nbnd_read( "%s" )
		f = open( "sele_EX.pk", "rb" )
		exc = pickle.load( f )
		f.close()
		self.e%02d = %s( self.mole, s_qm, s_mm, exc )
"""%( obj["typ"], who, key ) )

			# DISTANCE / ANGLE
			if( knd in [ "distance", "angle" ] ):
				f.write( """
		self.e%02d = %s( %s, %s, %s )
"""%( who, key, obj["kmb"], obj["ref"], str( [ int( i ) for i in obj["idx"].split() ] ) ) )

			# MULTIPLE_DISTANCE
			if( knd == "multiple_distance" ):
				f.write( """
		self.e%02d = %s( %s, %s, %s, %s )
"""%( who, key, obj["kmb"], obj["ref"], str( [ int( i ) for i in obj["idx"].split() ] ),
		str( [ float( i ) for i in obj["wei"].split() ] ) ) )
	
			who += 1

		# end constructor
		f.write( """
		self.mass = self.mole.mass

		self.flog = open( "log", "wt" )

		self.fdcd = qm3.io.dcd.dcd()
		self.fdcd.write( "dcd", self.mole.natm, self.sele )
		self.qdcd = False


	def log( self, buf ):
		self.flog.write( buf + "\\n" )
		self.flog.flush()


	def current_step( self, istep ):
		if( self.qdcd and istep % 1 == 0 ):
			self.mole.dcd_write( self.fdcd )


	def stop( self ):
		self.flog.close()
		self.fdcd.close()
""" )
		who = 0
		for key in iter( self.__eng ):
			if( key.split( "." )[-1] in [ "namd_pipe", "charmm_pipe", "dynamo_pipe", "py_lammps", "lammps_pipe" ] ):
				f.write( """\t\tself.e%02d.stop()\n"""%( who ) )
			who += 1

		f.write( """

	def update_coor( self ):
		for i in range( len( self.sele ) ):
			ii = 3 * self.sele[i]
			jj = 3 * i
			for j in [ 0, 1, 2 ]:
				self.mole.coor[ii+j] = self.coor[jj+j]


	def get_func( self ):
		self.update_coor()
""" )
		who = 0
		for key in iter( self.__eng ):
			if( not key.split( "." )[-1] in [ "Int_QMLJ", "Int_QMLJ_MMEL" ] ):
				f.write( """\t\tself.mole.func = 0.0
		self.e%02d.get_func( self.mole )
		self.deco["e%02d"] = self.mole.func
"""%( who, who ) )
			who += 1
		f.write( """\t\tself.func = sum( self.deco.values() )


	def get_grad( self ):
		self.update_coor()
		self.mole.grad = [ 0.0 for i in range( 3 * self.mole.natm ) ]
""" )
		who = 0
		for key in iter( self.__eng ):
			f.write( """\t\tself.mole.func = 0.0
		self.e%02d.get_grad( self.mole )
		self.deco["e%02d"] = self.mole.func
"""%( who, who ) )
			who += 1
		f.write( """\t\tself.func = sum( self.deco.values() )
		if( self.sele == [] ):
			self.grad = self.mole.grad[:]
		else:
			self.grad = []
			for i in self.sele:
				i3 = i * 3
				self.grad += self.mole.grad[i3:i3+3][:]




obj = my_problem()

""" )

		for key in iter( self.__job ):
			obj = { p: e.get() for p,e in self.__job[key][1:] }
			knd = key.split( "." )[-1]

			if( knd in [ "steepest_descent", "fire", "l_bfgs" ] ):
				f.write( """%s( obj, step_number = %s, step_size = %s, print_frequency = %s, gradient_tolerance = %s, log_function = obj.log )\n"""%(
					key, obj["stp"], obj["siz"], obj["prt"], obj["tol"] ) )

			elif( knd == "conjugate_gradient_plus" ):
				f.write( """%s( obj, step_number = %s, print_frequency = %s, gradient_tolerance = %s, log_function = obj.log )\n"""%(
					key, obj["stp"], obj["siz"], obj["prt"], obj["tol"] ) )

			elif( knd == "velocity_verlet" ):
				f.write( """qm3.actions.dynamics.assign_velocities( obj, temperature = %s )
%s( obj, step_number = %s, step_size = %s, print_frequency = %s, scale_frequency = %s, temperature_coupling = %s, temperature = %s, log_function = obj.log ) 
"""%( obj["tmp"], key, obj["stp"], obj["siz"], obj["prt"], obj["scl"], obj["cpl"], obj["tmp"] ) )

			elif( knd == "langevin_verlet" ):
				f.write( """qm3.actions.dynamics.assign_velocities( obj, temperature = %s )
%s( obj, step_number = %s, step_size = %s, print_frequency = %s, gamma_factor = %s, temperature = %s, log_function = obj.log ) 
"""%( obj["tmp"], key, obj["stp"], obj["siz"], obj["prt"], obj["gam"], obj["tmp"] ) )

		f.write( """
obj.mole.pdb_write( "last.pdb" )
f = open( "last.pk", "wb" )
pickle.dump( obj.mole, f )
f.close()
obj.stop()
""" )
		f.close()
		self.app.destroy()


	def __init__( self, pdb_name ):
		global	ENG, JOB
		self.__col = "#E6E6E6"
		self.__loc = "#CCCCCC"
		self.__fnt = ( "Courier New", "14" )
		self.__wid = {}
		# -------------------------------------------------
		self.__eng = collections.OrderedDict()
		self.__job = collections.OrderedDict()
		self.__pdb = pdb_name
		self.__mol = qm3.mol.molecule( self.__pdb )
		# -------------------------------------------------

		self.app = tkinter.Tk()
		self.app.title( "-- Helper --" )
		self.app.maxsize( width = 560, height = 2048 )
		self.app.minsize( width = 560, height = 480 )

		self.cnv = tkinter.Canvas( self.app, bg = self.__col )
		self.frm = tkinter.Frame( self.cnv, bg = self.__col )
		self.cnv.create_window( ( 0, 0 ), window = self.frm, anchor = tkinter.NW )
		self.vsb = tkinter.Scrollbar( self.app, orient = tkinter.VERTICAL, command = self.cnv.yview )
		self.vsb.pack( side = tkinter.RIGHT, fill = tkinter.Y )
		self.cnv.config( yscrollcommand = self.vsb.set )
		self.cnv.pack( side = tkinter.LEFT, expand = True, fill = tkinter.BOTH )
		self.frm.bind( "<Configure>", self.__update_canvas )

		f01 = tkinter.Frame( self.frm, bg = self.__col )
		f01.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 4, pady = 4 )
		l01 = tkinter.Label( f01, text = " box ", font = self.__fnt, bg = self.__col, justify = tkinter.CENTER )
		l01.pack( side = tkinter.LEFT, expand = 0, fill = None )
		Hint( l01, "Box size: fill only the first one (X) for a CUBIC symmetry" )
		self.__wid["box_x"] = tkinter.Entry( f01, font = self.__fnt, justify = tkinter.LEFT )
		self.__wid["box_x"].pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )
		self.__wid["box_y"] = tkinter.Entry( f01, font = self.__fnt, justify = tkinter.LEFT, width = 10 )
		self.__wid["box_y"].pack( side = tkinter.LEFT, expand = 0, fill = tkinter.X )
		self.__wid["box_z"] = tkinter.Entry( f01, font = self.__fnt, justify = tkinter.LEFT, width = 10 )
		self.__wid["box_z"].pack( side = tkinter.LEFT, expand = 0, fill = tkinter.X )

		f03 = tkinter.Frame( self.frm, bg = self.__col )
		f03.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 4, pady = 4 )
		l03 = tkinter.Label( f03, text = " sel ", font = self.__fnt, bg = self.__col, justify = tkinter.CENTER )
		l03.pack( side = tkinter.LEFT, expand = 0, fill = None )
		Hint( l03, """Selection of active atoms:

*            all atoms
int          atom number (C-indexing)
int-int      range of atom numbers
str:int      residue number (int) from chain (str)
str:int-int  range of residue numbers from chain (str)
str/int/lbl  atom label (lbl) from residue number (int) of chain (str)
str:int@rad  selection by residue around radii (rad) of residue number (int) of chain (str)
backbone     atoms from protein backbone: atom labels C, O, N, CA
not          found at any position negates the resulting selection""" )
		self.__wid["sel"] = tkinter.Text( f03, font = self.__fnt, height = 5, width = 50 )
		self.__wid["sel"].pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )

		f04 = tkinter.Frame( self.frm, bg = self.__col )
		f04.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 4, pady = 4 )
		l04 = tkinter.Label( f04, text = " sqm ", font = self.__fnt, bg = self.__col, justify = tkinter.CENTER )
		l04.pack( side = tkinter.LEFT, expand = 0, fill = None )
		Hint( l04, """Selection of QM atoms:

*            all atoms
int          atom number (C-indexing)
int-int      range of atom numbers
str:int      residue number (int) from chain (str)
str:int-int  range of residue numbers from chain (str)
str/int/lbl  atom label (lbl) from residue number (int) of chain (str)
str:int@rad  selection by residue around radii (rad) of residue number (int) of chain (str)
backbone     atoms from protein backbone: atom labels C, O, N, CA
not          found at any position negates the resulting selection""" )
		self.__wid["sqm"] = tkinter.Text( f04, font = self.__fnt, height = 4, width = 50 )
		self.__wid["sqm"].pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )

		f05 = tkinter.Frame( self.frm, bg = self.__col )
		f05.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 4, pady = 4 )
		l05 = tkinter.Label( f05, text = " lnk ", font = self.__fnt, bg = self.__col, justify = tkinter.CENTER )
		l05.pack( side = tkinter.LEFT, expand = 0, fill = None )
		Hint( l05, """Selection of Link Atoms, in the form of 'QM MM' pairs list with format:

str/int/lbl  atom label (lbl) from residue number (int) of chain (str)""", 800 )
		self.__wid["lnk"] = tkinter.Text( f05, font = self.__fnt, height = 3, width = 50 )
		self.__wid["lnk"].pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )

		f06 = tkinter.Frame( self.frm, bg = self.__col )
		f06.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 4, pady = 4 )
		l06 = tkinter.Label( f06, text = " nbn ", font = self.__fnt, bg = self.__col, justify = tkinter.CENTER )
		l06.pack( side = tkinter.LEFT, expand = 0, fill = None )
		Hint( l06, """Selection of QM-interacting (non-bonded) atoms:

*            all atoms
int          atom number (C-indexing)
int-int      range of atom numbers
str:int      residue number (int) from chain (str)
str:int-int  range of residue numbers from chain (str)
str/int/lbl  atom label (lbl) from residue number (int) of chain (str)
str:int@rad  selection by residue around radii (rad) of residue number (int) of chain (str)
backbone     atoms from protein backbone: atom labels C, O, N, CA
not          found at any position negates the resulting selection""" )
		self.__wid["nbn"] = tkinter.Text( f06, font = self.__fnt, height = 4, width = 50 )
		self.__wid["nbn"].pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )

		self.e_frm = tkinter.Frame( self.frm, bg = self.__col )
		self.e_frm.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 4, pady = 4 )
		f09 = tkinter.Frame( self.e_frm, bg = self.__col )
		f09.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH )
		l09 = tkinter.Label( f09, text = " eng ", font = self.__fnt, bg = self.__col, justify = tkinter.CENTER )
		l09.pack( side = tkinter.LEFT, expand = 0, fill = None )
		Hint( l09, """Select engines to work on the molecule:

press [+] button to add it with specific options
press [-] button for removing it""" )
		self.__wid["eng"] = ttk.Combobox( f09, font = self.__fnt, justify = tkinter.LEFT, state = "readonly", width = 40 )
		self.__wid["eng"]["values"] = list( ENG )
		self.__wid["eng"].pack( side = tkinter.LEFT, expand = 0, fill = tkinter.X )
		b01 = tkinter.Button( f09, text = "[+]", font = self.__fnt, justify = tkinter.CENTER, width = 4, command = self.__add_eng )
		b01.pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )
		b02 = tkinter.Button( f09, text = "[-]", font = self.__fnt, justify = tkinter.CENTER, width = 4, command = self.__rem_eng )
		b02.pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )

		self.j_frm = tkinter.Frame( self.frm, bg = self.__col )
		self.j_frm.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 4, pady = 4 )
		f10 = tkinter.Frame( self.j_frm, bg = self.__col )
		f10.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH )
		l10 = tkinter.Label( f10, text = " job ", font = self.__fnt, bg = self.__col, justify = tkinter.CENTER )
		l10.pack( side = tkinter.LEFT, expand = 0, fill = None )
		Hint( l10, """Select jobs to work on the molecule + engines:

press [+] button to add it with specific options
press [-] button for removing it""" )
		self.__wid["job"] = ttk.Combobox( f10, font = self.__fnt, justify = tkinter.LEFT, state = "readonly", width = 40 )
		self.__wid["job"]["values"] = list( JOB )
		self.__wid["job"].pack( side = tkinter.LEFT, expand = 0, fill = tkinter.X )
		b03 = tkinter.Button( f10, text = "[+]", font = self.__fnt, justify = tkinter.CENTER, width = 4, command = self.__add_job )
		b03.pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )
		b04 = tkinter.Button( f10, text = "[-]", font = self.__fnt, justify = tkinter.CENTER, width = 4, command = self.__rem_job )
		b04.pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )

		bbb = tkinter.Button( self.frm, text = " build! ", font = self.__fnt, justify = tkinter.CENTER, width = 4, command = self.__mkscript )
		bbb.pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )

		self.app.mainloop()






gui_application( sys.argv[1] )

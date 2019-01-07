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
	[ "qm3.engines.namd.namd", [ 
		[ "exe", "", "path of the NAMD executable, without options: str" ],
		[ "psf", "", "name of the PSF file used to define the model: str" ],
		[ "cut", "", "force-switching (floats): cut-on cut-off cut-list" ] ]
	],

	[ "qm3.engines.namd.namd_pipe", [ 
		[ "exe", "", "path of the NAMD executable, without options: str" ],
		[ "psf", "", "name of the PSF file used to define the model: str" ],
		[ "cut", "", "force-switching (floats): cut-on cut-off cut-list" ] ]
	],

	[ "qm3.engines.sander.sander", [
		[ "dir", "", "AMBERHOME path: str" ],
		[ "cut", "", "cut-off: float" ],
		[ "pbc", "", "periodic boundary conditions (ntb): on/off" ],
		[ "qm_sel", "", "selection mask of sander to set up the QM atoms: str" ],
		[ "qm_met", "", "semi-empirical hamiltionan inplemented in sander: str" ],
		[ "qm_chg", "", "charge of the QM atoms: int" ] ]
	],

	[ "qm3.engines.sander.py_sander", [
		[ "cut", "", "cut-off: float" ],
		[ "pbc", "", "periodic boundary conditions (ntb): on/off" ],
		[ "qm_met", "", "semi-empirical hamiltionan inplemented in sander: str" ],
		[ "qm_chg", "", "charge of the QM atoms: int" ] ]
	],

	[ "qm3.engines.dftb.dftb", [ 
		[ "exe", "", "path of the DFTB+ executable: str" ],
		[ "chg", "", "charge of the QM atoms: int" ],
		[ "prm", "", "folder of DFTB+ parameters set: str" ] ]
	],

	[ "qm3.engines.sqm.sqm", [ 
		[ "dir", "", "AMBERHOME path: str" ],
		[ "ini", "", """SQM(QM) flags: qm_theory = "DFTB", qmcharge = 0""" ] ]
	],

	[ "qm3.engines.gamess.gamess", [ 
		[ "exe", "", "path of the GAMESS-US launching script: str" ],
		[ "ini", "", "name of file containing the GAMESS-US input template: str" ] ]
	],

	[ "qm3.engines.nwchem.nwchem", [ 
		[ "exe", "", "path of the NWCHEM launching script: str" ],
		[ "ini", "", "name of file containing the NWCHEM input template: str" ] ]
	],

	[ "qm3.engines.demon.demon", [ 
		[ "exe", "", "path of the deMon2k launching script: str" ],
		[ "ini", "", "name of file containing the deMon2k input template: str" ] ]
	],

	[ "qm3.engines.orca.orca", [ 
		[ "exe", "", "path of the ORCA lunching script: str" ],
		[ "ini", "", "name of file containing the ORCA input template: str" ] ]
	],
                                 
    [ "qm3.engines.tchem.tchem", [
        [ "ini", "", "name of file containing the TeraChem input template: str" ] ]
    ],

	[ "qm3.engines.tchem.tchem_sckt", [ 
		[ "sck", "sckt_", "name of the UNIX socket linked to TeraChem: str" ] ]
	],

	[ "qm3.engines.gaussian.gaussian", [ 
		[ "exe", "", "path of the G09 launching script: str" ],
		[ "ini", "", "name of file containing the G09 head input template: str" ],
		[ "mid", "", "name of file containing the G09 middle input template: str" ],
		[ "end", "", "name of file containing the G09 last input template: str" ] ]
	],

	[ "qm3.engines.gaussian.gaussian_MMEL", [ 
		[ "exe", "", "path of the G09 launching script: str" ],
		[ "ini", "", "name of file containing the G09 head input template: str" ],
		[ "mid", "", "name of file containing the G09 middle input template: str" ],
		[ "end", "", "name of file containing the G09 last input template: str" ] ]
	],

	[ "qm3.engines.xpsi4.Psi4", [ 
		[ "opt", "", "python dictonary with Psi-4 options: { key: val }" ] ]
	],

	[ "qm3.engines._qmmm.Int_QMLJ", [
		[ "exc", "", "exclusion [QM (in sqm), MM (in nbn), scale] interaction triplets: C-index C-index float" ] ]
	],

	[ "qm3.engines._qmmm.Int_QMLJ_MMEL", [
		[ "exc", "", "exclusion [QM (in sqm), MM (in nbn), scale) interaction triplets: C-index C-index float" ] ]
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
		[ "scl", "100", """temperature scaling frequency in steps: int
    if scl <= 0 performs NVE""" ],
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


##################################################################################################
# Selections
#
SP0 = re.compile( "^([0-9]+)$" )
SP1 = re.compile( "^([0-9]+)-([0-9]+)$" )
SP2 = re.compile( "^([A-Z0-9]+):([0-9]+)$" )
SP3 = re.compile( "^([A-Z0-9]+):([0-9]+)-([0-9]+)$" )
SP4 = re.compile( "^([A-Z0-9]+)/([0-9]+)/(.+)$" )
SP5 = re.compile( "^([A-Z0-9]+):([0-9]+)@([0-9\.]+)$" )


def apply_selection( mol, sel, lst ):
	global	SP0, SP1, SP2, SP3, SP4, SP5
	_not = False
	for itm in lst:
		# -- all
		if( itm == "*" ):
			for i in range( mol.natm ):
				sel[i] = True
		# -- negate ALL selection (at the end...)
		elif( itm == "not" ):
			_not = True
		# -- atom number (C-indexing)
		elif( SP0.match( itm ) ):
			print( "SP0:", int(itm) )
			sel[int(itm)] = True
		# -- range of atom numbers
		elif( SP1.match( itm ) ):
			rng = [ int( i ) for i in SP1.findall( itm )[0] ]
			print( "SP1:", rng )
			for i in range( rng[0], rng[1] + 1 ):
				sel[i] = True
		# -- residue number (by chain)
		elif( SP2.match( itm ) ):
			rng = SP2.findall( itm )[0]
			print( "SP2:", rng )
			for i in list( mol.indx[rng[0]][int(rng[1])].values() ):
				sel[i] = True
		# -- range of residue numbers (by chain)
		elif( SP3.match( itm ) ):
			rng = SP3.findall( itm )[0]
			print( "SP3:", rng )
			for i in range( int( rng[1] ), int( rng[2] ) + 1 ):
				for j in list( mol.indx[rng[0]][i].values() ):
					sel[j] = True
		# -- chain / residue_number / atom_label
		elif( SP4.match( itm ) ):
			rng = SP4.findall( itm )[0]
			print( "SP4:", rng )
			sel[mol.indx[rng[0]][int(rng[1])][rng[2]]] = True
		# -- radial selection by residue around chain / residue_number
		elif( SP5.match( itm ) ):
			rng = SP5.findall( itm )[0]
			print( "SP5:", rng )
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


def flush_selections( conf ):
	f = open( "sele.pk", "wb" )
	pickle.dump( [ i for i in range( conf["MOL"].natm ) if conf["sel"][i] ], f )
	f.close()
	if( sum( conf["sqm"] ) > 0 ):
		f = open( "sele_QM.pk", "wb" )
		pickle.dump( [ i for i in range( conf["MOL"].natm ) if conf["sqm"][i] ], f )
		f.close()
		f = open( "sele_LA.pk", "wb" )
		pickle.dump( conf["lnk"], f )
		f.close()
		f = open( "sele_MM.pk", "wb" )
		pickle.dump( [ i for i in range( conf["MOL"].natm ) if conf["nbn"][i] ], f )
		f.close()



##################################################################################################
# Supporting templates...
#
def __common_namd( conf, opt ):
	g = open( "namd.inp", "wt" )
	g.write( """structure           %s
coordinates         %s
bincoordinates      namd.coor

paraTypeCharmm      on
parameters          PARAMETERS

fixedatoms          on
fixedatomsfile      namd.fix

extrabonds          off
extrabondsfile      namd.ic
"""%( opt["psf"], conf["mol"]["fname"] ) )
	if( conf["box"] != None ):
		g.write( """
cellBasisVector1    %.2lf   0.00   0.00
cellBasisVector2     0.00  %.2lf   0.00
cellBasisVector3     0.00   0.00  %.2lf

PME                 on
PMETolerance        0.000001
PMEGridSpacing      0.5
"""%( conf["box"][0], conf["box"][1], conf["box"][2] ) )
	g.write( """
exclude             scaled1-4
1-4scaling          0.5
switching           on
switchdist          %.1lf
cutoff              %.1lf
pairlistdist        %.1lf

wrapAll             off
wrapWater           off
nonbondedFreq       1
fullElectFrequency  1
stepspercycle       1
temperature         0.0
outputEnergies      1
outputname          namd.out
"""%( float( opt["cut"][0] ), float( opt["cut"][1] ), float( opt["cut"][2] ) ) )
	g.close()


def input_namd( conf, opt ):
	__common_namd( conf, opt )
	g = open( "namd.inp", "at" )
	g.write( """forcedcdfile        namd.force
forcedcdfreq        1
run                 0
output onlyforces   namd
""" )
	g.close()


def input_namd_pipe( conf, opt ):
	__common_namd( conf, opt )
	g = open( "namd.inp", "at" )
	g.write( """startup

set fd [ open "namd.pipe" r ]
while { [ gets $fd cmd ] >= 0 } {
    switch $cmd {
        "energy"      { run 0 }
        "gradient"    { run 0; output onlyforces namd }
        "charges"     { reloadCharges namd.chrg }
        "coordinates" { coorfile binread namd.coor }
        "exit"        { close $fd; exit }
    }
}
""" )
	g.close()


def input_sander( conf, opt ):
	g = open( "mdin", "wt" )
	g.write( """MM_slave (check: ntb & ibelly/bellymask)
&cntrl
 imin   = 0, ntx    = 1, irest  = 0,
 ntxo   = 1, ntpr   = 1, ntave  = 0,
 ntwr   = 1, iwrap  = 0, ntwx   = 0,
 ntwv   = 0, ntwf   = 1, ntwe   = 1,
 ioutfm = 1, ntwprt = 0, ntt    = 0,
 ntr    = 0, nstlim = 0, nscm   = 0,
 ntp    = 0, ifqnt  = %d,
 cut       = %.1lf,
 ntb       = %d,
 ibelly    = 1,
 bellymask = 'SANDER_MOVILE_SELECTION',
/
"""%( opt["qm_sel"] != None, float( opt["cut"] ), opt["pbc"] ) )
	if( opt["qm_sel"] != None ):
		g.write( """&qmmm
 qmmask    = '%s',
 qm_theory = '%s',
 qmcharge  = %d
/
"""%( opt["qm_sel"], opt["qm_met"], int( opt["qm_chg"] ) ) )
	g.close()



###################################################################################################
## Parsing...
##
#def __parse_options( cur, opt, siz, lin ):
#	i = cur + 1
#	if( i < siz ):
#		t = lin[i].strip().split()
#	while( i < siz and len( t ) > 0 ):
#		opt[t[0]] = " ".join( t[1:] )
#		i += 1
#		if( i < siz ):
#			t = lin[i].strip().split()
#	return( i )
#
#
#
#def parse_input( conf, lines ):
#	global	SP0, SP1, SP2, SP3, SP4, SP5, ENG, JOB
#	cur = 0
#	siz = len( lines )
#	while( cur < siz ):
#		t = lines[cur].strip().split()
#		if( len( t ) > 0 ):
#
#			# -- box size
#			if( t[0] == "box" and len( t ) == 4 ):
#				conf["box"] = [ float( j ) for j in t[1:] ]
#		
#			# -- molecule object
#			if( t[0] == "mol" and len( t ) ==  3 ):
#				tt = t[1].split( "." )
#				conf["imp"][tt[0]] = None
#				conf["mol"]["module"] = tt[0]
#				conf["mol"]["method"] = tt[1]
#				conf["mol"]["fname"]  = t[2]
#				if( conf["mol"]["module"] == "mol_io" ):
#					import mol_io
#					conf["MOL"] = mol_io.molecule()
#				elif( conf["mol"]["module"] == "x_namd" ):
#					import x_namd
#					conf["MOL"] = x_namd.molecule()
#				elif( conf["mol"]["module"] == "x_dynamo" ):
#					import x_dynamo
#					conf["MOL"] = x_dynamo.molecule()
#				elif( conf["mol"]["module"] == "x_sander" ):
#					import x_sander
#					conf["MOL"] = x_sander.molecule()
#				eval( "conf[\"MOL\"]." + conf["mol"]["method"] + "( \"" + conf["mol"]["fname"] + "\" )" )
#				if( conf["box"] != None ):
#					conf["MOL"].boxl = conf["box"][:]
#				conf["sel"] = [ True  for i in range( conf["MOL"].natm ) ]
#				conf["sqm"] = [ False for i in range( conf["MOL"].natm ) ]
#				conf["nbn"] = [ False for i in range( conf["MOL"].natm ) ]
#
#			# -- molecule output
#			if( t[0] == "out" and len( t ) ==  3 ):
#				conf["out"] = { "method": t[1], "fname": t[2] }
#
#			# -- trajectory
#			if( t[0] == "dcd" and len( t ) ==  3 ):
#				conf["imp"]["mol_io"] = None
#				conf["dcd"] = { "fname": t[1], "freq": int( t[2] ) }
#
#			# -- active atoms
#			if( t[0] == "sel" ):
#				conf["sel"] = [ False for i in range( conf["MOL"].natm ) ]
#				apply_selection( conf["MOL"], conf["sel"], t[1:] )
#	
#			# -- QM atoms: define PREVIOUS to LNK and NBN
#			if( t[0] == "sqm" ):
#				conf["sqm"] = [ False for i in range( conf["MOL"].natm ) ]
#				apply_selection( conf["MOL"], conf["sqm"], t[1:] )
#	
#			# -- Link-atoms bonds: define PREVIOUS to NBN
#			if( t[0] == "lnk" ):
#				tmp = []
#				for atm in t[1:]:
#					i,j,k = atm.split( "/" )
#					tmp.append( conf["MOL"].indx[i][int(j)][k] )
#				conf["lnk"] = [ [ tmp[2*i], tmp[2*i+1] ] for i in range( len( tmp ) // 2 ) ]
#	
#			# -- Non-bonded atoms
#			if( t[0] == "nbn" ):
#				conf["nbn"] = [ False for i in range( conf["MOL"].natm ) ]
#				apply_selection( conf["MOL"], conf["nbn"], t[1:] )
#				for i in range( conf["MOL"].natm ):
#					if( conf["sqm"][i] ):
#						conf["nbn"][i] = False
#				for i,j in conf["lnk"]:
#					conf["nbn"][j] = False
#
#			# -- engines
#			if( t[0] == "eng" and len( t ) == 2 ):
#	
#				if( t[1] == "x_namd.namd" or t[1] == "x_namd.namd_pipe" ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in ENG[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					obj["exe"] += " +setcpuaffinity +isomalloc_sync +idlepoll namd.inp > namd.out"
#					if( obj["class"] == "x_namd.namd_pipe" ):
#						obj["exe"] += " &"
#					obj["cut"] = [ float( i ) for i in obj["cut"].split() ]
#					conf["eng"].append( obj.copy() )
#					conf["MOL"]._occ = []
#					for i in range( conf["MOL"].natm ):
#						if( not conf["sel"][i] or conf["sqm"][i] ):
#							conf["MOL"]._occ.append( 1 )
#						else:
#							conf["MOL"]._occ.append( 0 )
#					conf["MOL"].pdb_write( "namd.fix" )
#					try:
#						conf["MOL"].namd_write( "namd.coor" )
#					except Exception as e:
#						print( "-- ERROR --\n" )
#						print( e )
#						print( "\n[try] mol x_namd.%s %s\n"%( conf["mol"]["method"], conf["mol"]["fname"] ) )
#						sys.exit( 1 )
#	
#				elif( t[1] == "x_sander.Sander" ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in ENG[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					obj["dir"] = "export AMBERHOME=%s; $AMBERHOME/bin/sander -O"%( obj["dir"] )
#					obj["pbc"] = obj["pbc"] == "on"
#					if( conf["box"] == None ):
#						conf["box"] = [ 1000, 1000, 1000 ]
#					conf["eng"].append( obj.copy() )
#					try:
#						conf["MOL"].inpcrd_write( "inpcrd" )
#					except Exception as e:
#						print( "-- ERROR --\n" )
#						print( e )
#						print( "\n[try] mol x_sander.%s %s\n"%( conf["mol"]["method"], conf["mol"]["fname"] ) )
#						sys.exit( 1 )
#	
#				elif( t[1] == "x_sander.py_sander" ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in ENG[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					obj["pbc"] = obj["pbc"] == "on"
#					if( conf["box"] == None ):
#						conf["box"] = [ 100, 100, 100 ]
#					conf["eng"].append( obj.copy() )
#	
#				elif( t[1] == "x_dftb.dftb" ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in ENG[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					obj["exe"] += " > dftb_in.log"
#					if( obj["prm"][-1] != "/" ):
#						obj["prm"] += "/"
#					conf["eng"].append( obj.copy() )
#	
#				elif( t[1] == "x_sqm.sqm" ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in ENG[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					obj["dir"] = "export AMBERHOME=%s; $AMBERHOME/bin/sqm"%( obj["dir"] )
#					conf["eng"].append( obj.copy() )
#	
#				elif( t[1] in [ "x_gamess.gamess", "x_nwchem.nwchem", "x_demon.demon", "x_orca.orca",
#						"x_gaussian.gaussian", "x_gaussian.gaussian_MMEL", "x_psi4.Psi4",
#						"x_tchem.tchem_sckt", "x_tchem.tchem" ] ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in ENG[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					conf["eng"].append( obj.copy() )
#	
#				elif( t[1] in [ "_qmmm.Int_QMLJ", "_qmmm.Int_QMLJ_MMEL" ] ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in ENG[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					obj["par"] = obj["par"].split()
#					obj["exc"] = obj["exc"].split()
#					tmp = []
#					for i in range( len( obj["exc"] ) // 3 ):
#						i3 = i * 3
#						buf = []
#						for j in [ 0, 1 ]:
#							if( SP0.match( obj["exc"][i3+j] ) ):
#								buf.append( obj["exc"][i3+j] )
#							elif( SP4.match( obj["exc"][i3+j] ) ):
#								k = SP4.findall( obj["exc"][i3+j] )[0]
#								buf.append( str( conf["MOL"].indx[k[0]][int(k[1])][k[2]] ) )
#						if( buf != [] ):
#							tmp.append( "[ " + buf[0] + ", " + buf[1] + ", " + obj["exc"][i3+2] + " ]" )
#					obj["exc"] = "[ " + ",".join( tmp ) + " ]"
#					conf["eng"].append( obj.copy() )
#	
#				elif( t[1] in [ "restraints.distance", "restraints.angle", "restraints.multiple_distance" ] ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in ENG[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					tmp = []
#					for itm in obj["idx"].split():
#						if( SP0.match( itm ) ):
#							tmp.append( itm )
#						elif( SP4.match( itm ) ):
#							k = SP4.findall( itm )[0]
#							tmp.append( str( conf["MOL"].indx[k[0]][int(k[1])][k[2]] ) )
#					obj["idx"] = "[ " + ",".join( tmp ) + " ]"
#					if( t[1] == "restraints.multiple_distance" ):
#						obj["wei"] = "[ " + ",".join( obj["wei"].split() ) + " ]"
#					conf["eng"].append( obj.copy() )
#	
#			# -- jobs
#			if( t[0] == "job" and len( t ) == 2 ):
#	
#				if( t[1] in [ "minimize.steepest_descent", "minimize.fire", "minimize.l_bfgs",
#					"minimize.conjugate_gradient_plus" ] ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in JOB[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					conf["job"].append( obj.copy() )
#	
#				if( t[1] == "dynamics.velocity_verlet" ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in JOB[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					conf["job"].append( obj.copy() )
#	
#				if( t[1] == "dynamics.langevin_verlet" ):
#					conf["imp"][t[1].split( "." )[0]] = None
#					obj = { p: d for p,d,h in JOB[t[1]] }; obj["class"] = t[1]
#					cur = __parse_options( cur, obj, siz, lines )
#					conf["job"].append( obj.copy() )
#
#		# -- cycle...
#		cur += 1
#
#	# -- summary
#	print( 80*"-" )
#	print( "box", conf["box"] )
#	print( "imp", conf["imp"] )
#	print( "mol", conf["mol"] )
#	print( "MOL", conf["MOL"] )
#	print( "dcd", conf["dcd"] )
#	print( "out", conf["out"] )
#	print( "sel", sum( conf["sel"] ) )
#	print( "sqm", sum( conf["sqm"] ) )
#	print( "lnk", conf["lnk"] )
#	print( "nbn", sum( conf["nbn"] ) )
#	print( "eng", conf["eng"] )
#	print( "job", conf["job"] )
#	print( 80*"-" )



##################################################################################################
# Translate...
#
def translate( conf ):
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

	for m in list( conf["imp"] ):
		f.write( "import %s\n"%( m ) )
	f.write( """import problem

class my_problem( problem.template ):
	def __init__( self ):
		problem.template.__init__( self )
		self.mol = %s.molecule()
		self.mol.%s( "%s" )
"""%( conf["mol"]["module"], conf["mol"]["method"], conf["mol"]["fname"] ) )

	if( conf["box"] != None ):
		f.write( '\t\tself.mol.boxl = [ %.2lf, %.2lf, %.2lf ]\n'%( conf["box"][0], conf["box"][1], conf["box"][2] ) )

	f.write( """		f = open( "sele.pk", "rb" )
		self.sel = pickle.load( f )
		f.close()
		self.size = len( self.sel ) * 3
""" )
	if( sum( conf["sel"] ) == conf["MOL"].natm ):
		f.write( '\t\tself.coor = self.mol.coor\n' )
	else:
		f.write( """		self.coor = []
		for i in self.sel:
			self.coor += self.mol.coor[3*i:3*i+3][:]
""" )
	f.write( """		self.func = 0.0
		self.grad = []
		self.hess = None
		self.deco = {}
""" )

	# -- CONSTRUCTOR
	clr = ""
	who = 0
	for obj in conf["eng"]:
		if( obj["class"] == "x_namd.namd" ):
			f.write( """		self.e%02d = x_namd.namd()
		self.e%02d.exe = "%s"
		self.mol.psf_read( "%s" )
"""%( who, who, obj["exe"], obj["psf"] ) )
			input_namd( conf, obj )
			clr = '\t\tself.e%02d.update_charges( self.mol, "%s" )\n'%( who, obj["psf"] )
	
		if( obj["class"] == "x_namd.namd_pipe" ):
			f.write( """		try:
			os.unlink( "namd.pipe" )
		except Exception as e:
			print( e )
		os.mkfifo( "namd.pipe" )
		os.system( "%s" )
		time.sleep( 10 )
		self.e%02d = x_namd.namd_pipe()
		self.mol.psf_read( "%s" )
"""%( obj["exe"], who, obj["psf"] ) )
			input_namd_pipe( conf, obj )
			clr = '\t\tself.e%02d.update_charges( self.mol )\n'%( who )
	
		if( obj["class"] == "x_sander.Sander" ):
			f.write( """		self.e%02d = x_sander.Sander()
		self.e%02d.exe = "%s"
		self.mol.prmtop_read( "prmtop" )
"""%( who, who, obj["dir"] ) )
			input_sander( box, mol, obj )
			clr = '\t\tself.e%02d.update_charges( self.mol )\n'%( who )
	
		if( obj["class"] == "x_sander.py_sander" ):
			if( sum( conf["sqm"] ) == 0 ):
				f.write( '\t\tself.e%02d = x_sander.py_sander( self.mol, prmtop = "prmtop", cutoff = %.1lf, PBC = %s )\n'%( who, obj["cut"], obj["pbc"] ) )
			else:
				f.write( """		f = open( "sele_QM.pk", "rb" )
		s_qm = pickle.load( f )
		f.close()
		self.e%02d = x_sander.py_sander( self.mol, prmtop = "prmtop", cutoff = %.1lf,
			qmsel = s_qm, method = "%s", charge = %d, PBC = %s )
		self.mol.prmtop_read( "prmtop" )
"""%( who, float( obj["cut"] ), obj["qm_met"], int( obj["qm_chg"] ), obj["pbc"] ) )
	
		if( obj["class"] in [ "x_dftb.dftb", "x_sqm.sqm", "x_gamess.gamess", "x_nwchem.nwchem", "x_demon.demon",
			"x_orca.orca", "x_psi4.Psi4", "x_gaussian.gaussian", "x_gaussian.gaussian_MMEL",
			"x_tchem.tchem_sckt", "x_tchem.tchem" ] ): 
			f.write( """		if( self.mol.chrg == [] ):
			self.mol.chrg = [ 0.0 for i in range( self.mol.natm ) ]
        f = open( "sele_QM.pk", "rb" )
        s_qm = pickle.load( f )
        f.close()
        f = open( "sele_LA.pk", "rb" )
        s_la = pickle.load( f )
        f.close()
        f = open( "sele_MM.pk", "rb" )
        s_mm = pickle.load( f )
        f.close()
        for i in s_qm:
            self.mol.chrg[i] = 0.0
""" )
			if( clr != "" ):
				f.write( clr )
	
	
		if( obj["class"] == "x_dftb.dftb" ):
			f.write( """		self.e%02d = x_dftb.dftb( self.mol, s_qm, s_mm, s_la )
		self.e%02d.exe = "%s"
		self.e%02d.prm = "%s"
"""%( who, who, obj["exe"], who, obj["prm"] ) )
	
		if( obj["class"] == "x_sqm.sqm" ):
			f.write( """		self.e%02d = x_sqm.sqm( self.mol, \"\"\"%s, \"\"\", s_qm, s_mm, s_la )
		self.e%02d.exe = \"%s\"
"""%( who, obj["ini"], who, obj["dir"] ) )
	
		if( obj["class"] == "x_tchem.tchem_sckt" ):
			f.write( """		self.e%02d = x_tchem.tchem_sckt( self.mol, "%s", s_qm, s_mm, s_la )
"""%( who, obj["sck"] ) )
	
		if( obj["class"] in [ "x_gamess.gamess", "x_nwchem.nwchem", "x_demon.demon", "x_orca.orca", "x_tchem.tchem" ] ):
			f.write( """		f = open( "%s", "rt" )
		i_qm = f.read()
		f.close()
		self.e%02d = %s( self.mol, i_qm, s_qm, s_mm, s_la )
		self.e%02d.exe = "%s"
"""%( obj["ini"], who, obj["class"], who, obj["exe"] ) )
		if( obj["class"] == "x_dftb.dftb" ):
			f.write( """		self.e%02d = x_dftb.dftb( self.mol, s_qm, s_mm, s_la )
		self.e%02d.exe = "%s"
		self.e%02d.prm = "%s"
"""%( who, who, obj["exe"], who, obj["prm"] ) )
	
	
		if( obj["class"] in [ "x_gaussian.gaussian", "x_gaussian.gaussian_MMEL" ] ):
			f.write( """		f = open( "%s", "rt" )
		i_qm = f.read()
		f.close()
		f = open( "%s", "rt" )
		m_qm = f.read()
		f.close()
		f = open( "%s", "rt" )
		e_qm = f.read()
		f.close()
		self.e%02d = %s( self.mol, i_qm, m_qm, e_qm, s_qm, s_mm, s_la )
		self.e%02d.exe = "%s"
"""%( obj["ini"], obj["mid"], obj["end"], who, obj["class"], who, obj["exe"] ) )
	
		if( obj["class"] == "x_psi4.Psi4" ):
			f.write( """		o_qm = %s
		self.e%02d = x_psi4.Psi4( self.mol, s_qm, o_qm, s_mm, s_la )
"""%( obj["opt"], who ) )
	
		if( obj["class"] in [ "_qmmm.Int_QMLJ", "_qmmm.Int_QMLJ_MMEL" ] ):
			f.write( """		self.mol.%s( "%s" )
		exc = %s
		self.e%02d = _qmmm.Int_QMLJ( self.mol, s_qm, s_mm, exc )
"""%( obj["par"][0], obj["par"][1], obj["exc"], who ) )
	
		if( obj["class"] in [ "restraints.distance", "restraints.angle" ] ):
			f.write( '\t\tself.e%02d = %s( %s, %s, %s )\n'%( who, obj["class"], obj["kmb"], obj["ref"], obj["idx"] ) )
	
		if( obj["class"] == "restraints.multiple_distance" ):
			f.write( '\t\tself.e%02d = %s( %s, %s, %s, %s )\n'%( who, obj["class"],
				 obj["kmb"], obj["ref"], obj["idx"], obj["wei"] ) )
	
		who += 1
	
	for obj in conf["job"]:
		if( obj["class"] in [ "dynamics.velocity_verlet", "dynamics.langevin_verlet" ] ):
			if( sum( conf["sel"] ) == conf["MOL"].natm ):
				f.write( '\t\tself.mass = self.mol.mass\n' )
			else:
				f.write( """		self.mass = []
		for i in self.sel:
			self.mass.append( self.mol.mass[i] )
""" )
	
	f.write( '\t\tself.flog = open( "log", "wt" )\n' )
	if( conf["dcd"] != None ):
		f.write( """		self.dcd = mol_io.dcd()
		self.dcd.write( "%s", self.mol.natm, self.sel )
"""%( conf["dcd"]["fname"] ) )

	# -- STOP method
	f.write( """

	def stop( self ):
		self.flog.close()
""" )
	if( conf["dcd"] != None ):
		f.write( '\t\tself.dcd.close()\n' )
	who = 0
	for obj in conf["eng"]:
		if( obj["class"] in [ "x_namd.namd_pipe", "x_sander.py_sander" ] ):
			f.write( '\t\tself.e%02d.stop()\n'%( who ) )
		who += 1
	
	# -- LOG & CURRENT_STEP method
	f.write( """

	def log( self, str ):
		self.flog.write( str + "\\n" )
		self.flog.flush()
""" )
	if( conf["dcd"] != None ):
		f.write( """

	def current_step( self, step ):
		if( step %% %d == 0 ):
			self.mol.dcd_write( self.dcd )
"""%( conf["dcd"]["freq"] ) )
	else:
		f.write( """

	def current_step( self, step ):
		pass
""" )

	# -- MM updating method
	f.write( """

	def update_coor( self ):
		for i in range( self.size // 3 ):
			ii = 3*self.sel[i]
			jj = 3*i
			for j in [0, 1, 2]:
				self.mol.coor[ii+j] = self.coor[jj+j]
""" )

	# -- FUNC method
	f.write( """

	def get_func( self ):
		self.update_coor()
""" )
	who = 0
	for obj in conf["eng"]:
		if( obj["class"][0:6] != "_qmmm." ):
			if( obj["class"] in [ "x_namd.namd", "x_namd.namd_pipe", "x_sander.Sander", "x_sander.py_sander" ] ):
				f.write( '\t\tself.e%02d.update_coor( self.mol )\n'%( who ) )
			f.write( """		self.mol.func = 0.0
		self.e%02d.get_func( self.mol )
		self.deco["e%02d"] = self.mol.func
"""%( who, who ) )
		who += 1
	
	# -- GRAD method
	f.write( """		self.func = sum( self.deco.values() )


	def get_grad( self ):
		self.update_coor()
		self.mol.grad = [ 0.0 for i in range( 3 * self.mol.natm ) ]
""" )
	who = 0
	for obj in conf["eng"]:
		if( obj["class"] in [ "x_namd.namd", "x_namd.namd_pipe", "x_sander.Sander", "x_sander.py_sander" ] ):
			f.write( '\t\tself.e%02d.update_coor( self.mol )\n'%( who ) )
			f.write( """		self.mol.func = 0.0
		self.e%02d.get_grad( self.mol )
		self.deco["e%02d"] = self.mol.func
"""%( who, who ) )
		who += 1
	
	f.write( """		self.func = sum( self.deco.values() )
		self.grad = []
		for i in self.sel:
			self.grad += self.mol.grad[3*i:3*i+3][:]




obj = my_problem()
""" )

	# -- JOBS
	for obj in conf["job"]:
		if( obj["class"] in [ "minimize.steepest_descent", "minimize.fire", "minimize.l_bfgs" ] ):
			f.write( '%s( obj, step_number = %s, step_size = %s, print_frequency = %s, gradient_tolerance = %s, log_function = obj.log )\n' %(
				obj["class"], obj["stp"], obj["siz"], obj["prt"], obj["tol"] ) )
	
		if( obj["class"] == "minimize.conjugate_gradient_plus" ):
			f.write( '%s( obj, step_number = %s, print_frequency = %s, gradient_tolerance = %s, log_function = obj.log )\n' %(
				obj["class"], obj["stp"], obj["prt"], obj["tol"] ) )
	
		if( obj["class"] == "dynamics.velocity_verlet" ):
			f.write( 'dynamics.assign_velocities( obj, temperature = %s )\n' %( obj["tmp"] ) )
			f.write( '%s( obj, step_number = %s, step_size = %s, print_frequency = %s, scale_frequency = %s, temperature = %s, log_function = obj.log )\n' %(
				obj["class"], obj["stp"], obj["siz"], obj["prt"], obj["scl"], obj["tmp"] ) )
	
		if( obj["class"] == "dynamics.langevin_verlet" ):
			f.write( 'dynamics.assign_velocities( obj, temperature = %s )\n' %( obj["tmp"] ) )
			f.write( '%s( obj, step_number = %s, step_size = %s, print_frequency = %s, gamma_factor = %s, temperature = %s, log_function = obj.log )\n' %(
				obj["class"], obj["stp"], obj["siz"], obj["prt"], obj["gam"], obj["tmp"] ) )

	# -- END
	if( conf["out"] != None ):
		if( conf["out"]["method"] == "pickled_molecule" ):
			f.write( """f = open( "%s", "wb" )
pickle.dump( obj.mol, f )
f.close()
"""%( conf["out"]["fname"] ) )
		else:
			f.write( 'obj.mol.%s( "%s" )\n'%( conf["out"]["method"], conf["out"]["fname"] ) )
	f.write( 'obj.stop()\n' )
	f.close()



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
				Hint( l, h, background = tbg )
				e = tkinter.Entry( f, font = self.__fnt, justify = tkinter.LEFT )
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
				Hint( l, h, background = tbg )
				e = tkinter.Entry( f, font = self.__fnt, justify = tkinter.LEFT )
				e.insert( 0, d )
				e.pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )
				self.__job[key].append( ( p, e ) )


	def __mkinput( self ):
#		f = open( "helper.inp", "wt" )
#		if( self.__wid["box_x"].get() != "" ):
#			f.write( "box " + self.__wid["box_x"].get() )
#			if( self.__wid["box_y"].get() != "" ):
#				f.write( " " + self.__wid["box_y"].get() )
#			else:
#				f.write( " " + self.__wid["box_x"].get() )
#			if( self.__wid["box_z"].get() != "" ):
#				f.write( " " + self.__wid["box_z"].get() )
#			else:
#				f.write( " " + self.__wid["box_x"].get() )
#			f.write( "\n" )
#		f.write( "mol " + self.__wid["mol_cb"].get() + " " + self.__wid["mol"].get() + "\n" )
#		tmp = self.__wid["sel"].get( 1.0, tkinter.END ).strip().replace( "\n", " " )
#		if( tmp != "" ):
#			f.write( "sel " + tmp + "\n" )
#		tmp = self.__wid["sqm"].get( 1.0, tkinter.END ).strip().replace( "\n", " " )
#		if( tmp != "" ):
#			f.write( "sqm " + tmp + "\n" )
#		tmp = self.__wid["lnk"].get( 1.0, tkinter.END ).strip().replace( "\n", " " )
#		if( tmp != "" ):
#			f.write( "lnk " + tmp + "\n" )
#		tmp = self.__wid["nbn"].get( 1.0, tkinter.END ).strip().replace( "\n", " " )
#		if( tmp != "" ):
#			f.write( "nbn " + tmp + "\n" )
#		if( self.__wid["dcd_s"].get() != "" ):
#			f.write( "dcd " + self.__wid["dcd_s"].get() + " " + self.__wid["dcd_n"].get() + "\n" )
#		f.write( "out " + self.__wid["out_cb"].get() + " " + self.__wid["out"].get() + "\n" )
#		for key in iter( self.__eng ):
#			f.write( "eng " + key + "\n" )
#			for p,e in self.__eng[key][1:]:
#				f.write( "\t" + p + " " + e.get() + "\n" )
#			f.write( "\n" )
#		for key in iter( self.__job ):
#			f.write( "job " + key + "\n" )
#			for p,e in self.__job[key][1:]:
#				f.write( "\t" + p + " " + e.get() + "\n" )
#			f.write( "\n" )
#		f.close()
		print( ">> 'run' script generated!" )
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
		self.__box = None
		self.__sel = [ False for i in range( self.__mol.natm ) ]
		self.__sqm = self.__sel[:]
		self.__nbn = self.__sel[:]
		self.__lnk = []
		self.__mod = []
		# -------------------------------------------------

		self.app = tkinter.Tk()
		self.app.title( "-- Helper --" )
		self.app.maxsize( width = 540, height = 2048 )
		self.app.minsize( width = 540, height = 480 )

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
		Hint( l01, "Fill only the first one (X) for a CUBIC symmetry" )
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
		Hint( l05, """Selection of Link Atoms, QM MM pairs list:

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

#		f07 = tkinter.Frame( self.frm, bg = self.__col )
#		f07.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH, padx = 4, pady = 4 )
#		l07 = tkinter.Label( f07, text = " dcd ", font = self.__fnt, bg = self.__col, justify = tkinter.CENTER )
#		l07.pack( side = tkinter.LEFT, expand = 0, fill = None )
#		Hint( l07, "Frequency for writing (>0) the trajectory file ('dcd')" )
#		self.__wid["dcd"] = tkinter.Entry( f07, font = self.__fnt, justify = tkinter.LEFT, width = 10 )
#		self.__wid["dcd"].insert( 0, "0" )
#		self.__wid["dcd"].pack( side = tkinter.LEFT, expand = 0, fill = tkinter.X )

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

		bbb = tkinter.Button( self.frm, text = " build! ", font = self.__fnt, justify = tkinter.CENTER, width = 4, command = self.__mkinput )
		bbb.pack( side = tkinter.LEFT, expand = 1, fill = tkinter.X )

		self.app.mainloop()






gui_application( sys.argv[1] )

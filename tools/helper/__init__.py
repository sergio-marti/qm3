# -*- coding: iso-8859-1 -*-
from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

import	bottle
import	os
import	re
import	config

import	qm3.mol
import	qm3.elements
import	qm3.utils
if( "NAMD" in config.MM_engines ):
	import	qm3.engines.namd


try:
	cwd = os.environ["HELPR_CWD"]
except:
	cwd = os.dirname( sys.argv[0] )
try:
	uwd = os.environ["HELPR_UWD"]
except:
	uwd = os.environ["HOME"]
try:
	pdb = os.path.abspath( os.environ["HELPR_COR"] )
	top = os.path.abspath( os.environ["HELPR_TOP"] )
except:
	pdb = "None"
	top = "None"


@bottle.route( "/" )
def __index():
	qm = ""
	for ee in [ j for i,j in sorted( [ k.lower(), k ] for k in config.QM_engines ) ]:
		qm += "<option>%s</option>"%( ee )
	mm = ""
	for ee in [ j for i,j in sorted( [ k.lower(), k ] for k in config.MM_engines ) ]:
		mm += "<option>%s</option>"%( ee )
	pdb_view = "None"
	if( pdb != "None" ):
		pdb_view = "<a href=\"/loadpdb\" target=\"_blank\">%s</a>"%( os.path.basename( pdb ) )
	qm_sel = ""
	qm_cor = "none"
	qm_chg = "block"
	qm_env = "block"
	if( top == "None" and pdb == "None" ):
		qm_sel = "*"
		qm_cor = "block"
		qm_chg = "none"
		qm_env = "none"
	return( bottle.template( "index.html", 
		pdb = pdb_view,
		top = os.path.basename( top ),
		qm_eng = qm,
		qm_sel = qm_sel,
		qm_chg = qm_chg,
		qm_cor = qm_cor,
		qm_env = qm_env,
		mm_eng = mm ) )


@bottle.route( "/jsmol/<sub:path>" )
def __jsmol( sub ):
	if( sub[-3:] == "png" ):
		return( bottle.static_file( sub, root = "jsmol", mimetype = "image/png" ) )
	else:
		return( bottle.static_file( sub, root = "jsmol" ) )


@bottle.route( "/editor" )
def __editor():
	return( bottle.template( "editor.html" ) )


@bottle.route( "/loadpdb" )
def __loadpdb():
	return( bottle.template( "loadpdb.html" ) )


@bottle.route( "/flushpdb" )
def __flushpdb():
	if( os.path.isfile( pdb ) ):
		job = os.path.basename( pdb )
		fld = os.path.dirname( pdb )
		return( bottle.static_file( job, root = fld, mimetype = "application/octet-stream" ) )


SP0 = re.compile( "^([0-9]+)$" )
SP1 = re.compile( "^([0-9]+)-([0-9]+)$" )
SP2 = re.compile( "^([A-Z0-9]+):([0-9]+)$" )
SP3 = re.compile( "^([A-Z0-9]+):([0-9]+)-([0-9]+)$" )
SP4 = re.compile( "^([A-Z0-9]+)/([0-9]+)/(.+)$" )
SP5 = re.compile( "^([A-Z0-9]+):([0-9]+)@([0-9\.]+)$" )
def __selection( mol, buf ):
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
		# -- atom number (NOT C-indexing)
		elif( SP0.match( itm ) ):
			sel[int(itm)-1] = True
		# -- range of atom numbers (NOT C-indexing)
		elif( SP1.match( itm ) ):
			rng = [ int( i ) for i in SP1.findall( itm )[0] ]
			for i in range( rng[0], rng[1] + 1 ):
				sel[i-1] = True
		# -- residue number (by chain)
		elif( SP2.match( itm ) ):
			rng = SP2.findall( itm )[0]
			for i in list( mol.indx[rng[0]][int(rng[1])].values() ):
				sel[i] = True
		# -- range of residue numbers (by chain)
		elif( SP3.match( itm ) ):
			rng = SP3.findall( itm )[0]
			for i in range( int( rng[1] ), int( rng[2] ) + 1 ):
				for j in list( mol.indx[rng[0]][i].values() ):
					sel[j] = True
		# -- chain / residue_number / atom_label
		elif( SP4.match( itm ) ):
			rng = SP4.findall( itm )[0]
			sel[mol.indx[rng[0]][int(rng[1])][rng[2]]] = True
		# -- radial selection by residue around chain / residue_number
		elif( SP5.match( itm ) ):
			rng = SP5.findall( itm )[0]
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


@bottle.route( "/mkinput", method = "POST" )
def __mkinput():
	job = os.path.join( uwd, "run_qm3" )
	try:
		mol = qm3.mol.molecule( pdb )
	except:
		mol = None
	sel = []
	# -- selections: QMsel, QMenv, Selec (possible link_atoms and exclussions...)
	_sqm = bottle.request.forms.get( "QMsel" ).strip()
	_eqm = bottle.request.forms.get( "QMenv" ).strip()
	_sel = bottle.request.forms.get( "Selec" ).strip()

	"""
casos posibles:

_sel == "" || _sel == "*":
	entra todo en la definición del sistema: sel = []
	
_sel != "" && _sel != "*":

	mol == None:
		solo una parte de los átomos QM es activa de cara a las acciones
		veo dos posibilidades:
			1) restringir las posibles selecciones: *, int, int-int
			2) cargar las coordenadas cartesianas en la molécula con un pickle
				y aplicar parseo de selecciones igualmente...

	mol != None:
		hay mezcla de QM y MM: aplicar parseo de selecciones
		es responsabilidad del usuario congelar átomos en el engine MM
	"""

	# -- let's go!
	f = open( job, "wt" )
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
try:
	import cStringIO
except:
	import io as cStringIO

import	qm3.mol
import	qm3.problem
""" )
	# -- imports: QMeng, MMeng, Todo[min_fire pes_scan min_quad nor_mode mep_path dyn_lang]
	# -- molecule definition (based on pdb or cartesian input provided: QMcor)
	f.write( """

class my_problem( qm3.problem.template ):
	def __init__( self ):
		qm3.problem.template.__init__( self )
""" )
	xyz = bottle.request.forms.get( "QMcor" ).strip()
	if( mol != None and top != "None" ):
		f.write( """		self.mole = qm3.mol.molecule( "%s" )\n"""%( pdb ) )
		if( top != "None" ):
			f.write( """		self.mole.psf_read( "%s" )
		self.mole.nbnd_read( "@@@@" )
"""%( top ) )
	else:
		f.write( """		f = cStringIO.StringIO( \"\"\"%d\n\n%s\"\"\" )
		f.seek( 0 )
		self.mole = qm3.mol.molecule()
		self.mole.xyz_read( f )
		f.close()
		self.mole.guess_atomic_numbers()
		self.mole.fill_masses()
		"""%( len( xyz.split() ) // 4, xyz.replace( "\r", "" ) ) )
	# -- boxl: box_X, box_Y, box_Z
	tmp = bottle.request.forms.get( "box_X" )
	try:
		box = [ float( tmp ), float( tmp ), float( tmp ) ]
	except:
		box = []
	try:
		tmp = float( bottle.request.forms.get( "box_Y" ) )
		box[1] = tmp
	except:
		pass
	try:
		tmp = float( bottle.request.forms.get( "box_Z" ) )
		box[2] = tmp
	except:
		pass
	if( box != [] ):
		f.write( """		self.mole.boxl = [ %.2lf, %.2lf, %.2lf ]\n"""%( box[0], box[1], box[2] ) )
	# -- definition of the active system (full or sele.pk based)
	if( sel == [] ):
		f.write( """
		self.size = 3 * self.mole.natm
		self.coor = self.mole.coor
		self.mass = self.mole.mass
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
	# -- generic
	f.write( """
		self.func = 0.0
		self.grad = []
		self.hess = []
""" )
	# -- engines: QMeng (QMchg QMmul QMmet QMbas), MMeng, Other (like restraints...)
# >> todos los restraints en una lista, de modo que los llama en bucle si hay alguno...
# >> parsear y definir lo que sea
	# -- zero QM charges on the MM engine (if any...)
	# -- update_coor method
	f.write( """

	def update_coor( self ):
		if( self.sele != [] ):
			for i in range( len( self.sele ) ):
				ii = 3 * self.sele[i]
				jj = 3 * i
				for j in [ 0, 1, 2 ]:
					self.mole.coor[ii+j] = self.coor[jj+j]
""" )
	# -- get_func method
	f.write( """

	def get_func( self ):
		self.update_coor()
		self.mole.func = 0.0
""" )
	# -- get_grad method
	f.write( """

	def get_grad( self ):
		self.update_coor()
		self.mole.func = 0.0
		self.mole.grad = [ 0.0 for i in range( 3 * self.mole.natm ) ]
""" )
	# -- get_hess method
	f.write( """

	def get_hess( self ):
		self.update_coor()
		self.mole.func = 0.0
		self.mole.grad = [ 0.0 for i in range( 3 * self.mole.natm ) ]
""" )
	# -- object generation
	f.write( """




obj = my_problem()

""" )
	# -- actions
	# -- stop engines (if needed...)
	f.close()
	return( "<h1>Check the script generated:</h1><h1>%s</h1>"%( job ) )


bottle.run( host = "localhost", port = 8080 )

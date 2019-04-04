# -*- coding: iso-8859-1 -*-
from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	bottle
import	os
import	config


try:
	cwd = os.environ["HELPR_CWD"]
except:
	cwd = ""
try:
	uwd = os.environ["HELPR_UWD"]
except:
	uwd = ""
try:
	pdb = os.path.abspath( os.environ["HELPR_COR"] )
except:
	pdb = "None"
try:
	top = os.path.abspath( os.environ["HELPR_TOP"] )
except:
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
	return( bottle.template( "index.html", 
		pdb = pdb_view,
		top = os.path.basename( top ),
		qm_eng = qm,
		mm_eng = mm ) )


@bottle.route( "/jsmol/<sub:path>" )
def __jsmol( sub ):
	if( sub[-3:] == "png" ):
		return( bottle.static_file( sub, root = "jsmol", mimetype = "image/png" ) )
	else:
		f = open( os.path.join( cwd, "jsmol", sub ), "rt" )
		r = f.read()
		f.close()
		return( r )


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


@bottle.route( "/mkinput", method = "POST" )
def __mkinput():
	print( list( bottle.request.forms.allitems() ) )
	return( "<h1>Input processed!</h1>" )


bottle.run( host = "localhost", port = 8080 )

# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

try:
	import	urllib2 as ulib
except:
	import	urllib.request as ulib

import	qm3.mol
import	qm3.io


###################################################################################################
# SDF
#
def sdf_read( fname = None ):
	mol = qm3.mol.molecule()
	f = qm3.io.open_r( fname )
	for i in range( 4 ):
		l = f.readline()
	mol.natm = int( l.strip().split()[0] )
	for i in range( mol.natm ):
		t = f.readline().split()
		mol.labl.append( t[3] )
		mol.resi.append( 1 )
		mol.resn.append( "XXX" )
		mol.segn.append( "X" )
		mol.coor += [ float( j ) for j in t[0:3] ]
	qm3.io.close( f, fname )
	mol.settle()
	return( mol )


def db_download( code ):
	svr, mid = code.split( ":" )
	if( svr == "pubchem" or svr == "" ):
		r = ulib.Request( "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/record/SDF/?record_type=3d&response_type=save"%( mid ) )
	elif( svr == "chemspider" ):
		r = ulib.Request( "http://www.chemspider.com/FilesHandler.ashx?type=str&3d=yes&id=" + mid )
#	r = ulib.Request( "https://cactus.nci.nih.gov/chemical/structure/%s/file?format=sdf&get3d=True"%( mol_id ) )
	f = ulib.urlopen( r )
	try:
		mol = sdf_read( f )
	except:
		print( "No 3D-structure found..." )
	f.close()
	return( mol )


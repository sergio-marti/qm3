# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

import	qm3.mol
import	qm3.io


###################################################################################################
# MOL2: obabel -i @@@ @@@ -o mol2
#
def mol2_read( fname = None ):
	mol = qm3.mol.molecule()
	f = qm3.io.open_r( fname )
	l = f.readline()
	while( l ):
		if( l.strip() == "@<TRIPOS>MOLECULE" ):
			f.readline()
			t = f.readline().strip().split()
			mol.natm = int( t[0] )
			nbnd = int( t[1] )
			bond = []
		if( l.strip() == "@<TRIPOS>ATOM" ):
			for i in range( mol.natm ):
				t = f.readline().strip().split()
				mol.labl.append( t[1] )
				mol.coor += [ float( t[2] ), float( t[3] ), float( t[4] ) ]
				mol.type.append( t[5] )
				mol.resi.append( int( t[6] ) )
				mol.resn.append( t[7][0:3] )
				mol.chrg.append( float( t[8] ) )
				mol.segn.append( "X" )
		if( l.strip() == "@<TRIPOS>BOND" ):
			for i in range( nbnd ):
				t  = f.readline().split()
				ii = int( t[1] ) - 1
				jj = int( t[2] ) - 1
				bond.append( [ ii, jj ] )
		l = f.readline()
	qm3.io.close( f, fname )
	mol.settle()
	return( mol, bond )

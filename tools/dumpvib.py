#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

import	struct
import	qm3.mol
import	qm3.utils
import	qm3.elements
import	qm3.maths.matrix
import	os, stat

if( len( sys.argv ) == 1 ):
	print( "%s NMODE_int SMB DMP [AMP=10.0]"%( sys.argv[0] ) )
	sys.exit( 1 )

sele = int( sys.argv[1] )

f = file( sys.argv[2], "rt" )
symb = [ i.title() for i in f.read().split() ]
f.close()
#m = qm3.mol.molecule()
#m.xyz_read( sys.argv[2] )
#symb = m.labl

mass = [ qm3.elements.mass[qm3.elements.rsymbol[i]] for i in symb ]
print( zip( symb, mass ) )
size = len( symb ) * 3
flen = os.stat( sys.argv[3] )[stat.ST_SIZE]
if( flen == size*(2+size)*8 ):
	f = file( sys.argv[3], "rb" )
	coor = struct.unpack( "%dd"%( size ), f.read( size * 8 )  )
	f.read( size * 8 )
	hess = struct.unpack( "%dd"%( size * size ), f.read( size * size * 8 )  )
	f.close()
# ------------------------------------------------------------
#	f = open( "hess", "wb" )
#	for i in range( size * size ):
#		f.write( struct.pack( "d", hess[i] ) )
#	f.close()
# ------------------------------------------------------------
	val, vec = qm3.utils.hessian_frequencies( mass, coor, hess, True )
	print( val[0:min(10,max(10,sele+1))] )
	ampl = 10.
	if( len( sys.argv ) == 5 ):
		ampl = float( sys.argv[4] )
	qm3.utils.normal_mode_view( coor, val, vec, symb, sele, afac = ampl )

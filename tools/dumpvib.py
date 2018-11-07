#!/usr/bin/env python

import	struct
import	qm3.utils
import	qm3.elements
import	qm3.maths.matrix
import	sys
import	os, stat


f = file( "symbols", "rt" )
symb = f.read().upper().split()
f.close()
mass = [ qm3.elements.mass[qm3.elements.rsymbol[i]] for i in symb ]
size = len( symb ) * 3
flen = os.stat( "update.dump" )[stat.ST_SIZE]
coor = []
hess = []
f = file( "update.dump", "rb" )
if( flen == 2*size*8+4*size*(size+1) ):
	print ">> C format"
	coor = struct.unpack( "%dd"%( size / 3 ), f.read( size * 8 )  )
	f.read( size * 8 )
	hess = f.read()
elif( flen == 2*size*8+4*size*(size+1) + 24 ):
	print ">> Fortran format"
	f.read( 4 )
	coor = struct.unpack( "%dd"%( size ), f.read( size * 8 )  )
	f.read( 12 + size * 8 )
	n = struct.unpack( "i", f.read(4) )[0]
	hess = struct.unpack( "%dd"%( n / 8 ), f.read( n ) )
f.close()
if( coor and hess ):
	hess = qm3.maths.matrix.swap_upper_diagonal_to_rows( size, hess )
	val, vec = qm3.utils.hessian_frequencies( mass, coor, hess, True )
	print val[0:10]
	fc = 10.
	if( len( sys.argv ) == 3 ):
		fc = float( sys.argv[2] )
	qm3.utils.normal_mode_view( coor, val, vec, symb, int( sys.argv[1] ), afac = fc )

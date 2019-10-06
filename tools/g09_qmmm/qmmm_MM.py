#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os


if( len( sys.argv ) != 4 ):
	print( "%s selection pdb psf"%( sys.argv[0] ) )
	sys.exit( 1 )
if( not os.access( "qmmm.fchk", os.R_OK ) ):
	print( "- 'qmmm.fchk' not found!" )
	sys.exit( 1 )


_ATM = 0.52917721067
coor = []
chrg = []
qmat = []


# -- load selections
f = open( sys.argv[1], "rt" )
for i in range( int( f.readline().strip() ) ):
	qmat.append( int( f.readline().strip() ) )
f.close()
qmat.sort()


# -- load coordinates and charges from qmmm.fchk
f = open( "qmmm.fchk", "rt" )
l = f.readline()
while( l != "" ):
	if( l[0:29] == "Current cartesian coordinates" ):
		i = int( l.strip().split()[-1] )
		j = int( i // 5 ) + ( i % 5 != 0 )
		i = 0
		while( i < j ):
			coor += [ float( k ) * _ATM for k in f.readline().strip().split() ]
			i += 1
	elif( l[0:11] == "ESP Charges" ):
		i = int( l.strip().split()[-1] )
		j = int( i // 5 ) + ( i % 5 != 0 )
		i = 0
		while( i < j ):
			chrg += [ float( k ) for k in f.readline().strip().split() ]
			i += 1
	l = f.readline()
f.close()


# -- modify PDB file
f = open( sys.argv[2], "rt" )
g = open( "__tmp__", "wt" )
i = 0
s = len( qmat )
w = 0
for l in f:
	if( l[0:4] == "ATOM" or l[0:4] == "HETA" ):
		if( w < s and i == qmat[w] ):
			j = w * 3
			g.write( l[0:30] + "%8.3lf%8.3lf%8.3lf"%( coor[j], coor[j+1], coor[j+2] ) + l[54:] )
			w += 1
		else:
			g.write( l )
		i += 1
	else:
		g.write( l )
g.close()
f.close()
os.rename( "__tmp__", sys.argv[2] )


# -- modify PSF
f = open( sys.argv[3], "rt" )
g = open( "__tmp__", "wt" )
g.write( f.readline() ); g.write( f.readline() )
l = f.readline(); g.write( l )
for i in range( int( l.strip().split()[0] ) + 1 ):
	g.write( f.readline() )
l = f.readline(); g.write( l )
n = int( l.strip().split()[0] )
s = len( qmat )
w = 0
for i in range( n ):
	l = f.readline();
	if( w < s and i == qmat[w] ):
		g.write( l[0:34] + "%10.6lf"%( chrg[w] ) + l[44:] )
		w += 1
	else:
		g.write( l )
g.write( f.read() )
g.close()
f.close()
os.rename( "__tmp__", sys.argv[3] )

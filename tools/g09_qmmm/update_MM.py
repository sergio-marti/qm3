#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os


for fn in [ "qmmm.sel", "qmmm.crd", "qmmm.psf", "qmmm.fchk" ]:
	if( not os.path.isfile( fn ) ):
		print( "-- missing: '%s'"%( fn ) )
		sys.exit( 1 )


_ATM = 0.52917721067
coor = []
chrg = []
qmat = []


# -- load selections
f = open( "qmmm.sel", "rt" )
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
f = open( "qmmm.crd", "rt" )
g = open( "__qmmm__", "wt" )
l = f.readline()
while( l[0] == "*" ):
	g.write( l )
	l = f.readline()
g.write( l )
w = 0
for i in range( int( l.strip() ) ):
	l = f.readline()
	if( i in qmat ):
		j = w * 3
#0123456789.123456789.123456789.123456789.123456789.
# 4338 1442 HOH  OH2    2.19700  18.51600 -16.13200 W    1441   0.00000
		g.write( l[0:20] + "%10.5lf%10.5lf%10.5lf"%( coor[j], coor[j+1], coor[j+2] ) + l[50:] )
		w += 1
	else:
		g.write( l )
g.close()
f.close()
os.rename( "__qmmm__", "qmmm.crd" )


# -- modify PSF
f = open( "qmmm.psf", "rt" )
g = open( "__qmmm__", "wt" )
g.write( f.readline() ); g.write( f.readline() )
l = f.readline(); g.write( l )
for i in range( int( l.strip().split()[0] ) + 1 ):
	g.write( f.readline() )
l = f.readline(); g.write( l )
for i in range( int( l.strip().split()[0] ) ):
	l = f.readline();
	chrg.append( float( l[34:44] ) )
	if( i in qmat ):
#0123456789.123456789.123456789.123456789.1234567890.
#       1 A    1    COP  C1   C2    -0.157000       12.0100           0   0.00000     -0.301140E-02
		g.write( l[0:34] + "%10.6lf"%( chrg[w] ) + "    " + l[48:] )
		w += 1
	else:
		g.write( l )
g.write( f.read() )
g.close()
f.close()
os.rename( "__qmmm__", "qmmm.psf" )

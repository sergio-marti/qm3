#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import  re

if( len( sys.argv ) == 1 ):
	print( "%s deps top prm" )
	sys.exit(1)

f = open( sys.argv[1], "rt" )
RES = f.readline().strip()
BND = []
ANG = []
DIE = []
for l in f:
	t = [ i.strip() for i in l.strip().split( "-" ) ]
	if( len( t ) == 2 ):
		BND.append( t[:] )
	elif( len( t ) == 3 ):
		ANG.append( t[:] )
	else:
		DIE.append( t[:] )
f.close()

TBL = {}
f = open( sys.argv[2], "rt" )
l = f.readline()
while( l != "" ):
	if( l[0:4] == "RESI" and l.find( RES ) > -1 ):
		t = l.strip().split()
		while( len( t ) != 0 ):
			if( t[0] == "ATOM" ):
				TBL[t[1]] = t[2]
			t = f.readline().strip().split()
	l = f.readline()
f.close()
 
f = open( sys.argv[3], "rt" )
b = f.readlines()
f.close()

for ai, aj in BND:
	pij = re.compile( "^%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( TBL[ai], TBL[aj] ) )
	pji = re.compile( "^%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( TBL[aj], TBL[ai] ) )
	for l in b:
		t = []
		if( pij.match( l ) ):
			t = l.strip().split()
		elif( pji.match( l ) ):
			t = l.strip().split()
		if( t ):
			print( ai, aj, float( t[2] ) * 4.184, float( t[3] ) )

for ai, aj, ak in ANG:
	pijk = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( TBL[ai], TBL[aj], TBL[ak] ) )
	pkji = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( TBL[ak], TBL[aj], TBL[ai] ) )
	for l in b:
		t = []
		if( pijk.match( l ) ):
			t = l.strip().split()
		elif( pkji.match( l ) ):
			t = l.strip().split()
		if( t ):
			print( ai, aj, ak, float( t[3] ) * 4.184, float( t[4] ) )

for ai, aj, ak, al in DIE:
	pijkl = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( TBL[ai], TBL[aj], TBL[ak], TBL[al] ) )
	plkji = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( TBL[al], TBL[ak], TBL[aj], TBL[ai] ) )
	pxjkx = re.compile( "^X[\ ]+%s[\ ]+%s[\ ]+X[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( TBL[aj], TBL[ak] ) )
	pxkjx = re.compile( "^X[\ ]+%s[\ ]+%s[\ ]+X[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( TBL[ak], TBL[aj] ) )
	for l in b:
		t = []
		if( pijkl.match( l ) ):
			t = l.strip().split()
		elif( plkji.match( l ) ):
			t = l.strip().split()
		elif( pxjkx.match( l ) ):
			t = l.strip().split()
		elif( pxkjx.match( l ) ):
			t = l.strip().split()
		if( t ):
			print( ai, aj, ak, al, float( t[4] ) * 4.184, int( t[5] ), float( t[6] ) )



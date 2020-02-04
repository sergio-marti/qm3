#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import  re

if( len( sys.argv ) == 1 ):
	print( "%s prm          [>> works on exclusions.src]"%( sys.argv[0] ) )
	sys.exit(1)

f = open( "exclusions.src", "rt" )
BND = []
ANG = []
DIE = []
for l in f:
	if( l.find( "||" ) > -1 ):
		t = [ i.strip() for i in l.strip().split( "||" )[1].split( "-" ) ]
		if( len( t ) == 2 ):
		    BND.append( t[:] )
		elif( len( t ) == 3 ):
		    ANG.append( t[:] )
		else:
		    DIE.append( t[:] )
f.close()

f = open( sys.argv[1], "rt" )
b = f.readlines()
f.close()

for ai, aj in BND:
	pij = re.compile( "^%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( ai, aj ) )
	pji = re.compile( "^%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( aj, ai ) )
	for l in b:
	    t = []
	    if( pij.match( l ) ):
	        t = l.strip().split()
	    elif( pji.match( l ) ):
	        t = l.strip().split()
	    if( t ):
	        print( ai, aj, float( t[2] ) * 4.184 * 2.0, float( t[3] ) )

for ai, aj, ak in ANG:
	pijk = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( ai, aj, ak ) )
	pkji = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( ak, aj, ai ) )
	for l in b:
	    t = []
	    if( pijk.match( l ) ):
	        t = l.strip().split()
	    elif( pkji.match( l ) ):
	        t = l.strip().split()
	    if( t ):
	        print( ai, aj, ak, float( t[3] ) * 4.184 * 2.0, float( t[4] ) )

for ai, aj, ak, al in DIE:
	pijkl = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+%s[\ ]+([0-9\-\.]+)[\ ]+([1-6])[\ ]+([0-9\.]+)"%( ai, aj, ak, al ) )
	plkji = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+%s[\ ]+([0-9\-\.]+)[\ ]+([1-6])[\ ]+([0-9\.]+)"%( al, ak, aj, ai ) )
	pxjkx = re.compile( "^X[\ ]+%s[\ ]+%s[\ ]+X[\ ]+([0-9\-\.]+)[\ ]+([1-6])[\ ]+([0-9\.]+)"%( aj, ak ) )
	pxkjx = re.compile( "^X[\ ]+%s[\ ]+%s[\ ]+X[\ ]+([0-9\-\.]+)[\ ]+([1-6])[\ ]+([0-9\.]+)"%( ak, aj ) )
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
	        if( float( t[4] ) != 0.0 ):
	            print( ai, aj, ak, al, float( t[4] ) * 4.184 * 2.0, int( t[5] ), float( t[6] ) )



#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import re

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

_BND = {}
for ai, aj in BND:
    pij = re.compile( "^%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( ai, aj ) )
    pji = re.compile( "^%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( aj, ai ) )
    for l in b:
        L = l.upper()
        t = []
        if( pij.match( L ) ):
            t = L.strip().split()
        elif( pji.match( L ) ):
            t = L.strip().split()
        if( t ):
            _BND["%s-%s"%( ai, aj )] = [ float( t[2] ) * 4.184 * 2.0, float( t[3] ) ]

_ANG = {}
for ai, aj, ak in ANG:
    pijk = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( ai, aj, ak ) )
    pkji = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)"%( ak, aj, ai ) )
    for l in b:
        L = l.upper()
        t = []
        if( pijk.match( L ) ):
            t = L.strip().split()
        elif( pkji.match( L ) ):
            t = L.strip().split()
        if( t ):
            _ANG["%s-%s-%s"%( ai, aj, ak )] = [ float( t[3] ) * 4.184 * 2.0, float( t[4] ) ]

_DIE = {}
for ai, aj, ak, al in DIE:
    pijkl = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+%s[\ ]+([0-9\-\.]+)[\ ]+([1-6])[\ ]+([0-9\.]+)"%( ai, aj, ak, al ) )
    plkji = re.compile( "^%s[\ ]+%s[\ ]+%s[\ ]+%s[\ ]+([0-9\-\.]+)[\ ]+([1-6])[\ ]+([0-9\.]+)"%( al, ak, aj, ai ) )
    pxjkx = re.compile( "^X[\ ]+%s[\ ]+%s[\ ]+X[\ ]+([0-9\-\.]+)[\ ]+([1-6])[\ ]+([0-9\.]+)"%( aj, ak ) )
    pxkjx = re.compile( "^X[\ ]+%s[\ ]+%s[\ ]+X[\ ]+([0-9\-\.]+)[\ ]+([1-6])[\ ]+([0-9\.]+)"%( ak, aj ) )
    for l in b:
        L = l.upper()
        t = []
        if( pijkl.match( L ) ):
            t = L.strip().split()
        elif( plkji.match( L ) ):
            t = L.strip().split()
        elif( pxjkx.match( L ) ):
            t = L.strip().split()
        elif( pxkjx.match( L ) ):
            t = L.strip().split()
        if( t ):
            if( float( t[4] ) != 0.0 ):
                _DIE["%s-%s-%s-%s"%( ai, aj, ak, al )] = [ float( t[4] ) * 4.184 * 2.0, int( t[5] ), float( t[6] ) ]


f = open( "exclusions.src", "rt" )
g = open( "exclusions.py", "wt" )
l = f.readline()
while( l != "" ):
    if( l.find( "||" ) > -1 ):
        g.write( l )
        t = [ i.strip() for i in l.strip().split( "||" )[1].split( "-" ) ]
        k = "-".join( t )
        if( len( t ) == 2 ):
            if( k in _BND ):
                t = f.readline().split()
                g.write( "        %s\n"%( " ".join( t[0:2] + [ "%.1lf,"%( _BND[k][0] ), "%.3lf,"%( _BND[k][1] ) ] + t[4:] ) ) )
            else:
                f.readline(); f.readline(); f.readline(); f.readline()
        elif( len( t ) == 3 ):
            if( k in _ANG ):
                t = f.readline().split()
                g.write( "        %s\n"%( " ".join( t[0:2] + [ "%.1lf,"%( _ANG[k][0] ), "%.1lf,"%( _ANG[k][1] ) ] + t[4:] ) ) )
            else:
                f.readline(); f.readline(); f.readline(); f.readline()
        else:
            if( k in _DIE ):
                t = f.readline().split()
                g.write( "        %s\n"%( " ".join( t[0:3] + [ "%d: ["%( _DIE[k][1] ), "%.3lf, %.1lf ],"%( _DIE[k][0], _DIE[k][2] ) ] + t[8:] ) ) )
            else:
                f.readline(); f.readline(); f.readline(); f.readline()
    else:
        g.write( l )
    l = f.readline()
g.close()
f.close()

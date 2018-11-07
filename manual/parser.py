#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

f = open( sys.argv[1], "rt" )
src = f.readlines()
f.close()

for i in range( len( src ) - 1, -1, -1 ):
	if( src[i][0:8] == "#SOURCE@" ):
		who = src[i][8:].strip()
		del src[i]
		f = open( who, "rt" )
		tmp = f.readlines()
		f.close()
		for j in range( len( tmp ) - 1, -1, -1 ):
			src.insert( i, tmp[j] )

f = open( sys.argv[1][4:] + ".tex", "wt" )
f.write( "".join( src ) )
f.close()

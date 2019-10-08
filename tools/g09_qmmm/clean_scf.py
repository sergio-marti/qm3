#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

f = open( "qmmm.log", "rt" )
for l in f:
	if( l[0:10] == " SCF Done:" ):
		e = float( l.strip().split()[-5] )
	elif( l[0:29] == " Self energy of the charges =" ):
		q = float( l.strip().split()[-2] )
f.close()
print( "%27.15lE"%( e - q ) )

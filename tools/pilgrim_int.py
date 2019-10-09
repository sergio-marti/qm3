#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

import	matplotlib.pyplot as plt
import	qm3.maths.interpolation

s = []
e = []
f = open( sys.argv[1], "rt" )
n = int( f.readline().strip() )
l = f.readline()
while( l != "" ):
	t = l.strip().split()
	_e = float( t[3] )
	_s = float( t[8] )
	if( _s != 0.0 ):
		e.append( _e )
		s.append( _s )
	for i in range( n + 1 ):
		f.readline()
	l = f.readline()
f.close()

o = qm3.maths.interpolation.lagrange( s, e )
i1 = o.calc( 0.0 )[0]
o = qm3.maths.interpolation.cubic_spline( s, e )
i2 = o.calc( 0.0 )[0]
o = qm3.maths.interpolation.hermite_spline( s, e, "akima" )
i3 = o.calc( 0.0 )[0]
print( "Lagrange: ", i1 )
print( "CSpline:  ", i2 )
print( "Akima:    ", i3 )
i4 = ( i1 + i2 + i3 ) / 3.0
print( "Average:  ", i4 )

plt.clf()
plt.grid( True )
plt.plot( s, e, 'o' )
plt.plot( [ 0.0 ] , [ i4 ], 'o' )
plt.show()

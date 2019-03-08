#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = range

import	re
import	math
import	qm3.actions.fitting
import	qm3.actions.minimize


def min_func( obj ):
	qm3.actions.minimize.steepest_descent( obj, step_number = 1000, print_frequency = 1000, gradient_tolerance = 100., step_size = 0.1 )
	qm3.actions.minimize.l_bfgs( obj, step_number = 10000, print_frequency = 1000, gradient_tolerance = 0.01, step_size = 0.1 )
	qm3.actions.minimize.conjugate_gradient_plus( obj, step_number = 1000, print_frequency = 100, gradient_tolerance = 0.01 )

def dynamo_dihedral( x, V ):
	return( 0.5 * sum( [ j*(1.+k*math.cos(i*x)) for i,j,k in zip( [0., 1., 2., 3.], V, [0., 1., -1., 1.] ) ] ) )


N = 1000
D = 2. * math.pi / float( N - 1 )
X = []
for i in range( N ):
	rx = - math.pi + float( i ) * D
	X.append( rx )

data = [ {}, {}, {}, {}, {} ]
keyw = None
nbpt = re.compile( "([^\ ]+)[\ ]+[0-9\.]+[\ ]+([0-9\.\-]+)[\ ]+([0-9\.]+)" )
f = file( sys.argv[1], "rt" )
l = f.readline()
while( l != "" ):
	if( l[0:4].lower() == "bond" ):
		keyw = 0
	if( l[0:4].lower() == "angl" ):
		keyw = 1
	if( l[0:4].lower() == "dihe" ):
		keyw = 2
	if( l[0:4].lower() == "impr" or l[0:4].lower() == "imph" ):
		keyw = 3
	if( l[0:4].lower() == "nonb" ):
		keyw = 4
	t = l.strip().split()
	if( keyw == 0 and len( t ) >= 4 ):
		data[keyw][" ".join( t[0:2] )] = t[2:4]
	if( keyw == 1 and len( t ) >= 5 ):
		data[keyw][" ".join( t[0:3] )] = t[3:5]
	if( ( keyw == 2 or keyw == 3 ) and len( t ) >= 7 ):
		k = " ".join( t[0:4] )
		if( data[keyw].has_key( k ) ):
			data[keyw][k] += t[4:7]
		else:
			data[keyw][k] = t[4:7]
	if( keyw == 4 and len( t ) >= 3 ):
		t = nbpt.findall( l )
		if( len( t ) > 0 ):
			t = t[0]
			data[keyw][t[0]] = t[1:]
	l = f.readline()
f.close()

f = file( "converted", "wt" )
f.write( "Bonds\n" )
for k in data[0].keys():
	t = k.upper().strip().split()
	f.write( "%-10s%-10s%20.3lf%8.3lf\n"%( t[0], t[1], float( data[0][k][0] ), float( data[0][k][1] ) ) )
f.write( "End\n\n" )

f.write( "Angles\n" )
for k in data[1].keys():
	t = k.upper().strip().split()
	f.write( "%-10s%-10s%-10s%20.3lf%8.3lf\n"%( t[0], t[1], t[2], float( data[1][k][0] ), float( data[1][k][1] ) ) )
f.write( "End\n\n" )

f.write( "Dihedrals\n" )
for k in data[2].keys():
	T = k.upper().strip().split()
	Y = [ .0 for i in range( N ) ]
	t = [ float( i ) for i in data[2][k] ]
	n = len( t )
	for i in range( 0, n, 3 ):
		r_k = float( t[i] )
		r_n = float( t[i+1] )
		r_d = float( t[i+2] ) / 180. * math.pi
		for j in range( N ):
			Y[j] += r_k * ( 1. + math.cos( r_n * X[j] - r_d ) )
	rX = []
	rY = []
	for j in range( N ):
		if( Y[j] <= 200. ):
			rX.append( X[j] )
			rY.append( Y[j] )
	try:
		o = qm3.actions.fitting.problem( rX, rY, dynamo_dihedral, [ .0, .0, .0, .0 ] )
		o.fit( min_func )
		o.table( "%s-%s-%s-%s.log"%( T[0], T[1], T[2], T[3] ) )
		print( "@ %s-%s-%s-%s: %.6lf "%( T[0], T[1], T[2], T[3], o.Rsq ), t )
		if( o.Rsq > .8 ):
			f.write( "%-10s%-10s%-10s%-10s%20.4lf%9.4lf%9.4lf%9.4lf\n"%( T[0], T[1], T[2], T[3], o.coor[0], o.coor[1], o.coor[2], o.coor[3] ) )
		else:
			f.write( "!= %-10s%-10s%-10s%-10s%20.4lf%9.4lf%9.4lf%9.4lf\n"%( T[0], T[1], T[2], T[3], o.coor[0], o.coor[1], o.coor[2], o.coor[3] ) )
	except:
		f.write( "%-10s%-10s%-10s%-10s%20.4lf%9.4lf%9.4lf%9.4lf\n"%( T[0], T[1], T[2], T[3], .0, .0, .0, .0 ) )
f.write( "End\n\n" )

f.write( "Impropers\n" )
for k in data[3].keys():
	T = k.upper().strip().split()
	Y = [ .0 for i in range( N ) ]
	t = [ float( i ) for i in data[3][k] ]
	n = len( t )
	for i in range( 0, n, 3 ):
		r_k = float( t[i] )
		r_d = float( t[i+2] )
		for j in range( N ):
			Y[j] += r_k * math.pow( X[j] * 180. / math.pi - r_d, 2 )
	rX = []
	rY = []
	for j in range( N ):
		if( Y[j] <= 200. ):
			rX.append( X[j] )
			rY.append( Y[j] )
	try:
		o = qm3.actions.fitting.problem( rX, rY, dynamo_dihedral, [ .0, .0, .0, .0 ] )
		o.fit( min_func )
		o.table( "%s-%s-%s-%s.log"%( T[0], T[1], T[2], T[3] ) )
		print( "@ %s-%s-%s-%s: %.6lf "%( T[0], T[1], T[2], T[3], o.Rsq ), t )
		if( o.Rsq > .8 ):
			f.write( "%-10s%-10s%-10s%-10s%20.4lf%12.4lf%12.4lf%12.4lf\n"%( T[0], T[1], T[2], T[3], o.coor[0], o.coor[1], o.coor[2], o.coor[3] ) )
		else:
			f.write( "!= %-10s%-10s%-10s%-10s%20.4lf%12.4lf%12.4lf%12.4lf\n"%( T[0], T[1], T[2], T[3], o.coor[0], o.coor[1], o.coor[2], o.coor[3] ) )
	except:
		f.write( "%-5s%-5s%-5s%-5s%20.4lf%12.4lf%12.4lf%12.4lf\n"%( T[0], T[1], T[2], T[3], .0, .0, .0, .0 ) )
f.write( "End\n\n" )

f.write( "Types\n" )
c = 2. / math.pow( 2., 1. / 6. )
for k in data[4].keys():
	f.write( "%-5s%20.4lf%9.4lf\n"%( k.upper(), float( data[4][k][1] ) * c, - float( data[4][k][0] ) ) )
f.write( "End\n" )
f.close()


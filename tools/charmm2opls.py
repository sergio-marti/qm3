#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = range

import	re
import	math
import	collections


if( len( sys.argv ) != 3 ):
	print( "%s Topology(RTF) Parameters(PRM)"%( sys.argv[0] ) )
	sys.exit( 1 )


opls_type = """MM_Definitions OPLS_AA 1.0

Types
! ---------------------------
"""

opls_resi = """
Electrostatics Scale 0.5
Lennard_Jones  Scale 0.5
Units kcal/mole

Residues
!-------------------------------------------------------------------------------
"""

opls_bond = """
Bonds
! ------------------------
"""

opls_angl = """
Angles
! -----------------------------
"""

opls_dihe = """
Dihedrals
! -----------------------------------------------------
"""

opls_impr = """
Impropers
! -----------------------------------------------------
"""


mass = { 1 : 1.00794, 2 : 4.00260, 3 : 6.94100, 4 : 9.01218, 5 : 10.8110, 6 : 12.0107, 7 : 14.0067, 8 : 15.9994, 9 : 18.9984, 10 : 20.1797,
	11 : 22.9898, 12 : 24.3050, 13 : 26.9815, 14 : 28.0855, 15 : 30.9738, 16 : 32.0650, 17 : 35.4530, 18 : 39.9480, 19 : 39.0983, 20 : 40.0780,
	21 : 44.9559, 22 : 47.8670, 23 : 50.9415, 24 : 51.9961, 25 : 54.9380, 26 : 55.8450, 27 : 58.9332, 28 : 58.6934, 29 : 63.5460, 30 : 65.3900,
	31 : 69.7230, 32 : 72.6400, 33 : 74.9216, 34 : 78.9600, 35 : 79.9040, 36 : 83.8000, 37 : 85.4678, 38 : 87.6200, 39 : 88.9059, 40 : 91.2240,
	41 : 92.9064, 42 : 95.9400, 43 : 98.9063, 44 : 101.0700, 45 : 102.9060, 46 : 106.4200, 47 : 107.8680, 48 : 112.4110, 49 : 114.8180, 50 : 118.7100,
	51 : 121.7600, 52 : 127.6000, 53 : 126.9040, 54 : 131.2930, 55 : 132.9050, 56 : 137.2370, 57 : 138.9050, 58 : 140.1160, 59 : 140.9080, 60 : 144.2400,
	61 : 146.9150, 62 : 150.3600, 63 : 151.9640, 64 : 157.2500, 65 : 158.9250, 66 : 162.5000, 67 : 164.9300, 68 : 167.2590, 69 : 168.9340, 70 : 173.0400,
	71 : 174.9670, 72 : 178.4900, 73 : 180.9480, 74 : 183.8400, 75 : 186.2070, 76 : 190.2300, 77 : 192.2170, 78 : 195.0780, 79 : 196.9670, 80 : 200.5900,
	81 : 204.3830, 82 : 207.2000, 83 : 208.9800, 84 : 208.9820, 85 : 209.9870, 86 : 222.0180, 87 : 223.0200, 88 : 226.0250, 89 : 227.0280, 90 : 232.0380,
	91 : 231.0360, 92 : 238.0290, 93 : 237.0480, 94 : 244.0640, 95 : 243.0610, 96 : 247.0700, 97 : 247.0700, 98 : 251.0800, 99 : 252.0830, 100 : 257.0950,
	101 : 258.0990, 102 : 259.1010, 103 : 262.1100, 104 : 261.1090, 105 : 262.1140, 106 : 263.1190, 107 : 262.1230, 108 : 265.1310, 109 : 266.1380 }


atom_type = collections.OrderedDict()
resi_name = None
resi_atom = []
resi_type = []
resi_chrg = []
resi_bond = []
resi_impr = []
f = open( sys.argv[1], "rt" )
l = f.readline()
while( l != "" ):
	t = l.upper().strip().split()
	if( len( t ) == 4 and t[0] == "MASS" ):
		m = float( t[3] )
		atom_type[t[2]] = sorted( [ ( math.fabs( m - mass[i] ), i ) for i in list( mass ) ] )[0][1]
	elif( len( t ) == 3 and t[0] == "RESI" ):
		resi_name = t[1]
	elif( len( t ) == 4 and t[0] == "ATOM" ):
		resi_atom.append( t[1] )
		resi_type.append( t[2] )
		resi_chrg.append( float( t[3] ) )
	elif( len( t ) >= 3 and t[0] == "BOND" ):
		resi_bond.append( " ".join( t[1:3] ) )
	elif( len( t ) >= 5 and t[0] == "IMPH" ):
		resi_impr.append( " ".join( t[1:5] ) ) #( t[1], t[2], t[3], t[4] ) )
	l = f.readline()
f.close()


ffld = { "bond": collections.OrderedDict(),
		"angl": collections.OrderedDict(),
		"dihe": collections.OrderedDict(),
		"impr": collections.OrderedDict(),
		"nbnd": collections.OrderedDict() }
f = open( sys.argv[2], "rt" )
key = None
pat = re.compile( "([^\ ]+)[\ ]+[0-9\.]+[\ ]+([0-9\.\-]+)[\ ]+([0-9\.]+)" )
cte = 2. / math.pow( 2., 1. / 6. )
l = f.readline()
while( l != "" ):
	if( l[0:4].upper() == "BOND" ):
		key = "bond"
	elif( l[0:4].upper() == "ANGL" ):
		key = "angl"
	elif( l[0:4].upper() == "DIHE" ):
		key = "dihe"
	elif( l[0:4].upper() in [ "IMPH", "IMPR" ] ):
		key = "impr"
	elif( l[0:4].upper() == "NONB" ):
		key = "nbnd"
	# --------------------------------------
	t = l.upper().strip().split()
	if( key == "bond" and len( t ) >= 4 ):
		ffld[key][" ".join( t[0:2] )] = ( float( t[2] ), float( t[3] ) )
	elif( key == "angl" and len( t ) >= 5 ):
		ffld[key][" ".join( t[0:3] )] = ( float( t[3] ), float( t[4] ) )
	elif( key in [ "dihe", "impr" ] and len( t ) >= 7 ):
		k = " ".join( t[0:4] )
		if( ffld[key].has_key( k ) ):
			ffld[key][k] += [ float( i ) for i in t[4:7] ]
		else:
			ffld[key][k] = [ float( i ) for i in t[4:7] ]
	elif( key == "nbnd" and len( t ) >= 3 ):
		t = pat.findall( l.upper() )
		if( len( t ) > 0 ):
			ffld[key][t[0][0]] = ( cte * float( t[0][2] ), math.fabs( float( t[0][1] ) ) )
	l = f.readline()
f.close()



def __dihe( v, x, y ):
	f = 0
	g = [ .0, .0, .0, .0 ]
	n = len( x )
	for i in range( n ):
		t_1   = ( 1.0 + math.cos( x[i] ) )
		t_2   = ( 1.0 - math.cos( 2.0 * x[i] ) )
		t_3   = ( 1.0 + math.cos( 3.0 * x[i] ) )
		t     = 0.5 * ( v[0] + v[1] * t_1 + v[2] * t_2 + v[3] * t_3 )
		d     = ( t - y[i] )
		f    += d * d
		g[0] += d
		g[1] += d * t_1
		g[2] += d * t_2
		g[3] += d * t_3
	return( f, g )


def __impr( v, x, y ):
	f = 0
	g = [ .0, .0, .0, .0 ]
	n = len( x )
	for i in range( n ):
		t_1   = ( 1.0 + math.cos( x[i] ) )
		t_2   = ( 1.0 - math.cos( 2.0 * x[i] ) )
		t     = 0.5 * ( v[0] + v[1] * t_1 + v[2] * t_2 )
		d     = ( t - y[i] )
		f    += d * d
		g[0] += d
		g[1] += d * t_1
		g[2] += d * t_2
	return( f, g )


def __fire( function, x, y ):
	msiz = 0.1
	coor = [ .0, .0, .0, .0 ]
	nstp = 0
	alph = 0.1
	ssiz = msiz
	velo = [ .0, .0, .0, .0 ]
	step = [ .0, .0, .0, .0 ]
	func, grad = function( coor, x, y )
	norm = math.sqrt( grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2] + grad[3] * grad[3] )
	it   = 0
	while( it < 5000 and norm > 1.e-8 ):
		vsiz = math.sqrt( velo[0] * velo[0] + velo[1] * velo[1] + velo[2] * velo[2] + velo[3] * velo[3] )
		vfac = - ( velo[0] * grad[0] + velo[1] * grad[1] + velo[2] * grad[2] + velo[3] * grad[3] )
		if( vfac > 0.0 ):
			for i in [0, 1, 2, 3]:
				velo[i] = ( 1.0 - alph ) * velo[i] - alph * grad[i] / norm * vsiz
			if( nstp > 5 ):
				ssiz  = min( ssiz * 1.1, msiz )
				alph *= 0.99
			nstp += 1
		else:
			velo  = [ .0, .0, .0, .0 ]
			alph  = 0.1
			nstp  = 0
			ssiz *= 0.5
		for i in [0, 1, 2, 3]:
			velo[i] -= ssiz * grad[i]
			step[i]  = ssiz * velo[i]
		disp = math.sqrt( step[0] * step[0] + step[1] * step[1] + step[2] * step[2] + step[3] * step[3] )
		if( disp > ssiz ):
			for i in [0, 1, 2, 3]:
				step[i] *= ssiz / disp
		for i in [0, 1, 2, 3]:
			coor[i] += step[i]
		func, grad = function( coor, x, y )
		norm = math.sqrt( grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2] + grad[3] * grad[3] )
		it  += 1
	print( "%20.10lf%20.10lf"%( func, norm ) )
	return( coor, func, norm )
	



N  = 100
_d = 2. * math.pi / float( N - 1 )
X  = [ - math.pi + float( i ) * _d for i in range( N ) ]

f = open( "opls", "wt" )
f.write( opls_type )
for k in list( atom_type ):
	f.write( "x%-4s%4d%10.5lf%10.5lf\n"%( k, atom_type[k], ffld["nbnd"][k][0], ffld["nbnd"][k][1] ) )
f.write( opls_resi )
f.write( "Residue %s\n"%( resi_name ) )
n = len( resi_atom )
f.write( "%d %d %d\n"%( n, len( resi_bond ), len( resi_impr ) ) )
for i in range( n ):
	f.write( "%-9sx%-5s%8.3lf\n"%( resi_atom[i], resi_type[i], resi_chrg[i] ) )
f.write( "\n" )
f.write( "\n".join( resi_bond ) )
f.write( "\n\n" )
f.write( "\n".join( resi_impr ) )
f.write( "\n" )
f.write( opls_bond )
for k in list( ffld["bond"] ):
	t = k.split()
	f.write( "x%-4sx%-4s%7.1lf%9.3lf\n"%( t[0], t[1], ffld["bond"][k][0], ffld["bond"][k][1] ) )
f.write( opls_angl )
for k in list( ffld["angl"] ):
	t = k.split()
	f.write( "x%-4sx%-4sx%-4s%6.1lf%10.2lf\n"%( t[0], t[1], t[2], ffld["angl"][k][0], ffld["angl"][k][1] ) )
f.write( opls_dihe )
for k in list( ffld["dihe"] ):
	y = [ 0.0 for i in range( N ) ]
	for i in range( 0, len( ffld["dihe"][k] ), 3 ):
		d = ffld["dihe"][k][i+2] / 180.0 * math.pi
		for j in range( N ):
			y[j] += ffld["dihe"][k][i] * ( 1 + math.cos( ffld["dihe"][k][i+1] * X[j] + d ) )
	rx = []
	ry = []
	for i in range( N ):
		if( math.fabs( y[i] ) <= 200.0 ):
			rx.append( X[i] )
			ry.append( y[i] )
	print( "%-30s"%( k ), end = "" )
	coor, func, norm = __fire( __dihe, rx, ry )
	if( func >= 10.0 ):
		f.write( "!-[fit] func/norm: %20.10lf / %20.10lf\n"%( func, norm ) )
	t = k.split()
	for i in [0, 1, 2, 3]:
		if( t[i] != "X" ):
			f.write( "x%-4s"%( t[i] ) )
		else:
			f.write( "X    " )
	f.write( "%8.3lf%9.3lf%9.3lf%9.3lf\n"%( coor[0], coor[1], coor[2], coor[3] ) )	
f.write( opls_impr )
for k in list( ffld["impr"] ):
	y = [ 0.0 for i in range( N ) ]
	for i in range( 0, len( ffld["impr"][k] ), 3 ):
		d = ffld["impr"][k][i+2] / 180.0 * math.pi
		for j in range( N ):
			y[j] += ffld["impr"][k][i] * ( 1 + math.cos( ffld["impr"][k][i+1] * X[j] + d ) )
	rx = []
	ry = []
	for i in range( N ):
		if( math.fabs( y[i] ) <= 200.0 ):
			rx.append( X[i] )
			ry.append( y[i] )
	print( "%-30s"%( k ), end = "" )
	coor, func, norm = __fire( __impr, rx, ry )
	if( func >= 10.0 ):
		f.write( "!-[fit] func/norm: %20.10lf / %20.10lf\n"%( func, norm ) )
	t = k.split()
	for i in [0, 1, 2, 3]:
		if( t[i] != "X" ):
			f.write( "x%-4s"%( t[i] ) )
		else:
			f.write( "X    " )
	f.write( "%8.3lf%9.3lf%9.3lf%9.3lf\n"%( coor[0], coor[1], coor[2], coor[3] ) )	
f.close()

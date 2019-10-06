#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	re
import	os
import	struct



_GAUS = """%nproc=1
%mem=2048mb
#p external=qmmm opt=(cart,ts,noeigentest) nosymm

.

0 1
"""



if( len( sys.argv ) != 5 ):
	print( "%s selection pdb psf non_bonded"%( sys.argv[0] ) )
	sys.exit( 1 )


_MAS = { 1 : 1.00794, 2 : 4.00260, 3 : 6.94100, 4 : 9.01218, 5 : 10.8110, 6 : 12.0107, 7 : 14.0067, 8 : 15.9994, 9 : 18.9984, 10 : 20.1797,
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

_SMB  = { 1 : "H", 2 : "He", 3 : "Li", 4 : "Be", 5 : "B", 6 : "C", 7 : "N", 8 : "O", 9 : "F", 10 : "Ne",
	11 : "Na", 12 : "Mg", 13 : "Al", 14 : "Si", 15 : "P", 16 : "S", 17 : "Cl", 18 : "Ar", 19 : "K", 20 : "Ca",
	21 : "Sc", 22 : "Ti", 23 : "V", 24 : "Cr", 25 : "Mn", 26 : "Fe", 27 : "Co", 28 : "Ni", 29 : "Cu", 30 : "Zn",
	31 : "Ga", 32 : "Ge", 33 : "As", 34 : "Se", 35 : "Br", 36 : "Kr", 37 : "Rb", 38 : "Sr", 39 : "Y", 40 : "Zr",
	41 : "Nb", 42 : "Mo", 43 : "Tc", 44 : "Ru", 45 : "Rh", 46 : "Pd", 47 : "Ag", 48 : "Cd", 49 : "In", 50 : "Sn",
	51 : "Sb", 52 : "Te", 53 : "I", 54 : "Xe", 55 : "Cs", 56 : "Ba", 57 : "La", 58 : "Ce", 59 : "Pr", 60 : "Nd",
	61 : "Pm", 62 : "Sm", 63 : "Eu", 64 : "Gg", 65 : "Tb", 66 : "Dy", 67 : "Ho", 68 : "Er", 69 : "Tm", 70 : "Yb",
	71 : "Lu", 72 : "Hf", 73 : "Ta", 74 : "W", 75 : "Re", 76 : "Os", 77 : "Ir", 78 : "Pt", 79 : "Au", 80 : "Hg",
	81 : "Tl", 82 : "Pb", 83 : "Bi", 84 : "Po", 85 : "At", 86 : "Rn", 87 : "Fr", 88 : "Ra", 89 : "Ac", 90 : "Th",
	91 : "Pa", 92 : "U", 93 : "Np", 94 : "Pu", 95 : "Am", 96 : "Cm", 97 : "Bk", 98 : "Cf", 99 : "Es", 100 : "Fm",
	101 : "Md", 102 : "No", 103 : "Lr", 104 : "Rf", 105 : "Db", 106 : "Sg", 107 : "Bh", 108 : "Hs", 109 : "Mt" }


labl = []
kind = []
coor = []
chrg = []
mass = []
nbnd = {}
qmat = []
mmat = []


# -- load selections
f = open( sys.argv[1], "rt" )
for i in range( int( f.readline().strip() ) ):
	qmat.append( int( f.readline().strip() ) )
for i in range( int( f.readline().strip() ) ):
	mmat.append( int( f.readline().strip() ) )
f.close()
qmat.sort()
mmat.sort()
mmat = list( set( mmat ).difference( set( qmat ) ) )


# -- load PDB file
f = open( sys.argv[2], "rt" )
for l in f:
	if( l[0:4] == "ATOM" or l[0:4] == "HETA" ):
		labl.append( l[12:17].strip() )
		coor += [ float( l[30:38] ), float( l[38:46] ), float( l[46:54] ) ]
f.close()


# -- load PSF
bnds = []
f = open( sys.argv[3], "rt" )
f.readline(); f.readline()
for i in range( int( f.readline().split()[0] ) + 1 ):
	f.readline()
n = int( f.readline().split()[0] )
if( len( labl ) != n ):
	print( "- Invalid number of atoms in PSF!" )
	sys.exit( 2 )
for i in range( n ):
	t = f.readline().strip().split()
	if( labl[i] != t[4] ):
		print( "- Wrong data (%d): %s/%s"%( i+1, labl[i], t[4] ) )
		sys.exit( 3 )
	kind.append( t[5] )
	chrg.append( float( t[6] ) )
	mass.append( float( t[7] ) )
f.readline()
n = int( f.readline().strip().split()[0] )
i = 0
while( i < n ):
	t = [ int( j ) - 1 for j in f.readline().strip().split() ]
	m = len( t ) // 2
	i += m
	for j in range( m ):
		bnds.append( [ t[2*j], t[2*j+1] ] )
f.close()


# -- calculate exclussions...
remo = []
excl = []
atmm = []
conn = []
laqm = []
lamm = []
nx12 = 0
nx13 = 0
nx14 = 0
for i in range( len( labl ) ):
	atmm.append( True )
	conn.append( [] )
for i in qmat:
	atmm[i] = False
for i,j in bnds:
	conn[i].append( j )
	conn[j].append( i )
for i in qmat:
	for j in conn[i]:
		if( j != i and atmm[j] ):
			laqm.append( i )
			lamm.append( j )
			remo.append( [ i, j ] )
			nx12 += 1
		for k in conn[j]:
			if( k != i and atmm[k] ):
				remo.append( [ i, k ] )
				nx13 += 1
			for l in conn[k]:
				if( k != i and l != j and l != i and atmm[l] ):
					excl.append( [ i, l ] )
					nx14 += 1
print( ">> %d exclussions generated (1-2:%d, 1-3:%d, 1-4:%d)"%( nx12 + nx13 + nx14, nx12, nx13, nx14 ) )


# -- build "qmmm.chrg"
mmat = list( set( mmat ).difference( set( lamm ) ) )
f = open( "qmmm.chrg", "wt" )
for i in range( len( mmat ) ):
	i3 = mmat[i] * 3
	f.write( "%20.10lf%20.10lf%20.10lf%8.3lf\n"%( coor[i3], coor[i3+1], coor[i3+2], chrg[mmat[i]] ) )
f.close()


# -- build "input.com" template
f = open( "input.com", "wt" )
f.write( _GAUS )
for i in qmat:
	i3  = i * 3
	sym = _SMB[sorted( [ ( math.fabs( _MAS[j] - mass[i] ), j ) for j in iter( _MAS ) ] )[0][1]]
	fix = 0
	if( i in laqm ):
		fix = -1
	f.write( "%-4s%4d%20.10lf%20.10lf%20.10lf\n"%( sym, fix, coor[i3], coor[i3+1], coor[i3+2] ) )
for i,j in zip( laqm, lamm ):
	i3 = i * 3
	j3 = j * 3
	vv = [ ii - jj for ii,jj in zip( coor[j3:j3+3], coor[i3:i3+3] ) ]
	mm = 1.1 / math.sqrt( sum( [ ii * ii for ii in vv ] ) )
	f.write( "%-4s%4d%20.10lf%20.10lf%20.10lf\n"%( " H", 0, coor[i3] + vv[0] * mm, coor[i3+1] + vv[1] * mm, coor[i3+2] + vv[2] * mm ) )
f.write( "\n\n\n" )
f.close()


# -- build binary "qmmm.nbnd": load NON_BONDED (kind, epsi _kcal/mol, rmin/2 _ang)
if( not os.access( "qmmm.nbnd", os.R_OK ) ):
	f = open( sys.argv[4], "rt" )
	for i,j,k in re.compile( "([A-Z0-9]+)[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)" ).findall( f.read() ):
		nbnd[i] = [ math.sqrt( float( j ) ), float( k ) ]
	f.close()
	f = open( "qmmm.nbnd", "wb" )
	f.write( struct.pack( "l", len( qmat ) * len( mmat ) - len( remo ) ) )
	nqm = len( qmat )
	nmm = len( mmat )
	for i in range( nqm ):
		for j in range( nmm ):
			# -- search for an excluded interaction (remo)
			q = True
			t = 0
			while( t < len( remo ) and q ):
				q = ( remo[t][0] == qmat[i] and remo[t][1] == mmat[j] )
				t += 1
			if( q ):
				x = 1.0
				# -- search for a scaled interaction (excl)
				q = True
				t = 0
				while( t < len( excl ) and q ):
					q = ( excl[t][0] == qmat[i] and excl[t][1] == mmat[j] )
					t += 1
				if( not q ):
					x = 0.5
				# -- register interaction
				f.write( struct.pack( "l", i ) )
				f.write( struct.pack( "l", j ) )
				f.write( struct.pack( "d", nbnd[kind[qmat[i]]][0] * nbnd[kind[mmat[j]]][0] ) )
				f.write( struct.pack( "d", nbnd[kind[qmat[i]]][1] + nbnd[kind[mmat[j]]][1] ) )
				f.write( struct.pack( "d", x ) )
	f.close()

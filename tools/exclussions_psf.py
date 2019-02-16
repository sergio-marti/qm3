#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
try:
	import	cPickle as pickle
except:
	import	pickle


#arguments: system.PSF QMSelection.pickled

f = open( sys.argv[1], "rt" )
# header
f.readline(); f.readline()
# title
for i in range( int( f.readline().strip().split()[0] ) ):
	f.readline()
# atoms
f.readline()
nat = int( f.readline().strip().split()[0] )
seg = []
res = []
lbl = []
for i in range( nat ):
	t = f.readline().strip().split()
	seg.append( t[1] )
	res.append( int( t[2] ) )
	lbl.append( t[4] )
# >> bonds <<
f.readline()
bnd = []
n = int( f.readline().strip().split()[0] )
i = 0
while( i < n ):
	t = [ int( j ) - 1 for j in f.readline().strip().split() ]
	m = len( t ) // 2
	i += m
	for j in range( m ):
		bnd.append( [ t[2*j], t[2*j+1] ] )
f.close()

smm = []
atm = []
for i in range( nat ):
	smm.append( True )
	atm.append( [] )

for i,j in bnd:
	atm[i].append( j )
	atm[j].append( i )

f = open( sys.argv[2], "rb" )
sqm = pickle.load( f )
f.close()
for i in sqm:
	smm[i] = False

# link atoms + 1-2 exclussions
ela = []
exc = []
for i in sqm:
	for j in atm[i]:
		if( j != i and smm[j] ):
			ela.append( [ i, j ] )
			exc.append( [ i, j, 0.0 ] )
			print( "%4s %4d %-6s -- %4s %4d %-6s    1-2"%( seg[i], res[i], lbl[i], seg[j], res[j], lbl[j] ) )

f = open( "sele_LA.pk", "wb" )
pickle.dump( ela, f )
f.close()

# 1-3 exclussions
for i in sqm:
	for j in atm[i]:
		for k in atm[j]:
			if( k != i and smm[k] ):
				exc.append( [ i, k, 0.0 ] )
				print( "%4s %4d %-6s -- %4s %4d %-6s    1-3"%( seg[i], res[i], lbl[i], seg[k], res[k], lbl[k] ) )

# 1-4 exclussions
for i in sqm:
	for j in atm[i]:
		for k in atm[j]:
			for l in atm[k]:
				if( k != i and l != j and l != i and smm[l] ):
					exc.append( [ i, l, 0.5 ] )
					print( "%4s %4d %-6s -- %4s %4d %-6s    1-4"%( seg[i], res[i], lbl[i], seg[l], res[l], lbl[l] ) )

exc.sort()
f = open( "sele_EX.pk", "wb" )
pickle.dump( exc, f )
f.close()
print( ">>", len( exc ), "exclussions generated!" )

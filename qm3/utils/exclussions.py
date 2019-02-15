# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

import	os
import	qm3.elements
try:
	import cPickle as pickle
except:
	import pickle



def __guess_bond( mole ):
	bond = []
	for i in range( mole.natm - 1 ):
		i3 = i * 3
		ri = qm3.elements.r_cov[mole.anum[i]] + 0.05
		for j in range( i + 1, mole.natm ):
			if( mole.anum[i] == 1 and mole.anum[j] == 1 ):
				continue
			rj = qm3.elements.r_cov[mole.anum[j]] + 0.05
			t  = ( ri + rj ) * ( ri + rj )
			if( qm3.utils.distanceSQ( mole.coor[i3:i3+3], mole.coor[j*3:j*3+3] ) <= t ):
				bond.append( [ i, j ] )
	return( bond )



def exclussions( data, sele_QM, ncpu = os.sysconf( 'SC_NPROCESSORS_ONLN' ) ):
	if( hasattr( data, "anum" ) and data.anum != [] ):
		try:
			import qm3.engines._mol_mech
			bond = qm3.engines._mol_mech.connectivity( ncpu, data )
		except:
			bond = __guess_bond( data )
		natm = data.natm
	elif( type( data ) == list and len( data[0] ) == 2 ):
		bond = data
		# -- dirty --
		natm = max( sum( bond, [] ) ) + 1
	else:
		return
	atmm = [ True ] * natm
	for i in sele_QM:
		atmm[i] = False
	conn = [[]] * natm
	for i,j in bond:
		conn[i].append( j )
		conn[j].append( i )
	latm = []
	excl = []
	nx12 = 0
	nx13 = 0
	nx14 = 0
	for i in sele_QM:
		for j in conn[i]:
			if( j != i and atmm[j] ):
				latm.append( [ i, j ] )
				excl.append( [ i, j, 0.0 ] )
				nx12 += 1
			for k in conn[j]:
				if( k != i and atmm[k] ):
					excl.append( [ i, k, 0.0 ] )
					nx13 += 1
				for l in conn[k]:
					if( k != i and l != j and l != i and atmm[l] ):
						exc.append( [ i, l, 0.5 ] )
						nx14 += 1
	f = open( "sele_LA.pk", "wb" )
	pickle.dump( latm, f )
	f.close()
	excl.sort()
	f = open( "sele_EX.pk", "wb" )
	pickle.dump( excl, f )
	f.close()
	print( ">> %d exclussions generated (1-2:%d, 1-3:%d, 1-4:%d)"%( nx12 + nx13 + nx14, nx12, nx13, nx14 ) )

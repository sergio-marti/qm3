# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.mol
import	qm3.setup._cions
import	qm3.utils
import	qm3.elements
import	qm3.constants
import	multiprocessing, qm3.utils.queue
import	qm3.maths.rand
import	time



def counter_ions( molec, fname = "ions.pdb", num = 1, chrg = 1.0, d_grd = 0.5, d_ion = 11.0, d_prt = 6.5 ):
	out = qm3.mol.molecule()
	out.natm = num
	out.coor = qm3.setup._cions.counter_ions( molec, num, chrg, d_grd, d_ion, d_prt, multiprocessing.cpu_count() )
	for i in range( num ):
		out.labl.append( "ION" )
		out.resi.append( i + 1 )
		out.resn.append( "ION" )
		out.segn.append( "I" )
	out.settle()
	return( out )



def solvate( molec, solvnt, radii = qm3.elements.r_vdw, transform = True ):
	t0 = time.time()
	sys.stderr.write( "+ Computing van-der-walls radii... " )
	r1m = []
	r2m = []
	for i in range( molec.natm ):
		r1m.append( radii[molec.anum[i]] )
		r2m.append( radii[molec.anum[i]] * radii[molec.anum[i]] )
	r1s = []
	r2s = []
	for i in range( solvnt.natm ):
		r1s.append( radii[solvnt.anum[i]] )
		r2s.append( radii[solvnt.anum[i]] * radii[solvnt.anum[i]] )
	sys.stderr.write( "done!\n" )
	# -------------------------------------------------------------------------------------------------------
	if( transform ):
		sys.stderr.write( "+ Transforming coordinates... " )
		qm3.utils.moments_of_inertia( molec.mass, molec.coor )
		qm3.utils.center( solvnt.mass, solvnt.coor )
		sys.stderr.write( "done!\n" )
	# -------------------------------------------------------------------------------------------------------
	sys.stderr.write( "+ Computing geometrical center (by residue)... " )
	gcm = []
	for i in range( len( molec.res_lim ) - 1 ):
		gc = [ 0.0, 0.0, 0.0 ]
		mi = molec.coor[3*molec.res_lim[i]:3*molec.res_lim[i]+3]
		ma = mi[:]
		nn = 0.0
		for j in range( molec.res_lim[i], molec.res_lim[i+1] ):
			nn += 1.0
			for k in [0, 1, 2]:
				mi[k] = min( mi[k], molec.coor[3*j+k] )
				ma[k] = max( ma[k], molec.coor[3*j+k] )
				gc[k] += molec.coor[3*j+k]
		gc = [ i / nn for i in gc ]
		gcm.append( ( gc[0], gc[1], gc[2], max( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip(mi,gc) ] ), sum( [ (ii-jj)*(ii-jj) for ii,jj in zip(ma,gc) ] ) ) ) )
	gcs = []
	for i in range( len( solvnt.res_lim ) - 1 ):
		gc = [ 0.0, 0.0, 0.0 ]
		mi = solvnt.coor[3*solvnt.res_lim[i]:3*solvnt.res_lim[i]+3]
		ma = mi[:]
		nn = 0.0
		for j in range( solvnt.res_lim[i], solvnt.res_lim[i+1] ):
			nn += 1.0
			for k in [0, 1, 2]:
				mi[k] = min( mi[k], solvnt.coor[3*j+k] )
				ma[k] = max( ma[k], solvnt.coor[3*j+k] )
				gc[k] += solvnt.coor[3*j+k]
		gc = [ i / nn for i in gc ]
		gcs.append( ( gc[0], gc[1], gc[2], max( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip(mi,gc) ] ), sum( [ (ii-jj)*(ii-jj) for ii,jj in zip(ma,gc) ] ) ) ) )
	sys.stderr.write( "done!\n" )
	# -------------------------------------------------------------------------------------------------------
	sel = []
	N = multiprocessing.cpu_count()
	Q = qm3.utils.queue.Queue( N )
	sys.stderr.write( "+ Num CPUs    : %d\n"%( N ) )
	for i in range( N ):
		multiprocessing.Process( target = __SolParSel, args = ( Q, N, i, molec, r1m, r2m, solvnt, r1s, r2s, gcm, gcs ) ).start()
	Q.serve()
	Q.data.sort()
	for i in Q.data:
		sel += range( solvnt.res_lim[i], solvnt.res_lim[i+1] )
	# ---------------------------------------------------------------------------------------
	out = solvnt.prune( sel )
	out.norm_resid()
	out.settle()
	sys.stderr.write( "+ Time        : %ld sec\n"%( time.time() - t0 ) )
	return( out )



def ionic_strength( molec, temp = 300.0, conc = 0.15, chains = None, mdst = 5.0 ):
	t0 = time.time()
	sys.stderr.write( "+ Computing geometrical center (by residue)... " )
	gcm = []
	wat = []
	mol = []
	tmp = []
	if( chains != None ):
		tmp = chains.strip().split()
	for i in range( len( molec.res_lim ) - 1 ):
		if( molec.segn[molec.res_lim[i]] in tmp or ( tmp == [] and molec.resn[molec.res_lim[i]][0:3].upper() in [ "HOH", "TIP", "WAT" ] ) ):
			wat.append( i )
		else:
			mol.append( i )
		gc = [ 0.0, 0.0, 0.0 ]
		mi = molec.coor[3*molec.res_lim[i]:3*molec.res_lim[i]+3]
		ma = mi[:]
		nn = 0.0
		for j in range( molec.res_lim[i], molec.res_lim[i+1] ):
			nn += 1.0
			for k in [0, 1, 2]:
				mi[k] = min( mi[k], molec.coor[3*j+k] )
				ma[k] = max( ma[k], molec.coor[3*j+k] )
				gc[k] += molec.coor[3*j+k]
		gc = [ i / nn for i in gc ]
		gcm.append( ( gc[0], gc[1], gc[2], max( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip(mi,gc) ] ), sum( [ (ii-jj)*(ii-jj) for ii,jj in zip(ma,gc) ] ) ) ) )
	sys.stderr.write( "done!\n" )
	# -------------------------------------------------------------------------------------------------------
	dens = qm3.constants.water_density( temp )
	mass = len( wat ) * ( qm3.elements.mass[8] + 2.0 * qm3.elements.mass[1] )
	NION = 2 * int( round( conc * mass / ( dens * 1000.0 ), 0 ) )
	sys.stderr.write( "+ Num NaCl    : %d\n"%( NION // 2 ) )
	# -------------------------------------------------------------------------------------------------------
	N = multiprocessing.cpu_count()
	sys.stderr.write( "+ Num CPUs    : %d\n"%( N ) )
	sys.stderr.write( "+ Pruning water molecules... " )
	Q = qm3.utils.queue.Queue( N )
	for i in range( N ):
		multiprocessing.Process( target = __IonParSel, args = ( Q, N, i, molec, gcm, wat, mol, mdst ) ).start()
	Q.serve()
	sys.stderr.write( "done!\n" )
	# -------------------------------------------------------------------------------------------------------
	sys.stderr.write( "+ Randomly placing salt ions (Na+/Cl-)... " )
	iml = qm3.mol.molecule()
	iml.natm = NION
	t = [ ( "SOD", "SOD" ), ( "CLA", "CLA" ) ]
	S = []
	for i in range( NION ):
		s = Q.data[:]
		w = qm3.maths.rand.randint( 0, len( s ) - 1 )
		iml.labl.append( t[i%2][0] )
		iml.resi.append( i+1 )
		iml.resn.append( t[i%2][1] )
		iml.segn.append( "IS" )
		iml.coor += molec.coor[3*molec.res_lim[s[w]]:3*molec.res_lim[s[w]]+3][:]
		S.append( s[w] )
		m = [ s[w] ]
		del s[w]
		Q = qm3.utils.queue.Queue( N )
		for j in range( N ):
			multiprocessing.Process( target = __IonParSel, args = ( Q, N, j, molec, gcm, s, m, mdst ) ).start()
		Q.serve()
	sys.stderr.write( "done!\n" )
	sel = []
	for i in S:
		sel += range( molec.res_lim[i], molec.res_lim[i+1] )
	iml.norm_resid()
	iml.settle()
	out = molec.prune( list( set( list( range( molec.natm ) ) ).difference( set( sel ) ) ) )
	out.append( iml )
	sys.stderr.write( "+ Time        : %ld sec\n"%( time.time() - t0 ) )
	return( out )



def __SolParSel( que, cpu, cur, mol, r1m, r2m, sol, r1s, r2s, gcm, gcs ):
	o = []
	n  = len( sol.res_lim ) - 1
	t0 = cur * n // cpu
	if( cur == cpu - 1 ):
		tf = n
	else:
		tf = ( cur + 1 ) * n // cpu
	sys.stderr.write( "> Process (%2d): %6d - %6d / %6d\n"%( cur, t0, tf, n ) )
	for i in range( t0 , tf ):
		f = False
		l = 0
		while( not f and l < len( mol.res_lim ) - 1 ):
			if( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip(gcs[i][0:3],gcm[l][0:3]) ] ) <= gcs[i][3] + gcm[l][3] + 25.0 ):
				j = sol.res_lim[i]
				while( not f and j < sol.res_lim[i+1] ):
					if( int( sol.anum[j] ) > 1 ):
						k = mol.res_lim[l]
						while( not f and k < mol.res_lim[l+1] ):
							if( int( mol.anum[k] ) > 1 ):
								r = 0.5625 * ( r2s[j] + r2m[k] + 2.0 * r1s[j] * r1m[k] )
								f |= ( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip( sol.coor[3*j:3*j+3], mol.coor[3*k:3*k+3] ) ] ) <= r )
							k += 1
					j += 1
			l += 1
		if( not f ):
			o.append( i )
	que.client( o )
	sys.stderr.write( "< Process (%2d): +%d molecules...\n"%( cur, len( o ) ) )



def __IonParSel( que, cpu, cur, mol, gcm, s_wat, s_mol, mdst ):
	r  = ( mdst + 5.0 ) * ( mdst + 5.0 )
	x  = mdst * mdst
	o  = []
	n  = len( s_wat )
	t0 = cur * n // cpu
	if( cur == cpu - 1 ):
		tf = n
	else:
		tf = ( cur + 1 ) * n // cpu
#	sys.stderr.write( "> Process (%2d): %6d - %6d / %6d\n"%( cur, t0, tf, n ) )
	for i in range( t0 , tf ):
		f = False
		l = 0
		while( not f and l < len( s_mol ) ):
			if( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip(gcm[s_wat[i]][0:3],gcm[s_mol[l]][0:3]) ] ) <= gcm[s_wat[i]][3] + gcm[s_mol[l]][3] + r ):
				j = mol.res_lim[s_wat[i]]
				while( not f and j < mol.res_lim[s_wat[i]+1] ):
					k = mol.res_lim[s_mol[l]]
					while( not f and k < mol.res_lim[s_mol[l]+1] ):
						f |= ( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip( mol.coor[3*j:3*j+3], mol.coor[3*k:3*k+3] ) ] ) <= x )
						k += 1
					j += 1
			l += 1
		if( not f ):
			o.append( s_wat[i] )
	que.client( o )


# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.constants



def __maphes4( molec, hess, hind ):
	n  = int( math.sqrt( len( molec.hess ) ) )
	n2 = n * 2
	if( hind[0] > -1 ):
		jj = 3 * ( n * hind[0] + hind[0] )
		molec.hess[jj]      += hess[0]
		molec.hess[jj+1]    += hess[1]
		molec.hess[jj+2]    += hess[2]
		molec.hess[jj+n]    += hess[12]
		molec.hess[jj+n+1]  += hess[13]
		molec.hess[jj+n+2]  += hess[14]
		molec.hess[jj+n2]   += hess[24]
		molec.hess[jj+n2+1] += hess[25]
		molec.hess[jj+n2+2] += hess[26]
	if( hind[1] > -1 ):
		jj = 3 * ( n * hind[1] + hind[1] )
		molec.hess[jj]      += hess[39]
		molec.hess[jj+1]    += hess[40]
		molec.hess[jj+2]    += hess[41]
		molec.hess[jj+n]    += hess[51]
		molec.hess[jj+n+1]  += hess[52]
		molec.hess[jj+n+2]  += hess[53]
		molec.hess[jj+n2]   += hess[63]
		molec.hess[jj+n2+1] += hess[64]
		molec.hess[jj+n2+2] += hess[65]
	if( hind[2] > -1 ):
		jj = 3 * ( n * hind[2] + hind[2] )
		molec.hess[jj]      += hess[78]
		molec.hess[jj+1]    += hess[79]
		molec.hess[jj+2]    += hess[80]
		molec.hess[jj+n]    += hess[90]
		molec.hess[jj+n+1]  += hess[91]
		molec.hess[jj+n+2]  += hess[92]
		molec.hess[jj+n2]   += hess[102]
		molec.hess[jj+n2+1] += hess[103]
		molec.hess[jj+n2+2] += hess[104]
	if( hind[3] > -1 ):
		jj = 3 * ( n * hind[3] + hind[3] )
		molec.hess[jj]      += hess[117]
		molec.hess[jj+1]    += hess[118]
		molec.hess[jj+2]    += hess[119]
		molec.hess[jj+n]    += hess[129]
		molec.hess[jj+n+1]  += hess[130]
		molec.hess[jj+n+2]  += hess[131]
		molec.hess[jj+n2]   += hess[141]
		molec.hess[jj+n2+1] += hess[142]
		molec.hess[jj+n2+2] += hess[143]
	if( hind[0] > -1 and hind[1] > -1 ):
		jj = 3 * ( n * hind[0] + hind[1] )
		molec.hess[jj]      += hess[3]
		molec.hess[jj+1]    += hess[4]
		molec.hess[jj+2]    += hess[5]
		molec.hess[jj+n]    += hess[15]
		molec.hess[jj+n+1]  += hess[16]
		molec.hess[jj+n+2]  += hess[17]
		molec.hess[jj+n2]   += hess[27]
		molec.hess[jj+n2+1] += hess[28]
		molec.hess[jj+n2+2] += hess[29]
		jj = 3 * ( n * hind[1] + hind[0] )
		molec.hess[jj]      += hess[36]
		molec.hess[jj+1]    += hess[37]
		molec.hess[jj+2]    += hess[38]
		molec.hess[jj+n]    += hess[48]
		molec.hess[jj+n+1]  += hess[49]
		molec.hess[jj+n+2]  += hess[50]
		molec.hess[jj+n2]   += hess[60]
		molec.hess[jj+n2+1] += hess[61]
		molec.hess[jj+n2+2] += hess[62]
	if( hind[0] > -1 and hind[2] > -1 ):
		jj = 3 * ( n * hind[0] + hind[2] )
		molec.hess[jj]      += hess[6]
		molec.hess[jj+1]    += hess[7]
		molec.hess[jj+2]    += hess[8]
		molec.hess[jj+n]    += hess[18]
		molec.hess[jj+n+1]  += hess[19]
		molec.hess[jj+n+2]  += hess[20]
		molec.hess[jj+n2]   += hess[30]
		molec.hess[jj+n2+1] += hess[31]
		molec.hess[jj+n2+2] += hess[32]
		jj = 3 * ( n * hind[2] + hind[0] )
		molec.hess[jj]      += hess[72]
		molec.hess[jj+1]    += hess[73]
		molec.hess[jj+2]    += hess[74]
		molec.hess[jj+n]    += hess[84]
		molec.hess[jj+n+1]  += hess[85]
		molec.hess[jj+n+2]  += hess[86]
		molec.hess[jj+n2]   += hess[96]
		molec.hess[jj+n2+1] += hess[97]
		molec.hess[jj+n2+2] += hess[98]
	if( hind[0] > -1 and hind[3] > -1 ):
		jj = 3 * ( n * hind[0] + hind[3] )
		molec.hess[jj]      += hess[9]
		molec.hess[jj+1]    += hess[10]
		molec.hess[jj+2]    += hess[11]
		molec.hess[jj+n]    += hess[21]
		molec.hess[jj+n+1]  += hess[22]
		molec.hess[jj+n+2]  += hess[23]
		molec.hess[jj+n2]   += hess[33]
		molec.hess[jj+n2+1] += hess[34]
		molec.hess[jj+n2+2] += hess[35]
		jj = 3 * ( n * hind[3] + hind[0] )
		molec.hess[jj]      += hess[108]
		molec.hess[jj+1]    += hess[109]
		molec.hess[jj+2]    += hess[110]
		molec.hess[jj+n]    += hess[120]
		molec.hess[jj+n+1]  += hess[121]
		molec.hess[jj+n+2]  += hess[122]
		molec.hess[jj+n2]   += hess[132]
		molec.hess[jj+n2+1] += hess[133]
		molec.hess[jj+n2+2] += hess[134]
	if( hind[1] > -1 and hind[2] > -1 ):
		jj = 3 * ( n * hind[1] + hind[2] )
		molec.hess[jj]      += hess[42]
		molec.hess[jj+1]    += hess[43]
		molec.hess[jj+2]    += hess[44]
		molec.hess[jj+n]    += hess[54]
		molec.hess[jj+n+1]  += hess[55]
		molec.hess[jj+n+2]  += hess[56]
		molec.hess[jj+n2]   += hess[66]
		molec.hess[jj+n2+1] += hess[67]
		molec.hess[jj+n2+2] += hess[68]
		jj = 3 * ( n * hind[2] + hind[1] )
		molec.hess[jj]      += hess[75]
		molec.hess[jj+1]    += hess[76]
		molec.hess[jj+2]    += hess[77]
		molec.hess[jj+n]    += hess[87]
		molec.hess[jj+n+1]  += hess[88]
		molec.hess[jj+n+2]  += hess[89]
		molec.hess[jj+n2]   += hess[99]
		molec.hess[jj+n2+1] += hess[100]
		molec.hess[jj+n2+2] += hess[101]
	if( hind[1] > -1 and hind[3] > -1 ):
		jj = 3 * ( n * hind[1] + hind[3] )
		molec.hess[jj]      += hess[45]
		molec.hess[jj+1]    += hess[46]
		molec.hess[jj+2]    += hess[47]
		molec.hess[jj+n]    += hess[57]
		molec.hess[jj+n+1]  += hess[58]
		molec.hess[jj+n+2]  += hess[59]
		molec.hess[jj+n2]   += hess[69]
		molec.hess[jj+n2+1] += hess[70]
		molec.hess[jj+n2+2] += hess[71]
		jj = 3 * ( n * hind[3] + hind[1] )
		molec.hess[jj]      += hess[111]
		molec.hess[jj+1]    += hess[112]
		molec.hess[jj+2]    += hess[113]
		molec.hess[jj+n]    += hess[123]
		molec.hess[jj+n+1]  += hess[124]
		molec.hess[jj+n+2]  += hess[125]
		molec.hess[jj+n2]   += hess[135]
		molec.hess[jj+n2+1] += hess[136]
		molec.hess[jj+n2+2] += hess[137]
	if( hind[2] > -1 and hind[3] > -1 ):
		jj = 3 * ( n * hind[2] + hind[3] )
		molec.hess[jj]      += hess[81]
		molec.hess[jj+1]    += hess[82]
		molec.hess[jj+2]    += hess[83]
		molec.hess[jj+n]    += hess[93]
		molec.hess[jj+n+1]  += hess[94]
		molec.hess[jj+n+2]  += hess[95]
		molec.hess[jj+n2]   += hess[105]
		molec.hess[jj+n2+1] += hess[106]
		molec.hess[jj+n2+2] += hess[107]
		jj = 3 * ( n * hind[3] + hind[2] )
		molec.hess[jj]      += hess[114]
		molec.hess[jj+1]    += hess[115]
		molec.hess[jj+2]    += hess[116]
		molec.hess[jj+n]    += hess[126]
		molec.hess[jj+n+1]  += hess[127]
		molec.hess[jj+n+2]  += hess[128]
		molec.hess[jj+n2]   += hess[138]
		molec.hess[jj+n2+1] += hess[139]
		molec.hess[jj+n2+2] += hess[140]




def mm_bond( molec, kumb, xref, a_i, a_j, skip_LE = 0.0, skip_BE = 9.e99,
			ffac = 1.0, grad = False, gfac = [ 1.0, 1.0 ], hess = False, hind = [ -1, -1 ] ):
	"""
	bond = force_constant * ( distance - reference )^2

	force_constant [kJ/mol.A^2]
	reference [A]
	"""
	ai = 3 * a_i
	aj = 3 * a_j
	dr = [ i-j for i,j in zip( molec.coor[ai:ai+3], molec.coor[aj:aj+3] ) ]
	r2 = sum( [ i * i for i in dr ] )
	vv = math.sqrt( r2 )
	df = kumb * ( vv - xref )
	if( vv >= skip_LE and vv <= skip_BE ):
		molec.func += df * ( vv - xref ) * ffac
		if( grad ):
			df *= 2.0 / vv
			for i in [0, 1, 2]:
				molec.grad[ai+i] += df * dr[i] * gfac[0]
				molec.grad[aj+i] -= df * dr[i] * gfac[1]
		if( hess ):
			tt  = ( 2.0 * kumb - df ) / r2
			hxx = ( tt * dr[0] * dr[0] + df )
			hxy =   tt * dr[0] * dr[1]
			hxz =   tt * dr[0] * dr[2]
			hyy = ( tt * dr[1] * dr[1] + df )
			hyz =   tt * dr[1] * dr[2]
			hzz = ( tt * dr[2] * dr[2] + df )
			# ii & jj -- hessian should have been previously initialized...
			n  = int( math.sqrt( len( molec.hess ) ) )
			n2 = n * 2
			for ii in hind:
				if( ii > -1 ):
					jj = 3 * ( n * ii + ii )
					molec.hess[jj]      += hxx
					molec.hess[jj+1]    += hxy
					molec.hess[jj+2]    += hxz
					molec.hess[jj+n]    += hxy
					molec.hess[jj+n+1]  += hyy
					molec.hess[jj+n+2]  += hyz
					molec.hess[jj+n2]   += hxz
					molec.hess[jj+n2+1] += hyz
					molec.hess[jj+n2+2] += hzz
			# ij & ji (only if both atoms are involved...)
			if( hind[0] > -1 and hind[1] > -1 ):
				for ii,jj in [ ( hind[0], hind[1] ), ( hind[1], hind[0] ) ]:
					kk = 3 * ( n * ii + jj )
					molec.hess[kk]      -= hxx
					molec.hess[kk+1]    -= hxy
					molec.hess[kk+2]    -= hxz
					molec.hess[kk+n]    -= hxy
					molec.hess[kk+n+1]  -= hyy
					molec.hess[kk+n+2]  -= hyz
					molec.hess[kk+n2]   -= hxz
					molec.hess[kk+n2+1] -= hyz
					molec.hess[kk+n2+2] -= hzz
	return( vv )


def mm_angle( molec, kumb, xref, a_i, a_j, a_k,
			ffac = 1.0, grad = False, gfac = [ 1.0, 1.0, 1.0 ], hess = False, hind = [ -1, -1, -1 ] ):
	"""
	angle = force_constant * ( angle - reference )^2

	force_constant [kJ/mol.rad^2]
	reference [rad]
	return_value [deg]
	"""
	ai = 3 * a_i
	aj = 3 * a_j
	ak = 3 * a_k
	dij = [ i-j for i,j in zip( molec.coor[ai:ai+3], molec.coor[aj:aj+3] ) ]
	rij = math.sqrt( sum( [ i * i for i in dij ] ) )
	dij = [ i / rij for i in dij ]
	dkj = [ k-j for k,j in zip( molec.coor[ak:ak+3], molec.coor[aj:aj+3] ) ]
	rkj = math.sqrt( sum( [ i * i for i in dkj ] ) )
	dkj = [ i / rkj for i in dkj ]
	dot = sum( [ i * j for i,j in zip( dij, dkj ) ] )
	dot = min( 1.0, max( -1.0, dot ) )
	vv  = math.acos( dot )
	dv  = ( vv - xref )
	df  = kumb * dv
	molec.func += df * dv * ffac
	if( grad ):
		dx  = - 1.0 / math.sqrt( 1.0 - dot * dot )
		df *= 2.0 * dx
		dti = [ ( dkj[i] - dot * dij[i] ) / rij for i in [ 0, 1, 2 ] ]
		dtk = [ ( dij[i] - dot * dkj[i] ) / rkj for i in [ 0, 1, 2 ] ]
		dtj = [ - ( dti[i] + dtk[i] ) for i in [ 0, 1, 2 ] ]
		for i in [0, 1, 2]:
			molec.grad[ai+i] += df * dti[i] * gfac[0]
			molec.grad[aj+i] += df * dtj[i] * gfac[1]
			molec.grad[ak+i] += df * dtk[i] * gfac[2]
	if( hess ):
		disp = 1.e-4
		gbak = molec.grad[:]
		xhes = []
		for w in [ i for i,j in zip( [ ai, aj, ak ], hind ) if j > -1 ]:
			for j in [0, 1, 2]:
				cbak = molec.coor[w+j]
				molec.coor[w+j] = cbak + disp
				molec.grad = [ 0.0 for i in range( len( gbak ) ) ]
				mm_angle( molec, kumb, xref, a_i, a_j, a_k, ffac = 0.0, grad = True, gfac = [ 1.0, 1.0, 1.0 ], hess = False )
				fgrd = molec.grad[ai:ai+3] + molec.grad[aj:aj+3] + molec.grad[ak:ak+3]
				molec.coor[w+j] = cbak - disp
				molec.grad = [ 0.0 for i in range( len( gbak ) ) ]
				mm_angle( molec, kumb, xref, a_i, a_j, a_k, ffac = 0.0, grad = True, gfac = [ 1.0, 1.0, 1.0 ], hess = False )
				bgrd = molec.grad[ai:ai+3] + molec.grad[aj:aj+3] + molec.grad[ak:ak+3]
				molec.coor[w+j] = cbak
				xhes += [ ( ii - jj ) / ( 2.0 * disp ) for ii,jj in zip( fgrd, bgrd ) ]
		molec.grad = gbak[:]
		n  = int( math.sqrt( len( molec.hess ) ) )
		n2 = n * 2
		if( hind[0] > -1 ):
			jj = 3 * ( n * hind[0] + hind[0] )
			molec.hess[jj]      += xhes[0]
			molec.hess[jj+1]    += xhes[1]
			molec.hess[jj+2]    += xhes[2]
			molec.hess[jj+n]    += xhes[6]
			molec.hess[jj+n+1]  += xhes[7]
			molec.hess[jj+n+2]  += xhes[8]
			molec.hess[jj+n2]   += xhes[12]
			molec.hess[jj+n2+1] += xhes[13]
			molec.hess[jj+n2+2] += xhes[14]
		if( hind[1] > -1 ):
			jj = 3 * ( n * hind[1] + hind[1] )
			molec.hess[jj]      += xhes[30]
			molec.hess[jj+1]    += xhes[31]
			molec.hess[jj+2]    += xhes[32]
			molec.hess[jj+n]    += xhes[39]
			molec.hess[jj+n+1]  += xhes[40]
			molec.hess[jj+n+2]  += xhes[41]
			molec.hess[jj+n2]   += xhes[48]
			molec.hess[jj+n2+1] += xhes[49]
			molec.hess[jj+n2+2] += xhes[50]
		if( hind[2] > -1 ):
			jj = 3 * ( n * hind[2] + hind[2] )
			molec.hess[jj]      += xhes[60]
			molec.hess[jj+1]    += xhes[61]
			molec.hess[jj+2]    += xhes[62]
			molec.hess[jj+n]    += xhes[69]
			molec.hess[jj+n+1]  += xhes[70]
			molec.hess[jj+n+2]  += xhes[71]
			molec.hess[jj+n2]   += xhes[78]
			molec.hess[jj+n2+1] += xhes[79]
			molec.hess[jj+n2+2] += xhes[80]
		if( hind[0] > -1 and hind[1] > -1 ):
			jj = 3 * ( n * hind[0] + hind[1] )
			molec.hess[jj]      += xhes[3]
			molec.hess[jj+1]    += xhes[4]
			molec.hess[jj+2]    += xhes[5]
			molec.hess[jj+n]    += xhes[12]
			molec.hess[jj+n+1]  += xhes[13]
			molec.hess[jj+n+2]  += xhes[14]
			molec.hess[jj+n2]   += xhes[21]
			molec.hess[jj+n2+1] += xhes[22]
			molec.hess[jj+n2+2] += xhes[23]
			jj = 3 * ( n * hind[1] + hind[0] )
			molec.hess[jj]      += xhes[27]
			molec.hess[jj+1]    += xhes[28]
			molec.hess[jj+2]    += xhes[29]
			molec.hess[jj+n]    += xhes[36]
			molec.hess[jj+n+1]  += xhes[37]
			molec.hess[jj+n+2]  += xhes[38]
			molec.hess[jj+n2]   += xhes[45]
			molec.hess[jj+n2+1] += xhes[46]
			molec.hess[jj+n2+2] += xhes[47]
		if( hind[0] > -1 and hind[2] > -1 ):
			jj = 3 * ( n * hind[0] + hind[2] )
			molec.hess[jj]      += xhes[6]
			molec.hess[jj+1]    += xhes[7]
			molec.hess[jj+2]    += xhes[8]
			molec.hess[jj+n]    += xhes[15]
			molec.hess[jj+n+1]  += xhes[16]
			molec.hess[jj+n+2]  += xhes[17]
			molec.hess[jj+n2]   += xhes[24]
			molec.hess[jj+n2+1] += xhes[25]
			molec.hess[jj+n2+2] += xhes[26]
			jj = 3 * ( n * hind[2] + hind[0] )
			molec.hess[jj]      += xhes[54]
			molec.hess[jj+1]    += xhes[55]
			molec.hess[jj+2]    += xhes[56]
			molec.hess[jj+n]    += xhes[63]
			molec.hess[jj+n+1]  += xhes[64]
			molec.hess[jj+n+2]  += xhes[65]
			molec.hess[jj+n2]   += xhes[72]
			molec.hess[jj+n2+1] += xhes[73]
			molec.hess[jj+n2+2] += xhes[74]
		if( hind[1] > -1 and hind[2] > -1 ):
			jj = 3 * ( n * hind[1] + hind[2] )
			molec.hess[jj]      += xhes[33]
			molec.hess[jj+1]    += xhes[34]
			molec.hess[jj+2]    += xhes[35]
			molec.hess[jj+n]    += xhes[42]
			molec.hess[jj+n+1]  += xhes[43]
			molec.hess[jj+n+2]  += xhes[44]
			molec.hess[jj+n2]   += xhes[51]
			molec.hess[jj+n2+1] += xhes[52]
			molec.hess[jj+n2+2] += xhes[53]
			jj = 3 * ( n * hind[2] + hind[1] )
			molec.hess[jj]      += xhes[57]
			molec.hess[jj+1]    += xhes[58]
			molec.hess[jj+2]    += xhes[59]
			molec.hess[jj+n]    += xhes[66]
			molec.hess[jj+n+1]  += xhes[67]
			molec.hess[jj+n+2]  += xhes[68]
			molec.hess[jj+n2]   += xhes[75]
			molec.hess[jj+n2+1] += xhes[76]
			molec.hess[jj+n2+2] += xhes[77]
	return( vv * qm3.constants.R2D )


def mm_dihedral( molec, data, a_i, a_j, a_k, a_l,
				ffac = 1.0, grad = False, gfac = [ 1.0, 1.0, 1.0, 1.0 ], hess = False, hind = [ -1, -1, -1, -1 ] ):
	"""
	dihedral = force_constant * ( 1 + cos( periodicity * angle - displacement ) )

	force_constant [kJ/mol]
	displacement [rad]

	data = [ frc_per=1, dsp_per=1, frc_per=2, dsp_per=2, ..., frc_per=6, dsp_per=6 ]
	"""
	ai  = 3 * a_i
	aj  = 3 * a_j
	ak  = 3 * a_k
	al  = 3 * a_l
	dji = [ ii-jj for ii,jj in zip( molec.coor[aj:aj+3], molec.coor[ai:ai+3] ) ]
	dkj = [ ii-jj for ii,jj in zip( molec.coor[ak:ak+3], molec.coor[aj:aj+3] ) ]
	rkj = math.sqrt( sum( [ ii*ii for ii in dkj ] ) )
	dlk = [ ii-jj for ii,jj in zip( molec.coor[al:al+3], molec.coor[ak:ak+3] ) ]
	vt  = [ dji[1] * dkj[2] - dkj[1] * dji[2], dji[2] * dkj[0] - dkj[2] * dji[0], dji[0] * dkj[1] - dkj[0] * dji[1] ]
	rt2 = sum( [ ii*ii for ii in vt ] )
	vu  = [ dkj[1] * dlk[2] - dlk[1] * dkj[2], dkj[2] * dlk[0] - dlk[2] * dkj[0], dkj[0] * dlk[1] - dlk[0] * dkj[1] ]
	ru2 = sum( [ ii*ii for ii in vu ] )
	vtu = [ vt[1] * vu[2] - vu[1] * vt[2], vt[2] * vu[0] - vu[2] * vt[0], vt[0] * vu[1] - vu[0] * vt[1] ]
	rtu = math.sqrt( rt2 * ru2 )
	cs1 = sum( [ ii*jj for ii,jj in zip( vt, vu ) ] ) / rtu
	cs1 = min( 1.0, max( -1.0, cs1 ) )
	sn1 = sum( [ ii*jj for ii,jj in zip( dkj, vtu ) ] ) / ( rkj * rtu )
	cs2 = cs1 * cs1 - sn1 * sn1
	sn2 = 2.0 * cs1 * sn1
	cs3 = cs1 * cs2 - sn1 * sn2
	sn3 = cs1 * sn2 + sn1 * cs2
	cs4 = cs1 * cs3 - sn1 * sn3
	sn4 = cs1 * sn3 + sn1 * cs3
	cs5 = cs1 * cs4 - sn1 * sn4
	sn5 = cs1 * sn4 + sn1 * cs4
	cs6 = cs1 * cs5 - sn1 * sn5
	sn6 = cs1 * sn5 + sn1 * cs5
	dph = 0.0
	if( data[0] != 0.0 ):
		cd  = math.cos( data[1] )
		sd  = math.sin( data[1] )
		dph += data[0] * ( cs1 * sd - sn1 * cd )
		molec.func += data[0] * ( 1.0 + cs1 * cd + sn1 * sd ) * ffac
	if( data[2] != 0.0 ):
		cd  = math.cos( data[3] )
		sd  = math.sin( data[3] )
		dph += data[2] * 2.0 * ( cs2 * sd - sn2 * cd )
		molec.func += data[2] * ( 1.0 + cs2 * cd + sn2 * sd ) * ffac
	if( data[4] != 0.0 ):
		cd  = math.cos( data[5] )
		sd  = math.sin( data[5] )
		dph += data[4] * 3.0 * ( cs3 * sd - sn3 * cd )
		molec.func += data[4] * ( 1.0 + cs3 * cd + sn3 * sd ) * ffac
	if( data[6] != 0.0 ):
		cd  = math.cos( data[7] )
		sd  = math.sin( data[7] )
		dph += data[6] * 4.0 * ( cs4 * sd - sn4 * cd )
		molec.func += data[6] * ( 1.0 + cs4 * cd + sn4 * sd ) * ffac
	if( data[8] != 0.0 ):
		cd  = math.cos( data[9] )
		sd  = math.sin( data[9] )
		dph += data[8] * 5.0 * ( cs5 * sd - sn5 * cd )
		molec.func += data[8] * ( 1.0 + cs5 * cd + sn5 * sd ) * ffac
	if( data[10] != 0.0 ):
		cd  = math.cos( data[11] )
		sd  = math.sin( data[11] )
		dph += data[10] * 6.0 * ( cs6 * sd - sn6 * cd )
		molec.func += data[10] * ( 1.0 + cs6 * cd + sn6 * sd ) * ffac
	if( grad ):
		dki = [ ii-jj for ii,jj in zip( molec.coor[ak:ak+3], molec.coor[ai:ai+3] ) ]
		dlj = [ ii-jj for ii,jj in zip( molec.coor[al:al+3], molec.coor[aj:aj+3] ) ]
		dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
				( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
				( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
		dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
				( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
				( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
		molec.grad[ai]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph * gfac[0]
		molec.grad[ai+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph * gfac[0]
		molec.grad[ai+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph * gfac[0]
		molec.grad[aj]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph * gfac[1]
		molec.grad[aj+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph * gfac[1]
		molec.grad[aj+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph * gfac[1]
		molec.grad[ak]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph * gfac[2]
		molec.grad[ak+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph * gfac[2]
		molec.grad[ak+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph * gfac[2]
		molec.grad[al]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph * gfac[3]
		molec.grad[al+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph * gfac[3]
		molec.grad[al+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph * gfac[3]
	if( hess ):
		disp = 1.e-4
		gbak = molec.grad[:]
		xhes = []
		for w in [ i for i,j in zip( [ ai, aj, ak, al ], hind ) if j > -1 ]:
			for j in [0, 1, 2]:
				cbak = molec.coor[w+j]
				molec.coor[w+j] = cbak + disp
				molec.grad = [ 0.0 for i in range( len( gbak ) ) ]
				mm_dihedral( molec, data, a_i, a_j, a_k, a_l, ffac = 0.0, grad = True, gfac = [ 1.0, 1.0, 1.0, 1.0 ], hess = False )
				fgrd = molec.grad[ai:ai+3] + molec.grad[aj:aj+3] + molec.grad[ak:ak+3] + molec.grad[al:al+3]
				molec.coor[w+j] = cbak - disp
				molec.grad = [ 0.0 for i in range( len( gbak ) ) ]
				mm_dihedral( molec, data, a_i, a_j, a_k, a_l, ffac = 0.0, grad = True, gfac = [ 1.0, 1.0, 1.0, 1.0 ], hess = False )
				bgrd = molec.grad[ai:ai+3] + molec.grad[aj:aj+3] + molec.grad[ak:ak+3] + molec.grad[al:al+3]
				molec.coor[w+j] = cbak
				xhes += [ ( ii - jj ) / ( 2.0 * disp ) for ii,jj in zip( fgrd, bgrd ) ]
		molec.grad = gbak[:]
		__maphes4( molec, xhes, hind )
	ang = qm3.constants.R2D * math.acos( cs1 )
	if( sn1 <= 0.0 ):
		ang = -ang
	return( ang )


def mm_improper( molec, kumb, xref, a_i, a_j, a_k, a_l,
				ffac = 1.0, grad = False, gfac = [ 1.0, 1.0, 1.0, 1.0 ], hess = False, hind = [ -1, -1, -1, -1 ] ):
	"""
	improper = force_constant * ( angle - reference )^2

	force_constant [kJ/mol.rad^2]
	reference [deg]
	a_i should be central atom
	"""
	ai  = 3 * a_i
	aj  = 3 * a_j
	ak  = 3 * a_k
	al  = 3 * a_l
	dji = [ ii-jj for ii,jj in zip( molec.coor[aj:aj+3], molec.coor[ai:ai+3] ) ]
	dkj = [ ii-jj for ii,jj in zip( molec.coor[ak:ak+3], molec.coor[aj:aj+3] ) ]
	dlk = [ ii-jj for ii,jj in zip( molec.coor[al:al+3], molec.coor[ak:ak+3] ) ]
	vt  = [ dji[1] * dkj[2] - dkj[1] * dji[2], dji[2] * dkj[0] - dkj[2] * dji[0], dji[0] * dkj[1] - dkj[0] * dji[1] ]
	vu  = [ dkj[1] * dlk[2] - dlk[1] * dkj[2], dkj[2] * dlk[0] - dlk[2] * dkj[0], dkj[0] * dlk[1] - dlk[0] * dkj[1] ]
	vtu = [ vt[1] * vu[2] - vu[1] * vt[2], vt[2] * vu[0] - vu[2] * vt[0], vt[0] * vu[1] - vu[0] * vt[1] ]
	rt2 = sum( [ ii*ii for ii in vt ] )
	ru2 = sum( [ ii*ii for ii in vu ] )
	rtu = math.sqrt( rt2 * ru2 )
	rkj = math.sqrt( sum( [ ii*ii for ii in dkj ] ) )
	cos = sum( [ ii*jj for ii,jj in zip( vt, vu ) ] ) / rtu
	sin = sum( [ ii*jj for ii,jj in zip( dkj, vtu ) ] ) / ( rkj * rtu )
	cos = min( 1.0, max( -1.0, cos ) )
	ang = qm3.constants.R2D * math.acos( cos )
	if( sin <= 0.0 ):
		ang = -ang
	if( math.fabs( ang + xref ) < math.fabs( ang - xref ) ):
		xref = -xref
	dt  = ang - xref
	while( dt >  180.0 ):
		dt -= 360.0
	while( dt < -180.0 ):
		dt += 360.0
	dt /= qm3.constants.R2D
	molec.func += kumb * dt * dt * ffac
	if( grad ):
		dph = 2.0 * kumb * dt
		dki = [ ii-jj for ii,jj in zip( molec.coor[ak:ak+3], molec.coor[ai:ai+3] ) ]
		dlj = [ ii-jj for ii,jj in zip( molec.coor[al:al+3], molec.coor[aj:aj+3] ) ]
		dvt = [ ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj ), 
				( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj ), 
				( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj ) ]
		dvu = [ ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj ),
				( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj ),
				( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj ) ]
		molec.grad[ai]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph * gfac[0]
		molec.grad[ai+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph * gfac[0]
		molec.grad[ai+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph * gfac[0]
		molec.grad[aj]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph * gfac[1]
		molec.grad[aj+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph * gfac[1]
		molec.grad[aj+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph * gfac[1]
		molec.grad[ak]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph * gfac[2]
		molec.grad[ak+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph * gfac[2]
		molec.grad[ak+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph * gfac[2]
		molec.grad[al]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph * gfac[3]
		molec.grad[al+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph * gfac[3]
		molec.grad[al+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph * gfac[3]
	if( hess ):
		disp = 1.e-4
		gbak = molec.grad[:]
		xhes = []
		for w in [ i for i,j in zip( [ ai, aj, ak, al ], hind ) if j > -1 ]:
			for j in [0, 1, 2]:
				cbak = molec.coor[w+j]
				molec.coor[w+j] = cbak + disp
				molec.grad = [ 0.0 for i in range( len( gbak ) ) ]
				mm_improper( molec, kumb, xref, a_i, a_j, a_k, a_l, ffac = 0.0, grad = True, gfac = [ 1.0, 1.0, 1.0, 1.0 ], hess = False )
				fgrd = molec.grad[ai:ai+3] + molec.grad[aj:aj+3] + molec.grad[ak:ak+3] + molec.grad[al:al+3]
				molec.coor[w+j] = cbak - disp
				molec.grad = [ 0.0 for i in range( len( gbak ) ) ]
				mm_improper( molec, kumb, xref, a_i, a_j, a_k, a_l, ffac = 0.0, grad = True, gfac = [ 1.0, 1.0, 1.0, 1.0 ], hess = False )
				bgrd = molec.grad[ai:ai+3] + molec.grad[aj:aj+3] + molec.grad[ak:ak+3] + molec.grad[al:al+3]
				molec.coor[w+j] = cbak
				xhes += [ ( ii - jj ) / ( 2.0 * disp ) for ii,jj in zip( fgrd, bgrd ) ]
		molec.grad = gbak[:]
		__maphes4( molec, xhes, hind )
	return( ang )




class distance( object ):
	def __init__( self, kumb, xref, indx, skip_LE = 0.0, skip_BE = 9.e99 ):
		self.kumb = kumb
		self.xref = xref
		self.indx = indx[:]
		self.skpL = skip_LE
		self.skpB = skip_BE
		self.ffac = 1.0
		self.gfac = [ 1.0, 1.0 ]
		self.hind = [ -1, -1 ]


	def get_func( self, molec ):
		return( mm_bond( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.skpL, self.skpB,
				self.ffac ) )


	def get_grad( self, molec ):
		return( mm_bond( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.skpL, self.skpB,
				self.ffac, True, self.gfac ) )


	def get_hess( self, molec ):
		return( mm_bond( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.skpL, self.skpB,
				self.ffac, True, self.gfac, True, self.hind ) )



class angle( object ):
	def __init__( self, kumb, xref, indx ):
		self.kumb = kumb
		self.xref = xref / qm3.constants.R2D
		self.indx = indx[:]
		self.ffac = 1.0
		self.gfac = [ 1.0, 1.0, 1.0 ]
		self.hind = [ -1, -1, -1 ]


	def get_func( self, molec ):
		return( mm_angle( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2],
				self.ffac ) )


	def get_grad( self, molec ):
		return( mm_angle( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2],
				self.ffac, True, self.gfac ) )


	def get_hess( self, molec ):
		return( mm_angle( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2],
				self.ffac, True, self.gfac, True, self.hind ) )



class dihedral( object ):
	def __init__( self, data, indx ):
		"""
	data = {  periodicity: [ force_constant [kJ/mol], displacement [degrees] ], ... }

	X - C_sp3 - C_sp3 - X   =>  { 3: [ 0.8159, 0.0 ] }

	valid periodicities = [ 1 : 6 ]
		"""
		self.data = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
		for i in range( 6 ):
			if( i+1 in data ):
				self.data[2*i]   = data[i+1][0]
				self.data[2*i+1] = data[i+1][1] / qm3.constants.R2D
		self.indx = indx[:]
		self.ffac = 1.0
		self.gfac = [ 1.0, 1.0, 1.0, 1.0 ]
		self.hind = [ -1, -1, -1, -1 ]


	def get_func( self, molec ):
		return( mm_dihedral( molec, self.data, self.indx[0], self.indx[1], self.indx[2], self.indx[3],
				self.ffac ) )


	def get_grad( self, molec ):
		return( mm_dihedral( molec, self.data, self.indx[0], self.indx[1], self.indx[2], self.indx[3],
				self.ffac, True, self.gfac ) )


	def get_hess( self, molec ):
		return( mm_dihedral( molec, self.data, self.indx[0], self.indx[1], self.indx[2], self.indx[3],
				self.ffac, True, self.gfac, True, self.hind ) )



class improper( object ):
	def __init__( self, kumb, xref, indx ):
		self.kumb = kumb
		self.xref = xref
		self.indx = indx[:]
		self.ffac = 1.0
		self.gfac = [ 1.0, 1.0, 1.0, 1.0 ]
		self.hind = [ -1, -1, -1, -1 ]


	def get_func( self, molec ):
		return( mm_improper( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2], self.indx[3], self.ffac ) )


	def get_grad( self, molec ):
		return( mm_improper( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2], self.indx[3],
				self.ffac, True, self.gfac ) )


	def get_hess( self, molec ):
		return( mm_improper( molec, self.kumb, self.xref, self.indx[0], self.indx[1], self.indx[2], self.indx[3],
				self.ffac, True, self.gfac, True, self.hind ) )



class multiple_distance( object ):
	def __init__( self, kumb, xref, indx, weigth ):
		"""
	multiple_distance = force_constant * ( value - reference )^2

	value = SUM weigth_i * distance_i

	force_constant [kJ/mol.A^2]
	reference [A]
		"""
		if( len( weigh ) * 2 != len( indx ) ):
			print( "- restraints.multiple_distance: Number of ATOMS should be TWICE the number of WEIGHTS!" )
			return( None )
		self.kumb = kumb
		self.xref = xref
		self.indx = indx[:]
		self.weig = weigth[:]
		self.size = len( weigth )


	def get_func( self, molec ):
		dr = []
		r  = []
		for i in range( self.size ):
			i3 = i * 3
			ii = 3 * self.indx[2*i]
			jj = 3 * self.indx[2*i+1]
			dr += [ j-k for j,k in zip( molec.coor[ii:ii+3], molec.coor[jj:jj+3] ) ]
			r.append( math.sqrt( sum( [ j * j for j in dr[i3:i3+3] ] ) ) )
		vv = sum( [ i * j for i,j in zip( r, self.weig ) ] )
		df = self.kumb * ( vv - self.xref )
		molec.func += df * ( vv - self.xref )
		return( vv )


	def get_grad( self, molec ):
		dr = []
		r  = []
		for i in range( self.size ):
			i3 = i * 3
			ii = 3 * self.indx[2*i]
			jj = 3 * self.indx[2*i+1]
			dr += [ j-k for j,k in zip( molec.coor[ii:ii+3], molec.coor[jj:jj+3] ) ]
			r.append( math.sqrt( sum( [ j * j for j in dr[i3:i3+3] ] ) ) )
		vv = sum( [ i * j for i,j in zip( r, self.weig ) ] )
		df = self.kumb * ( vv - self.xref )
		molec.func += df * ( vv - self.xref )
		for i in range( self.size ):
			i3 = i * 3
			ii = 3 * self.indx[2*i]
			jj = 3 * self.indx[2*i+1]
			tt = 2.0 * self.weig[i] * df / r[i]
			for j in [0, 1, 2]:
				molec.grad[ii+j] += tt * dr[i3+j]
				molec.grad[jj+j] -= tt * dr[i3+j]
		return( vv )



class tether( object ):
	def __init__( self, molec, kumb, indx ):
		"""
	thether = force_constant * SUM ( cartesian - reference )^2

	force_constant [kJ/mol.A^2]
	reference [A]
		"""
		self.kumb = kumb
		self.indx = {}
		for i in indx:
			self.indx[i] = molec.coor[3*i:3*i+3][:]


	def get_func( self, molec ):
		for w in self.indx:
			w3 = w * 3
			dr = [ i-j for i,j in zip( molec.coor[w3:w3+3], self.indx[w] ) ]
			molec.func += self.kumb * sum( [ i * i for i in dr ] )
		return( None )


	def get_grad( self, molec ):
		for w in self.indx:
			w3 = w * 3
			dr = [ i-j for i,j in zip( molec.coor[w3:w3+3], self.indx[w] ) ]
			molec.func += self.kumb * sum( [ i * i for i in dr ] ) 
			for i in [0, 1, 2]:
				molec.grad[w3+i] += 2.0 * self.kumb * dr[i]
		return( None )


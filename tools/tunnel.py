#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

import	os
import	math
import	qm3.mol
import	qm3.elements
import	qm3.constants
import	qm3.maths.interpolation
import	qm3.actions.rate
import	matplotlib.pyplot as plt
import	matplotlib.backends.backend_pdf



def s_coor( mol, lst ):
	o = 0.0
	for i in range( mol.natm ):
		i3 = i * 3
		o += mol.mass[i] * sum( [ (i-j)*(i-j) for i,j in zip( mol.coor[i3:i3+3], lst[i3:i3+3] ) ] )
	return( math.sqrt( o ) )



if( len( sys.argv ) < 5 ):
	print( "%s xyz reac_ener_kJ/mol prod_path reac_path [temperature=298.15]"%( sys.argv[0] ) )
	sys.exit( -1 )


temp = 298.15
if( len( sys.argv ) == 6 ):
	temp = float( sys.argv[5] )


# -- parse molecule
mol = qm3.mol.molecule()
mol.xyz_read( sys.argv[1] )
mol.guess_atomic_numbers()
mol.fill_masses()
mas = []
for i in range( mol.natm ):
	t = math.sqrt( mol.mass[i] )
	mas += [ t, t, t ]

# -- "problemize"
mol.size = 3 * mol.natm
mol.func = 0
mol.grad = [ 0.0 for i in range( mol.size ) ]
mol.hess = [ 0.0 for i in range( mol.size * mol.size ) ]

# -- reactant potential energy
ref = float( sys.argv[2] )

irc = []
# -- parse irc towards products
f  = open( sys.argv[3], "rt" )
qm3.actions.rate.path_read( f, mol )
neg, val, zpe, eta, tau = qm3.actions.rate.curvature( mol )
wig = 1.0 + 1.0 / 24.0 * math.pow( val[0] * 100.0 * qm3.constants.C * qm3.constants.H / ( qm3.constants.KB * temp ), 2.0 )
acc = 0.0
lst = mol.coor[:]
irc.append( [ 0.0, zpe, mol.func - ref + zpe, val, None, None ] )
while( qm3.actions.rate.path_read( f, mol ) ):
	acc += s_coor( mol, lst )
	lst = mol.coor[:]
	neg, val, zpe, eta, tau = qm3.actions.rate.curvature( mol )
	irc.append( [ acc, zpe, mol.func - ref + zpe, val, eta, tau ] )
f.close()

# -- parse irc towards reactants (skip first structure, already parsed...)
f  = open( sys.argv[4], "rt" )
qm3.actions.rate.path_read( f, mol )
acc = 0.0
lst = mol.coor[:]
while( qm3.actions.rate.path_read( f, mol ) ):
	acc -= s_coor( mol, lst )
	lst = mol.coor[:]
	neg, val, zpe, eta, tau = qm3.actions.rate.curvature( mol )
	irc.append( [ acc, zpe, mol.func - ref + zpe, val, eta, tau ] )
f.close()

irc.sort()

# interpolate eta/tau for the saddle
who = [ irc[i][0] for i in range( len( irc ) ) ].index( 0.0 )
crd = []
eta = []
tau = []
for i in range( len( irc ) ):
	if( i != who ):
		crd.append( irc[i][0] )
		eta.append( irc[i][4] )
		tau.append( irc[i][5] )
i_eta = qm3.maths.interpolation.hermite_spline( crd, eta )
i_tau = qm3.maths.interpolation.hermite_spline( crd, tau )
irc[who][4] = i_eta.calc( 0.0 )[0]
irc[who][5] = i_tau.calc( 0.0 )[0]

pdf = matplotlib.backends.backend_pdf.PdfPages( "tunnel.pdf" )

plt.clf()
plt.grid( True )
plt.title( "ZPE" )
plt.plot( [ irc[i][0] for i in range( len( irc ) ) ], [ irc[i][1] for i in range( len( irc ) ) ], '-o' )
pdf.savefig()
plt.close()

plt.clf()
plt.grid( True )
plt.title( "Vadi" )
plt.plot( [ irc[i][0] for i in range( len( irc ) ) ], [ irc[i][2] for i in range( len( irc ) ) ], '-o' )
pdf.savefig()
plt.close()

plt.clf()
plt.grid( True )
plt.title( "eta" )
plt.plot( crd, eta, '-o' )
pdf.savefig()
plt.close()

plt.clf()
plt.grid( True )
plt.title( "tau" )
plt.plot( crd, tau, '-o' )
pdf.savefig()
plt.close()

print()
print( "%16s%64s"%( "s ", "Frequencies " ) )
print( 80 * "-" )
for i in range( len( irc ) ):
	print( "%16.4lf%64s"%( irc[i][0], "".join( [ "%8.1lf"%( j ) for j in irc[i][3][0:8] ] ) ) )
print()
print( "WIGNER: %.4lf (%.2lf)"%( wig, temp ) )
print()
print( "%16s%16s%16s%16s%16s%16s"%( "s ", "V ", "ZPE ", "V_adi ", "eta ", "tau " ) )
print( 96 * "-" )
for i in range( len( irc ) ):
	print( "%16.4lf%16.4lf%16.4lf%16.4lf%16.4lf%16.4lf"%( irc[i][0], irc[i][2] - irc[i][1], irc[i][1], irc[i][2], irc[i][4], irc[i][5] ) )

# -- effective reduced masses...
crd = []
adi = []
eta = []
tau = []
for i in range( len( irc ) ):
	crd.append( irc[i][0] )
	adi.append( irc[i][2] )
	eta.append( irc[i][4] )
	tau.append( irc[i][5] )
mef = qm3.actions.rate.effective_reduced_masses( crd, eta, tau )

#if( os.path.isfile( "tunnel.mef" ) ):
#	f = open( "tunnel.mef", "rt" )
#	mef = [ float( i ) for i in f.read().strip().split() ]
#	f.close()

plt.clf()
plt.grid( True )
plt.title( "mueff" )
plt.plot( crd, mef, '-o' )
pdf.savefig()
plt.close()

# -- kappa calculation
k, e0, p0 = qm3.actions.rate.transmission_coefficient( crd, adi, 1.0, temp )
k, em, pm = qm3.actions.rate.transmission_coefficient( crd, adi, mef, temp )

plt.clf()
plt.grid( True )
plt.title( "prob" )
plt.plot( e0, p0, '-o' )
plt.plot( em, pm, '-o' )
pdf.savefig()
plt.close()

pdf.close()

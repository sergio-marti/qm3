# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	qm3.constants
import	qm3.io
import	os
import	struct
import	multiprocessing


__prmtop = ""



def coordinates_read( mol, fname = None ):
	f = qm3.io.open_r( fname )
	f.readline()
	n = int( f.readline().strip() )
	if( mol.natm != n ):
		print( "- Wrong atom number... (or broken!)" )
		f.close()
		return
	n *= 3
	mol.coor = []
	while( len( mol.coor ) < n ):
		mol.coor += [ float( j ) for j in f.readline().split() ]
	mol.boxl = [ float( j ) for j in f.readline().split()[0:3] ]
	qm3.io.close( f, fname )


def coordinates_write( mol, fname = None ):
	f = qm3.io.open_w( fname )
	f.write( "comment\n%d\n"%( mol.natm ) )
	n = 3 * mol.natm
	for i in range( n ):
		f.write( "%12.7lf"%( mol.coor[i] ) )
		if( (i+1)%6 == 0 ):
			f.write( "\n" )
	if( n%6 != 0 ):
		f.write( "\n" )
	f.write( "%12.7lf%12.7lf%12.7lf%12.7lf%12.7lf%12.7lf\n"%( mol.boxl[0], mol.boxl[1], mol.boxl[2], 90.0, 90.0, 90.0 ) )
	qm3.io.close( f, fname )


def topology_read( mol, fname = None ):
	global	__prmtop
	__prmtop = ""
	f = qm3.io.open_r( fname )
	if( mol.natm > 0 ):
		l = f.readline()
		while( l ):
			if( l[0:12].upper() == "%FLAG CHARGE" ):
				f.readline()
				mol.chrg = []
				while( len( mol.chrg ) < mol.natm ):
					mol.chrg += [ float( i ) / 18.2223 for i in f.readline().strip().split() ]
			else:
				__prmtop += l
			l = f.readline()
	else:
		mol.initialize()
		nres = 0
		lres = []
		l = f.readline()
		while( l ):
			if( l[0:12].upper() == "%FLAG CHARGE" ):
				f.readline()
				while( len( mol.chrg ) < mol.natm ):
					mol.chrg += [ float( i ) / 18.2223 for i in f.readline().strip().split() ]
			elif( l[0:14].upper() == "%FLAG POINTERS" ):
				__prmtop += l
				__prmtop += f.readline()
				l = f.readline()
				__prmtop += l
				mol.natm = int( l.strip().split()[0] )
				l = f.readline()
				__prmtop += l
				nres = int( l.strip().split()[1] )
			elif( l[0:15].upper() == "%FLAG ATOM_NAME" ):
				__prmtop += l
				__prmtop += f.readline()
				while( len( mol.labl ) < mol.natm ):
					l = f.readline()
					__prmtop += l
					mol.labl += l.strip().split()
			elif( l[0:19].upper() == "%FLAG ATOMIC_NUMBER" ):
				__prmtop += l
				__prmtop += f.readline()
				while( len( mol.anum ) < mol.natm ):
					l = f.readline()
					__prmtop += l
					mol.anum += [ int( i ) for i in l.strip().split() ]
			elif( l[0:10].upper() == "%FLAG MASS" ):
				__prmtop += l
				__prmtop += f.readline()
				while( len( mol.mass ) < mol.natm ):
					l = f.readline()
					__prmtop += l
					mol.mass += [ float( i ) for i in l.strip().split() ]
			elif( l[0:19].upper() == "%FLAG RESIDUE_LABEL" ):
				__prmtop += l
				__prmtop += f.readline()
				while( len( lres ) < nres ):
					l = f.readline()
					__prmtop += l
					lres += l.strip().split()
			elif( l[0:21].upper() == "%FLAG RESIDUE_POINTER" ):
				__prmtop += l
				__prmtop += f.readline()
				while( len( mol.res_lim ) < nres ):
					l = f.readline()
					__prmtop += l
					mol.res_lim += [ int( i ) - 1 for i in l.strip().split() ]
				mol.res_lim.append( mol.natm )
			elif( l[0:20].upper() == "%FLAG BOX_DIMENSIONS" ):
				__prmtop += l
				__prmtop += f.readline()
				l = f.readline()
				__prmtop += l
				mol.boxl = [ float( i ) for i in l.strip().split()[-3:] ]
			else:
				__prmtop += l
			l = f.readline()
		mol.segn = [ "A" for i in range( mol.natm ) ]
		mol.seg_lim = [ 0, mol.natm ]
		for i in range( nres ):
			for j in range( mol.res_lim[i], mol.res_lim[i+1] ):
				mol.resi.append( i + 1 )
				mol.resn.append( lres[i] )
	qm3.io.close( f, fname )


def topology_write( mol, fname = None ):
	global	__prmtop
	if( __prmtop and mol.chrg ):
		f = qm3.io.open_w( fname )
		f.write( __prmtop )
		f.write( "%FLAG CHARGE                                                                    \n" )
		f.write( "%FORMAT(5E16.8)                                                                 \n" )
		for i in range( mol.natm ):
			f.write( "%16.8lf"%( mol.chrg[i] * 18.2223 ) )
			if( (i+1)%5 == 0 ):
				f.write( "\n" )
		qm3.io.close( f, fname )



try:
	import	qm3.engines._sander
	class py_sander( qm3.engines._sander.sander ):
		def __init__( self, mol, prmtop, cutoff = 10.0, PBC = True, qmsel = None, method = "AM1", charge = 0 ):
			qm3.engines._sander.sander.__init__( self, mol, prmtop, cutoff, PBC, qmsel, method, charge )
except:
	pass



class sander( object ):

	def __init__( self ):
		self.exe = "bash r.sander"


	def update_coor( self, mol ):
		coordinates_write( mol, "inpcrd" )


	def update_chrg( self, mol ):
		topology_write( mol, "prmtop" )


	def get_func( self, mol ):
		self.update_coor( mol )
		os.system( self.exe )
		f = open( "mden", "rt" )
		l = f.readline()
		k = True
		while( l and k ):
			if( l[0:2] == "L6" ):
				try:
					mol.func += qm3.constants.K2J * float( l.strip().split()[2] )
					k = False
				except:
					pass
			l = f.readline()
		f.close()


	def get_grad( self, mol ):
		self.update_coor( mol )
		self.get_func( mol )
		# dirty NETCDF wrapper...
		f = open( "mdfrc", "rb" )
		b = f.read()
		f.close()
		o_x = ord( "x" ); o_y = ord( "y" ); o_z = ord( "z" )
		n = len( b )
		w = 0
		f = True
		while( w < n-2 and f ):
			# bytes (python3) / str (python2) compliant...
			if( ( b[w] == o_x and b[w+1] == o_y and b[w+2] == o_z ) or ( b[w] == "x" and b[w+1] == "y" and b[w+2] == "z" ) ):
				w += 7
				f  = False
			w += 1
		g = struct.unpack( ">%df"%( mol.natm * 3 ), b[w:] )
		for i in range( mol.natm ):
			i3 = i * 3
			for j in [0, 1, 2]:
				mol.grad[i3+j] -= g[i3+j] * qm3.constants.K2J




# ------------------------------------------------------------------------------------
# - ORCA_FAKE - ("orca" in current folder, for UNMODIFIED sander versions...)
# ------------------------------------------------------------------------------------
##!/usr/bin/env python
#
#f = open( "inpfile.xyz", "rt" )
#nQM = int( f.readline().strip() )
#f.close()
#f = open( "orc_job.engrad", "wt" )
#f.write( """#
## Number of atoms
##
# %d
##
## The current total energy in Eh
##
#      0.000000000000
##
## The current gradient in Eh/bohr
##
#"""%( nQM ) )
#for i in range( nQM ):
#	f.write( "       0.000000000000\n       0.000000000000\n       0.000000000000\n" )
#f.close()
#
#f = open( "ptchrg.xyz", "rt" )
#nMM = int( f.readline().strip() )
#f.close()
#f = open( "orc_job.pcgrad", "wt" )
#f.write( "%d\n"%( nMM ) )
#for i in range( nMM ):
#	f.write( "   0.000000000000   0.000000000000   0.000000000000\n" )
#f.close()
#
#f = open( "orc_job.dat", "wt" )
#f.write( "FINAL SINGLE POINT ENERGY         0.000000000000\n" )
#f.write( "Total Dipole Moment    :      0.00000       0.00000       0.00000\n" )
#f.write( "                        -----------------------------------------\n" )
#f.write( "Magnitude (a.u.)       :      0.00000\n" )
#f.close()
# ------------------------------------------------------------------------------------

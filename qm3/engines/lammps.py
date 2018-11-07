# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.constants
import	qm3.elements
import	qm3.io
import	qm3.mol
import	os, stat
import	time
import	struct



def lammps_read( fname = None ):
	mol = qm3.mol.molecule()
	f = qm3.io.open_r( fname )
	l = f.readline().lower()
	q = True
	m = {}
	while( l != "" and q ):
		if( l.find( "xlo" ) > 0 and l.find( "xhi" ) > 0 ):
			t = l.strip().split()
			mol.boxl[0] = float( t[1] ) - float( t[0] )
		elif( l.find( "ylo" ) > 0 and l.find( "yhi" ) > 0 ):
			t = l.strip().split()
			mol.boxl[1] = float( t[1] ) - float( t[0] )
		elif( l.find( "zlo" ) > 0 and l.find( "zhi" ) > 0 ):
			t = l.strip().split()
			mol.boxl[2] = float( t[1] ) - float( t[0] )
		elif( l.find( "masses" ) >= 0 ):
			f.readline()
			l = f.readline().strip()
			while( l != "" ):
				t = l.split()
				m[t[0]] = float( t[1] )
				l = f.readline().strip()
		elif( l.find( "atoms" ) == 0 ):
			q = False
			f.readline()
			l = f.readline().strip()
			while( l != "" ):
				t = l.split()
				mol.segn.append( "XXX" )
				mol.resn.append( "X" )
				mol.resi.append( int( t[1] ) )
				mol.labl.append( t[2] )
				mol.mass.append( m[t[2]] )
				mol.chrg.append( float( t[3] ) )
				mol.coor += [ float( t[4] ), float( t[5] ), float( t[6] ) ]
				mol.natm += 1
				l = f.readline().strip()
		l = f.readline().lower()
	qm3.io.close( f, fname )
	mol.settle()
	return( mol )


def guess_labels( mol ):
	for i in range( mol.natm ):
		mol.labl[i] = qm3.elements.symbol[ sorted( [ ( math.fabs( qm3.elements.mass[j] - mol.mass[i] ), j ) for j in iter( qm3.elements.mass ) if j > 0 ] )[0][1] ]




sys.path.insert( 0, "/Users/smarti/Devel/lammps" )
try:
	import	lammps as _lammps
	class py_lammps( object ):
		def __init__( self, inp, name = "serial", cmdargs = [ "-sc", "none" ] ):
			self.lmp = _lammps.lammps( name, cmdargs )
			self.lmp.file( inp )
			self.chg = self.lmp.gather_atoms( "q", 1, 1 )
			self.crd = self.lmp.gather_atoms( "x", 1, 3 )


		def stop( self ):
			self.lmp.close()


		def update_charges( self, mol ):
			for i in range( mol.natm ):
				self.chg[i] = mol.chrg[i]
			self.lmp.scatter_atoms( "q", 1, 3, self.chg )


		def update_coor( self, mol ):
			for i in range( 3 * mol.natm ):
				self.crd[i] = mol.coor[i]
			self.lmp.scatter_atoms( "x", 1, 3, self.crd )


		def get_func( self, mol ):
			self.update_coor( mol )
			self.lmp.command( "run 0" )
			mol.func += self.lmp.get_thermo( "pe" ) * qm3.constants.K2J


		def get_grad( self, mol ):
			self.get_func( mol )
			frz = self.lmp.gather_atoms( "f", 1, 3 )
			for i in range( 3 * mol.natm ):
				mol.grad[i] -= frz[i] * qm3.constants.K2J


except:
	pass




class lammps_pipe( object ):

	def __init__( self, inp ):
		self.pfd = open( "lammps.pipe", "wt" )
		f = open( inp, "rt" )
		self.pfd.write( f.read() + "\n" )
		self.pfd.flush()
		f.close()


	def stop( self ):
		self.pfd.flush()
		self.pfd.close()


	def update_charges( self, mol ):
		f = open( "lammps.chrg", "wt" )
		f.write( "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\n"%( mol.natm ) )
		f.write( "ITEM: BOX BOUNDS pp pp pp\n0.0000000000000000e+00 0.0000000000000000e+00\n0.0000000000000000e+00 0.0000000000000000e+00\n0.0000000000000000e+00 0.0000000000000000e+00\n" )
		f.write( "ITEM: ATOMS id q\n" )
		for i in range( mol.natm ):
			f.write( "%d %.6lf\n"%( i+1, mol.chrg[i] ) )
		f.close()
		self.pfd.write( "read_dump lammps.chrg 0 q box no format native\n" )
		self.pfd.flush()


	def update_coor( self, mol ):
		f = open( "lammps.coor", "wt" )
		f.write( "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\n"%( mol.natm ) )
		f.write( "ITEM: BOX BOUNDS pp pp pp\n0.0000000000000000e+00 0.0000000000000000e+00\n0.0000000000000000e+00 0.0000000000000000e+00\n0.0000000000000000e+00 0.0000000000000000e+00\n" )
		f.write( "ITEM: ATOMS id x y z\n" )
		for i in range( mol.natm ):
			i3 = i * 3
			f.write( "%d %.10lf %.10lf %.10lf\n"%( i+1, mol.coor[i3], mol.coor[i3+1], mol.coor[i3+2] ) )
		f.close()
		self.pfd.write( "read_dump lammps.coor 0 x y z box no format native\n" )
		self.pfd.flush()


	def get_func( self, mol ):
# ------------------------------------------------
#		self.pfd.write( "reset_timestep 0\nrun 0\n" )
#		self.pfd.write( "print $(pe) file lammps.ener screen no\n" )
#		self.pfd.write( "print \"done\" file lammps.lock screen no\n" )
#		self.pfd.flush()
#		while( not os.access( "lammps.lock", os.R_OK ) ):
#			time.sleep( 0.01 )
#		while( os.stat( "lammps.lock" )[stat.ST_SIZE] < 5 ):
#			time.sleep( 0.01 )
#		os.unlink( "lammps.lock" )
#		f = open( "lammps.ener", "rt" )
#		mol.func += float( f.read() ) * qm3.constants.K2J
#		f.close()
# patched (verlet.cpp) -------------------------
		self.pfd.write( "reset_timestep 0\nrun 0\n" )
		self.pfd.flush()
		while( not os.access( "lammps.ener", os.R_OK ) ):
			time.sleep( 0.01 )
		while( os.stat( "lammps.ener" )[stat.ST_SIZE] < 8 ):
			time.sleep( 0.01 )
		f = open( "lammps.ener", "rb" )
		mol.func += struct.unpack( "d", f.read( 8 ) )[0] * qm3.constants.K2J
		f.close()
		os.unlink( "lammps.ener" )
# ------------------------------------------------


	def get_grad( self, mol ):
# ------------------------------------------------
#		self.pfd.write( "reset_timestep 0\nrun 0\n" )
#		self.pfd.write( "print $(pe) file lammps.ener screen no\n" )
#		self.pfd.write( "write_dump all custom lammps.force id fx fy fz modify sort id format line \"%d %.10lf %.10lf %.10lf\"\n" )
#		self.pfd.write( "print \"done\" file lammps.lock screen no\n" )
#		self.pfd.flush()
#		while( not os.access( "lammps.lock", os.R_OK ) ):
#			time.sleep( 0.01 )
#		while( os.stat( "lammps.lock" )[stat.ST_SIZE] < 5 ):
#			time.sleep( 0.01 )
#		os.unlink( "lammps.lock" )
#		f = open( "lammps.ener", "rt" )
#		mol.func += float( f.read() ) * qm3.constants.K2J
#		f.close()
#		f = open( "lammps.force", "rt" )
#		l = f.readline()
#		while( l.strip() != "ITEM: ATOMS id fx fy fz" ):
#			l = f.readline()
#		for i in range( mol.natm ):
#			t = f.readline().strip().split()
#			for j in [0, 1, 2]:
#				mol.grad[3*i+j] -= float( t[j+1] ) * qm3.constants.K2J
#		f.close()
# patched (verlet.cpp) -------------------------
		self.get_func( mol )
		while( not os.access( "lammps.force", os.R_OK ) ):
			time.sleep( 0.01 )
		while( os.stat( "lammps.force" )[stat.ST_SIZE] < mol.natm * 28 ):
			time.sleep( 0.01 )
		f = open( "lammps.force", "rb" )
		for i in range( mol.natm ):
			w = 3 * ( struct.unpack( "i", f.read( 4 ) )[0] - 1 )
			t = struct.unpack( "3d", f.read( 24 ) )
			for j in [0, 1, 2]:
				mol.grad[w+j] -= t[j] * qm3.constants.K2J
		f.close()
		os.unlink( "lammps.force" )
# ------------------------------------------------




class lammps( object ):

	def __init__( self ):
		self.exe = "mpirun -n 4 /Users/smarti/Devel/lammps/lmp_mpi-16Mar18 -in lammps.inp -sc none -log lammps.log"


	def update_coor( self, mol ):
		f = open( "lammps.xyzq", "wt" )
		f.write( "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\n"%( mol.natm ) )
		f.write( "ITEM: BOX BOUNDS pp pp pp\n0.0000000000000000e+00 0.0000000000000000e+00\n0.0000000000000000e+00 0.0000000000000000e+00\n0.0000000000000000e+00 0.0000000000000000e+00\n" )
		f.write( "ITEM: ATOMS id x y z q\n" )
		for i in range( mol.natm ):
			i3 = i * 3
			f.write( "%d %.10lf %.10lf %.10lf %.6lf\n"%( i+1, mol.coor[i3], mol.coor[i3+1], mol.coor[i3+2], mol.chrg[i] ) )
		f.close()


	def update_charges( self, mol ):
		self.update_coor( mol )


	def get_func( self, mol ):
		os.system( self.exe )
		f = open( "lammps.ener", "rt" )
		mol.func += float( f.read() ) * qm3.constants.K2J
		f.close()


	def get_grad( self, mol ):
		self.get_func( mol )
		f = open( "lammps.force", "rt" )
		l = f.readline()
		while( l.strip() != "ITEM: ATOMS id fx fy fz" ):
			l = f.readline()
		for i in range( mol.natm ):
			i3 = i * 3
			t = f.readline().strip().split()
			for j in [0, 1, 2]:
				mol.grad[i3+j] -= float( t[j+1] ) * qm3.constants.K2J
		f.close()





LAMMPS_INP = """############################################################
units           real
atom_style      full

pair_style      buck/coul/long 11.0
kspace_style    pppm 1.e-4
bond_style      harmonic
angle_style     harmonic
dihedral_style  charmm
improper_style  harmonic

read_data       data
read_dump       lammps.xyzq 0 x y z q box no format native

pair_modify     tail yes
neighbor        2.0 bin
neigh_modify    delay 5
timestep        1.0
thermo          1
thermo_style    multi

reset_timestep  0
run             0
print           $(pe) file lammps.ener screen no
write_dump      all custom lammps.force id fx fy fz modify sort id format line "%d %.10lf %.10lf %.10lf"
############################################################
"""


"""

QM atoms whould be defined within LAMMPS:
	  i) without charges in the data file
	 ii) group qmatm id/molecule @@@
	iii) neigh_modify exclude group qmatm qmatm

in this way there's no need for a _qmmm.* fix, since the energy and the gradient
arising from the lennard-jones are calcualted by LAMMPS

"""

#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	time
import	qm3.mol
import	qm3.problem
import	qm3.constants
import	socket
import	qm3.engines.sqm
import	qm3.engines.namd
import	qm3.engines._qmmm
import	qm3.actions.minimize
try:
	import	cStringIO as io
except:
	import	io


os.environ["QM3_LIBSQM"] = "/Users/smarti/Devel/amber/18/qm3/libsqm.so"


class envi_problem( qm3.problem.template ):
	def __init__( self, molec ):
		qm3.problem.template.__init__( self )

		self.mol = molec

		self.sel = list( range( 4, self.mol.natm ) )

		os.system( "bash r.namd &" )
		while( not os.path.isfile( "namd.shmid" ) ):
			time.sleep( 1 )
		time.sleep( 4 )

		self.eng = qm3.engines.namd.namd_shm()

		self.size = 3 * len( self.sel )
		self.coor = []
		for i in self.sel:
			self.coor += self.mol.coor[3*i:3*i+3]


	def update_coor( self ):
		for i in range( self.size // 3 ):
			for j in [0, 1, 2]:
				self.mol.coor[3*self.sel[i]+j] = self.coor[3*i+j]


	def get_func( self ):
		self.update_coor()
		self.mol.func = 0.0
		self.eng.get_func( self.mol )
		self.func = self.mol.func


	def get_grad( self ):
		self.update_coor()
		self.mol.func = 0.0
		self.mol.grad = [ 0.0 for i in range( 3 * self.mol.natm ) ]
		self.eng.get_grad( self.mol )
		self.func = self.mol.func
		self.grad = []
		for i in self.sel:
			self.grad += self.mol.grad[3*i:3*i+3]


class core_problem( qm3.problem.template ):
	def __init__( self ):
		qm3.problem.template.__init__( self )

		self.mol = qm3.mol.molecule( "pdb" )
		self.mol.psf_read( "psf" )
		self.mol.guess_atomic_numbers()
		self.mol.nbnd_read( "pes/non_bonded" )
		self.mol.boxl = [ 20., 20., 20. ]

		self.sel = [ 0, 1, 2, 3 ]

		f = io.StringIO( """_slave_
&qmmm
maxcyc    = 0,
qm_theory = "PM6",
qmcharge  = 0,
qmmm_int  = 1,
verbosity = 4
 /
qm3_atoms
""" )
		f.seek( 0 )
		s_mm = list( range( 4, self.mol.natm ) )
		self.eng = qm3.engines.sqm.dl_sqm( self.mol, f, self.sel, s_mm, [] )
		f.close()

		self.fix = qm3.engines._qmmm.Int_QMLJ( self.mol, self.sel, s_mm, [] )

		self.size = 3 * len( self.sel )
		self.coor = []
		self.mass = []
		self.anum = []
		for i in self.sel:
			self.coor += self.mol.coor[3*i:3*i+3]
			self.mass.append( self.mol.mass[i] )
			self.anum.append( self.mol.anum[i] )

		self.envi = envi_problem( self.mol )


	def update_coor( self ):
		for i in range( self.size // 3 ):
			for j in [0, 1, 2]:
				self.mol.coor[3*self.sel[i]+j] = self.coor[3*i+j]


	def get_func( self ):
		self.update_coor()
		self.mol.func = 0.0
		self.eng.get_func( self.mol )
		self.fix.get_func( self.mol )
		self.func = self.mol.func


	def get_grad( self ):
		self.update_coor()
		self.mol.func = 0.0
		self.mol.grad = [ 0.0 for i in range( 3 * self.mol.natm ) ]
		self.eng.get_grad( self.mol )
		self.fix.get_func( self.mol )
		self.fix.get_grad( self.mol )
		self.func = self.mol.func
		self.grad = []
		for i in self.sel:
			self.grad += self.mol.grad[3*i:3*i+3]


	def get_hess( self ):
# only when optimizing...
#		self.relax_envi()
		self.num_hess()


	def relax_envi( self ):
		self.get_func()
		self.envi.eng.update_chrg( self.mol )
		qm3.actions.minimize.fire( self.envi, print_frequency = 1, gradient_tolerance = 1, step_number = 1000 )
		self.mol.pdb_write( "mm_opt.pdb" )


	def fake_fchk( self, cmd, dst, update = True ):
		if( update ):
			# ------- Update coordiantes ----
			f = open( "%s.gjf"%( dst ), "rt" )
			f.readline()
			for i in range( self.size // 3 ):
				t = [ float( j ) for j in f.readline().strip().split()[1:] ]
				for j in [0, 1, 2]:
					self.coor[3*i+j] = t[j]
			f.close()
			self.relax_envi()
		# -----------------------------------------------------------------
		if( cmd == "force" ):
			self.get_grad()
		elif( cmd == "freq=noraman" ):
			self.get_hess()
		f = open( "%s.fchk"%( dst ), "wt" )
		f.write( "Number of atoms                            I%17d\n"%( self.size // 3 ) )
		f.write( "Charge                                     I%17d\n"%( 0 ) )
		f.write( "Multiplicity                               I%17d\n"%( 1 ) )
		e = self.func / qm3.constants.H2J
		f.write( "Total Energy                               R%27.15lE\n"%( e ) )
		f.write( "Atomic numbers                             I   N=%12d\n"%( self.size // 3 ) )
		t = []
		for i in range( self.size // 3 ):
			t.append( "%12d"%( self.anum[i] ) )
			if( (i+1)%6 == 0 ):
				t.append( "\n" )
		f.write( "".join( t ) )
		if( (self.size//3)%6 != 0 ):
			f.write( "\n" )
		f.write( "Current cartesian coordinates              R   N=%12d\n"%( self.size ) )
		t = []
		for i in range( self.size ):
			t.append( "%16.8lE"%( self.coor[i] / qm3.constants.A0 ) )
			if( (i+1)%5 == 0 ):
				t.append( "\n" )
		f.write( "".join( t ) )
		if( self.size%5 != 0 ):
			f.write( "\n" )
		f.write( "Real atomic weights                        R   N=%12d\n"%( self.size // 3 ) )
		t = []
		for i in range( self.size // 3 ):
			t.append( "%16.8lE"%( self.mass[i] ) )
			if( (i+1)%5 == 0 ):
				t.append( "\n" )
		f.write( "".join( t ) )
		if( (self.size//3)%5 != 0 ):
			f.write( "\n" )
		f.write( "Cartesian Gradient                         R   N=%12d\n"%( self.size ) )
		c = qm3.constants.A0 / qm3.constants.H2J
		t = []
		for i in range( self.size ):
			t.append( "%16.8lE"%( self.grad[i] * c ) )
			if( (i+1)%5 == 0 ):
				t.append( "\n" )
		f.write( "".join( t ) )
		if( self.size%5 != 0 ):
			f.write( "\n" )
		if( cmd == "freq=noraman" ):
			n = self.size * ( self.size + 1 ) // 2
			f.write( "Cartesian Force Constants                  R   N=%12d\n"%( n ) )
			c = qm3.constants.A0 * qm3.constants.A0 / qm3.constants.H2J
			t = []
			k = 0
			for i in range( self.size ):
				for j in range( i + 1 ):
					t.append( "%16.8lE"%( self.hess[i+self.size*j] * c ) )
					if( (k+1)%5 == 0 ):
						t.append( "\n" )
					k += 1
			f.write( "".join( t ) )
			if( n%5 != 0 ):
				f.write( "\n" )
		f.close()



obj = core_problem()

#qm3.actions.minimize.baker( obj, print_frequency = 1, gradient_tolerance = 1, step_number = 1000, follow_mode = 0 )
#qm3.engines.namd.pdb_write( obj.mol, "t.pdb", [ 0, 1, 2, 3 ] )
#obj.fake_fchk( "freq=noraman", "t", False )

#qm3.actions.minimize.baker( obj, print_frequency = 1, gradient_tolerance = 1, step_number = 1000 )
#qm3.engines.namd.pdb_write( obj.mol, "r.pdb", [ 0, 1, 2, 3 ] )
#obj.fake_fchk( "freq=noraman", "r", False )

if( True ):
	sck = socket.socket( socket.AF_UNIX, socket.SOCK_STREAM )
	unx = "/tmp/qm3_unix"
	os.system( "/bin/rm -vf " + unx )
	sck.bind( unx )
	while( True ):
		sck.listen( 1 )
		son, adr = sck.accept()
		try:
			buf = ""
			while( len( buf ) < 1024 ):
				buf += son.recv( 1024 )
			cmd, fld = buf.strip().split()
			obj.fake_fchk( cmd, fld, True )
			son.sendall( "done" )
		except:
			son.sendall( "fail" )
		son.close()
	sck.close()

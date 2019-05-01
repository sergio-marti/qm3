# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	qm3.utils
import	qm3.engines



def dftb_input( obj, mol, run, det = "Yes" ):
	f = open( "dftb_in.hsd", "wt" )
	f.write( """Driver = {}
Geometry = GenFormat {
  %d C
  %s
"""%( len( obj.sel ) + len( obj.lnk ), str.join( " ", obj.tbl ) ) )
	j = 0
	for i in obj.sel:
		i3 = i * 3
		f.write( "  %4d%4d%20.10lf%20.10lf%20.10lf\n"%( j + 1, obj.tbl.index( obj.smb[j] ) + 1,
			mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
			mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
			mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
		j += 1
	if( obj.lnk ):
		obj.vla = []
		k = len( obj.sel )
		w = obj.tbl.index( "H" ) + 1
		for i,j in obj.lnk:
			c, v = qm3.utils.LA_coordinates( i, j, mol )
			f.write( "  %4d%4d%20.10lf%20.10lf%20.10lf\n"%( k + 1, w, c[0], c[1], c[2] ) )
			obj.vla.append( ( obj.sel.index( i ), k, [-v[0], -v[1], -v[2]] ) )
			k += 1
	f.write( "}\n" )
	f.write( """Hamiltonian = DFTB {
  SCC = Yes
  MaxSCCIterations = 1000
  SlaterKosterFiles = Type2FileNames {
    Prefix = \"%s\"
    Separator = \"-\"
    Suffix = \".skf\"
  }
  MaxAngularMomentum {
"""%( obj.prm ) )
	for e in obj.tbl:
		f.write( "    %s = \"%s\"\n"%( e, obj.ang[e] ) )
	f.write( "  }\n  Charge = %d\n"%( obj.chg ) )
	if( obj.opt ):
		f.write( "\n" + obj.opt + "\n" )
	if( os.access( "charges.bin", os.R_OK ) ):
		f.write( "  ReadInitialCharges = Yes\n" )
	if( obj.nbn ):
		g = open( "charges.dat", "wt" )
		for i in obj.nbn:
			i3 = i * 3
			g.write( "%20.10lf%20.10lf%20.10lf%12.6lf\n"%(
				mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
				mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
				mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), mol.chrg[i] ) )
		g.close()
		f.write( """  ElectricField = {
    PointCharges = {
      CoordsAndCharges [Angstrom] = DirectRead {
        Records = %d
        File = \"charges.dat\"
      }
    }
  }
"""%( len( obj.nbn ) ) )
	f.write( """}
Options { WriteDetailedOut = %s }
Analysis {
  MullikenAnalysis = Yes
  WriteBandOut = No
"""%( det ) )
	if( run == "grad" ):
		f.write( "  CalculateForces = Yes\n" )
	else:
		f.write( "  CalculateForces = No\n" )
	f.write( """}
ParserOptions { WriteHSDInput = No }
""" )
	f.close()



class dftb( qm3.engines.qmbase ):

	def __init__( self, mol, sele, nbnd = [], link = [], chrg = 0, parm = "", xopt = None ):
		qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
		self.chg = chrg
		self.prm = parm
		# extra non-default options (such as ThirdOrderFull, HubbardDerivs, HCorrection, Dispersion & so on...)
		self.opt = xopt
		self.exe = "bash r.dftb"
		self.smb = [ i.title() for i in self.smb ]
		self.tbl = { i:None for i in self.smb }
		if( self.lnk ):
			self.tbl["H"] = None
		self.tbl = list( self.tbl )
		self.ang = {}
		tmp = { 1: "s", 2: "p", 3: "d", 4: "f", 5: "g" }
		for elm in self.tbl:
			f = open( parm + "%s-%s.skf"%( elm, elm ), "rt" )
			self.ang[elm] = tmp[int( f.readline().strip().split()[-1] )]
			f.close()


	def mk_input( self, mol, run ):
		dftb_input( self, mol, run )


	def parse_log( self, mol, run ):
		f = open( "detailed.out", "rt" )
		l = f.readline()
		while( l ):
			if( l.strip() == "Net atomic charges (e)" or l.strip() == "Atomic gross charges (e)" ):
				f.readline()
				for i in range( len( self.sel ) ):
					mol.chrg[self.sel[i]] = float( f.readline().split()[1] )
			if( l[0:20].strip() == "Total energy:" ):
				mol.func += float( l.split()[2] ) * self._ce
			if( l.strip() == "Total Forces" and run == "grad" ):
				g = []
				for i in range( len( self.sel ) + len( self.vla ) ):
					g += [ - float( j ) * self._cg for j in f.readline().split()[1:] ]
				qm3.utils.LA_gradient( self.vla, g )
				for i in range( len( self.sel ) ):
					i3 = i * 3
					for j in [0, 1, 2]:
						mol.grad[3*self.sel[i]+j] += g[i3+j]
			if( l.strip() == "Forces on external charges" and run == "grad" ):
				for i in range( len( self.nbn ) ):
					g = [ - float( j ) * self._cg for j in f.readline().split() ]
					for j in [0, 1, 2]:
						mol.grad[3*self.nbn[i]+j] += g[j]
			l = f.readline()
		f.close()



try:
	import	ctypes
	class dl_dftb( qm3.engines.qmbase ):
	
		def __init__( self, mol, sele, nbnd = [], link = [], chrg = 0, parm = "", xopt = None ):
			qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
			self.chg = chrg
			self.prm = parm
			# extra non-default options (such as ThirdOrderFull, HubbardDerivs, HCorrection, Dispersion & so on...)
			self.opt = xopt
			self.exe = "bash r.dftb"
			self.smb = [ i.title() for i in self.smb ]
			self.tbl = { i:None for i in self.smb }
			if( self.lnk ):
				self.tbl["H"] = None
			self.tbl = list( self.tbl )
			self.ang = {}
			tmp = { 1: "s", 2: "p", 3: "d", 4: "f", 5: "g" }
			for elm in self.tbl:
				f = open( parm + "%s-%s.skf"%( elm, elm ), "rt" )
				self.ang[elm] = tmp[int( f.readline().strip().split()[-1] )]
				f.close()

			self.nQM = len( self.sel ) + len( self.lnk )
			self.nMM = len( self.nbn )
			self.siz = 1 + 3 * ( self.nQM + self.nMM ) + self.nMM + self.nQM
			self.vec = ( ctypes.c_double * self.siz )()
			self.lib = ctypes.CDLL( os.getenv( "QM3_LIBDFTB" ) )
			self.lib.qm3_dftb_calc_.argtypes = [ ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_double ) ]
			self.lib.qm3_dftb_calc_.restype = None

			dftb_input( self, mol, "grad", "No" )
			self.lib.qm3_dftb_init_()
	
	
		def update_coor( self, mol ):
			for i in range( len( self.sel ) ):
				i3 = self.sel[i] * 3
				j3 = i * 3
				for j in [0, 1, 2]:
					self.vec[j3+j] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) 
			self.vla = []
			k = len( self.sel )
			for i in range( len( self.lnk ) ):
				j3 = k * 3
				c, v = qm3.utils.LA_coordinates( self.lnk[i][0], self.lnk[i][1], mol )
				for j in [0, 1, 2]:
					self.vec[j3+j] = c[j]
				self.vla.append( ( self.sel.index( self.lnk[i][0] ), k, v[:] ) )
				k += 1
			for i in range( self.nMM ):
				i3 = self.nbn[i] * 3
				j3 = ( self.nQM + i ) * 3
				for j in [0, 1, 2]:
					self.vec[j3+j] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) 
			for i in range( self.nMM ):
				self.vec[3*(self.nQM+self.nMM)+i] = mol.chrg[self.nbn[i]]


		def get_func( self, mol ):
			self.update_coor( mol )
			self.lib.qm3_dftb_calc_( ctypes.c_int( self.nQM ), ctypes.c_int( self.nMM ), ctypes.c_int( self.siz ), self.vec )
			mol.func += self.vec[0] * self._ce
			for i in range( len( self.sel ) ):
				mol.chrg[self.sel[i]] = self.vec[i+1]
	
	
		def get_grad( self, mol ):
			self.update_coor( mol )
			self.lib.qm3_dftb_calc_( ctypes.c_int( self.nQM ), ctypes.c_int( self.nMM ), ctypes.c_int( self.siz ), self.vec )
			mol.func += self.vec[0] * self._ce
			for i in range( len( self.sel ) ):
				mol.chrg[self.sel[i]] = self.vec[i+1]
			g = [ j * self._cg for j in self.vec[self.nQM+1:] ]
			qm3.utils.LA_gradient( self.vla, g )
			for i in range( len( self.sel ) ):
				i3 = i * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.sel[i]+j] += g[i3+j]
			for i in range( self.nMM ):
				i3 = ( self.nQM + i ) * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.nbn[i]+j] += g[i3+j]

except:
	pass

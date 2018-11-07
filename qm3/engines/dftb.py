# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	qm3.utils
import	qm3.engines



class dftb( qm3.engines.qmbase ):

	def __init__( self, mol, sele, nbnd = [], link = [], hami = None ):
		qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
		self.chg = 0
		self.prm = ""
		self.exe = "/Users/smarti/Devel/dftb+/dftb+_1.3.vecLib > dftb_in.log"
		self.smb = [ i.title() for i in self.smb ]
		self.tbl = { i:None for i in self.smb }
		if( self.lnk ):
			self.tbl["H"] = None
		self.tbl = list( self.tbl )
		self.ang = { "H": "s", "C": "p", "N": "p", "O": "p", "P": "d", "S": "d" }
		# extra non-default options (such as ThirdOrderFull, HubbardDerivs & so on...)
		self.ham = hami


	def mk_input( self, mol, run ):
		f = open( "dftb_in.hsd", "wt" )
		f.write( "Driver = {}\nGeometry = GenFormat {\n  %d C\n  %s\n"%( len( self.sel ) + len( self.lnk ), str.join( " ", self.tbl ) ) )
		j = 0
		for i in self.sel:
			i3 = i * 3
			f.write( "  %4d%4d%20.10lf%20.10lf%20.10lf\n"%( j + 1, self.tbl.index( self.smb[j] ) + 1,
				mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
				mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
				mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
			j += 1
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			w = self.tbl.index( "H" ) + 1
			for i,j in self.lnk:
				c, v = qm3.utils.LA_coordinates( i, j, mol )
				f.write( "  %4d%4d%20.10lf%20.10lf%20.10lf\n"%( k + 1, w, c[0], c[1], c[2] ) )
				self.vla.append( ( self.sel.index( i ), k, [-v[0], -v[1], -v[2]] ) )
				k += 1
		f.write( "}\n" )
		f.write( "Hamiltonian = DFTB {\n  SCC = Yes\n  MaxSCCIterations = 1000\n  SlaterKosterFiles = Type2FileNames {\n    Prefix = \"%s\"\n    Separator = \"-\"\n    Suffix = \".skf\"\n  }\n  MaxAngularMomentum {\n"%( self.prm ) )
		for e in self.tbl:
			f.write( "    %s = \"%s\"\n"%( e, self.ang[e] ) )
		f.write( "  }\n  Charge = %d\n"%( self.chg ) )
		if( self.ham ):
			f.write( "\n" + self.ham + "\n" )
		if( os.access( "charges.bin", os.R_OK ) ):
			f.write( "  ReadInitialCharges = Yes\n" )
		if( self.nbn ):
			g = open( "charges.dat", "wt" )
			for i in self.nbn:
				i3 = i * 3
				g.write( "%20.10lf%20.10lf%20.10lf%12.6lf\n"%(
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), mol.chrg[i] ) )
			g.close()
			f.write( "  ElectricField = {\n    PointCharges = {\n      CoordsAndCharges [Angstrom] = DirectRead {\n        Records = %d\n        File = \"charges.dat\"\n      }\n    }\n  }\n"%( len( self.nbn ) ) )
		f.write( "}\n" )
		f.write( "Options { WriteDetailedOut = Yes }\nAnalysis {\n  MullikenAnalysis = Yes\n  CalculateForces = Yes\n  WriteBandOut = No\n}\nParserOptions { ParserVersion = 5 }\n" )
		f.close()


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
					g += [ - float( j ) * self._cg for j in f.readline().split() ]
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

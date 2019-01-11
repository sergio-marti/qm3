# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	math
import	qm3.utils
import	qm3.engines



class nwchem( qm3.engines.qmbase ):

	def __init__( self, mol, ini, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
		self.ini = ini
		self.exe = "bash r.nwchem"


	def mk_input( self, mol, run ):
		f = open( "nwchem.nw", "wt" )
		f.write( "start nwchem\ngeometry units angstroms nocenter noautoz noautosym\n" )
		j = 0
		for i in self.sel:
			i3 = i * 3
			f.write( "%4s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j],
				mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
				mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
				mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
			j += 1
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			for i,j in self.lnk:
				c, v = qm3.utils.LA_coordinates( i, j, mol )
				f.write( "%-4s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] ) )
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1
		f.write( "end\n" )
		if( os.access( "nwchem.movecs", os.R_OK ) ):
			f.write( self.ini.replace( "@@@", "vectors input nwchem.movecs" ) )
		else:
			f.write( self.ini.replace( "@@@\n", "" ) )
		if( self.nbn ):
			f.write( "set bq:max_nbq %d\n"%( len( self.nbn ) + 1 ) )
			f.write( "bq units angstroms\n force nwchem.mmgrad\n load nwchem.mmchrg units angstroms format 1 2 3 4\nend\n" )
			g = open( "nwchem.mmchrg", "wt" )
			for i in self.nbn:
				i3 = i * 3
				g.write( "%20.10lf%20.10lf%20.10lf%12.4lf\n"%(
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), mol.chrg[i] ) )
			g.close()
		who = "scf"
		if( self.ini.lower().find( "dft" ) > -1 ):
			who = "dft"
		if( run == "grad" ):
			f.write( "task " + who + " gradient\n" )
		else:
			f.write( "task " + who + "\n" )
		f.close()


	def parse_log( self, mol, run ):
		f = open( "nwchem.log", "rt" )
		l = f.readline()
		while( l ):
			L = l.strip()
			# read energy
			if( L.find( "Total " ) > -1 and L.find( " energy " ) > -1 ):
				t = L.split()
				if( len( t ) == 5 ):
					mol.func += float( t[4] ) * self._ce
			# read gradient and LAs projection
			if( run == "grad" and L.find( "ENERGY GRADIENTS" ) > -1 ):
				f.readline(); f.readline(); f.readline()
				g = []
				for i in range( len( self.sel ) + len( self.lnk ) ):
					g += [ float( j ) * self._cg for j in f.readline().strip().split()[-3:] ]
				qm3.utils.LA_gradient( self.vla, g )
				for i in range( len( self.sel ) ):
					i3 = i * 3
					for j in [0, 1, 2]:
						mol.grad[3*self.sel[i]+j] += g[i3+j]
			l = f.readline()
		f.close()
		if( self.nbn and os.access( "nwchem.mmgrad", os.R_OK ) ):
			f = open( "nwchem.mmgrad", "rt" )
			f.readline()
			for i in range( len( self.nbn ) ):
				t = [ float( j ) for j in f.readline().strip().split() ]
				mol.grad[3*self.nbn[i]]   += t[0] * self._cg
				mol.grad[3*self.nbn[i]+1] += t[1] * self._cg
				mol.grad[3*self.nbn[i]+2] += t[2] * self._cg
			f.close()


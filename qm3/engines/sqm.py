# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	math
import	qm3.constants
import	qm3.elements
import	struct
import	qm3.utils
import	qm3.engines



class sqm( qm3.engines.qmbase ):

	def __init__( self, mol, ini, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
		self.ini = ini
		self.exe = "AMBERHOME=/Users/smarti/Devel/amber/16 /Users/smarti/Devel/amber/16/bin/sqm"


	def mk_input( self, mol, run ):
		f = open( "mdin", "wt" )
		f.write( "single point energy calculation\n" )
		f.write( "&qmmm\n" + self.ini + "maxcyc = 0,\nqmmm_int = 1,\nverbosity = 4\n /\n" )
		j = 0
		for i in self.sel:
			i3 = i * 3
			f.write( "%3d%4s%20.10lf%20.10lf%20.10lf\n"%( qm3.elements.rsymbol[self.smb[j]], self.smb[j],
				mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
				mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
				mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
			j += 1
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			for i,j in self.lnk:
				c, v = qm3.utils.LA_coordinates( i, j, mol )
				# To allow the interaction of the Link-Atom with the environment change atomic number to "1"
				f.write( "%3d%4s%20.10lf%20.10lf%20.10lf\n"%( -1, "H", c[0], c[1], c[2] ) )
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1
		f.write( "\n" )
		if( self.nbn ):
			f.write( "#EXCHARGES\n" )
			for i in self.nbn:
				i3 = i * 3
				f.write( "%3d%4s%20.10lf%20.10lf%20.10lf%12.4lf\n"%( 1, "H",
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), mol.chrg[i] ) )
			f.write( "#END\n" )
		f.close()


	def parse_log( self, mol, run ):
		f = open( "mm_output", "rb" )
		f.read( 4 )
		mol.func += struct.unpack( "d", f.read( 8 ) )[0] * qm3.constants.K2J
		f.read( 4 )
		if( run == "grad" ):
			n = struct.unpack( "i", f.read( 4 ) )[0] // 8
			g = [ i * qm3.constants.K2J for i in struct.unpack( "%dd"%( n ), f.read( 8 * n ) ) ]
			f.read( 4 )
			qm3.utils.LA_gradient( self.vla, g )
			for i in range( len( self.sel ) ):
				i3 = i * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.sel[i]+j] += g[i3+j]
			if( self.nbn ):
				n = struct.unpack( "i", f.read( 4 ) )[0] // 8
				g = [ i * qm3.constants.K2J for i in struct.unpack( "%dd"%( n ), f.read( 8 * n ) ) ]
				for i in range( len( self.nbn ) ):
					i3 = i * 3
					for j in [0, 1, 2]:
						mol.grad[3*self.nbn[i]+j] += g[i3+j]
		f.close()


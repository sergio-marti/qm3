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


def sqm_input( obj, mol, run ):
	f = open( "mdin", "wt" )
	f.write( "single point energy calculation\n" )
	f.write( "&qmmm\n" + obj.ini + "maxcyc = 0,\nqmmm_int = 1,\nverbosity = 4\n /\n" )
	j = 0
	for i in obj.sel:
		i3 = i * 3
		f.write( "%3d%4s%20.10lf%20.10lf%20.10lf\n"%( qm3.elements.rsymbol[obj.smb[j]], obj.smb[j],
			mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
			mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
			mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
		j += 1
	if( obj.lnk ):
		obj.vla = []
		k = len( obj.sel )
		for i,j in obj.lnk:
			c, v = qm3.utils.LA_coordinates( i, j, mol )
			# To allow the interaction of the Link-Atom with the environment change atomic number to "1"
			f.write( "%3d%4s%20.10lf%20.10lf%20.10lf\n"%( -1, "H", c[0], c[1], c[2] ) )
			obj.vla.append( ( obj.sel.index( i ), k, v[:] ) )
			k += 1
	f.write( "\n" )
	if( obj.nbn ):
		f.write( "#EXCHARGES\n" )
		for i in obj.nbn:
			i3 = i * 3
			f.write( "%3d%4s%20.10lf%20.10lf%20.10lf%12.4lf\n"%( 1, "H",
				mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
				mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
				mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), mol.chrg[i] ) )
		f.write( "#END\n" )
	f.close()



class sqm( qm3.engines.qmbase ):

	def __init__( self, mol, ini, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
		self.ini = ini
		self.exe = "bash r.sqm"


	def mk_input( self, mol, run ):
		sqm_input( self, mol, run )


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



try:
	import	ctypes
	class dl_sqm( qm3.engines.qmbase ):
	
		def __init__( self, mol, ini, sele, nbnd = [], link = [] ):
			qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )

			self.nQM = len( self.sel ) + len( self.lnk )
			self.siz = 1 + 3 * ( self.nQM + len( self.nbn ) ) + self.nQM
			self.vec = ( ctypes.c_double * self.siz )()
			self.lib = ctypes.CDLL( os.getenv( "QM3_LIBSQM" ) )
			self.lib.qm3_sqm_calc_.argtypes = [ ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( ctypes.c_double ) ]
			self.lib.qm3_sqm_calc_.restype = None

			sqm_input( self, mol, ini )
			self.lib.qm3_sqm_init_()
	
	
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
			for i in range( len( self.nbn ) ):
				i3 = self.nbn[i] * 3
				j3 = ( self.nQM + i ) * 3
				for j in [0, 1, 2]:
					self.vec[j3+j] = mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) 


		def get_func( self, mol ):
			self.update_coor( mol )
			self.lib.qm3_sqm_calc_( ctypes.c_int( self.siz ), self.vec )
			mol.func += self.vec[0] * qm3.constants.K2J
			for i in range( len( self.sel ) ):
				mol.chrg[self.sel[i]] = self.vec[i+1]
	
	
		def get_grad( self, mol ):
			self.update_coor( mol )
			self.lib.qm3_sqm_calc_( ctypes.c_int( self.siz ), self.vec )
			mol.func += self.vec[0] * qm3.constants.K2J
			for i in range( len( self.sel ) ):
				mol.chrg[self.sel[i]] = self.vec[i+1]
			g = [ j * qm3.constants.K2J for j in self.vec[self.nQM+1:] ]
			qm3.utils.LA_gradient( self.vla, g )
			for i in range( len( self.sel ) ):
				i3 = i * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.sel[i]+j] += g[i3+j]
			for i in range( len( self.nbn ) ):
				i3 = ( self.nQM + i ) * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.nbn[i]+j] += g[i3+j]

except:
	pass

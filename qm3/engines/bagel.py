# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os, os.path
import	re
import	math
import	json
import	qm3.utils
import	qm3.engines



class bagel( qm3.engines.qmbase ):

	def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
		self.exe = "bash r.bagel"


	def mk_input( self, mol, run ):
		atm = []
		j = 0
		for i in self.sel:
			i3 = i * 3
			atm.append( { "atom": self.smb[j], "xyz": [
				mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
				mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
				mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 )
				 ] } )
			j += 1
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			for i,j in self.lnk:
				c, v = qm3.utils.LA_coordinates( i, j, mol )
				atm.append( { "atom": "H", "xyz": [ c[0], c[1], c[2] ] } )
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1
		if( self.nbn ):
			for i in self.nbn:
				i3 = i * 3
				atm.append( { "atom": "Q", "xyz": [
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 )
					 ], "charge": mol.chrg[i] } )
		s_wf = ""
		if( os.path.isdir( "bagel.archive" ) ):
			s_wf = "{ \"title\" : \"save_ref\", \"file\" : \"bagel\" },"
		f = open( "bagel.json", "wt" )
		buf = self.inp.replace( "qm3_atoms", json.dumps( atm, indent = 0 ) )
		buf = buf.replace( "qm3_guess", s_wf )
		f.write( buf )
		f.close()


	def parse_log( self, mol, run ):
		if( os.access( "ENERGY.out", os.R_OK ) ):
			f = open( "ENERGY.out", "rt" )
			mol.func += float( f.read().strip() ) * self._ce
			f.close()
			os.unlink( "ENERGY.out" )
		if( run == "grad" and os.access( "FORCE_0.out", os.R_OK ) ):
			f = open( "FORCE_0.out", "rt" )
			f.readline()
			g = []
			for i in range( len( self.sel ) + len( self.lnk ) ):
				g += [ float( j ) * self._cg for j in f.readline().strip().split()[1:] ]
			if( self.nbn ):
				for i in self.nbn:
					i3 = i * 3
					t = f.readline().strip().split()
					for j in [0, 1, 2]:
						mol.grad[i3+j] += float( t[j+1] ) * self._cg
			f.close()
			qm3.utils.LA_gradient( self.vla, g )
			# copy array
			for i in range( len( self.sel ) ):
				i3 = i * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.sel[i]+j] += g[i3+j]
			os.unlink( "FORCE_0.out" )


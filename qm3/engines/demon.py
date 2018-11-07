# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	qm3.utils
import	qm3.engines



class demon( qm3.engines.qmbase ):

	def __init__( self, mol, ini, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
		self.ini = ini
		self.exe = "/Users/smarti/Devel/deMon2k/4/deMon_4.4.1"


	def mk_input( self, mol, run ):
		f = open( "deMon.inp", "wt" )
		f.write( "title slave\nsymmetry off\nqm/mm charmm\nvisualization off\n" )
		f.write( self.ini )
		if( os.access( "deMon.rst", os.R_OK ) ):
			f.write( "guess restart\n" )
		if( self.nbn ):
			f.write( "embed file\n" )
		f.write( "geometry cartesian angstrom\n" )
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
		f.close()
		if( self.nbn ):
			f = open( "deMon.cub", "wt" )
			for i in self.nbn:
				i3 = i * 3
				f.write( "%20.10lf%20.10lf%20.10lf%12.4lf\n"%(
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), mol.chrg[i] ) )
			f.close()


	def parse_log( self, mol, run ):
		f = open( "deMon.qmm", "rt" )
		mol.func += float( f.readline().split()[-1] ) * self._ce
		if( run == "grad" ):
			l = f.readline()
			while( l ):
				# read gradient and LAs projection
				if( l.strip() == "QMFORCES" ):
					g = []
					for i in range( len( self.sel ) + len( self.lnk ) ):
						g += [ float( j ) * self._cg for j in f.readline().strip().split()[-3:] ]
					qm3.utils.LA_gradient( self.vla, g )
					for i in range( len( self.sel ) ):
						i3 = i * 3
						for j in [0, 1, 2]:
							mol.grad[3*self.sel[i]+j] += g[i3+j]
				# read MM gradient
				if( self.nbn and l.strip() == "EMBEDFORCES" ):
					for i in range( len( self.nbn ) ):
						t = [ float( j ) * self._cg for j in f.readline().strip().split() ]
						mol.grad[3*self.nbn[i]]   += t[0]
						mol.grad[3*self.nbn[i]+1] += t[1]
						mol.grad[3*self.nbn[i]+2] += t[2]
				l = f.readline()
		f.close()

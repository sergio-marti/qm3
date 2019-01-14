# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	re
import	os
import	glob
import	math
import	qm3.utils
import	qm3.engines



class lsdalton( qm3.engines.qmbase ):

	def __init__( self, mol, ini, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
		self.ini = ini
		self.bas = "6-31G*"
		self.chg = 0
		self.exe = "bash r.lsdalton"


	def mk_input( self, mol, run ):
		f = open( "LSDALTON.INP", "wt" )
		if( run == "grad" ):
			f.write( self.ini.replace( "@@@", "**RESPONS\n*MOLGRA" ) )
		else:
			f.write( self.ini.replace( "@@@", "" ) )
		f.close()
		f = open( "MOLECULE.INP", "wt" )
		n = len( self.sel ) + len( self.lnk ) + len( self.nbn )
		q = self.chg + sum( [ mol.chrg[i] for i in self.nbn ] )
		f.write( "ATOMBASIS\n.\n.\nAtomtypes=%d Charge=%.1lf Angstrom Nosymmetry\n"%( n, q ) )
		j = 0
		for i in self.sel:
			i3 = i * 3
			f.write( "Charge=%.1lf Atoms=1 Basis=%s\n%-2s%20.10lf%20.10lf%20.10lf\n"%( 
				mol.anum[i], self.bas, self.smb[j], 
				mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
				mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
				mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
			j += 1
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			for i,j in self.lnk:
				c, v = qm3.utils.LA_coordinates( i, j, mol )
				f.write( "Charge=1.0 Atoms=1 Basis=%s\n%-2s%20.10lf%20.10lf%20.10lf\n"%( 
				self.bas, "H", c[0], c[1], c[2] ) )
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1
		if( self.nbn ):
			for i in self.nbn:
				i3 = i * 3
				f.write( "Charge=%.3lf Atoms=1 pointcharge\nbq%20.10lf%20.10lf%20.10lf\n"%( 
					mol.chrg[i],
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3] / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
			f.close()


	def parse_log( self, mol, run ):
		f = open( "LSDALTON.OUT", "rt" )
		e = re.compile( "[\ ]+Final[\ ][^\ ]*[\ ]energy:[\ ]+([0-9\.\-]+)" )
		l = f.readline()
		while( l != "" ):
			if( e.match( l ) ):
				mol.func += self._ce * float( e.findall( l )[0] )
			if( run == "grad" and l.find( "Molecular gradient (au)" ) > 0 ):
				f.readline(); f.readline()
				g = []
				for i in range( len( self.sel ) + len( self.lnk ) + len( self.nbn ) ):
					g += [ float( j ) * self._cg for j in f.readline().strip().split()[1:] ]
				qm3.utils.LA_gradient( self.vla, g )
				for i in range( len( self.sel ) ):
					i3 = i * 3
					for j in [0, 1, 2]:
						mol.grad[3*self.sel[i]+j] += g[i3+j]
				n = len( self.sel ) + len( self.lnk )
				for i in range( len( self.nbn ) ):
					i3 = ( n + i ) * 3
					for j in [0, 1, 2]:
						mol.grad[3*self.nbn[i]+j] += g[i3+j]
			l = f.readline()
		f.close()


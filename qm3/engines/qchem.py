# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os, os.path
import	re
import	math
import	qm3.utils
import	qm3.engines



class qchem( qm3.engines.qmbase ):

	def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
		self.exe = "bash r.qchem"


	def mk_input( self, mol, run ):
		s_qm = ""
		j = 0
		for i in self.sel:
			i3 = i * 3
			s_qm += "%2s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j],
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) )
			j += 1
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			for i,j in self.lnk:
				c, v = qm3.utils.LA_coordinates( i, j, mol )
				s_qm += "%-2s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] )
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1
		s_mm = ""
		if( self.nbn ):
			s_mm = "$external_charges\n"
			for i in self.nbn:
				i3 = i * 3
				s_mm += "%20.10lf%20.10lf%20.10lf%12.4lf\n"%( 
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), 
					mol.chrg[i] )
			s_mm += "$end"
		s_rn = "single_point"
		if( run == "grad" ):
			s_rn = "force"
		s_wf = ""
		if( os.path.isdir( "qchem.tmp" ) ):
			s_wf = "scf_guess read"
		f = open( "qchem.inp", "wt" )
		buf = self.inp.replace( "qm3_atoms", s_qm[:-1] )
		buf = buf.replace( "qm3_job", s_rn )
		buf = buf.replace( "qm3_guess", s_wf )
		buf = buf.replace( "qm3_charges", s_mm )
		f.write( buf )
		f.close()


	def parse_log( self, mol, run ):
		f = open( "qchem.log", "rt" )
		mol.func += float( re.compile( "The QM part of the Energy is[\ ]+([0-9\.\-]+)" ).findall( f.read() )[0] ) * self._ce
		f.close()
		if( run == "grad" ):
			f = open( "efield.dat", "rt" )
			if( self.nbn ):
				for i in self.nbn:
					i3 = i * 3
					t = f.readline().strip().split()
					for j in [0, 1, 2]:
						mol.grad[i3+j] += - self._cg * mol.chrg[i] * float( t[j] )
			g = []
			for i in range( len( self.sel ) + len( self.lnk ) ):
				g += [ float( j ) * self._cg for j in f.readline().strip().split() ]
			f.close()
			qm3.utils.LA_gradient( self.vla, g )
			# copy array
			for i in range( len( self.sel ) ):
				i3 = i * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.sel[i]+j] += g[i3+j]




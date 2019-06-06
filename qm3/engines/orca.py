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



class orca( qm3.engines.qmbase ):

	def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
		self.exe = "bash r.orca"


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
			s_mm = "%pointcharges \"orca.pc\""
			f = open( "orca.pc", "wt" )
			f.write( "%d\n"%( len( self.nbn ) ) )
			for i in self.nbn:
				i3 = i * 3
				f.write( "%12.4lf%20.10lf%20.10lf%20.10lf\n"%( mol.chrg[i],
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3] / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
			f.close()
		s_rn = ""
		if( run == "grad" ):
			s_rn = "engrad"
		f = open( "orca.inp", "wt" )
		buf = self.inp.replace( "qm3_atoms", s_qm[:-1] )
		buf = buf.replace( "qm3_job", s_rn )
		buf = buf.replace( "qm3_charges", s_mm )
		f.write( buf )
		f.close()


	def parse_log( self, mol, run ):
		if( run == "grad" ):
			f = open( "orca.engrad", "rt" )
			t = re.compile( "[0-9\.\-]+" ).findall( f.read() )
			f.close()
			n = int( t[0] )
			mol.func += float( t[1] ) * self._ce
			g = [ float( t[i] ) * self._cg for i in range( 2, 2 + n * 3 ) ]
			qm3.utils.LA_gradient( self.vla, g )
			for i in range( len( self.sel ) ):
				i3 = i * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.sel[i]+j] += g[i3+j]
			if( self.nbn and os.access( "orca.pcgrad", os.R_OK ) ):
				f = open( "orca.pcgrad", "rt" )
				t = re.compile( "[0-9\.\-]+" ).findall( f.read() )
				f.close()
				n = int( t[0] )
				g = [ float( t[i] ) * self._cg for i in range( 1, 1 + n * 3 ) ]
				for i in range( len( self.nbn ) ):
					i3 = i * 3
					for j in [0, 1, 2]:
						mol.grad[3*self.nbn[i]+j] += g[i3+j]
		else:
			f = open( "orca.out", "rt" )
			mol.func += self._ce * float( re.compile( "FINAL SINGLE POINT ENERGY[\ ]*([0-9\.\-]+)" ).findall( f.read() )[0] )
			f.close()
		for ff in glob.glob( "orca.*" ):
			if( ff != "orca.gbw" and ff != "orca.ges" ):
				os.unlink( ff )





"""
~/Devel/orca/3.0.2/bin/orca_chelpg orca.gbw		>> standard output


[...]

CHELPG Charges            
--------------------------------
  0   O   :      -0.751399
  1   H   :       0.402272
  2   H   :       0.394430
  3   O   :      -0.798244
  4   H   :       0.376724
  5   H   :       0.376218
--------------------------------
Total charge:    -0.000000
--------------------------------

CHELPG charges calculated...

[...]



Automatic calculation of CHELPG charges using the default values can also be achieved by specifying
! CHELPG
in the simple input section.

"""

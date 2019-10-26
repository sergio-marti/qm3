# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	math
import	struct
import	qm3.elements
import	qm3.engines



class gamess( qm3.engines.qmbase ):

	def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
		self.exe = "bash r.gamess"
		self.chk = None


	def mk_input( self, mol, run ):
		s_qm = ""
		j = 0
		for i in self.sel:
			i3 = i * 3
			s_qm += "%4s%6.1lf%20.10lf%20.10lf%20.10lf\n"%( self.smb[j], float( qm3.elements.rsymbol[self.smb[j]] ),
				mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
				mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
				mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) )
			j += 1
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			for i,j in self.lnk:
				c, v = qm3.engines.LA_coordinates( i, j, mol )
				s_qm += "%-4s%6.1lf%20.10lf%20.10lf%20.10lf\n"%( "H", 1., c[0], c[1], c[2] )
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1
		if( run == "grad" ):
			s_rn = "gradient"
		elif( run == "hess" ):
			s_rn = "hessian"
		else:
			s_rn = "energy"
		s_wf = ""
		if( self.chk ):
			s_wf = " $guess\nguess=moread\n $end\n%s"%( self.chk )
		f = open( "gamess.inp", "wt" )
		buf = self.inp.replace( "qm3_atoms", s_qm[:-1] )
		buf = buf.replace( "qm3_guess", s_wf )
		buf = buf.replace( "qm3_job", s_rn )
		f.write( buf )
		f.close()
		if( self.nbn ):
			f = open( "mm_charges", "wb" )
#			f.write( struct.pack( "i", 4 ) + struct.pack( "i", len( self.nbn ) ) + struct.pack( "i", 4 ) )
			# "-fdefault-integer-8" / "-i8" flag for Gamess-US compilation
			f.write( struct.pack( "i", 8 ) + struct.pack( "l", len( self.nbn ) ) + struct.pack( "i", 8 ) )
			for i in self.nbn:
				i3 = i * 3
				f.write( struct.pack( "i", 32 ) )
				f.write( struct.pack( "d", mol.coor[i3]   - mol.boxl[0] * round( ( mol.coor[i3]   ) / mol.boxl[0], 0 ) ) )
				f.write( struct.pack( "d", mol.coor[i3+1] - mol.boxl[1] * round( ( mol.coor[i3+1] ) / mol.boxl[1], 0 ) ) )
				f.write( struct.pack( "d", mol.coor[i3+2] - mol.boxl[2] * round( ( mol.coor[i3+2] ) / mol.boxl[2], 0 ) ) )
				f.write( struct.pack( "d", mol.chrg[i] ) )
				f.write( struct.pack( "i", 32 ) )
			f.close()


	def parse_log( self, mol, run ):
		self.chk = ""
		f = open( "gamess.data", "rt" )
		l = f.readline()
		while( l ):
			if( l[0:2] == "E(" ):
				mol.func += float( l.split()[-5][:-1] ) * self._ce
			if( l[0:5] == " $VEC" ):
				while( l[0:5] != " $END" ):
					self.chk += l
					l = f.readline()
				self.chk += l
			if( l[0:6] == " $GRAD" and run in [ "grad", "hess" ] ):
				f.readline()
				g = []
				for i in range( len( self.sel ) + len( self.lnk ) ):
					g += [ float( j ) * self._cg for j in f.readline().split()[2:5] ]
				qm3.engines.LA_gradient( self.vla, g )
				for i in range( len( self.sel ) ):
					i3 = i * 3
					for j in [0, 1, 2]:
						mol.grad[3*self.sel[i]+j] += g[i3+j]
			if( l[0:6] == " $HESS" and run == "hess" ):
				t = [ [] for i in range( 3 * ( len( self.sel ) + len( self.lnk ) ) ) ]
				f.readline()
				l = f.readline()
				while( l[0:5] != " $END" ):
					i = int( l[0:2].strip() ) - 1
					t[i] += [ float( l[j*15+5:j*15+20].strip() ) * self._ch for j in range( ( len( l ) - 6 ) // 15 ) ]
					l = f.readline()
				n = 3 * len( self.sel )
				k = 0
				for i in range( n ):
					for j in range( n ):
						mol.hess[i*n+j] += t[i][j]
			l = f.readline()
		f.close()
		f = open( "gamess.out", "rt" )
		l = f.readline()
		while( l and l[0:13] != " NET CHARGES:" ):
			l = f.readline()
		if( l[0:13] == " NET CHARGES:" ):
			f.readline(); f.readline(); f.readline()
			for i in self.sel:
				mol.chrg[i] = float( f.readline().split()[1] )
		f.close()
		if( self.nbn and os.access( "mm_output", os.R_OK ) and run in [ "grad", "hess" ] ):
			f = open( "mm_output", "rb" )
			f.read( struct.unpack( "i", f.read( 4 ) )[0] + 4 )
			f.read( struct.unpack( "i", f.read( 4 ) )[0] + 4 )
			n = len( self.nbn )
			f.read( 4 )
			tx = struct.unpack( "%dd"%( n ), f.read( 8 * n ) )
			f.read( 8 )
			ty = struct.unpack( "%dd"%( n ), f.read( 8 * n ) )
			f.read( 8 )
			tz = struct.unpack( "%dd"%( n ), f.read( 8 * n ) )
			f.close()
			for i in range( n ):
				mol.grad[3*self.nbn[i]]   += tx[i] * self._cg
				mol.grad[3*self.nbn[i]+1] += ty[i] * self._cg
				mol.grad[3*self.nbn[i]+2] += tz[i] * self._cg




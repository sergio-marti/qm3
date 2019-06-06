# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	math
import	socket
import	struct
import	qm3.utils
import	qm3.engines
try:
	import cStringIO as io
except:
	import io



class tchem( qm3.engines.qmbase ):

	def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, inp, sele, nbnd, link )
		self.exe = "bash r.tchem"


	def mk_input( self, mol, run ):
		f = open( "tchem_qm.xyz", "wt" )
		f.write( "%d\n\n"%( len( self.sel ) + len( self.lnk ) ) )
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
				f.write( "%4s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] ) )
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1
		f.close()
		s_mm = ""
		s_rn = "energy"
		if( run == "grad" ):
			s_rn = "gradient"
		if( self.nbn ):
			s_mm = "pointcharges  tchem_mm.xyz\namber         yes\n"
			f = open( "tchem_mm.xyz", "wt" )
			f.write( "%d\n\n"%( len( self.nbn ) ) )
			for i in self.nbn:
				i3 = i * 3
				f.write( "%12.4lf%20.10lf%20.10lf%20.10lf\n"%( mol.chrg[i],
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
			f.close()
		f = open( "tchem.inp", "wt" )
		buf = self.inp.replace( "qm3_job", s_rn )
		buf = buf.replace( "qm3_charges", s_mm )
		f.write( buf )
		f.close()


	def parse_log( self, mol, run ):
		f = open( "tchem.log", "rt" )
		l = f.readline()
		while( l ):
			L = l.strip()
			# read energy
			if( L[0:13] == "FINAL ENERGY:" ):
				mol.func += float( L.split()[2] ) * self._ce
			# read gradient and perform LAs projection
			if( run == "grad" and L[0:5] == "dE/dX" ):
				g = []
				for i in range( len( self.sel ) + len( self.lnk ) ):
					g += [ float( j ) * self._cg for j in f.readline().strip().split()[-3:] ]
				qm3.utils.LA_gradient( self.vla, g )
				for i in range( len( self.sel ) ):
					i3 = i * 3
					for j in [0, 1, 2]:
						mol.grad[3*self.sel[i]+j] += g[i3+j]
			# read MM gradient
			if( run == "grad" and self.nbn and L[0:30] == "------- MM / Point charge part" ):
				for i in range( len( self.nbn ) ):
					t = [ float( j ) * self._cg for j in f.readline().strip().split() ]
					mol.grad[3*self.nbn[i]]   += t[0]
					mol.grad[3*self.nbn[i]+1] += t[1]
					mol.grad[3*self.nbn[i]+2] += t[2]
			l = f.readline()
		f.close()




class tchem_sckt( qm3.engines.qmbase ):

	def __init__( self, unx, mol, sele, nbnd = [], link = [] ):
		f = io.StringIO( "" )
		qm3.engines.qmbase.__init__( self, mol, f, sele, nbnd, link )

		self.ssz = 1000 * 3 * 8
		self.unx = unx
		self.tot = len( sele ) + len( link )
		self.smb = "".join( [ "%-2s"%( i ) for i in self.smb ] + len( link ) * [ "H " ] )


	def connect( self, mol ):
		self.sck = socket.socket( socket.AF_UNIX, socket.SOCK_STREAM )
		self.sck.connect( self.unx )
		# n_qm + n_mm
		self.sck.send( struct.pack( "i", self.tot ) + struct.pack( "i", len( self.nbn ) ) )
		# smb
		if( sys.version_info[0] == 2 ):
			self.sck.send( self.smb )
		else:
			self.sck.send( bytes( self.smb, encoding = "ascii" ) )
		# mm_chg
		msg = b""
		for i in self.nbn:
			msg += struct.pack( "d", mol.chrg[i] )
			if( len( msg ) >= self.ssz ):
				self.sck.send( msg )
				msg = b""
		if( len( msg ) > 0 ):
			self.sck.send( msg )


	def get_grad( self, mol, gradient = True ):
		# setup QM
		crd = []
		for i in self.sel:
			i3 = i * 3
			for j in [0, 1, 2]:
				crd.append( mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) )
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			for i,j in self.lnk:
				c, v = qm3.utils.LA_coordinates( i, j, mol )
				crd += c[:]
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1

		self.connect( mol )

		# send QM coordinates
		msg = b""
		for i in range( self.tot ):
			i3 = i * 3
			for j in [ 0, 1, 2 ]:
				msg += struct.pack( "d", crd[i3+j] )
			if( len( msg ) >= self.ssz ):
				self.sck.send( msg )
				msg = b""
		if( len( msg ) > 0 ):
			self.sck.send( msg )

		# send MM coordinates
		msg = b""
		for i in self.nbn:
			i3 = i * 3
			for j in [ 0, 1, 2 ]:
				msg += struct.pack( "d", mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) )
			if( len( msg ) >= self.ssz ):
				self.sck.send( msg )
				msg = b""
		if( len( msg ) > 0 ):
			self.sck.send( msg )

		# recv energy
		mol.func += self._ce * float( struct.unpack( "d", self.sck.recv( 8 ) )[0] )

		# recv QM charges
		msg = b""
		tot = 8 * self.tot
		cur = tot
		while( cur > 0 ):
			msg += self.sck.recv( min( cur, self.ssz ) )
			cur = tot - len( msg )
		t = struct.unpack( "%dd"%( self.tot ), msg )
		for i in range( len( self.sel ) ):
			mol.chrg[self.sel[i]] = t[i]

		# recv QM + MM gradients
		msg = b""
		siz = 3 * ( self.tot + len( self.nbn ) )
		tot = 8 * siz
		cur = tot
		while( cur > 0 ):
			msg += self.sck.recv( min( cur, self.ssz ) )
			cur = tot - len( msg )
		tg = struct.unpack( "%dd"%( siz ), msg )

		self.sck.close()

		if( gradient ):
			# QM gradients (project out LAs)
			g = []
			for i in range( self.tot ):
				i3 = i * 3
				for j in [0, 1, 2]:
					g.append( tg[i3+j] * self._cg )
			qm3.utils.LA_gradient( self.vla, g )
			for i in range( len( self.sel ) ):
				i3 = i * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.sel[i]+j] += g[i3+j]
			# MM gradients
			k = 3 * self.tot
			for i in range( len( self.nbn ) ):
				i3 = i * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.nbn[i]+j] += tg[k+i3+j] * self._cg


	# mpi hook for terachem only properly calculates gradients... :(
	def get_func( self, mol ):
		self.get_grad( mol, gradient = False )

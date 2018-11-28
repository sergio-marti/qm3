# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	qm3.maths.matrix
import	qm3.utils
import	qm3.engines



#
# use with "prop=(field,read)" for calculating the electrostatic charges gradient
#
class gaussian( qm3.engines.qmbase ):

	def __init__( self, mol, ini, mid, end, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
		self.ini = ini
		self.mid = mid
		self.end = end
		self.exe = ". /Users/smarti/Devel/g09/pgi.imac64/g09.profile; g09 g09.com"


	def mk_input( self, mol, run ):
		f = open( "g09.com", "wt" )
		T = ""
		if( os.access( "g09.chk", os.R_OK ) ):
			T = "guess=(read)"
		if( run == "grad" ):
			T += " force"
		if( run == "hess" ):
			T += " freq=noraman"
		f.write( self.ini % ( T ) )
		j = 0
		for i in self.sel:
			i3 = i * 3
			f.write( "%2s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j],
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
			j += 1
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			for i,j in self.lnk:
				c, v = qm3.utils.LA_coordinates( i, j, mol )
				f.write( "%-2s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] ) )
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1
		if( self.nbn ):
			# "charge" keyword
			f.write( "\n" )
			for i in self.nbn:
				i3 = i * 3
				f.write( "%20.10lf%20.10lf%20.10lf%12.4lf\n"%( 
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), 
					mol.chrg[i] ) )
			if( self.mid != "" ):
				f.write( "\n" )
				f.write( self.mid )
			# "prop=(field,read)" keyword
			f.write( "\n" )
			for i in self.nbn:
				i3 = i * 3
				f.write( "%20.10lf%20.10lf%20.10lf\n"%( 
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
		else:
			if( self.mid != "" ):
				f.write( "\n" )
				f.write( self.mid )
		if( self.end != "" ):
			f.write( "\n" )
			f.write( self.end )
		f.write( "\n\n\n\n\n" )
		f.close()


	def parse_log( self, mol, run ):
		fd = open( "Test.FChk", "rt" )
		ln = fd.readline()
		while( ln ):
			# read energy
			if( ln[0:12] == "Total Energy" ):
				mol.func += float( ln.split()[3] ) * self._ce
			# read gadient
			if( run in [ "grad", "hess" ] and ln[0:18] == "Cartesian Gradient" ):
				i = int( ln.split()[-1] )
				j = int( i // 5 ) + ( i%5 != 0 )
				i = 0
				g = []
				while( i < j ):
					ln = fd.readline()
					for itm in ln.split():
						g.append( float( itm ) * self._cg )
					i += 1
				# remove LAs from gradient
				qm3.utils.LA_gradient( self.vla, g )
				# copy array
				for i in range( len( self.sel ) ):
					i3 = i * 3
					for j in [0, 1, 2]:
						mol.grad[3*self.sel[i]+j] += g[i3+j]
				# read hessian (columns)
				if( run == "hess" ):
					ln = fd.readline()
					i = int( ln.split()[-1] )
					j = int( i // 5 ) + ( i % 5 != 0 )
					i = 0
					h = []
					while( i < j ):
						ln = fd.readline()
						for itm in ln.split():
							h.append( float( itm ) * self._ch )
						i += 1
					# truncate LAs and swap hessian (cols>>rows)
					i = 3 * len( self.sel )
					j = i * ( i + 1 ) // 2
					mol.hess = qm3.maths.matrix.from_upper_diagonal_columns( h[0:j], i )
			# read charges
			if( ln[0:11] == "ESP Charges" ):
				i = int( ln.split()[-1] )
				j = int( i // 5 ) + ( i % 5 != 0 )
				i = 0
				k = 0
				while( i < j ):
					ln = fd.readline()
					for itm in ln.split():
						if( k < len( self.sel ) ):
							mol.chrg[self.sel[k]] = float( itm )
							k += 1
					i += 1
			ln = fd.readline()
		fd.close()
		# remove autoenergy of the charges, and calculate MM gradient
		if( self.nbn ):
			fd = open( "g09.log", "rt" )
			ln = fd.readline()
			fl = True
			while( ln and fl ):
				if( ln[0:29] == " Self energy of the charges =" ):
					fl = False
					mol.func -= float( ln.split()[-2] ) * self._ce
				ln = fd.readline()
			if( run in [ "grad", "hess" ] ):
				fl = True
				while( ln and fl ):
					if( ln.strip() == "Potential          X             Y             Z" ):
						fl = False
						for i in range( 1 + len( self.sel ) + len( self.lnk ) ):
							fd.readline()
						for i in self.nbn:
							i3 = i * 3
							t = fd.readline().split()[2:]
							for j in [0, 1, 2]:
								mol.grad[i3+j] += - self._cg * mol.chrg[i] * float( t[j] )
					ln = fd.readline()
			fd.close()
		# return
		os.unlink( "Test.FChk" )



#
# charges gradient MUST be calculated after using "qmmm.Int_QMLJ_MMEL"
#
class gaussian_MMEL( qm3.engines.qmbase ):

	def __init__( self, mol, ini, mid, end, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
		self.ini = ini
		self.mid = mid
		self.end = end
		self.exe = ". /Users/smarti/Devel/g09/pgi.imac64/g09.profile; g09 g09.com"


	def mk_input( self, mol, run ):
		f = open( "g09.com", "wt" )
		T = ""
		if( os.access( "g09.chk", os.R_OK ) ):
			T = "guess=(read)"
		if( run == "grad" ):
			T += " force"
		if( run == "hess" ):
			T += " freq=noraman"
		f.write( self.ini % ( T ) )
		j = 0
		for i in self.sel:
			i3 = i * 3
			f.write( "%2s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j],
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) )
			j += 1
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			for i,j in self.lnk:
				c, v = qm3.utils.LA_coordinates( i, j, mol )
				f.write( "%-2s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] ) )
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1
		if( self.mid != "" ):
			f.write( "\n" )
			f.write( self.mid )
		if( self.nbn ):
			# "charge" keyword
			f.write( "\n" )
			for i in self.nbn:
				i3 = i * 3
				f.write( "%20.10lf%20.10lf%20.10lf%12.4lf\n"%( 
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ), 
					mol.chrg[i] ) )
		if( self.end != "" ):
			f.write( "\n" )
			f.write( self.end )
		f.write( "\n\n\n\n\n" )
		f.close()


	def parse_log( self, mol, run ):
		fd = open( "Test.FChk", "rt" )
		ln = fd.readline()
		while( ln ):
			# read energy
			if( ln[0:12] == "Total Energy" ):
				mol.func += float( ln.split()[3] ) * self._ce
			# read gadient
			if( run in [ "grad", "hess" ] and ln[0:18] == "Cartesian Gradient" ):
				i = int( ln.split()[-1] )
				j = int( i // 5 ) + ( i % 5 != 0 )
				i = 0
				g = []
				while( i < j ):
					ln = fd.readline()
					for itm in ln.split():
						g.append( float( itm ) * self._cg )
					i += 1
				# remove LAs from gradient
				qm3.utils.LA_gradient( self.vla, g )
				# copy array
				for i in range( len( self.sel ) ):
					i3 = i * 3
					for j in [0, 1, 2]:
						mol.grad[3*self.sel[i]+j] += g[i3+j]
				# read hessian (columns)
				if( run == "hess" ):
					ln = fd.readline()
					i = int( ln.split()[-1] )
					j = int(i/5) + (i%5!=0)
					i = 0
					h = []
					while( i < j ):
						ln = fd.readline()
						for itm in ln.split():
							h.append( float( itm ) * self._ch )
						i += 1
					# truncate LAs and swap hessian (cols>>rows)
					i = 3 * len( self.sel )
					j = i * ( i + 1 ) / 2
					mol.hess = qm3.maths.matrix.from_upper_diagonal_columns( h[0:j], i )
			# read charges
			if( ln[0:11] == "ESP Charges" ):
				i = int( ln.split()[-1] )
				j = int( i // 5 ) + ( i % 5 != 0 )
				i = 0
				k = 0
				while( i < j ):
					ln = fd.readline()
					for itm in ln.split():
						if( k < len( self.sel ) ):
							mol.chrg[self.sel[k]] = float( itm )
							k += 1
					i += 1
			ln = fd.readline()
		fd.close()
		# remove autoenergy of the charges
		if( self.nbn ):
			fd = open( "g09.log", "rt" )
			ln = fd.readline()
			fl = True
			while( ln and fl ):
				if( ln[0:29] == " Self energy of the charges =" ):
					fl = False
					mol.func -= float( ln.split()[-2] ) * self._ce
				ln = fd.readline()
			fd.close()
		# return
		os.unlink( "Test.FChk" )



#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

import	matplotlib.pyplot as plt
import	math
import	os


DIME = 2


def muller_brown( coor ):
	A  = [ -200.0, -100.0, -170.0, 15.0 ]
	a  = [   -1.0,   -1.0,   -6.5,  0.7 ]
	b  = [    0.0,    0.0,   11.0,  0.6 ]
	c  = [  -10.0,  -10.0,   -6.5,  0.7 ]
	xo = [    1.0,    0.0,   -0.5, -1.0 ]
	yo = [    0.0,    0.5,    1.5,  1.0 ]
	func = 0.0
	grad = [ .0, .0 ]
	for i in range( 4 ):
		f = A[i] * math.exp( a[i] * math.pow( coor[0] - xo[i], 2.0 ) + b[i] * ( coor[0] - xo[i] ) * ( coor[1] - yo[i] ) + c[i] * math.pow( coor[1] - yo[i], 2.0 ) )
		func    += f
		grad[0] += f * ( 2.0 * a[i] * ( coor[0] - xo[i] ) + b[i] * ( coor[1] - yo[i] ) )
		grad[1] += f * ( b[i] * ( coor[0] - xo[i] ) + 2.0 * c[i] * ( coor[1] - yo[i] ) )
	return( func, grad )




class PEB:
	def __init__( self, surface, r0, rf, kumb, npt ):
		self.np = npt
		self.r0 = r0
		self.rf = rf
		self.ku = kumb

		self.size = npt * DIME
		self.coor = [ .0 ] * self.size
		diff = [ ( i - j ) / float( self.np + 1 ) for i,j in zip( self.rf, self.r0 ) ]
		self.coor = [ .0 ] * self.size
		for i in range( self.np ):
			for j in [ 0, 1 ]:
				self.coor[DIME*i+j] = self.r0[j] + diff[j] * float( i + 1 )
		self.surf = surface

		# >> precalculate MB points for plotting <<
		n  = 200
		dx = ( 1.5 - -2.0 ) / float( n )
		dy = ( 2.0 - -0.5 ) / float( n )
		self.plot_x = [ -2.0 + dx * i for i in range( n ) ]
		self.plot_y = [ -0.5 + dy * i for i in range( n ) ]
		self.plot_z = []
		for i in range( n ):
			self.plot_z.append( [] )
			for j in range( n ):
				self.plot_z[-1].append( self.surf( [ self.plot_x[j], self.plot_y[i] ] )[0] )
		min_z = -150.
		max_z = +100.
		n = 30
		self.plot_l = [ min_z + ( max_z - min_z ) / float( n ) * i for i in range( n ) ]
		self.plot_c = [ "black" ] * n


	def eval( self ):
		func = 0.0
		grad = [ .0 ] * self.size
		# -- potential
		for i in range( self.np ):
			f, g = self.surf( self.coor[DIME*i:DIME*i+DIME] )
			func += f
			for j in [0, 1]:
				grad[DIME*i+j] += g[j]
		# -- springs
		dr = [ i - j for i,j in zip( self.coor[0:DIME], self.r0 ) ]
		func += self.ku * sum( [ i * i for i in dr ] ) * 0.5
		for j in [0, 1]:
			grad[j] += self.ku * dr[j]
		for i in range( 1, self.np ):
			ii = DIME * ( i - 1 )
			jj = DIME * i
			dr = [ i - j for i,j in zip( self.coor[jj:jj+DIME], self.coor[ii:ii+DIME] ) ]
			func += self.ku * sum( [ i * i for i in dr ] ) * 0.5
			for j in [0, 1]:
				grad[ii+j] -= self.ku * dr[j]
				grad[jj+j] += self.ku * dr[j]
		ii = DIME * ( self.np - 1 )
		dr = [ i - j for i,j in zip( self.rf, self.coor[ii:ii+DIME] ) ]
		func += self.ku * sum( [ i * i for i in dr ] ) * 0.5
		for j in [0, 1]:
			grad[ii+j] -= self.ku * dr[j]
		return( func, grad )


	def plot( self, step = 0 ):
		plt.clf()
		plt.grid( True )
		plt.contour( self.plot_x, self.plot_y, self.plot_z, levels = self.plot_l, colors = self.plot_c )
		plt.plot( [ self.r0[0] ] + [ self.coor[DIME*i]   for i in range( self.np ) ] + [ self.rf[0] ],
				  [ self.r0[1] ] + [ self.coor[DIME*i+1] for i in range( self.np ) ] + [ self.rf[1] ], 'b-o' )
		plt.savefig( "snap.%04d.png"%( step ) )




class NEB( PEB ):
	"""
g_i = g_i^{PES} + g_i^s
~~~~
g_i^{ PES } = nabla V left( r_i right )
~~~~
g_i^s = - k_{ i+1 } ( r_{i+1} - r_i ) + k_i (  r_i - r_{i-1} )
newline
g_i^o = left( g_i^{PES} - g_i^{PES} cdot %tau_{{}parallel{}} %tau_{{}parallel{}} right) 
+ g_i^s cdot %tau_{{}parallel{}} %tau_{{}parallel{}}
~~~~
"where " %tau_{{}parallel{}} " is the unit tangent to the path"
newline
g_i^{NEB} = g_i^o+ f( %fi_i ) left( g_i^s - g_i^s cdot %tau_{{}parallel{}}%tau_{{}parallel{}} right)
~~~~
f( %fi_i ) = 1 over 2 left( 1 + cos left( %pi  cos( %fi_i) right) right)
~~~~
cos( %fi_i ) = { left( r_{i+1} - r_i right) cdot left( r_i - r_{i-1} right) } over{ left lline r_{i+1} - r_i right rline left lline r_i - r_{i-1} right rline }
	"""
	def __init__( self, surface, r0, rf, kumb, npt ):
		PEB.__init__( self, surface, r0, rf, kumb, npt )


	def eval( self ):
		func = 0.0
		grad = [ .0 ] * self.size
		# first node
		f, g = self.surf( self.coor[0:DIME] )
		t1 = [ i - j for i,j in zip( self.coor[0:DIME], self.r0 ) ]
		n1 = math.sqrt( sum( [ i * i for i in t1 ] ) )
		t2 = [ i - j for i,j in zip( self.coor[DIME:2*DIME], self.coor[0:DIME] ) ]
		n2 = math.sqrt( sum( [ i * i for i in t2 ] ) )
		tt = [ i / n1 + j / n2 for i,j in zip( t1, t2 ) ]
		mm = math.sqrt( sum( [ i * i for i in tt ] ) )
		tt = [ i / mm for i in tt ]
		func += f
		func += self.ku * sum( [ i * i for i in t1 ] ) * 0.5
		gt = sum( [ i * j for i,j in zip(  g, tt ) ] )
		s1 = sum( [ i * j for i,j in zip( t1, tt ) ] )
		s2 = sum( [ i * j for i,j in zip( t2, tt ) ] )
		for j in [0, 1]:
			grad[j] += g[j] - gt * tt[j] + self.ku * s1 * tt[j] - self.ku * s2 * tt[j]
		t1 = t2[:]
		n1 = n2
		# midd nodes
		for e in range( 1, self.np - 1 ):
			ii = DIME * e
			f, g = self.surf( self.coor[ii:ii+DIME] )
			jj = DIME * ( e + 1 )
			t2 = [ i - j for i,j in zip( self.coor[jj:jj+DIME], self.coor[ii:ii+DIME] ) ]
			n2 = math.sqrt( sum( [ i * i for i in t2 ] ) )
			tt = [ i / n1 + j / n2 for i,j in zip( t1, t2 ) ]
			mm = math.sqrt( sum( [ i * i for i in tt ] ) )
			tt = [ i / mm for i in tt ]
			func += f
			func += self.ku * sum( [ i * i for i in t1 ] ) * 0.5
			gt = sum( [ i * j for i,j in zip(  g, tt ) ] )
			s1 = sum( [ i * j for i,j in zip( t1, tt ) ] )
			s2 = sum( [ i * j for i,j in zip( t2, tt ) ] )
			for j in [0, 1]:
				grad[ii+j] += g[j] - gt * tt[j] + self.ku * s1 * tt[j] - self.ku * s2 * tt[j]
			t1 = t2[:]
			n1 = n2
		# last node
		ii = DIME * ( self.np - 1 )
		f, g = self.surf( self.coor[ii:ii+DIME] )
		t2 = [ i - j for i,j in zip( self.rf, self.coor[ii:ii+DIME] ) ]
		n2 = math.sqrt( sum( [ i * i for i in t1 ] ) )
		tt = [ i / n1 + j / n2 for i,j in zip( t1, t2 ) ]
		mm = math.sqrt( sum( [ i * i for i in tt ] ) )
		tt = [ i / mm for i in tt ]
		func += f
		func += self.ku * sum( [ i * i for i in t1 ] ) * 0.5
		gt = sum( [ i * j for i,j in zip(  g, tt ) ] )
		s1 = sum( [ i * j for i,j in zip( t1, tt ) ] )
		s2 = sum( [ i * j for i,j in zip( t2, tt ) ] )
		for j in [0, 1]:
			grad[ii+j] += g[j] - gt * tt[j] + self.ku * s1 * tt[j] - self.ku * s2 * tt[j]
		return( func, grad )




def fire( obj, step_size = 0.1, step_number = 1000, gradient_tolerance = 0.1 ):
	nstp = 0
	ssiz = step_size
	alph = 0.1
	velo = [ .0 ] * obj.size
	step = [ .0 ] * obj.size
	func, grad = obj.eval()
	norm = math.sqrt( sum( [ i * i for i in grad ] ) )
	it = 0
	obj.plot( it )
	print( "%10d%20.5lf%20.10lf%20.10lf"%( it, func, norm, ssiz ) )
	while( it < step_number and norm > gradient_tolerance ):
		vsiz = math.sqrt( sum( [ i * i for i in velo ] ) )
		vfac = - sum( [ i * j for i,j in zip(  velo, grad ) ] )
		if( vfac > 0.0 ):
			for i in range( obj.size ):
				velo[i] = ( 1.0 - alph ) * velo[i] - alph * grad[i] / norm * vsiz
			if( nstp > 5 ):
				ssiz  = min( ssiz * 1.1, step_size )
				alph *= 0.99
			nstp += 1
		else:
			velo = [ .0 ] * obj.size
			alph  = 0.1
			ssiz *= 0.5
			nstp  = 0
		for i in range( obj.size ):
			velo[i] -= ssiz * grad[i]
			step[i]  = ssiz * velo[i]
		tmp   = math.sqrt( sum( [ i * i for i in step ] ) )
		if( tmp > ssiz ):
			step = [ i * ssiz / tmp for i in step ]
		for i in range( obj.size ):
			obj.coor[i] += step[i]
		func, grad = obj.eval()
		norm = math.sqrt( sum( [ i * i for i in grad ] ) )
		it += 1
		obj.plot( it )
		print( "%10d%20.5lf%20.10lf%20.10lf"%( it, func, norm, ssiz ) )








r0 = [ -0.5582251892711165, 1.4417247583632917 ]
r1 = [ -0.0500285453045008, 0.4667145765710418 ]
r2 = [  0.6234864227409558, 0.0280395760047727 ]

os.system( "/bin/rm -f snap.????.png" )

#fire( PEB( muller_brown, r0, r2, 1500., 7 ) )
fire( NEB( muller_brown, r0, r2, 10., 20 ), step_number = 200 )

os.system( "ffmpeg -y -r 6 -pattern_type glob -i \"snap.????.png\" -q 0 -c:v mpeg4 -b:v 12M neb.avi > /dev/null 2> /dev/null; /bin/rm -f snap.????.png" )

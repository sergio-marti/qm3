# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.maths.matrix
import	qm3.maths.interpolation
import	qm3.utils._mpi



#
# Chem. Phys. Lett. v446, p182 (2007) [10.1016/j.cplett.2007.08.017]
# J. Comput. Chem. v35, p1672 (2014) [10.1002/jcc.23673]
# J. Phys. Chem. A v121, p9764 (2017) [10.1021/acs.jpca.7b10842]
#

def string_distribute( ncrd, nwin, rcrd, rmet, interpolant = qm3.maths.interpolation.hermite_spline ):
	nc2 = ncrd * ncrd
	arc = [ 0.0 ]
	for i in range( 1, nwin ):
# -----------------------------------------------------------------
		if( rmet == None ):
			arc.append( arc[i-1] + math.sqrt( sum( [ math.pow( rcrd[i*ncrd+j] - rcrd[(i-1)*ncrd+j], 2.0 ) for j in range( ncrd ) ] ) ) )
		else:
			# use the metric tensor for arc length calculation (eq 7 @ 10.1002/jcc.23673)
			tmp = [ rcrd[i*ncrd+j] - rcrd[(i-1)*ncrd+j] for j in range( ncrd ) ]
			mat = qm3.maths.matrix.inverse( [ 0.5 * ( rmet[i*nc2+j] + rmet[(i-1)*nc2+j] ) for j in range( nc2 ) ], ncrd, ncrd )
			mat = qm3.maths.matrix.mult( mat, ncrd, ncrd, tmp, ncrd, 1 )
			arc.append( arc[i-1] + math.sqrt( sum( [ tmp[j] * mat[j] for j in range( ncrd ) ] ) ) )
# -----------------------------------------------------------------
	dcv = arc[-1] / float( nwin - 1.0 )
	ind = [ j * dcv for j in range( nwin ) ]
	out = [ 0.0 for i in range( ncrd * nwin ) ]
	for i in range( ncrd ):
		tmp = [ rcrd[j*ncrd+i] for j in range( nwin ) ]
		out[i] = tmp[0]
		out[(nwin-1)*ncrd+i] = tmp[-1]
		eng = interpolant( arc, tmp )
		for j in range( 1, nwin - 1 ):
			out[j*ncrd+i] = eng.calc( ind[j] )[0]
	return( out, arc[-1] )



def string_integrate( ncrd, nwin, i_from = 0, i_to = -1, interpolant = qm3.maths.interpolation.hermite_spline, kumb = None, spline_derivatives = True ):
	# average metrics
	ncr2 = ncrd * ncrd
	rmet = [ 0.0 for i in range( ncr2 * nwin ) ]
	for i in range( nwin ):
		f = open( "string.%03d.met"%( i ), "rt" )
		k = 0
		n = 0.0
		for l in f:
			k += 1
			if( k >= i_from and ( k <= i_to or i_to == -1 ) ):
				t = [ float( j ) for j in l.strip().split() ]
				for j in range( ncr2 ):
					rmet[i*ncr2+j] += t[j]
				n += 1.0
		f.close()
	rmet = [ rmet[j] / n for j in range( ncr2 * nwin ) ]
	f = open( "string.metrics", "wt" )
	for i in range( nwin ):
		f.write( "".join( [ "%20.10lf"%( rmet[i*ncr2+j] ) for j in range( ncr2 ) ] ) + "\n" )
	f.close()
	# average string collective variables
	rcrd = [ 0.0 for i in range( ncrd * nwin ) ]
	f = open( "string.dat", "rt" )
	k = 0
	n = 0.0
	for l in f:
		k += 1
		if( k >= i_from and ( k <= i_to or i_to == -1 ) ):
			t = [ float( j ) for j in l.strip().split() ]
			for j in range( ncrd * nwin ):
				rcrd[j] += t[j]
			n += 1.0
	f.close()
	rcrd = [ rcrd[j] / n for j in range( ncrd * nwin ) ]
	f = open( "string.avevars", "wt" )
	for i in range( nwin ):
		f.write( "".join( [ "%20.10lf"%( rcrd[i*ncrd+j] ) for j in range( ncrd ) ] ) + "\n" )
	f.close()
	rcrd = string_distribute( ncrd, nwin, rcrd, rmet, interpolant )[0]
	f = open( "string.colvars", "wt" )
	for i in range( nwin ):
		f.write( "".join( [ "%20.10lf"%( rcrd[i*ncrd+j] ) for j in range( ncrd ) ] ) + "\n" )
	f.close()
	dFdz = [ 0.0 for i in range( ncrd * nwin ) ]
	if( kumb == None or len( kumb ) != ncrd ):
		# average forces from datas (eq. 21 @ 10.1016/j.cplett.2007.08.017)
		for i in range( nwin ):
			f = open( "string.%03d.frc"%( i ), "rt" )
			k = 0
			n = 0.0
			for l in f:
				k += 1
				if( k >= i_from and ( k <= i_to or i_to == -1 ) ):
					t = [ float( j ) for j in l.strip().split() ]
					for j in range( ncrd ):
						dFdz[i*ncrd+j] += t[j]
					n += 1.0
			f.close()
	else:
		# evaluate forces from datas (eq. 21 @ 10.1016/j.cplett.2007.08.017)
		for i in range( nwin ):
			f = open( "string.%03d.cvs"%( i ), "rt" )
			k = 0
			n = 0.0
			for l in f:
				k += 1
				if( k >= i_from and ( k <= i_to or i_to == -1 ) ):
					t = [ float( j ) for j in l.strip().split() ]
					for j in range( ncrd ):
						dFdz[i*ncrd+j] += kumb[j] * ( rcrd[i*ncrd+j] - t[j] )
					n += 1.0
			f.close()
	dFdz = [ dFdz[j] / n for j in range( ncrd * nwin ) ]
	f = open( "string.dFdz", "wt" )
	for i in range( nwin ):
		for j in range( ncrd ):
			f.write(  "%20.10lf"%( dFdz[i*ncrd+j] )  )
		f.write( "\n" )
	f.close()
# -------------------------------------------------------------------------------
	if( spline_derivatives ):
		# spline interpolated collective variables derivative
		dzds = [ 0.0 for j in range( ncrd * nwin ) ]
		for i in range( ncrd ):
			tmp = [ rcrd[j*ncrd+i] for j in range( nwin ) ]
			eng = interpolant( range( nwin ), tmp )
			for j in range( nwin ):
				dzds[j*ncrd+i] = eng.calc( j )[1]
# -------------------------------------------------------------------------------
	else:
		# linear interpolation of the collective variables derivative (ds = 1)
		dzds = [ 0.0 for j in range( ncrd * nwin ) ]
		for j in range( ncrd ):
			dzds[j] = ( rcrd[ncrd+j] - rcrd[j] )
			dzds[(nwin-1)*ncrd+j] = ( rcrd[(nwin-1)*ncrd+j] - rcrd[(nwin-2)*ncrd+j] )
		for i in range( 1, nwin - 1 ):
			for j in range( ncrd ):
				dzds[i*ncrd+j] = ( rcrd[(i+1)*ncrd+j] - rcrd[(i-1)*ncrd+j] ) * 0.5
# -------------------------------------------------------------------------------
	f = open( "string.dzds", "wt" )
	for i in range( nwin ):
		for j in range( ncrd ):
			f.write(  "%20.10lf"%( dzds[i*ncrd+j] )  )
		f.write( "\n" )
	f.close()
	# thermodynamic integration (eq. 20 @ 10.1016/j.cplett.2007.08.017)
	mfep = [ 0.0 for j in range( nwin ) ]
	for i in range( 1, nwin ):
		mfep[i] = mfep[i-1] + sum( [ dzds[i*ncrd+j] * dFdz[i*ncrd+j] for j in range( ncrd ) ] )
	f = open( "string.mfep", "wt" )
	f.write( "\n".join( [ "%20.10lf"%( mfep[j] ) for j in range( nwin ) ] ) )
	f.close()



class string( object ):

	def __init__( self, node, conf, molec ):
		"""
String config:
------------------------------------------------------------------------
ncrd      nwin      tstp|1e-5
dist      atom_i    atom_j    kumb      min_val|.0      max_val|9.e99
...
dist      atom_i    atom_j    kumb      min_val|.0      max_val|9.e99
ref_1,1   ...       ref_1,nc
...       ...       ...    
ref_nw,1  ...       ref_nw,nc
------------------------------------------------------------------------
tstp ~ 0.001 / 100.0  (dt / gamma) ~ 1.e-5 (dyn) / 1.e-4 (min)
kumb ~ 3000
------------------------------------------------------------------------
"""
		self.node = node
		f = open( conf, "rt" )
		t = f.readline().strip().split()
		self.ncrd = int( t[0] )
		self.nwin = int( t[1] )
		self.tstp = 0.001 / 100.0
		if( len( t ) == 3 ):
			self.tstp = float( t[2] )
		self.jidx = {}
		self.func = []
		self.atom = []
		self.bcrd = []
		self.rcrd = []
		self.kumb = []
		for i in range( self.ncrd ):
			t = f.readline().strip().split()
			if( t[0][0:4] == "dist" and ( len( t ) == 4 or len( t ) == 6 ) ):
				self.func.append( self.distance )
				a_i = int( t[1] )
				a_j = int( t[2] )
				self.atom.append( ( a_i, a_j ) )
				self.kumb.append( float( t[3] ) )
				if( len( t ) == 6 ):
					self.bcrd.append( ( float( t[4] ), float( t[5] ) ) )
				else:
					self.bcrd.append( ( 0.0, 9.0e99 ) )
				self.jidx[a_i] = True
				self.jidx[a_j] = True
		self.jidx = { jj: ii for ii,jj in zip( range( len( self.jidx ) ), sorted( self.jidx ) ) }
		self.idxj = { self.jidx[ii]: ii for ii in iter( self.jidx ) }
		tmp = []
		for i in range( self.nwin ):
			tmp += [ float( j ) for j in f.readline().strip().split() ]
		f.close()
		self.rcrd = tmp[self.node*self.ncrd:(self.node+1)*self.ncrd]
		if( self.node == 0 ):
			self.icrd = tmp[:]
		self.mass = []
		for i in range( len( self.jidx ) ):
			for j in [0, 1, 2]:
				self.mass.append( 1.0 / molec.mass[self.idxj[i]] )
		self.jcol = 3 * len( self.jidx )
		self.jaco = []
		self.ccrd = []
		self.cmet = []
		# file handlers
		self.fcvs = open( "string.%03d.cvs"%( self.node ), "wt" )
		self.ffrc = open( "string.%03d.frc"%( self.node ), "wt" )
		self.fmet = open( "string.%03d.met"%( self.node ), "wt" )
		if( self.node == 0 ):
			self.fcnv = open( "convergence.dat", "wt" )
			self.fstr = open( "string.dat", "wt" )


	def stop( self ):
		self.fcvs.close()
		self.ffrc.close()
		self.fmet.close()
		if( self.node == 0 ):
			self.fcnv.close()
			self.fstr.close()


	def s_grad( self, molec ):
		# calculate current CVs
		self.ccrd = []
		jaco = [ 0.0 for i in range( self.ncrd * self.jcol ) ]
		for i in range( self.ncrd ):
			self.func[i]( i, molec, jaco )
		# translate gradients into the molecule
		diff = [ self.kumb[i] * ( self.ccrd[i] - self.rcrd[i] ) for i in range( self.ncrd ) ]
		grad = qm3.maths.matrix.mult( diff, 1, self.ncrd, jaco, self.ncrd, self.jcol )
		for i in range( len( self.jidx ) ):
			i3 = i * 3
			for j in [0, 1, 2]:
				molec.grad[3*self.idxj[i]+j] += grad[i3+j]
		# flush current colective variables
		self.fcvs.write( "".join( [ "%20.10lf"%( i ) for i in self.ccrd ] ) + "\n" )
		self.fcvs.flush()
		# flush current forces
		self.ffrc.write( "".join( [ "%20.10lf"%( -i ) for i in diff ] ) + "\n" )
		self.ffrc.flush()
		# calculate current metric tensor M (eq. 7 @ 10.1016/j.cplett.2007.08.017)
		self.cmet = [ 0.0 for i in range( self.ncrd * self.ncrd ) ]
		for i in range( self.ncrd ):
			for j in range( i, self.ncrd ):
				self.cmet[i*self.ncrd+j] = sum( [ jaco[i*self.jcol+k] * self.mass[k] * jaco[j*self.jcol+k] for k in range( self.jcol ) ] )
				self.cmet[j*self.ncrd+i] = self.cmet[i*self.ncrd+j]
		# flush current metric
		self.fmet.write( "".join( [ "%20.10lf"%( i ) for i in self.cmet ] ) + "\n" )
		self.fmet.flush()
		# perform dynamics on the reference CVs and box'em (eq. 17 @ 10.1016/j.cplett.2007.08.017)
		grad = qm3.maths.matrix.mult( diff, 1, self.ncrd, self.cmet, self.ncrd, self.ncrd )
		for i in range( self.ncrd ):
			self.rcrd[i] += grad[i] * self.tstp
			self.rcrd[i] = min( max( self.rcrd[i], self.bcrd[i][0] ), self.bcrd[i][1] )


	def s_dist( self ):
		# redistribute reference CVs by arc-length every step
		qm3.utils._mpi.barrier()
		if( self.node == 0 ):
			# get current string from nodes
			tmp_c = self.rcrd[:]
			tmp_m = self.cmet[:]
			for i in range( 1, self.nwin ):
				tmp_c += qm3.utils._mpi.recv_r8( i, self.ncrd )
				tmp_m += qm3.utils._mpi.recv_r8( i, self.ncrd * self.ncrd )
			# re-parametrize string
			tmp_c = string_distribute( self.ncrd, self.nwin, tmp_c, tmp_m )[0]
			# send back new string to nodes
			for i in range( 1, self.nwin ):
				qm3.utils._mpi.send_r8( i, tmp_c[i*self.ncrd:(i+1)*self.ncrd] )
			self.rcrd = tmp_c[self.node*self.ncrd:(self.node+1)*self.ncrd][:]
			# store re-parametrized string
			self.fstr.write( "".join( [ "%20.10lf"%( tmp_c[j] ) for j in range( self.ncrd * self.nwin ) ] ) + "\n" )
			self.fstr.flush()
			# store current convergence
			ncrd2 = self.ncrd * self.ncrd
			tmp_a = []
			tmp_b = []
			for i in range( self.nwin ):
				tmp_a += [ tmp_c[i*self.ncrd+j] - self.icrd[i*self.ncrd+j] for j in range( self.ncrd ) ]
				tmp_b += qm3.maths.matrix.mult( tmp_i[i*ncrd2:(i+1)*ncrd2], self.ncrd, self.ncrd, tmp_a[i*self.ncrd:(i+1)*self.ncrd], self.ncrd, 1 )
			self.fcnv.write( "%20.10lf\n"%( math.sqrt( sum( [ tmp_a[i] * tmp_b[i] for i in range( self.ncrd * self.nwin ) ] ) / float( self.nwin ) ) ) )
			self.fcnv.flush()
		else:
			qm3.utils._mpi.send_r8( 0, self.rcrd )
			qm3.utils._mpi.send_r8( 0, self.cmet )
			self.rcrd = qm3.utils._mpi.recv_r8( 0, self.ncrd )


	def get_grad( self, molec ):
		self.s_grad( molec )
		self.s_dist()


	def distance( self, icrd, molec, jacob ):
		ai = self.atom[icrd][0]
		aj = self.atom[icrd][1]
		dd = [ (jj-ii) for ii,jj in zip( molec.coor[3*ai:3*ai+3], molec.coor[3*aj:3*aj+3] ) ]
		vv = math.sqrt( sum( [ ii*ii for ii in dd ] ) )
		self.ccrd.append( vv )
		for k in [0, 1, 2]:
			jacob[icrd*self.jcol+3*self.jidx[ai]+k] -= dd[k] / vv
			jacob[icrd*self.jcol+3*self.jidx[aj]+k] += dd[k] / vv



###############################################################################
# Iterable version of the string
#
class stepped_string( string ):

	def __init__( self, node, conf, molec ):
		string.__init__( self, node, conf, molec )


	def get_grad( self, molec ):
		self.s_grad( molec )

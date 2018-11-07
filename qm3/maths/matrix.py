# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math



version = None
try:
	import qm3.maths._matrix
	version = qm3.maths._matrix.version
except:
	try:
		import numpy
		version = "Numpy-" + numpy.version.version
	except:
		pass



def norm( vec ):
	m = math.sqrt( sum( [ i*i for i in vec ] ) )
	return( [ i/m for i in vec ] )



def dot_product( va, vb ):
	return( sum( [ i*j for i,j in zip( va, vb ) ] ) )



def cross_product( *vec ):
	siz = len( vec )
	dim = siz + 1
	o = [ 0.0 for i in range( dim ) ]
	if( sum( [ len( vec[i] ) == dim for i in range( siz ) ] ) == siz ):
		for i in range( dim ):
			t = []
			for j in range( siz ):
				for k in range( dim ):
					if( k != i ):
						t.append( vec[j][k] )
			o[i] = math.pow( -1, i ) * det( t, siz )
			del t
	return( o )



#def ud_ij( mat, n , i, j ):
#	if( i < 0 or j < 0 or i >= n or j >= n ):
#		return( None )
#	elif( j >= i ):
#		return( mat[(i*n)-((i-1)*i)/2+(j-i)] )
#	else:
#		return( mat[(j*n)-((j-1)*j)/2+(i-j)] )
#
#
#def swap_upper_diagonal_to_columns( n, v ):
#	o = range( n*(n+1)/2 )
#	for i in range( n ):
#		o[i*(i+1)/2+i] = v[i*n-(i-1)*i/2]
#		for j in range( i + 1, n ):
#			o[(j-1)*j/2+i+j] = v[i*n-(i-1)*i/2+(j-i)]
#	return( o )
#
#
#def swap_upper_diagonal_to_rows( n, v ):
#	o = range( n*(n+1)/2 )
#	for i in range( n ):
#		o[i*n-(i-1)*i/2] = v[i*(i+1)/2+i]
#		for j in range( i + 1, n ):
#			o[i*n-(i-1)*i/2+(j-i)] = v[(j-1)*j/2+i+j]
#	return( o )



def from_diagonal( lst, row ):
	out = [ 0.0 for i in range( row * row ) ]
	for i in range( row ):
		out[i * row + i] = float( lst[i] )
	return( out )



def from_upper_diagonal_rows( lst, row ):
	if( len( lst ) == row * (row + 1) / 2 ):
		out = [ 0.0 for i in range( row * row ) ]
		k = 0
		for i in range( row ):
			for j in range( i, row ):
				out[i * row + j] = float( lst[k] )
				k += 1
		for i in range( row - 1 ):
			for j in range( i + 1, row ):
				out[j * row + i] = out[i * row + j]
		return( out )
	else:
		raise Exception( "matrix.from_upper_diagonal_rows: invalid dimensions" )



# Import from fortran column-packed symmetric matrix upper diagonal
def from_upper_diagonal_columns( lst, row ):
	if( len( lst ) == row * (row + 1) / 2 ):
		out = [ 0.0 for i in range( row * row ) ]
		k = 0
		for j in range( row ):
			for i in range( j+1 ):
				out[i * row + j] = float( lst[k] )
				k += 1
		for i in range( row - 1 ):
			for j in range( i + 1, row ):
				out[j * row + i] = out[i * row + j]
		return( out )
	else:
		raise Exception( "matrix.from_upper_diagonal_columns: invalid dimensions" )



def get_row( mat, row, col, i ):
	return( mat[i * col : (i + 1) * col] )



def get_column( mat, row, col, j ):
	return( [ mat[i * col + j] for i in range( row ) ] )



def mprint( mat, row, col, fmt = "%10.3lf" ):
	print( "%d rows x %d columns"%( row, col ) )
	for i in range( row ):
		for j in range( col ):
			print( fmt%( mat[i* col + j] ), end = "" )
		print()



def det( mat, row ):
	try:
		return( qm3.maths._matrix.det( mat, row ) )
	except:
		try:
			return( numpy.linalg.det( numpy.array( mat, dtype = float ).reshape( ( row, row ) ) ) )
		except:
			raise Exception( "matrix.det: neither qm3.maths._matrix.so nor numpy available..." )



def T( mat, row, col ):
	out = []
	for j in range( col ):
		for i in range( row ):
			out.append( mat[i * col + j] )
	return( out )



def add_row( mat, row, col, lst ):
	if( len( lst ) == col ):
		mat += [ float( lst[i] ) for i in range( col ) ]



def add_column( mat, row, col, lst ):
	if( len( lst ) == row ):
		for i in range( row-1, -1, -1 ):
			mat.insert( i * col + row - 1, float( lst[i] ) )



def mult( a_mat, a_row, a_col, b_mat, b_row, b_col ):
	if( a_col != b_row ):
		raise Exception( "matrix.mult: wrong dimensions: A_{0..%d,0..%d} * B_{0..%d,0..%d}"%( a_row-1, a_col-1, b_row-1, b_col-1 ) )
	else:
		try:
			return( qm3.maths._matrix.matmul( a_mat, a_row, a_col, b_mat, b_row, b_col ) )
		except:
			try:
				return( numpy.dot( numpy.array( a_mat, dtype = float ).reshape( ( a_row, a_col ) ), numpy.array( b_mat, dtype = float ).reshape( ( b_row, b_col ) ) ).flatten().tolist() )
			except:
				raise Exception( "matrix.mult: neither qm3.maths._matrix.so nor numpy available..." )



def diag( mat, row ):
	if( len( mat ) != row * row ):
		raise Exception( "matrix.diag: non square matrix: %d vs %d^2=%d"%( len( mat ), row, row * row ) )
	else:
		try:
			return( qm3.maths._matrix.diag( mat, row ) )
		except:
			try:
				val, vec = numpy.linalg.eigh( numpy.array( mat, dtype = float ).reshape( ( row, row ) ) )
				return( val.flatten().tolist(), vec.flatten().tolist() )
			except:
				raise Exception( "matrix.diag: neither qm3.maths._matrix.so nor numpy available..." )



def inverse( mat, row, col ):
	if( row == 1 or col == 1 ):
		# V^{{}+{}}_{n,1} = V^T_{n,1} over { V^T_{n,1} * V_{1,n} }
		# V^{{}+{}}_{1,n} = V^T_{1,n} over { V^T_{1,n} * V_{n,1} }
		t = sum( [ mat[i] * mat[i] for i in range( col * row ) ] )
		return( [ mat[i] / t for i in range( col * row ) ] )
	else:
		try:
			return( qm3.maths._matrix.invert( mat, row, col ) )
		except:
			if( row == col ):
				try:
					return( numpy.linalg.inv( numpy.array( mat, dtype = float ).reshape( ( row, col ) ) ).flatten().tolist() )
				except:
					raise Exception( "matrix.inverse: neither qm3.maths._matrix.so nor numpy available..." )
			else:
				try:
					return( numpy.linalg.pinv( numpy.array( mat, dtype = float ).reshape( ( row, col ) ) ).flatten().tolist() )
				except:
					raise Exception( "matrix.inverse: neither qm3.maths._matrix.so nor numpy available..." )



# TODO ============================================================================================
#
#	.to_upper_diagonal_rows()					return upper-diagonal row-packed
#	.to_upper_diagonal_columns()				return upper-diagonal column-packed (fortran)
#
# =================================================================================================

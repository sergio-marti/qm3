import qm3.maths.stats
import qm3.maths.matrix

print( qm3.maths.stats.stats( [ 0.388, 0.389, 0.410 ] ) )
print( qm3.maths.stats.autocorrelation( range( 10 ) ) )
print( qm3.maths.stats.autocorrelation( range( 10 ), 2 ) )
print( qm3.maths.stats.sampling_ratio( range( 10 ) ) )

print( 80*"-" )
print( qm3.maths.stats.k_means( [ 7, 4, 10, 16, 13, 7, 3, 5, 7, 3, 13, 14, 12, 11, 10, 7, 7, 5, 3, 3 ], 2 )[0] )
x = [ [ 7, 4, 10, 16, 13 ], [ 7, 7, 5, 3, 3 ], [ 7, 3, 5, 7, 3 ], [ 13, 14, 12, 11, 10 ] ]
qm3.maths.matrix.mprint( sum( x, [] ), 4, 5 )
o = qm3.maths.stats.PCA( x )
print( o.val )
qm3.maths.matrix.mprint( o.select( [ 0,2,3 ] ), 3, 5 )
qm3.maths.matrix.mprint( o.select( [ 0,2,3 ], False ), 4, 5 )

print( 80*"-" )
import numpy
print( qm3.maths.stats.npk_means( numpy.array( [ 7, 4, 10, 16, 13, 7, 3, 5, 7, 3, 13, 14, 12, 11, 10, 7, 7, 5, 3, 3 ] ), 2 )[0] )
x = numpy.array( [ 7, 4, 10, 16, 13, 7, 7, 5, 3, 3, 7, 3, 5, 7, 3, 13, 14, 12, 11, 10 ] ).reshape( (4,5) )
o = qm3.maths.stats.npPCA( x )
print( o.val )
print( o.select( [ 1,2,3 ] ) )
print( o.select( [ 1,2,3 ], False ) )

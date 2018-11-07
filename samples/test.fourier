import math
import qm3.maths.fourier

n = 100
x = [ -math.pi+2.*math.pi*float(i)/float( n ) for i in range( n+1 ) ]
y = [ math.sin( i ) for i in x ]
y = [ 1.0 for i in range( n // 2 ) ] + [ 0.0 ] + [ -1.0 for i in range( n // 2 ) ]

o = qm3.maths.fourier.series( 100, x, y )
o.integrate()
z = [ o.calc( x[i] ) for i in range( n+1 ) ]
for i in range( n+1 ):
    print( "%20.10lf%20.10lf%20.10lf"%( x[i], y[i], z[i] ) )

#import matplotlib.pyplot as plt
#plt.grid( True )
#plt.ylim( -1.1, 1.1 )
#plt.plot( x, y, "ro", markersize = 8.0 )
#plt.plot( x, z, "g-", linewidth = 3.0 )
#plt.tight_layout()
#plt.savefig( "11.png" )


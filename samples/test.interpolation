import math
import qm3.maths.interpolation

n = 100
x = []
y = []
d = ( 10.0 - 0.0 ) / n
t = 0.0
for i in range( n ):
    x.append( t )
    y.append( math.sin( t ) )
    t += d
gg = qm3.maths.interpolation.gaussian( x, y, 0.005 )
ll = qm3.maths.interpolation.lagrange( x, y, points = 3 )
cs = qm3.maths.interpolation.cubic_spline( x, y )
ss = qm3.maths.interpolation.hermite_spline( x, y, "steffen" )
aa = qm3.maths.interpolation.hermite_spline( x, y, "akima" )
fc = qm3.maths.interpolation.hermite_spline( x, y, "fritsch_carlson" )

r = 3.05
print( "Real:     ", r, math.sin( r ), math.cos( r ) )
print( "Lagrange: ", ll.calc( r ) )
print( "Gaussian: ", gg.calc( r ) )
print( "CSpline:  ", cs.calc( r ) )
print( "Steffen:  ", ss.calc( r ) )
print( "Akima:    ", aa.calc( r ) )
print( "Fritsch.C:", fc.calc( r ) )
print( 80*'-' )

r = ( x[0] + x[1] ) * 0.5
print( "Real:     ", r, math.sin( r ), math.cos( r ) )
print( "Lagrange: ", ll.calc( r ) )
print( "Gaussian: ", gg.calc( r ) )
print( "CSpline:  ", cs.calc( r ) )
print( "Steffen:  ", ss.calc( r ) )
print( "Akima:    ", aa.calc( r ) )
print( "Fritsch.C:", fc.calc( r ) )
print( 80*'-' )

r = ( x[-1] + x[-2] ) * 0.5
print( "Real:     ", r, math.sin( r ), math.cos( r ) )
print( "Lagrange: ", ll.calc( r ) )
print( "Gaussian: ", gg.calc( r ) )
print( "CSpline:  ", cs.calc( r ) )
print( "Steffen:  ", ss.calc( r ) )
print( "Akima:    ", aa.calc( r ) )
print( "Fritsch.C:", fc.calc( r ) )

n = 10
x = [ -1 + i * 2.0/(n-1.0) for i in range( n ) ]
y = x[:]
z = []
for i in range( n ):
    for j in range( n ):
        z.append( x[i] * x[i] + y[j] * y[j] )
i2d = qm3.maths.interpolation.interpolate_2d( x, y, z, 
    interpolant = lambda x,y : qm3.maths.interpolation.hermite_spline( x, y, method = "akima" ) )

print( 80*'-' )
x = [ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. ]
y = [ 10., 10., 10., 10., 10., 10., 10.5, 15., 50., 60., 85. ]
print( "CSpline: ", qm3.maths.interpolation.cubic_spline( x, y ).calc( 9.5 ) )
print( "Steffen: ", qm3.maths.interpolation.hermite_spline( x, y ).calc( 9.5 ) )
print( "Fitpack: ", qm3.maths.interpolation.tensioned_spline( x, y, .0 ).calc( 9.5 ) )
print( "TFitpack:", qm3.maths.interpolation.tensioned_spline( x, y, 10. ).calc( 9.5 ) )

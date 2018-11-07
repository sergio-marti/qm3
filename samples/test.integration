import qm3.maths.interpolation
import qm3.maths.integration

x = []
y = []
a = -10.0
b =  10.0
h = 0.01
t = a
for i in range( int( ( b - a ) / h ) + 2 ):
    x.append( t )
    y.append( t * t )
    t += h
print( "x^2 (-10,10) = 2/3 * 1000 = 666.66666666666663" )
print( "Simpson(discrete):   ", qm3.maths.integration.Simpson(x,y) )
oc = qm3.maths.interpolation.cubic_spline( x, y )
oh = qm3.maths.interpolation.hermite_spline( x, y )
ol = qm3.maths.interpolation.lagrange( x, y )
print( "Simpson(cspln,100):  ", qm3.maths.integration.Simpson_f( lambda x: oc.calc(x)[0], a, b, n = 100 ) )
print( "Simpson(cspln,auto): ", qm3.maths.integration.Simpson_f( lambda x: oc.calc(x)[0], a, b ) )
print( "Simpson(herm.,auto): ", qm3.maths.integration.Simpson_f( lambda x: oh.calc(x)[0], a, b ) )
print( "Simpson(lagra,auto): ", qm3.maths.integration.Simpson_f( lambda x: ol.calc(x)[0], a, b ) )
print( "Simpson(func,auto):  ", qm3.maths.integration.Simpson_f( lambda x: x * x, a, b ) )
print( "Romberg(cspln):      ", qm3.maths.integration.Romberg( lambda x: oc.calc(x)[0], a, b ) )
print( "Romberg(hermite):    ", qm3.maths.integration.Romberg( lambda x: oh.calc(x)[0], a, b ) )
print( "Romberg(lagra):      ", qm3.maths.integration.Romberg( lambda x: ol.calc(x)[0], a, b ) )
print( "Romberg(func):       ", qm3.maths.integration.Romberg( lambda x: x * x, a, b ) )
print( "Gauss-Legendre(func) ", qm3.maths.integration.Gauss( lambda x: x * x, a, b ) )
print( "Trapecios(discrete): ", sum( ( y[i+1] + y[i] ) * ( x[i+1] - x[i] ) * 0.5 for i in range( len( y ) - 1 ) ) )

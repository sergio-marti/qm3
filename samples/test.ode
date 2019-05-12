import math
import qm3.maths.ode

# --- single equation (A --> B, 1st order: depends on f)
print( "\nKinetics: A  --> B\n" )
print( "%12s%12s%12s%12s"%( "Time", "[A](Euler)", "[A](RK4)", "[A](anal)" ) )
a = 0.1
def f( t, x ):
    return( [ - a * x[0] ] )
x0 = [ 10. ]
cx = x0[:]
ex = x0[:]
t0 = 0.0
tf = 10.0
dt = ( tf - t0 ) / 100.0
for i in range( int( ( tf - t0 ) / dt ) ): 
    ct = t0 + i * dt
    if( i%10 == 0 ):
        print( "%12.6lf%12.6lf%12.6lf%12.6lf"%( ct, ex[0], cx[0], math.exp( - a * ct ) * x0[0] ) )
    cx = qm3.maths.ode.RK4_step( f, dt, ct, cx )
    if( i == 0 ):
        e0 = ex[:]
        ex = qm3.maths.ode.Euler_init( f, dt, ct, e0 )
    else:
        tt = qm3.maths.ode.Euler_step( f, dt, ct, ex, e0 )
        e0 = ex[:]
        ex = tt[:]
print( "%12.6lf%12.6lf%12.6lf%12.6lf"%( tf, ex[0], cx[0], math.exp( - a * tf ) * x0[0] ) )


# --- multiple equations (A --> B --> C, 1st order)
print( "\nKinetics: A  --> B --> C\n" )
print( "%12s%12s%12s%12s%12s%12s%12s"%( "Time", "[A](Euler)", "[A](RK4)", "[A](anal)",
    "[B](Euler)", "[B](RK4)", "[B](anal)" ) )
a = 0.1
b = 0.2
def f( t, x ):
    return( [ - a * x[0], a * x[0] - b * x[1] ] )
x0 = [ 10.0, 0.0 ]
cx = x0[:]
ex = x0[:]
t0 = 0.0
tf = 10.0
dt = ( tf - t0 ) / 100.0
for i in range( int( ( tf - t0 ) / dt ) ): 
    ct = t0 + i * dt
    if( i%10 == 0 ):
        print( "%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf"%( ct, ex[0], cx[0],
            x0[0]*math.exp(-a*ct), ex[1], cx[1],
            a/(b-a)*x0[0]*(math.exp(-a*ct)-math.exp(-b*ct)) ) )
    cx = qm3.maths.ode.RK4_step( f, dt, ct, cx )
    if( i == 0 ):
        e0 = ex[:]
        ex = qm3.maths.ode.Euler_init( f, dt, ct, e0 )
    else:
        tt = qm3.maths.ode.Euler_step( f, dt, ct, ex, e0 )
        e0 = ex[:]
        ex = tt[:]
print( "%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf"%( tf, ex[0], cx[0], x0[0]*math.exp(-a*tf), ex[1], cx[1],
    a/(b-a)*x0[0]*(math.exp(-a*tf)-math.exp(-b*tf)) ) )


# -- second order: dynamics
print( "\nDynamics: d^2 q / dt^2 = -1/m dV / dq    on V = x^2+y^2\n" )
print( "%12s%12s%12s%12s%12s%12s%12s%12s"%( "Time", "x", "y", "dx/dt", "dy/dt", "V", "K", "U" ) )
def f( t, x, v ):
    return( [ - 2 * x[0], - 2 * x[1] ] )
cx = [ 10.0, 10.0 ]
cv = [ 0.0, 0.0 ]
t0 = 0.0
tf = 10.0
dt = ( tf - t0 ) / 100.0
for i in range( int( ( tf - t0 ) / dt ) ): 
    ct = t0 + i * dt
    if( i%10 == 0 or True ):
        ek = 0.5 * ( cv[0] * cv[0] + cv[1] * cv[1] )
        ep = cx[0] * cx[0] + cx[1] * cx[1]
        print( "%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf"%( ct, cx[0], cx[1], cv[0], cv[1],
            ep, ek, ep+ek ) )
    cx, cv = qm3.maths.ode.RK42_step( f, dt, ct, cx, cv )
print( "%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf%12.6lf"%( tf, cx[0], cx[1], cv[0], cv[1], ep, ek, ep+ek ) )


# Least-Squares Finite Elements integration (2D)
print( "\n" )
nx  = 50
ny  = 50
x0  = -0.5
y0  = -0.5
dx  = math.fabs( 2.0 * x0 ) / float( nx - 1.0 )
dy  = math.fabs( 2.0 * y0 ) / float( ny - 1.0 )
crx = [ x0 + i * dx for i in range( nx ) ]
cry = [ y0 + i * dy for i in range( ny ) ]
grd = []
for i in range( nx ):
    for j in range( ny ):
        grd.append( [ 2 * crx[i], 2 * cry[j] ] )
F = qm3.maths.ode.least_squares_finite_elements_2d( nx, dx, ny, dy, grd )

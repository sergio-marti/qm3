import qm3.maths.roots

def f( x ):
    return( x*x-2 )

def g( x ):
    return( 2*x )

print( "bisect:        ", qm3.maths.roots.bisect( f, .0, 2. ) )
print( "ridders:       ", qm3.maths.roots.ridders( f, .0, 2. ) )
print( "newton_raphson:", qm3.maths.roots.newton_raphson( f, .0 ) )
print( "halley:        ", qm3.maths.roots.halley( f, g, .1 ) )
print( "newton_raphson:", qm3.maths.roots.multi_newton_raphson( [ lambda x: x[0]*x[0]-2 ], [ .0 ] )[0] )

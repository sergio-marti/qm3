py_exe=python3

rm -f 04.points
$py_exe << EOD >> 04.points
mov = [ (-1,0), (+1,0), (0,-1), (0,+1) ]
lst = []
f = open( "03.path", "rt" )
for l in f:
    t = l.split()
    if( len( t ) == 3 ):
        ci = int( round( ( float( t[0] ) - 1.4 ) / 0.1, 0 ) )
        cj = int( round( ( float( t[1] ) - 1.4 ) / 0.1, 0 ) )
        k  = "%d:%d"%( ci, cj )
        if( not k in lst ):
            lst.append( k )
        for i,j in mov:
            ti = min( 22, max( 0, ci + i ) )
            tj = min( 22, max( 0, cj + j ) )
            k  = "%d:%d"%( ti, tj )
            if( not k in lst ):
                lst.append( k )
x = []
y = []
for k in lst:
    t = k.split( ":" )    
    x.append( 1.4 + float( t[0] ) * 0.1 )
    y.append( 1.4 + float( t[1] ) * 0.1 )
    print( " ".join( t ) )

import qm3.maths.grids
import matplotlib.pyplot as pyplot
g = qm3.maths.grids.grid( "03.grid" )
pyplot.clf()
pyplot.grid( True )
pyplot.plot( x, y, 'og' )
g.plot2d()
EOD

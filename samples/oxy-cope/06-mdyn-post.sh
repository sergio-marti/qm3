py_exe=python3

$py_exe << EOD
import os
import matplotlib.pyplot as pyplot

d1 = []
d2 = []
tt = []
f = open( "colvar", "rt" )
for l in f:
    t = l.split()
    if( len( t ) == 5 ):
        tt.append( float( t[0] ) )
        d1.append( float( t[1] ) * 10 )
        d2.append( float( t[2] ) * 10 )
f.close()
pyplot.clf()
pyplot.grid( True )
pyplot.plot( tt, d1, '-g', linewidth = 2 )
pyplot.plot( tt, d2, '-b', linewidth = 2 )
pyplot.savefig( "06.mdyn_coor.pdf" )
pyplot.show()

pyplot.clf()
pyplot.grid( True )
pyplot.xlim( ( -3.0, 3.0 ) )
pyplot.ylim( ( -35.0, 5.0 ) )
f = open( "hills", "rt" )
dat = f.readlines()
f.close()
l = []
for r,c in [ (2,'r'), (4,'b'), (6,'g'), (8,'y'), (10,'k') ]:
    f = open( "tmp", "wt" )
    f.write( "".join( dat[0:r*1000+1] ) )
    f.close()
    os.system( "./bin/plumed sum_hills --hills tmp" )
    x = []
    y = []
    m = -9e99
    f = open( "fes.dat", "rt" )
    b = f.readlines()
    f.close()
    os.unlink( "fes.dat" )
    for i in range( 5, len( b ) ):
        t = b[i].split()
        x.append( float( t[0] ) * 10 )
        y.append( float( t[1] ) / 4.184 )
        if( x[-1] > -0.1 and x[-1] < 0.1 ):
            m = max( m, y[-1] )
    l.append( pyplot.plot( x, [ i - m for i in y ], '-'+c, linewidth = 2 ) )
    l[-1][0].set_label( "%d ps"%( r * 10 ) )
pyplot.legend()
pyplot.savefig( "06.mdyn_fes.pdf" )
pyplot.show()
EOD

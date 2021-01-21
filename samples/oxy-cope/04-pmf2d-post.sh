py_exe=python3

$py_exe << EOD
import glob
import qm3.utils.free_energy
import qm3.maths.grids

pmf = qm3.utils.free_energy.umbint_2d( glob.glob( "04.dat.??.??" ) )
pmf.setup( ( 22 * 2, 22 * 2 ) )
pmf.integrate()
grd = qm3.maths.grids.grid( "umbint_2d" )
grd.z = [ min( i, 220 ) / 4.184 for i in grd.z ]
grd.save( "04.pmf" )
grd.plot2d( levels = 40, fname = "04.pmf2d.pdf" )
grd.plot3d()
EOD


$py_exe << EOD
import glob
import matplotlib.pyplot as pyplot
pyplot.clf()
pyplot.grid( True )
for fn in glob.glob( "04.dat.??.??" ):
    print( fn )
    x = []
    y = []
    f = open( fn, "rt" )
    f.readline()
    for l in f:
        t = l.split()
        x.append( float( t[0] ) )
        y.append( float( t[1] ) )
    f.close()
    pyplot.plot( x, y, '-' )
pyplot.savefig( "04.pmf2d.overlap.pdf" )
pyplot.show()
EOD


rm -f 04.summary
../../tools/s_loca.py 04.pmf 1.54 2.99 -1 >> 04.summary
echo >> 04.summary
echo "===================================" >> 04.summary
echo >> 04.summary
../../tools/s_loca.py 04.pmf 1.66 1.89  0 >> 04.summary
echo >> 04.summary
echo "===================================" >> 04.summary
echo >> 04.summary
../../tools/s_loca.py 04.pmf 3.16 1.55 -1 >> 04.summary
echo >> 04.summary
echo "===================================" >> 04.summary


$py_exe << EOD >> 04.summary
import  re
f = open( "04.summary", "rt" )
p = re.compile( "\[([0-9\.\-]+)\, ([0-9\.\-]+)\] ([0-9\.\-]+)" ).findall( f.read() )
f.close()
for k in range( 1, len( p ), 2 ):
    ri = float( p[k][0] )
    rj = float( p[k][1] )
    ee = round( float( p[k][2] ), 2 )
    ci = int( round( ( ri - 1.4 ) / 0.1, 0 ) )
    cj = int( round( ( rj - 1.4 ) / 0.1, 0 ) )
    print( round( ri, 4 ), round( rj, 4 ), round( ee, 2 ), ci, cj )
EOD



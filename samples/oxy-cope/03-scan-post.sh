py_exe=python3

$py_exe << EOD
import  qm3.maths.grids
import  glob
import  io
f = io.StringIO()
for grd in glob.glob( "03.ene.*" ):
    g = open( grd, "rt" )
    for l in g:
        f.write( l )
    g.close()
f.seek( 0 )
grd = qm3.maths.grids.grid()
grd.regular( f, ( 31, 31 ), (.1, .1) )
t = min( grd.z )
grd.z = [ ( i - t ) / 4.184 for i in grd.z ]
grd.save( "03.grid" )
grd.plot2d( levels = 40, fname = "03.scan.pdf" )
grd.plot3d()
EOD

rm -f 03.summary
../../tools/s_loca.py 03.grid 1.54 2.82 -1 >> 03.summary
echo >> 03.summary
echo "===================================" >> 03.summary
echo >> 03.summary
../../tools/s_loca.py 03.grid 1.66 1.93 0  >> 03.summary
echo >> 03.summary
echo "===================================" >> 03.summary
echo >> 03.summary
../../tools/s_loca.py 03.grid 2.87 1.54 -1 >> 03.summary

$py_exe << EOD >> 03.summary
import  re
f = open( "03.summary", "rt" )
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

../../tools/s_path.py 03.grid 1.6792 1.8627 >& /dev/null
mv path.log 03.path

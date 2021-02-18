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


$py_exe << EOD
import  os
import  qm3.utils
import  qm3.mol

import  matplotlib.pyplot as plt
from    matplotlib.backends.backend_pdf import PdfPages
import  qm3.maths.grids


m = qm3.mol.molecule( "03.pes.00.00.pdb" )
bnd = [ ( m.indx["A"][1]["C2"], m.indx["A"][1]["C3"] ),
        ( m.indx["A"][1]["C5"], m.indx["A"][1]["C6"] ),
        ( m.indx["A"][1]["C3"], m.indx["A"][1]["O7"] ) ]


if( not os.path.isfile( "03.info" ) ):
    f = open( "03.info", "wt" )
    f.write( "%-10s"%( "#" ) )
    for ii,jj in bnd:
        f.write( "%12s"%( m.labl[ii] + "-" + m.labl[jj] ) )
    f.write( "\n" )
    for i in range( 0, 22 ):
        for j in range( 0, 22 ):
            if( os.path.isfile( "03.pes.%02d.%02d.pdb"%( i, j ) ) ):
                m.pdb_read( "03.pes.%02d.%02d.pdb"%( i, j ) )
                print( i, j )
                b = []
                for ii,jj in bnd:
                    b.append( "%12.6lf"%( qm3.utils.distance( m.coor[3*ii:3*ii+3], m.coor[3*jj:3*jj+3] ) ) )
                f.write( "%4d%4d  %s\n"%( i, j, "".join( b ) ) )
                f.flush()
        f.write( "\n" )
    f.close()


pdf = PdfPages( "03.info.pdf" )
grd = qm3.maths.grids.grid()
for k in range( len( bnd ) ):
    grd.parse( "03.info", "0:1:%d"%( 2 + k ) )
    plt.title( m.labl[bnd[k][0]] + "-" + m.labl[bnd[k][1]] )
    plt.grid( True )
    rz  = grd.rotate()
    nx  = len( grd.x )
    ny  = len( grd.y )
    lx  = []
    ly  = []
    lz  = []
    for i in range( ny ):
        lx.append( grd.x[:] )
        ly.append( nx * [ grd.y[i] ] )
        lz.append( rz[i*nx:(i+1)*nx][:] )
    cnt = plt.contourf( lx, ly, lz, cmap = "coolwarm" )
    plt.colorbar( cnt )
    pdf.savefig()
    plt.close()
pdf.close()
EOD

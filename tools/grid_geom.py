#!/usr/bin/env python3
import  os
import  qm3.utils
import  qm3.engines.dynamo

import  matplotlib.pyplot as plt
from    matplotlib.backends.backend_pdf import PdfPages
import  qm3.maths.grids


m = qm3.engines.dynamo.coordinates_read( "crd.0.0" )
bnd = [ ( m.indx["A"][145]["SG"], m.indx["A"][613]["C36"] ),
        ( m.indx["A"][145]["SG"], m.indx["A"][145]["HG"] ),
        ( m.indx["A"][145]["HG"], m.indx["A"][41]["NE2"] ),
        ( m.indx["A"][145]["HG"], m.indx["A"][613]["C72"] ) ]


if( not os.path.isfile( "info" ) ):
    f = open( "info", "wt" )
    f.write( "%-10s"%( "#" ) )
    for ii,jj in bnd:
        f.write( "%12s"%( m.labl[ii] + "-" + m.labl[jj] ) )
    f.write( "\n" )
    for i in range( 0, 51 ):
        for j in range( 0, 51 ):
            if( os.path.isfile( "crd.%d.%d"%( i, j ) ) ):
                m = qm3.engines.dynamo.coordinates_read( "crd.%d.%d"%( i, j ) )
                print( i, j )
                b = []
                for ii,jj in bnd:
                    b.append( "%12.6lf"%( qm3.utils.distance( m.coor[3*ii:3*ii+3], m.coor[3*jj:3*jj+3] ) ) )
                f.write( "%4d%4d  %s\n"%( i, j, "".join( b ) ) )
                f.flush()
        f.write( "\n" )
    f.close()


pdf = PdfPages( "info.pdf" )
grd = qm3.maths.grids.grid()
for k in range( len( bnd ) ):
    grd.parse( "info", "0:1:%d"%( 2 + k ) )
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
#    plt.imshow( lz, cmap = "coolwarm" )
#    plt.colorbar()
    pdf.savefig()
    plt.close()
pdf.close()

#!/usr/bin/env python3

import  qm3.mol
import  qm3.engines.colvar_s
import  qm3.utils
import  qm3.elements
import  qm3.maths.matrix
import  qm3.maths.interpolation
import  matplotlib.pyplot as plt
import  matplotlib.backends.backend_pdf
import  pickle
import  math
import  os


pdf = matplotlib.backends.backend_pdf.PdfPages( "cvs_dist.pdf" )

bnd = [ (2240,9403), (2240,2241), (616,2241), (9407,9404) ]
win = 60
dim = len( bnd )
dat = [ [] for i in range( dim ) ]


g = open( "pmf_s.cnf", "wt" )
g.write( "%d %d\n"%( dim, win ) )
for a1,a2 in bnd:
    g.write( "dist%10d%10d\n"%( a1, a2 ) )
g.close()


for i in range( win ):
    print( i )
    m = qm3.mol.molecule( "../node.%02d"%( i ) )
    for j in range( dim ):
        dat[j].append( qm3.utils.distance( m.coor[3*bnd[j][0]:3*bnd[j][0]+3], m.coor[3*bnd[j][1]:3*bnd[j][1]+3] ) )

g = open( "string", "wt" )
for i in range( win ):
    for j in range( dim ):
        g.write( "%12.4lf"%( dat[j][i] ) )
    g.write( "\n" )
g.close()



# ----------------------------------------------------------------------------------------------------
class metric_tensor( object ):
    def __init__( self, conf ):
        f = open( conf, "rt" )
        t = f.readline().strip().split()
        self.ncrd = int( t[0] )
        self.nwin = int( t[1] )
        self.jidx = {}
        self.func = []
        self.atom = []
        for i in range( self.ncrd ):
            t = f.readline().strip().split()
            if( t[0][0:4] == "dist" and len( t ) == 3 ):
                self.func.append( self.distance )
                a_i = int( t[1] )
                a_j = int( t[2] )
                self.atom.append( ( a_i, a_j ) )
                self.jidx[a_i] = True
                self.jidx[a_j] = True
        f.close()
        self.jidx = { jj: ii for ii,jj in zip( range( len( self.jidx ) ), sorted( self.jidx ) ) }
        self.idxj = { self.jidx[ii]: ii for ii in iter( self.jidx ) }
        self.jcol = 3 * len( self.jidx )


    def build( self, molec ):
        jaco = [ 0.0 for i in range( self.ncrd * self.jcol ) ]
        for i in range( self.ncrd ):
            self.func[i]( i, molec, jaco )
        cmet = [ 0.0 for i in range( self.ncrd * self.ncrd ) ]
        for i in range( self.ncrd ):
            for j in range( i, self.ncrd ):
                cmet[i*self.ncrd+j] = sum( [ jaco[i*self.jcol+k] * jaco[j*self.jcol+k] for k in range( self.jcol ) ] )
                cmet[j*self.ncrd+i] = cmet[i*self.ncrd+j]
        return( cmet )


    def distance( self, icrd, molec, jacob ):
        ai = self.atom[icrd][0]
        aj = self.atom[icrd][1]
        dd = [ (jj-ii) for ii,jj in zip( molec.coor[3*ai:3*ai+3], molec.coor[3*aj:3*aj+3] ) ]
        vv = math.sqrt( sum( [ ii*ii for ii in dd ] ) )
        for k in [0, 1, 2]:
            jacob[icrd*self.jcol+3*self.jidx[ai]+k] -= dd[k] / vv
            jacob[icrd*self.jcol+3*self.jidx[aj]+k] += dd[k] / vv
# ----------------------------------------------------------------------------------------------------


if( True ):
    print( "* calculating metrics..." )
    x = metric_tensor( "pmf_s.cnf" )
    f = open( "pmf_s.met", "wt" )
    for i in range( win ):
        c = x.build( qm3.mol.molecule( "../../node.%02d"%( i ) ) )
        f.write( "".join( [ "%20.10lf"%( i ) for i in c ] ) + "\n" )
    f.close()
    qm3.maths.matrix.mprint( c, dim, dim )
else:
    print( "* identity matrix for metrics..." )
    m = qm3.mol.molecule( "node.30" )
    c = [ 0.0 for i in range( dim * dim ) ]
    for i in range( dim ):
        c[i*dim+i] = 1.0
    f = open( "pmf_s.met", "wt" )
    for i in range( win ):
        f.write( "".join( [ "%20.10lf"%( i ) for i in c ] ) + "\n" )
    f.close()
    qm3.maths.matrix.mprint( c, dim, dim )


eng = qm3.engines.colvar_s.colvar_s( 0.0, 0.0, "pmf_s.cnf", "../pmf_s.str", "pmf_s.met" )
plt.clf()
plt.grid( True )
plt.plot( eng.arcl, '-o' )
pdf.savefig()
plt.show()

plt.clf()
plt.grid( True )
arc = [ 0.0 ]
for i in range( win - 1 ):
    arc.append( arc[-1] + eng.arcl[i] )
tad = [ [] for i in range( dim ) ]
x = [ arc[-1] / ( win - 1 ) * i for i in range( win ) ]
for j in range( dim ):
    fix = qm3.maths.interpolation.hermite_spline( arc, dat[j] )
    for i in range( win ):
        tad[j].append( fix.calc( x[i] )[0] )
    plt.plot( dat[j], '-' )
    plt.plot( tad[j], 'o' )
pdf.savefig()
plt.show()

g = open( "pmf_s.str", "wt" )
for i in range( win ):
    for j in range( dim ):
        g.write( "%12.4lf"%( tad[j][i] ) )
    g.write( "\n" )
g.close()

eng = qm3.engines.colvar_s.colvar_s( 0.0, 0.0, "pmf_s.cnf", "pmf_s.str", "pmf_s.met" )
plt.clf()
plt.grid( True )
plt.plot( eng.arcl, '-o' )
pdf.savefig()
plt.show()

pdf.close()

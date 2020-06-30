# ----------------------------------------------------------------------------
#morse(d,a,x)=d*(1-exp(-a*(x-1.5)))**2
#cross(d,a,b,c,x,y)=d*exp(-(a*(x-1.5)**2+b*(x-1.5)*(y-1.5)+c*(y-1.5)**2))
#unset key
#unset colorbox
#set grid
#set contour
#set cntrparam levels auto 25
#set range[1:6]
#set yrange[1:6]
#set surface
#set isosamples 40
#set pm3d at s
#set hidden3d trianglepattern 7
#splot morse(10,1,x) + morse(10,1,y) + cross(20,1,1,1,x,y)
# ----------------------------------------------------------------------------
# mpirun -n 4 python test.pneb
# ----------------------------------------------------------------------------
import qm3.problem
import qm3.engines.mmres
import math
import qm3.actions.neb
import qm3.actions.minimize
import qm3.maths.matrix
import qm3.utils._mpi
import sys


class my_pes( qm3.problem.template ):
    def __init__( self, coor ):
        self.size = 9
        self.func = 0.0
        self.coor = coor[:]
        self.grad = []
        self.hess = []
        self.mass = [ 1.0, 1.0, 1.0 ]
        self.fix  = qm3.engines.mmres.angle( 2., 179., [ 0, 1, 2 ] )


    def neb_data( self, node ):
        f = open( "../node.%02d"%( node ), "wt" )
        f.write( "3\n" )
        d1, d2 = self.distances()
        f.write( "@%12.6lf%12.6lf%12.6lf\n"%( d1, d2, self.func ) )
        f.write( "X%20.10lf%20.10lf%20.10lf\n"%( self.coor[0], self.coor[1], self.coor[2] ) )
        f.write( "X%20.10lf%20.10lf%20.10lf\n"%( self.coor[3], self.coor[4], self.coor[5] ) )
        f.write( "X%20.10lf%20.10lf%20.10lf\n"%( self.coor[6], self.coor[7], self.coor[8] ) )
        f.close()


    def distances( self ):
        return(
            math.sqrt( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip( self.coor[0:3], self.coor[3:6] ) ] ) ),
            math.sqrt( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip( self.coor[3:6], self.coor[6:9] ) ] ) ) )


    def get_func( self ):
        r1, r2 = self.distances()
        self.func  = 10.0 * math.pow( 1.0 - math.exp( - ( r1 - 1.5 ) ), 2.0 )
        self.func += 10.0 * math.pow( 1.0 - math.exp( - ( r2 - 1.5 ) ), 2.0 )
        self.func += 20.0 * math.exp( - ( math.pow( r1 - 1.5, 2.0 ) + ( r1 - 1.5 ) * ( r2 - 1.5 )
            + math.pow( r2 - 1.5, 2.0 ) ) )
        self.fix.get_func( self )


    def get_grad( self ):
        self.grad = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
        d1 = [ (ii-jj) for ii,jj in zip( self.coor[0:3], self.coor[3:6] ) ]
        r1 = math.sqrt( sum( [ ii*ii for ii in d1 ] ) )
        e1 = math.exp( - ( r1 - 1.5 ) )
        self.func = 10.0 * math.pow( 1.0 - e1, 2.0 )
        for i in [ 0, 1, 2 ]:
            self.grad[i]   += 20.0 * e1 * ( 1.0 - e1 ) * d1[i] / r1
            self.grad[3+i] -= 20.0 * e1 * ( 1.0 - e1 ) * d1[i] / r1
        d2 = [ (ii-jj) for ii,jj in zip( self.coor[3:6], self.coor[6:9] ) ]
        r2 = math.sqrt( sum( [ ii*ii for ii in d2 ] ) )
        e2 = math.exp( - ( r2 - 1.5 ) )
        self.func += 10.0 * math.pow( 1.0 - e2, 2.0 )
        for i in [ 0, 1, 2 ]:
            self.grad[3+i] += 20.0 * e2 * ( 1.0 - e2 ) * d2[i] / r2
            self.grad[6+i] -= 20.0 * e2 * ( 1.0 - e2 ) * d2[i] / r2
        c1 = ( r1 - 1.5 )
        c2 = ( r2 - 1.5 )
        cx = ( r1 - 1.5 ) * ( r2 - 1.5 )
        cf = 20.0 * math.exp( - ( c1 * c1 + cx + c2 * c2 ) )
        self.func += cf
        for i in [ 0, 1, 2 ]:
            self.grad[i]   -= cf * ( 2.0 * c1 + c2 ) * d1[i] / r1
            self.grad[3+i] += cf * ( 2.0 * c1 + c2 ) * d1[i] / r1
            self.grad[3+i] -= cf * ( c1 + 2.0 * c2 ) * d2[i] / r2
            self.grad[6+i] += cf * ( c1 + 2.0 * c2 ) * d2[i] / r2
        self.fix.get_grad( self )



reac = [ -1.286, -1.264, -1.235, -0.375, -0.363, -0.296, 1.482, 1.521, 1.586 ]
prod = [ -1.559, -1.465, -1.552,  0.331,  0.344,  0.372, 1.270, 1.236, 1.291 ]

nwin = 24

node, ncpu = qm3.utils._mpi.init()

if( node == 0 ):
    def std_log( txt ):
        sys.stdout.write( txt + "\n" )
        sys.stdout.flush()
else:
    def std_log( txt ):
        pass

os.mkdir( "%02d"%( node ) )
os.chdir( "%02d"%( node ) )

oneb = qm3.actions.neb.parall_neb( my_pes( reac ), [0, 1, 2], 1., node, ncpu,
    qm3.actions.neb.distribute( nwin, [ reac, prod ] ) )
qm3.actions.minimize.fire( oneb, step_number = 200, gradient_tolerance = 0.001, step_size = 0.1,
    print_frequency = 1, log_function = std_log )
qm3.utils._mpi.stop()

# for i in {0..39}; do python test.colvar_s $i; done; python test.colvar_s.integrate

import sys
import math
import qm3.problem
import qm3.engines.mmres
import qm3.engines.colvar_s
import qm3.actions.dynamics
try:
    import cStringIO as io
except:
    import io


class my_pes( qm3.problem.template ):
    def __init__( self, node, coor ):
        self.size = 9
        self.func = 0.0
        self.coor = coor[:]
        self.grad = []
        self.hess = []
        self.mass = [ 1.0, 1.0, 1.0 ]
        self.fix  = qm3.engines.mmres.angle( 2., 179., [ 0, 1, 2 ] )
        self.node = node
        f = io.StringIO( """2 24
dist      0     1
dist      1     2
""" )
        f.seek( 0 ) 
        kumb      = 100.
        # colvar_s: 0 ... 1.4  >> 40 windows / i: 0 ... 39 / 0.0 + i * 0.035
        xref      = 0.0 + node * 0.035
        self.cvs  = qm3.engines.colvar_s.colvar_s( kumb, xref, f, "string.colvars", "string.metrics", self )
        f.close()
        self.dat  = open( "pmf.%03d.dat"%( self.node ), "wt" )
        self.dat.write( "%20.10lf%20.10lf\n"%( kumb, xref ) )
        self.trj  = open( "pmf.%03d.xyz"%( self.node ), "wt" )
        self.flg  = open( "pmf.%03d.log"%( self.node ), "wt" )
        self.isok = False
        self.curr = None


    def current_step( self, istep ):
        if( self.isok ):
            if( istep % 100 == 0 ):
                self.trj.write( "3\n" )
                self.trj.write( "%12.6lf%12.6lf\n"%( self.distances() ) )
                self.trj.write( "X%20.10lf%20.10lf%20.10lf\n"%( self.coor[0], self.coor[1], self.coor[2] ) )
                self.trj.write( "X%20.10lf%20.10lf%20.10lf\n"%( self.coor[3], self.coor[4], self.coor[5] ) )
                self.trj.write( "X%20.10lf%20.10lf%20.10lf\n"%( self.coor[6], self.coor[7], self.coor[8] ) )
                self.trj.flush()
            self.dat.write( "%20.10lf\n"%( self.curr ) )


    def log( self, txt ):
        self.flg.write( txt + "\n" )
        self.flg.flush()


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
        self.curr = self.cvs.get_grad( self )



tran = [ -1.4294229501, -1.3751538183, -1.3866947262,
         -0.0140125282,  0.0355803387,  0.0421168135,
          1.3610977165,  1.3502627200,  1.4284399852 ]

obj = my_pes( int( sys.argv[1] ), tran )
qm3.actions.dynamics.assign_velocities( obj, temperature = 10., project = True )
obj.isok = False
qm3.actions.dynamics.langevin_verlet( obj, step_size = 0.001,
    temperature = 10., gamma_factor = 50., print_frequency = 100, project = True,
    log_function = obj.log, step_number = 1000 )
obj.isok = True
qm3.actions.dynamics.langevin_verlet( obj, step_size = 0.001,
    temperature = 10., gamma_factor = 50., print_frequency = 100, project = True,
    log_function = obj.log, step_number = 5000 )
obj.dat.close()
obj.flg.close()
obj.trj.close()

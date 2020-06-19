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
# mpirun -n 4 python test.string
# ----------------------------------------------------------------------------
import qm3.problem
import qm3.engines.mmres
import math
import qm3.actions.dynamics
import qm3.actions.string
import qm3.maths.matrix
import qm3.utils._mpi
import struct
import os


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
        self.cvs  = qm3.actions.string.string( self.node, "string.def", self )
        self.trj  = open( "string.%03d.xyz"%( self.node ), "wt" )
        self.flg  = open( "string.%03d.log"%( self.node ), "wt" )
        self.log( "Starting point: %12.6lf%12.6lf\n"%( self.distances() ) )


    def current_step( self, istep ):
        if( istep % 100 == 0 ):
            self.trj.write( "3\n" )
            self.trj.write( "%12.6lf%12.6lf\n"%( self.distances() ) )
            self.trj.write( "X%20.10lf%20.10lf%20.10lf\n"%( self.coor[0], self.coor[1], self.coor[2] ) )
            self.trj.write( "X%20.10lf%20.10lf%20.10lf\n"%( self.coor[3], self.coor[4], self.coor[5] ) )
            self.trj.write( "X%20.10lf%20.10lf%20.10lf\n"%( self.coor[6], self.coor[7], self.coor[8] ) )
            self.trj.flush()


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
        self.cvs.get_grad( self )



ncrd = 2
nwin = 24

reac = [ -1.286, -1.264, -1.235, -0.375, -0.363, -0.296, 1.482, 1.521, 1.586 ]
prod = [ -1.559, -1.465, -1.552,  0.331,  0.344,  0.372, 1.270, 1.236, 1.291 ]
diff = [ (i-j)/(nwin-1) for i,j in zip( prod, reac ) ]

node, ncpu = qm3.utils._mpi.init()
size = nwin // ncpu

if( node == 0 ):
    f = open( "string.def", "wt" )
    f.write( "%d %d 1.e-5\n"%( ncrd, nwin ) )
    f.write( "dist %6d%6d%8.1lf%6.1lf%6.1lf\n"%( 0, 1, 300., 1.0, 6.0 ) )
    f.write( "dist %6d%6d%8.1lf%6.1lf%6.1lf\n"%( 1, 2, 300., 1.0, 6.0 ) )
    for i in range( nwin ):
        coor = [ ii+i*jj for ii,jj in zip( reac, diff ) ]
        r1 = math.sqrt( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip( coor[0:3], coor[3:6] ) ] ) )
        r2 = math.sqrt( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip( coor[3:6], coor[6:9] ) ] ) )
        f.write( "%12.6lf%12.6lf\n"%( r1, r2 ) )
    f.close()
qm3.utils._mpi.barrier()

cwd = os.getcwd()
obj_list = []
for j in range( size ):
#    os.mkdir( cwd + os.sep + str( size*node+j ) )
#    os.chdir( cwd + os.sep + str( size*node+j ) )
    coor = [ ii+(size*node+j)*jj for ii,jj in zip( reac, diff ) ]
    obj_list.append( my_pes( size*node+j, coor ) )

for j in range( size ):
#    os.chdir( cwd + os.sep + str( size*node+j ) )
    qm3.actions.dynamics.assign_velocities( obj_list[j], temperature = 10., project = True )

dyn_list = []
for j in range( size ):
#    os.chdir( cwd + os.sep + str( size*node+j ) )
    dyn_list.append( qm3.actions.dynamics.langevin_verlet( obj_list[j], step_size = 0.001,
        temperature = 10., gamma_factor = 50., print_frequency = 100, project = True,
        log_function = obj_list[j].log ) )

for i in range( 20000 ):
    for j in range( size ):
#        os.chdir( cwd + os.sep + str( size*node+j ) )
        dyn_list[j].integrate()
    qm3.utils._mpi.barrier()
    if( node == 0 ):
#        os.chdir( cwd + os.sep + "0" )
        tmp_c = []
        tmp_m = []
        for j in range( size ):
            tmp_c += obj_list[j].cvs.rcrd[:]
            tmp_m += obj_list[j].cvs.cmet[:]
        for k in range( 1, ncpu ):
            for j in range( size ):
                tmp_c += qm3.utils._mpi.recv_r8( k, ncrd )
                tmp_m += qm3.utils._mpi.recv_r8( k, ncrd * ncrd )
        tmp_c = qm3.actions.string.string_distribute( ncrd, nwin, tmp_c, tmp_m )[0]
        for j in range( size ):
            obj_list[j].cvs.rcrd = tmp_c[j*ncrd:(j+1)*ncrd][:]
        for k in range( 1, ncpu ):
            for j in range( size ):
                qm3.utils._mpi.send_r8( k, tmp_c[(size*k+j)*ncrd:(size*k+j+1)*ncrd] )
        obj_list[0].cvs.fstr.write( "".join( [ "%20.10lf"%( tmp_c[j] ) for j in range( ncrd * nwin ) ] ) + "\n" )
        obj_list[0].cvs.fstr.flush()
        ncrd2 = ncrd * ncrd
        tmp_a = []
        tmp_b = []
        for i in range( nwin ):
            tmp_i = qm3.maths.matrix.inverse( [ tmp_m[i*ncrd2+j] for j in range( ncrd2 ) ], ncrd, ncrd )
            tmp_a += [ tmp_c[i*ncrd+j] - obj_list[0].cvs.icrd[i*ncrd+j] for j in range( ncrd ) ]
            tmp_b += qm3.maths.matrix.mult( tmp_i, ncrd, ncrd, tmp_a[i*ncrd:(i+1)*ncrd], ncrd, 1 )
        obj_list[0].cvs.fcnv.write( "%20.10lf\n"%( math.sqrt( sum( [ tmp_a[i] * tmp_b[i] for i in range( ncrd * nwin ) ] )
             / float( nwin ) ) ) )
        obj_list[0].cvs.fcnv.flush()
    else:
        for j in range( size ):
            qm3.utils._mpi.send_r8( 0, obj_list[j].cvs.rcrd )
            qm3.utils._mpi.send_r8( 0, obj_list[j].cvs.cmet )
        for j in range( size ):
            obj_list[j].cvs.rcrd = qm3.utils._mpi.recv_r8( 0, ncrd )

for j in range( size ):
#    os.chdir( cwd + os.sep + str( size*node+j ) )
    dyn_list[j].stats()
    obj_list[j].cvs.stop()
    obj_list[j].flg.close()
    obj_list[j].trj.close()

    f = open( "string.%03d.rst"%( size*node+j ), "wb" )
    for k in range( obj_list[j].size ):
        f.write( struct.pack( "d", obj_list[j].coor[k] ) )
    for k in range( obj_list[j].size ):
        f.write( struct.pack( "d", obj_list[j].velo[k] ) )
    f.close()

qm3.utils._mpi.stop()

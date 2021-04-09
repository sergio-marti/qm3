#
# mpirun -n 8 python 04-string.py
#
import  qm3.mol
import  qm3.problem
import  qm3.elements
import  qm3.engines.mol_mech
import  qm3.engines.xtb
import  qm3.engines.mmres
import  qm3.actions.string
import  qm3.actions.dynamics
import  qm3.utils._mpi
import  math
import  os
import  time
import  pickle


# ---------------------------------------------------
# >> problem config

class my_problem( qm3.problem.template ):
    def __init__( self, node ):
        qm3.problem.template.__init__( self )

        self.mol = qm3.mol.molecule( "../04.str_seed.%02d.pdb"%( node ) )
        self.mol.boxl = [ 40., 40., 40. ]

        self.emm = qm3.engines.mol_mech.simple_force_field( self.mol, 8 )
        self.emm.cut_on   = 12.0
        self.emm.cut_off  = 14.0
        self.emm.cut_list = 18.0
        self.emm.system_read( "../02.mm_data.pk", self.mol )
        self.mol.mass = [ qm3.elements.mass[i] for i in self.mol.anum ]

        f = open( "../01.sele_QM.pk", "rb" )
        sqm = pickle.load( f )
        f.close()
        f = open( "../02.sele_MM.pk", "rb" )
        smm = pickle.load( f )
        f.close()
        f = open( "../02.sele.pk", "rb" )
        self.sele = pickle.load( f )
        f.close()
        f = open( "../02.fixed.pk", "rb" )
        for i in pickle.load( f ):
            self.emm.free[i] = False
        f.close()

        self.emm.qm_atoms( sqm )

        self.eqm = qm3.engines.xtb.dl_xtb( self.mol, 0, 0, sqm, smm )

        self.cvs = qm3.actions.string.string( node, "../04.string.def", 0.0 )

        self.size = len( self.sele ) * 3
        self.coor = []
        self.mass = []
        for i in self.sele:
            i3 = i * 3
            self.coor += self.mol.coor[i3:i3+3]
            self.mass.append( self.mol.mass[i] )

        self.log = open( "log", "wt" )


    def logging( self, txt ):
        self.log.write( txt + "\n" )
        self.log.flush()


    def update_coor( self ):
        for i in range( len( self.sele ) ):
            i3 = i * 3
            I3 = self.sele[i] * 3
            for j in [0, 1, 2]:
                self.mol.coor[I3+j] = self.coor[i3+j]


    def get_func( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.emm.get_func( self.mol )
        self.eqm.get_func( self.mol )
        self.dat = []
        for itm in self.umb:
            self.dat.append( itm.get_func( self.mol ) )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.mol.natm * 3 ) ]
        self.emm.get_grad( self.mol )
        self.eqm.get_grad( self.mol )
        self.cvs.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = []
        for i in self.sele:
            i3 = i * 3
            self.grad += self.mol.grad[i3:i3+3]



ncrd = 2
nwin = 64

node, ncpu = qm3.utils._mpi.init()
size = nwin // ncpu

cwd = os.getcwd()

obj = []
dyn = []
for j in range( size ):
    os.mkdir( "%s/scratch_%02d"%( cwd, size*node+j ) )
    os.chdir( "%s/scratch_%02d"%( cwd, size*node+j ) )
    obj.append( my_problem( size*node+j ) )
    qm3.actions.dynamics.assign_velocities( obj[j], temperature = 300., project = True )
    dyn.append( qm3.actions.dynamics.langevin_verlet( obj[j], step_size = 0.001,
        temperature = 300., gamma_factor = 100., print_frequency = 10,
        project = True, log_function = obj[j].logging ) )


def take_step( distrib = True ):
    for j in range( size ):
        os.chdir( "%s/scratch_%02d"%( cwd, size*node+j ) )
        dyn[j].integrate()
    qm3.utils._mpi.barrier()
    if( node == 0 ):
        ncrd2 = ncrd * ncrd
        tmp_c = []
        tmp_m = []
        for j in range( size ):
            tmp_c += obj[j].cvs.rcrd[:]
            tmp_m += obj[j].cvs.cmet[:]
        for i in range( 1, ncpu ):
            for j in range( size ):
                tmp_c += qm3.utils._mpi.recv_r8( i, ncrd )
                tmp_m += qm3.utils._mpi.recv_r8( i, ncrd2 )
        if( distrib ):
            tmp_c = qm3.actions.string.string_distribute( ncrd, nwin, tmp_c, tmp_m )[0]
        for j in range( size ):
            obj[j].cvs.rcrd = tmp_c[j*ncrd:(j+1)*ncrd][:]
        for i in range( 1, ncpu ):
            for j in range( size ):
                qm3.utils._mpi.send_r8( i, tmp_c[(size*i+j)*ncrd:(size*i+j+1)*ncrd] )
        obj[0].cvs.fstr.write( "".join( [ "%20.10lf"%( tmp_c[j] ) for j in range( ncrd * nwin ) ] ) + "\n" )
        obj[0].cvs.fstr.flush()
        tmp_a = []
        tmp_b = []
        for i in range( nwin ):
            tmp_i = qm3.maths.matrix.inverse( [ tmp_m[i*ncrd2+j] for j in range( ncrd2 ) ], ncrd, ncrd )
            tmp_a += [ tmp_c[i*ncrd+j] - obj[0].cvs.icrd[i*ncrd+j] for j in range( ncrd ) ]
            tmp_b += qm3.maths.matrix.mult( tmp_i, ncrd, ncrd, tmp_a[i*ncrd:(i+1)*ncrd], ncrd, 1 )
        obj[0].cvs.fcnv.write( "%20.10lf\n"%( math.sqrt( sum( [ tmp_a[i] * tmp_b[i]
            for i in range( ncrd * nwin ) ] ) / float( nwin ) ) ) )
        obj[0].cvs.fcnv.flush()
    else:
        for j in range( size ):
            qm3.utils._mpi.send_r8( 0, obj[j].cvs.rcrd )
            qm3.utils._mpi.send_r8( 0, obj[j].cvs.cmet )
        for j in range( size ):
            obj[j].cvs.rcrd = qm3.utils._mpi.recv_r8( 0, ncrd )


for _ in range( 500 ):
    take_step( False )
for j in range( size ):
    if( size*node+j > 0 and size*node+j < ncpu - 1 ):
        obj[j].cvs.tstp = 0.001
for _ in range( 1000 ):
    take_step()


for j in range( size ):
    dyn[j].stats()
    os.chdir( "%s/scratch_%02d"%( cwd, size*node+j ) )
    obj[j].mol.pdb_write( "../04.string.%02d.pdb"%( size*node+j ) )
    obj[j].cvs.stop()

qm3.utils._mpi.barrier()
qm3.utils._mpi.stop()

#
# mpirun -n 8 python 03-pneb.py
#
import  qm3.mol
import  qm3.elements
import  qm3.problem
import  qm3.engines.mol_mech
import  qm3.engines.xtb
import  qm3.actions.minimize
import  qm3.actions.neb
import  qm3.utils._mpi
import  os
import  sys
import  pickle


class envi_problem( qm3.problem.template ):
    def __init__( self, molec ):
        qm3.problem.template.__init__( self )

        self.mol = molec

        self.emm = qm3.engines.mol_mech.simple_force_field( self.mol )
        self.emm.cut_on   = 12.0
        self.emm.cut_off  = 14.0
        self.emm.cut_list = 18.0
        self.emm.system_read( "../02.mm_data.pk", self.mol )

        f = open( "../02.sele_MM.pk", "rb" )
        self.sele = pickle.load( f )
        f.close()

        f = open( "../01.sele_QM.pk", "rb" )
        for i in pickle.load( f ):
            self.emm.free[i] = False
        f.close()
        f = open( "../02.fixed.pk", "rb" )
        for i in pickle.load( f ):
            self.emm.free[i] = False
        f.close()

        self.size = len( self.sele ) * 3
        self.coor = []
        for i in self.sele:
            i3 = i * 3
            self.coor += self.mol.coor[i3:i3+3]

        self.log = open( "log", "wt" )


    def feed_log( self, txt ):
        self.log.write( txt + "\n" )


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
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.mol.natm * 3 ) ]
        self.emm.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = []
        for i in self.sele:
            i3 = i * 3
            self.grad += self.mol.grad[i3:i3+3]



class core_problem( qm3.problem.template ):
    def __init__( self ):
        qm3.problem.template.__init__( self )

        self.mol = qm3.mol.molecule( "../03.pes.02.15.pdb" )
        self.mol.boxl = [ 40., 40., 40. ]

        self.envi = envi_problem( self.mol )

        f = open( "../01.sele_QM.pk", "rb" )
        self.sele = pickle.load( f )
        f.close()
        f = open( "../02.sele_MM.pk", "rb" )
        smm = pickle.load( f )
        f.close()

        self.eqm = qm3.engines.xtb.dl_xtb( self.mol, 0, 0, self.sele, smm )

        self.size = len( self.sele ) * 3
        self.coor = []
        self.mass = []
        for i in self.sele:
            i3 = i * 3
            self.coor += self.mol.coor[i3:i3+3]
            self.mass.append( qm3.elements.mass[self.mol.anum[i]] )

        self.mm_flg = False
        self.get_func()
        self.mm_flg = True


    def neb_data( self, who ):
        f = open( "../03.neb.%02d.pdb"%( who ), "wt" )
        f.write( "REMARK: %.6lf\n"%( self.func ) )
        self.mol.pdb_write( f )
        f.close()


    def update_coor( self ):
        for i in range( len( self.sele ) ):
            i3 = i * 3
            I3 = self.sele[i] * 3
            for j in [0, 1, 2]:
                self.mol.coor[I3+j] = self.coor[i3+j]
        if( self.mm_flg ):
            self.envi.feed_log( str( self.mol.chrg[0:17] ) )
            qm3.actions.minimize.fire( self.envi, print_frequency = 10, gradient_tolerance = 0.2,
                step_size = 0.1, step_number = 200, log_function = self.envi.feed_log )


    def get_func( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.eqm.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.mol.natm * 3 ) ]
        self.eqm.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = []
        for i in self.sele:
            i3 = i * 3
            self.grad += self.mol.grad[i3:i3+3]





node, ncpu = qm3.utils._mpi.init()
if( node == 0 ):
    def gen_log( txt ):
        sys.stdout.write( txt + "\n" )
        sys.stdout.flush()
else:
    def gen_log( txt ):
        pass

os.mkdir( "scratch_%02d"%( node ) )
os.chdir( "scratch_%02d"%( node ) )

f = open( "../01.sele_QM.pk", "rb" )
sel = pickle.load( f )
f.close()

crd = []
mol = qm3.mol.molecule( "../03.pes.02.15.pdb" )
tmp = []
for i in sel:
    tmp += mol.coor[3*i:3*i+3][:]
crd.append( tmp[:] )

mol.pdb_read( "../03.pes.16.01.pdb" )
tmp = []
for i in sel:
    tmp += mol.coor[3*i:3*i+3][:]
crd.append( tmp[:] )

crd = qm3.actions.neb.distribute( 24, crd )

obj = qm3.actions.neb.parall_neb( core_problem(), sel, 200., node, ncpu, crd )
qm3.actions.minimize.fire( obj, step_number = 1000, print_frequency = 1, step_size = 0.1,
        gradient_tolerance = 0.1, log_function = gen_log, fire2 = True, exit_uphill = True )

qm3.utils._mpi.barrier()

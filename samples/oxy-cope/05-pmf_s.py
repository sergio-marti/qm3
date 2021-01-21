import  sys
import  qm3.mol
import  qm3.problem
import  qm3.elements
import  qm3.engines.mol_mech
import  qm3.engines.xtb
import  qm3.engines.colvar_s
import  qm3.actions.dynamics
import  os
import  pickle


# ---------------------------------------------------
# >> problem config

class my_problem( qm3.problem.template ):
    def __init__( self, node ):
        qm3.problem.template.__init__( self )

        self.mol = qm3.mol.molecule( "../05.pmf_s.seed.%02d.pdb"%( node ) )
        self.mol.boxl = [ 40., 40., 40. ]

        self.emm = qm3.engines.mol_mech.simple_force_field( self.mol, 8 )
        self.emm.cut_on   = 12.0
        self.emm.cut_off  = 14.0
        self.emm.cut_list = 18.0
        self.emm.system_read( "../02.mm_data.pk", self.mol )

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

        self.umb = []
        self.dat = []

        self.size = len( self.sele ) * 3
        self.coor = []
        self.mass = []
        for i in self.sele:
            i3 = i * 3
            self.coor += self.mol.coor[i3:i3+3]
            self.mass.append( qm3.elements.mass[self.mol.anum[i]] )

        kb = 2000
        rx = node * 0.110703125  # s: 7.602 / 0-63 windows
        self.umb = qm3.engines.colvar_s.colvar_s( kb, rx,
                "../05.pmf_s.cnf", "../05.pmf_s.str", "../05.pmf_s.met", self.mol )

        self.fdat = open( "../05.dat_s.%02d"%( node ), "wt" )
        self.fdat.write( "%20.10lf%20.10lf\n"%( kb, rx ) )
        self.fdat.flush()


    def current_step( self, istep ):
        if( istep > 2000 ):
            self.fdat.write( "%20.10lf\n"%( self.dat ) )
            self.fdat.flush()
        if( istep % 1000 == 0 ):
            self.emm.update_non_bonded( self.mol )


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
        self.dat = self.umb.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.mol.natm * 3 ) ]
        self.emm.get_grad( self.mol )
        self.eqm.get_grad( self.mol )
        self.dat = self.umb.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = []
        for i in self.sele:
            i3 = i * 3
            self.grad += self.mol.grad[i3:i3+3]




node = int( sys.argv[1] )

os.mkdir( "scratch_%02d"%( node ) )
os.chdir( "scratch_%02d"%( node ) )

ff = open( "log", "wt" )

def feed_log( txt ):
    ff.write( txt + "\n" )


obj = my_problem( node )
qm3.actions.dynamics.assign_velocities( obj, 300. )
qm3.actions.dynamics.langevin_verlet( obj, step_size = 0.001, temperature = 300.0,
    gamma_factor = 100.0, print_frequency = 100, step_number = 6000, log_function = feed_log )
obj.mol.pdb_write( "../05.pmf_s.%02d.pdb"%( node ) )

ff.close()

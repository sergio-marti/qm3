import  qm3.mol
import  qm3.elements
import  qm3.problem
import  qm3.engines.mol_mech
import  qm3.engines.xtb
import  qm3.actions.minimize
import  os
import  pickle


def nolog( txt ):
    pass


class envi_problem( qm3.problem.template ):
    def __init__( self, molec ):
        qm3.problem.template.__init__( self )

        self.mol = molec

        self.emm = qm3.engines.mol_mech.simple_force_field( self.mol )
        self.emm.cut_on   = 12.0
        self.emm.cut_off  = 14.0
        self.emm.cut_list = 18.0
        self.emm.system_read( "02.mm_data.pk", self.mol )

        f = open( "02.sele_MM.pk", "rb" )
        self.sele = pickle.load( f )
        f.close()

        f = open( "01.sele_QM.pk", "rb" )
        for i in pickle.load( f ):
            self.emm.free[i] = False
        f.close()
        f = open( "02.fixed.pk", "rb" )
        for i in pickle.load( f ):
            self.emm.free[i] = False
        f.close()

        self.size = len( self.sele ) * 3
        self.coor = []
        for i in self.sele:
            i3 = i * 3
            self.coor += self.mol.coor[i3:i3+3]


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

        self.mol = qm3.mol.molecule( "03.neb.11.pdb" )
        self.mol.boxl = [ 40., 40., 40. ]

        self.envi = envi_problem( self.mol )

        f = open( "01.sele_QM.pk", "rb" )
        self.sele = pickle.load( f )
        f.close()
        f = open( "02.sele_MM.pk", "rb" )
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

        self.mm_flg = True
        self.cnt = 0


    def update_coor( self ):
        for i in range( len( self.sele ) ):
            i3 = i * 3
            I3 = self.sele[i] * 3
            for j in [0, 1, 2]:
                self.mol.coor[I3+j] = self.coor[i3+j]
        if( self.mm_flg ):
            qm3.actions.minimize.fire( self.envi, print_frequency = 10, gradient_tolerance = 0.5,
                step_size = 0.1, step_number = 200, log_function = nolog )


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


    def get_hess( self ):
        if( not os.path.isfile( "update.dump" ) or self.cnt % 50 == 0 ):
            self.update_coor()
            self.mm_flg = False
            self.num_hess( central = True )
            self.mm_flg = True
            qm3.utils.manage_hessian( self.coor, self.grad, self.hess, should_update = False )
        else:
            self.get_grad()
            self.hess = [ .0 for i in range( self.size * self.size ) ]
            qm3.utils.manage_hessian( self.coor, self.grad, self.hess, should_update = True )
        qm3.utils.raise_hessian_RT( self.mass, self.coor, self.hess )
        self.cnt += 1




obj = core_problem()

obj.mm_flg = False
obj.get_func()
obj.mm_flg = True

qm3.actions.minimize.baker( obj, print_frequency = 1,
        step_size = 0.1, gradient_tolerance = 1., follow_mode = 0 )
obj.mol.pdb_write( "03.saddle.pdb" )

smb = [ qm3.elements.symbol[obj.mol.anum[i]] for i in obj.sele ]
val, vec = qm3.utils.hessian_frequencies( obj.mass, obj.coor, obj.hess )
print( val )
for i in range( min( 10, obj.size ) ):
    if( val[i] < .0 ):
        qm3.utils.normal_mode_view( obj.coor, val, vec, smb, i, afac = 4 )
    else:
        qm3.utils.normal_mode_view( obj.coor, val, vec, smb, i )

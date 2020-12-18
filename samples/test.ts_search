import qm3.problem
import qm3.fio.dcd
import qm3.fio.xplor
import qm3.mol
import qm3.utils
import qm3.elements
import qm3.actions.minimize
import qm3.engines.namd
import qm3.engines.dftb
import qm3.engines.mmint
import time
import os


class MM_problem( qm3.problem.template ):

    def __init__( self, molec, nbn ):
        qm3.problem.template.__init__( self )

        self.mol = molec

        try:
            os.unlink( "namd.pipe" )
        except:
            pass
        os.mkfifo( "namd.pipe" )
        os.system( "bash r.namd &" )
        time.sleep( 10 )
        self.eng = qm3.engines.namd.namd_pipe()

        self.sel  = nbn
        self.coor = sum( [ self.mol.coor[3*i:3*i+3] for i in self.sel ], [] )
        self.size = len( self.coor )
        self.grad = []
        self.hess = []


    def update_coor( self ):
        for i in range( len( self.sel ) ):
            for j in [0, 1, 2]:
                self.mol.coor[3*self.sel[i]+j] = self.coor[3*i+j]


    def get_func( self ):
        self.update_coor()
        self.eng.update_coor( self.mol )
        self.mol.func = .0
        self.eng.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.eng.update_coor( self.mol )
        self.mol.func = .0
        self.mol.grad = [ .0 for i in range( 3 * self.mol.natm ) ]
        self.eng.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = sum( [ self.mol.grad[3*i:3*i+3] for i in self.sel ], [] )



class QM_problem( qm3.problem.template ):

    def __init__( self, molec ):
        qm3.problem.template.__init__( self )

        self.mol = molec

        self.sel = sorted( self.mol.indx["A"][1].values() )
        self.nbn = list( set( list( range( self.mol.natm ) ) ).difference( set( self.sel ) ) )
        for i in self.sel:
            self.mol.chrg[i] = .0
        self.coor = sum( [ self.mol.coor[3*i:3*i+3] for i in self.sel ], [] )
        self.mass = [ self.mol.mass[i] for i in self.sel ]
        self.size = len( self.coor )

        self.eng = qm3.engines.dftb.dftb( self.mol, self.sel, self.nbn )
        self.eng.chg = 0
        self.eng.prm = "./bin/3ob-3-1/"
        self.eng.exe = "bash r.dftb"

        self.fix = qm3.engines.mmint.QMLJ( self.mol, self.sel, self.nbn, [] )

        self.grad = []
        self.hess = []

        self.pMM = MM_problem( molec, self.nbn )

        self.cnt = 0

        self.dcd = qm3.fio.dcd.dcd()
        self.dcd.open_write( "dcd", self.mol.natm )


    def current_step( self, istep ):
        self.dcd.append( self.mol )


    def update_coor( self ):
        for i in range( len( self.sel ) ):
            for j in [0, 1, 2]:
                self.mol.coor[3*self.sel[i]+j] = self.coor[3*i+j]


    def get_func( self ):
        self.update_coor()
        self.mol.func = .0
        self.eng.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = .0
        self.mol.grad = [ .0 for i in range( 3 * self.mol.natm ) ]
        self.eng.get_grad( self.mol )
        self.fix.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = sum( [ self.mol.grad[3*i:3*i+3] for i in self.sel ], [] )


    def get_hess( self ):
        self.update_coor()
        print( 100*"#" )
        print( [ self.mol.chrg[i] for i in self.sel ] )
        self.pMM.eng.update_chrg( self.mol )
        qm3.actions.minimize.fire( self.pMM, step_number = 1000, print_frequency = 100, gradient_tolerance = 1. )
        print( 100*"#" )
        # ---------------------------------------------
        if( self.cnt % 10 == 0 ):
            self.num_hess( central = True )
            qm3.utils.manage_hessian( self.coor, self.grad, self.hess, should_update = False )
        else:
            self.mol.func = .0
            self.mol.grad = [ .0 for i in range( 3 * self.mol.natm ) ]
            self.eng.get_grad( self.mol )
            self.fix.get_grad( self.mol )
            self.func = self.mol.func
            self.grad = sum( [ self.mol.grad[3*i:3*i+3] for i in self.sel ], [] )
            self.hess = [ .0 for i in range( self.size * self.size ) ]
            qm3.utils.manage_hessian( self.coor, self.grad, self.hess, should_update = True )
        self.cnt += 1
        # ---------------------------------------------
        qm3.utils.raise_hessian_RT( self.mass, self.coor, self.hess )



mol = qm3.mol.molecule( "pdb" )
mol.boxl = [ 40.0, 40.0, 40.0 ]
qm3.fio.xplor.psf_read( self.mol, "psf" )
qm3.engines.mmint.non_bonded( self.mol, "non_bonded" )


obj = QM_problem( mol )

obj.get_func()
qm3.actions.minimize.steepest_descent( obj.pMM, step_number = 1000, step_size = 1.,
        print_frequency = 100, gradient_tolerance = 100 )
qm3.actions.minimize.fire( obj.pMM, step_number = 1000, step_size = 1.0,
        print_frequency = 100, gradient_tolerance = 1.0 )

qm3.actions.minimize.baker( obj, step_number = 100, step_size = 0.1,
        print_frequency = 1, gradient_tolerance = 2.0, follow_mode = 0 )
mol.pdb_write( "out.pdb" )
obj.dcd.close()
obj.pMM.eng.stop()

val, vec = qm3.utils.hessian_frequencies( obj.mass, obj.coor, obj.hess )
print( val )
for i in range( min( 10, obj.size ) ):
    if( val[i] < .0 ):
        qm3.utils.normal_mode_view( obj.coor, val, vec, obj.eng.smb, i, afac = 4 )
    else:
        qm3.utils.normal_mode_view( obj.coor, val, vec, obj.eng.smb, i )

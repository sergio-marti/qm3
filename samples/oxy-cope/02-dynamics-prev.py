import  qm3.mol
import  qm3.problem
import  qm3.elements
import  qm3.engines.mol_mech
import  qm3.engines.xtb
import  qm3.actions.minimize
import  os
import  time
import  pickle


class my_problem( qm3.problem.template ):
    def __init__( self ):
        qm3.problem.template.__init__( self )

        self.mol = qm3.mol.molecule( "01.pdb" )
        self.mol.boxl = [ 40., 40., 40. ]
        self.mol.anum = [ qm3.elements.rsymbol[self.mol.labl[i][0]] for i in range( self.mol.natm ) ]
        self.mol.chrg = [ 0.0 ] * 17 + [ -0.834, 0.417, 0.417 ] * ( ( self.mol.natm - 17 ) // 3 )

        self.emm = qm3.engines.mol_mech.simple_force_field( self.mol, 32 )
        print( "ncpu:", self.emm.ncpu )
        self.emm.cut_on   = 12.0
        self.emm.cut_off  = 14.0
        self.emm.cut_list = 18.0
        self.emm.topology( self.mol, qchg = False )
        self.emm.parameters( self.mol )

        self.emm.system_write( "02.mm_data.pk", self.mol )

        f = open( "01.sele_QM.pk", "rb" )
        sqm = pickle.load( f )
        f.close()
        f = open( "01.sele_MM.pk", "rb" )
        smm = pickle.load( f )
        f.close()

        self.emm.qm_atoms( sqm )

        self.eqm = qm3.engines.xtb.dl_xtb( self.mol, 0, 0, sqm, smm )

        self.size = self.mol.natm * 3
        self.coor = self.mol.coor


    def current_step( self, istep ):
        if( ( istep + 1 ) % 100 == 0 ):
            self.emm.update_non_bonded( self.mol )


    def get_func( self ):
        self.mol.func = 0.0
        self.emm.get_func( self.mol )
        self.eqm.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.mol.natm * 3 ) ]
        self.emm.get_grad( self.mol )
        self.eqm.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = self.mol.grad[:]




obj = my_problem()

qm3.actions.minimize.fire( obj, print_frequency = 100, 
    gradient_tolerance = 10.0, step_number = 1000, fire2 = True )

obj.mol.pdb_write( "02.prev.pdb" )

sqm = list( sorted( obj.mol.indx["A"][1].values() ) )

sel = obj.mol.sph_sel( sqm, 14.0 )
f = open( "02.sele.pk", "wb" )
pickle.dump( sel, f )
f.close()

smm = list( set( sel ).difference( set( sqm ) ) )
f = open( "02.sele_MM.pk", "wb" )
pickle.dump( smm, f )
f.close()

fix = list( set( range( obj.mol.natm ) ).difference( set( sel ) ) )
f = open( "02.fixed.pk", "wb" )
pickle.dump( fix, f )
f.close()

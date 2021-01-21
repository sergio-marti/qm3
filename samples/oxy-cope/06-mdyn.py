import  sys
import  qm3.mol
import  qm3.fio.dcd
import  qm3.problem
import  qm3.elements
import  qm3.engines.mol_mech
import  qm3.engines.xtb
import  qm3.engines.mmres
import  qm3.engines._plumed
import  qm3.actions.dynamics
import  os
import  pickle


class my_problem( qm3.problem.template ):
    def __init__( self ):
        qm3.problem.template.__init__( self )

        self.mol = qm3.mol.molecule( "03.pes.02.14.pdb" )
        self.mol.boxl = [ 40., 40., 40. ]

        self.emm = qm3.engines.mol_mech.simple_force_field( self.mol, 8 )
        self.emm.cut_on   = 12.0
        self.emm.cut_off  = 14.0
        self.emm.cut_list = 18.0
        self.emm.system_read( "02.mm_data.pk", self.mol )
        self.mol.fill_masses()

        f = open( "01.sele_QM.pk", "rb" )
        sqm = pickle.load( f )
        f.close()
        f = open( "02.sele_MM.pk", "rb" )
        smm = pickle.load( f )
        f.close()
        f = open( "02.sele.pk", "rb" )
        self.sele = pickle.load( f )
        f.close()
        f = open( "02.fixed.pk", "rb" )
        for i in pickle.load( f ):
            self.emm.free[i] = False
        f.close()

        self.emm.qm_atoms( sqm )

        self.eqm = qm3.engines.xtb.dl_xtb( self.mol, 0, 0, sqm, smm )

        self.rst = qm3.engines._plumed.Plumed( self.mol, 0.001, False )
        # -- walls
        self.umb = []
        self.umb.append( qm3.engines.mmres.distance( 10000., 1.4, [ 1, 2 ], skip_BE = 1.4 ) )
        self.umb.append( qm3.engines.mmres.distance( 10000., 4.2, [ 1, 2 ], skip_LE = 4.2 ) )
        self.umb.append( qm3.engines.mmres.distance( 10000., 1.4, [ 4, 5 ], skip_BE = 1.4 ) )
        self.umb.append( qm3.engines.mmres.distance( 10000., 4.2, [ 4, 5 ], skip_LE = 4.2 ) )

        self.size = len( self.sele ) * 3
        self.coor = []
        self.mass = []
        for i in self.sele:
            i3 = i * 3
            self.coor += self.mol.coor[i3:i3+3]
            self.mass.append( self.mol.mass[i] )

        self.dcd = qm3.fio.dcd.dcd()
        self.dcd.open_write( "06.dcd", self.mol.natm )

        self.time = -1



    def current_step( self, istep ):
        if( istep % 100 == 0 ):
            self.dcd.append( self.mol )
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
        for itm in self.umb:
            itm.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.mol.natm * 3 ) ]
        self.emm.get_grad( self.mol )
        self.eqm.get_grad( self.mol )
        self.time += 1
        self.rst.get_grad( self.mol, self.time )
        for itm in self.umb:
            itm.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = []
        for i in self.sele:
            i3 = i * 3
            self.grad += self.mol.grad[i3:i3+3]



f = open( "plumed.dat", "wt" )
f.write( """
d1:  DISTANCE ATOMS=2,3 NOPBC
d2:  DISTANCE ATOMS=5,6 NOPBC
ant: COMBINE ARG=d1,d2 COEFFICIENTS=1.0,-1.0 PERIODIC=NO
METAD ...
LABEL=met
ARG=ant
SIGMA=0.005
HEIGHT=1.2
PACE=10
FILE=hills
... METAD
PRINT ARG=d1,d2,ant,met.bias FILE=colvar STRIDE=1
""" )
f.close()


obj = my_problem()
qm3.actions.dynamics.assign_velocities( obj, 300. )
qm3.actions.dynamics.langevin_verlet( obj, step_size = 0.001, temperature = 300.0,
    gamma_factor = 100.0, print_frequency = 100, step_number = 100000 )

obj.mol.pdb_write( "06.pdb" )
obj.dcd.close()

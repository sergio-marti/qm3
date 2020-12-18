import sys
import qm3.actions.dynamics
import qm3.engines.mmres
import qm3.engines._plumed
import qm3.engines.mmint
import qm3.problem
import qm3.mol
import qm3.fio.dcd
import qm3.engines.namd
import qm3.engines.dftb
import time
import os
try:
    import cPickle as pickle
except:
    import pickle


class my_problem( qm3.problem.template ):

    def __init__( self ):
        qm3.problem.template.__init__( self )

        os.mkfifo( "namd.pipe" )
        os.system( "bash r.namd &" )
        time.sleep( 10 )

        self.mol = qm3.mol.molecule( "pdb" )
        self.mol.boxl = [ 40., 40., 40. ]
        self.mol.nbnd_read( "non_bonded" )
        qm3.engines.mmint.non_bonded( self.mol, "non_bonded" )

        self.emm = qm3.engines.namd.namd_pipe()
        s_qm = sorted( self.mol.indx["A"][1].values() )
        for i in s_qm:
            self.mol.chrg[i] = 0.0
        f = open( "in16", "rb" )
        s_mm = pickle.load( f )
        f.close()

        self.eqm = qm3.engines.dftb.dftb( self.mol, s_qm, s_mm )
        self.eqm.chg = 0
        self.eqm.prm = "./bin/3ob-3-1/"
        self.eqm.exe = "bash r.dftb"

        self.fix = qm3.engines.mmint.QMLJ( self.mol, s_qm, s_mm, [] )

        self.rst = qm3.engines._plumed.Plumed( self.mol, 0.001, False )

        # -- walls
        self.ld1 = qm3.engines.mmres.distance( 10000., 1.4, [ 1, 2 ], skip_BE = 1.4 )
        self.ud1 = qm3.engines.mmres.distance( 10000., 4.2, [ 1, 2 ], skip_LE = 4.2 )
        self.ld2 = qm3.engines.mmres.distance( 10000., 1.4, [ 4, 5 ], skip_BE = 1.4 )
        self.ud2 = qm3.engines.mmres.distance( 10000., 4.2, [ 4, 5 ], skip_LE = 4.2 )

        self.size = self.mol.natm * 3
        self.coor = self.mol.coor[:]
        self.grad = []
        self.mass = self.mol.mass

        self.time = -1
        self.dcd = qm3.fio.dcd.dcd()
        self.dcd.open_write( "dcd", self.mol.natm )
        self.mol.dcd_write( self.dcd )

        self.emm.update_chrg( self.mol )


    def current_step( self, istep ):
        if( istep%100 == 0 ):
            self.dcd.append( self.mol )


    def get_func( self ):
        self.mol.coor = self.coor[:]
        self.emm.update_coor( self.mol )
        self.mol.func = .0
        self.emm.get_func( self.mol )
        self.eqm.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.mol.coor = self.coor[:]
        self.emm.update_coor( self.mol )
        self.mol.func = .0
        self.mol.grad = [ .0 for i in range( 3 * self.mol.natm ) ]
        self.emm.get_grad( self.mol )
        self.eqm.get_grad( self.mol )
        self.fix.get_grad( self.mol )
        self.time += 1
        self.rst.get_grad( self.mol, self.time )
        self.ld1.get_grad( self.mol )
        self.ud1.get_grad( self.mol )
        self.ld2.get_grad( self.mol )
        self.ud2.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = self.mol.grad[:]



f = open( "namd.inp", "wt" )
f.write( """
structure           psf
coordinates         pdb
paraTypeCharmm      on
parameters          wat.prm
parameters          acs.prm
fixedatoms          on
fixedatomsfile      pdb
cellBasisVector1    40. .0 .0
cellBasisVector2    .0 40. .0
cellBasisVector3    .0 .0 40.
PME                 on
PMETolerance        0.000001
PMEGridSpacing      0.5
exclude             scaled1-4
1-4scaling          0.5
switching           on
switchdist          12.0
cutoff              14.0
pairlistdist        16.0
wrapAll             off
wrapWater           off
nonbondedFreq       1
fullElectFrequency  1
stepspercycle       1
temperature         0.0
outputEnergies      1
outputname          namd.out
startup
############################################################
set fd [ open "namd.pipe" r ]
while { [ gets $fd cmd ] >= 0 } {
    switch $cmd {
        "energy"      { run 0 }
        "gradient"    { run 0; output onlyforces namd }
        "charges"     { reloadCharges namd.chrg }
        "coordinates" { coorfile binread namd.coor }
        "exit"        { close $fd; exit }
    }
}
""" )
f.close()


f = open( "r.namd", "wt" )
f.write( "./bin/namd2 +setcpuaffinity +isomalloc_sync +idlepoll namd.inp > namd.out\n" )
f.close()


f = open( "r.dftb", "wt" )
f.write( """
export OMP_NUM_THREADS=4
./bin/dftb+ > dftb_in.log
""" )
f.close()


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

t0 = time.time()
obj = my_problem()
qm3.actions.dynamics.assign_velocities( obj, 300.0, project = True )
dyn = qm3.actions.dynamics.langevin_verlet( obj, step_size = 0.001, print_frequency = 10,
    step_number = 100000, project = True )
obj.dcd.close()
obj.emm.stop()
print( time.time() - t0 )

import qm3.mol
import qm3.fio.xplor
import qm3.problem
import qm3.utils._mpi
import qm3.actions.dynamics
import qm3.engines.namd
import qm3.engines.dftb
import qm3.engines.mmint
import qm3.engines.mmres
import os
import time
try:
    import cPickle as pickle
except:
    import pickle



class my_problem( qm3.problem.template ):

    def __init__( self ):
        qm3.problem.template.__init__( self )

        self.node, self.ncpu = qm3.utils._mpi.init()

        os.mkdir( str( self.node ) )
        os.chdir( str( self.node ) )

        f = open( "namd.inp", "wt" )
        f.write( """structure           ../psf
coordinates         ../seed
paraTypeCharmm      on
parameters          ../wat.prm
parameters          ../acs.prm
fixedatoms          on
fixedatomsfile      ../pdb
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
        os.mkfifo( "namd.pipe" )
        os.system( "bash ../r.namd &" )

        self.mol = qm3.mol.molecule( "../seed" )
        self.mol.boxl = [ 40., 40., 40. ]
        qm3.fio.xplor.psf_read( self.mol, "../psf" )
        qm3.engines.mmint.non_bonded( self.mol, "../non_bonded" )

        self.mass = self.mol.mass

        self.emm = qm3.engines.namd.namd_pipe()
        s_qm = sorted( self.mol.indx["A"][1].values() )
        for i in s_qm:
            self.mol.chrg[i] = 0.0
        f = open( "../in16", "rb" )
        s_mm = pickle.load( f )
        f.close()

        self.eqm = qm3.engines.dftb.dftb( self.mol, s_qm, s_mm )
        self.eqm.chg = 0
        self.eqm.prm = "./bin/3ob-3-1/"
        self.eqm.exe = "bash ../r.dftb"

        self.fix = qm3.engines.mmint.QMLJ( self.mol, s_qm, s_mm, [] ) 

        self.size = self.mol.natm * 3
        self.coor = self.mol.coor[:]
        self.grad = []
        self.hess = []

        kmb = 2400.
        # -1.8 ... 1.7 / 0 ... 35 windows
        ref = -1.8 + 0.1 * float( self.node )
        idx = [ self.mol.indx["A"][1]["C2"], self.mol.indx["A"][1]["C3"],
                self.mol.indx["A"][1]["C5"], self.mol.indx["A"][1]["C6"] ]
        wei = [ 1., -1. ]
        self.dat = open( "pmf.%03d.dat"%( self.node ), "wt" )
        self.dat.write( "%20.10lf%20.10lf\n"%( kmb, ref ) )
        self.rst = qm3.engines.mmres.multiple_distance( kmb, ref, idx, wei )
        self.cnt = 0

        self.emm.update_chrg( self.mol )

        self.flg = open( "pmf.%03d.log"%( self.node ), "wt" )
        self.log( str( sum( self.mass ) ) + " g/mol" )
        self.log( str( self.node ) + " / " + str( kmb ) + " / " + str( ref ) + " / " + str( wei ) + " / " + str( idx ) )


    def log( self, txt ):
        self.flg.write( txt + "\n" )
        self.flg.flush()


    def get_func( self ):
        self.mol.coor = self.coor[:]
        self.emm.update_coor( self.mol )
        self.mol.func = 0.0
        self.emm.get_func( self.mol )
        self.eqm.get_func( self.mol )
        self.rst.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.mol.coor = self.coor[:]
        self.emm.update_coor( self.mol )
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( 3 * self.mol.natm ) ]
        self.emm.get_grad( self.mol )
        self.eqm.get_grad( self.mol )
        self.fix.get_grad( self.mol )
        if( self.cnt > 2000 ):
            self.dat.write( "%20.10lf\n"%( self.rst.get_grad( self.mol ) ) )
        else:
            self.rst.get_grad( self.mol )
        self.cnt += 1
        self.func = self.mol.func
        self.grad = self.mol.grad[:]



obj = my_problem()
qm3.actions.dynamics.assign_velocities( obj, 300. )
qm3.actions.dynamics.langevin_verlet( obj, temperature = 300., gamma_factor = 50.,
    print_frequency = 10, project = True, log_function = obj.log, step_number = 4000 )
obj.dat.close()
obj.emm.stop()
obj.flg.close()
qm3.utils._mpi.barrier()
qm3.utils._mpi.stop()

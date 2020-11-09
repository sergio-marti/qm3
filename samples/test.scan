import os
import sys
import qm3.fio.xplor
import qm3.actions.minimize
import qm3.engines.mmres
import qm3.engines.mmint
import qm3.problem
import qm3.mol
import qm3.engines.namd
import qm3.engines.dftb
import qm3.utils
import time
try:
    import cPickle as pickle
except:
    import pickle


class my_problem( qm3.problem.template ):

    def __init__( self, old_i, old_j, new_i, new_j ):
        qm3.problem.template.__init__( self )

        os.mkdir( str( new_i ) + "." + str( new_j ) )
        os.chdir( str( new_i ) + "." + str( new_j ) )

        f = open( "namd.inp", "wt" )
        f.write( """structure           ../psf
coordinates         ../pes.%d.%d
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
"""%( old_i, old_j ) )
        f.close()
        os.mkfifo( "namd.pipe" )
        os.system( "bash ../r.namd &" )
        time.sleep( 10 )

        self.mol = qm3.mol.molecule( "../pes.%d.%d"%( old_i, old_j ) )
        self.mol.boxl = [ 40., 40., 40. ]
        qm3.fio.xplor.psf_read( self.mol, "../psf" )
        qm3.engines.mmint.non_bonded( self.mol, "../non_bonded" )

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

        self.rst_I = qm3.engines.mmres.distance( 5000., 1.50 + 0.075 * float( new_i ),
            [ self.mol.indx["A"][1]["C2"], self.mol.indx["A"][1]["C3"] ] )
        self.rst_J = qm3.engines.mmres.distance( 5000., 3.00 - 0.075 * float( new_j ),
            [ self.mol.indx["A"][1]["C5"], self.mol.indx["A"][1]["C6"] ] )

        self.size = self.mol.natm * 3
        self.coor = self.mol.coor[:]
        self.grad = []
        self.hess = []

        self.emm.update_chrg( self.mol )

        self.flg = open( "../log.%d.%d"%( new_i, new_j ), "wt" )


    def log( self, txt ):
        self.flg.write( txt + "\n" )
        self.flg.flush()


    def get_func( self ):
        self.mol.coor = self.coor[:]
        self.emm.update_coor( self.mol )
        self.mol.func = .0
        self.emm.get_func( self.mol )
        self.last_QM = self.mol.func
        self.eqm.get_func( self.mol )
        self.last_QM = self.mol.func - self.last_QM
        self.rst_I.get_func( self.mol )
        self.rst_J.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.mol.coor = self.coor[:]
        self.emm.update_coor( self.mol )
        self.mol.func = .0
        self.mol.grad = [ .0 for i in range( 3 * self.mol.natm ) ]
        self.emm.get_grad( self.mol )
        self.last_QM = self.mol.func
        self.eqm.get_grad( self.mol )
        self.last_QM = self.mol.func - self.last_QM
        self.fix.get_grad( self.mol )
        self.rst_I.get_grad( self.mol )
        self.rst_J.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = self.mol.grad[:]



old_i = int( sys.argv[1] )
old_j = int( sys.argv[2] )
new_i = int( sys.argv[3] )
new_j = int( sys.argv[4] )

obj = my_problem( old_i, old_j, new_i, new_j )

if( new_i == 0 and new_j == 0 ):
    qm3.actions.minimize.steepest_descent( obj, step_number = 100, print_frequency = 10,
        gradient_tolerance = 100., step_size = 1.0, log_function = obj.log )

qm3.actions.minimize.conjugate_gradient_plus( obj, step_number = 2000,
    print_frequency = 100, gradient_tolerance = 0.1, log_function = obj.log )

obj.mol.pdb_write( "../pes.%d.%d"%( new_i, new_j ) )

obj.log( "\n\n%20.10lf%20.10lf%40.10lf%20.10lf"%( 
    qm3.utils.distance( obj.coor[obj.rst_I.indx[0]*3:obj.rst_I.indx[0]*3+3],
        obj.coor[obj.rst_I.indx[1]*3:obj.rst_I.indx[1]*3+3] ),
    qm3.utils.distance( obj.coor[obj.rst_J.indx[0]*3:obj.rst_J.indx[0]*3+3],
        obj.coor[obj.rst_J.indx[1]*3:obj.rst_J.indx[1]*3+3] ), obj.func, obj.last_QM ) )

obj.flg.close()
obj.emm.stop()

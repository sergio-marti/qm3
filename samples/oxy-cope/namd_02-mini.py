import  qm3.mol
import  qm3.fio.xplor
import  qm3.problem
import  qm3.engines.namd
import  qm3.engines.xtb
import  qm3.engines.mmint
import  qm3.actions.minimize
import  os
import  time
import  pickle


class my_problem( qm3.problem.template ):
    def __init__( self ):
        qm3.problem.template.__init__( self )

        self.mol = qm3.mol.molecule( "01.pdb" )
        self.mol.boxl = [ 40., 40., 40. ]
        qm3.fio.xplor.psf_read( self.mol, "01.psf" )
        self.mol.guess_atomic_numbers()
        qm3.engines.mmint.non_bonded( self.mol, "01.non_bonded" )

        os.system( "rm -vf namd.*" )
        f = open( "namd.inp", "wt" )
        f.write( """
structure           01.psf
coordinates         01.pdb
paraTypeCharmm      on
parameters          01.wat.prm
parameters          01.acs.prm
fixedatoms          on
fixedatomsfile      01.pdb
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
pairlistdist        18.0
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
        "gradient"    { run 0; output onlyforces shm }
        "charges"     { reloadCharges shm }
        "coordinates" { coorfile shmread }
        "exit"        { close $fd; exit }
    }
}
""" )
        f.close()
        os.mkfifo( "namd.pipe" )
        os.system( "NAMD_SHM=1 ./bin/namd2 +ppn 1 +setcpuaffinity +isomalloc_sync +idlepoll namd.inp > namd.out &" )
        while( not os.path.isfile( "namd.shmid" ) ):
            time.sleep( 1 )
        time.sleep( 1 )
        self.emm = qm3.engines.namd.namd_shm()

        f = open( "01.sele_QM.pk", "rb" )
        sqm = pickle.load( f )
        f.close()
        f = open( "01.sele_MM.pk", "rb" )
        smm = pickle.load( f )
        f.close()

        self.eqm = qm3.engines.xtb.dl_xtb( self.mol, 0, 0, sqm, smm )

        self.fix = qm3.engines.mmint.QMLJ( self.mol, sqm, smm, [] )

        self.size = self.mol.natm * 3
        self.coor = self.mol.coor
        for i in sqm:
            self.mol.chrg[i] = 0.0
        self.emm.update_chrg( self.mol )


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
        self.fix.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = self.mol.grad[:]




obj = my_problem()

qm3.actions.minimize.fire( obj, print_frequency = 10, 
    gradient_tolerance = 10.0, step_number = 1000, fire2 = True )

sqm = list( sorted( obj.mol.indx["A"][1].values() ) )

sel = obj.mol.sph_sel( sqm, 14.0 )
f = open( "02.sele.pk", "wb" )
pickle.dump( sel, f )
f.close()

smm = list( set( sel ).difference( set( sqm ) ) )
f = open( "02.sele_MM.pk", "wb" )
pickle.dump( smm, f )
f.close()

fix = list( set( range( obj.mol.natm ) ).difference( set( smm ) ) )
f = open( "02.fixed.pk", "wb" )
pickle.dump( fix, f )
f.close()

qm3.engines.namd.pdb_write( obj.mol, "02.prev.pdb", fix )

obj.emm.stop()

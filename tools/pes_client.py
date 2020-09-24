#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-

import    sys
import    os
os.environ["QM3_LIBXTB"] = "/Users/smarti/Devel/qm3/xtb.so"

import    qm3.mol
import    qm3.engines.namd
import    qm3.engines.xtb
import    qm3.engines._mmint
import    qm3.engines.mmres
import    qm3.actions.minimize
import    qm3.problem
import    pickle
import    time
import    socket
import    struct


class my_problem( qm3.problem.template ):
    def __init__( self, who, cx, cy ):
        qm3.problem.template.__init__( self )

        self.mol = qm3.mol.molecule( who )
        self.mol.boxl = [ 80.30, 69.03, 75.00 ]
        self.mol.psf_read( "../psf" )
        self.mol.guess_atomic_numbers()
        self.mol.nbnd_read( "../non_bonded" )

        f = open( "namd.inp", "wt" )
        f.write( """
structure           ../psf
coordinates         %s
paraTypeCharmm      on
parameters          ../par_all22_prot.inp
parameters          ../9RQ.prm
fixedatoms          on
fixedatomsfile      ../pdb.fixed
cellBasisVector1    80.30   0.00   0.00
cellBasisVector2     0.00  69.03   0.00
cellBasisVector3     0.00   0.00  75.00
exclude             scaled1-4
1-4scaling          1.0
switching           on
switchdist          14.0
cutoff              16.0
pairlistdist        18.0
wrapAll             off
wrapWater           off
nonbondedFreq       1
fullElectFrequency  1
stepspercycle       1
temperature         0.0
outputEnergies      1
outputname          namd.out
run 0
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
"""%( who ) )
        f.close()
        os.mkfifo( "namd.pipe" )
        os.system( "NAMD_SHM=1 namd2 +setcpuaffinity +isomalloc_sync +idlepoll +ppn 1 namd.inp > namd.out &" )
        while( not os.path.isfile( "namd.shmid" ) ):
            time.sleep( 1 )
        time.sleep( 2 )

        self.emm = qm3.engines.namd.namd_shm()

        f = open( "../sele.pk", "rb" )
        self.sel = pickle.load( f )
        f.close()

        f = open( "../sele_QM.pk", "rb" )
        sqm = pickle.load( f )
        f.close()
        f = open( "../sele_LA.pk", "rb" )
        sla = pickle.load( f )
        f.close()
        f = open( "../sele_MM.pk", "rb" )
        smm = pickle.load( f )
        f.close()

        self.eqm = qm3.engines.xtb.dl_xtb( self.mol, 1, 0, sqm, smm, sla )

        f = open( "../sele_EX.pk", "rb" )
        exc = pickle.load( f )
        f.close()
        self.fix = qm3.engines._mmint.QMLJ( self.mol, sqm, smm, exc )

        self.umb = []
        # 0 -- 29
        self.umb.append( qm3.engines.mmres.multiple_distance( 5000., - 1.5 + 0.1 * cx, [ self.mol.indx["X"][1]["O2"], self.mol.indx["X"][1]["CX"], self.mol.indx["X"][1]["CX"], self.mol.indx["X"][2]["O"] ], [ 1.0, -1.0 ] ) )
        # 0 -- 43
        self.umb.append( qm3.engines.mmres.multiple_distance( 5000., - 2.5 + 0.1 * cy, [ self.mol.indx["X"][2]["O"], self.mol.indx["X"][2]["H2"], self.mol.indx["X"][2]["H2"], self.mol.indx["X"][1]["O2"] ], [ 1.0, -1.0 ] ) )

        self.exc = []
        # CB - CA || CT2 - CT1
        self.exc.append( qm3.engines.mmres.distance( 1861.9, 1.538, [ 525, 523 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 525 ), -1 ]
        # CB - CA || CT2 - CT1
        self.exc.append( qm3.engines.mmres.distance( 1861.9, 1.538, [ 912, 910 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 912 ), -1 ]
        # CB - CA || CT2 - CT1
        self.exc.append( qm3.engines.mmres.distance( 1861.9, 1.538, [ 983, 981 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 983 ), -1 ]
        #------------------------------------------------------------------
        # CB - CA - N || CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.angle( 585.8, 113.5, [ 525, 523, 521 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 525 ), -1, -1 ]
        # CB - CA - HA || CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.angle( 292.9, 111.0, [ 525, 523, 524 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 525 ), -1, -1 ]
        # CB - CA - C || CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.angle( 435.1, 108.0, [ 525, 523, 536 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 525 ), -1, -1 ]
        # CB - CA - N || CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.angle( 585.8, 113.5, [ 912, 910, 908 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 912 ), -1, -1 ]
        # CB - CA - HA || CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.angle( 292.9, 111.0, [ 912, 910, 911 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 912 ), -1, -1 ]
        # CB - CA - C || CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.angle( 435.1, 108.0, [ 912, 910, 923 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 912 ), -1, -1 ]
        # CB - CA - N || CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.angle( 585.8, 113.5, [ 983, 981, 979 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 983 ), -1, -1 ]
        # CB - CA - HA || CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.angle( 292.9, 111.0, [ 983, 981, 982 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 983 ), -1, -1 ]
        # CB - CA - C || CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.angle( 435.1, 108.0, [ 983, 981, 994 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 983 ), -1, -1 ]
        #------------------------------------------------------------------
        # CB - CA - N - C || CT2 - CT1 - NH1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 1: [ 15.062, 0.0 ], }, [ 525, 523, 521, 519 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 525 ), -1, -1, -1 ]
        # CB - CA - N - HN || CT2 - CT1 - NH1 - H
        # CB - CA - C - N || CT2 - CT1 - C - N
        # CB - CA - C - O || CT2 - CT1 - C - O
        self.exc.append( qm3.engines.mmres.dihedral( { 1: [ 11.715, 0.0 ], }, [ 525, 523, 536, 537 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 525 ), -1, -1, -1 ]
        # HB1 - CB - CA - N || HA - CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 526, 525, 523, 521 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 526 ), self.sel.index( 525 ), -1, -1 ]
        # HB1 - CB - CA - HA || HA - CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 526, 525, 523, 524 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 526 ), self.sel.index( 525 ), -1, -1 ]
        # HB1 - CB - CA - C || HA - CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 526, 525, 523, 536 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 526 ), self.sel.index( 525 ), -1, -1 ]
        # HB2 - CB - CA - N || HA - CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 527, 525, 523, 521 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 527 ), self.sel.index( 525 ), -1, -1 ]
        # HB2 - CB - CA - HA || HA - CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 527, 525, 523, 524 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 527 ), self.sel.index( 525 ), -1, -1 ]
        # HB2 - CB - CA - C || HA - CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 527, 525, 523, 536 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 527 ), self.sel.index( 525 ), -1, -1 ]
        # CG - CB - CA - N || CPH1 - CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 529, 525, 523, 521 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 529 ), self.sel.index( 525 ), -1, -1 ]
        # CG - CB - CA - HA || CPH1 - CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 529, 525, 523, 524 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 529 ), self.sel.index( 525 ), -1, -1 ]
        # CG - CB - CA - C || CPH1 - CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 529, 525, 523, 536 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 529 ), self.sel.index( 525 ), -1, -1 ]
        # CB - CA - N - C || CT2 - CT1 - NH1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 1: [ 15.062, 0.0 ], }, [ 912, 910, 908, 906 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 912 ), -1, -1, -1 ]
        # CB - CA - N - HN || CT2 - CT1 - NH1 - H
        # CB - CA - C - N || CT2 - CT1 - C - NH1
        # CB - CA - C - O || CT2 - CT1 - C - O
        self.exc.append( qm3.engines.mmres.dihedral( { 1: [ 11.715, 0.0 ], }, [ 912, 910, 923, 924 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 912 ), -1, -1, -1 ]
        # HB1 - CB - CA - N || HA - CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 913, 912, 910, 908 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 913 ), self.sel.index( 912 ), -1, -1 ]
        # HB1 - CB - CA - HA || HA - CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 913, 912, 910, 911 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 913 ), self.sel.index( 912 ), -1, -1 ]
        # HB1 - CB - CA - C || HA - CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 913, 912, 910, 923 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 913 ), self.sel.index( 912 ), -1, -1 ]
        # HB2 - CB - CA - N || HA - CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 914, 912, 910, 908 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 914 ), self.sel.index( 912 ), -1, -1 ]
        # HB2 - CB - CA - HA || HA - CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 914, 912, 910, 911 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 914 ), self.sel.index( 912 ), -1, -1 ]
        # HB2 - CB - CA - C || HA - CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 914, 912, 910, 923 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 914 ), self.sel.index( 912 ), -1, -1 ]
        # CG - CB - CA - N || CPH1 - CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 917, 912, 910, 908 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 917 ), self.sel.index( 912 ), -1, -1 ]
        # CG - CB - CA - HA || CPH1 - CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 917, 912, 910, 911 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 917 ), self.sel.index( 912 ), -1, -1 ]
        # CG - CB - CA - C || CPH1 - CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 917, 912, 910, 923 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 917 ), self.sel.index( 912 ), -1, -1 ]
        # CB - CA - N - C || CT2 - CT1 - NH1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 1: [ 15.062, 0.0 ], }, [ 983, 981, 979, 977 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 983 ), -1, -1, -1 ]
        # CB - CA - N - HN || CT2 - CT1 - NH1 - H
        # CB - CA - C - N || CT2 - CT1 - C - NH1
        # CB - CA - C - O || CT2 - CT1 - C - O
        self.exc.append( qm3.engines.mmres.dihedral( { 1: [ 11.715, 0.0 ], }, [ 983, 981, 994, 995 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 0.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 983 ), -1, -1, -1 ]
        # HB1 - CB - CA - N || HA - CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 984, 983, 981, 979 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 984 ), self.sel.index( 983 ), -1, -1 ]
        # HB1 - CB - CA - HA || HA - CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 984, 983, 981, 982 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 984 ), self.sel.index( 983 ), -1, -1 ]
        # HB1 - CB - CA - C || HA - CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 984, 983, 981, 994 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 984 ), self.sel.index( 983 ), -1, -1 ]
        # HB2 - CB - CA - N || HA - CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 985, 983, 981, 979 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 985 ), self.sel.index( 983 ), -1, -1 ]
        # HB2 - CB - CA - HA || HA - CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 985, 983, 981, 982 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 985 ), self.sel.index( 983 ), -1, -1 ]
        # HB2 - CB - CA - C || HA - CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 985, 983, 981, 994 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 985 ), self.sel.index( 983 ), -1, -1 ]
        # CG - CB - CA - N || CPH1 - CT2 - CT1 - NH1
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 988, 983, 981, 979 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 988 ), self.sel.index( 983 ), -1, -1 ]
        # CG - CB - CA - HA || CPH1 - CT2 - CT1 - HB
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 988, 983, 981, 982 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 988 ), self.sel.index( 983 ), -1, -1 ]
        # CG - CB - CA - C || CPH1 - CT2 - CT1 - C
        self.exc.append( qm3.engines.mmres.dihedral( { 3: [ 1.674, 0.0 ], }, [ 988, 983, 981, 994 ] ) )
        self.exc[-1].ffac = 0.0
        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]
        self.exc[-1].hind = [ self.sel.index( 988 ), self.sel.index( 983 ), -1, -1 ]

        self.size = 3 * len( self.sel )
        self.coor = []
        self.mass = []
        for i in self.sel:
            self.coor += self.mol.coor[3*i:3*i+3]
            self.mass.append( self.mol.mass[i] )
        self.func = 0
        self.grad = []
        self.xumb = []

        for i in sqm:
            self.mol.chrg[i] = 0.0
        self.emm.update_chrg( self.mol )


    def update_coor( self ):
        for i in range( len( self.sel ) ):
            i3 = self.sel[i] * 3
            j3 = i * 3
            self.mol.coor[i3:i3+3] = self.coor[j3:j3+3]


    def get_func( self ):
        self.update_coor()
        self.mol.func = 0
        self.emm.get_func( self.mol )
        self.eqm.get_func( self.mol )
        self.xumb = []
        for umb in self.umb:
            self.xumb.append( umb.get_func( self.mol ) )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0
        self.mol.grad = [ 0 for i in range( 3 * self.mol.natm ) ]
        self.emm.get_grad( self.mol )
        self.eqm.get_grad( self.mol )
        self.fix.get_grad( self.mol )
        self.xumb = []
        for umb in self.umb:
            self.xumb.append( umb.get_grad( self.mol ) )
        for exc in self.exc:
            exc.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = sum( [ self.mol.grad[3*i:3*i+3] for i in self.sel ], [] )



pnt = 0
while( os.path.isdir( "_f%06d"%( pnt ) ) ):
    pnt += 1
os.mkdir( "_f%06d"%( pnt ) )
os.chdir( "_f%06d"%( pnt ) )

sck = socket.socket( socket.AF_INET, socket.SOCK_STREAM ) 
sck.connect( ( "127.0.0.1", 6969 ) )
ox, oy, cx, cy = struct.unpack( "4i", sck.recv( 16 ) )
sck.close()

if( ox >= 0 and oy >= 0 and cx >= 0 and cy >= 0 ):
    obj = my_problem( "../pes.%02d.%02d"%( ox, oy ), cx, cy )
    qm3.actions.minimize.fire( obj, step_number = 1000, gradient_tolerance = 1.0, step_size = 0.1, print_frequency = 1 )
    f = open( "../fixed.pk", "rb" )
    t = pickle.load( f )
    f.close()
    f = open( "../pes.%02d.%02d"%( cx, cy ), "wt" )
    f.write( "REMARK: %12.6lf%12.6lf%18.6lf\n"%( obj.xumb[0], obj.xumb[1], obj.func ) )
    qm3.engines.namd.pdb_write( obj.mol, f, fixed = t )
    f.close()
    obj.emm.stop()
    f = open( "namd.shmid", "rt" )
    os.system( "ipcrm -m " + f.readline().strip() )
    f.close()

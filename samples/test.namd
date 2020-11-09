import qm3.mol
import qm3.fio.xplor
import qm3.engines.namd
import qm3.engines.mmint
import qm3.problem
import qm3.maths.matrix

class my_eng( qm3.problem.template ):

    def __init__( self ):
        self.mol = qm3.mol.molecule( "pdb" )
        qm3.fio.xplor.psf_read( self.mol, "psf" )
        qm3.engines.mmint.non_bonded( self.mol, "nbn" )
        self.mol.chrg[0] = 0.0
        self.mol.chrg[1] = 0.0
        self.mol.chrg[2] = 0.0

        self.eng = qm3.engines.namd.namd()
        self.eng.exe = "./bin/namd2 +ppn 2 namd.inp > namd.out"
        qm3.fio.xplor.psf_write( self.mol, "psf" )

        self.fix = qm3.engines.mmint.QMLJ( self.mol, [ 0, 1, 2 ], [ 3, 4, 5 ], [] )
        
        self.size = 3 * self.mol.natm
        self.coor = self.mol.coor
        self.func = 0.0
        self.grad = []


    def get_func( self ):
        self.mol.func = 0.0
        self.eng.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.size ) ]
        self.eng.get_grad( self.mol )
        self.fix.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = self.mol.grad




f = open( "pdb", "wt" )
f.write( """ATOM      1  OH2 HOH     1       0.121   0.069  -0.225  1.00  0.00      W   
ATOM      2  H1  HOH     1      -0.527   0.165  -0.946  1.00  0.00      W   
ATOM      3  H2  HOH     1       0.702  -0.636  -0.547  1.00  0.00      W   
ATOM      4  OH2 HOH     2      -0.451   1.127   2.211  0.00  0.00      W   
ATOM      5  H1  HOH     2      -0.292   0.595   1.399  0.00  0.00      W   
ATOM      6  H2  HOH     2       0.058   1.927   2.010  0.00  0.00      W   
END
""" )
f.close()

f = open( "par", "wt" )
f.write( """* Force Field Parameter File.
* 

BOND
OT   HT    450.000     0.9572
HT   HT      0.000     1.5139

ANGLE
HT   OT   HT     55.000   104.5200

NONBONDED  NBXMOD 5  GROUP SWITCH CDIEL -
CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0  EPS 1.0  E14FAC 0.83333333  WMIN 1.4
OT      0.00   -0.1521    1.7682
HT      0.00   -0.0460    0.224
""" )
f.close()

f = open( "nbn", "wt" )
f.write( """OT      0.1521    1.7682
HT      0.0460    0.224
""" )
f.close()

f = open( "psf", "wt" )
f.write( """PSF

       1 !NTITLE
 REMARKS original generated structure x-plor psf file

       6 !NATOM
       1 W    1    HOH  OH2  OT    -0.834000       15.9994           0
       2 W    1    HOH  H1   HT     0.417000        1.0080           0
       3 W    1    HOH  H2   HT     0.417000        1.0080           0
       4 W    2    HOH  OH2  OT    -0.834000       15.9994           0
       5 W    2    HOH  H1   HT     0.417000        1.0080           0
       6 W    2    HOH  H2   HT     0.417000        1.0080           0

       6 !NBOND: bonds
       1       2       1       3       2       3       4       5
       4       6       5       6

       2 !NTHETA: angles
       2       1       3       5       4       6

       0 !NPHI: dihedrals


       0 !NIMPHI: impropers


       0 !NDON: donors


       0 !NACC: acceptors


       0 !NNB

       0       0       0       0       0       0

       1       0 !NGRP
       0       0       0

""" )
f.close()

f = open( "namd.inp", "wt" )
f.write( """structure           psf
coordinates         pdb
bincoordinates      namd.coor
paraTypeCharmm      on
parameters          par
fixedatoms          on
fixedatomsfile      pdb
exclude             scaled1-4
1-4scaling          0.5
switching           on
switchdist          4.5
cutoff              6.0
pairlistdist        8.0
wrapAll             off
wrapWater           off
nonbondedFreq       1
fullElectFrequency  1
stepspercycle       1
temperature         0.0
outputname          namd.out
outputEnergies      1
forcedcdfile        namd.force
forcedcdfreq        1
run                 0
output onlyforces   namd
""" )
f.close()


obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )

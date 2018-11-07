import sys
sys.path.append( "/Users/smarti/Devel/ParmEd/lxml/build/lib.macosx-10.13-intel-2.7" )
sys.path.append( "/Users/smarti/Devel/ParmEd/build/lib.macosx-10.13-intel-2.7" )
import parmed

f = open( "pdb", "wt" )
f.write( """ATOM      1  OH2 WAT     1       0.121   0.069  -0.225  1.00  0.00      W   
ATOM      2  H1  WAT     1      -0.527   0.165  -0.946  1.00  0.00      W   
ATOM      3  H2  WAT     1       0.702  -0.636  -0.547  1.00  0.00      W   
ATOM      4  OH2 HOH     2      -0.451   1.127   2.211  0.00  0.00      W   
ATOM      5  H1  HOH     2      -0.292   0.595   1.399  0.00  0.00      W   
ATOM      6  H2  HOH     2       0.058   1.927   2.010  0.00  0.00      W   
END
""" )
f.close()

pdb = parmed.load_file( "pdb" )
pdb.box = [ 100, 100, 100, 90., 90., 90. ]
parmed.write_PDB( pdb, "pdb" )

f = open( "top.rtf", "wt" )
f.write( """* Topology File.
* 
   22   1
MASS     1 OT     15.99940  O
MASS     2 HT      1.00800  H

DEFA FIRS NONE LAST NONE   

RESI HOH 0.000
GROUP   
ATOM OH2  OT     -0.834
ATOM H1   HT      0.417
ATOM H2   HT      0.417
BOND OH2 H1 OH2 H2 H1 H2 
ANGLE H1 OH2 H2
ACCEPTOR OH2   

RESI WAT 0.000
GROUP   
ATOM OH2  OT     -0.834
ATOM H1   HT      0.417
ATOM H2   HT      0.417

END
""" )
f.close()

f = open( "par.prm", "wt" )
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

prm = parmed.charmm.CharmmParameterSet( "top.rtf", "par.prm" )
print( prm.residues )
xml = parmed.openmm.OpenMMParameterSet.from_parameterset( prm )
xml.write( "prm.xml" )

f = open( "psf", "wt" )
f.write( """PSF

       1 !NTITLE
 REMARKS original generated structure x-plor psf file

       6 !NATOM
       1 W    1    WAT  OH2  OT    -0.834000       15.9994           0
       2 W    1    WAT  H1   HT     0.417000        1.0080           0
       3 W    1    WAT  H2   HT     0.417000        1.0080           0
       4 W    2    HOH  OH2  OT    -0.834000       15.9994           0
       5 W    2    HOH  H1   HT     0.417000        1.0080           0
       6 W    2    HOH  H2   HT     0.417000        1.0080           0

       2 !NBOND: bonds
       4       5       4       6

       1 !NTHETA: angles
       5       4       6

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

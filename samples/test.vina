import qm3.mol
import qm3.actions.vina
try:
    import cStringIO as io
except:
    import io

f = io.StringIO( """
ATOM      1  OH2 HOH     1       0.121   0.069  -0.225  0.00  0.00      W   
ATOM      2  H1  HOH     1      -0.527   0.165  -0.946  0.00  0.00      W   
ATOM      3  H2  HOH     1       0.702  -0.636  -0.547  0.00  0.00      W   
ATOM      4  OH2 HOH     2      -0.451   1.127   2.211  0.00  0.00      W   
ATOM      5  H1  HOH     2      -0.292   0.595   1.399  0.00  0.00      W   
ATOM      6  H2  HOH     2       0.058   1.927   2.010  0.00  0.00      W   
END
""" )
m = qm3.mol.molecule( f )
f.close()

f = io.StringIO( """PSF

       1 !NTITLE
 REMARKS original generated structure x-plor psf file

       6 !NATOM
       1 W    1    HOH  OH2  OT    -0.834000       15.9994           0
       2 W    1    HOH  H1   HT     0.417000        1.0080           0
       3 W    1    HOH  H2   HT     0.417000        1.0080           0
       4 W    2    HOH  OH2  OT    -0.834000       15.9994           0
       5 W    2    HOH  H1   HT     0.417000        1.0080           0
       6 W    2    HOH  H2   HT     0.417000        1.0080           0

       4 !NBOND: bonds
       1       2       1       3       4       5       4       6
""" )
qm3.actions.vina.psf_read( m, f )
f.close()

qm3.actions.vina.receptor( m, fname = "x.rec", sele = [ 3, 4, 5 ] )
qm3.actions.vina.ligand( m, fname = "x.lig", sele = [ 0, 1, 2 ] )
cen = [ 0.0, 0.0, 0.0 ]
for i in [ 0, 1, 2 ]:
    for j in [0, 1, 2]:
        cen[j] += m.coor[3*i+j]
cen = [ i / 3.0 for i in cen ]
f = open( "vina.conf", "wt" )
f.write( """
receptor  = x.rec
ligand    = x.lig
center_x  =  0.099
center_y  = -0.134
center_z  = -0.573
size_x    = 30.
size_y    = 30.
size_z    = 30.
num_modes = 10
energy_range = 5
""" )
f.close()

import os
os.system( "./bin/vina --cpu 2 --config vina.conf" )

qm3.actions.vina.parse_dock( m, fname = "x.lig_out.pdbqt", sele = [ 0, 1, 2 ], model = 10 )
m.pdb_write()

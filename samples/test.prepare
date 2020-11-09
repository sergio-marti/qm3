import qm3.mol
import qm3.elements
import qm3.utils.prepare
import qm3.constants
import qm3.utils
import os
import sys

try:
    import cStringIO as io
except:
    import io

f = io.StringIO( """ATOM      1  N   NH3     1       0.040   0.029  -0.171  0.00  0.00      A   
ATOM      2  H1  NH3     1       0.917  -0.069  -0.589  0.00  0.00      A   
ATOM      3  H2  NH3     1      -0.443   0.874  -0.242  0.00  0.00      A   
ATOM      4  H3  NH3     1      -0.354  -0.719   0.318  0.00  0.00      A   
END
""" )
o = qm3.mol.molecule( f )
f.close()
o.chrg = [ -0.729296, 0.243101, 0.243119, 0.243076 ]
t = qm3.utils.prepare.counter_ions( o, num = 1, d_prt = 6., d_ion = 6. )
for i in range( t.natm ):
    t.labl[i] = "SOD"
    t.resn[i] = "SOD"
o.append( t )
o.anum = [ 7, 1, 1, 1, 11 ]
o.mass = [ qm3.elements.mass[i] for i in o.anum ]
qm3.utils.moments_of_inertia( o.mass, o.coor )

d = qm3.constants.water_density()
m = qm3.elements.calc_mass( "H2O1" )
v = 20
n = int( round( d * qm3.constants.NA * ( v * v * v * 1.e-24 ) / m, 0 ) )
f = open( "wat", "wt" )
f.write( """
HETATM    1  OH2 HOH     1      10.203   7.604  12.673
HETATM    2  H1  HOH     1       9.626   6.787  12.673
HETATM    3  H2  HOH     1       9.626   8.420  12.673
CONECT    1    2
CONECT    1    3
CONECT    1    2    3
END
""" )
f.close()
f = open( "inp", "wt" )
f.write( """tolerance 2.0
output pdb
filetype pdb
structure wat
  number %d
  inside cube 0.5 0.5 0.5 %.1lf
end structure
"""%( n, v - 0.5 ) )
f.close()
os.system( "./bin/packmol < inp" )

s = qm3.mol.molecule( "pdb" )
f.close()
s.anum = [ qm3.elements.rsymbol[s.labl[i][0]] for i in range( s.natm ) ]
s.mass = [ qm3.elements.mass[s.anum[i]] for i in range( s.natm ) ]
o.append( qm3.utils.prepare.solvate( o, s ) )
o.pdb_write( "pdb" )

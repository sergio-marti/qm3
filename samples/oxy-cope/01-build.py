import  os
import  re
import  qm3.constants
import  qm3.elements
import  qm3.utils.prepare
import  qm3.mol
import  pickle
import  qm3.engines.namd


f = open( "01.acs", "wt" )
f.write( """ATOM      1  C1  COP     1      -1.039   1.239  -0.102  1.00  0.00      A
ATOM      2  C2  COP     1       0.454   1.279   0.023  0.00  0.00      A
ATOM      3  C3  COP     1       1.031  -0.118   0.450  0.00  0.00      A
ATOM      4  C4  COP     1       0.626  -1.160  -0.570  0.00  0.00      A
ATOM      5  C5  COP     1      -0.543  -1.803  -0.531  0.00  0.00      A
ATOM      6  C6  COP     1      -1.688   0.831  -1.194  0.00  0.00      A
ATOM      7  O7  COP     1       2.448  -0.036   0.595  0.00  0.00      A
ATOM      8  H8  COP     1       2.840   0.139  -0.277  0.00  0.00      A
ATOM      9  H9  COP     1      -1.613   1.486   0.793  0.00  0.00      A
ATOM     10  H10 COP     1       0.762   2.019   0.768  0.00  0.00      A
ATOM     11  H11 COP     1       0.917   1.565  -0.927  0.00  0.00      A
ATOM     12  H12 COP     1       0.620  -0.382   1.443  0.00  0.00      A
ATOM     13  H13 COP     1       1.279  -1.261  -1.441  0.00  0.00      A
ATOM     14  H14 COP     1      -0.871  -2.457  -1.331  0.00  0.00      A
ATOM     15  H15 COP     1      -1.224  -1.681   0.304  0.00  0.00      A
ATOM     16  H16 COP     1      -1.169   0.554  -2.106  0.00  0.00      A
ATOM     17  H17 COP     1      -2.766   0.759  -1.193  0.00  0.00      A
""" )
f.close()


d = qm3.constants.water_density()
m = qm3.elements.calc_mass( "H2O1" )
v = 40
n = int( round( d * qm3.constants.NA * ( v * v * v * 1.e-24 ) / m, 0 ) )
f = open( "01.hoh", "wt" )
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
f = open( "01.inp", "wt" )
f.write( """tolerance 2.0
output 01.wat
filetype pdb
structure 01.hoh
  number %d
  inside cube 0.5 0.5 0.5 %.1lf
end structure
"""%( n, v - 0.5 ) )
f.close()
os.system( "./bin/packmol < 01.inp" )


m = qm3.mol.molecule( "01.acs" )
m.anum = [ qm3.elements.rsymbol[m.labl[i][0]] for i in range( m.natm ) ]
m.fill_masses()
s = qm3.mol.molecule( "01.wat" )
s.anum = [ qm3.elements.rsymbol[s.labl[i][0]] for i in range( s.natm ) ]
s.fill_masses()
m.append( qm3.utils.prepare.solvate( m, s ) )
for i in range( m.res_lim[1], m.natm ):
    m.segn[i] = "W"
m.settle()
m.norm_resid()
m.pdb_write( "01.pdb" )

s = list( sorted( m.indx["A"][1].values() ) )
f = open( "01.sele_QM.pk", "wb" )
pickle.dump( s, f )
f.close()

f = open( "01.sele_MM.pk", "wb" )
pickle.dump( list( set( m.sph_sel( s, 18.0 ) ).difference( set( s ) ) ), f )
f.close()

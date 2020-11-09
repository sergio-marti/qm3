import os
import qm3.mol
import qm3.engines.gromacs
import qm3.problem
import qm3.maths.matrix
try:
    import cStringIO as io
except:
    import io

class my_eng( qm3.problem.template ):

    def __init__( self ):
        self.mol = qm3.engines.gromacs.coordinates_read( "gromacs.g96" )

        self.eng = qm3.engines.gromacs.gromacs()
        self.eng.exe  = "source ./bin/gmx_rc; "
        self.eng.exe += "gmx mdrun -s gromacs.tpr -e gromacs.edr -o gromacs.trr -g gromacs.log -rerun gromacs.g96 >& gromacs.out2"

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
        self.func = self.mol.func
        self.grad = self.mol.grad



os.system( "rm -f gromacs.*" )

f = io.StringIO( """ATOM      1  OH2 QWT     1       0.121   0.069  -0.225  1.00  0.00      W   
ATOM      2  H1  QWT     1      -0.527   0.165  -0.946  1.00  0.00      W   
ATOM      3  H2  QWT     1       0.702  -0.636  -0.547  1.00  0.00      W   
ATOM      4  OH2 MWT     2      -0.451   1.127   2.211  0.00  0.00      W   
ATOM      5  H1  MWT     2      -0.292   0.595   1.399  0.00  0.00      W   
ATOM      6  H2  MWT     2       0.058   1.927   2.010  0.00  0.00      W   
END
""" )
f.seek( 0 )
qm3.engines.gromacs.coordinates_write( qm3.mol.molecule( f ), "gromacs.g96" )


f = open( "gromacs.itp", "wt" )
f.write( """
[ moleculetype ]
MWT    2

[ atoms ]
1     opls_111  1       MWT             OH2             1       -0.834
2     opls_112  1       MWT              H1             1        0.417
3     opls_112  1       MWT              H2             1        0.417

[ bonds ]
1   2   1   0.09572   376560.0
1   3   1   0.09572   376560.0


[ angles ]
2   1   3   1   104.52   460.24


[ moleculetype ]
QWT    1

[ atoms ]
1     opls_111  1       QWT             OH2             1        0.000
2     opls_112  1       QWT              H1             1        0.000
3     opls_112  1       QWT              H2             1        0.000

[ bonds ]
1   2   1   0.1   0.00000001
1   3   1   0.1   0.00000001
2   3   1   0.1   0.00000001
""" )
f.close()

f = open( "gromacs.top", "wt" )
f.write( """#include "oplsaa.ff/forcefield.itp"
#include "gromacs.itp"

[ system ]
QM/MM water dimer

[ molecules ]
QWT 1
MWT 1
""" )
f.close()

f = open( "gromacs.mdp", "wt" )
f.write( """
define = -DFLEXIBLE
integrator = md-vv
dt = 0.001
nsteps = 1
nstfout = 1
nstenergy = 1

cutoff-scheme = group
pbc = no
ns-type = simple
rlist = 0.8
coulombtype = reaction-field-zero
rcoulomb-switch = 0.45
rcoulomb = 0.6
vdwtype = switch
rvdw-switch = 0.45
rvdw = 0.6
""" )
f.close()

cmd  = "source ./bin/gmx_rc; "
cmd += "gmx grompp -f gromacs.mdp -c gromacs.g96 -p gromacs.top -o gromacs.tpr -maxwarn 100 >& gromacs.out1"
os.system( cmd )


obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )

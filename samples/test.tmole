import qm3.mol
import qm3.problem
import qm3.engines.tmole
import qm3.maths.matrix
import os

class my_eng( qm3.problem.template ):

    def __init__( self ):
        self.mol = qm3.mol.molecule()
        self.mol.xyz_read( "xyz" )
        self.mol.chrg = [ -0.834, 0.417, 0.417, -0.834, 0.417, 0.417 ]
        self.mol.anum = [ 8, 1, 1, 8, 1, 1 ]

        self.eng = qm3.engines.tmole.tmole( self.mol, [ 0, 1, 2 ], [ 3, 4, 5 ] )

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
    


f = open( "tmole.rc", "wt" )
f.write( """
export TURBODIR=/usr/local/chem/tmole_7.1
export PATH=$TURBODIR/bin:$PATH

# PARALLEL ------------------------------------------ (*_omp)
export LD_LIBRARY_PATH=$TURBODIR/bin:$LD_LIBRARY_PATH
unset OMP_DYNAMIC
export OMP_NUM_THREADS=4
""" )
f.close()


f = open( "xyz", "wt" )
f.write( """3

O       6.778   2.205  -2.327
H       7.469   2.103  -1.644
H       7.332   2.143  -3.133
""" )
f.close()
f = open( "setup", "wt" )
f.write( """
source tmole.rc

x2t xyz > coord
define << EOD

slave
a coord
*
no
b all def-SVP
*
eht
y
0
y
dft
func b3-lyp
on
*
*
EOD
""" )
f.close()
os.system( "rm -vf basis control coord mos; bash setup" )


f = open( "control", "rt" )
buf = f.readlines()
f.close()
f = open( "control", "wt" )
for l in buf:
    if( l.strip() == "$drvopt" ):
        f.write( l )
        f.write( "   point charges\n" )
    elif( l.strip() == "$end" ):
        f.write( "$point_charges file=charges\n" )
        f.write( "$point_charge_gradients file=charges.gradient\n" )
        f.write( l )
    else:
        f.write( l )
f.close()


f = open( "xyz", "wt" )
f.write( """6

O       0.12109      0.06944     -0.22458
H      -0.52694      0.16499     -0.94583
H       0.70159     -0.63565     -0.54677
O      -0.45114      1.12675      2.21102
H      -0.29157      0.59483      1.39876
H       0.05804      1.92714      2.01036
""" )
f.close()


obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )

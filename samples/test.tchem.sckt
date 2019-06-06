import qm3.mol
import qm3.problem
import qm3.engines.tchem
import qm3.maths.matrix
import os
import time


class my_eng( qm3.problem.template ):

    def __init__( self ):
        self.mol = qm3.mol.molecule()
        self.mol.xyz_read( "xyz" )
        self.mol.chrg = [ -0.834, 0.417, 0.417, -0.834, 0.417, 0.417 ]
        self.mol.anum = [ 8, 1, 1, 8, 1, 1 ]
        self.eng = qm3.engines.tchem.tchem_sckt( "sckt_", self.mol, [ 0, 1, 2 ], [ 3, 4, 5 ] )
        self.eng.exe = "bash r.tchem"

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


f = open( "r.tchem", "wt" )
f.write( """
source /usr/local/chem/mpich2-1.4.1p1/rc

export TeraChem=/usr/local/chem/tchem_1.93p
export NBOEXE=$TeraChem/bin/nbo6.i4.exe
export OMP_NUM_THREADS=16
export LD_LIBRARY_PATH=$TeraChem/lib:$LD_LIBRARY_PATH
export PATH=$TeraChem/bin:$PATH

rm -rf tchem_; mkdir tchem_; cd tchem_; mkdir scr
mpirun -n 1 $TeraChem/bin/terachem -U1 -Mtcport_ >& tchem.log
""" )
f.close()
os.system( "bash r.tchem &" )
time.sleep( 1 )


f = open( "input", "wt" )
f.write( """basis         def2-svp
charge        0
spinmult      1
method        b3lyp
dftd          no
gpus          1 0
dftgrid       1
threall       1.e-11
convthre      3.e-5
guess         scr/c0
run           gradient
coordinates   tchem_qm.xyz
pointcharges  tchem_mm.xyz
amber         yes
""" )
f.close()


f = open( "r.sckt", "wt" )
f.write( """
source /usr/local/chem/mpich2-1.4.1p1/rc

rm -f sckt_ sckt_.log
mpirun -n 1 ./a.out tcport_ sckt_ input >& sckt_.log
""" )
f.close()
os.system( "bash r.sckt &" )
time.sleep( 1 )


obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )

import qm3.mol
import qm3.problem
import qm3.engines.dftd4
import qm3.maths.matrix
try:
  import cStringIO as io
except:
  import io
import os

os.environ["QM3_LIBDFTD4"] = "./bin/dftd4.so"


class my_eng( qm3.problem.template ):

    def __init__( self ):
        self.mol = qm3.mol.molecule()
        f = io.StringIO( """6

O       0.12109      0.06944     -0.22458
H      -0.52694      0.16499     -0.94583
H       0.70159     -0.63565     -0.54677
O      -0.45114      1.12675      2.21102
H      -0.29157      0.59483      1.39876
H       0.05804      1.92714      2.01036
""" )
        self.mol.xyz_read( f )
        f.close()
        self.mol.chrg = [ -0.834, 0.417, 0.417, -0.834, 0.417, 0.417 ]
        self.mol.anum = [ 8, 1, 1, 8, 1, 1 ]
#        self.eng = qm3.engines.dftd4.dftd4( self.mol, [ 0, 1, 2 ] )
#        self.eng.exe = "bash r.dftd4"
        prm = { "chrg": 0.0, "s6": 1.00, "s8": 2.02929367, "a1": 0.40868035, "a2": 4.53807137  }
        self.eng = qm3.engines.dftd4.dl_dftd4( self.mol, prm, [ 0, 1, 2 ] )

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
    


f = open( "r.dftd4", "wt" )
f.write( "./bin/dftd4 dftd4.xyz --chrg 0 --func b3-lyp --grad >& dftd4.log" )
f.close()


obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3, fmt = "%20.10lf" )

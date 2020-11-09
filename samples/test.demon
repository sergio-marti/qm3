import qm3.mol
import qm3.problem
import qm3.engines.demon
import qm3.maths.matrix
try:
  import cStringIO as io
except:
  import io


class my_eng( qm3.problem.template ):

    def __init__( self ):
        self.mol = qm3.mol.molecule()
        self.mol.xyz_read( "xyz" )
        self.mol.chrg = [ -0.834, 0.417, 0.417, -0.834, 0.417, 0.417 ]
        self.mol.anum = [ 8, 1, 1, 8, 1, 1 ]
        f = io.StringIO( """
title _slave_
symmetry off
visualization off
charge 0
multiplicity 1
basis (def2-svp)
scftype rks tol=1.0e-6
vxctype b3lyp
eris mixed
qm3_guess
embed file
qm/mm charmm
geometry cartesian angstrom
qm3_atoms
""" )
        f.seek( 0 )
        self.eng = qm3.engines.demon.demon( self.mol, f, [ 0, 1, 2 ], [ 3, 4, 5 ] )
        self.eng.exe = "./bin/deMon_4.4.1"

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


# https://bse.pnl.gov/bse/portal
f = open( "BASIS", "wt" )
f.write( """O-HYDROGEN H (Def2-SVP)
    3
    1    0    3
     13.0107010              0.19682158E-01   
      1.9622572              0.13796524       
      0.44453796             0.47831935       
    2    0    1
      0.12194962             1.0000000        
    2    1    1
      0.8000000              1.0000000        
O-OXYGEN O (Def2-SVP)
    6
    1    0    5
   2266.1767785             -0.53431809926E-02      
    340.87010191            -0.39890039230E-01      
     77.363135167           -0.17853911985    
     21.479644940           -0.46427684959    
      6.6589433124          -0.44309745172    
    2    0    1
      0.80975975668          1.0000000        
    3    0    1
      0.25530772234          1.0000000        
    2    1    3
     17.721504317            0.43394573193E-01      
      3.8635505440           0.23094120765    
      1.0480920883           0.51375311064    
    3    1    1
      0.27641544411          1.0000000        
    3    2    1
      1.2000000              1.0000000        
""" )
f.close()

obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )

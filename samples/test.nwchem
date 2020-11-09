import qm3.mol
import qm3.problem
import qm3.engines.nwchem
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
start nwchem
geometry units angstroms nocenter noautoz noautosym
qm3_atoms
end
basis
 * library def2-svp
end
charge 0
dft
 qm3_guess
 mulliken
 mult 1
 iterations 100
 direct
 xc b3lyp 
end
qm3_charges
qm3_job
""" )
        f.seek( 0 )
        self.eng = qm3.engines.nwchem.nwchem( self.mol, f, [ 0, 1, 2 ], [ 3, 4, 5 ] )
        self.eng.exe = "NWCHEM_BASIS_LIBRARY=./bin/nwchemlib/"
        self.eng.exe += " mpirun -n 2 ./bin/nwchem nwchem.nw > nwchem.log"

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


obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )

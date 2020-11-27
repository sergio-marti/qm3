import os
import sys
import qm3.mol
import qm3.problem
import qm3.engines.lio
import qm3.maths.matrix
try:
  import cStringIO as io
except:
  import io


os.environ["LIOHOME"] = "./bin/slio"
os.environ["OMP_NUM_THREADS"] = "2"
os.environ["QM3_LIBLIO"] = "./bin/lio.so"


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
        self.mol.chrg = [ -0.834, 0.417, 0.417, -0.834, 0.417, 0.417 ]
        self.mol.anum = [ 8, 1, 1, 8, 1, 1 ]
        if( len( sys.argv ) == 3 ):
            f = io.StringIO( """&lio
 natom       = qm3_natm
 nsol        = qm3_nchg
 charge      = 0
 writeforces = t
 mulliken    = t
 nmax        = 100
 verbose     = 4
 basis_set   = "6-31G**"
 told        = 1.0d-9
 use_libxc   = .true.
ex_functional_id = %s
ec_functional_id = %s
&end
"""%( sys.argv[1], sys.argv[2] ) )
            os.system( "xc-info %s 2>&1"%( sys.argv[1] ) )
            os.system( "xc-info %s 2>&1"%( sys.argv[2] ) )
        else:
            f = io.StringIO( """&lio
 natom       = qm3_natm
 nsol        = qm3_nchg
 charge      = 0
 writeforces = t
 mulliken    = t
 nmax        = 100
 verbose     = 4
 basis_set   = "6-31G**"
 told        = 1.0d-9
&end
""" )
#        self.eng = qm3.engines.lio.lio( self.mol, f, [ 0, 1, 2 ], [ 3, 4, 5 ] )
#        self.eng.exe = "bash r.lio"
        self.eng = qm3.engines.lio.dl_lio( self.mol, f, [ 0, 1, 2 ], [ 3, 4, 5 ] )

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
    

f = open( "r.lio", "wt" )
f.write( """
export LIOHOME=./bin/slio
export PATH=$LIOHOME:$PATH
export LD_LIBRARY_PATH=$LIOHOME:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=2
lio.x -i lio.inp -c lio.xyz > lio.log
""" )
f.close()


print( 80 * "=" )
obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )

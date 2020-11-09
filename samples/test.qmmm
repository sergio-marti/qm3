import    qm3.mol
import    qm3.problem
import    qm3.engines.mol_mech
import    qm3.engines.sqm
import    qm3.actions.minimize
import    pickle
import    os
try:
    import    cStringIO as io
except:
    import    io


#qm3.engines.mol_mech.mol_mech_so = False
os.environ["QM3_LIBSQM"] = "./bin/sqm.so"


class my_problem( qm3.problem.template ):
    def __init__( self ):
        qm3.problem.template.__init__( self )

        self.mol = qm3.mol.molecule()
        f = io.StringIO( """22

N    0.07218  -0.00593   0.03395
H   -0.34526  -0.46245   0.82019
C    1.50632   0.07581  -0.03848
H    1.91105  -0.30487   0.88829
C    2.04479  -0.73937  -1.24235
H    1.65228  -1.77726  -1.14982
H    1.62683  -0.31520  -2.18212
C    3.58048  -0.79949  -1.33174
H    3.98508   0.23579  -1.41135
H    3.97548  -1.23992  -0.38766
C    4.06792  -1.62441  -2.53577
H    3.65199  -2.65336  -2.44080
H    3.64058  -1.16493  -3.45645
C    5.59652  -1.68158  -2.63884
H    6.01772  -0.65489  -2.71249
H    6.02714  -2.18646  -1.74652
N    6.00363  -2.43769  -3.84584
H    7.04153  -2.46918  -3.91356
H    5.63318  -3.40807  -3.79236
H    5.61916  -1.97186  -4.69280
C    1.90740   1.52565  -0.15031
O    2.90065   1.95032   0.44298
""" )
        self.mol.xyz_read( f )
        f.close()
        self.mol.guess_atomic_numbers()
        self.mol.type = []
        self.mol.chrg = []
        f = io.StringIO( """N.3     -0.4700
H        0.3100
C.3      0.0700
Hn       0.0900
C.3     -0.1800
Hn       0.0900
Hn       0.0900
C.3     -0.1800
Hn       0.0900
Hn       0.0900
C.3     -0.1800
Hn       0.0900
Hn       0.0900
C.3      0.2100
Hn       0.0500
Hn       0.0500
N.3     -0.3000
H        0.3300
H        0.3300
H        0.3300
C.co     0.5100
O.2     -0.5100
       5       3       8       5      11       8      14      11
      17      14       1       2       1       3      21       3
       3       4       5       6       5       7       8       9
       8      10      11      12      11      13      14      15
      14      16      22      21      17      18      17      19
      17      20
""" )
        for i in range( self.mol.natm ):
            t = f.readline().strip().split()
            self.mol.type.append( t[0] )
            self.mol.chrg.append( float( t[1] ) )
        bnd = []
        tmp = [ int( i )-1 for i in f.read().split() ]
        for i in range( len( tmp ) // 2 ):
            bnd.append( [ tmp[2*i], tmp[2*i+1] ] )
        f.close()

        self.eMM = qm3.engines.mol_mech.simple_force_field( self.mol )
        self.eMM.topology( self.mol, bond = bnd, qtyp = False, qchg = False )
        self.eMM.parameters( self.mol )
        self.eMM.update_non_bonded( self.mol )

        sqm = [ 11-1, 12-1, 13-1, 14-1, 15-1, 16-1, 17-1, 18-1, 19-1, 20-1 ]
        smm = list( set( list( range( self.mol.natm ) ) ).difference( set( sqm ) ) )

        f = io.StringIO( """
_slave_
&qmmm 
maxcyc    = 0,
qm_theory = "AM1",
qmcharge  = +1,
qmmm_int  = 1,
verbosity = 4
 /
qm3_atoms
qm3_charges
""" )
        self.eQM = qm3.engines.sqm.dl_sqm( self.mol, f, sqm, smm, [ [ 11-1, 8-1 ] ] )
        f.close()

        self.eMM.qm_atoms( sqm )

        self.sel = list( range( self.mol.natm ) )

        self.size = 3 * len( self.sel )
        self.coor = []
        self.mass = []
        for i in self.sel:
            self.coor += self.mol.coor[3*i:3*i+3]
            self.mass.append( self.mol.mass[i] )

        self.fd = open( "xyz", "wt" )


    def current_step( self, step ):
        self.mol.xyz_write( self.fd )
        self.fd.flush()


    def update_coor( self ):
        for i in range( len( self.sel ) ):
            for j in [0, 1, 2]:
                self.mol.coor[3*self.sel[i]+j] = self.coor[3*i+j]


    def get_func( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.eQM.get_func( self.mol )
        self.eMM.get_func( self.mol, qprint = True )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( 3 * self.mol.natm ) ]
        self.eQM.get_grad( self.mol )
        self.eMM.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = []
        for i in self.sel:
            self.grad += self.mol.grad[3*i:3*i+3]
        qm3.utils.project_RT_modes( self.mass, self.coor, self.grad, [] )



obj = my_problem()
qm3.actions.minimize.l_bfgs( obj, print_frequency = 1, step_number = 1000, gradient_tolerance = 1.0, step_size = 0.1 )
obj.fd.close()
obj.get_grad()
print( obj.func )
import qm3.maths.matrix
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )

import qm3.mol
import qm3.engines.mol_mech
import qm3.maths.matrix

try:
    import cStringIO as io
except:
    import io

f = io.StringIO( """21

C   -0.01476   0.12442  -0.20255
C    0.82769   0.56556   0.74695
C    1.66489  -0.32106   1.63054
C    2.79301   0.40129   2.43845
C    2.38361  -1.23996   0.68764
N    0.81099  -1.18271   2.45228
H    3.38479  -0.35878   2.95699
C    2.22092   1.37970   3.53401
O    3.72665   1.08583   1.54760
H   -0.14531  -0.93232  -0.39818
H   -0.59221   0.82684  -0.78640
H    0.92057   1.63685   0.89636
H    3.93090   0.46971   0.82467
O    1.50556   2.34387   3.14131
O    2.51366   1.07455   4.72941
H    0.36554  -0.63370   3.15950
H    1.36912  -1.89187   2.88308
H    0.11546  -1.61071   1.87517
N    3.28838  -0.43777  -0.19999
H    2.26545  -2.33303   0.65253
H    3.33456   0.55611  -0.09966""" )

m = qm3.mol.molecule()
m.xyz_read( f )
f.close()
m.guess_atomic_numbers()
m.chrg = [ .0 for i in range( m.natm ) ]
m.chrg[5]  = +1.0
m.chrg[13] = -0.5
m.chrg[14] = -0.5
print( m.boxl )

#qm3.engines.mol_mech.mol_mech_so = False
print( "Binary version available:", qm3.engines.mol_mech.mol_mech_so )

o = qm3.engines.mol_mech.simple_force_field( m )
o.topology( m, impr = [ [ 8-1, 14-1, 15-1, 4-1, 300.0, 0.0 ] ] )
o.parameters( m )
o.cut_on   = -1
o.cut_off  = -1
o.cut_list = -1

print( 80*"=" )
print( "%-4s%-4s%8s%10s "%( "", "Atm", "Chrg", "Typ" ), "Con" )
for i in range( m.natm ):
    print( "%-4d%-4s%8.3lf%10s "%( i, m.labl[i], m.chrg[i], m.type[i] ), o.conn[i] )
print( 80*"-" )
print( "Charge: %8.3lf"%( sum( m.chrg ) ) )

m.func = 0.0
m.grad = [ 0.0 for i in range( 3 * m.natm ) ]
o.get_grad( m, qprint = True )
qm3.maths.matrix.mprint( m.grad, m.natm, 3 )

#################################################################################
import qm3.actions.minimize
import qm3.problem
class my_obj( qm3.problem.template ):
    def __init__( self, mol, eng ):
        self.mol = mol
        self.eng = eng
        self.size = 3 * self.mol.natm
        self.coor = self.mol.coor
        self.func = 0.0
        self.grad = []

    def get_func( self ):
        self.mol.func = 0.0
        self.eng.get_func( self.mol, True )
        self.func = self.mol.func

    def get_grad( self ):
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.size ) ]
        self.eng.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = self.mol.grad
    
x = my_obj( m, o )
qm3.actions.minimize.fire( x, step_number = 1000, print_frequency = 100, gradient_tolerance = 1.0 )
x.mol.xyz_write( "xyz" )
x.get_func()

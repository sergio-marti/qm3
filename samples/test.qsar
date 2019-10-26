import qm3.mol
import qm3.utils.qsar

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
m.chrg = [ -0.371, -0.116, 0.182, -0.024, -0.051, -0.878, 0.247, 0.634, -0.649, 0.255, 0.262, 0.269, 0.402, -0.570, -0.573, 0.360, 0.336, 0.335, -0.629, 0.221, 0.358 ]

m = qm3.utils.qsar.remove_non_polar_H( m )
o = qm3.utils.qsar.topological_index( m )

print( 80*"=" )
print( "%-4s%-4s%8s "%( "", "Atm", "Chrg" ), "Con" )
for i in range( m.natm ):
    print( "%-4d%-4s%8.3lf "%( i, m.labl[i], m.chrg[i] ), o.conn[i] )
print( 80*"-" )
print( "Charge: %8.3lf"%( sum( m.chrg ) ) )

print( 80*"=" )
print( "Top_index[Randic]:", o.Randic() )
print( "Top_index[Hosoya]:", o.Hosoya() )
if( len( o.bond ) <= 20 ):
    print( "Top_index[Wiener]:", o.Wiener() )
print( "Top_index[QNA(1,0)]:", o.QNA( 1, 0 ) )
print( "Top_index[QNA(0,1)]:", o.QNA( 0, 1 ) )
print( "Top_index[QNA(1,1)]:", o.QNA( 1, 1 ) )

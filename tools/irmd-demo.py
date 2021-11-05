#!/usr/bin/env python3
import  sys
import  math
import  io
import  os
import  struct
import  numpy as np
import  matplotlib.pyplot as plt
import  qm3.problem
import  qm3.mol
import  qm3.constants
import  qm3.utils
import  qm3.fio.dcd
import  qm3.engines.xtb
import  qm3.actions.dynamics
import  qm3.actions.minimize


class my_problem( qm3.problem.template ):
    def __init__( self, flag = False ):
        qm3.problem.template.__init__( self )
        self.mol = qm3.mol.molecule()
        if( not os.path.isfile( "optimized" ) ):
            f = io.StringIO( """30

O       1.15100      2.40460     -0.00590
O      -1.30410     -2.38240     -0.02630
N       1.51530     -1.09090      0.00191
N      -1.43260      1.09110     -0.01080
C       0.73300      0.03270     -0.00121
C       2.90680      0.65281      0.00070
C       2.84700     -0.74240      0.00320
C      -0.71860     -0.00400     -0.00440
C       1.56550      1.12820     -0.00191
C      -2.97420     -0.67090      0.00430
C      -1.63270     -1.22380     -0.00351
C      -2.78790      0.70710     -0.00340
C       4.16970      1.27760      0.00130
C       3.99180     -1.54730      0.00610
C      -4.21060     -1.27990      0.01170
C      -3.89020      1.54340     -0.00191
C       5.32340      0.48510      0.00400
C       5.23350     -0.90530      0.00660
C      -5.32200     -0.43760      0.01320
C      -5.16420      0.96170      0.00630
H       1.16990     -2.04090      0.00290
H       4.25710      2.36030     -0.00050
H       3.92280     -2.62990      0.00801
H      -4.32160     -2.35730      0.01600
H      -3.77840      2.62240     -0.00720
H       6.30010      0.96150      0.00440
H       6.14250     -1.50140      0.00880
H      -6.32200     -0.86310      0.01910
H      -6.04510      1.59850      0.00740
H       1.91420      3.00820     -0.00610
""" )
        else:
            f = open( "optimized", "rt" )
        self.mol.xyz_read( f )
        self.mol.guess_atomic_numbers()
        self.mol.fill_masses()
        self.eng = qm3.engines.xtb.dl_xtb( self.mol, 0, 0, list( range( self.mol.natm ) ) )

        self.size = 3 * self.mol.natm
        self.coor = self.mol.coor
        self.mass = self.mol.mass

        self.fvel = None
        if( flag ):
            self.fvel = open( "velo.bin", "wb" )
            self.fdcd = qm3.fio.dcd.dcd()
            self.fdcd.open_write( "dcd", self.mol.natm )


    def current_step( self, istep ):
        if( self.fvel != None ):
            for itm in self.velo:
                self.fvel.write( struct.pack( "d", itm ) )
            self.fvel.flush()
            self.fdcd.append( self.mol )


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



stp = 20000
dlt = 0.2 # fs

if( not os.path.isfile( "optimized" ) ):
    obj = my_problem()
    qm3.actions.minimize.fire( obj, step_number = 1000 )
    obj.mol.xyz_write( "optimized" )
    sys.exit( 1 )

if( not os.path.isfile( "velo.bin" ) ):
    obj = my_problem( True )
    qm3.actions.dynamics.assign_velocities( obj, 300.0 )
    qm3.actions.dynamics.velocity_verlet( obj, temperature = 300.0, scale_frequency = 10, step_number = stp, step_size = dlt / 1000.0, print_frequency = 1 )
    obj.fvel.close()
    obj.fdcd.close()
    sys.exit( 1 )

obj = my_problem()
velo = np.zeros( ( obj.size, stp ) )
f = open( "velo.bin", "rb" )
for i in range( stp ):
    for j in range( obj.size ):
        velo[j,i] = struct.unpack( "d", f.read( 8 ) )[0]
f.close()

#] ir-md: http://quantum-machine.org/gdml/doc/applications.html
ac2 = [ np.correlate( v, v, "full" ) for v in velo ]
ac2 /= np.linalg.norm( ac2, axis = 1 )[:, None]
ac2 = np.mean( ac2, axis = 0 )
dos = np.abs( np.fft.fft( ac2 ) ) ** 2
dos /= np.linalg.norm( dos ) / 2
dos = dos[:stp]
frq = np.fft.fftfreq( 2 * stp - 1, dlt )[:stp]
frq /= 1.0e-15 * qm3.constants.C * 100.0 # 1/fs >> 1/cm
import scipy.ndimage.filters
dos = scipy.ndimage.filters.gaussian_filter1d( dos, sigma = 10 )
llm = frq >= 500
ulm = frq <= 4000
dos /= np.max( dos[[ i for i in range( stp ) if llm[i] and ulm[i] ]] )

#] hessian-based spectrum
obj.get_func()
chg = obj.mol.chrg[:]
obj.num_hess()
val, vec = qm3.utils.hessian_frequencies( obj.mass, obj.coor, obj.hess )
print( val )
iri = qm3.utils.intensities( chg, vec )
sx, sy = qm3.utils.spectrum( val, iri, sigm = 50.0, minf = 500.0 )

#] experimental data: http://www.irug.org/jcamp-details?id=1855
dx = []
dy = []
f = open( "data", "rt" )
for l in f:
    t = [ float( i ) for i in l.split() ]
    dx.append( t[0] )
    dy.append( t[1] )
f.close()

#] plot!
plt.clf()
plt.grid( True )
plt.yticks( [] )
plt.xlim( 4000, 500 )
plt.ylim( -0.05, 1.05 )
plt.plot( frq, dos, '-r', linewidth = 2 )
plt.plot( sx, sy, '-b' )
plt.plot( dx, dy, '-k' )
plt.tight_layout()
plt.savefig( "spectrum.pdf" )
plt.show()

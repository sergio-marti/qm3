import  os
import  io
import  math
import  qm3.mol
import  qm3.engines.xtb
import  qm3.actions.ring_polymer


os.environ["QM3_LIBXTB"] = "./bin/xtb.so"


class my_model:
    def __init__( self, molec ):
        self.eqm = qm3.engines.xtb.dl_xtb( molec, 0, 0, list( range( molec.natm ) ) )

    def get_grad( self, molec ):
        self.eqm.get_grad( molec )



class my_rpmd( qm3.actions.ring_polymer.md_template ):
    def __init__( self ):
        qm3.actions.ring_polymer.md_template.__init__( self )

        self.mole = qm3.mol.molecule()
        f = io.StringIO( """6

C       0.02707      0.00232      0.02707
O       0.83930      0.83494      0.83930
H      -0.60585     -0.62664      0.66390
H      -0.60416      0.62975     -0.60416
H       0.66390     -0.62664     -0.60585
H       1.39155      0.25024      1.39155
""" )
        self.mole.xyz_read( f )
        self.mole.guess_atomic_numbers() 
        self.mole.fill_masses()

        self.engn = my_model( self.mole )

        self.sele = list( range( self.mole.natm ) )
        self.size = 3 * len( self.sele )
        self.coor = []
        self.mass = []
        for i in self.sele:
            self.coor += self.mole.coor[3*i:3*i+3]
            self.mass.append( self.mole.mass[i] )

        self.fd = open( "trj", "wt" )


    def current_step( self, istep ):
        if( istep % 10 == 0 ):
            self.fd.write( "%d\n\n"%( self.mole.natm + ( self.rp_bead - 1 ) * self.rp_dime ) )
            for i in range( self.mole.natm ):
                if( not i in self.rp_atom ):
                    self.fd.write( "%-2s%20.10lf%20.10lf%20.10lf\n"%( self.mole.labl[i],
                        self.coor[3*i], self.coor[3*i+1], self.coor[3*i+2] ) )
            for i in self.rp_atom:
                for j in range( self.rp_bead ):
                    self.fd.write( "%-2s%20.10lf%20.10lf%20.10lf\n"%( self.mole.labl[i],
                        self.rp_coor[i][3*j], self.rp_coor[i][3*j+1], self.rp_coor[i][3*j+2] ) )
            





obj = my_rpmd()
obj.setup( [ 5 ], 8 )
import  time
t = time.time()
obj.integrate( step_number = 10000, print_frequency = 1 )
print( time.time() - t )
obj.fd.close()

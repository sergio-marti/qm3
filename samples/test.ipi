import qm3.actions.ipi
import qm3.mol
import qm3.engines.dftb
import qm3.utils._mpi
import os


# i-pi input.xml & sleep 5; mpirun -n 4 python test.ipi


class problem( qm3.actions.ipi.ipi_problem ):
    def __init__( self, node ):
        qm3.actions.ipi.ipi_problem.__init__( self )

        os.mkdir( "%02d"%( node ) )
        os.chdir( "%02d"%( node ) )

        self.mol = qm3.mol.molecule()
        self.mol.xyz_read( "../init.xyz" )
        self.mol.chrg = [ 0.0 for i in range( self.mol.natm ) ]

        self.eng = qm3.engines.dftb.dftb( self.mol, list( range( self.mol.natm ) ) )
        self.eng.exe = "./bin/dftb+ > dftb_in.log"
        self.eng.chg = +1.0
        self.eng.prm = "./bin/3ob-3-1/"

        self.size = 3 * self.mol.natm
        self.coor = self.mol.coor[:]


    def get_func( self ):
        self.mol.func = 0.0
        self.mol.coor = self.coor[:]
        self.eng.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.mol.func = 0.0
        self.mol.coor = self.coor[:]
        self.mol.grad = [ 0.0 for i in range( 3 * self.mol.natm ) ]
        self.eng.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = self.mol.grad


node, ncpu = qm3.utils._mpi.init()
obj = problem( node )
obj.connect( "/tmp/ipi_zundel" )
qm3.utils._mpi.stop()

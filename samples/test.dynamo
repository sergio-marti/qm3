import qm3.engines.dynamo
import qm3.engines.mmint
import qm3.problem
import qm3.maths.matrix
import os
import time

class my_eng( qm3.problem.template ):

    def __init__( self ):
        self.mol = qm3.engines.dynamo.coordinates_read( "crd" )
        self.mol.chrg = [ 0.000, 0.000, 0.000, -0.834, 0.417, 0.417 ]
        self.mol.epsi = [ 0.9374, 0.0000, 0.0000, 0.9374, 0.000, 0.000 ] # Sqrt( eps * 4.184 ) >> Sqrt(kJ/mol)
        self.mol.rmin = [ 1.6610, 0.0000, 0.0000, 1.6610, 0.000, 0.000 ] # rmin/2 >> Ang

        try:
            os.unlink( "dynamo.pipe" )
        except Exception as e:
            print( e )
        os.mkfifo( "dynamo.pipe" )
        os.system( "./a.out dynamo.pipe > dynamo.log &" )
        time.sleep( 10 )
        self.eng = qm3.engines.dynamo.dynamo_pipe()
        self.eng.update_chrg( self.mol )

        self.fix = qm3.engines.mmint.QMLJ( self.mol, [ 0, 1, 2 ], [ 3, 4, 5 ], [] )
        
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
        self.fix.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = self.mol.grad




f = open( "crd", "wt" )
f.write( """6 2 1
Subsystem     1  W
     2
Residue     1  HOH
     3
     1   OH2           8        0.1210000000        0.0690000000       -0.2250000000
     2   H1            1       -0.5270000000        0.1650000000       -0.9460000000
     3   H2            1        0.7020000000       -0.6360000000       -0.5470000000
Residue     2  HOH
     3
     4   OH2           8       -0.4510000000        1.1270000000        2.2110000000
     5   H1            1       -0.2920000000        0.5950000000        1.3990000000
     6   H2            1        0.0580000000        1.9270000000        2.0100000000
""" )
f.close()

f = open( "seq", "wt" )
f.write( """sequence
1
subsystem W
2
HOH 2
end
end
""" )
f.close()

f = open( "opls", "wt" )
f.write( """MM_Definitions OPLS_AA 1.0
Types
HW  1   0.00000  0.00000
OW  8   3.15061  0.15210 
End
Electrostatics Scale 0.5
Lennard_Jones  Scale 0.5
Units kcal/mole
Residues
Residue HOH
  3  3  0
OH2      OW      -0.834
H1       HW       0.417
H2       HW       0.417
OH2 H1 ; OH2 H2 ; H1 H2
End
Parameters
Bonds
HW OW      529.6    0.9572
HW HW       38.25   1.5139
End
Angles
HW OW HW   34.05    104.52
HW HW OW    0.0      37.74
End
Dihedrals
End
End
End
""" )
f.close()


obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )
obj.eng.stop()

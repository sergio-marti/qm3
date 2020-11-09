import qm3.mol
import qm3.engines.sander
import qm3.engines.mmint
import qm3.problem
import qm3.maths.matrix

f = open( "pdb", "wt" )
f.write( """ATOM      1  OH2 HOH     1       0.121   0.069  -0.225  1.00  0.00      W   
ATOM      2  H1  HOH     1      -0.527   0.165  -0.946  1.00  0.00      W   
ATOM      3  H2  HOH     1       0.702  -0.636  -0.547  1.00  0.00      W   
ATOM      4  OH2 HOH     2      -0.451   1.127   2.211  0.00  0.00      W   
ATOM      5  H1  HOH     2      -0.292   0.595   1.399  0.00  0.00      W   
ATOM      6  H2  HOH     2       0.058   1.927   2.010  0.00  0.00      W   
END
""" )
f.close()

f = open( "prmtop", "wt" )
f.write( """%VERSION  VERSION_STAMP = V0001.000  DATE = 05/04/18  16:32:20                  
%FLAG TITLE                                                                     
%FORMAT(20a4)                                                                   
default_name                                                                    
%FLAG POINTERS                                                                  
%FORMAT(10I8)                                                                   
       6       2       6       0       0       0       0       0       0       0
       8       2       0       0       0       2       0       0       2       0
       0       0       0       0       0       0       0       1       3       0
       0
%FLAG ATOM_NAME                                                                 
%FORMAT(20a4)                                                                   
OH2 H1  H2  OH2 H1  H2  
%FLAG CHARGE                                                                    
%FORMAT(5E16.8)                                                                 
 -1.51973982E+01  7.59869910E+00  7.59869910E+00 -1.51973982E+01  7.59869910E+00
  7.59869910E+00
%FLAG ATOMIC_NUMBER                                                             
%FORMAT(10I8)                                                                   
       8       1       1       8       1       1
%FLAG MASS                                                                      
%FORMAT(5E16.8)                                                                 
  1.60000000E+01  1.00800000E+00  1.00800000E+00  1.60000000E+01  1.00800000E+00
  1.00800000E+00
%FLAG ATOM_TYPE_INDEX                                                           
%FORMAT(10I8)                                                                   
       1       2       2       1       2       2
%FLAG NUMBER_EXCLUDED_ATOMS                                                     
%FORMAT(10I8)                                                                   
       2       1       1       2       1       1
%FLAG NONBONDED_PARM_INDEX                                                      
%FORMAT(10I8)                                                                   
       1       2       2       3
%FLAG RESIDUE_LABEL                                                             
%FORMAT(20a4)                                                                   
HOH HOH 
%FLAG RESIDUE_POINTER                                                           
%FORMAT(10I8)                                                                   
       1       4
%FLAG BOND_FORCE_CONSTANT                                                       
%FORMAT(5E16.8)                                                                 
  5.53000000E+02  5.53000000E+02
%FLAG BOND_EQUIL_VALUE                                                          
%FORMAT(5E16.8)                                                                 
  1.51360000E+00  9.57200000E-01
%FLAG ANGLE_FORCE_CONSTANT                                                      
%FORMAT(5E16.8)                                                                 

%FLAG ANGLE_EQUIL_VALUE                                                         
%FORMAT(5E16.8)                                                                 

%FLAG DIHEDRAL_FORCE_CONSTANT                                                   
%FORMAT(5E16.8)                                                                 

%FLAG DIHEDRAL_PERIODICITY                                                      
%FORMAT(5E16.8)                                                                 

%FLAG DIHEDRAL_PHASE                                                            
%FORMAT(5E16.8)                                                                 

%FLAG SCEE_SCALE_FACTOR                                                         
%FORMAT(5E16.8)                                                                 

%FLAG SCNB_SCALE_FACTOR                                                         
%FORMAT(5E16.8)                                                                 

%FLAG SOLTY                                                                     
%FORMAT(5E16.8)                                                                 
  0.00000000E+00  0.00000000E+00
%FLAG LENNARD_JONES_ACOEF                                                       
%FORMAT(5E16.8)                                                                 
  5.81935564E+05  0.00000000E+00  0.00000000E+00
%FLAG LENNARD_JONES_BCOEF                                                       
%FORMAT(5E16.8)                                                                 
  5.94825035E+02  0.00000000E+00  0.00000000E+00
%FLAG BONDS_INC_HYDROGEN                                                        
%FORMAT(10I8)                                                                   
       3       6       1       0       3       2       0       6       2      12
      15       1       9      12       2       9      15       2
%FLAG BONDS_WITHOUT_HYDROGEN                                                    
%FORMAT(10I8)                                                                   

%FLAG ANGLES_INC_HYDROGEN                                                       
%FORMAT(10I8)                                                                   

%FLAG ANGLES_WITHOUT_HYDROGEN                                                   
%FORMAT(10I8)                                                                   

%FLAG DIHEDRALS_INC_HYDROGEN                                                    
%FORMAT(10I8)                                                                   

%FLAG DIHEDRALS_WITHOUT_HYDROGEN                                                
%FORMAT(10I8)                                                                   

%FLAG EXCLUDED_ATOMS_LIST                                                       
%FORMAT(10I8)                                                                   
       2       3       3       0       5       6       6       0
%FLAG HBOND_ACOEF                                                               
%FORMAT(5E16.8)                                                                 

%FLAG HBOND_BCOEF                                                               
%FORMAT(5E16.8)                                                                 

%FLAG HBCUT                                                                     
%FORMAT(5E16.8)                                                                 

%FLAG AMBER_ATOM_TYPE                                                           
%FORMAT(20a4)                                                                   
OW  HW  HW  OW  HW  HW  
%FLAG TREE_CHAIN_CLASSIFICATION                                                 
%FORMAT(20a4)                                                                   
BLA BLA BLA BLA BLA BLA 
%FLAG JOIN_ARRAY                                                                
%FORMAT(10I8)                                                                   
       0       0       0       0       0       0
%FLAG IROTAT                                                                    
%FORMAT(10I8)                                                                   
       0       0       0       0       0       0
%FLAG SOLVENT_POINTERS                                                          
%FORMAT(3I8)                                                                    
       0       2       1
%FLAG ATOMS_PER_MOLECULE                                                        
%FORMAT(10I8)                                                                   
       3       3
%FLAG BOX_DIMENSIONS                                                            
%FORMAT(5E16.8)                                                                 
  9.00000000E+01  1.00000000E+02  1.00000000E+02  1.00000000E+02
%FLAG RADIUS_SET                                                                
%FORMAT(1a80)                                                                   
modified Bondi radii (mbondi)                                                   
%FLAG RADII                                                                     
%FORMAT(5E16.8)                                                                 
  1.50000000E+00  1.50000000E+00  1.50000000E+00  1.50000000E+00  1.50000000E+00
  1.50000000E+00
%FLAG SCREEN                                                                    
%FORMAT(5E16.8)                                                                 
  8.00000000E-01  8.00000000E-01  8.00000000E-01  8.00000000E-01  8.00000000E-01
  8.00000000E-01
%FLAG IPOL                                                                      
%FORMAT(1I8)                                                                    
       0
""" )
f.close()


##############################################################################################
# sander: MM
class my_eng( qm3.problem.template ):
    def __init__( self ):
        self.mol = qm3.mol.molecule( "pdb" )
        qm3.engines.sander.topology_read( self.mol, "prmtop" )
        self.mol.boxl = [ 100, 100, 100 ]
        self.mol.chrg = [ 0.000, 0.000, 0.000, -0.834, 0.417, 0.417 ]
        self.mol.anum = [ 8, 1, 1, 8, 1, 1 ]
        self.mol.epsi = [ 0.9374, 0.0000, 0.0000, 0.9374, 0.000, 0.000 ] # Sqrt( eps * 4.184 ) >> Sqrt(kJ/mol)
        self.mol.rmin = [ 1.6610, 0.0000, 0.0000, 1.6610, 0.000, 0.000 ] # rmin/2 >> Ang

        self.eng = qm3.engines.sander.sander()
        self.eng.exe = "./bin/sander -O"
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


f = open( "mdin", "wt" )
f.write( """slave_frozen_QM
&cntrl
 imin   = 0, ntx    = 1, irest  = 0,
 ntxo   = 1, ntpr   = 1, ntave  = 0,
 ntwr   = 1, iwrap  = 0, ntwx   = 0,
 ntwv   = 0, ntwf   = 1, ntwe   = 1,
 ioutfm = 1, ntwprt = 0, ntt    = 0,
 ntr    = 0, nstlim = 0, nscm   = 0,
 ntp    = 0, ntb    = 0, ifqnt  = 0,
 cut       = 16.0,
 ibelly    = 1, 
 bellymask = '!:1',
/

""" )
f.close()

obj = my_eng()
obj.get_grad()
print( obj.func )
import qm3.maths.matrix
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )


##############################################################################################
# Python: MM
class my_eng( qm3.problem.template ):

    def __init__( self ):
        self.mol = qm3.mol.molecule( "pdb" )
        self.mol.boxl = [ 100, 100, 100 ]

        self.eng = qm3.engines.sander.py_sander( self.mol, "prmtop",
             cutoff = 16.0, PBC = False, qmsel = [ 0, 1, 2 ], method = "EXTERN" )

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



obj = my_eng()
obj.get_grad()
print( obj.func )
import qm3.maths.matrix
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )


##############################################################################################
# Extern(QM) MM
class my_eng( qm3.problem.template ):
    def __init__( self ):
        self.mol = qm3.mol.molecule( "pdb" )
        qm3.engines.sander.topology_read( self.mol, "prmtop" )
        self.mol.boxl = [ 100, 100, 100 ]
        self.mol.chrg = [ 0.000, 0.000, 0.000, -0.834, 0.417, 0.417 ]
        self.mol.anum = [ 8, 1, 1, 8, 1, 1 ]
        self.mol.epsi = [ 0.9374, 0.0000, 0.0000, 0.9374, 0.000, 0.000 ] # Sqrt( eps * 4.184 ) >> Sqrt(kJ/mol)
        self.mol.rmin = [ 1.6610, 0.0000, 0.0000, 1.6610, 0.000, 0.000 ] # rmin/2 >> Ang

        self.eng = qm3.engines.sander.sander()
        self.eng.exe = "./bin/sander -O"
        self.eng.update_chrg( self.mol )

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


f = open( "mdin", "wt" )
f.write( """slave_frozen_QM
&cntrl
 imin   = 0, ntx    = 1, irest  = 0,
 ntxo   = 1, ntpr   = 1, ntave  = 0,
 ntwr   = 1, iwrap  = 0, ntwx   = 0,
 ntwv   = 0, ntwf   = 1, ntwe   = 1,
 ioutfm = 1, ntwprt = 0, ntt    = 0,
 ntr    = 0, nstlim = 0, nscm   = 0,
 ntp    = 0, ntb    = 0, ifqnt  = 1,
 cut       = 16.0,
/
&qmmm
 qmcut     = 16.0,
 qmmask    = ':1',
 qmcharge  = 0,
 spin      = 1,
 qm_ewald  = 0,
 qm_theory = 'EXTERN',
/
""" )
f.close()

obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )

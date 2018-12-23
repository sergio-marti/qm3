import qm3.mol
import qm3.engines.openmm
import qm3.maths.matrix

import simtk.openmm
import simtk.openmm.app
import simtk.unit

mol = qm3.mol.molecule( "pdb" )
prm = simtk.openmm.app.ForceField( "prm.xml" )
psf = simtk.openmm.app.CharmmPsfFile( "psf" )
omm = prm.createSystem( psf.topology, nonbondedMethod = simtk.openmm.app.CutoffNonPeriodic, nonbondedCutoff = 6. * simtk.unit.angstrom, rigidWater = False )
eng = qm3.engines.openmm.py_openmm( omm, psf.topology, [0, 1, 2] )
mol.chrg = [ .0, .0, .0, -0.834, 0.417, 0.417 ] 
eng.update_chrg( mol )
mol.func = 0.0
mol.grad = [ 0.0 for i in range( 3 * mol.natm ) ]
eng.get_grad( mol )
print( mol.func )
qm3.maths.matrix.mprint( mol.grad, mol.natm, 3 )

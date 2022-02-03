# -*- coding: iso-8859-1 -*-
import sys
import os


sys.path.insert( 0, os.getenv( "QM3_OPENMM" ) )
try:
    import openmm
    import openmm.app
    import openmm.unit
    class py_openmm( object ):
        def __simulation( self ):
            self.sim = openmm.app.Simulation( self.top, self.sys,
                openmm.CustomIntegrator( 0.001 ),
                openmm.Platform.getPlatformByName( self.knd ) )
            
        def __init__( self, omm_system, topology, qm_excl = [], platform = "CPU" ):
            self.sys = omm_system
            self.top = topology
            self.knd = platform
            self.nbn = None
            self.sim = None
            for i in range( self.sys.getNumForces() ):
                if( type( self.sys.getForce( i ) ) == openmm.NonbondedForce ):
                    self.nbn = self.sys.getForce( i )
            try:
                n = len( qm_excl )
                for i in range( 0, n - 1 ):
                    for j in range( i + 1, n ):
                        self.nbn.addException( qm_excl[i], qm_excl[j], 0.0, 0.0, 0.0 )
            except:
                pass
            self.__simulation()


        def update_chrg( self, mol ):
            for i in range( mol.natm ):
                t = self.nbn.getParticleParameters( i )
                self.nbn.setParticleParameters( i, mol.chrg[i], t[1], t[2] )
            self.__simulation()
            self.update_coor( mol )


        def update_coor( self, mol ):
            tmp = []
            for i in range( mol.natm ):
                i3 = i * 3
                tmp.append( openmm.Vec3( mol.coor[i3], mol.coor[i3+1], mol.coor[i3+2] ) * openmm.unit.angstrom )
            self.sim.context.setPositions( tmp )


        def get_func( self, mol ):
            self.update_coor( mol )
            stt = self.sim.context.getState( getEnergy = True, getForces = False )
            mol.func += stt.getPotentialEnergy().value_in_unit( openmm.unit.kilojoule/openmm.unit.mole )


        def get_grad( self, mol ):
            self.update_coor( mol )
            stt = self.sim.context.getState( getEnergy = True, getForces = True )
            mol.func += stt.getPotentialEnergy().value_in_unit( openmm.unit.kilojoule/openmm.unit.mole )
            frc = stt.getForces()
            for i in range( mol.natm ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    mol.grad[i3+j] -= frc[i][j].value_in_unit( openmm.unit.kilojoule/(openmm.unit.angstrom*openmm.unit.mole) )

except:
    pass

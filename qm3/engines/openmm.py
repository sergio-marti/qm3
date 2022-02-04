# -*- coding: iso-8859-1 -*-
import sys
import os


sys.path.insert( 0, os.getenv( "QM3_OPENMM" ) )
try:
    import openmm
    import openmm.app
    import openmm.unit

    class py_openmm( object ):
        def __init__( self, omm_system, topology, qm_atom = [], platform = "CPU" ):
            self.sys = omm_system
            self.top = topology
            self.knd = platform
            self.nbn = None
            self.sim = None

            nqm = len( qm_atom )
            if( nqm > 0 ):
                msk = [ 0 for i in range( self.sys.getNumParticles() ) ]
                for i in qm_atom:
                    msk[i] = 1
                for ii in range( self.sys.getNumForces() ):
                    cur = self.sys.getForce( ii )
                    if( type( cur ) == openmm.NonbondedForce ):
                        for i in range( 0, nqm - 1 ):
                            for j in range( i + 1, nqm ):
                                cur.addException( qm_atom[i], qm_atom[j], 0.0, 0.0, 0.0, replace = True )
                        self.nbn = cur
                    elif( type( cur ) == openmm.HarmonicBondForce ):
                        for i in range( cur.getNumBonds() ):
                            tmp = cur.getBondParameters( i )
                            if( msk[tmp[0]] == 1 and msk[tmp[1]] == 1 ):
                                cur.setBondParameters( i, tmp[0], tmp[1], 0.0, 0.0 )
                    elif( type( cur ) == openmm.HarmonicAngleForce ):
                        for i in range( cur.getNumAngles() ):
                            tmp = cur.getAngleParameters( i )
                            if( msk[tmp[0]] + msk[tmp[1]] + msk[tmp[2]] >= 2 ):
                                cur.setAngleParameters( i, tmp[0], tmp[1], tmp[2], 0.0, 0.0 )
                    elif( type( cur ) == openmm.PeriodicTorsionForce ):
                        for i in range( cur.getNumTorsions() ):
                            tmp = cur.getTorsionParameters( i )
                            if( msk[tmp[0]] + msk[tmp[1]] + msk[tmp[2]] + msk[tmp[3]] >= 3 ):
                                cur.setTorsionParameters( i, tmp[0], tmp[1], tmp[2], tmp[3], 1, 0.0, 0.0 )
                    elif( type( cur ) == openmm.CMMotionRemover ):
                        pass
                    else:
                        print( ">> Unhandled QM atoms at: %s [%d]"%( type( cur ), ii ) )

            if( self.nbn == None ):
                for i in range( self.sys.getNumForces() ):
                    if( type( self.sys.getForce( i ) ) == openmm.NonbondedForce ):
                        self.nbn = self.sys.getForce( i )

            self.sim = openmm.app.Simulation( self.top, self.sys,
                openmm.CustomIntegrator( 0.001 ),
                openmm.Platform.getPlatformByName( self.knd ) )


        def update_chrg( self, mol ):
            for i in range( mol.natm ):
                t = self.nbn.getParticleParameters( i )
                self.nbn.setParticleParameters( i, mol.chrg[i], t[1], t[2] )
            self.nbn.updateParametersInContext( self.sim.context )


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

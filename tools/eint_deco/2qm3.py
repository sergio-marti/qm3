#!/usr/bin/env python3
import  pickle
import  qm3.mol
import  qm3.engines.sander

mol = qm3.mol.molecule()
bnd = qm3.engines.sander.topology_read( mol, "prmtop" )
f = open( "bonds.pk", "wb" )
pickle.dump( bnd, f )
f.close()
f = open( "molec.pk", "wb" )
pickle.dump( mol, f )
f.close()

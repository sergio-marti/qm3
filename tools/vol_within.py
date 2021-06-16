#!/usr/bin/env python3
import  qm3.mol
import  qm3.engines.dynamo
import  qm3.elements
import  struct

mol = qm3.engines.dynamo.coordinates_read( "../00.r" )

sus = list( mol.indx["X"][1].values() )

f = open( "within.acs", "wb" )
f.write( struct.pack( "i", len( sus ) ) )
for i in sus:
    i3 = i * 3
    for j in [0, 1, 2]:
        f.write( struct.pack( "d", mol.coor[i3+j] ) )
    f.write( struct.pack( "d", qm3.elements.r_vdw[mol.anum[i]] ) )
f.close()

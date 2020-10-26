# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import re
import qm3.fio
from qm3.engines._mmint import *


# >> add a parser function for charmm/par, dynamo/opls?

def non_bonded( mol, fname ):
    out = True
    if( mol.type == [] ):
        print( "- Molecule types undefined!" )
        return( False )
    mol.epsi = []
    mol.rmin = []
    nbd = {}
    f = qm3.fio.open_r( fname )
    for i,j,k in re.compile( "([A-Z0-9]+)[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)" ).findall( f.read() ):
        nbd[i] = [ math.sqrt( float( j ) * qm3.constants.K2J ), float( k ) ]
    qm3.fio.close( f, fname )
    for i in range( mol.natm ):
        if( mol.type[i] in nbd ):
            mol.epsi.append( nbd[mol.type[i]][0] )
            mol.rmin.append( nbd[mol.type[i]][1] )
        else:
            mol.epsi.append( None )
            mol.rmin.append( None )
            print( "- Atom index %d misses Non-Bonded data..."%( i+1 ) )
            out = False
    return( out )



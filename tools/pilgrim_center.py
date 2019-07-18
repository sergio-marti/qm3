#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange

import	qm3.mol
import	qm3.utils

import	qm3.engines.dynamo


mol = qm3.engines.dynamo.coordinates_read( "sadd" )
mol.fill_masses()


sel = []
for a in [ "C9", "H13", "H14", "C10", "C11", "H15", "H16", "C12", "H17", "H18", "C19", "H22", "H23", "H24", "C18", "C13", "C14", "H19", "C15", "O2", "C16", "H20", "C17", "H21", "H26" ]:
	sel.append( mol.indx["MOL"][1][a] )
for a in [ "C1", "H11", "H12", "N1X", "C2X", "O2X", "N3X", "H3X", "C4Y", "O4X", "N10", "C10A", "C4A", "N5", "C5A", "C9A", "C9", "H9", "C8X", "C7", "C6X", "H6", "C7M", "H7M1", "H7M2", "H7M3", "C8M", "H8M1", "H8M2", "H8M3" ]:
	sel.append( mol.indx["FAD"][1][a] )
sel.sort()


mas = []
crd = []
for i in sel:
	mas.append( mol.mass[i] )
	crd += mol.coor[3*i:3*i+3]
cen = qm3.utils.center( mas, crd )
for i in range( mol.natm ):
	i3 = i * 3
	for j in [0, 1, 2]:
		mol.coor[i3+j] -= cen[j]


qm3.engines.dynamo.coordinates_write( mol, "sadd.cen" )
mol.xyz_write( "acs", sele = sel )

#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
from __future__ import print_function, division
import  sys
if( sys.version_info[0] == 2 ):
    range = xrange

import  qm3.mol
import  qm3.engines.namd
import  qm3.elements
import  pickle

m = qm3.mol.molecule( sys.argv[1] )
m.boxl = [ float( sys.argv[2] ), float( sys.argv[3] ), float( sys.argv[4] ) ]

# fix each axis by step
for W in [0, 1, 2]:
    ori = [ .0, .0, .0 ]
    npt = 0
    # segment_1: protein, segment_2: substrate
    for i in range( m.res_lim[m.seg_lim[2]] ):
        for j in [0, 1, 2]:
            ori[j] += m.coor[3*i+j]
        npt += 1
    ori = [ i / npt for i in ori ]
    for i in range( m.natm ):
        m.coor[3*i+W] -= ori[W]
    # segnemt_3...: counterions and waters
    for i in range( m.seg_lim[2], m.seg_lim[-1] ):
        cen = [ .0, .0, .0 ]
        npt = 0
        for j in range( m.res_lim[i], m.res_lim[i+1] ):
            for k in [0, 1, 2]:
                cen[k] += m.coor[3*j+k]
            npt += 1
        cen = [ ( cen[j] / npt ) for j in [0, 1, 2] ]
        if( round( cen[W] / m.boxl[W], 0 ) ):
            for j in range( m.res_lim[i], m.res_lim[i+1] ):
                m.coor[3*j+W] -= m.boxl[W] * round( m.coor[3*j+W] / m.boxl[W], 0 )

m.pdb_write( "pdb" )

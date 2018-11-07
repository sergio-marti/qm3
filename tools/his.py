#!/usr/bin/env python

import	sys
import	qm3.mol
import	qm3.utils


m = qm3.mol.molecule( sys.argv[1] )

his = {}
for i in xrange( m.natm ):
	if( m.resn[i] == "HIS" and m.labl[i] == "ND1" ):
		if( not m.segn[i] in his ):
			his[m.segn[i]] = {}
		if( not m.resi[i] in his[m.segn[i]] ):
			his[m.segn[i]][m.resi[i]] = {}
		if( not "ND1" in his[m.segn[i]][m.resi[i]] ):
			his[m.segn[i]][m.resi[i]]["ND1"] = i
	elif( m.resn[i] == "HIS" and m.labl[i] == "NE2" ):
		if( not m.segn[i] in his ):
			his[m.segn[i]] = {}
		if( not m.resi[i] in his[m.segn[i]] ):
			his[m.segn[i]][m.resi[i]] = {}
		if( not "NE2" in his[m.segn[i]][m.resi[i]] ):
			his[m.segn[i]][m.resi[i]]["NE2"] = i

cut = 3.2
for chain in sorted( his ):
	for resi in sorted( his[chain] ):
		jd = m.indx[chain][resi]["ND1"]
		je = m.indx[chain][resi]["NE2"]
		print chain, resi
		sd = m.sph_sel( [ his[chain][resi]["ND1"] ], 4. )
		se = m.sph_sel( [ his[chain][resi]["NE2"] ], 4. )
		print "\t- HSD:"
		for i in sd:
			if( m.labl[i] in [ "O", "SD" ] ):
				t = round( qm3.utils.distance( m.coor[3*jd:3*jd+3], m.coor[3*i:3*i+3] ), 3 )
				if( t <= cut  ):
					print "\t\t[Hd] %4s %4d %4s %4s"%( m.segn[i], m.resi[i], m.resn[i], m.labl[i] ), t
			elif( m.labl[i][0:2] in [ "OE", "OD" ] ):
				t = round( qm3.utils.distance( m.coor[3*jd:3*jd+3], m.coor[3*i:3*i+3] ), 3 )
				if( t <= cut  ):
					print "\t\t[Hd] %4s %4d %4s %4s"%( m.segn[i], m.resi[i], m.resn[i], m.labl[i] ), t
		for i in se:
			if( m.labl[i] in [ "N", "SG", "OG", "NZ", "NE", "NH1", "NH2", "ND2", "NE1" ] ):
				t = round( qm3.utils.distance( m.coor[3*je:3*je+3], m.coor[3*i:3*i+3] ), 3 )
				if( t <= cut  ):
					print "\t\t[ e] %4s %4d %4s %4s"%( m.segn[i], m.resi[i], m.resn[i], m.labl[i] ), t
			elif( m.labl[i] == "NE2" and m.resn[i] == "TRP" ):
				t = round( qm3.utils.distance( m.coor[3*je:3*je+3], m.coor[3*i:3*i+3] ), 3 )
				if( t <= cut  ):
					print "\t\t[ e] %4s %4d %4s %4s"%( m.segn[i], m.resi[i], m.resn[i], m.labl[i] ), t
		print "\n\t- HSE:"
		for i in sd:
			if( m.labl[i] in [ "N", "SG", "OG", "NZ", "NE", "NH1", "NH2", "ND2", "NE1" ] ):
				t = round( qm3.utils.distance( m.coor[3*jd:3*jd+3], m.coor[3*i:3*i+3] ), 3 )
				if( t <= cut  ):
					print "\t\t[ d] %4s %4d %4s %4s"%( m.segn[i], m.resi[i], m.resn[i], m.labl[i] ), t
			elif( m.labl[i] == "NE2" and m.resn[i] == "TRP" ):
				t = round( qm3.utils.distance( m.coor[3*jd:3*jd+3], m.coor[3*i:3*i+3] ), 3 )
				if( t <= cut  ):
					print "\t\t[ d] %4s %4d %4s %4s"%( m.segn[i], m.resi[i], m.resn[i], m.labl[i] ), t
		for i in se:
			if( m.labl[i] in [ "O", "SD" ] ):
				t = round( qm3.utils.distance( m.coor[3*je:3*je+3], m.coor[3*i:3*i+3] ), 3 )
				if( t <= cut  ):
					print "\t\t[He] %4s %4d %4s %4s"%( m.segn[i], m.resi[i], m.resn[i], m.labl[i] ), t
			elif( m.labl[i][0:2] in [ "OE", "OD" ] ):
				t = round( qm3.utils.distance( m.coor[3*je:3*je+3], m.coor[3*i:3*i+3] ), 3 )
				if( t <= cut  ):
					print "\t\t[He] %4s %4d %4s %4s"%( m.segn[i], m.resi[i], m.resn[i], m.labl[i] ), t
		print 80*"#"

# -*- coding: iso-8859-1 -*-
from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange


QM_engines = {
	"Gaussian": "qm3.engines.gaussian.gaussian",
	"Amber SQM": "qm3.engines.sqm.dl_sqm",
	"deMon2k": "qm3.engines.demon.demon",
	"fDynamo": "qm3.engines.dynamo.py_dynamo",
	"Orca": "qm3.engines.orca.orca",
	"NWChem": "qm3.engines.nwchem.nwchem",
	"DFTB+": "qm3.engines.dftb.dftb" }


MM_engines = {
	"fDynamo": "qm3.engines.dynamo.py_dynamo",
	"CHARMM": "qm3.engines.charmm.charmm_shm",
	"NAMD": "qm3.engines.namd.namd_shm" }

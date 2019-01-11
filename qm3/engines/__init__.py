# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	qm3.constants



class qmbase( object ):

	def __init__( self, mol, sele, nbnd = [], link = [] ):
		self._cx = qm3.constants.A0
		self._ce = qm3.constants.H2J
		self._cg = self._ce / qm3.constants.A0
		self._ch = self._cg / qm3.constants.A0

		self.sel = sele[:]
		self.lnk = link[:]
#		t = [ j for i,j in self.lnk ]
#		self.nbn = [ i for i in nbnd if not i in t ]
		self.nbn = sorted( set( nbnd ).difference( set( sele + sum( link, [] ) ) ) )

		self.exe = ""
		self.vla = []
#		self.smb = [ i.title() for i in mol.guess_symbols( sele ) ]
		self.smb = mol.guess_symbols( sele )

		if( not mol.chrg ):
			mol.chrg = [ 0.0 for i in range( mol.natm ) ]


	def mk_input( self, mol, run ):
		pass


	def parse_log( self, mol, run ):
		pass


	def get_func( self, mol ):
		self.mk_input( mol, "ener" )
		os.system( self.exe )
		self.parse_log( mol, "ener" )


	def get_grad( self, mol ):
		self.mk_input( mol, "grad" )
		os.system( self.exe )
		self.parse_log( mol, "grad" )


	def get_hess( self, mol ):
		self.mk_input( mol, "hess" )
		os.system( self.exe )
		self.parse_log( mol, "hess" )

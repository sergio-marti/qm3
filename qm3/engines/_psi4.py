# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	math
import	qm3.constants
import	qm3.utils



#
# Set environment variable: PSI_SCRATCH
#
sys.path.insert( 0, os.getenv( "QM3_PSI4" ) )
try:
	import	psi4

	class py_psi4( object ):

		def __init__( self, mol, sele, opts = { "reference": "rks", "basis": "6-31g*", "d_convergence": 6, "scf_type": "direct", "guess": "read",
												"output": False, "charge": 0, "method": "b3lyp", "ncpus": 1, "memory": "1024 MB" }, nbnd = [], link = [] ):
			self._ce = qm3.constants.H2J
			self._cg = self._ce / qm3.constants.A0
			self.sel = sele[:]
			self.lnk = link[:]
			t = [ j for i,j in self.lnk ]
			self.nbn = [ i for i in nbnd if not i in t ]
	
			self.smb = mol.guess_symbols( sele )
			buf = "\n%d 1\n"%( opts.pop( "charge" ) )
			j = 0
			for i in sele:
				i3 = i * 3
				buf += "%-2s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j], 
					mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
					mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
					mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) )
				j += 1
			self.vla = []
			if( self.lnk ):
				k = len( self.sel )
				for i,j in self.lnk:
					c, v = qm3.utils.LA_coordinates( i, j, mol )
					buf += "%-2s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] )
					self.vla.append( ( self.sel.index( i ), k, v[:] ) )
					k += 1
	
			# -- psi4 --
			if( opts.pop( "output" ) ):
				psi4.core.set_output_file( "psi4.out", False )
			else:
				psi4.core.be_quiet()
			psi4.set_memory( opts.pop( "memory" ) )
			psi4.set_num_threads( opts.pop( "ncpus" ) )
			buf += "symmetry c1\nno_reorient\nno_com\n"
			self.QMatm = psi4.geometry( buf )
			psi4.activate( self.QMatm )
			if( self.nbn ):
				self.MMatm = psi4.QMMM()
				self.MMatm.charges = []
				for i in self.nbn:
					i3 = i * 3
					self.MMatm.charges.append( [ mol.chrg[i],
						mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
						mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
						mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ] )
				self.MMatm.populateExtern()
				psi4.core.set_global_option_python( "EXTERN", self.MMatm.extern )
			self.met = opts.pop( "method" )
			psi4.set_options( opts )
	
	
		def update_coor( self, mol ):
			crd = []
			for i in self.sel:
				i3 = i * 3
				crd.append( [ ( mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) ) / qm3.constants.A0 for j in [0, 1, 2] ] )
			self.vla = []
			if( self.lnk ):
				k = len( self.sel )
				for i,j in self.lnk:
					c, v = qm3.utils.LA_coordinates( i, j, mol )
					crd.append( [ c[j] / qm3.constants.A0 for j in [0, 1, 2] ] )
					self.vla.append( ( self.sel.index( i ), k, v[:] ) )
					k += 1
			self.QMatm.set_geometry( psi4.core.Matrix.from_list( crd ) )
			self.QMatm.update_geometry()
			if( self.nbn ):
				self.MMatm = psi4.QMMM()
				self.MMatm.charges = []
				f = open( "grid.dat", "wt" )
				for i in self.nbn:
					i3 = i * 3
					t = [ mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) for j in [0, 1, 2] ]
					self.MMatm.charges.append( [ mol.chrg[i], t[0], t[1], t[2] ] )
					f.write( "%20.10lf%20.10lf%20.10lf\n"%( t[0], t[1], t[2] ) )
				f.close()
				self.MMatm.populateExtern()
				psi4.core.set_global_option_python( "EXTERN", self.MMatm.extern )
	
	
		def get_func( self, mol ):
			self.update_coor( mol )
			mol.func += psi4.energy( self.met, resturn_wfn = False ) * self._ce
	
	
		def get_grad( self, mol ):
			self.update_coor( mol )
			g, wfn = psi4.gradient( self.met, return_wfn = True )
			mol.func += psi4.get_variable( 'CURRENT ENERGY' ) * self._ce
			g = sum( g.to_array().tolist(), [] )
			qm3.utils.LA_gradient( self.vla, g )
			for i in range( len( self.sel ) ):
				i3 = i * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.sel[i]+j] += g[i3+j] * self._cg
			if( self.nbn ):
				ef = psi4.core.OEProp( wfn )
				ef.add( "GRID_FIELD" )
				ef.compute()
				efx = ef.Exvals()
				efy = ef.Eyvals()
				efz = ef.Ezvals()
				for i in range( len( self.nbn ) ):
					mol.grad[3*self.nbn[i]]   -= self._cg * mol.chrg[self.nbn[i]] * efx[i]
					mol.grad[3*self.nbn[i]+1] -= self._cg * mol.chrg[self.nbn[i]] * efy[i]
					mol.grad[3*self.nbn[i]+2] -= self._cg * mol.chrg[self.nbn[i]] * efz[i]
				

except:
	pass


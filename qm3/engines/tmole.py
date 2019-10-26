# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	os
import	math
import	qm3.engines



class tmole( qm3.engines.qmbase ):

	def __init__( self, mol, sele, nbnd = [], link = [] ):
		qm3.engines.qmbase.__init__( self, mol, sele, nbnd, link )
		self.exe_ene = ". ./tmole.rc; dscf 1>  tmole.log 2>> tmole.log"
#		self.exe_ene = ". ./tmole.rc; dscf_omp 1>  tmole.log 2>> tmole.log"
		self.exe_grd = ". ./tmole.rc; grad 1>> tmole.log 2>> tmole.log"


	def mk_input( self, mol, run ):
		f = open( "coord", "wt" )
		f.write( "$coord\n" )
		j = 0
		for i in self.sel:
			i3 = i * 3
			f.write( "%20.10lf%20.10lf%20.10lf%4s\n"%(
				( mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ) ) / self._cx, 
				( mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ) ) / self._cx,
				( mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) / self._cx, self.smb[j] ) )
			j += 1
		if( self.lnk ):
			self.vla = []
			k = len( self.sel )
			for i,j in self.lnk:
				c, v = qm3.engines.LA_coordinates( i, j, mol )
				f.write( "%20.10lf%20.10lf%20.10lf   H\n"%( c[0] / self._cx, c[1] / self._cx, c[2] / self._cx ) )
				self.vla.append( ( self.sel.index( i ), k, v[:] ) )
				k += 1
		f.write( "$user-defined bonds\n$end" )
		f.close()
		if( self.nbn ):
			f = open( "charges", "wt" )
			f.write( "$point_charges nocheck\n" )
			for i in self.nbn:
				i3 = i * 3
				f.write( "%20.10lf%20.10lf%20.10lf%12.4lf\n"%(
					( mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ) ) / self._cx, 
					( mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ) ) / self._cx,
					( mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ) / self._cx, mol.chrg[i] ) )
			f.write( "$end" )
			f.close()


	def parse_log( self, mol, run ):
		f = open( "energy", "rt" )
		f.readline()
		mol.func += float( f.readline().split()[1] ) * self._ce
		f.close()
		os.unlink( "energy" )
		if( run == "grad" ):
			f = open( "gradient", "rt" )
			for i in range( 2 + len( self.sel ) + len( self.lnk ) ):
				f.readline()
			g = []
			for i in range( len( self.sel ) + len( self.lnk ) ):
				g += [ float( j.replace( "D", "E" ) ) * self._cg for j in f.readline().strip().split() ]
			qm3.engines.LA_gradient( self.vla, g )
			for i in range( len( self.sel ) ):
				i3 = i * 3
				for j in [0, 1, 2]:
					mol.grad[3*self.sel[i]+j] += g[i3+j]
			f.close()
			os.unlink( "gradient" )
			# read MM gradient
			if( self.nbn ):
				f = open( "charges.gradient", "rt" )
				f.readline()
				for i in range( len( self.nbn ) ):
					t = [ float( j.replace( "D", "E" ) ) * self._cg for j in f.readline().strip().split() ]
					mol.grad[3*self.nbn[i]]   += t[0]
					mol.grad[3*self.nbn[i]+1] += t[1]
					mol.grad[3*self.nbn[i]+2] += t[2]
				f.close()
				os.unlink( "charges.gradient" )


	def get_grad( self, mol ):
		self.mk_input( mol, "grad" )
		os.system( self.exe_ene )
		os.system( self.exe_grd )
		self.parse_log( mol, "grad" )






"""
# tmole.rc  ------------------------------------------------

export TURBODIR=/usr/local/chem/tmole_7.1
export PATH=$TURBODIR/bin:$PATH

# PARALLEL ------------------------------------------ (*_omp)
export LD_LIBRARY_PATH=$TURBODIR/bin:$LD_LIBRARY_PATH
unset OMP_DYNAMIC
export OMP_NUM_THREADS=4

# ----------------------------------------------------------

# system setup (previous to any calculation)  --------------

x2t xyz > coord
define << EOD

slave
a coord
*
no
b all def-SVP       << BASIS SET
*
eht
y
0                   << MOLECULAR CHARGE
y
dft
func b3-lyp         << FUNCTIONAL
on
*
*
EOD


# edit & change "control" file for QM/MM calculations  -----

$drvopt
   point charges

$point_charges file=charges
$point_charge_gradients file=charges.gradient
$end

# ----------------------------------------------------------
"""



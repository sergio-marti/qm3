# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	math
import	qm3.constants
import	qm3.io
import	qm3.mol
import	os, stat
import	struct




def coordinates_read( fname = None ):
	"""
    1 COP   C1         1    2.122800112    1.917699933    2.044047117
||||| ----- |||||       |||||||||||||||---------------|||||||||||||||
.123456789.123456789.123456789.123456789.123456789.123456789.123456789
	"""
	mol = qm3.mol.molecule()
	f = qm3.io.open_r( fname )
	l = f.readline()
	while( l.strip() != "POSITION" ):
		l = f.readline()
	l = f.readline()
	while( l.strip() != "END" ):
		mol.natm += 1
		mol.segn.append( "A" )
		mol.resi.append( int( l[0:5] ) )
		mol.resn.append( l[6:11].strip() )
		mol.labl.append( l[12:17].strip() )
		mol.coor.append( float( l[24:39].strip() ) * 10.0 )
		mol.coor.append( float( l[39:54].strip() ) * 10.0 )
		mol.coor.append( float( l[54:69].strip() ) * 10.0 )
		l = f.readline()
	l = f.readline()
	if( l.strip() == "BOX" ):
		mol.boxl = [ float( i ) * 10.0 for i in f.readline().strip().split() ]
	qm3.io.close( f, fname )
	mol.settle()
	return( mol )
	

def coordinates_write( mol, fname = None ):
	f = qm3.io.open_w( fname )
	f.write( "TITLE\nGenerated with QM3\nEND\nPOSITION\n" )
	for i in range( len( mol.res_lim ) - 1 ):
		for j in range( mol.res_lim[i], mol.res_lim[i+1] ):
			f.write( "%5d %-5s %-5s%7d%15.9lf%15.9lf%15.9lf\n"% ( ( i + 1 ) % 100000,
				mol.resn[j], mol.labl[j], ( j + 1 ) % 10000000,
				# GROMCAS uses nanometers...
				mol.coor[3*j] / 10.0, mol.coor[3*j+1] / 10.0, mol.coor[3*j+2] / 10.0 ) )
	f.write( "END\n" )
	if( mol.boxl != [ qm3.mol.MXLAT, qm3.mol.MXLAT, qm3.mol.MXLAT ]  ):
		f.write( "BOX\n%15.9lf%15.9lf%15.9lf\nEND\n"%( mol.boxl[0] / 10.0, mol.boxl[1] / 10.0, mol.boxl[2] / 10.0 ) )
	qm3.io.close( f, fname )



class gromacs( object ):
	"""
	rm -f \#* borra.*

	gmx editconf -f cope.pdb -box 4.0 4.0 4.0 -c -o borra.g96

	gmx solvate -cp borra.g96 -box 4.0 4.0 4.0 -o borra.g96

	vi -o borra.g96 cope.top

	gmx grompp -f cope.mdp -c borra.g96 -p cope.top -n cope.ndx -o borra.tpr -po borra.mdp
	"""

	def __init__( self ):
		self.exe  = "source /Users/smarti/Devel/gromacs/2019.b1/bin/GMXRC.bash;"
		self.exe += "gmx mdrun -s gromacs.tpr -e gromacs.edr -o gromacs.trr -g gromacs.log -rerun gromacs.g96 >& /dev/null"


	def update_coor( self, mol ):
		coordinates_write( mol, "gromacs.g96" )


#	def update_chrg( self, mol, fname ):
#	       ...NOT AVAILABLE...


	@staticmethod
	def __parse_edr( fname ):
		f = open( fname, "rb" )
		f.read( 8 )
		n = struct.unpack( ">i", f.read( 4 ) )[0]
		k = []
		for i in range( n ):
			t = struct.unpack( ">i", f.read( 4 ) )[0]
			if( t%4 != 0 ):
				k.append( f.read( 4 * ( ( t // 4 ) + 1 ) )[0:t] )
			else:
				k.append( f.read( t ) )
			t = struct.unpack( ">i", f.read( 4 ) )[0]
			if( t%4 != 0 ):
				t = 4 * ( ( t // 4 ) + 1 )
			f.read( t )
		f.seek( os.stat( fname )[stat.ST_SIZE] - n * 4 )
		return( dict( zip( k, struct.unpack( ">%df"%( n ), f.read( n * 4 ) ) ) ) )


	@staticmethod
	def __parse_trr( fname ):
		f = open( fname, "rb" )
		f.read( 32 )
		box_siz = struct.unpack( ">i", f.read( 4 ) )[0]
		vir_siz = struct.unpack( ">i", f.read( 4 ) )[0]
		prs_siz = struct.unpack( ">i", f.read( 4 ) )[0]
		f.read( 8 )
		crd_siz = struct.unpack( ">i", f.read( 4 ) )[0]
		vel_siz = struct.unpack( ">i", f.read( 4 ) )[0]
		frz_siz = struct.unpack( ">i", f.read( 4 ) )[0]
		natoms  = struct.unpack( ">i", f.read( 4 ) )[0]
		f.read( 16 )
		f.read( box_siz + vir_siz + prs_siz + crd_siz + vel_siz )
#		g = [ - i / 10.0 for i in struct.unpack( ">%df"%( 3 * natoms ), f.read( frz_siz ) ) ]
		g = []
		for i in range( 3 * natoms ):
			g.append( - struct.unpack( ">f", f.read( 4 ) )[0] / 10.0 )
		f.close()
		return( g )


	def get_func( self, mol ):
		self.update_coor( mol )
		try:
			os.unlink( "gromacs.log" )
			os.unlink( "gromacs.edr" )
			os.unlink( "gromacs.trr" )
		except:
			pass
		os.system( self.exe )
		mol.func += gromacs.__parse_edr( "gromacs.edr" )["Potential"]


	def get_grad( self, mol ):
		self.get_func( mol )
		g = gromacs.__parse_trr( "gromacs.trr" )
		for i in range( mol.natm ):
			i3 = i * 3
			for j in [0, 1, 2]:
				mol.grad[i3+j] += g[i3+j]








# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange
import	struct
import	math
import	qm3.mol
import	qm3.constants
import	qm3.io
import	os, stat
import	time



def coordinates_read( fname = None ):
	"""
* TITLE
*  DATE:    11/ 8/18     12:47:12      CREATED BY USER: smarti
*
 3092
    1    1 COPE C1    -0.98294   1.05507   0.12061 ACS  1      0.00000
           |||| ||||----------||||||||||---------- |||| ||||
.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789
                      ||||||||  ||||||||--------------------||||||||||||||||||||--------------------  ||||||||  ||||||||
         1         1  COPE      C1             -0.9830000000        1.0550000000        0.1210000000  ACS       1               0.0000000000
	"""
	mol = qm3.mol.molecule()
	f = qm3.io.open_r( fname )
	l = f.readline()
	while( l[0] == "*" ):
		l = f.readline()
	if( l.strip().split()[-1].upper() == "EXT" ):
		mol.natm = int( l.strip().split()[0] )
		for i in range( mol.natm ):
			l = f.readline()
			mol.resn.append( l[22:30].strip() )
			mol.labl.append( l[32:40].strip() )
			mol.coor.append( float( l[40:60].strip() ) )
			mol.coor.append( float( l[60:80].strip() ) )
			mol.coor.append( float( l[80:100].strip() ) )
			mol.segn.append( l[102:110].strip() )
			mol.resi.append( int( l[112:120] .strip() ) )
	else:
		mol.natm = int( l.strip() )
		for i in range( mol.natm ):
			l = f.readline()
			mol.resn.append( l[11:15].strip() )
			mol.labl.append( l[16:20].strip() )
			mol.coor.append( float( l[20:30].strip() ) )
			mol.coor.append( float( l[30:40].strip() ) )
			mol.coor.append( float( l[40:50].strip() ) )
			mol.segn.append( l[51:55].strip() )
			mol.resi.append( int( l[56:60] .strip() ) )
	qm3.io.close( f, fname )
	mol.settle()
	return( mol )
	

def coordinates_write( mol, fname = None ):
	f = qm3.io.open_w( fname )
	fmt = "%5d%5d %-4s %-4s%10.5lf%10.5lf%10.5lf %-4s %-4d%10.5lf\n"
	flg = ( mol.natm >= 100000 )
	if( not flg ):
		i = 0
		while( i < mol.natm and not flg ):
			flg |= ( len( mol.labl[i] ) > 4 ) or ( len( mol.resn[i] ) > 4 ) or ( len( mol.segn[i] ) > 4 )
			i += 1
	if( flg ):
		fmt = "%10d%10d  %-8s  %-8s%20.10lf%20.10lf%20.10lf  %-8s  %-8d%20.10lf\n"
		f.write( "* QMCube\n*\n%10d  EXT\n"%( mol.natm ) )
	else:
		f.write( "* QMCube\n*\n%5d\n"%( mol.natm ) )
	k = 1
	for i in range( len( mol.res_lim ) - 1 ):
		for j in range( mol.res_lim[i], mol.res_lim[i+1] ):
			f.write( fmt%( k, i+1, mol.resn[j], mol.labl[j],mol.coor[3*j],
				mol.coor[3*j+1], mol.coor[3*j+2], mol.segn[j], mol.resi[j], 0.0 ) )
			k += 1
	qm3.io.close( f, fname )


def selection( mol, sele, fname = None ):
	f = qm3.io.open_w( fname )
	t = [ False for i in range( mol.natm ) ]
	for i in sele:
		t[i] = True
	k = 0
	c = 0
	f.write( "defi s%02d sele -\n"%( k ) )
	for i in range( len( mol.res_lim ) - 1 ):
		if( sum( t[mol.res_lim[i]:mol.res_lim[i+1]] ) == mol.res_lim[i+1] - mol.res_lim[i] ):
			f.write( "    ( segid %s .and. resi %d )"%( mol.segn[mol.res_lim[i]], mol.resi[mol.res_lim[i]] ) )
			c += 1
			if( c > 50 ):
				k += 1
				f.write( " -\nend\ndefi s%02d sele -\n"%( k ) )
				c = 0
			else:
				f.write( " .or. -\n" )
		else:
			for j in range( mol.res_lim[i], mol.res_lim[i+1] ):
				if( t[j] ):
					f.write( "    ( segid %s .and. resi %d .and. type %s )"%( mol.segn[j], mol.resi[j], mol.labl[j] ) )
					c += 1
					if( c > 50 ):
						k += 1
						f.write( " -\nend\ndefi s%02d sele -\n"%( k ) )
						c = 0
					else:
						f.write( " .or. -\n" )
	f.write( "end\n" )
	qm3.io.close( f, fname )



class charmm_pipe( object ):

	def __init__( self, charmm_path, fname, sleep = 10 ):
		os.mkfifo( "charmm.pipe" )
		os.system( charmm_path + " < charmm.pipe > charmm.log &" )
		self.pfd = open( "charmm.pipe", "wt" )
		f = open( fname, "rt" )
		self.pfd.write( f.read() )
		self.pfd.flush()
		f.close()
		time.sleep( sleep )
		self.pfd.write( "format (G20.14)\n" )
		self.pfd.flush()


	def stop( self ):
		self.pfd.write( "stop\n" )
		self.pfd.flush()
		self.pfd.close()


	def update_charges( self, mol, sele ):
		f = open( "charmm.chrg", "wb" )
		f.write( struct.pack( "i", 84 ) + b"COOR" )
		for i in [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]:
			f.write( struct.pack( "i", i ) )
		f.write( struct.pack( "i", 84 ) )
		f.write( struct.pack( "i", 164 ) )
		f.write( struct.pack( "i", 2 ) )
		f.write( b"* Updating selected charges...                                                  " )
		f.write( b"*                                                                               " )
		f.write( struct.pack( "i", 164 ) )
		n = len( sele )
		f.write( struct.pack( "i", 4 ) + struct.pack( "i", n ) + struct.pack( "i", 4 ) )
		for i in [0, 1, 2]:
			f.write( struct.pack( "i", n * 4 ) )
			for j in range( n ):
				f.write( struct.pack( "f", 0.0 ) )
			f.write( struct.pack( "i", n * 4 ) )
		f.write( struct.pack( "i", n * 4 ) )
		for i in sele:
			f.write( struct.pack( "f", mol.chrg[i] ) )
		f.write( struct.pack( "i", n * 4 ) )
		f.close()
		self.pfd.write( "open unit 99 read file name charmm.chrg\n" )
		self.pfd.write( "read coor comp sele core end file unit 99\n" )
		self.pfd.write( "close unit 99\n" )
		self.pfd.write( "scal char copy wmai sele core show end\n" )
		self.pfd.flush()


	def update_coor( self, mol ):
		f = open( "charmm.coor", "wb" )
		f.write( struct.pack( "i", 84 ) + b"COOR" )
		for i in [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]:
			f.write( struct.pack( "i", i ) )
		f.write( struct.pack( "i", 84 ) )
		f.write( struct.pack( "i", 164 ) )
		f.write( struct.pack( "i", 2 ) )
		f.write( b"* Updating coordinates...                                                       " )
		f.write( b"*                                                                               " )
		f.write( struct.pack( "i", 164 ) )
		f.write( struct.pack( "i", 4 ) + struct.pack( "i", mol.natm ) + struct.pack( "i", 4 ) )
		# here should come crystal information, but charmm doesn't seem to miss it... :d
		for i in [0, 1, 2]:
			f.write( struct.pack( "i", mol.natm * 4 ) )
			for j in range( mol.natm ):
				f.write( struct.pack( "f", mol.coor[3*j+i] ) )
			f.write( struct.pack( "i", mol.natm * 4 ) )
		f.write( struct.pack( "i", mol.natm * 4 ) )
		for i in range( mol.natm):
			f.write( struct.pack( "f", 0.0 ) )
		f.write( struct.pack( "i", mol.natm * 4 ) )
		f.close()
		self.pfd.write( "open unit 99 read file name charmm.coor\n" )
		self.pfd.write( "read coor file unit 99\n" )
		self.pfd.write( "close unit 99\n" )
		self.pfd.flush()


	def get_func( self, mol ):
		self.pfd.write( "energy\n" )
		self.pfd.write( "open unit 99 write card name charmm.ener\n" )
		self.pfd.write( "write ener unit 99\n" )
		self.pfd.write( "* ?ENER\n*\n" )
		self.pfd.write( "close unit 99\n" )
		self.pfd.flush()
		while( not os.access( "charmm.ener", os.R_OK ) ):
			time.sleep( 0.01 )
		while( os.stat( "charmm.ener" )[stat.ST_SIZE] < 78 ):
			time.sleep( 0.01 )
		f = open( "charmm.ener", "rt" )
		mol.func += float( f.readline().strip().split()[-1] ) * qm3.constants.K2J
		f.close()
		os.unlink( "charmm.ener" )


	def get_grad( self, mol ):
		self.pfd.write( "energy\n" )
		self.pfd.write( "coor force comp\n" )
		self.pfd.write( "open unit 99 write file name charmm.force\n" )
		self.pfd.write( "write coor comp unit 99\n" )
		self.pfd.write( "* ?ENER\n*\n" )
		self.pfd.write( "close unit 99\n" )
		self.pfd.flush()
		sz = ( 4 + 84 + 4 ) + ( 4 + 164 + 4 ) + ( 4 + 4 + 4 ) + 4 * ( 4 + 4 * mol.natm + 4 )
		while( not os.access( "charmm.force", os.R_OK ) ):
			time.sleep( 0.01 )
		while( os.stat( "charmm.force" )[stat.ST_SIZE] < sz ):
			time.sleep( 0.01 )
		f = open( "charmm.force", "rb" )
		f.read( 100 )
		mol.func += float( str( f.read( 160 ) ).split( "*" )[1].strip() ) * qm3.constants.K2J
		f.read( 16 )
		if( os.stat( "charmm.force" )[stat.ST_SIZE] > sz ):
			f.read( 56 )
		for i in [0, 1, 2]:
			f.read( 4 )
			for j in range( mol.natm ):
				mol.grad[3*j+i] += struct.unpack( "f", f.read( 4 ) )[0] * qm3.constants.K2J
			f.read( 4 )
		f.close()
		os.unlink( "charmm.force" )





CHARMM_INP = """
* title
*

prnl 6
wrnl 6
bomblvl -1

open read form unit 10 name top
read rtf card unit 10
close unit 10

open read form unit 10 name prm
read parameter card unit 10
close unit 10

!read sequence COPE 1
!generate ACS setup noangle nodihedral
!read sequence HOH 1025
!generate WAT setup noangle nodihedral
!open unit 10 write form name psf
!write psf card unit 10
!close unit 10

open unit 10 read form name psf
read psf card unit 10
close unit 10

open unit 10 read form name card
read coor card unit 10
close unit 10

faster on

set box 31.4
calc fft = int( @box ) + 1
set con 12.0
set cof 14.0
set cls 16.0

crys defi cubic @box @box @box 90.0 90.0 90.0
crys build cutoff @cof noper 0
imag fixed xcen 0.0 ycen 0.0 zcen 0.0 sele all end

defi core sele ( segi ACS .and. resi 1 ) show end
cons fix sele core end

!nbonds -
!	elec fswitch vdw vswitch -
!	cutnb @cls ctofnb @cof ctonnb @con

nbonds -
	atom ewald pmewal kappa 0.32 fftx @fft ffty @fft fftz @fft -
	order 6 cutim @cof cutnb @cls ctofnb @cof ctonnb @con
"""

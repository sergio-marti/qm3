#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-
import    sys
import    os
import    xnamd
import    qm3.mol
import    qm3.problem
import    qm3.constants
import    qm3.engines.sqm
import    qm3.engines._qmmm
import    qm3.actions.minimize
import    io


# ======================================================================
"""
export GauExe=`pwd`/pilgrim_qm3.py
export GauFchk=/usr/bin/true

rm -rf 1-GTS 2-PLG_DATA 3-PLG_OUTPUT 4-PLG_RST 5-MOLDEN 6-PLOTFILES TMP pif.* tracking
pilgrim.py --gather
pilgrim.py --input << EOD
add temp 300
add chem reac: r --> t --> p
add path t
sbw=-1.05
sfw=+1.05
ds=0.01
hsteps=10
..
end
EOD

cat > pif.calcs << EOD
start_meppoint t gaussian
[Pilgrim_gradhess] [Pilgrim_name]
[Pilgrim_geometry]           
end_meppoint
EOD

pilgrim.py --pfn  | tee log
pilgrim.py --path | tee -a log
pilgrim.py --plot
"""
# ======================================================================


# set the right folder...
os.chdir( "/Users/smarti/Devel/Pilgrim/tests/qm3" )


# -- define selections (avoid overlap and sort'em...)
sQM = [ 0, 1, 2, 3 ]
sMM = list( range( 4, 796 ) )

#sMM = list( set( sMM ).difference( set( sQM ) ) )
#sQM.sort()
#sMM.sort()


# -- parse Pilgrim input
cmd, dst = sys.stdin.readline().strip().split()
tmp = []
for l in sys.stdin:
    tmp += [ float( i ) for i in l.strip().split()[1:] ]


# -- define molecule
mol = qm3.mol.molecule( "setup/pdb" )
who = dst.split( "/" )[-1]
if( who == "bw1" or who == "fw1" ):
    print( ">> loading TS geometry: " + who )
    mol.xyz_read( "0xyz", replace = True )
else:
    mol.xyz_read( "xyz", replace = True )
mol.boxl = [ 20., 20., 20., ]
mol.psf_read( "setup/psf" )
mol.nbnd_read( "non_bonded" )
mol.guess_atomic_numbers()


# -- update QM atoms (>> mol)
for i in range( len( sQM ) ):
    for j in [0, 1, 2]:
        mol.coor[3*sQM[i]+j] = tmp[3*i+j]


# -- define QM engine
os.environ["QM3_LIBSQM"] = "/Users/smarti/Devel/amber/18/qm3/libsqm.so"
f = io.StringIO( """_slave_
&qmmm
maxcyc    = 0,
qm_theory = "PM6",
qmcharge  = 0,
qmmm_int  = 1,
verbosity = 4
 /
qm3_atoms
qm3_charges
""" )
f.seek( 0 )
#eQM = qm3.engines.sqm.sqm( mol, f, sQM, sMM, [] )
eQM = qm3.engines.sqm.dl_sqm( mol, f, sQM, sMM, [] )
f.close()


# -- define problem object (recycle: mol, eQM and sMM )
class my_problem( qm3.problem.template ):
    def __init__( self, mol, eQM, sMM ):
        qm3.problem.template.__init__( self )
        self.mol  = mol
        self.sel  = sMM
        self.eMM  = xnamd.namd()
        self.eQM  = eQM
        self.size = 3 * len( self.sel )
        self.coor = []
        for i in self.sel:
            self.coor += self.mol.coor[3*i:3*i+3]


    def update_coor( self ):
        for i in range( self.size // 3 ):
            for j in [0, 1, 2]:
                self.mol.coor[3*self.sel[i]+j] = self.coor[3*i+j]


    def get_func( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.eQM.get_func( self.mol )
        self.eMM.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( 3 * self.mol.natm ) ]
        self.eQM.get_grad( self.mol )
        self.eMM.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = []
        for i in self.sel:
            self.grad += self.mol.grad[3*i:3*i+3]

        
# -- minize environment (and store result for future...)
obj = my_problem( mol, eQM, sMM )
qm3.actions.minimize.fire( obj, print_frequency = 100, gradient_tolerance = 0.1, step_number = 2000 )
mol.xyz_write( "xyz" )


# -- satisfy Pilgrim...
mol.func = 0.0
mol.grad = [ 0.0 for i in range( 3 * mol.natm ) ]
eQM.get_grad( mol )
fQM = qm3.engines._qmmm.Int_QMLJ( mol, sQM, sMM, [] )
fQM.get_func( mol )
fQM.get_grad( mol )
size = len( sQM ) * 3
func = mol.func
grad = []
for i in sQM:
    grad += mol.grad[3*i:3*i+3]


# -- numerical hessian
if( cmd == "freq=noraman" ):
    dsp = 1.e-4
    tmp = []
    for i in sQM:
        for j in [0, 1, 2]:
            t = mol.coor[3*i+j]
            mol.coor[3*i+j] = t + dsp
            mol.grad = [ 0.0 for k in range( 3 * mol.natm ) ]
            eQM.get_grad( mol )
            fQM.get_grad( mol )
            gp = []
            for k in sQM:
                gp += mol.grad[3*k:3*k+3]
            mol.coor[3*i+j] = t - dsp
            mol.grad = [ 0.0 for k in range( 3 * mol.natm ) ]
            eQM.get_grad( mol )
            fQM.get_grad( mol )
            gm = []
            for k in sQM:
                gm += mol.grad[3*k:3*k+3]
            tmp.append( [ ( gp[k] - gm[k] ) / ( 2. * dsp ) for k in range( size ) ] )    
            mol.coor[3*i+j] = t
    mol.hess = []
    for i in range( size ):
        for j in range( size ):
            mol.hess.append( ( tmp[i][j] + tmp[j][i] ) * 0.5 )
    fQM.get_hess( mol )


# -- fake fchk
f = open( "%s.fchk"%( dst ), "wt" )
f.write( "Number of atoms                            I%17d\n"%( size // 3 ) )
f.write( "Charge                                     I%17d\n"%( 0 ) )
f.write( "Multiplicity                               I%17d\n"%( 1 ) )
f.write( "Total Energy                               R%27.15lE\n"%( func / qm3.constants.H2J ) )
f.write( "Atomic numbers                             I   N=%12d\n"%( size // 3 ) )
t = []
for i in range( size // 3 ):
    t.append( "%12d"%( mol.anum[sQM[i]] ) )
    if( (i+1)%6 == 0 ):
        t.append( "\n" )
f.write( "".join( t ) )
if( (size//3)%6 != 0 ):
    f.write( "\n" )
f.write( "Current cartesian coordinates              R   N=%12d\n"%( size ) )
t = []
k = 0
for i in range( size // 3 ):
    for j in [0, 1, 2]:
        t.append( "%16.8lE"%( mol.coor[sQM[i]*3+j] / qm3.constants.A0 ) )
        k += 1
        if( k%5 == 0 ):
            t.append( "\n" )
f.write( "".join( t ) )
if( size%5 != 0 ):
    f.write( "\n" )
f.write( "Real atomic weights                        R   N=%12d\n"%( size // 3 ) )
t = []
for i in range( size // 3 ):
    t.append( "%16.8lE"%( mol.mass[sQM[i]] ) )
    if( (i+1)%5 == 0 ):
        t.append( "\n" )
f.write( "".join( t ) )
if( (size//3)%5 != 0 ):
    f.write( "\n" )
f.write( "Cartesian Gradient                         R   N=%12d\n"%( size ) )
c = qm3.constants.A0 / qm3.constants.H2J
t = []
for i in range( size ):
    t.append( "%16.8lE"%( grad[i] * c ) )
    if( (i+1)%5 == 0 ):
        t.append( "\n" )
f.write( "".join( t ) )
if( size%5 != 0 ):
    f.write( "\n" )
if( cmd == "freq=noraman" ):
    n = size * ( size + 1 ) // 2
    f.write( "Cartesian Force Constants                  R   N=%12d\n"%( n ) )
    c = qm3.constants.A0 * qm3.constants.A0 / qm3.constants.H2J
    t = []
    k = 0
    for i in range( size ):
        for j in range( i + 1 ):
            t.append( "%16.8lE"%( mol.hess[i+size*j] * c ) )
            if( (k+1)%5 == 0 ):
                t.append( "\n" )
            k += 1
    f.write( "".join( t ) )
    if( n%5 != 0 ):
        f.write( "\n" )
f.close()

print( "normal termination" )

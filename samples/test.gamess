import qm3.mol
import qm3.problem
import qm3.engines.gamess
import qm3.maths.matrix
try:
  import cStringIO as io
except:
  import io


class my_eng( qm3.problem.template ):

    def __init__( self ):
        self.mol = qm3.mol.molecule()
        self.mol.xyz_read( "xyz" )
        self.mol.chrg = [ -0.834, 0.417, 0.417, -0.834, 0.417, 0.417 ]
        self.mol.anum = [ 8, 1, 1, 8, 1, 1 ]
        f = io.StringIO( """
 $contrl
runtyp=qm3_job
coord=cart
nosym=1
nprint=7
scftyp=rhf
units=angs
icharg=0
mult=1
maxit=200
dfttyp=b3lyp
 $end

 $scf
dirscf=.true.
conv=1.d-6
 $end

 $system
mwords=10
memddi=1024
 $end

 $basis
extfil=.true.
gbasis=AABBCCDD
 $end

 $data
c1
qm3_atoms
 $end

 $elpot
iepot=1
where=pdc
 $end

 $pdc
ptsel=geodesic
constr=charge
 $end

qm3_guess
""" )
        f.seek( 0 )
        self.eng = qm3.engines.gamess.gamess( self.mol, f, [ 0, 1, 2 ], [ 3, 4, 5 ] )
        self.eng.exe = "bash r.gamess"

        self.size = 3 * self.mol.natm
        self.coor = self.mol.coor
        self.func = 0.0
        self.grad = []
    

    def get_func( self ):
        self.mol.func = 0.0
        self.eng.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.size ) ]
        self.eng.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = self.mol.grad
    


f = open( "xyz", "wt" )
f.write( """6

O       0.12109      0.06944     -0.22458
H      -0.52694      0.16499     -0.94583
H       0.70159     -0.63565     -0.54677
O      -0.45114      1.12675      2.21102
H      -0.29157      0.59483      1.39876
H       0.05804      1.92714      2.01036
""" )
f.close()


# https://bse.pnl.gov/bse/portal
f = open( "basis", "wt" )
f.write( """H AABBCCDD
S   3
  1     13.0107010              0.19682158E-01   
  2      1.9622572              0.13796524       
  3      0.44453796             0.47831935       
S   1
  1      0.12194962             1.0000000        
P   1
  1      0.8000000              1.0000000        

O AABBCCDD
S   5
  1   2266.1767785             -0.53431809926E-02      
  2    340.87010191            -0.39890039230E-01      
  3     77.363135167           -0.17853911985    
  4     21.479644940           -0.46427684959    
  5      6.6589433124          -0.44309745172    
S   1
  1      0.80975975668          1.0000000        
S   1
  1      0.25530772234          1.0000000        
P   3
  1     17.721504317            0.43394573193E-01      
  2      3.8635505440           0.23094120765    
  3      1.0480920883           0.51375311064    
P   1
  1      0.27641544411          1.0000000        
D   1
  1      1.2000000              1.0000000
""" )
f.close()


f = open( "r.gamess", "wt" )
f.write( """GWD=./bin
SCR=`pwd`
DDI=$GWD/ddikick.x
EXE=$GWD/gamess.00.x
JOB=gamess.temp

export   INPUT=$SCR/gamess.inp
export  EXTBAS=basis
export AUXDATA=$GWD/auxdata
export ERICFMT=$AUXDATA/ericfmt.dat
export BASPATH=$AUXDATA/BASES
export  DASORT=$SCR/$JOB.F20
export   PUNCH=$SCR/gamess.data
export DICTNRY=$SCR/gamess.guess
export DFTGRID=$SCR/$JOB.F22

ulimit -c 0
rm -f $JOB.* $PUNCH
export DDI_RSH=ssh
export DDI_VER=new
export NNODES=1
export NCPUS=1
export HOSTLIST="`hostname`:cpus=$NCPUS"
$DDI $EXE $JOB -ddi $NNODES $NCPUS $HOSTLIST -scr $SCR < /dev/null >& gamess.out
""" )
f.close()

obj = my_eng()
obj.get_grad()
print( obj.func )
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )

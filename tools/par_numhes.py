#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
from __future__ import print_function, division
import  sys
if( sys.version_info[0] == 2 ):
    range = xrange

import  os
import  time

try:
    import cPickle as pickle
except:
    import pickle
try:
    import cStringIO as io
except:
    import io

import  qm3.mol
import  qm3.problem
import  qm3.engines.sqm
import  qm3.utils
import  qm3.utils._mpi


os.environ["AMBERHOME"] = "/Users/smarti/Devel/amber/18"
os.environ["QM3_LIBSQM"] = "/Users/smarti/Devel/amber/18/qm3/libsqm.so"


class my_problem( qm3.problem.template ):
    def __init__( self, node ):
        qm3.problem.template.__init__( self )

        try:
            os.mkdir( str( node ) )
        except:
            pass
        try:
            os.chdir( str( node ) )
        except:
            sys.exit( 1 )

        self.mol = qm3.mol.molecule()
        f = io.StringIO( """27

Si         10.6199068403        7.9323943986        8.7879816551
O           9.0734254523        7.3970643400        9.1540684162
O          10.6177609394        9.5696549634        9.0141306065
Si          9.4401539974       10.7086517093        9.2002571964
Si          7.6727218048        8.2660128921        9.1891363626
O           8.0640397314        9.8341895057        9.5334914065
O           9.2522735204       11.5409246246        7.7857485914
O          10.8007593346        7.5945973616        7.1857167985
O           6.9913323546        8.2178917047        7.6814690520
Si         10.2568996490        8.2676141922        5.7924984688
Si          7.3563651278        8.6659954997        6.1596098990
O           8.6524762731        7.8106950351        5.6099664753
O          10.2812221248        9.8935632439        5.9625090567
O           7.6704150853       10.2677915282        6.0427041527
Si          9.1392540558       11.0542165469        6.2175778533
O          11.2135513663        7.8984798620        4.5364151126
O           9.7067004971       11.7439705063       10.4241556893
O          11.7769178601        7.1259807877        9.5767943035
H          11.5629083312        6.9052025139       10.5156864008
O           6.6852146939        7.7485558757       10.3659838673
H           5.7798037170        8.1544495512       10.4082897510
O           4.5946860366       10.2056912618        7.4547464206
H           5.3004055631        9.7433781200        7.9231302939
O           9.2330905755       12.2979176587        5.1992514639
H           9.9100534235       12.2046467779        4.4868343017
H          11.1361268091        7.0285123854        4.0736203370
H           8.9274672451       12.3062547556       10.6221200055
""" )
        self.mol.xyz_read( f )
        f.close()
        self.mol.guess_atomic_numbers()

        f = io.StringIO( """
_slave_
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
        self.eqm = qm3.engines.sqm.dl_sqm( self.mol, f, list( range( self.mol.natm ) ) )
        f.close()

        self.size = 3 * self.mol.natm
        self.coor = self.mol.coor
        self.mass = self.mol.mass
        self.hess = []


    def hess_step( self, who, dsp = 1.e-4 ):
        self.mol.func = 0.0
        bak = self.mol.coor[who]
        self.mol.coor[who] = bak + dsp
        self.mol.grad = [ 0.0 for i in range( self.size ) ]
        self.eqm.get_grad( self.mol )
        gf = self.mol.grad[:]
        self.mol.coor[who] = bak - dsp
        self.mol.grad = [ 0.0 for i in range( self.size ) ]
        self.eqm.get_grad( self.mol )
        self.mol.coor[who] = bak
        self.hess.append( [ ( i - j ) / ( 2.0 * dsp ) for i,j in zip( gf, self.mol.grad ) ] )




node, ncpu = qm3.utils._mpi.init()
obj = my_problem( node )

t = time.time()

for ite in range( obj.size ):
    if( ite % ncpu == node ):
        obj.hess_step( ite )

qm3.utils._mpi.barrier()
print( node, len( obj.hess ) )

if( node == 0 ):
    tmp = []
    cur = 0
    for ite in range( obj.size ):
        who = ite % ncpu
        if( who == 0 ):
            tmp.append( obj.hess[cur] )
            cur += 1
        else:
            tmp.append( qm3.utils._mpi.recv_r8( who, obj.size ) )
    hes = []
    for i in range( obj.size ):
        for j in range( obj.size ):
            hes.append( 0.5 * ( tmp[i][j] + tmp[j][i] ) )
    print( time.time() - t )
    val, vec = qm3.utils.hessian_frequencies( obj.mass, obj.coor, hes )
    print( val )
    f = open( "../parall.dump", "wb" )
    pickle.dump( hes, f )
    f.close()
else:
    for itm in obj.hess:
        qm3.utils._mpi.send_r8( 0, itm )

qm3.utils._mpi.barrier()
qm3.utils._mpi.stop()

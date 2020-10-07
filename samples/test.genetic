#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import    sys
if( sys.version_info[0] == 2 ):
    range = xrange

import    math
import    qm3.maths.rand
import    qm3.actions.genetic
import    qm3.actions.fitting
import    time
import    os
import    sys
try:
    import    cPickle as pickle
except:
    import    pickle

try:
    import    qm3.utils._mpi
    node, ncpu = qm3.utils._mpi.init()
except:
    node, ncpu = ( 0, 1 )


n = 100
if( node == 0 ):
    if( not os.access( "data", os.R_OK ) ):
        x0 = 0.0
        xf = 10.0
        d = ( xf - x0 ) / ( n - 1 )
        X = [ x0 + i * d for i in range( n ) ]
        Y = [ math.cos( X[i] ) + qm3.maths.rand.gauss( 0.0, 0.2 ) for i in range( n ) ]
        r, C, t = qm3.actions.fitting.poly_fit( X, Y, 5 )
        f = open( "data", "wb" )
        pickle.dump( X, f )
        pickle.dump( Y, f )
        pickle.dump( C, f )
        f.close()
    else:
        f = open( "data", "rb" )
        X = pickle.load( f )
        Y = pickle.load( f )
        C = pickle.load( f )
        f.close()

if( ncpu > 1 ):
    if( node == 0 ):
        for i in range( 1, ncpu ):
            qm3.utils._mpi.send_r8( i, X )
            qm3.utils._mpi.send_r8( i, Y )
    else:
        X = qm3.utils._mpi.recv_r8( 0, n )
        Y = qm3.utils._mpi.recv_r8( 0, n )


class p5_fitting:
    def __init__( self, x, y ):
        self.n = len( x )
        self.x = x
        self.y = y
        self.size = 6
        self.coor = []
        self.func = 0.0
        self.grad = []


    def poly( self, w ):
        return( sum( [ self.coor[i] * math.pow( self.x[w], i ) for i in range( self.size ) ] ) )


    def get_func( self ):
        self.func = 0.0
        for i in range( self.n ):
            self.func += math.pow( self.y[i] - self.poly( i ), 2.0 )


    def get_grad( self ):
        self.get_func()
        fbak = self.func
        self.grad = []
        for i in range( self.size ):
            tmp = self.coor[i]
            self.coor[i] = tmp + 1.0e-6
            self.get_func()
            f_for = self.func
            self.coor[i] = tmp - 1.0e-6
            self.get_func()
            f_bak = self.func
            self.coor[i] = tmp
            self.grad.append( ( f_for - f_bak ) / 2.0e-6 )
        self.func = fbak


    def current_step( self, istep ):
        pass


obj = p5_fitting( X, Y )

t0 = time.time()
if( ncpu > 1 ):
    f = open( "log.mpi", "at" )
    qm3.actions.genetic.mpi_diffevo( obj,
        [ ( -5, 5 ) ] * 6,
        mpi_node = node, mpi_ncpu = ncpu, 
        population_size = 10 )
else:
    if( len( sys.argv ) == 3 ):
       f = open( "log.smp", "at" )
       a = [ obj ]
       for i in range( 1, int( sys.argv[1] ) ):
           a.append( p5_fitting( X, Y ) )
       qm3.actions.genetic.smp_diffevo( a,
           [ ( -5, 5 ) ] * 6,
           population_size = int( sys.argv[2] ) )
    else:
        f = open( "log.serial", "at" )
        qm3.actions.genetic.diffevo( obj,
            [ ( -5, 5 ) ] * 6,
            population_size = 40 )
t0 = time.time() - t0

if( node == 0 ):
    f.write( "%12.6lf%20.3lf%20.10lf\n"%( t0, obj.func, math.sqrt( sum( [ math.pow( obj.coor[i] - C[i], 2.0 ) for i in range( obj.size ) ] ) ) ) )
    f.close()

if( ncpu > 1 ):
    qm3.utils._mpi.barrier()
    qm3.utils._mpi.stop()

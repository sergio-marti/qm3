# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.problem



class muller_brown( qm3.problem.template ):

    def __init__( self ):
        self.A  = [ -200.0, -100.0, -170.0, 15.0 ]
        self.a  = [ -1.0, -1.0, -6.5, 0.7 ]
        self.b  = [ 0.0, 0.0, 11.0, 0.6 ]
        self.c  = [ -10.0, -10.0, -6.5, 0.7 ]
        self.xo = [ 1.0, 0.0, -0.5, -1.0 ]
        self.yo = [ 0.0, 0.5, 1.5, 1.0 ]
        self.mass = [ 1.0, 1.0 ]
        self.size = 2
        self.coor = [ 0.0, 0.0 ]
        self.func = 0.0
        self.grad = []
        self.hess = []


    def get_func( self ):
        self.func = 0.0
        for i in range( 4 ):
            self.func += self.A[i] * math.exp( self.a[i] * math.pow( self.coor[0] - self.xo[i], 2.0 ) + self.b[i] * ( self.coor[0] - self.xo[i] ) * ( self.coor[1] - self.yo[i] ) + self.c[i] * math.pow( self.coor[1] - self.yo[i], 2.0 ) )


    def get_grad( self ):
        self.func = 0.0
        self.grad = [ 0.0, 0.0 ]
        for i in range( 4 ):
            f = self.A[i] * math.exp( self.a[i] * math.pow( self.coor[0] - self.xo[i], 2.0 ) + self.b[i] * ( self.coor[0] - self.xo[i] ) * ( self.coor[1] - self.yo[i] ) + self.c[i] * math.pow( self.coor[1] - self.yo[i], 2.0 ) )
            self.func    += f
            self.grad[0] += f * ( 2.0 * self.a[i] * ( self.coor[0] - self.xo[i] ) + self.b[i] * ( self.coor[1] - self.yo[i] ) )
            self.grad[1] += f * ( self.b[i] * ( self.coor[0] - self.xo[i] ) + 2.0 * self.c[i] * ( self.coor[1] - self.yo[i] ) )


    def get_hess( self ):
        self.func = 0.0
        self.grad = [ 0.0, 0.0 ]
        self.hess = [ 0.0, 0.0, 0.0, 0.0 ]
        for i in range( 4 ):
            f = self.A[i] * math.exp( self.a[i] * math.pow( self.coor[0] - self.xo[i], 2.0 ) + self.b[i] * ( self.coor[0] - self.xo[i] ) * ( self.coor[1] - self.yo[i] ) + self.c[i] * math.pow( self.coor[1] - self.yo[i], 2.0 ) )
            self.func    += f
            self.grad[0] += f * ( 2.0 * self.a[i] * ( self.coor[0] - self.xo[i] ) + self.b[i] * ( self.coor[1] - self.yo[i] ) )
            self.grad[1] += f * ( self.b[i] * ( self.coor[0] - self.xo[i] ) + 2.0 * self.c[i] * ( self.coor[1] - self.yo[i] ) )
            self.hess[0] += f * ( 2.0 * self.a[i] + math.pow( 2.0 * self.a[i] * ( self.coor[0] - self.xo[i] ) + self.b[i] * ( self.coor[1] - self.yo[i] ), 2.0 ) )
            self.hess[3] += f * ( 2.0 * self.c[i] + math.pow( 2.0 * self.c[i] * ( self.coor[1] - self.yo[i] ) + self.b[i] * ( self.coor[0] - self.xo[i] ), 2.0 ) )
            t = f * ( self.b[i] + ( 2.0 * self.a[i] * ( self.coor[0] - self.xo[i] ) + self.b[i] *( self.coor[1] - self.yo[i] ) ) * ( self.b[i] * ( self.coor[0] - self.xo[i] ) + 2.0 * self.c[i] * ( self.coor[1] - self.yo[i] ) ) )
            self.hess[1] += t
            self.hess[2] += t



class cerjan_miller( qm3.problem.template ):

    def __init__( self ):
        self.a = 1.0
        self.b = 1.2
        self.c = 1.0
        self.mass = [ 1.0, 1.0 ]
        self.size = 2
        self.coor = [ 0.0, 0.0 ]
        self.func = 0.0
        self.grad = []
        self.hess = []


    def get_func( self ):
        x2 = self.coor[0] * self.coor[0]
        gx = math.exp( - x2 )
        y2 = self.coor[1] * self.coor[1]
        self.func = ( self.a - self.b * y2 ) * x2 * gx + self.c * 0.5 * y2


    def get_grad( self ):
        x2 = self.coor[0] * self.coor[0]
        gx = math.exp( - x2 )
        y2 = self.coor[1] * self.coor[1]
        self.func = ( self.a - self.b * y2 ) * x2 * gx + self.c * 0.5 * y2
        self.grad = [ 0.0, 0.0 ]
        self.grad[0] = - 2.0 * gx * self.coor[0] * ( x2 - 1.0 ) * ( self.a - self.b * y2 )
        self.grad[1] = self.coor[1] * ( self.c - 2.0 * self.b * x2 * gx )


    def get_hess( self ):
        x2 = self.coor[0] * self.coor[0]
        gx = math.exp( - x2 )
        y2 = self.coor[1] * self.coor[1]
        self.func = ( self.a - self.b * y2 ) * x2 * gx + self.c * 0.5 * y2
        self.grad = [ 0.0, 0.0 ]
        self.grad[0] = - 2.0 * gx * self.coor[0] * ( x2 - 1.0 ) * ( self.a - self.b * y2 )
        self.grad[1] = self.coor[1] * ( self.c - 2.0 * self.b * x2 * gx )
        self.hess = [ 0.0, 0.0, 0.0, 0.0 ]
        self.hess[0] = 2.0 * gx * ( 1.0 - 5.0 * x2 + 2.0 * x2 * x2 ) * ( self.a - self.b * y2 )
        self.hess[1] = 4.0 * self.b * gx * self.coor[0] * self.coor[1] * ( x2 - 1.0 )
        self.hess[2] = self.hess[1]
        self.hess[3] = self.c - 2.0 * self.b * gx * x2 





# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange



class template( object ):
    """
    Interface class for generic problem (energy/gradient/hessian)
    """
    def __init__( self ):
        self.size = 0
        self.func = 0.0
        self.coor = []
        self.grad = []
        self.hess = []
        self.mass = []


    def current_step( self, istep ):
        pass


    def get_func( self ):
        pass


    def get_grad( self ):
        pass


    def get_hess( self ):
        pass


    # >> paralelizar esto empleando subdirectorios (para evitar problemas al usar programas externos...)
    def num_grad( self, dsp = 1.e-4, central = True ):
        self.grad = []
        if( central ):
            for i in range( self.size ):
                t = self.coor[i]
                self.coor[i] = t + dsp
                self.get_func()
                fp = self.func
                self.coor[i] = t - dsp
                self.get_func()
                self.grad.append( ( fp - self.func ) / ( 2. * dsp ) )
                self.coor[i] = t
            self.get_func()
        else:
            self.get_func()
            f_bak = self.func
            for i in range( self.size ):
                t = self.coor[i]
                self.coor[i] = t + dsp
                self.get_func()
                self.grad.append( ( self.func - f_bak ) / dsp  )
                self.coor[i] = t
            self.func = f_bak


    def num_hess( self, dsp = 1.e-4, central = True ):
        hess = []
        if( central ):
            for i in range( self.size ):
                t = self.coor[i]
                self.coor[i] = t + dsp
                self.get_grad()
                gp = self.grad[:]
                self.coor[i] = t - dsp
                self.get_grad()
                hess.append( [ ( gp[j] - self.grad[j] ) / ( 2. * dsp ) for j in range( self.size ) ] )
                self.coor[i] = t
            self.get_grad()
        else:
            self.get_grad()
            f_bak = self.func
            g_bak = self.grad[:]
            for i in range( self.size ):
                t = self.coor[i]
                self.coor[i] = t + dsp
                self.get_grad()
                hess.append( [ ( self.grad[j] - g_bak[j] ) / dsp for j in range( self.size ) ] )
                self.coor[i] = t
            self.func = f_bak
            self.grad = g_bak[:]
        # symmetrize
        self.hess = []
        for i in range( self.size ):
            for j in range( self.size ):
                self.hess.append( ( hess[i][j] + hess[j][i] ) * 0.5 )

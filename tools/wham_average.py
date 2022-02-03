#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-
import    sys
import    math
import    qm3.constants
import    qm3.maths.rand



class wham( object ):

    def __init__( self, data ):
        self.__dd = []
        self.__dv = []
        self.__nw = 0
        self.__fc = []
        self.__rf = []
        self.__NN = []
        self.__nb = None
        self.__nn = []
        self.__nv = []
        self.__mx = None
        self.__Mx = None
        self.crd  = []
        self.pmf  = []
        self.var  = []
        # ---------------------------------------------
        for fn in data:
            f = open( fn, "rt" )
            t = [ float( i ) for i in f.readline().strip().split() ]
            self.__fc.append( t[0] )
            self.__rf.append( t[1] )
            if( not self.__mx ): self.__mx = t[1]
            if( not self.__Mx ): self.__Mx = t[1]
            n = 0.0
            for l in f:
                t = [ float( i ) for i in l.strip().split() ]
                self.__dd.append( t[0] )
                self.__dv.append( t[1] )
                self.__mx = min( self.__mx, t[0] )
                self.__Mx = max( self.__Mx, t[0] )
                n += 1.0
            f.close()
            self.__NN.append( n )
            self.__nw += 1


    def setup( self, nbins = None ):
        self.__nb = nbins
        if( not self.__nb ):
            self.__nb = 2 * self.__nw
        __db      = ( self.__Mx - self.__mx ) / float( self.__nb )
        self.__nn = [ 0.0 for i in range( self.__nb ) ]
        self.__nv = [ 0.0 for i in range( self.__nb ) ]
        self.crd  = [ self.__mx + __db * ( float( i ) + 0.5 ) for i in range( self.__nb ) ]
        for j in range( len( self.__dd ) ):
            w = sorted( [ ( math.fabs( self.crd[i] - self.__dd[j] ), i ) for i in range( self.__nb ) ] )[0][1]
            self.__nn[w] += 1.0
            self.__nv[w] += self.__dv[j]
        for j in range( self.__nb ):
            self.__nv[j] /= self.__nn[j]
            print( self.__nn[j], self.crd[j], self.__nv[j] )


    def integrate( self, temperature = 300.0, maxit = 20000, toler = 1.0e-3 ):
        __rt     = temperature * 1.0e-3 * qm3.constants.R
        __ff     = [ 0.0 for i in range( self.__nw ) ]
        __rr     = [ 0.0 for i in range( self.__nb ) ]
        __uu     = []
        __eu     = []
        self.pmf = []
        for i in range( self.__nb ):
            __uu.append( [] )
            __eu.append( [] )
            for j in range( self.__nw ):
                x = 0.5 * self.__fc[j] * math.pow( self.crd[i] - self.__rf[j], 2.0 )
                __uu[-1].append( x )
                __eu[-1].append( math.exp( - x / __rt ) )
        f = False
        I = 0
        while( I < maxit and not f ):
            __fo = __ff[:]
            for i in range( self.__nb ):
                __rr[i] = self.__nn[i] / sum( [ self.__NN[j] * math.exp( - ( __uu[i][j] - __ff[j] ) / __rt ) for j in range( self.__nw ) ] )
            for j in range( self.__nw ):
                __ff[j] = - __rt * math.log( sum( [ __eu[i][j] * __rr[i] for i in range( self.__nb ) ] ) )
            f = max( [ math.fabs( __fo[i] - __ff[i] ) for i in range( self.__nw ) ] ) < toler
            I += 1
        print( "#%9s%20d"%( f, I ) )
        if( f ):
            x = sum( __rr )
            for i in range( self.__nb ):
                self.var.append( self.__nv[i] * sum( [ math.exp( - ( __uu[i][j] - __ff[j] ) / __rt ) for j in range( self.__nw ) ] ) * __rr[i] )
                __rr[i] /= x
                if( __rr[i] > 0.0 ):
                    self.pmf.append( - __rt * math.log( __rr[i] ) )
                else:
                    self.pmf.append( 0.0 )
            x = max( self.pmf )
            print( "%20s%20s%20s"%( "Coordinate", "PMF", "<Variable>" ) )
            for i in range( self.__nb ):
                self.pmf[i] -= x
                print( "%20.10lf%20.10lf%20.10lf"%( self.crd[i], self.pmf[i], self.var[i] ) )
        return( f )




x = wham( sys.argv[1:] )
x.setup()
x.integrate()

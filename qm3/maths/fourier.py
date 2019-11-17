# -*- coding: iso-8859-1 -*-
#
#    http://en.wikipedia.org/wiki/Fourier_series
#
"""
http://rosettacode.org/wiki/Fast_Fourier_transform
http://en.wikipedia.org/wiki/Cooley–Tukey_FFT_algorithm
http://en.wikipedia.org/wiki/Discrete_Fourier_transform
http://es.wikipedia.org/wiki/Transformada_de_Fourier_discreta
"""

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.maths.integration
import qm3.maths.interpolation



class series( object ):

    def __init__( self, size, x, y, interpolant = qm3.maths.interpolation.hermite_spline ):
        self.size = size
        self.base = 0.0
        self.bcof = []
        self.acof = []
        self.x    = x[:]
        self.y    = y[:]
        self.peri = 0.5 * ( self.x[-1] - self.x[0] )
        self.inte = interpolant( x, y )


    def integrate( self ):
        self.base = 0.5 * qm3.maths.integration.Simpson_f( lambda r: self.inte.calc( r )[0], self.x[0], self.x[-1], n = 1000 )[0]
        for i in range( self.size ):
            f = float( i + 1 ) * math.pi / self.peri
            self.acof.append( qm3.maths.integration.Simpson_f( lambda r: self.inte.calc( r )[0] * math.cos( r * f ), self.x[0], self.x[-1], n = 1000 )[0] )
            self.bcof.append( qm3.maths.integration.Simpson_f( lambda r: self.inte.calc( r )[0] * math.sin( r * f ), self.x[0], self.x[-1], n = 1000 )[0] )


    def calc( self, r, items = None ):
        o = self.base
        n = self.size
        if( items ):
            if( items > 0 and items <= self.size ):
                n = items
        for i in range( n ):
            f = float( i + 1 ) * math.pi / self.peri
            o += self.acof[i] * math.cos( f * r ) + self.bcof[i] * math.sin( f * r )
        return( o / self.peri )





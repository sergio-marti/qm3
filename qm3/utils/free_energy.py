# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.constants
import qm3.maths.ode
import qm3.maths.stats
import qm3.maths.rand
import time



#
# Phys. Rev. Lett. v91, p140601 (2003) [10.1103/PhysRevLett.91.140601]
# J. Comp. Phys. v22, p245 (1976) [10.1016/0021-9991(76)90078-4]
# https://github.com/chqm3.maths.oderalab/pymbar.git
#
def bennett_acceptance_ratio( forward, reverse, temperature = 300. ):
    def __diff( size, exp_f, exp_r, dF ):
        edf = math.exp( - dF )
        f_f = 0.0
        f_r = 0.0
        for i in range( size ):
            f_f += 1.0 / ( 1.0 + edf / exp_f[i] )
            f_r += 1.0 / ( 1.0 + 1.0 / ( exp_r[i] * edf ) )
        return( math.log( f_f ) - math.log( f_r ) )

    rt = 1.0e-3 * qm3.constants.R * temperature
    exp_f = []
    exp_r = []
    size  = min( len( forward ), len( reverse ) )
    for i in range( size ):
        exp_f.append( math.exp( - forward[i] / rt ) )
        exp_r.append( math.exp( - reverse[i] / rt ) )
    x0 = - math.log( sum( exp_f ) / size )
    xf =   math.log( sum( exp_r ) / size )
    f0 = __diff( size, exp_f, exp_r, x0 )
    ff = __diff( size, exp_f, exp_r, xf )
    while( f0 * ff > 0.0 ):
        x0 -= 0.1
        f0 = __diff( size, exp_f, exp_r, x0 )
        xf += 0.1
        ff = __diff( size, exp_f, exp_r, xf )
    fm = 1.0
    i  = 0
    ii = 1000
    while( i < ii and math.fabs( fm ) > 1.e-8 ):
        xm = ( x0 + xf ) * 0.5
        fm = __diff( size, exp_f, exp_r, xm )
        if( f0 * fm <= 0.0 ):
            xf = xm
            ff = fm
        elif( ff * fm <= 0.0 ):
            x0 = xm
            f0 = fm
        else:
            return( None )
        i += 1
    if( i >= ii ):
        return( None )
    edf = math.exp( - xm )
    af  = 0.0
    af2 = 0.0
    ar  = 0.0
    ar2 = 0.0
    for i in range( size ):
        t    = 1.0 / ( 1.0 + edf / exp_f[i] )
        af  += t
        af2 += t * t
        t    = 1.0 / ( 1.0 + 1.0 / ( exp_r[i] * edf ) )
        ar  += t
        ar2 += t * t
    af  = af  / size
    ar  = ar  / size
    af2 = af2 / size
    ar2 = ar2 / size
    return( { "dF": rt * xm, "Error": rt * math.sqrt( ( af2 / ( af * af ) + ar2 / ( ar * ar ) - 2.0 ) / size ) } )



#
# First order expansion + k_means clustering
# ----------------------------------------------------------------------------------------------
# Chipot, C. Free Energy Calculations in Biological Systems. How Useful Are They in Practice?
#        New Algorithms for Macromolecular Simulation, 2006, Springer) 3-540-25542-7
# Straatsma, T.P.; Berendsen, H.J.C.; Stam, A.J. Molecular Physics, 1986, 57, 89-95 [10.1080/00268978600100071]
# Smith, E.B.; Wells, B.H. Molecular Physics, 1984, 52, 701-704 [10.1080/00268978400101481]
# ----------------------------------------------------------------------------------------------
#%DELTA F = - RT ln left langle e^{-{{U_j - U_i} over RT}} right rangle_i  +- RT {  { %delta %varepsilon } over { left langle e^{-{{U_j - U_i} over RT}} right rangle_i} }
#~~~~~~~~
#%delta %varepsilon^2 = { { 1 + 2 %tau } over N } left( left langle e^{-2{{U_j - U_i} over RT}} right rangle_i - left langle e^{-{{U_j - U_i} over RT}} right rangle_i^2 right)
#~~~~~~~~
#1 + 2 %tau = { 1 + r_1 } over { 1 - r_1 } 
#~~~~~~~~
#r_1 = { sum_{i=2}^N{left( x_i - bar x right) left( x_{i-1} - bar x right) } } over { sum_{i=1}^N{left( x_i - bar x right) ^2} }
# ----------------------------------------------------------------------------------------------
#
def fep_integrate( dene, temperature = 300.0, clusters = 1, tries = 10 ):
    def __integrate( dene, temperature ):
        rt  = temperature * 1.0e-3 * qm3.constants.R
        nn  = len( dene )
        exp = [ math.exp( - dene[i] / rt ) for i in range( nn ) ]
        mm  = sum( exp ) / float( nn )
        ex2 = sum( [ math.exp( - 2.0 * dene[i] / rt ) for i in range( nn ) ] ) / float( nn )
        try:
            r1  = sum( [ ( exp[i] - mm ) * ( exp[i-1] - mm ) for i in range( 1, nn ) ] )
            r1 /= sum( [ math.pow( exp[i] - mm, 2.0 ) for i in range( nn ) ] )
        except:
            r1 = 0.0
        sr  = ( 1.0 + r1 ) / ( 1.0 - r1 )
        err = rt * math.sqrt( sr / float( nn ) * math.fabs( ex2 - mm * mm ) ) / mm
        out = { "Samples": nn, "dF": - rt * math.log( mm ), "Error": err, "Sampling Ratio": sr, "Autocorrelation": r1 }
        return( out )

    if( clusters == 1 ):
        out = [ __integrate( dene, temperature ) ]
    else:
        out = []
        rt  = temperature * 1.0e-3 * qm3.constants.R
        exp = [ math.exp( - dene[i] / rt ) for i in range( len( dene ) ) ]
        dat = {}
        for i in range( tries ):
            grp = qm3.maths.stats.k_means( exp, clusters )[0]
            for j in range( clusters ):
                n = len( grp[j] )
                e = [ dene[exp.index(k)] for k in grp[j] ]
                k = "%d:%.3lf"%( n, sum( e ) / float( n ) )
                dat[k] = e[:]
            time.sleep( qm3.maths.rand.randint( 1, 3 ) )
        for k in iter( dat ):
            if( len( dat[k] ) > 1 ):
                out.append( __integrate( dat[k], temperature ) )
    return( out )



#
# J. Chem. Phys. v123, p144104 (2005) [10.1063/1.2052648]
# J. Chem. Phys. v124, p234106 (2006) [10.1063/1.2206775]
#
class umbint( object ):

    def __init__( self, data ):
        self.__gs = 1.0 / math.sqrt( 2.0 * math.pi )
        self.__fc = []
        self.__rf = []
        self.__mm = []
        self.__ss = []
        self.__ff = []
        self.__nw = 0
        self.__nb = None
        self.__db = 0.0
        self.__mx = None
        self.__Mx = None
        self.crd  = []
        self.pmf  = []
        # ---------------------------------------------
        for fn in data:
            f = open( fn, "rt" )
            t = [ float( i ) for i in f.readline().strip().split() ]
            self.__fc.append( t[0] )
            self.__rf.append( t[1] )
            if( not self.__mx ): self.__mx = t[1]
            if( not self.__Mx ): self.__Mx = t[1]
            n = 0.0; m = 0.0; s = 0.0
            for l in f:
                t = float( l.strip() )
                n += 1.0
                m += t
                s += t * t
                self.__mx = min( self.__mx, t )
                self.__Mx = max( self.__Mx, t )
            f.close()
            self.__ff.append( n )
            m /= n
            self.__mm.append( m )
            self.__ss.append( math.sqrt( math.fabs( s / n - m * m ) ) )
            self.__nw += 1


    def setup( self, nbins = None ):
        self.__nb = nbins
        if( not self.__nb ):
            self.__nb = 2 * self.__nw
        self.__db = ( self.__Mx - self.__mx ) / float( self.__nb )
        self.crd  = [ self.__mx + self.__db * ( float( i ) + 0.5 ) for i in range( self.__nb ) ]


    def __prob( self, x, m, s ):
        return( self.__gs / s * math.exp( - 0.5 * math.pow( ( x - m ) / s, 2.0 ) ) )


    def __dAdx( self, x, RT ):
        pt = 0.0
        px = 0.0
        for j in range( self.__nw ):
            p   = self.__ff[j] * self.__prob( x, self.__mm[j], self.__ss[j] )
            pt += p
            px += p * ( RT * ( x - self.__mm[j] ) / ( self.__ss[j] * self.__ss[j] ) - self.__fc[j] * ( x - self.__rf[j] ) )
        return( px / pt )


    def integrate( self, temperature = 300.0 ):
        self.__rt = temperature * 1.0e-3 * qm3.constants.R
        self.pmf  = [ 0.0 ]
        e = 0.0
        l = self.__dAdx( self.crd[0], self.__rt )
        for i in range( 1, self.__nb ):
            e = self.__dAdx( self.crd[i], self.__rt )
            self.pmf.append( self.pmf[-1] + 0.5 * self.__db * ( l + e ) )
            l = e
        print( "#%19s%20s"%( "Reference", "PMF" ) )
        x = max( self.pmf )
        for i in range( self.__nb ):
            print( "%20.10lf%20.10lf"%( self.crd[i], self.pmf[i] - x ) )


    # Eq 6: using all (correlated) data
    def __vdAdx( self, x, RT ):
        pt = 0.0
        px = 0.0
        for j in range( self.__nw ):
            p   = self.__ff[j] * self.__prob( x, self.__mm[j], self.__ss[j] )
            s2  = self.__ss[j] * self.__ss[j]
            pt += p
            px += p * p * RT * RT * ( 2.0 * math.pow( x - self.__mm[j], 2.0 ) + s2 ) / ( self.__ff[j] * s2 * s2 ) 
        return( px / ( pt * pt ) )


    def error( self, bin_a, bin_b = None ):
        if( bin_b == None ):
            bin_b = self.pmf.index( 0.0 )
        t_a = min( bin_a, bin_b )
        t_b = max( bin_a, bin_b )
        bin_a = t_a
        bin_b = t_b
        pmf = self.pmf[bin_b] - self.pmf[bin_a]
        # references are not sorted... (average all dispersions)
        ssm = sum( self.__ss ) / self.__nw
        err = 0.0
        for j in range( bin_a, bin_b + 1 ):
            err += self.__vdAdx( self.crd[j], self.__rt )
        err *= ( ( self.crd[bin_b] - self.crd[bin_a] ) * ssm / self.__gs - 2.0 * ssm * ssm )
        return( pmf, err )



#
# wham.F90 fDynamo module
# Comput. Phys. Communications v135, p40 (2001) [10.1016/S0010-4655(00)00215-0]
# J. Chem. Theory Comput. v6, p3713 (2010) [10.1021/ct100494z]
#
class wham( object ):

    def __init__( self, data ):
        self.__dd = []
        self.__nw = 0
        self.__fc = []
        self.__rf = []
        self.__ss = []
        self.__gm = []
        self.__gr = []
        self.__nb = None
        self.__nn = []
        self.__mx = None
        self.__Mx = None
        self.crd  = []
        self.pmf  = []
        # ---------------------------------------------
        for fn in data:
            f = open( fn, "rt" )
            t = [ float( i ) for i in f.readline().strip().split() ]
            self.__fc.append( t[0] )
            self.__rf.append( t[1] )
            if( not self.__mx ): self.__mx = t[1]
            if( not self.__Mx ): self.__Mx = t[1]
            n = 0.0
            gm = 0.0
            gr = 0.0
            for l in f:
                t = float( l.strip() )
                self.__dd.append( t )
                self.__mx = min( self.__mx, t )
                self.__Mx = max( self.__Mx, t )
                n += 1.0
                gm += t
                gr += t * t
            f.close()
            self.__ss.append( n )
            gm /= n
            self.__gm.append( gm )
            self.__gr.append( math.sqrt( math.fabs( gr / n - gm * gm ) ) )
            self.__nw += 1


    def setup( self, nbins = None ):
        self.__nb = nbins
        if( not self.__nb ):
            self.__nb = 2 * self.__nw
        self.__db = ( self.__Mx - self.__mx ) / float( self.__nb )
        self.__nn = [ 0.0 for i in range( self.__nb ) ]
        self.crd  = [ self.__mx + self.__db * ( float( i ) + 0.5 ) for i in range( self.__nb ) ]
        for t in self.__dd:
#            w = sorted( [ ( math.fabs( self.crd[i] - t ), i ) for i in range( self.__nb ) ] )[0][1]
            try:
                w = int( ( t - self.__mx ) / self.__db )
                self.__nn[w] += 1.0
            except:
                pass
        self.__ff = [ 0.0 for i in range( self.__nw ) ]
        del( self.__dd )


    def integrate( self, temperature = 300.0, maxit = 10000, toler = 1.0e-3, qprint = True ):
        __rt     = temperature * 1.0e-3 * qm3.constants.R
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
            __fo = self.__ff[:]
            for i in range( self.__nb ):
                __rr[i] = self.__nn[i] / sum( [ self.__ss[j] * math.exp( - ( __uu[i][j] - self.__ff[j] ) / __rt ) for j in range( self.__nw ) ] )
            for j in range( self.__nw ):
                self.__ff[j] = - __rt * math.log( sum( [ __eu[i][j] * __rr[i] for i in range( self.__nb ) ] ) )
            f = max( [ math.fabs( __fo[i] - self.__ff[i] ) for i in range( self.__nw ) ] ) < toler
            I += 1
        print( "#%9s%20d"%( f, I ) )
        if( f ):
            x = sum( __rr )
            for i in range( self.__nb ):
                __rr[i] /= x
                if( __rr[i] > 0.0 ):
                    self.pmf.append( - __rt * math.log( __rr[i] ) )
                else:
                    self.pmf.append( 0.0 )
            if( qprint ):
                x = max( self.pmf )
                print( "#%9s%20s%20s%20s"%( "Points", "Reference", "Density", "PMF" ) )
                for i in range( self.__nb ):
                    print( "%10.0lf%20.10lf%20.10lg%20.10lf"%( self.__nn[i], self.crd[i], __rr[i], self.pmf[i] - x ) )
        return( f )


    # mostly instructive (*very* slow): use tools/wham_bootstrapping.c instead
    def bootstrap( self, samples = 200, by_window = True, temperature = 300., maxit = 10000, toler = 1.0e-3 ):
        self.integrate( temperature, maxit, toler, False )
        m_pmf = self.pmf[:]
        self.rms = [ i * i for i in m_pmf ]
        if( by_window ):
            for ii in range( 1, samples ):
                self.__nn = [ .0 for i in range( self.__nb ) ]
                for i in range( self.__nw ):
                    j = 0
                    while( j < int( self.__ss[i] ) ):
                        t = qm3.maths.rand.gauss( self.__gm[i], self.__gr[i] )
                        try:
                            w = int( ( t - self.__mx ) / self.__db )
                            self.__nn[w] += 1.0
                            j += 1
                        except:
                            pass
                self.integrate( temperature, maxit, toler, False )
                for i in range( self.__nb ):
                    m_pmf[i] += self.pmf[i]
                    self.rms[i] += self.pmf[i] * self.pmf[i]
        else:
            pop = sum( self.__nn )
            den = [ i / pop for i in self.__nn ]
            xdn = max( den )
            nbn = self.__nb - 1
            for ii in range( 1, samples ):
                self.__nn = [ .0 for i in range( self.__nb ) ]
                c = 0
                while( c < pop ):
                    w = qm3.maths.rand.randint( 0, nbn )
                    if( den[w] >= qm3.maths.rand.random() * xdn ):
                        self.__nn[w] += 1.0
                        c += 1
                self.integrate( temperature, maxit, toler, False )
                for i in range( self.__nb ):
                    m_pmf[i] += self.pmf[i]
                    self.rms[i] += self.pmf[i] * self.pmf[i]
        print( "#%19s%20s%20s"%( "Reference", "<PMF>", "RMS" ) )
        for i in range( self.__nb ):
            self.pmf[i] = m_pmf[i] / float( samples )
            self.rms[i] = math.sqrt( math.fabs( self.rms[i] / float( samples ) - self.pmf[i] * self.pmf[i] ) )
            print( "%20.10lf%20.10lf%20.10lf"%( self.crd[i], self.pmf[i], self.rms[i] ) )



#
# J. Chem. Phys. v131, p34109 (2009) [10.1063/1.3175798]
#
class umbint_2d( object ):

    def __init__( self, data ):
        self.__gs = 1.0 / ( 2.0 * math.pi )
        self.__fc = []
        self.__rf = []
        self.__mm = []
        self.__ss = []
        self.__cx = []
        self.__ff = []
        self.__nw = 0
        self.__nb = None
        self.__db = 0.0
        self.__mx = None
        self.__Mx = None
        self.crd  = []
        self.pmf  = []
        # ---------------------------------------------
        for fn in data:
            f = open( fn, "rt" )
            t = [ float( i ) for i in f.readline().strip().split() ]
            self.__fc.append( t[0] ); self.__fc.append( t[2] )
            self.__rf.append( t[1] ); self.__rf.append( t[3] )
            if( not self.__mx ): self.__mx = [ t[1], t[3] ]
            if( not self.__Mx ): self.__Mx = [ t[1], t[3] ]
            n = 0.0; m1 = 0.0; s1 = 0.0; m2 = 0.0; s2 = 0.0; cx = 0.0
            for l in f:
                t   = [ float( i ) for i in l.strip().split() ]
                n  += 1.0
                m1 += t[0]
                s1 += t[0] * t[0]
                m2 += t[1]
                s2 += t[1] * t[1]
                cx += t[0] * t[1]
                self.__mx[0] = min( self.__mx[0], t[0] )
                self.__Mx[0] = max( self.__Mx[0], t[0] )
                self.__mx[1] = min( self.__mx[1], t[1] )
                self.__Mx[1] = max( self.__Mx[1], t[1] )
            f.close()
            self.__ff.append( n )
            m1 /= n
            self.__mm.append( m1 )
            self.__ss.append( math.sqrt( math.fabs( s1 / n - m1 * m1 ) ) )
            m2 /= n
            self.__mm.append( m2 )
            self.__ss.append( math.sqrt( math.fabs( s2 / n - m2 * m2 ) ) )
            self.__cx.append( cx / n - m1 * m2 )
            self.__nw += 1


    def setup( self, nbins ):
        self.__nb = nbins[:]
        self.__db = [ ( self.__Mx[i] - self.__mx[i] ) / float( self.__nb[i] ) for i in [0, 1] ]
        self.crx  = [ self.__mx[0] + self.__db[0] * ( float( i ) + 0.5 ) for i in range( self.__nb[0] ) ]
        self.cry  = [ self.__mx[1] + self.__db[1] * ( float( i ) + 0.5 ) for i in range( self.__nb[1] ) ]


    #P( x,y ) = 1 over { 2 %pi s_x s_y sqrt{1 - %rho^2} } exp left[ - {z  over {2 (1- %rho^2)} } right] newline
    #%rho = V_{x,y} over { s_x s_y } ~~~~ 
    #z = left( {x - m_x} over s_x right)^2 + left( {y - m_y} over s_x right)^2 -{ 2 %rho (x-m_x)(y-m_y)  } over{ s_x s_y }
    def __prob( self, x, mx, sx, y, my, sy, c ):
        r = c / ( sx * sy )
        t = 1.0 - r * r
        z = math.pow( ( x - mx ) / sx, 2.0 ) + math.pow( ( y - my ) / sy, 2.0 ) - 2.0 * r * ( x - mx ) * ( y - my ) / ( sx * sy )
        return( self.__gs / ( sx * sy * math.sqrt( t ) ) * math.exp( - 0.5 * z / t ) )


    #dAdz_{ x,y/nx·ny } ( %xi_{x,y} ) = left[ 1 over %beta left[ matrix{ s_x^2 # V_{x,y} ## V_{x,y} # s_y^2 } right]^{-1}  ( %xi - langle %xi rangle )_{x,y} - left[ matrix{ K_x# 0## 0# K_y} right] ( %xi -  %xi^r )_{x,y} right] { N_i P_i(%xi) } over { sum_j N_j P_j(%xi) }
    def integrate( self, temperature = 300.0 ):
        __rt = temperature * 1.0e-3 * qm3.constants.R
        dAdz = []
        for i in range( self.__nb[0] ):
            for j in range( self.__nb[1] ):
                pt = 0.0
                pk = [ 0.0, 0.0 ]
                for k in range( self.__nw ):
                    p      = self.__ff[k] * self.__prob( self.crx[i], self.__mm[2*k], self.__ss[2*k], self.cry[j], self.__mm[2*k+1], self.__ss[2*k+1], self.__cx[k] )
                    pt    += p
                    dt     = math.pow( self.__ss[2*k] * self.__ss[2*k+1], 2.0 ) - self.__cx[k] * self.__cx[k]
                    pk[0] += p * ( __rt / dt * ( self.__ss[2*k+1] * self.__ss[2*k+1] * ( self.crx[i] - self.__mm[2*k] ) - self.__cx[k] * ( self.cry[j] - self.__mm[2*k+1] ) ) - self.__fc[2*k] * ( self.crx[i] - self.__rf[2*k] ) )
                    pk[1] += p * ( __rt / dt * ( self.__ss[2*k] * self.__ss[2*k] * ( self.cry[j] - self.__mm[2*k+1] ) - self.__cx[k] * ( self.crx[i] - self.__mm[2*k] ) ) - self.__fc[2*k+1] * ( self.cry[j] - self.__rf[2*k+1] ) )
                dAdz.append( [ pk[0] / pt, pk[1] / pt ] )
        self.pmf = qm3.maths.ode.least_squares_finite_elements_2d( self.__nb[0], self.__db[0], self.__nb[1], self.__db[1], dAdz )
        x = min( self.pmf )
        k = 0
        fd = open( "umbint_2d", "wt" )
        for i in range( self.__nb[0] ):
            for j in range( self.__nb[1] ):
                self.pmf[k] -= x
                fd.write( "%20.10lf%20.10lf%20.10lf\n"%( self.crx[i], self.cry[j], self.pmf[k] ) )
                k += 1
            fd.write( "\n" )
        fd.close()



try:
    import matplotlib.pyplot
    def plot_data( data, dsigma = 2.0 ):
        __mm = []
        __ss = []
        __dd = []
        __Mx = 0
        for fn in data:
            f = open( fn, "rt" )
            t = [ float( i ) for i in f.readline().strip().split() ]
            n = 0.0
            m = 0.0
            s = 0.0
            __dd.append( [] )
            for l in f:
                t = float( l.strip() )
                m += t
                s += t * t
                n += 1.0
                __dd[-1].append( t )
            f.close()
            m /= n
            __mm.append( m )
            __ss.append( math.sqrt( math.fabs( s / n - m * m ) ) )
            __Mx = max( __Mx, n )
        t = sorted( __mm )
        __mw = __mm.index( t[0]  )
        __Mw = __mm.index( t[-1] )
        matplotlib.pyplot.grid( True )
        matplotlib.pyplot.xlim( 0.0, __Mx )
        matplotlib.pyplot.ylim( __mm[__mw] - 2.0 * __ss[__mw], __mm[__Mw] + 2.0 * __ss[__Mw] )
        for i in range( len( __mm ) ):
            x = [ float( j ) for j in range( len( __dd[i] ) ) ]
            matplotlib.pyplot.plot( x, __dd[i] )
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig( "overlap_coor.pdf" )
        matplotlib.pyplot.clf()
        matplotlib.pyplot.xlim( __mm[__mw] - 2.0 * __ss[__mw], __mm[__Mw] + 2.0 * __ss[__Mw] )
        f = 1.0 / math.sqrt( 2.0 * math.pi )
        for i in range( len( __mm ) ):
            nx = 100
            mx = __mm[i] - dsigma * __ss[i]
            dx = 2.0 * ( dsigma * __ss[i] ) / float( nx )
            x = [ mx + float( j ) * dx for j in range( nx + 1 ) ]
            y = [ f / __ss[i] * math.exp( - 0.5 * math.pow( ( x[j] - __mm[i] ) / __ss[i], 2.0 ) ) for j in range( nx + 1 ) ]
            matplotlib.pyplot.plot( x, y )
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig( "overlap_gaus.pdf" )
except:
        pass




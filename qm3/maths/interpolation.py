# -*- coding: iso-8859-1 -*-
import math
        


def find_center( rx, x ):
    lo = 0
    hi = len( x ) - 2
    if( rx <= x[0] ):
        return( lo )
    elif( rx >= x[hi] ):
        return( hi )
    else:
        mm = ( lo + hi ) // 2
        while( rx < x[mm] or rx >= x[mm+1] ):
            if( rx < x[mm] ):
                hi = mm
            else:
                lo = mm
            mm = ( lo + hi ) // 2
        return( mm )



class gaussian( object ):

    def __init__( self, x, y, sigma = 0.1 ):
        self.n = len( x )
        self.x = []
        self.y = []
        for i,j in sorted( [ ( x[k], y[k] ) for k in range( self.n ) ] ):
            self.x.append( i )
            self.y.append( j )
        self.sigma = sigma


    def calc( self, rx ):
        ry = 0.0
        dy = 0.0
        s  = 0.0
        for i in range( self.n ):
            w   = math.exp( - math.pow( ( rx - self.x[i] ) / self.sigma, 2.0 ) )
            ry += self.y[i] * w
            dy -= self.y[i] * w * 2.0 * ( rx - self.x[i] ) / self.sigma
            s  += w
        return( ry / s, dy / s )



class bezier( object ):

    def __init__( self, x, y ):
        self.n = len( x )
        self.x = []
        self.y = []
        for i,j in sorted( [ ( x[k], y[k] ) for k in range( self.n ) ] ):
            self.x.append( i )
            self.y.append( j )


    # iterative deCasteljau
    def calc( self, rx ):
        u = max( 0.0, min( 1.0, ( rx - self.x[0] ) / ( self.x[-1] - self.x[0] ) ) )
        t = self.y[:]
        for k in range( 1, self.n ):
            for i in range( self.n - k ):
                t[i] = ( 1.0 - u ) * t[i] + u * t[i+1]
        return( t[0], None )



class bspline( object ):

    def __init__( self, x, y, order = 2 ):
        self.n = len( x )
        self.x = []
        self.y = []
        for i,j in sorted( [ ( x[k], y[k] ) for k in range( self.n ) ] ):
            self.x.append( i )
            self.y.append( j )
        self.order = order


    def calc( self, rx ):
        i = min( max( 0, find_center( rx, self.x ) ), self.n - 1 )
        p = min( min( i, self.order ), self.n - 1 - i )
        N = [ 1.0 ] + [ 0.0 for j in range( p ) ]
        L = [ 0.0 for j in range( p + 1 ) ]
        R = [ 0.0 for j in range( p + 1 ) ]
        N[0] = 1.0
        for j in range( 1, p + 1 ):
            L[j] = rx - self.x[i+1-j]
            R[j] = self.x[i+j] - rx
            s    = 0.0
            for r in range( j ):
                t    = N[r] / ( R[r+1] + L[j-r] )
                N[r] = s + R[r+1] * t
                s    = L[j-r] * t
            N[j] = s
        ry = 0.0
        for j in range( p + 1 ):
            ry += N[j] * self.y[i-p+j+1]
        return( ry, None )



class linear( object ):

    def __init__( self, x, y ):
        self.n = len( x )
        self.x = []
        self.y = []
        for i,j in sorted( [ ( x[k], y[k] ) for k in range( self.n ) ] ):
            self.x.append( i )
            self.y.append( j )


    def calc( self, rx ):
        i  = find_center( rx, self.x )
        ry = self.y[i] + ( self.y[i+1] - self.y[i] ) * ( rx - self.x[i] ) / ( self.x[i+1] - self.x[i] )
        return( ry, None )



class lagrange( object ):

    def __init__( self, x, y, points = 2 ):
        self.n = len( x )
        self.x = []
        self.y = []
        for i,j in sorted( [ ( x[k], y[k] ) for k in range( self.n ) ] ):
            self.x.append( i )
            self.y.append( j )
        self.points = points


    def calc( self, rx ):
        i  = find_center( rx, self.x ) + 1
        i0 = max( 1, i - self.points )
        ix = min( self.n - 1, i + self.points )
        if( i0 == 0 and ix+self.points < self.n ):
            ix += self.points
        if( ix == self.n-1 and i0-self.points > 0 ):
            i0 -= self.points
        ry = 0.0
        dy = 0.0
        for i in range( i0, ix + 1 ):
            t = 1.0
            v = 0.0
            for j in range( i0, ix + 1 ):
                if( i != j ):
                    t *= ( rx - self.x[j] ) / ( self.x[i] - self.x[j] )
                    # - derivative ------------------------------------
                    u = 1.0 / ( self.x[i] - self.x[j] )
                    for k in range( i0, ix + 1 ):
                        if( k != i and k != j ):
                            u *= ( rx - self.x[k] ) / ( self.x[i] - self.x[k] )
                    v += u
                    # -------------------------------------------------
            ry += self.y[i] * t
            dy += self.y[i] * v
        return( ry, dy )



class aitken( object ):

    def __init__( self, x, y ):
        self.n = len( x )
        self.x = []
        self.y = []
        for i,j in sorted( [ ( x[k], y[k] ) for k in range( self.n ) ] ):
            self.x.append( i )
            self.y.append( j )


    def calc( self, rx ):
        i  = find_center( rx, self.x )
        # -----------------------------------------------
        if( i < 1 ):
            x = self.x[0:4]
            y = self.y[0:4]
        elif( i > self.n - 3 ):
            x = self.x[-4:]
            y = self.y[-4:]
        else:
            x = self.x[i-1:i+3]
            y = self.y[i-1:i+3]
        # -----------------------------------------------
        x0s  = x[0] * x[0]
        x1s  = x[1] * x[1]
        x2s  = x[2] * x[2]
        x3s  = x[3] * x[3]
        dx01 = x[0] - x[1]
        dx02 = x[0] - x[2]
        dx03 = x[0] - x[3]
        dx12 = x[1] - x[2]
        dx13 = x[1] - x[3]
        dx23 = x[2] - x[3]
        den  = dx01 * dx02 * dx03 * dx12 * dx13 * dx23
        dy01 = y[0] - y[1]
        dy02 = y[0] - y[2]
        dy03 = y[0] - y[3]
        dy12 = y[1] - y[2]
        dy13 = y[1] - y[3]
        dy23 = y[2] - y[3]
        # -----------------------------------------------
        a0  = x[0] * dx02 * x[2] * dx03 * dx23 * x[3] * y[1] + x1s * x[1] * ( x[0] * dx03 * x[3] * y[2] + x2s * ( x[0] * y[3] - x[3] * y[0] ) + x[2] * ( x3s * y[0] - x0s * y[3] ) ) + x[1] * ( x0s * dx03 * x3s *y[2] + x2s * x[2] * ( x0s * y[3] - x3s * y[0] ) + x2s * ( x3s * x[3] * y[0] - x0s * x[0] * y[3] ) ) + x1s * ( x[0] * x[3] *( x3s - x0s ) * y[2] + x2s * x[2] * ( x[3] * y[0] - x[0] * y[3] ) + x[2] * ( x0s * x[0] * y[3] - x3s * x[3] * y[0] ) )
        a1 = x0s * dx03 * x3s * dy12 + x2s * x[2] * ( x3s * dy01 + x0s * dy13 ) + x1s * ( x3s * x[3] * dy02 + x0s * x[0] * dy23 - x2s * x[2] * dy03 ) - x2s * ( x3s * x[3] * dy01 + x0s * x[0] * dy13 ) + x1s * x[1] * ( - x3s * dy02 + x2s * dy03 - x0s * dy23 )
        a2 = - x[0] * x[3] * ( x0s - x3s ) * dy12 + x[2] * ( x3s * x[3] * dy01 + x0s * x[0] * dy13 ) + x1s * x[1] * ( x[3] * dy02 + x[0] * dy23 - x[2] * dy03 ) - x2s * x[2] * ( x[3]* dy01 + x[0] * dy13 ) + x[1] * ( - x3s * x[3] * dy02 + x2s * x[2] * dy03 - x0s * x[0] * dy23 )
        
        a3 = x[0] * dx03 * x[3] * dy12 + x2s * ( x[3] * dy01 + x[0] * dy13 ) + x[1] * ( x3s * dy02 + x0s * dy23 - x2s * dy03 ) - x[2] * ( x3s * dy01 + x0s * dy13 ) + x1s * ( - x[3] * dy02 + x[2] * dy03 - x[0] * dy23 )
        # -----------------------------------------------
        ry = ( a0 + rx * ( a1 + rx * ( a2 + rx * a3 ) ) ) / den
        dy = ( a1 + rx * ( 2.0 * a2 + 3.0 * rx * a3 ) ) / den
        return( ry, dy )



class cubic_spline( object ):

    def __init__( self, x, y ):
        self.n = len( x )
        self.x = []
        self.y = []
        for i,j in sorted( [ ( x[k], y[k] ) for k in range( self.n ) ] ):
            self.x.append( i )
            self.y.append( j )
        self.y2 = [ 0.0 for i in range( self.n ) ]
        u = [ 0.0 for i in range( self.n ) ]
        u[0] = 0.0
        u[self.n-1] = 0.0
        self.y2[0] = 0.0
        self.y2[self.n-1] = 0.0
        for i in range( 1, self.n - 1 ):
            s = ( self.x[i] - self.x[i-1] ) / ( self.x[i+1] - self.x[i-1] )
            p = s * self.y2[i-1] + 2.0
            self.y2[i] = ( s - 1.0 ) / p
            u[i]=( 6.0 * ( ( self.y[i+1] - self.y[i] ) / ( self.x[i+1] - self.x[i] ) - ( self.y[i] - self.y[i-1] ) / ( self.x[i] - self.x[i-1] ) ) / ( self.x[i+1] - self.x[i-1] ) - s * u[i-1] ) / p
        for i in range( self.n-2, -1, -1 ):
            self.y2[i] = self.y2[i] * self.y2[i+1] + u[i]

    
    def calc( self, rx ):
        klo = max( 0, find_center( rx, self.x ) )
        khi = min( self.n - 1, klo + 1 )
        h   = self.x[khi] - self.x[klo]
        a   = ( self.x[khi] - rx ) / h
        b   = ( rx - self.x[klo] ) / h
        ry  = a * self.y[klo] + b * self.y[khi] + ( ( a * a * a - a ) * self.y2[klo] + ( b * b * b - b ) * self.y2[khi] ) * ( h * h ) / 6.0
        dy  = ( self.y[khi] - self.y[klo] ) / h + h * ( ( 3.0 * b * b - 1.0 ) * self.y2[khi] - ( 3.0 * a * a - 1.0 ) * self.y2[klo] ) / 6.0
        return( ry, dy )



class hermite_spline( object ):
    """
    Available methods: steffen  /  akima  /  [fritsch_carlson]
    """
    def __init__( self, x, y, method = "fritsch_carlson" ):
        self.n = len( x )
        self.x = []
        self.y = []
        for i,j in sorted( [ ( x[k], y[k] ) for k in range( self.n ) ] ):
            self.x.append( i )
            self.y.append( j )
        dx = []
        dy = []
        m  = []
        self.c1 = []
        self.c2 = []
        self.c3 = []
        for i in range( self.n - 1 ):
            dx.append( self.x[i+1] - self.x[i] )
            dy.append( self.y[i+1] - self.y[i] )
            m.append( dy[i] / dx[i] )
        # -------------------------------------------------------------------
        # Steffen
        if( method == "steffen" ):
            self.c1.append( m[0] )
            for i in range( self.n - 2 ):
                self.c1.append( ( math.copysign( 1.0, m[i] ) + math.copysign( 1.0, m[i+1] ) ) * min( math.fabs( m[i] ), math.fabs( m[i+1] ), 0.5 * math.fabs( ( dx[i] * m[i+1] + dx[i+1] * m[i] ) / ( dx[i] + dx[i+1] ) ) ) )
            self.c1.append( m[-1] )
        # -------------------------------------------------------------------
        # Akima
        elif( method == "akima" ):
            M  = [ 2.0 * m[0] - m[1], 2.0 * m[0] - m[1] ] + m[:]
            M += [ 2.0 * m[-1] - m[-2], 2.0 * ( 2.0 * m[-1] - m[-2] ) - m[-1] ]
            for i in range( self.n ):
                a = math.fabs( M[i+3] - M[i+2] )
                b = math.fabs( M[i+1] - M[i] )
                if( a+b > 0.0 ):
                    self.c1.append( ( b * M[i+2] + a * M[i+1] ) / ( a + b ) )
                else:
                    self.c1.append( ( M[i+2] + M[i+1] ) / 2.0 )
        # -------------------------------------------------------------------
        # Fritsch-Carlson
        else:
            self.c1.append( m[0] )
            for i in range( self.n - 2 ):
                if( m[i] * m[i+1] <= 0.0 ):
                    self.c1.append( 0.0 )
                else:
                    t = dx[i] + dx[i+1]
                    self.c1.append( 3.0 * t / ( ( t + dx[i+1] ) / m[i] + ( t + dx[i] ) / m[i+1] ) )
            self.c1.append( m[-1] )
        # -------------------------------------------------------------------
        for i in range( self.n - 1 ):
            t = self.c1[i] + self.c1[i+1] - m[i] - m[i]
            self.c2.append( ( m[i] - self.c1[i] - t ) / dx[i] )
            self.c3.append( t / ( dx[i] * dx[i] ) )


    def calc( self, rx ):
        i  = max( 0, find_center( rx, self.x ) )
        h  = rx - self.x[i]
        h2 = h * h
        ry = self.y[i] + ( self.c1[i] + self.c2[i] * h + self.c3[i] * h2 ) * h
        dy = self.c1[i] + ( 2.0 * self.c2[i] + 3.0 * self.c3[i] * h ) * h
        return( ry, dy )



class interpolate_2d( object ):
    """
    Z Data should be porperly sorted as fixed_X, changing_Y (see grids.parse for sorting)
    
    Interpolant can make use of a lambda function:

            interpolant = lambda x,y: interpolation.hermite_spline( x, y, method = "akima" )
    """
    def __init__( self, x, y, z, interpolant = cubic_spline ):
        self.nx = len( x )
        self.ny = len( y )
        self.x  = x[:]
        self.y  = y[:]
        self.z  = z[:]
        self.II = interpolant
        self.Ix = []
        for i in range( self.nx ):
            t = [ self.z[self.ny*i+j] for j in range( self.ny ) ]
            self.Ix.append( self.II( self.y, t ) )
        self.Iy = []
        for j in range( self.ny ):
            t = [ self.z[self.ny*i+j] for i in range( self.nx ) ]
            self.Iy.append( self.II( self.x, t ) )


    def calc( self, rx, ry ):
        t  = [ self.Ix[i].calc( ry )[0] for i in range( self.nx ) ]
        ox = self.II( self.x, t ).calc( rx )
        t  = [ self.Iy[i].calc( rx )[0] for i in range( self.ny ) ]
        oy = self.II( self.y, t ).calc( ry )
        return( ( ox[0] + oy[0] ) * 0.5, ox[1], oy[1] )


try:
    import qm3.maths._fitpack
    class tensioned_spline:

        def __init__( self, x, y, tension = 0.0 ):
            self.n = len( x )
            self.s = tension
            self.x = []
            self.y = []
            for i,j in sorted( [ ( x[k], y[k] ) for k in range( self.n ) ] ):
                self.x.append( i )
                self.y.append( j )
            self.y2 = qm3.maths._fitpack.init( self.x, self.y, self.s )
    
        
        def calc( self, rx ):
            return( qm3.maths._fitpack.calc( rx, self.x, self.y, self.y2, self.s ) )

except:
    pass





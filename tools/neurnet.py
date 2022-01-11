#!/usr/bin/env python3
import  math
import  numpy



class model( object ):
    def __init__( self, isize, layers ):
        self.isiz = isize
        self.coef = [ numpy.random.randn( isize, layers[0] ), numpy.zeros( ( 1, layers[0] ) ) ]
        for i in range( 1, len( layers ) ):
            self.coef += [ numpy.random.randn( layers[i-1], layers[i] ), numpy.zeros( ( 1, layers[i] ) ) ]
        self.coef += [ numpy.random.randn( layers[-1], 1 ), numpy.zeros( ( 1, 1 ) ) ]
        self.oref = 0.0


    def __random( self ):
        for i in range( 0, len( self.coef ), 2 ):
            self.coef[i]   = numpy.random.randn( self.coef[i].shape[0], self.coef[i].shape[1] )
            self.coef[i+1] = numpy.zeros( self.coef[i+1].shape )


    @staticmethod
    def f_actv( vec ):
        return( numpy.log( 1.0 + numpy.exp( numpy.where( vec > 709.7827, 709.7827, vec ) ) ) )

    @staticmethod
    def g_actv( vec ):
        return( 1.0 / ( 1.0 + numpy.exp( numpy.where( vec < -709.7827, -709.7827, - vec ) ) ) )
    

    def get_func( self, inp ):
        out = numpy.array( inp ).reshape( ( 1, self.isiz ) )
        for i in range( 0, len( self.coef ), 2 ):
            out = self.f_actv( numpy.dot( out, self.coef[i] ) + self.coef[i+1] )
        return( float( out ) + self.oref )


    def get_grad( self, inp ):
        f_tmp = numpy.dot( numpy.array( inp ).reshape( ( 1, self.isiz ) ), self.coef[0] ) + self.coef[1]
        g_tmp = self.coef[0] * self.g_actv( f_tmp )
        for i in range( 2, len( self.coef ), 2 ):
            f_tmp = numpy.dot( self.f_actv( f_tmp ), self.coef[i] ) + self.coef[i+1]
            g_tmp = numpy.dot( g_tmp, self.coef[i] * self.g_actv( f_tmp ) )
        return( float( self.f_actv( f_tmp ) ), g_tmp.flatten().tolist() )


    def __grad( self, iset, oset ):
        grd = []
        for i in range( len( self.coef ) ):
            grd.append( numpy.zeros( self.coef[i].shape ) )
        los = 0.0
        for i in range( len( oset ) ):
            buf = []
            inp = numpy.array( iset[i] ).reshape( ( 1, self.isiz ) ) 
            tmp = inp.copy()
            for j in range( 0, len( self.coef ), 2 ):
                tmp = numpy.dot( tmp, self.coef[j] ) + self.coef[j+1]
                buf.append( tmp.copy() )
                tmp = self.f_actv( tmp )
                buf.append( tmp.copy() )
            dif = 2.0 * ( buf[-1] - numpy.array( oset[i] ).reshape( ( 1, 1 ) ) )
            los += 0.25 * float( dif ) * float( dif )
            f_idx = -2
            g_idx = -1
            g_tmp = [ dif * self.g_actv( buf[f_idx] ) ]
            for j in range( len( self.coef ) - 2, 0, -2 ):
                f_idx -= 1
                g_tmp.insert( 0, numpy.dot( buf[f_idx].T, g_tmp[g_idx] ) )
                f_idx -= 1
                g_tmp.insert( 0, self.g_actv( buf[f_idx] ) * numpy.dot( g_tmp[g_idx], self.coef[j].T ) )
                g_idx -= 2
            g_tmp.insert( 0, numpy.dot( inp.T, g_tmp[g_idx] ) )
            for j in range( len( self.coef ) ):
                grd[j] += g_tmp[j]
        out = []
        for m in grd:
            out += m.flatten().tolist()
        return( los, numpy.array( out ) )


    def __fire( self, iset, oset, step_size = 0.0001, step_number = 1000, max_tries = 2 ):
        coor = []
        shap = []
        for c in self.coef:
            coor += c.flatten().tolist()
            shap.append( c.shape )
        size = len( coor )
        coor = numpy.array( coor )
        nstp = 0
        ssiz = step_size
        alph = 0.1
        dstp = 5
        velo = numpy.zeros( size )
        step = numpy.zeros( size )
        loss, grad = self.__grad( iset, oset )
        last = loss
        xlss = loss
        xcof = [ m.copy() for m in self.coef ]
        norm = math.sqrt( numpy.dot( grad, grad ) )
        print( "          %30.1lf%30.1lf"%( xlss, norm ) )
        nfun = 0
        it = 0
        while( it < step_number and nfun < max_tries ):
            if( - numpy.dot( velo, grad ) > 0.0 ):
                vsiz = math.sqrt( numpy.dot( velo, velo ) )
                velo = ( 1.0 - alph ) * velo - alph * grad / norm * vsiz
                if( nstp > dstp ):
                    ssiz = min( ssiz * 1.1, step_size )
                    alph *= 0.99
                nstp += 1
            else:
                alph = 0.1
                ssiz *= 0.5
                nstp = 0
                velo = numpy.zeros( size )
            velo -= ssiz * grad
            step = ssiz * velo
            tmp  = math.sqrt( numpy.dot( step, step ) )
            if( tmp > ssiz ):
                step *= ssiz / tmp
            coor += step
            tmp = 0
            for j in range( len( shap ) ):
                dsp = shap[j][0] * shap[j][1]
                self.coef[j] = coor[tmp:tmp+dsp].reshape( shap[j] )
                tmp += dsp
            loss, grad = self.__grad( iset, oset )
            norm = math.sqrt( numpy.dot( grad, grad ) )
            nfun += loss >= last
            last = loss
#            print( "%10d%30.1lf%30.1lf%30.1lf"%( it, loss, norm, loss - xlss ) )
            if( loss < xlss ):
                xlss = loss
                xcof = [ m.copy() for m in self.coef ]
            it += 1
        print( "%10d%30.1lf%30.1lf%30.1lf"%( it, loss, norm, loss - xlss ) )
        self.coef = xcof
        return( xlss )


    def train( self, iset, oset, loss_tolerance = 1.0, population = 10 ):
        dim = len( iset )
        tol = loss_tolerance * dim
        print( "cutoff    : %.1lf x %d = %.1lf"%( loss_tolerance, dim, tol ) )
        print( "population: %d"%( population ) )
        acc = 0
        x_loss = self.__fire( iset, oset, step_size = 0.01, step_number = 1000 )
        x_coef = [ m.copy() for m in self.coef ]
        while( acc < population ):
            self.__random()
            loss = self.__fire( iset, oset, step_size = 0.01, step_number = 1000 )
            if( loss < tol ):
                acc += 1
                if( loss < x_loss ):
                    x_loss = loss
                    x_coef = [ m.copy() for m in self.coef ]
                print( "try [%4d/%d] %.6lf"%( acc, population, loss ) )
        self.coef = x_coef
        print( "restoring best fit (%.6lf)"%( x_loss ) )





if( __name__ == "__main__" ):
    import  os
    import  sys
    import  pickle
    import  matplotlib.pyplot as plt
    
    x   = numpy.linspace( -2, 2, 11 )
    inp = []
    out = []
    for i in range( len( x ) ):
        for j in range( len( x ) ):
            inp.append( [ x[i], x[j] ] )
            out.append( x[i] * x[i] + x[j] * x[j] )

    obj = model( 2, [ 12, 4 ] )
    obj.train( inp, out, loss_tolerance = 0.1 )

    x  = numpy.linspace( -2, 2, 21 )
    nx = []
    ny = []
    nz = []
    for i in range( 21 ):
        for j in range( 21 ):
            nx.append( x[i] )
            ny.append( x[j] )
            nz.append( obj.get_func( [ x[i], x[j] ] ) )
    plt.clf()
    plt.grid( True )
    fig = plt.figure()
    ax3 = fig.add_subplot( 111, projection = "3d" )
    ax3.scatter( [ i[0] for i in inp ], [ i[1] for i in inp ], out )
    ax3.scatter( nx, ny, nz )
    plt.show()

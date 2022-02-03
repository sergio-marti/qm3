# -*- coding: iso-8859-1 -*-
import  math
import  qm3.maths.matrix
import  qm3.utils

try:
    import  qm3.engines._ml_info
    ml_info_so = True
except:
    ml_info_so = False



class model( object ):
    """
    ===  Sequential Softplus-activated Dense Layers  ===

    import  tensorflow as tf
    import  numpy as np

    ref = np.min( out )
    out -= ref
    mod = tf.keras.models.Sequential( [ 
        tf.keras.layers.Dense( 256, activation = "softplus" ),
        tf.keras.layers.Dense( 128, activation = "softplus" ),
        tf.keras.layers.Dense(  64, activation = "softplus" ),
        tf.keras.layers.Dense(  32, activation = "softplus" ),
        tf.keras.layers.Dense(   1, activation = "softplus" ) ] )
    mod.compile( optimizer = tf.keras.optimizers.Adam( learning_rate = 0.0001 ),
        loss = "MeanSquaredError" )
    stp = tf.keras.callbacks.ModelCheckpoint( "model.h5", monitor = "loss",
        verbose = 1, mode = "min", save_best_only = True )
    mod.fit( inp, out, batch_size = 2048, epochs = 10000, callbacks = [ stp ], verbose = 0 )

    mod = tf.keras.models.load_model( "model.h5" )
    obj = qm3.engines.ml_model.model()
    for cof in mod.trainable_weights:
        tmp = cof.numpy()
        obj.coef.append( tmp.flatten().tolist() )
        obj.dime.append( tmp.shape )
    obj.oref = float( ref )
    """
    def __init__( self ):
        self.coef = []
        self.dime = []
        self.oref = 0.0


    @staticmethod
    def f_actv( vec ):
        # math.log( sys.float_info.max )
        return( [ math.log( 1.0 + math.exp( min( 709.7827, x ) ) ) for x in vec ] )


    @staticmethod
    def g_actv( vec ):
        # math.log( sys.float_info.max )
        return( [ 1.0 / ( 1.0 + math.exp( min( 709.7827, - x ) ) ) for x in vec ] )


    # -- scalar
    def get_func( self, data ):
        siz = len( data )
        out = data[:]
        for i in range( 0, len( self.dime ), 2 ):
            out = qm3.maths.matrix.mult( out, 1, siz, self.coef[i], self.dime[i][0], self.dime[i][1] )
            for j in range( self.dime[i][1] ):
                out[j] += self.coef[i+1][j]
            out = self.f_actv( out )
            siz = self.dime[i][1]
        return( out[0] + self.oref )


    # -- same dimension as input data
    def get_grad( self, data ):
        f_tmp = qm3.maths.matrix.mult( data, 1, len( data ),
            self.coef[0], self.dime[0][0], self.dime[0][1] )
        for j in range( self.dime[0][1] ):
            f_tmp[j] += self.coef[1][j]
        f_siz = self.dime[0][1]
        g_tmp = self.coef[0][:]
        x_tmp = self.g_actv( f_tmp )
        for i in range( self.dime[0][0] ):
            for j in range( self.dime[0][1] ):
                g_tmp[i*self.dime[0][1]+j] *= x_tmp[j]
        g_siz = [ self.dime[0][0], self.dime[0][1] ]
        for k in range( 2, len( self.dime ), 2 ):
            f_tmp = qm3.maths.matrix.mult( self.f_actv( f_tmp ), 1, f_siz,
                        self.coef[k], self.dime[k][0], self.dime[k][1] )
            for j in range( self.dime[k][1] ):
                f_tmp[j] += self.coef[k+1][j]
            f_siz = self.dime[k][1]
            y_tmp = self.coef[k][:]
            x_tmp = self.g_actv( f_tmp )
            for i in range( self.dime[k][0] ):
                for j in range( self.dime[k][1] ):
                    y_tmp[i*self.dime[k][1]+j] *= x_tmp[j]
            g_tmp = qm3.maths.matrix.mult( g_tmp, g_siz[0], g_siz[1],
                        y_tmp, self.dime[k][0], self.dime[k][1] )
            g_siz[1] = self.dime[k][1]
        return( self.f_actv( f_tmp )[0] + self.oref, g_tmp )



try:
    import  numpy
    
    class np_model( object ):
        """
        ===  Sequential Softplus-activated Dense Layers  ===
    
        import  tensorflow as tf
    
        ref = np.min( out )
        out -= ref
        mod = tf.keras.models.Sequential( [ 
            tf.keras.layers.Dense( 256, activation = "softplus" ),
            tf.keras.layers.Dense( 128, activation = "softplus" ),
            tf.keras.layers.Dense(  64, activation = "softplus" ),
            tf.keras.layers.Dense(  32, activation = "softplus" ),
            tf.keras.layers.Dense(   1, activation = "softplus" ) ] )
        mod.compile( optimizer = tf.keras.optimizers.Adam( learning_rate = 0.0001 ),
            loss = "MeanSquaredError" )
        stp = tf.keras.callbacks.ModelCheckpoint( "model.h5", monitor = "loss",
            verbose = 1, mode = "min", save_best_only = True )
        mod.fit( inp, out, batch_size = 2048, epochs = 10000, callbacks = [ stp ], verbose = 0 )

        mod = tf.keras.models.load_model( "model.h5" )
        obj = qm3.engines.ml_model.np_model()
        obj.coef = [ m.numpy() for m in mod.trainable_weights ]
        obj.oref = float( ref )
        """
        def __init__( self ):
            self.coef = []
            self.oref = 0.0
    
    
        @staticmethod
        def f_actv( vec ):
            # math.log( sys.float_info.max )
            return( numpy.log( 1.0 + numpy.exp( numpy.where( vec > 709.7827, 709.7827, vec ) ) ) )
    
    
        @staticmethod
        def g_actv( vec ):
            # math.log( sys.float_info.max )
            return( 1.0 / ( 1.0 + numpy.exp( numpy.where( vec < -709.7827, -709.7827, - vec ) ) ) )
    
    
        # -- scalar
        def get_func( self, data ):
            out = numpy.array( data ).reshape( ( 1, len( data ) ) )
            for i in range( 0, len( self.coef ), 2 ):
                out = self.f_actv( numpy.dot( out, self.coef[i] ) + self.coef[i+1] )
            return( float( out ) + self.oref )
    
    
        # -- same dimension as input data
        def get_grad( self, data ):
            inp = numpy.array( data ).reshape( ( 1, len( data ) ) )
            f_tmp = numpy.dot( inp, self.coef[0] ) + self.coef[1]
            g_tmp = self.coef[0] * self.g_actv( f_tmp )
            for i in range( 2, len( self.coef ), 2 ):
                f_tmp = numpy.dot( self.f_actv( f_tmp ), self.coef[i] ) + self.coef[i+1]
                g_tmp = numpy.dot( g_tmp, self.coef[i] * self.g_actv( f_tmp ) )
            return( float( self.f_actv( f_tmp ) ) + self.oref, g_tmp.flatten().tolist() )
except:
    pass


try:
    import  numpy
    import  tensorflow

    class tf_model( object ):
        def __init__( self, model ):
            self.mode = tensorflow.keras.models.load_model( model )
            self.oref = 0.0


        def get_func( self, data ):
            return( float( self.mode( numpy.array( data ).reshape( ( 1, len( data ) ) ), training = False ) ) + self.oref )


        def get_grad( self, data ):
            inp = tensorflow.convert_to_tensor( numpy.array( data ).reshape( ( 1, len( data ) ) ) )
            with tensorflow.GradientTape() as grd:
                grd.watch( inp )
                lss = self.mode( inp )
            ene = float( self.mode( inp, training = False ) ) + self.oref
            return( ene, grd.gradient( lss, inp ).numpy().flatten().tolist() )
except:
    pass



class delta_template( object ):
    def __init__( self, model, molec, sele ):
        self.mode = model
        self.sele = [ i * 3 for i in sele ]


    def get_info( self, molec ):
        raise( NotImplementedError )


    def get_jaco( self, molec ):
        raise( NotImplementedError )


    def get_func( self, molec ):
        molec.func += self.mode.get_func( self.get_info( molec ) )


    def get_grad( self, molec ):
         etmp, gtmp = self.mode.get_grad( self.get_info( molec ) )
         jaco = self.get_jaco( molec )
         siz = len( self.sele ) * 3
         dim = len( gtmp )
         grd = qm3.maths.matrix.mult( gtmp, 1, dim, jaco, dim, siz )
         molec.func += etmp
         k = 0
         for i in self.sele:
             for j in [0, 1, 2]:
                 molec.grad[i+j] += grd[k]
                 k += 1


    def num_grad( self, molec, disp = 1.e-3 ):
        molec.func += self.mode.get_func( self.get_info( molec ) )
        for i in self.sele:
            for j in [0, 1, 2]:
                bak = molec.coor[i+j]
                molec.coor[i+j] = bak + disp
                ffw = self.mode.get_func( self.get_info( molec ) )
                molec.coor[i+j] = bak - disp
                bbw = self.mode.get_func( self.get_info( molec ) )
                molec.grad[i+j] += ( ffw - bbw ) / ( 2.0 * disp )
                molec.coor[i+j] = bak



# [10.1103/PhysRevLett.108.058301]
class delta_coul( delta_template ):
    def __init__( self, model, molec, sele, anum = False ):
        delta_template.__init__( self, model, molec, sele )
        self.anum = []
        for i in range( len( sele ) ):
            if( anum ):
                self.anum.append( molec.anum[i] )
            else:
                self.anum.append( 1.0 )


    def get_info( self, molec ):
        if( ml_info_so ):
            crd = []
            for i in self.sele:
                crd += molec.coor[i:i+3][:]
            return( qm3.engines._ml_info.coul_info( self.anum, crd ) )
        else:
            siz = len( self.sele )
            dim = siz * ( siz + 1 ) // 2
            out = [ 0.0 for i in range( dim ) ]
            for i in range( siz ):
                i3 = self.sele[i]
                ww = i * siz - ( ( i - 1 ) * i ) // 2
                out[ww] = 0.5 * math.pow( self.anum[i], 2.4 )
                for j in range( i + 1, siz ):
                    j3 = self.sele[j]
                    out[ww+(j-i)] = self.anum[i] * self.anum[j] / math.sqrt( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip( molec.coor[i3:i3+3], molec.coor[j3:j3+3] ) ] ) )
            return( out )


    def get_jaco( self, molec ):
        if( ml_info_so ):
            crd = []
            for i in self.sele:
                crd += molec.coor[i:i+3][:]
            return( qm3.engines._ml_info.coul_jaco( self.anum, crd ) )
        else:
            siz = len( self.sele )
            dim = siz * ( siz + 1 ) // 2
            out = [ 0.0 for i in range( dim * siz * 3 ) ]
            row = 0
            for i in range( siz ):
                i3 = self.sele[i]
                for j in range( i, siz ):
                    if( j != i ):
                        j3 = self.sele[j]
                        dr = [ (ii-jj) for ii,jj in zip( molec.coor[i3:i3+3], molec.coor[j3:j3+3] ) ]
                        r2 = sum( [ ii * ii for ii in dr ] )
                        zz = self.anum[i] * self.anum[j] / ( r2 * sqrt( r2 ) )
                        for k in [0, 1, 2]:
                            out[row+i3+k] = - dr[k] * zz
                            out[row+j3+k] =   dr[k] * zz
                    row += siz * 3
            return( out )



# [10.1016/j.cpc.2018.03.016]
class delta_rada( delta_template ):
    def __init__( self, model, molec, sele ):
        delta_template.__init__( self, model, molec, sele )


    def setup( self, molec ):
        siz = len( self.sele )
        self.refe = []
        for i in range( siz ):
            dst = []
            for j in range( siz ):
                rr = [ ii - jj for ii,jj in zip( molec.coor[self.sele[i]:self.sele[i]+3], molec.coor[self.sele[j]:self.sele[j]+3] ) ]
                dst.append( ( math.sqrt( sum( [ ii * ii for ii in rr ] ) ), j * 3, rr ) )
            dst.sort()
            self.refe.append( dst[1][1] )
            vec = [ ii / dst[1][0] for ii in dst[1][2] ]
            nrm = []
            for j in range( 2, len( self.sele ) ):
                nrm.append( ( math.fabs( sum( [ ii * jj for ii,jj in zip( dst[j][2], vec ) ] ) ), dst[j][1] ) )
            nrm.sort()
            self.refe.append( nrm[1][1] )


    def get_info( self, molec, coor = None ):
        crd = []
        for i in self.sele:
            crd += molec.coor[i:i+3][:]
        if( ml_info_so ):
            return( qm3.engines._ml_info.rada_info( self.refe, crd ) )
        else:
            siz = len( self.sele )
            out = []
            for i in range( siz ):
                i3 = i * 3
                a3 = self.refe[2*i]
                b3 = self.refe[2*i+1]
                e1 = qm3.maths.matrix.norm( [ jj - ii for ii,jj in zip( crd[i3:i3+3], crd[a3:a3+3] ) ] )
                e2 = [ jj - ii for ii,jj in zip( crd[i3:i3+3], crd[b3:b3+3] ) ]
                dd = sum( [ ii * jj for ii,jj in zip( e1, e2 ) ] )
                e2 = qm3.maths.matrix.norm( [ ii - dd * jj for ii,jj in zip( e2, e1 ) ] )
                e3 = [ e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0] ]
                for j in range( siz ):
                    if( j != i ):
                        j3 = j * 3
                        rr = [ jj - ii for ii,jj in zip( crd[i3:i3+3], crd[j3:j3+3] ) ]
                        mm = math.sqrt( sum( [ ii * ii for ii in rr ] ) )
                        out.append( 1. / mm )
                        out.append( ( rr[0] * e1[0] + rr[1] * e1[1] + rr[2] * e1[2] ) / mm )
                        out.append( ( rr[0] * e2[0] + rr[1] * e2[1] + rr[2] * e2[2] ) / mm )
                        out.append( ( rr[0] * e3[0] + rr[1] * e3[1] + rr[2] * e3[2] ) / mm )
            return( out )



# [10.1063/1.3553717]
class delta_acsf( delta_template ):
    def __init__( self, model, molec, sele ):
        delta_template.__init__( self, model, molec, sele )


    def setup( self, cutx = 8.0, eta2 = [ 3.57, 1.43, 0.71, 0.36, 0.21, 0.11 ], dse5 = 1.0, eta5 = [ 1.07, 0.71, 0.36, 0.18 ] ):
        self.cutx = cutx
        self.eta2 = eta2
        self.dse5 = dse5
        self.pre5 = math.pow( 2.0, 1.0 - dse5 )
        self.eta5 = eta5
        self.size = len( eta2 ) + len( eta5 )


    @staticmethod
    def __fcut( dst, cut ):
        if( dst > cut ):
            return( 0.0 )
        else:
            return( 0.5 * ( math.cos( math.pi * dst / cut ) + 1.0 ) )
    
    @staticmethod
    def __dist( ci, cj ):
        return( math.sqrt( sum( [ (i-j)*(i-j) for i,j in zip( ci, cj ) ] ) ) )
    
    @staticmethod
    def __cosang( ci, cj, ck ):
        vij = [ (j-i) for i,j in zip( ci, cj ) ]
        mij = sum( [ i*i for i in vij ] )
        vik = [ (k-i) for i,k in zip( ci, ck ) ]
        mik = sum( [ i*i for i in vik ] )
        return( sum( [ i*j for i,j in zip( vij, vik ) ] ) / math.sqrt( mij * mik ) )

    def get_info( self, molec ):
        crd = []
        for i in self.sele:
            crd += molec.coor[i:i+3][:]
        if( ml_info_so ):
            return( qm3.engines._ml_info.acsf_info( self.cutx, self.eta2, self.dse5, self.pre5, self.eta5, crd ) )
        else:
            siz = len( self.sele )
            out = [ 0.0 for i in range( siz * self.size ) ]
            for i in range( siz ):
                i3 = self.sele[i]
                ci = molec.coor[i3:i3+3][:]
                dd = i * self.size
                for j in range( siz ):
                    if( j != i ):
                        j3 = self.sele[j]
                        cj = molec.coor[j3:j3+3][:]
                        dij = self.__dist( ci, cj )
                        fij = self.__fcut( dij, self.cutx )
                        if( fij > 0.0 ):
                            ww = 0
                            for ee in self.eta2:
                                out[dd+ww] += fij * math.exp( - ee * dij * dij )
                                ww += 1
                            for k in range( siz ):
                                if( k != j and k != i ):
                                    k3 = self.sele[k]
                                    ck = molec.coor[k3:k3+3][:]
                                    dik = self.__dist( ci, ck )
                                    fik = self.__fcut( dik, self.cutx )
                                    if( fik > 0.0 ):
                                        ww = len( self.eta2 )
                                        for ee in self.eta5:
                                            out[dd+ww] += self.pre5 * fij * fik * \
                                                math.pow( 1.0 + self.__cosang( ci, cj, ck ), self.dse5 ) * \
                                                math.exp( - ee * ( dij * dij + dik * dik ) )
                                            ww += 1
            return( out )

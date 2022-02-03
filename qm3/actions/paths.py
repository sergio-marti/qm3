# -*- coding: iso-8859-1 -*-
import sys
import math
import qm3.maths.matrix



__vcut = 0.00035481432270250985



def default_log( txt ):
    sys.stdout.write( txt + "\n" )
    sys.stdout.flush()



# mass weighted:  xyz * sqrt(m)  ;  grd / sqrt(m)  ;  hes / sqrt(mi * mj)
def __project_RT_modes( mas, crd, grd, hes = [] ):
    siz = len( crd )
    mtt = 0.0
    cen = [ 0.0, 0.0, 0.0 ]
    for i in range( siz // 3 ):
        k = i * 3
        mtt += mas[k] * mas[k]
        for j in [0, 1, 2]:
            cen[j] += crd[k+j] * mas[k]
    cen = [ cen[i] / mtt for i in [0, 1, 2] ]
    rtm = [ 0.0 for i in range( 6 * siz ) ]
    for i in range( siz // 3 ):
        k              = i * 3
        rtm[k]         = mas[k]
        rtm[siz+k+1]   = mas[k]
        rtm[2*siz+k+2] = mas[k]
        rtm[3*siz+k+1] = - ( crd[k+2] - cen[2] * mas[k] )
        rtm[3*siz+k+2] =   ( crd[k+1] - cen[1] * mas[k] )
        rtm[4*siz+k]   =   ( crd[k+2] - cen[2] * mas[k] )
        rtm[4*siz+k+2] = - ( crd[k  ] - cen[0] * mas[k] )
        rtm[5*siz+k]   = - ( crd[k+1] - cen[1] * mas[k] )
        rtm[5*siz+k+1] =   ( crd[k  ] - cen[0] * mas[k] )
    for i in range( 6 ):
        for j in range( i ):
            tmp = sum( [ rtm[i*siz+k] * rtm[j*siz+k] for k in range( siz ) ] )
            for k in range( siz ):
                rtm[i*siz+k] -= tmp * rtm[j*siz+k]
        tmp = math.sqrt( sum( [ rtm[i*siz+k] * rtm[i*siz+k] for k in range( siz ) ] ) )
        for k in range( siz ):
            rtm[i*siz+k] /= tmp
    # gradient
    for i in range( 6 ):
        tmp = sum( [ grd[k] * rtm[i*siz+k] for k in range( siz ) ] )
        for k in range( siz ):
            grd[k] -= tmp * rtm[i*siz+k]
    # hessian
    if( hes != [] ):
        ixx = [ 0.0 for i in range( siz * siz ) ]
        for i in range( siz ):
            ixx[siz*i+i] += 1.
            for j in range( siz ):
                for k in range( 6 ):
                    ixx[siz*i+j] -= rtm[k*siz+i] * rtm[k*siz+j]
        tmp = qm3.maths.matrix.mult( ixx, siz, siz, qm3.maths.matrix.mult( hes, siz, siz, ixx, siz, siz ), siz, siz )
        for i in range( siz * siz ):
            hes[i] = tmp[i]



# use positive/forward or negative/reverse
def initial_step( obj, step_size = 0.0053, project_RT = True ):
    s = min( len( obj.mass ), obj.size )
    k = obj.size // s
    w = [ 0.0 for i in range( obj.size ) ]
    for i in range( s ):
        w[i*k] = math.sqrt( obj.mass[i] )
        for j in range( 1, k ):
            w[i*k+j] = w[i*k]
    x = [ obj.coor[i] * w[i] for i in range( obj.size ) ]
    obj.get_hess()
    h = []
    k = 0
    for i in range( obj.size ):
        for j in range( obj.size ):
            h.append( obj.hess[k] / ( w[i] * w[j] ) )
            k += 1
    g = [ obj.grad[i] / w[i] for i in range( obj.size ) ]
    if( project_RT ):
        __project_RT_modes( w, x, g, h )
    val, vec = qm3.maths.matrix.diag( h, obj.size )
    nskp = sum( [ 1 for i in range( obj.size ) if val[i] < __vcut ] )
    print( "initial_step:", val[0:10] )
    tmp  = step_size / math.sqrt( sum( [ vec[i*obj.size] * vec[i*obj.size] for i in range( obj.size ) ] ) )
    dx   = [ vec[i*obj.size] * tmp for i in range( obj.size ) ]
    return( nskp, dx, [ dx[i] / w[i] for i in range( obj.size ) ] )



def euler( obj, 
            step_number = 100,
            step_size = 0.0053,            # use positive/forward or negative/reverse
            gradient_tolerance = 1.0,
            print_frequency = 10,
            project_RT = True,
            from_saddle = True,
            estabilize = True,
            log_function = default_log ):
    log_function( "\n---------------------------------------- Minimum Path (Euler)\n" )
    log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
    log_function( "Step Number:        %20d"%( step_number ) )
    log_function( "Step Size:          %20.10lg"%( step_size ) )
    log_function( "Print Frequency:    %20d"%( print_frequency ) )
    log_function( "Gradient Tolerance: %20.10lg"%( gradient_tolerance ) )
    log_function( "Project RT modes:   %20s"%( project_RT ) )
    log_function( "From Saddle:        %20s"%( from_saddle ) )
    log_function( "%10s%25s%25s"%( "Step", "Function", "Gradient" ) )
    log_function( "-" * 60 )
    s    = min( len( obj.mass ), obj.size )
    k    = obj.size // s
    w    = [ 0.0 for i in range( obj.size ) ]
    for i in range( s ):
        w[i*k] = math.sqrt( obj.mass[i] )
        for j in range( 1, k ):
            w[i*k+j] = w[i*k]
    dx = [ 0.0 for i in range( obj.size ) ]
    x  = [ obj.coor[i] * w[i] for i in range( obj.size ) ]
    if( from_saddle ):
        nskp, dx, tt = initial_step( obj, step_size, project_RT )
    step_size = math.fabs( step_size )
    grms      = gradient_tolerance * 2.0
    it1       = 0
    it2       = step_number // 10
    # -- the minimal amount of iterations sholud be tunned (no information about the topology: nskp)
    while( it1 < step_number and ( grms > gradient_tolerance or it1 < it2 ) ):
        for i in range( obj.size ):
            x[i] += dx[i]
            obj.coor[i] = x[i] / w[i]
        obj.current_step( it1 )
        obj.get_grad()
        gn = [ obj.grad[i] / w[i] for i in range( obj.size ) ]
        if( project_RT ):
            __project_RT_modes( w, x, gn )
        gg   = math.sqrt( sum( [ i * i for i in gn ] ) )
        dx   = [ - step_size * gn[i] / gg for i in range( obj.size ) ]
        grms = math.sqrt( sum( [ i * i for i in obj.grad ] ) / float( obj.size ) )
        it1 += 1
        if( it1%print_frequency == 0 ):
            log_function( "%10ld%25.5lf%25.10lf"%( it1, obj.func, grms ) )
    if( it1%print_frequency != 0 ):
        log_function( "%10ld%25.5lf%25.10lf"%( it1, obj.func, grms ) )
    log_function( "-" * 60 + "\n" )



def taylor( obj, 
            step_number = 100,
            step_size = 0.0053,            # use positive/forward or negative/reverse
            gradient_tolerance = 1,
            print_frequency = 10,
            project_RT = True,
            from_saddle = True,
            avoid_recrossing = 50,
            log_function = default_log ):
    log_function( "\n---------------------------------------- Minimum Path (Taylor)\n" )
    log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
    log_function( "Step Number:        %20d"%( step_number ) )
    log_function( "Step Size:          %20.10lg"%( step_size ) )
    log_function( "Print Frequency:    %20d"%( print_frequency ) )
    log_function( "Gradient Tolerance: %20.10lg"%( gradient_tolerance ) )
    log_function( "Project RT modes:   %20s"%( project_RT ) )
    log_function( "From Saddle:        %20s"%( from_saddle ) )
    log_function( "%10s%25s%25s"%( "Step", "Function", "Gradient" ) )
    log_function( "-" * 60 )
    s = min( len( obj.mass ), obj.size )
    k = obj.size // s
    w = [ 0.0 for i in range( obj.size ) ]
    for i in range( s ):
        w[i*k] = math.sqrt( obj.mass[i] )
        for j in range( 1, k ):
            w[i*k+j] = w[i*k]
    dx = [ 0.0 for i in range( obj.size ) ]
    v  = [ 0.0 for i in range( obj.size ) ]
    x  = [ obj.coor[i] * w[i] for i in range( obj.size ) ]
    if( from_saddle ):
        nskp, dx, tt = initial_step( obj, step_size, project_RT )
        if( avoid_recrossing > 0 ):
            ox   = dx[:]
    step_size = math.fabs( step_size )
    grms      = gradient_tolerance * 2.0
    it1       = 0
    # -- the minimal amount of iterations sholud be tunned (no information about the topology: nskp)
    itm       = step_number // 10
    while( it1 < step_number and ( grms > gradient_tolerance or it1 < itm ) ):
        for i in range( obj.size ):
            x[i] += dx[i]
            obj.coor[i] = x[i] / w[i]
        obj.current_step( it1 )
        obj.get_hess()
        h = []
        k = 0
        for i in range( obj.size ):
            for j in range( obj.size ):
                h.append( obj.hess[k] / ( w[i] * w[j] ) )
                k += 1
        g = [ obj.grad[i] / w[i] for i in range( obj.size ) ]
        if( project_RT ):
            __project_RT_modes( w, x, g, h )
        gg = math.sqrt( sum( [ i * i for i in g ] ) )
        # Eqs 2, 4, 7 & 13 of J. Chem. Phys. v88, p922 (1988) [10.1063/1.454172]
        v0 = [ - g[i] / gg for i in range( obj.size ) ]
        tt = qm3.maths.matrix.mult( h, obj.size, obj.size, v0, obj.size, 1 )
        pp = sum( [ i * j for i,j in zip( v0, tt ) ] )
        dx = []
        for i in range( obj.size ):
            v1 = ( tt[i] - pp * v0[i] ) / gg
            dx.append( step_size * ( v0[i] + 0.5 * step_size * v1 ) )
        # avoid recrossing
        if( from_saddle and it1 <= avoid_recrossing and avoid_recrossing > 0 ):
            tmp = sum( [ ox[i] * dx[i] for i in range( obj.size ) ] )
            if( tmp < 0.0 ):
                dx = [ -dx[i] for i in range( obj.size ) ]
        grms = math.sqrt( sum( [ i * i for i in obj.grad ] ) / float( obj.size ) )
        it1 += 1
        if( it1%print_frequency == 0 ):
            log_function( "%10ld%25.5lf%25.10lf"%( it1, obj.func, grms ) )
    if( it1%print_frequency != 0 ):
        log_function( "%10ld%25.5lf%25.10lf"%( it1, obj.func, grms ) )
    log_function( "-" * 60 + "\n" )



def baker( obj, 
            step_number = 100,
            step_size = 0.0053,            # use positive/forward or negative/reverse
            gradient_tolerance = 1,
            print_frequency = 10,
            project_RT = True,
            from_saddle = True,
            avoid_recrossing = 50,
            log_function = default_log ):
    log_function( "\n---------------------------------------- Minimum Path (Baker)\n" )
    log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
    log_function( "Step Number:        %20d"%( step_number ) )
    log_function( "Step Size:          %20.10lg"%( step_size ) )
    log_function( "Print Frequency:    %20d"%( print_frequency ) )
    log_function( "Gradient Tolerance: %20.10lg"%( gradient_tolerance ) )
    log_function( "Project RT modes:   %20s"%( project_RT ) )
    log_function( "From Saddle:        %20s"%( from_saddle ) )
    log_function( "Avoid Recrossing:   %20d\n"%( avoid_recrossing ) )
    log_function( "%10s%20s%20s%10s"%( "Step", "Function", "Gradient", "Nskip" ) )
    log_function( "-" * 60 )
    lrge = 1.0e+6
    step = 50.0
    tol2 = 1.0e-8
    mxit = 999
    s    = min( len( obj.mass ), obj.size )
    k    = obj.size // s
    w    = [ 0.0 for i in range( obj.size ) ]
    for i in range( s ):
        w[i*k] = math.sqrt( obj.mass[i] )
        for j in range( 1, k ):
            w[i*k+j] = w[i*k]
    dx = [ 0.0 for i in range( obj.size ) ]
    gx = [ 0.0 for i in range( obj.size ) ]
    x  = [ obj.coor[i] * w[i] for i in range( obj.size ) ]
    if( from_saddle ):
        nskp, dx, tt = initial_step( obj, step_size, project_RT )
        if( avoid_recrossing > 0 ):
            ox   = dx[:]
    else:
        nskp = 7
    step_size = math.fabs( step_size )
    mskp      = 6 * project_RT
    grms      = gradient_tolerance * 2.0
    it1       = 0
    flg       = True
    while( it1 < step_number and ( grms > gradient_tolerance or nskp > mskp ) and flg ):
        for i in range( obj.size ):
            x[i] += dx[i]
            obj.coor[i] = x[i] / w[i]
        obj.current_step( it1 )
        obj.get_hess()
        h = []
        k = 0
        for i in range( obj.size ):
            for j in range( obj.size ):
                h.append( obj.hess[k] / ( w[i] * w[j] ) )
                k += 1
        g = [ obj.grad[i] / w[i] for i in range( obj.size ) ]
        if( project_RT ):
            __project_RT_modes( w, x, g, h )
        val, vec = qm3.maths.matrix.diag( h, obj.size )
        nskp = sum( [ 1 for i in range( obj.size ) if val[i] < __vcut ] )
        grms = math.sqrt( sum( [ i * i for i in g ] ) )
        # transform gradient vector to the local hessian modes
        for i in range( obj.size ):
            dx[i]  = 0.0
            val[i] /= grms
            gx[i]  = sum( [ g[j] * vec[j*obj.size+i] for j in range( obj.size ) ] )
        # minimize along the selected modes (skip first)
        lmbd = 0.0
        if( val[0] < 0.0 ):
            lmbd = val[0] - step
            l1   = val[0]
            l2   = - lrge
        ovr = sum( [ gx[j] * gx[j] / ( lmbd - val[j] ) for j in range( obj.size ) ] )
        i   = 0
        while( i < mxit and math.fabs( lmbd - ovr ) >= tol2 ):
            if( val[0] > 0.0 ):
                lmbd = ovr;
            else:
                if( ovr < lmbd ):
                    l1 = lmbd;
                if( ovr > lmbd ):
                    l2 = lmbd;
                if( l2 > - lrge ):
                    lmbd = 0.5 * ( l1 + l2 )
                elif( l2 == - lrge ):
                    lmbd -= step;
            ovr = sum( [ gx[j] * gx[j] / ( lmbd - val[j] ) for j in range( obj.size ) ] )
            i += 1
        if( i > mxit ):
            log_function( "\n -- Too much lambda iterations..." )
            flg = False
        for i in range( obj.size ):
            dx[i] += sum( [ vec[i*obj.size+j] * gx[j] / ( lmbd - val[j] ) for j in range( obj.size ) ] )
        # check final step (too small or large...)
        ovr = math.sqrt( sum( [ dx[i] * dx[i] for i in range( obj.size ) ] ) )
        if( ovr < tol2 ):
            log_function( "\n -- The step size is *very* small..." )
            flg = False
        # scale long steps...
        if( ovr > step_size ):
            for i in range( obj.size ):
                dx[i] *= step_size / ovr
        # avoid recrossing
        if( from_saddle and it1 <= avoid_recrossing and nskp > mskp and avoid_recrossing > 0 ):
            tmp = sum( [ ox[i] * dx[i] for i in range( obj.size ) ] )
            if( tmp < 0.0 ):
                dx = [ -dx[i] for i in range( obj.size ) ]
        grms = math.sqrt( sum( [ i * i for i in obj.grad ] ) / float( obj.size ) )
        it1 += 1
        if( it1%print_frequency == 0 ):
            log_function( "%10ld%20.5lf%20.10lf%10ld"%( it1, obj.func, grms, nskp ) )
    if( it1%print_frequency != 0 ):
        log_function( "%10ld%20.5lf%20.10lf%10ld"%( it1, obj.func, grms, nskp ) )
    log_function( "-" * 60 + "\n" )



def page_mciver( obj, 
            step_number = 100,
            step_size = 0.0053,            # use positive/forward or negative/reverse
            gradient_tolerance = 1,
            print_frequency = 10,
            project_RT = True,
            from_saddle = True,
            avoid_recrossing = 50,
            log_function = default_log ):
    log_function( "\n---------------------------------------- Minimum Path (Page-McIver:LQA)\n" )
    log_function( "Degrees of Freedom: %20ld"%( obj.size ) )
    log_function( "Step Number:        %20d"%( step_number ) )
    log_function( "Step Size:          %20.10lg"%( step_size ) )
    log_function( "Print Frequency:    %20d"%( print_frequency ) )
    log_function( "Gradient Tolerance: %20.10lg"%( gradient_tolerance ) )
    log_function( "Project RT modes:   %20s"%( project_RT ) )
    log_function( "From Saddle:        %20s"%( from_saddle ) )
    log_function( "Avoid Recrossing:   %20d\n"%( avoid_recrossing ) )
    log_function( "%10s%20s%20s%10s"%( "Step", "Function", "Gradient", "Nskip" ) )
    log_function( "-" * 60 )
    it2m = 1000
    it3m = 100000
    s    = min( len( obj.mass ), obj.size )
    k    = obj.size // s
    w    = [ 0.0 for i in range( obj.size ) ]
    for i in range( s ):
        w[i*k] = math.sqrt( obj.mass[i] )
        for j in range( 1, k ):
            w[i*k+j] = w[i*k]
    dx = [ 0.0 for i in range( obj.size ) ]
    v  = [ 0.0 for i in range( obj.size ) ]
    x  = [ obj.coor[i] * w[i] for i in range( obj.size ) ]
    if( from_saddle ):
        nskp, dx, tt = initial_step( obj, step_size, project_RT )
        if( avoid_recrossing > 0 ):
            ox   = dx[:]
    else:
        nskp = 7
    step_size = math.fabs( step_size )
    mskp      = 6 * project_RT
    grms      = gradient_tolerance * 2.0
    it1       = 0
    flg       = True
    while( it1 < step_number and ( grms > gradient_tolerance or nskp > mskp ) and flg ):
        for i in range( obj.size ):
            x[i] += dx[i]
            obj.coor[i] = x[i] / w[i]
        obj.current_step( it1 )
        obj.get_hess()
        h = []
        k = 0
        for i in range( obj.size ):
            for j in range( obj.size ):
                h.append( obj.hess[k] / ( w[i] * w[j] ) )
                k += 1
        g = [ obj.grad[i] / w[i] for i in range( obj.size ) ]
        if( project_RT ):
            __project_RT_modes( w, x, g, h )
        val, vec = qm3.maths.matrix.diag( h, obj.size )
        nskp = sum( [ 1 for i in range( obj.size ) if val[i] < __vcut ] )
        grms = math.sqrt( sum( [ i * i for i in g ] ) )
        # transform gradient vector to the local hessian modes
        for i in range( obj.size ):
            g[i]   /= grms
            val[i] /= grms
        for i in range( obj.size ):
            v[i]  = sum( [ g[j] * vec[j*obj.size+i] for j in range( obj.size ) ] )
        # search for the PM step
        pm_dt = 0.2 * step_size
        pm_t  = 1.e10
        pm_ot = 0.0
        it2   = 0
        it3   = 0
        while( math.fabs( 1.0 - pm_ot / pm_t ) > 1.0e-6 and it2 < it2m and it3 < it3m ):
            it2  += 1
            pm_ot = pm_t
            pm_dt = 0.5 * pm_dt
            pm_t  = 0.0
            pm_ft = math.sqrt( sum( [ math.pow( v[i] * math.exp( - val[i] * pm_t ), 2.0 ) for i in range( obj.size ) ] ) )
            pm_s  = 0.0;
            it3   = 0;
            while( pm_s < step_size and it3 < it3m ):
                it3  += 1
                pm_os = pm_s
                pm_of = pm_ft
                pm_t += pm_dt
                pm_ft = math.sqrt( sum( [ math.pow( v[i] * math.exp( - val[i] * pm_t ), 2.0 ) for i in range( obj.size ) ] ) )
                pm_s += 0.5 * pm_dt * ( pm_ft + pm_of )
            # does not converge if it cannot reach the minimum with the current step-size (we have ended...)
            if( pm_os != pm_s ):
                pm_t -= ( step_size - pm_s ) * pm_dt / ( pm_os - pm_s )
            else:
                log_function( "\n -- The current step-size did not converge..." )
                flg = False
        if( math.fabs( 1.0 - pm_t / pm_ot ) <= 1.0e-6 and flg ):
            for i in range( obj.size ):
                tmp = val[i] * pm_t
                if( math.fabs( tmp ) < 1.e-8 ):
                    v[i] *= - pm_t * ( 1.0 - tmp  * ( 0.5 - tmp / 6.0 ) )
                else:
                    v[i] *= ( math.exp( - tmp ) - 1.0 ) / val[i]
            for i in range( obj.size ):
                dx[i] = sum( [ v[j] * vec[i*obj.size+j] for j in range( obj.size ) ] )
            # avoid recrossing
            if( from_saddle and it1 <= avoid_recrossing and nskp > mskp and avoid_recrossing > 0 ):
                tmp = sum( [ ox[i] * dx[i] for i in range( obj.size ) ] )
                if( tmp < 0.0 ):
                    dx = [ -dx[i] for i in range( obj.size ) ]
        grms = math.sqrt( sum( [ i * i for i in obj.grad ] ) / float( obj.size ) )
        it1 += 1
        if( it1%print_frequency == 0 ):
            log_function( "%10ld%20.5lf%20.10lf%10ld"%( it1, obj.func, grms, nskp ) )
    if( it1%print_frequency != 0 ):
        log_function( "%10ld%20.5lf%20.10lf%10ld"%( it1, obj.func, grms, nskp ) )
    log_function( "-" * 60 + "\n" )




#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdio.h>


long xilaenv_( long* );
long ilaenv_( long*, char*, char*, long*, long*, long*, long* );
void dgetrf_( long*, long*, double*, long*, long*, long* );
void dgetri_( long*, double*, long*, long*, double*, long*, long* );
void dspev_ ( char*, char*, long*, double*, double*, double*, long*, double*, long* );
void dgemm_ ( char*, char*, long*, long*, long*, double*, double*, long*, double*, long*, double*, double*, long* );


double __dot_product( long n, double *v1, double *v2 ) {
    double	o = 0.0;
    long	i;
    for( i = 0; i < n; i++ ) o += v1[i] * v2[i];
    return( o );
}



static void __transpose( double *mat, long n, long m ) {
    double	*t, s;
    long	i, j, k;

    if( n == m ) {
    	for( i=0; i<n-1; i++ )
    		for( j=i+1; j<n; j++ ) {
    			s = mat[i*n+j];
    			mat[i*n+j] = mat[j*n+i];
    			mat[j*n+i] = s;
    		}
    } else {
    	t = (double*)malloc((n*m)*sizeof(double));
    	k = 0;
    	for( j=0; j<m; j++) for( i=0; i<n; i++ ) t[k++] = mat[i*m+j];
    	for( i=0; i<n*m; i++ ) mat[i] = t[i];
    	free( t );
    }
}



static void __mult( double *a_mat, long a_n, long a_m, double *b_mat, long b_n, long b_m, double *c_mat ) {
/*
    double	t;
    long	i, j, k;

    if( a_m == b_n ) {

    	for( i=0; i<a_n; i++ )
    		for( j=0; j<b_m; j++ ) {
    			t = .0;
    			for( k=0; k<a_m; k++ )
    				t += a_mat[i*a_m+k] * b_mat[k*b_m+j];
    			c_mat[i*b_m+j] = t;
    		}
    }

*/
/*
SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
M    	specifies the number of rows of the matrix op( A ) and of the  matrix  C
N    	specifies the number of columns of the matrix op( B ) and the number of columns of the matrix C
K    	specifies the number of columns of the matrix op( A ) and the number of rows of the matrix op( B )
LDA    	specifies the first dimension of A as declared in the calling (sub) program. When TRANSA = 'N' or 'n' then
    	LDA must be at least max( 1, m ), otherwise LDA must be at least max( 1, k ).
LDB    	specifies the first dimension of B as declared in the calling (sub) program. When TRANSB = 'N' or 'n' then
    	LDB must be at least max( 1, k ), otherwise  LDB must be at least  max( 1, n ).
LDC    	specifies the first dimension of C as declared in the calling (sub) program. LDC must be at least max( 1, m ).
*/

    char	ca, cb;
    double	ta, tb;
    long	i, m, n, k, la, lb, lc;

    if( a_m == b_n ) {
    	for( i = 0; i < a_n * b_m; i++ ) c_mat[i] = .0;
    	// B matrix (C >> Fortran): b_n x b_m
    	ca = 'N'; m = b_m;
    	// A matrix (C >> Fortran): a_n x a_m
    	cb = 'N'; n = a_n;
    	// Common: a_m == b_n
        k  = a_m; // = b_n;
    	la = m; lb = k; lc = m;
    	ta = 1.0; tb = 0.0;

    	dgemm_( &ca, &cb, &m, &n, &k, &ta, b_mat, &la, a_mat, &lb, &tb, c_mat, &lc );
    }

}



double __det( double *mat, long n ) {
    double	o, sgn, *tmp;
    long	i, j, k, l;

    o = .0;
    if( n == 2 ) {
    	o = mat[0] * mat[3] - mat[1] * mat[2];
    } else if( n == 3 ) {
    	o = mat[0] * mat[4] * mat[8] + mat[2] * mat[3] * mat[7] + mat[1] * mat[5] * mat[6]-
    		mat[2] * mat[4] * mat[6] - mat[1] * mat[3] * mat[8] - mat[0] * mat[5] * mat[7];
    } else {
    	sgn = 1.;
    	for( i=0; i<n; i++ ) {
    		if( mat[i] != .0 ) {
    			l = 0;
    			tmp = (double*)malloc(((n-1)*(n-1))*sizeof(double));
    			for( j=1; j<n; j++ )
    				for( k=0; k<n; k++ ) 
    					if( k != i ) 
    						tmp[l++] = mat[j*n+k];
    			o += sgn * mat[i] * __det( tmp, n-1 );
    			free( tmp );
    		}
    		sgn *= -1.;
    	}
    }
    return( o );
}



static void __invert( double *mat, long n ) {
    double	*dge_work;
    long	dge_n, dge_lwork, *dge_ipiv, dge_info;
//    long	nb, ila_n1;
    long	nb, ila_ispec, ila_n1, ila_n2, ila_n3, ila_n4;
    char	*ila_name, *ila_opts;

    __transpose( mat, n, n );
//    ila_n1 = n;
    ila_ispec = 1;
    ila_n2 = -1;
    ila_n3 = -1;
    ila_n4 = -1;
    ila_name = (char*) malloc( 10 * sizeof( char ) );
    ila_opts = (char*) malloc( 10 * sizeof( char ) );
    snprintf( ila_name, 10, "DGETRI" );
    snprintf( ila_opts, 10, " " );
//    nb = ilaenv_( &ila_ispec, ila_name, ila_opts, &ila_n1, &ila_n2, &ila_n3, &ila_n4 );
    nb = xilaenv_( &ila_n1 );
    dge_n = n;
    dge_lwork = dge_n * nb;
    dge_work = (double*) malloc(dge_lwork*sizeof(double));
    dge_ipiv = (long*) malloc(dge_n*sizeof(long));
    dgetrf_( &dge_n, &dge_n, mat, &dge_n, dge_ipiv, &dge_info );
    dgetri_( &dge_n, mat, &dge_n, dge_ipiv, dge_work, &dge_lwork, &dge_info );
    free( dge_ipiv ); 
    free( dge_work ); 
	free( ila_name );
	free( ila_opts );
    __transpose( mat, n, n );
}



static void __pseudoinvert( double *mat, long n, long m ) {
    double	*Tmat, *tmp;
    long	i;

    if( n > m ) {
    	// A^{{}+{}}_{m,n} = left( A^T_{m,n} * A_{n,m} right)^{-1}_{m,m} * A^T_{m,n}
    	// A^{{}+{}}_{m,n} * A_{n,m} = I_{m,m}
    	Tmat=(double*)malloc((n*m)*sizeof(double));
    	tmp=(double*)malloc((m*m)*sizeof(double));
    	for( i=0; i<n*m; i++ ) Tmat[i] = mat[i];
    	__transpose( Tmat, n, m );
    	__mult( Tmat, m, n, mat, n, m, tmp );
    	__invert( tmp, m );
    	__mult( tmp, m, m, Tmat, m, n, mat );
    } else {
    	// A^{{}+{}}_{m,n} = A^T_{m,n} * left( A_{n,m} * A^T_{m,n} right)^{-1}_{n,n}
    	// A_{n,m} * A^{{}+{}}_{m,n} = I_{n,n}
    	Tmat=(double*)malloc((n*m)*sizeof(double));
    	tmp=(double*)malloc((n*n)*sizeof(double));
    	for( i=0; i<n*m; i++ ) Tmat[i] = mat[i];
    	__transpose( Tmat, n, m );
    	__mult( mat, n, m, Tmat, m, n, tmp );
    	__invert( tmp, n );
    	__mult( Tmat, m, n, tmp, n, n, mat );
    }
    free( tmp ); free( Tmat );
}



static void __diag( double *mat, double *evec, double *eval, long n ) {
    long	f, i, j, k;
    double	*work, *udcp;

    udcp = (double*)malloc(n*(n+1)/2*sizeof(double));
    k = 0;
    for( i=0; i<n; i++ )
    	for( j=0; j<=i; j++ )
    		udcp[k++] = mat[i+j*n];
    work = (double*)malloc(3*n*sizeof(double));
    dspev_( "V", "U", &n, udcp, eval, evec, &n, work, &f );
    free( udcp ); free( work );
    __transpose( evec, n, n );
}



long __jacobi( double *mat, double *vec, long n ) {
#define JACOBI_MAX  10000
#define JACOBI_CUT  1.e-10
    double	*t, mx, mk, vtet, vsen, vcos, ds;
    long	i, j, k, nit, mi, mj, *idx;

    for( i=0; i<n; i++ ) {
    	for( j=0; j<n; j++ )
    		vec[i*n+j] = .0;
    	vec[i*n+i] = 1.;
    }

    t = (double*)malloc(2*n*sizeof(double));
    idx = (long*)malloc(2*n*sizeof(long));

    nit = 0;
    while( nit < JACOBI_MAX ) {
    	//--------------------------------------------------------------- Piv
    	mi = 0; 
    	mj = 1; 
    	mx = fabs( mat[1] );
    	ds = .0;
    	for( i=0; i<n-1; i++ )
    		for( j=i+1; j<n; j++ ) {
    			mk = fabs( mat[i*n+j] );
    			ds += mk;
    			if( mk > mx ) {
    				mi = i;
    				mj = j;
    				mx = mk;
    			}
    		}
    	if( ds < JACOBI_CUT ) break;
    	//--------------------------------------------------------------- Ang
    	if( mat[mj*(n+1)] != mat[mi*(n+1)] ) {
    		vtet = .5 * atan( -2. * mat[mi*n+mj] / ( mat[mj*(n+1)] - mat[mi*(n+1)] ) );
    	} else { vtet = 0.250 * M_PI; }
    	vsen = sin( vtet );
    	vcos = cos( vtet );
    	//--------------------------------------------------------------- Rot
    	for( k=0; k<n; k++ ) {
    		if( k <= mi ) idx[k] = k * n + mi;
    		else idx[k] = mi * n + k;

    		if( k <= mj ) idx[n+k] = k * n + mj;
    		else idx[n+k] = mj * n + k;

    		if( k != mi && k != mj ) {
    			t[k] = vcos * mat[k*n+mi] + vsen * mat[k*n+mj];
    			t[n+k] = vcos * mat[k*n+mj] - vsen * mat[k*n+mi];
    		}
    	}
    	t[mi] = vcos * vcos * mat[mi*n+mi]
    		+ vsen * vsen * mat[mj*n+mj]
    		+ 2. * vsen * vcos * mat[mi*n+mj];
    	t[n+mj] = vsen * vsen * mat[mi*n+mi] 
    		+ vcos * vcos * mat[mj*n+mj]
    		- 2. * vsen * vcos * mat[mi*n+mj];
    	t[n+mi] = ( vcos * vcos - vsen * vsen ) * mat[mi*n+mj]
    		+ vsen * vcos * ( mat[mj*n+mj] - mat[mi*n+mi] );
    	t[mj] = t[n+mi];
    	for( k=0; k<2*n; k++) mat[idx[k]] = t[k];
    	for( i=0; i<n-1; i++ )
    		for( j=i+1; j<n; j++ )
    			mat[j*n+i] = mat[i*n+j];
    	//--------------------------------------------------------------- Eig
    	for( k=0; k<n; k++ ) {
    		i = k * n + mi;
    		j = k * n + mj;
    		t[0] = vcos * vec[i] + vsen * vec[j];
    		t[1] = vcos * vec[j] - vsen * vec[i];
    		vec[i] = t[0];
    		vec[j] = t[1];
    	}
    	nit++;
    }
    free( idx ); free( t );
    return( ds < JACOBI_CUT );
}



void __fjacobi_find_max( double *mat, long n, long *p, long *q ) {
    double	t = 0.0;
    long	i, j;

    for( i=0; i<n; i++ )
    	for( j=0; j < n; j++ )
    		if( i != j && fabs( mat[i*n+j] ) > t ) {
    			t = fabs( mat[i*n+j] );
    			*p = i; *q = j;
    		}
}
void __fjacobi_fast_rotation( double *mat, long n, double si, double co, long p, long q ) {
    double	ap, aq;
    long	i;

    for( i=0; i<n; i++ ) {
    	ap = mat[i*n+p];
    	aq = mat[i*n+q];
    	mat[i*n+p] = co * ap - si * aq;
    	mat[i*n+q] = si * ap + co * aq;
    }
    for( i=0; i<n; i++ ) {
    	ap = mat[p*n+i];
    	aq = mat[q*n+i];
    	mat[p*n+i] = co * ap - si * aq;
    	mat[q*n+i] = si * ap + co * aq;
    }
}
long __fjacobi( double *mat, double *vec, long n ) {
#define FJACOBI_MAX  100000
#define FJACOBI_CUT  1.e-10

    double		si, co, *w, *v, dpq, t, x;
    long		i, j, p, q, c;

    w = (double*)malloc( n*n * sizeof( double ) );
    v = (double*)malloc( n*n * sizeof( double ) );
    for( i=0; i<n; i++ ) for( j=0; j<n; j++ ) if( i == j ) { vec[i*n+j] = 1.0; } else { vec[i*n+j] = 0.0; }
    
    x = FJACOBI_CUT;
    c = 0;
    while( c < FJACOBI_MAX && x >= FJACOBI_CUT ) {
    	__fjacobi_find_max( mat, n, &p, &q );
    	dpq = mat[q*n+q] - mat[p*n+p];
    	if( dpq == 0.0 ) { t = 1.0; } else {
    		x = dpq / mat[p*n+q] * 0.50;
    		t = ( x / fabs( x ) ) / ( fabs( x ) + sqrt( x * x + 1.0 ) );
    	}
    	co = 1.0 / sqrt( t * t + 1.0 );
    	si = t * co;
    	__fjacobi_fast_rotation( mat, n, si, co, p, q );
    	for( i=0; i<n; i++ ) for( j=0; j<n; j++ ) if( i == j ) { w[i*n+j] = 1.0; } else { w[i*n+j] = 0.0; }
    	w[p*n+p] = co;
    	w[p*n+q] = si;
    	w[q*n+p] = - si;
    	w[q*n+q] = co;
    	__mult( vec, n, n, w, n, n, v );
    	for( i=0; i<n*n; i++ ) vec[i] = v[i];
    	x = 0.0; for( i=0; i<n-1; i++ ) for( j=i+1; j<n; j++ ) x+= fabs( mat[i*n+j] );
    	c++;
    }
    free( w ); free( v );
    return( c < FJACOBI_MAX );
}



/* ===========================================================================================================================================

//gfortran -c -O1 ilaenv_fix.f
//gfortran -c -O2 lapack_deps.f
//(gcc -c matrix_core.c)&&(gcc matrix_core.o lapack_deps.o ilaenv_fix.o -L/opt/local/Cellar/gcc/5.2.0/lib/gcc/5 -lm -lgfortran)&&(./a.out)

int main() {

    double		*a_dat, *b_dat, *vec, *val;
    long		a_row, a_col, b_row, b_col;
    long		i, j;

    a_row = 3; a_col = 3;
    a_dat = (double*) malloc( (a_row*a_col) * sizeof( double ) );

    a_dat[0] = 1.; 		a_dat[1] = 2.; 		a_dat[2] = 3.;
    a_dat[3] = 4.; 		a_dat[4] = 5.; 		a_dat[5] = 6.;
    a_dat[6] = 7.; 		a_dat[7] = 8.; 		a_dat[8] = 9.;
printf( "\n[]\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", a_dat[i*a_col+j] ); } printf( "\n" ); }
    __transpose( a_dat, a_row, a_col );
printf( "\n[]^T\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", a_dat[i*a_col+j] ); } printf( "\n" ); }

    a_dat[0] = 1.; 		a_dat[1] = 2.; 		a_dat[2] = 7.;
    a_dat[3] = 2.; 		a_dat[4] = 5.; 		a_dat[5] = 8.;
    a_dat[6] = 7.; 		a_dat[7] = 8.; 		a_dat[8] = 15.;

printf( "\n[]\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", a_dat[i*a_col+j] ); } printf( "\n" ); }
    printf( "\ndet: %8.3lf\n", __det( a_dat, a_row ) );

    vec = (double*) malloc( (a_row*a_col) * sizeof( double ) );
    __jacobi( a_dat, vec, a_row );
printf( "\n(Jacobi) Eig-val\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", a_dat[i*a_col+j] ); } printf( "\n" ); }
printf( "Eig-vec\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", vec[i*a_col+j] ); } printf( "\n" ); }

    a_dat[0] = 1.; 		a_dat[1] = 2.; 		a_dat[2] = 7.;
    a_dat[3] = 2.; 		a_dat[4] = 5.; 		a_dat[5] = 8.;
    a_dat[6] = 7.; 		a_dat[7] = 8.; 		a_dat[8] = 15.;
    __fjacobi( a_dat, vec, a_row );
printf( "\n(fJacobi) Eig-val\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", a_dat[i*a_col+j] ); } printf( "\n" ); }
printf( "Eig-vec\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", vec[i*a_col+j] ); } printf( "\n" ); }

    a_dat[0] = 1.; 		a_dat[1] = 2.; 		a_dat[2] = 7.;
    a_dat[3] = 2.; 		a_dat[4] = 5.; 		a_dat[5] = 8.;
    a_dat[6] = 7.; 		a_dat[7] = 8.; 		a_dat[8] = 15.;
    val = (double*) malloc( a_row * sizeof( double ) );
    __diag( a_dat, vec, val, a_row );
printf( "\n(Lapack) Eig-val\n" ); for( i = 0; i < a_row; i++ ) { printf( "%8.3lf", val[i] ); }
printf( "\nEig-vec\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", vec[i*a_col+j] ); } printf( "\n" ); }

    a_dat[0] = 1.; 		a_dat[1] = 2.; 		a_dat[2] = 7.;
    a_dat[3] = 2.; 		a_dat[4] = 5.; 		a_dat[5] = 8.;
    a_dat[6] = 7.; 		a_dat[7] = 8.; 		a_dat[8] = 15.;

    __invert( a_dat, a_row );
printf( "\n[]^-1\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", a_dat[i*a_col+j] ); } printf( "\n" ); }

    a_dat[0] = 1.; 		a_dat[1] = 2.; 		a_dat[2] = 7.;
    a_dat[3] = 2.; 		a_dat[4] = 5.; 		a_dat[5] = 8.;
    a_dat[6] = 7.; 		a_dat[7] = 8.; 		a_dat[8] = 15.;

    __mult( a_dat, a_row, a_col, a_dat, a_row, a_col, vec );
printf( "\n[]*[]\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", vec[i*a_col+j] ); } printf( "\n" ); }

    free( a_dat );
    free( vec );
    free( val );

    a_row = 3; a_col = 4;
    a_dat = (double*) malloc( (a_row*a_col) * sizeof( double ) );
    a_dat[0] = -9.417; 		a_dat[1] = -0.812; 		a_dat[2] = 2.968;		a_dat[3] = 3.915;
    a_dat[4] = 1.128; 		a_dat[5] = 3.908; 		a_dat[6] = 0.675;		a_dat[7] = 0.801;
    a_dat[8] = 0.217; 		a_dat[9] = -0.265; 		a_dat[10] = 0.858;		a_dat[11] =-0.184;

    __pseudoinvert( a_dat, a_row, a_col );
printf( "\n[]+\n" ); for( i = 0; i < a_row; i++ ) { for( j = 0; j < a_col; j++ ) { printf( "%8.3lf", a_dat[i*a_col+j] ); } printf( "\n" ); }

    free( a_dat );

    return( 0 );
}
=========================================================================================================================================== */

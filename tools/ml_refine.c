#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdio.h>
#include<pthread.h>
#include<time.h>


#define STEP_NUMBER	1000
#define	MAX_TRIES	4


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

int		NCPU;
double	STEP_SIZE;

typedef struct { int row; int col; double *mat; } Tmat;


void dgemm_ ( char*, char*, long*, long*, long*, double*, double*, long*, double*, long*, double*, double*, long* );

void mult( Tmat ma, Tmat mb, Tmat mc ) {
    char	transa, transb;
    double	alpha, beta;
    long	m, n, k, lda, ldb, ldc;
    int		i;

   	for( i = 0; i < mc.row * mc.col; i++ ) mc.mat[i] = 0.0;
	m      = mb.col;
	k      = mb.row; //ma.col;
	n      = ma.row;
	transa = 'N';
	transb = 'N';
	lda    = m;
	ldb    = k;
	ldc    = m;
	alpha  = 1.0;
	beta   = 0.0;
   	dgemm_( &transa, &transb, &m, &n, &k, &alpha, mb.mat, &lda, ma.mat, &ldb, &beta, mc.mat, &ldc );
}

void mult_TN( Tmat ma, Tmat mb, Tmat mc ) {
    char	transa, transb;
    double	alpha, beta;
    long	m, n, k, lda, ldb, ldc;
    int		i;

   	for( i = 0; i < mc.row * mc.col; i++ ) mc.mat[i] = 0.0;
	m      = mb.col;
	k      = mb.row; //ma.row;
	n      = ma.col;
	transa = 'N';
	transb = 'T';
	lda    = m;
	ldb    = n;
	ldc    = m;
	alpha  = 1.0;
	beta   = 0.0;
   	dgemm_( &transa, &transb, &m, &n, &k, &alpha, mb.mat, &lda, ma.mat, &ldb, &beta, mc.mat, &ldc );
}

void mult_NT( Tmat ma, Tmat mb, Tmat mc ) {
    char	transa, transb;
    double	alpha, beta;
    long	m, n, k, lda, ldb, ldc;
    int		i;

   	for( i = 0; i < mc.row * mc.col; i++ ) mc.mat[i] = 0.0;
	m      = mb.row;
	k      = mb.col; //ma.col;
	n      = ma.row;
	transa = 'T';
	transb = 'N';
	lda    = k;
	ldb    = k;
	ldc    = m;
	alpha  = 1.0;
	beta   = 0.0;
   	dgemm_( &transa, &transb, &m, &n, &k, &alpha, mb.mat, &lda, ma.mat, &ldb, &beta, mc.mat, &ldc );
}


void f_actv( Tmat mat ) {
	int		i;
	for( i = 0; i < mat.row * mat.col; i++ ) mat.mat[i] = log( 1.0 + exp( min( 709.7827, mat.mat[i] ) ) );
}

void g_actv( Tmat mat ) {
	int		i;
	for( i = 0; i < mat.row * mat.col; i++ ) mat.mat[i] = 1.0 / ( 1.0 + exp( min( 709.7827, - mat.mat[i] ) ) );
}


double predict( int dim, Tmat *mod, Tmat inp ) {
	int		i, j;
	Tmat	out, tmp;
	double	ret;

	out.row = inp.row;
	out.col = inp.col;
	out.mat = (double *) malloc( out.row * out.col * sizeof( double ) );
	for( i = 0; i < out.row * out.col; i++ ) out.mat[i] = inp.mat[i];
	for( i = 0; i < dim; i += 2 ) {
		tmp.row = out.row;
		tmp.col = mod[i].col;
		tmp.mat = (double*) malloc( tmp.row * tmp.col * sizeof( double ) );
		mult( out, mod[i], tmp );
		free( out.mat );
		for( j = 0; j < tmp.col; j++ ) tmp.mat[j] += mod[i+1].mat[j];
		f_actv( tmp );
		out.row = tmp.row;
		out.col = tmp.col;
		out.mat = tmp.mat;
	}
	ret = out.mat[0];
	free( out.mat );
	return( ret );
}


void set_copy( int dim, Tmat *ma, Tmat *mb ) {
    int		i, j;
    for( i = 0; i < dim; i++ ) for( j = 0; j < ma[i].row * ma[i].col; j++ ) mb[i].mat[j] = ma[i].mat[j];
}


void set_zeros( int dim, Tmat *ma ) {
    int		i, j;
    for( i = 0; i < dim; i++ ) for( j = 0; j < ma[i].row * ma[i].col; j++ ) ma[i].mat[j] = 0.0;
}


typedef struct { int dim; Tmat *mod; int ii; int jj; Tmat *inp; double *out; double los; Tmat *grd; } job_arg;

void* worker( void *args ) {
    job_arg		*arg = (job_arg*) args;
	int			i, j, k, l;
	Tmat		*buf, tmp, *fub;
	double		dif;

	buf = (Tmat*) malloc( arg->dim * sizeof( Tmat ) );
	fub = (Tmat*) malloc( arg->dim * sizeof( Tmat ) );
	for( i = 0; i < arg->dim; i++ ) {
		buf[i].row = 1;
		buf[i].col = arg->mod[i].col;
		buf[i].mat = (double*) malloc( buf[i].col * sizeof( double ) );
		fub[i].row = arg->mod[i].row;
		fub[i].col = arg->mod[i].col;
		fub[i].mat = (double*) malloc( fub[i].row * fub[i].col * sizeof( double ) );
	}
	for( i = arg->ii; i < arg->jj; i++ ) {
		// forward
		tmp.row = 1;
		tmp.col = arg->inp[i].col;
		tmp.mat = arg->inp[i].mat;
		for( j = 0; j < arg->dim; j += 2 )  {
			mult( tmp, arg->mod[j], buf[j] );
			for( k = 0; k < buf[j].col; k++ ) {
				buf[j].mat[k] += arg->mod[j+1].mat[k];
				buf[j+1].mat[k] = buf[j].mat[k];
			}
			f_actv( buf[j+1] );
			tmp.col = buf[j+1].col;
			tmp.mat = buf[j+1].mat;
		}
		dif = ( tmp.mat[0] - arg->out[i] );
		arg->los += dif * dif;
		// reverse
		l = arg->dim - 1;
		for( j = 0; j < fub[l].row * fub[l].col; j++ ) fub[l].mat[j] = buf[l-1].mat[j];
		g_actv( fub[l] );
		for( j = 0; j < fub[l].row * fub[l].col; j++ ) fub[l].mat[j] *= 2.0 * dif;
		for( k = arg->dim - 2; k > 0; k -= 2 ) {
			mult_TN( buf[k-1], fub[l], fub[k] );
			mult_NT( fub[l], arg->mod[k], fub[k-1] );
			g_actv( buf[k-2] );
			for( j = 0; j < buf[k-2].col; j++ ) fub[k-1].mat[j] *= buf[k-2].mat[j];
			l -= 2;
		}
		mult_TN( arg->inp[i], fub[1], fub[0] );
		// accumulate
		for( j = 0; j < arg->dim; j++ )
			for( k = 0; k < arg->grd[j].row * arg->grd[j].col; k++ )
				arg->grd[j].mat[k] += fub[j].mat[k];

	}
	for( i = 0; i < arg->dim; i++ ) { free( buf[i].mat ); free( fub[i].mat ); }
	free( buf ); free( fub );
	pthread_exit( NULL );
}


void prm_grad( int dim, Tmat *mod, int nel, Tmat *inp, double *out, double *los, Tmat *grd ) {
	int			i, j, k, *num;
	Tmat		**buf;
    pthread_t	*pid;
    job_arg		*arg;

	buf = (Tmat**) malloc( NCPU * sizeof( Tmat* ) );
	for( i = 0; i < NCPU; i++ ) {
		buf[i] = (Tmat*) malloc( dim * sizeof( Tmat ) );
		for( j = 0; j < dim; j++ ) {
			buf[i][j].row = mod[j].row;
			buf[i][j].col = mod[j].col;
			buf[i][j].mat = (double*) malloc( mod[j].row * mod[j].col * sizeof( double ) );
		}
		set_zeros( dim, buf[i] );
	}
	num = (int*) calloc( NCPU, sizeof( int ) );
	for( i = 0; i < nel; i++ ) num[i%NCPU]++;
   	pid = (pthread_t*) malloc( NCPU * sizeof( pthread_t ) );
   	arg = (job_arg*) malloc( NCPU * sizeof( job_arg ) );
	j   = 0;
    for( i = 0; i < NCPU; i++ ) {
		arg[i].dim = dim;
		arg[i].mod = mod;
		arg[i].inp = inp;
		arg[i].out = out;
		arg[i].los = 0.0;
		arg[i].grd = buf[i];
		arg[i].ii  = j;
		j += num[i];
		arg[i].jj  = j;
    	pthread_create( &pid[i], NULL, worker, (void*) &arg[i] );
	}
   	for( i = 0; i < NCPU; i++ ) pthread_join( pid[i], NULL );
	set_zeros( dim, grd );
	*los = 0.0;
   	for( i = 0; i < NCPU; i++ ) {
		for( j = 0; j < dim; j++ )
			for( k = 0; k < grd[j].row * grd[j].col; k++ )
				grd[j].mat[k] += buf[i][j].mat[k];
		*los += arg[i].los;
	}
	free( num ); free( pid ); free( arg );
	for( i = 0; i < NCPU; i++ ) { for( j = 0; j < dim; j++ ) free( buf[i][j].mat ); free( buf[i] ); } free( buf );
}


double dot_product( int dim, Tmat *ma, Tmat *mb ) {
    int		i, j;
    double	o = 0.0;
    for( i = 0; i < dim; i++ ) 
		for( j = 0; j < ma[i].row * ma[i].col; j++ )
			o += ma[i].mat[j] * mb[i].mat[j];
    return( o );
}

int fire( int dim, Tmat *mod, int nel, Tmat *inp, double *out ) {
	int		i, j, nstp, dstp, nfun, it;
	double	loss, norm, ssiz, alph, xlss, vsiz, last, tmp;
	Tmat	*grad, *velo, *step, *xcof;
	time_t	t0;

	grad = (Tmat*) malloc( dim * sizeof( Tmat ) );
	velo = (Tmat*) malloc( dim * sizeof( Tmat ) );
	step = (Tmat*) malloc( dim * sizeof( Tmat ) );
	xcof = (Tmat*) malloc( dim * sizeof( Tmat ) );
	for( i = 0; i < dim; i++ ) {
		grad[i].row = mod[i].row;
		grad[i].col = mod[i].col;
		grad[i].mat = (double*) malloc( mod[i].row * mod[i].col * sizeof( double ) );
		velo[i].row = mod[i].row;
		velo[i].col = mod[i].col;
		velo[i].mat = (double*) malloc( mod[i].row * mod[i].col * sizeof( double ) );
		step[i].row = mod[i].row;
		step[i].col = mod[i].col;
		step[i].mat = (double*) malloc( mod[i].row * mod[i].col * sizeof( double ) );
		xcof[i].row = mod[i].row;
		xcof[i].col = mod[i].col;
		xcof[i].mat = (double*) malloc( mod[i].row * mod[i].col * sizeof( double ) );
	}

	printf( "fire step_size:   %.6lf\n", STEP_SIZE );
	nstp = 0;
	ssiz = STEP_SIZE;
	alph = 0.1;
	dstp = 5;
	set_zeros( dim, velo );
	set_zeros( dim, step );
	t0   = time( NULL );
	prm_grad( dim, mod, nel, inp, out, &loss, grad );
	printf( ">> %ld seconds / gradient\n", time( NULL ) - t0 );
	norm = dot_product( dim, grad, grad );
	xlss = loss;
	last = loss;
	set_copy( dim, mod, xcof );
	printf( "          %30.1lf%30.1lf\n", loss, loss / nel );
	fflush( stdout );
	nfun = 0;
	it   = 0;
	while( it < STEP_NUMBER && nfun < MAX_TRIES ) {
		if( - dot_product( dim, velo, grad ) > 0.0 ) {
			vsiz = dot_product( dim, velo, velo );
			for( i = 0; i < dim; i++ )
				for( j = 0; j < mod[i].row * mod[i].col; j++ )
					velo[i].mat[j] = ( 1.0 - alph ) * velo[i].mat[j] - alph * grad[i].mat[j] / norm * vsiz;
			if( nstp > dstp ) { ssiz = min( ssiz * 1.1, STEP_SIZE ); alph *= 0.99; }
			nstp++;
		} else {
			alph = 0.1;
			ssiz *= 0.5;
			nstp = 0;
			set_zeros( dim, velo );
		}
		tmp = 0.0;
		for( i = 0; i < dim; i++ )
			for( j = 0; j < mod[i].row * mod[i].col; j++ ) {
				velo[i].mat[j] -= ssiz * grad[i].mat[j];
				step[i].mat[j] = ssiz * velo[i].mat[j];
				tmp += step[i].mat[j] * step[i].mat[j];
			}
		tmp = sqrt( tmp );
		if( tmp > ssiz )
			for( i = 0; i < dim; i++ )
				for( j = 0; j < mod[i].row * mod[i].col; j++ )
					step[i].mat[j] *= ssiz / tmp;
		for( i = 0; i < dim; i++ )
			for( j = 0; j < mod[i].row * mod[i].col; j++ )
				mod[i].mat[j] += step[i].mat[j];
		prm_grad( dim, mod, nel, inp, out, &loss, grad );
		norm = dot_product( dim, grad, grad );
		if( loss >= last ) nfun++;
		last = loss;
		printf( "%10d%30.1lf%30.1lf\n", it, loss, loss - xlss );
		if( loss < xlss ) {
			xlss = loss;
			nfun = 0;
			set_copy( dim, mod, xcof );
		}
		it++;
	}
	printf( "%10d%30.1lf%30.1lf\n", it, xlss, xlss / nel );
	set_copy( dim, xcof, mod );

	for( i = 0; i < dim; i++ ) {
		free( grad[i].mat ); free( velo[i].mat ); free( step[i].mat ); free( xcof[i].mat );
	}
	free( grad ); free( velo ); free( step ); free( xcof );
	return( it );
}


int  main( int argc, char **argv ) {
	FILE	*fd;
	char	buf[8];
	int		i, j, k, dim, nel;
	Tmat	*mod, *inp;
	double	*out;

	NCPU = atoi( argv[3] );
	STEP_SIZE = atof( argv[4] );
	if( STEP_SIZE < 0 ) STEP_SIZE = 0.01;
	printf( "------------------------------------------------------------------\n" );
	fd = fopen( argv[1], "rb" );
	fread( &buf, 4, 1, fd );
	memcpy( &dim, buf, 4 );
	printf( "number of layers: %d\n", (int)(dim/2) );
	mod = (Tmat*) malloc( dim * sizeof( Tmat ) );
	for( i = 0; i < dim; i++ ) {
		fread( &buf, 4, 1, fd );
		memcpy( &(mod[i].row), buf, 4 );
		fread( &buf, 4, 1, fd );
		memcpy( &(mod[i].col), buf, 4 );
		printf( "    %2d: %5d x %5d\n", i, mod[i].row, mod[i].col );
		mod[i].mat = (double*) malloc( mod[i].row * mod[i].col * sizeof( double ) );
		for( j = 0; j < mod[i].row * mod[i].col; j++ ) {
			fread( &buf, 8, 1, fd );
			memcpy( &(mod[i].mat[j]), buf, 8 );
		}
	}
	fclose( fd );
	printf( "------------------------------------------------------------------\n" );
	fd = fopen( argv[2], "rb" );
	fread( &buf, 4, 1, fd );
	memcpy( &k, buf, 4 );
	fread( &buf, 4, 1, fd );
	memcpy( &nel, buf, 4 );
	printf( "input dimension:  %d\n", k );
	printf( "number of inputs: %d\n", nel );
	printf( "number of cpus:   %d\n", NCPU );
	inp = (Tmat*) malloc( nel * sizeof( Tmat ) );
	for( i = 0; i < nel; i++ ) {
		inp[i].row = 1;
		inp[i].col = k;
		inp[i].mat = (double*) malloc( k * sizeof( double ) );
		for( j = 0; j < k; j++ ) {
			fread( &buf, 8, 1, fd );
			memcpy( &(inp[i].mat[j]), buf, 8 );
		}
	}
	out = (double*) malloc( nel * sizeof( double ) );
	for( i = 0; i < nel; i++ ) {
		fread( &buf, 8, 1, fd );
		memcpy( &(out[i]), buf, 8 );
	}
	fclose( fd );
	printf( "------------------------------------------------------------------\n" );

	for( k = 0; k < 10; k++ ) {
		if( fire( dim, mod, nel, inp, out ) < STEP_NUMBER ) STEP_SIZE *= 0.1;
		printf( "------------------------------------------------------------------\n" );
		fd = fopen( "last.raw", "wb" );
		memcpy( buf, &dim, 4 );
		fwrite( &buf, 4, 1, fd );
		for( i = 0; i < dim; i++ ) {
			memcpy( buf, &(mod[i].row), 4 );
			fwrite( &buf, 4, 1, fd );
			memcpy( buf, &(mod[i].col), 4 );
			fwrite( &buf, 4, 1, fd );
			for( j = 0; j < mod[i].row * mod[i].col; j++ ) {
				memcpy( buf, &(mod[i].mat[j]), 8 );
				fwrite( &buf, 8, 1, fd );
			}
		}
		fclose( fd );
	}

	for( i = 0; i < dim; i++ ) free( mod[i].mat ); free( mod );
	for( i = 0; i < nel; i++ ) free( inp[i].mat ); free( inp );
	free( out );
	return( 0 );
}

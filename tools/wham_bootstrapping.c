#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>


#define max(a,b) (((a)>(b))?(a):(b))


double dev_random( FILE *fd ) {
	char buf[8];
	uint64_t rnd;
	fread( buf, 1, 8, fd );
	memcpy( &rnd, buf, 8 );
	return( ( (double) rnd ) / 18446744073709551616.0 ); // 2^64
}

double gauss( FILE *fd, double mean, double stdv ) {
	return( mean + stdv * sqrt( -2.0 * log( dev_random( fd ) ) ) * cos( 2.0 * M_PI * dev_random( fd ) ) );
}

long randint( FILE *fd, long a, long b ) {
	return( a + (long)( dev_random( fd ) * (double)( b - a - 1.0 ) ) );
}


void wham( long try, long n_win, long n_bin, double kbt, double *_ss, double *_nn,
			double *_uu, double *_eu, double *_ff, double *_fo, double *_rr, double *pmf, double *rms ) {
	double t, r, f_tol = 1.e-3;
	long i, j, k, f, mxit = 10000;

	for( i = 0; i < n_bin; i++ ) { _rr[i] = 0.0; }
	f = 0;
	k = 0;
	while( k < mxit && ! f ) {
		for( j = 0; j < n_win; j++ ) _fo[j] = _ff[j];
		for( i = 0; i < n_bin; i++ ) {
			for( t = 0.0, j = 0; j < n_win; j++ ) { t += _ss[j] * exp( - ( _uu[i*n_win+j] - _ff[j] ) / kbt ); }
			_rr[i] = _nn[i] / t;
		}
		for( j = 0; j < n_win; j++ ) {
			for( t = 0.0, i = 0; i < n_bin; i++ ) { t += _eu[i*n_win+j] * _rr[i]; }
			_ff[j] = - kbt * log( t );
		}
		for( t = fabs( _fo[0] - _ff[0] ), j = 1; j < n_win; j++ ) { t = max( t, fabs( _fo[j] - _ff[j] ) ); }
		f = t < f_tol;
		k += 1;
	}
	printf( "#%10ld%4ld%10ld\n", try, f, k );
	for( t = 0.0, i = 0; i < n_bin; i++ ) { t += _rr[i]; }
	for( i = 0; i < n_bin; i++ ) {
		_rr[i] /= t;
		if( _rr[i] > 0.0 ) {
			r = - kbt * log( _rr[i] );
			pmf[i] += r;
			rms[i] += r * r;
		}
	}
}

int main( int argc, char **argv ){
	FILE *fd;
	double x_min, x_max, b_wid, f_tol, kbt;
	double *_fc = NULL, *_rf = NULL, *_gm = NULL, *_gr = NULL, *_nn = NULL, *_ss = NULL, *_ff = NULL;
   	double *_rr = NULL, *_uu = NULL, *_eu = NULL, *_fo = NULL, *crd = NULL, *pmf = NULL, *rms = NULL, r1, r2;
	long i, j, k, c, n_bin, n_bts, n_win, f_dec;

	n_win = argc - 7;
	x_min = atof( argv[1] );
	x_max = atof( argv[2] );
	n_bin = atol( argv[3] ); if( n_bin <= 0 ) { n_bin = 2 * n_win; }
	b_wid = ( x_max - x_min ) / ( (double) n_bin );
	kbt   = atof( argv[4] ) * 0.0083144621;
	n_bts = atol( argv[5] ); if( n_bts < 1 ) { n_bts = 1; }
	f_dec = atol( argv[6] );

	printf( "#n_win: %ld\n", n_win );
	printf( "#x_min: %.4lf\n", x_min );
	printf( "#x_max: %.4lf\n", x_max );
	printf( "#n_bin: %ld\n", n_bin );
	printf( "#d_bin: %lf\n", b_wid );
	printf( "#temp.: %.2lf\n", atof( argv[4] ) );
	printf( "#n_bts: %ld\n", n_bts );
	printf( "#f_dec: %ld\n", f_dec );

	_fc = (double*) malloc( n_win * sizeof( double ) );
	_rf = (double*) malloc( n_win * sizeof( double ) );
	_gm = (double*) malloc( n_win * sizeof( double ) );
	_gr = (double*) malloc( n_win * sizeof( double ) );
	_ss = (double*) malloc( n_win * sizeof( double ) );
	_nn = (double*) malloc( n_bin * sizeof( double ) );
	crd = (double*) malloc( n_bin * sizeof( double ) );
	pmf = (double*) malloc( n_bin * sizeof( double ) );
	rms = (double*) malloc( n_bin * sizeof( double ) );
	for( i = 0; i < n_bin; i++ ) {
		_nn[i] = 0.0;
		pmf[i] = 0.0; 
		rms[i] = 0.0;
		crd[i] = x_min + b_wid * ( (double) i + 0.5 );
	}

	for( i = 7; i < argc; i++ ) {
		j = i - 7;
		fd = fopen( argv[i], "rt" );
		fscanf( fd, "%lf%lf", &r1, &r2 );
		_fc[j] = r1;
		_rf[j] = r2;
		_ss[j] = 0.0;
		_gm[j] = 0.0;
		_gr[j] = 0.0;
		while( fscanf( fd, "%lf", &r1 ) == 1 ) {
			k = (long)( ( r1 - x_min ) / b_wid );
			if( k >= 0 && k < n_bin ) {
				_nn[k] += 1.0;
				_ss[j] += 1.0;
				_gm[j] += r1;
				_gr[j] += r1 * r1;
			}
		}
		_gm[j] /= _ss[j];
		_gr[j] = sqrt( fabs( _gr[j] / _ss[j] - _gm[j] * _gm[j] ) );
		fclose( fd );
	}

	_uu = (double*) malloc( n_bin * n_win * sizeof( double ) );
	_eu = (double*) malloc( n_bin * n_win * sizeof( double ) );
	for( i = 0; i < n_bin; i++ ) {
		for( j = 0; j < n_win; j++ ) {
			r1 = 0.5 * _fc[j] * ( crd[i] - _rf[j] ) * ( crd[i] - _rf[j] );
			_uu[i*n_win+j] = r1;
			_eu[i*n_win+j] = exp( - r1 / kbt );
		}
	}

	_ff = (double*) malloc( n_win * sizeof( double ) );
	_fo = (double*) malloc( n_win * sizeof( double ) );
	_rr = (double*) malloc( n_bin * sizeof( double ) );
	for( i = 0; i < n_win; i++ ) { _ff[i] = 0.0; }

	wham( 0, n_win, n_bin, kbt, _ss, _nn, _uu, _eu, _ff, _fo, _rr, pmf, rms );

	fd = fopen( "/dev/urandom", "rb" );

	if( f_dec == 1 ) {
		// window-based gaussian distribution
		for( c = 1; c < n_bts; c++ ) {
			for( i = 0; i < n_bin; i++ ) _nn[i] = 0.0;
			for( j = 0; j < n_win; j++ ) {
				i = 0;
				while( i < (long) _ss[j] ) {
					k  = (long)( ( gauss( fd, _gm[j], _gr[j] ) - x_min ) / b_wid );
					if( k >= 0 && k < n_bin ) {
						_nn[k] += 1.0;
						i += 1;
					}
				}
			}
			wham( c, n_win, n_bin, kbt, _ss, _nn, _uu, _eu, _ff, _fo, _rr, pmf, rms );
		}
	} else {
		// global-based MC density
		double popu, mxdn, *dens = NULL;
		dens = (double*) malloc( n_bin * sizeof( double ) );
		for( popu =  _nn[0], i = 1; i < n_bin; i++ ) { popu += _nn[i]; }
		for( mxdn = dens[0], i = 0; i < n_bin; i++ ) { dens[i] = _nn[i] / popu; mxdn = max( mxdn, dens[i] ); }
		for( c = 1; c < n_bts; c++ ) {
			for( i = 0; i < n_bin; i++ ) _nn[i] = 0.0;
			i = 0;
			while( i < (long) popu ) {
				k = randint( fd, 0, n_bin );
				if( dens[k] >= dev_random( fd ) * mxdn ) {
					_nn[k] += 1.0;
					i += 1;
				}
			}
			wham( c, n_win, n_bin, kbt, _ss, _nn, _uu, _eu, _ff, _fo, _rr, pmf, rms );
		}
		free( dens );
	}

	fclose( fd );

	printf( "#%19s%20s%20s\n", "Coordinate", "<PMF>", "Stdev" );
	for( i = 0; i < n_bin; i++ ) {
		pmf[i] /= (double) n_bts;
		rms[i] = sqrt( fabs( rms[i] / n_bts - pmf[i] * pmf[i] ) );
		printf( "%20.10lf%20.10lf%20.10lf\n", crd[i], pmf[i], rms[i] );
	}

	free( _fc ); free( _rf ); free( _gm ); free( _gr ); free( _ss ); free( _nn );
	free( _uu ); free( _eu ); free( _rr ); free( _fo ); free( _ff );
	free( crd ); free( pmf ); free( rms );
	return( 0 );
}

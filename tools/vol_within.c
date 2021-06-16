#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
	FILE	*fd, *hd;
	char	buf[256];
	double	t, dx, dy, dz, lx, ly, lz, vv;
	double	*xyzr;
	int		i, j, k, l, m, na, nx, ny, nz, ss, cr, cx, cy, cz, c2, tx, ty, tz;
	double	***grid;

	fd = fopen( "within.acs", "rb" );
	if( fd == NULL ) return( -1 );
	fread( buf, 1, 4, fd );
	memcpy( &ss, &buf[0], 4 );
	xyzr = (double*) malloc( ss * 4 * sizeof( double ) );
	for( i = 0; i < ss * 4; i++ ) {
		fread( buf, 1, 8, fd );
		memcpy( &xyzr[i], &buf[0], 8 );
		xyzr[i] /= 0.52917721092;
	}
	fclose( fd );

	// tune vdw radii (x3)
	for( i = 3; i < ss * 4; i+=4 ) xyzr[i] *= 3.0;

	fd = fopen( "volume.cube", "rt" );
	if( fd == NULL ) return( -1 );

	hd = fopen( "within.cube", "wt" );
	if( hd == NULL ) { fclose( fd ); return( -1 ); }

	// head (dgrid has a bug and there's only one header line...)
	fgets( buf, 256, fd );
	fputs( buf, hd );
	fgets( buf, 256, fd );
	fputs( buf, hd );

	// number of atoms, and lower corner
	fgets( buf, 256, fd );
	sscanf( buf, "%d%lf%lf%lf", &na, &lx, &ly, &lz );
	fprintf( hd, "%5d%12.6lf%12.6lf%12.6lf\n", ss, lx, ly, lz );

	// X: grid points and vector displacement
	fgets( buf, 256, fd );
	sscanf( buf, "%d%lf%lf%lf", &nx, &dx, &dy, &dz );
	fprintf( hd, "%5d%12.6lf%12.6lf%12.6lf\n", nx, dx, dy, dz );

	// Y: grid points and vector displacement
	fgets( buf, 256, fd );
	sscanf( buf, "%d%lf%lf%lf", &ny, &dx, &dy, &dz );
	fprintf( hd, "%5d%12.6lf%12.6lf%12.6lf\n", ny, dx, dy, dz );
	
	// Z: grid points and vector displacement
	fgets( buf, 256, fd );
	sscanf( buf, "%d%lf%lf%lf", &nz, &dx, &dy, &dz );
	fprintf( hd, "%5d%12.6lf%12.6lf%12.6lf\n", nz, dx, dy, dz );

	grid = malloc( nx * sizeof( double* ) );
	for( i = 0; i < nx; i++ ) {
		grid[i] = malloc( ny * sizeof( double* ) );
		for( j = 0 ; j < ny; j++ ) 
			grid[i][j] = malloc( nz * sizeof( double* ) );
	}

	// Coordinates
	for( i = 0; i < ss; i++ ) {
		j = i * 4;
		fprintf( hd, "%5d%12.6lf%12.6lf%12.6lf%12.6lf\n", 0, 0.0, xyzr[j], xyzr[j+1], xyzr[j+2] );
	}

	for( i = 0; i < na; i++ ) fgets( buf, 256, fd );
	
	// Cube formatting: 6 columns (max) of %13.5le
	vv = 0.0;
	for( i = 0; i < nx; i++ ) {
		for( j = 0; j < ny; j++ ) {
			k = 0;
			while( k < nz ) {
				fgets( buf, 256, fd );
				for( l = 0; l < (int)( ( strnlen( buf, 256 ) - 1 ) / 13 ); l++ ) {
					sscanf( &buf[13*l], "%lf", &t );
					if( t > 0.0 ) vv += 1.0;
					grid[i][j][k] = t;
					k++;
				}
			}
		}
	}

	printf( "VOL: %.0lf\n", vv );

	fclose( fd );

	vv = 0.0;
	for( l = 0; l < ss; l ++ ) {
		m  = l * 4;
		cr = (int)( xyzr[m+3] / dz );
		c2 = cr * cr;
		cx = (int)( ( xyzr[m]   - lx ) / dz );
		cy = (int)( ( xyzr[m+1] - ly ) / dz );
		cz = (int)( ( xyzr[m+2] - lz ) / dz );
    	for( i = -cr; i <= cr; i++ ) {
			tx = cx + i;
    		for( j = -cr; j <= cr; j++ ) {
				ty = cy + j;
    			for( k = -cr; k <= cr; k++ ) {
					tz = cz + k;
    				if( i*i + j*j + k*k <= c2 && tx >=0 && tx < nx && ty >= 0 && ty < ny && tz >= 0 && tz < nz && grid[tx][ty][tz] == 1.0 ) {
						grid[tx][ty][tz] = 2.0;
						vv += 1.0;
					}
				}
			}
		}
	}

	printf( "VOL: %.0lf\n", vv );

	for( i = 0; i < nx; i++ ) {
		for( j = 0; j < ny; j++ ) {
			for( k = 0; k < nz; k++ ) {
				fprintf( hd, "%13.5le", grid[i][j][k] );
				if( k%6 == 5 ) fprintf( hd, "\n" );
			}
			if( k%6 != 5 ) fprintf( hd, "\n" );
		}
	}

	fclose( hd );

	for( i = 0; i < nx; i++ ) {
		for( j = 0 ; j < ny; j++ ) free( grid[i][j] );
		free( grid[i] );
	}
	free( grid );

	return( 0 );
}

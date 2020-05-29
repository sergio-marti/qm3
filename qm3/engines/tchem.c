#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <mpi.h>
#include <signal.h>


/* compile mpich2-1.4.1p1 with the following flags:

./configure --prefix=@@@@ --enable-fc  --enable-shared  --with-namepublisher=file

*/


// use an integer multiple of 3*8 = 24 bytes
#define SCKTSIZE    24000


#define min(a,b) (((a)<(b))?(a):(b))


MPI_Comm    com = -1;


void sig_handler( int sig ) {
    int		fake = 0;

    MPI_Send( &fake, 1, MPI_INT, 0, 0, com );
    MPI_Comm_disconnect( &com );
    MPI_Finalize();
    exit( 0 );
}


int main( int argc, char **argv ) {
    int					sckt, chld, tmp, flg, i, j, k, r, a, n_qm, n_mm, dim = 2 * 32 * 128;
    char				prt[MPI_MAX_PORT_NAME], *o, *p, *buf, *smb, inp[dim+1];
    double				dip[4], ene, *qm_crd, *qm_chg, *mm_crd, *mm_chg, *grd;
    char				s_ene[8], *s_qm_crd, *s_mm_crd, *s_mm_chg;
    struct sockaddr_un	srvr;
    FILE				*fd;
    MPI_Status			sts;
    struct sigaction	sig;


    if( argc != 4 ) return( 2 );

    sig.sa_handler = sig_handler;
    sigemptyset( &sig.sa_mask );
    sig.sa_flags = SA_RESTART;
    if( sigaction( SIGINT, &sig, NULL ) == -1 ) { return( 3 ); }

    MPI_Init( NULL, NULL );

    // input:parse
    if( ( fd = fopen( argv[3], "rt" ) ) == NULL ) { 
    	MPI_Finalize();
    	return( 1 );
    }
    buf = (char*) malloc( 256 * sizeof( char ) );
    memset( &inp[0], 32, dim );
    i = 0;
    while( fgets( buf, 256, fd ) != NULL ) {
    	j = 0; o = &buf[0]; p = &buf[0]; while( ! isblank( *p ) ) { p++; j++; } *p = 0;
    	memcpy( &inp[i*128], o, j ); i++;
    	p++; while( isblank( *p ) ) p++;
    	j = 0; o = p; while( *p != 0 ) { p++; j++; } p--; *p = 0;
    	memcpy( &inp[i*128], o, j - 1 ); i++;
    }
/*
    memcpy( &inp[i*128], "run",      3 ); i++; memcpy( &inp[i*128], "gradient", 8 ); i++;
    memcpy( &inp[i*128], "amber",    5 ); i++; memcpy( &inp[i*128], "yes",      3 ); i++;
    memcpy( &inp[i*128], "dftgrid",  7 ); i++; memcpy( &inp[i*128], "1",        1 ); i++;
    memcpy( &inp[i*128], "threall",  7 ); i++; memcpy( &inp[i*128], "1.e-11",   6 ); i++;
    memcpy( &inp[i*128], "convthre", 8 ); i++; memcpy( &inp[i*128], "3.e-5",    5 ); i++;
*/
    memcpy( &inp[i*128], "end",     3 ); i++;
    dim = ( i + 1 ) * 128; inp[dim] = 0;
fprintf( stderr, "input_buff:\n%s\n-------------\n", inp );
    free( buf );
    fclose( fd );

    // mpi:connect
/*
    if( ( fd = fopen( argv[1], "rt" ) ) == NULL ) { 
    	close( sckt );
    	MPI_Finalize();
    	return( 1 );
    }
    buf = (char*) malloc( 256 * sizeof( char ) );
    while( fgets( buf, 256, fd ) != NULL ) {
    	k = strlen( buf );
    	i = 0; j = 0;
    	o = &buf[0]; 
    	while( *o != 0 && j == 0 ) {
    		if( *o == 'p' && i+10 < k && strncmp( o, "port_name:", 10 ) == 0 ) {
    			o += 10; while( *o == ' ' ) o++;
    			p = o; while( *p != 0 ) p++; p--; *p = 0;
    			j = 1;
    			snprintf( prt, MPI_MAX_PORT_NAME, "%s", o );
    		}
    		o++; i++;
    	}
    }
    free( buf );
    fclose( fd );
*/

    MPI_Lookup_name( argv[1], MPI_INFO_NULL, prt );
fprintf( stderr, "terachem_port: %s\n", prt );
    MPI_Comm_connect( prt, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &com );

    // unix:initialize
    sckt = socket( AF_UNIX, SOCK_STREAM, 0 );
    srvr.sun_family = AF_UNIX;
    snprintf( srvr.sun_path, 80, "%s", argv[2] );
    unlink( srvr.sun_path );

    if( bind( sckt, (struct sockaddr *) &srvr, sizeof( struct sockaddr_un ) ) == -1 ) {
    	MPI_Finalize();
    	return( 1 );
    }

    if( listen( sckt, 4 ) == -1 ) {
    	MPI_Finalize();
    	return( 1 );
    }
fprintf( stderr, "unix: listening!\n" );

    flg = 0;

    while( 1 ){

    	// unix:accept
    	if( ( chld = accept( sckt, NULL, NULL ) ) != -1 ) {
fprintf( stderr, "\n>> child accepted!\n" );
    
    		// unix:read n_qm + n_mm && smb
    		bzero( s_ene, 8 );
    		recv( chld, &s_ene[0], 8, 0 );
    		memcpy( &n_qm, &s_ene[0], 4 );
    		memcpy( &n_mm, &s_ene[4], 4 );
//fprintf( stderr, "n_qm = %d / n_mm = %d\n", n_qm, n_mm );

    		if( flg == 0 ) { smb = (char*) malloc( ( 2 * n_qm + 1 ) * sizeof( char ) ); }
    		bzero( smb, 2 * n_qm + 1 );
    		recv( chld, &smb[0], 2 * n_qm, 0 );
//fprintf( stderr, "smb = [%s]\n", smb );

    		if( flg == 0 ) {

    			// mpi:init
    			MPI_Send( inp, dim, MPI_CHAR, 0, 2, com );
    			MPI_Send( &n_qm, 1, MPI_INT, 0, 2, com );
    			MPI_Send( smb, 2 * n_qm, MPI_CHAR, 0, 2, com );

    			// local arrays
    			qm_crd = (double*) malloc( 3 * n_qm * sizeof( double ) );
    			qm_chg = (double*) malloc( n_qm * sizeof( double ) );
    			mm_crd = (double*) malloc( 3 * n_mm * sizeof( double ) );
    			mm_chg = (double*) malloc( n_mm * sizeof( double ) );
    			grd    = (double*) malloc( 3 * ( n_qm + n_mm ) * sizeof( double ) );
    
    			s_qm_crd = (char*) malloc( 8 * 3 * n_qm * sizeof( char ) );
    			s_mm_chg = (char*) malloc( 8 * n_mm * sizeof( char ) );
    			s_mm_crd = (char*) malloc( 8 * 3 * n_mm * sizeof( char ) );
    			buf      = (char*) malloc( SCKTSIZE * sizeof( char ) );
    		}
    
    		// unix:read MM charges
    		k = 0; r = 8 * n_mm;
    		if( r > 0 ) {
//fprintf( stderr, "reading MM-charges\n" );
    			bzero( s_mm_chg, r );
    			a = recv( chld, &s_mm_chg[0], min( r, SCKTSIZE ), 0 );
    			k += a; r -= a;
//fprintf( stderr, "\t%d bytes remaining...\n", r );
    			while( r > 0 ) { 
    				a = recv( chld, &s_mm_chg[k], min( r, SCKTSIZE ), 0 );
    				k += a; r -= a;
//fprintf( stderr, "\t%d bytes remaining...\n", r );
    			}
    			for( i = 0; i < n_mm; i++ ) memcpy( &mm_chg[i], &s_mm_chg[8*i], 8 );
    		}

    		// unix: read QM coordinates
//fprintf( stderr, "reading QM-coordinates\n" );
    		k = 0; r = 8 * 3 * n_qm;
    		bzero( s_qm_crd, r );
    		a = recv( chld, &s_qm_crd[0], min( r, SCKTSIZE ), 0 );
    		k += a; r -= a;
//fprintf( stderr, "\t%d bytes remaining...\n", r );
    		while( r > 0 ) {
    			a = recv( chld, &s_qm_crd[k], min( r, SCKTSIZE ), 0 );
    			k += a; r -= a;
//fprintf( stderr, "\t%d bytes remaining...\n", r );
    		}
    		for( i = 0; i < 3 * n_qm; i++ ) memcpy( &qm_crd[i], &s_qm_crd[8*i], 8 );

    		// unix: read MM coordinates
    		k = 0; r = 8 * 3 * n_mm;
    		if( r > 0 ) {
//fprintf( stderr, "reading MM-coordinates\n" );
    			bzero( s_mm_crd, r );
    			a = recv( chld, &s_mm_crd[0], min( r, SCKTSIZE ), 0 );
    			k += a; r -= a;
//fprintf( stderr, "\t%d bytes remaining...\n", r );
    			while( r > 0 ) {
    				a = recv( chld, &s_mm_crd[k], min( r, SCKTSIZE ), 0 );
    				k += a; r -= a;
//fprintf( stderr, "\t%d bytes remaining...\n", r );
    			}
    			for( i = 0; i < 3 * n_mm; i++ ) memcpy( &mm_crd[i], &s_mm_crd[8*i], 8 );
    		}

    		// mpi:write info
    		MPI_Send( &n_qm, 1, MPI_INT, 0, 2, com );
    		MPI_Send( smb, 2 * n_qm, MPI_CHAR, 0, 2, com );
    		MPI_Send( qm_crd, 3 * n_qm, MPI_DOUBLE, 0, 2, com );
    		MPI_Send( &n_mm, 1, MPI_INT, 0, 2, com );
    		MPI_Send( mm_chg, n_mm, MPI_DOUBLE, 0, 2, com );
    		MPI_Send( mm_crd, 3 * n_mm, MPI_DOUBLE, 0, 2, com );

    		// mpi:read info
    		MPI_Recv( &ene, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, com, &sts );
    		MPI_Recv( qm_chg, n_qm, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, com, &sts );
    		MPI_Recv( dip, 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, com, &sts );
    		MPI_Recv( dip, 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, com, &sts );
    		MPI_Recv( dip, 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, com, &sts );
    		MPI_Recv( grd, 3 * ( n_qm + n_mm ), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, com, &sts );

    		// unix:write energy
fprintf( stderr, "sendind energy: %.6lf\n", ene );
    		bzero( s_ene, 8 );
    		memcpy( &s_ene[0], &ene, 8 );
    		send( chld, s_ene, 8, 0 );

    		// unix:write QM charges
//fprintf( stderr, "sending QM-charges:\n" );
    		bzero( buf, SCKTSIZE );
    		j = 0;
    		for( i = 0; i < n_qm; i++ ) {
    			memcpy( &buf[j], &qm_chg[i], 8 ); j += 8;
    			if( j >= SCKTSIZE ) {
    				send( chld, buf, SCKTSIZE, 0 );
    				bzero( buf, SCKTSIZE );
//fprintf( stderr, "\t%d/%d bytes sent...\n", j, 8*n_qm );
    				j = 0;
    			}
    		}
    		if( j > 0 ) { 
    			send( chld, buf, j, 0 );
//fprintf( stderr, "\t%d/%d bytes sent...\n", j, 8*n_qm );
    		}

    		// unix:write QM+MM gradients
//fprintf( stderr, "sending QM&MM-gradients\n" );
    		bzero( buf, SCKTSIZE );
    		j = 0;
    		for( i = 0; i < ( n_qm + n_mm ); i++ ) {
    			memcpy( &buf[j],   &grd[3*i], 8 ); j += 8;
    			memcpy( &buf[j], &grd[3*i+1], 8 ); j += 8;
    			memcpy( &buf[j], &grd[3*i+2], 8 ); j += 8;
    			if( j >= SCKTSIZE ) { 
    				send( chld, buf, SCKTSIZE, 0 );
    				bzero( buf, SCKTSIZE );
//fprintf( stderr, "\t%d/%d bytes sent...\n", j, 24*(n_qm+n_mm) );
    				j = 0;
    			 }
    		}
    		if( j > 0 ) { 
    			send( chld, buf, j, 0 );
//fprintf( stderr, "\t%d/%d bytes sent...\n", j, 24*(n_qm+n_mm) );
    		}

//    		free( smb ); free( qm_crd ); free( qm_chg ); free( mm_crd ); free( mm_chg ); free( grd );
//    		free( s_qm_crd ); free( s_mm_crd ); free( s_mm_chg ); free( buf );

    		if( flg == 0 ) { flg = 1; }
    		close( chld );
fprintf( stderr, "done!\n" );
    	}

    }

    // never reached...
    close( sckt );

    return( 0 );
}

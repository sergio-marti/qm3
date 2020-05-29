#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


typedef struct ptn_node { long i, j, k; struct ptn_node *n; } ptn_lst;
typedef struct sph_node { long r; ptn_lst *p; struct sph_node *n; } sph_lst;
typedef struct { long _i0, _if, *npt; char *grd; double dsp, *bmn, *crd; sph_lst **hsh; } targ;


//pthread_mutex_t    lock;



void* __fill_atoms( void *args ) {
    targ		*a = (targ*) args;
    long		l, l3, i, j, k, n01;
    ptn_lst		*c;

    n01 = a->npt[0] * a->npt[1];
    for( l = a->_i0; l < a->_if; l++ ) {
    	l3 = l * 3;
    	i  = (long)( ( a->crd[l3]   - a->bmn[0] ) / a->dsp );
    	j  = (long)( ( a->crd[l3+1] - a->bmn[1] ) / a->dsp );
    	k  = (long)( ( a->crd[l3+2] - a->bmn[2] ) / a->dsp );
    	c  = a->hsh[l]->p;
    	while( c != NULL ) {
//pthread_mutex_lock( &lock );
    		a->grd[ (i + c->i) + (j + c->j) * a->npt[0] + (k + c->k) * n01 ] += 1;
//pthread_mutex_unlock( &lock );
    		c = c->n;
    	}
    }
    return( NULL );
}


static PyObject* __molecular_grid( PyObject *self, PyObject *args ){
    PyObject	*o_xyz, *o_rad, *o_flg;
    long		i, j, k, l, l3, cpu, siz, rad, rd2, *rng, nit, knd = 0;
    double		mad, tmp, cnt;
    double		cte = 4.1887902047863905;
    double		vol = 0.0, err = 0.0, exc = 0.0;
    double		bmax[3] = { -9999., -9999., -9999. };
    double		bmin[3] = { +9999., +9999., +9999. };;
    long		npt[3];
    char		*grd;
    double		*crd, dsp, d3;
    sph_lst		*sph, **hsh;
    sph_lst		*sph_cur;
    ptn_lst		*ptn_cur;
    pthread_t	*pid;
    targ		*arg;
    time_t		t0;

    o_flg = Py_False;
    if( PyArg_ParseTuple( args, "lOOd|O", &cpu, &o_rad, &o_xyz, &dsp, &o_flg ) ) {

    	t0 = time( NULL );
    	d3  = dsp * dsp * dsp;
    	siz = (long) PyList_Size( o_rad );
    	crd = (double*) malloc( 3 * siz * sizeof( double ) );
    	if( siz < cpu ) { cpu = 1; }
fprintf(stderr,"CPU: %ld\n",cpu);
fprintf(stderr,"DSP: %lf\n",dsp);

    	hsh = (sph_lst**) malloc( siz * sizeof( sph_lst* ) );
    	sph = (sph_lst*)  malloc( sizeof( sph_lst ) );
    	sph->r = -1;
    	sph->p = NULL;
    	sph->n = NULL;

    	mad = .0;
    	for( l = 0; l < siz; l++ ) {
    		tmp = PyFloat_AsDouble( PyList_GetItem( o_rad, l ) );
    		mad = max( mad, tmp );
    		rad = (long)( tmp / dsp );

    		sph_cur = sph;
    		while( sph_cur != NULL && sph_cur->r != rad ) { sph_cur = sph_cur->n; }
    		if( sph_cur == NULL ) {
    			for( sph_cur = sph; sph_cur->n != NULL; sph_cur = sph_cur->n );

    			sph_cur->n = (sph_lst*) malloc( sizeof( sph_lst ) );
    			sph_cur->n->r = rad;
    			sph_cur->n->n = NULL;
    			sph_cur->n->p = (ptn_lst*) malloc( sizeof( ptn_lst ) );
    			sph_cur->n->p->i = 9999;
    			sph_cur->n->p->j = 9999;
    			sph_cur->n->p->k = 9999;
    			sph_cur->n->p->n = NULL;

    			ptn_cur = sph_cur->n->p;
    			rd2 = rad * rad;
    			cnt = 0.0;
    			for( i = -rad; i <= rad; i++ ) {
    				for( j = -rad; j <= rad; j++ ) {
    					for( k = -rad; k <= rad; k++ ) {
    						if( i*i + j*j + k*k <= rd2 ) {
    							ptn_cur->n = (ptn_lst*) malloc( sizeof( ptn_lst ) );
    							ptn_cur->n->i = i;
    							ptn_cur->n->j = j;
    							ptn_cur->n->k = k;
    							ptn_cur->n->n = NULL;
    							ptn_cur = ptn_cur->n;	
    							cnt += 1.0;
    						}
    					}
    				}
    			}
    			ptn_cur = sph_cur->n->p;
    			sph_cur->n->p = sph_cur->n->p->n;
    			free( ptn_cur );
    			sph_cur = sph_cur->n;
    			err = max( err, fabs( cte * tmp * tmp * tmp - cnt * d3 ) );
    			knd++;
    		}
    		hsh[l] = sph_cur;
    		l3     = l * 3;
    		for( i = 0; i < 3; i++ ) {
    			crd[l3+i] = PyFloat_AsDouble( PyList_GetItem( o_xyz, l3+i ) );
    			bmin[i] = min( bmin[i], crd[l3+i] );
    			bmax[i] = max( bmax[i], crd[l3+i] );
    		}
    	}
    	mad *= 1.1;
    	for( i = 0; i < 3; i++ ) { 
    		bmin[i] -= mad; bmax[i] += mad;
    		npt[i] = (long)( ( bmax[i] - bmin[i] ) / dsp ) + 1;
    	}

fprintf(stderr,"MIN: %8.3lf%8.3lf%8.3lf\n",bmin[0],bmin[1],bmin[2]);
fprintf(stderr,"MAX: %8.3lf%8.3lf%8.3lf\n",bmax[0],bmax[1],bmax[2]);
fprintf(stderr,"NPT: %8ld%8ld%8ld\n",npt[0],npt[1],npt[2]);
fprintf(stderr,"KND: %8ld\n",knd);
fprintf(stderr,"ERR: %.3lf _A^3\n",err*siz);

    	grd = (char*) malloc( npt[0] * npt[1] * npt[2] * sizeof( char ) );
    	for( l = 0; l < npt[0] * npt[1] * npt[2]; l++ ) grd[l] = 0;

    	pid = (pthread_t*) malloc( cpu * sizeof( pthread_t ) );
    	arg = (targ*) malloc( cpu * sizeof( targ ) );
    	rng = (long*) malloc( ( cpu + 1 ) * sizeof( long ) );
    	for( l = 0; l <= cpu; l++ ) rng[l] = 0;

    	nit = (long)round( (float)siz / (float)cpu );
    	if( nit != 0 ) for( l = 0; l < cpu; l++ ) rng[l] = l * nit;
    	rng[cpu] = siz;
    	for( l = 0; l < cpu; l++ ) {
    		arg[l]._i0 = rng[l];
    		arg[l]._if = rng[l+1];
    		arg[l].dsp = dsp;
    		arg[l].npt = npt;
    		arg[l].bmn = bmin;
    		arg[l].crd = crd;
    		arg[l].grd = grd;
    		arg[l].hsh = hsh;
    		pthread_create( &pid[l], NULL, __fill_atoms, (void*) &arg[l] );
    	}
    	for( l = 0; l < cpu; l++ ) pthread_join( pid[l], NULL );

    if( o_flg == Py_True ) {
    	FILE* fd = fopen( "volume.cube", "wt" );
    	fprintf( fd, "QM3:\n-- Molecular Volume --\n" );
    	tmp = 0.52917721092;
    	fprintf( fd, "%5ld%12.6lf%12.6lf%12.6lf\n", siz, bmin[0] / tmp, bmin[1] / tmp, bmin[2] / tmp );
    	fprintf( fd, "%5ld%12.6lf%12.6lf%12.6lf\n", npt[0], dsp / tmp, 0.0, 0.0 );
    	fprintf( fd, "%5ld%12.6lf%12.6lf%12.6lf\n", npt[1], 0.0, dsp / tmp, 0.0 );
    	fprintf( fd, "%5ld%12.6lf%12.6lf%12.6lf\n", npt[2], 0.0, 0.0, dsp / tmp );
    	for( l = 0; l < siz; l++ )
    		fprintf( fd, "%5ld%12.6lf%12.6lf%12.6lf%12.6lf\n", (long)0, 0.0,
    			crd[3*l] / tmp, crd[3*l+1] / tmp, crd[3*l+2] / tmp );
    	nit = npt[0] * npt[1];
    	l = 0;
    	for( i = 0; i < npt[0]; i++ ) {
    		for( j = 0; j < npt[1]; j++ ) {
    			for( k = 0; k < npt[2]; k++ ) {
    				if( grd[i+j*npt[0]+k*nit] > 0 ) { vol += 1.0; }
//    				if( grd[i+j*npt[0]+k*nit] > 1 ) { exc += 1.0; }
    				fprintf( fd, "%13.5le", (float) grd[i + j * npt[0] + k * nit] );
    				l++;
    				if( l%6 == 0 ) { fprintf( fd, "\n" ); }
    			}
    		}
    	}
    	fclose( fd );
    } else {
    	for( l = 0; l < npt[0] * npt[1] * npt[2]; l++ ) {
    		if( grd[l] > 0 ) { vol += 1.0; }
//    		if( grd[l] > 1 ) { exc += 1.0; }
    	}
    }

//fprintf(stderr,"EXC: %lf _A^3\n",exc*d3);
fprintf(stderr,"VOL: %lf _A^3\n",vol*d3);

    	free( pid ); free( rng ); free( arg );
    	free( crd ); 
    	sph_cur = sph;
    	while( sph_cur != NULL ) {
    		ptn_cur = sph_cur->p;
    		while( ptn_cur != NULL ) {
    			ptn_cur = ptn_cur->n;
    			free( sph_cur->p );
    			sph_cur->p = ptn_cur;
    		}
    		sph_cur = sph_cur->n;
    		free( sph );
    		sph = sph_cur;
    	}
    	free( hsh );
    	free( grd );

fprintf( stderr, "TIM: %ld _sec\n", time( NULL ) - t0 );
    	return( Py_BuildValue( "d", vol*d3 ) );
    }
    Py_INCREF( Py_None ); return( Py_None );
}


/*
This version is not useful for "real" problems in which intersections take
    place among more than 2 spheres...
The result is that the multiple exclussion volumes are considered more than one
    time, thus providing smaller volumes (compared to the grid method...)

Follow (someday) these analytical formulas:

    Molecular physics, 1991, 72, 1313-1345
>>    Journal of Physical Chemistry, 1987, 91, 4121-4122
    Mathematical Notes, 2001, 69, 421–428
    Nano Studies, 2015, 11, 111-126
*/
static PyObject* __molecular_spheres( PyObject *self, PyObject *args ){
    PyObject	*o_xyz, *o_rad;
    long		i, j, i3, j3, siz;
    double		cte = 1.0471975511965976;
    double		vol = 0.0, exc = 0.0;
    double		*crd, *r, *r2, d, d2, hi, hj;
    time_t		t0;

    if( PyArg_ParseTuple( args, "OO", &o_rad, &o_xyz ) ) {

    	t0 = time( NULL );
    	siz = (long) PyList_Size( o_rad );
    	crd = (double*) malloc( 3 * siz * sizeof( double ) );
    	r   = (double*) malloc( siz * sizeof( double ) );
    	r2  = (double*) malloc( siz * sizeof( double ) );

    	for( i = 0; i < siz; i++ ) {
    		r[i]  = PyFloat_AsDouble( PyList_GetItem( o_rad, i ) );
    		r2[i] = r[i] * r[i];
    		vol  += r2[i] * r[i];
    		i3    = i * 3;
    		for( j = 0; j < 3; j++ ) {
    			crd[i3+j] = PyFloat_AsDouble( PyList_GetItem( o_xyz, i3+j ) );
    		}
    	}
    	for( i = 0; i < siz - 1; i++ ) {
    		i3 = i * 3;
    		for( j = i + 1; j < siz; j++ ) {
    			j3  = j * 3;
    			d2  = ( crd[i3]   - crd[j3]   ) * ( crd[i3]   - crd[j3]   );
    			d2 += ( crd[i3+1] - crd[j3+1] ) * ( crd[i3+1] - crd[j3+1] );
    			d2 += ( crd[i3+2] - crd[j3+2] ) * ( crd[i3+2] - crd[j3+2] );
    			if( d2 < r2[i] + r2[j] + 2.0 * r[i] * r[j] ) {
    				d    = sqrt( d2 );
    				hi   = r[i] * ( 1.0 - ( r2[i] + d2 - r2[j] ) / ( 2.0 * r[i] * d ) );
    				hj   = r[j] * ( 1.0 - ( r2[j] + d2 - r2[i] ) / ( 2.0 * r[j] * d ) );
    				exc += hi * hi * ( 3.0 * r[i] - hi ) + hj * hj * ( 3.0 * r[j] - hj );
    			}
    		}
    	}
    	vol *= 4.0 * cte;
    	exc *= cte;
fprintf(stderr,"SPH: %lf _A^3\n",vol);
fprintf(stderr,"EXC: %lf _A^3\n",exc);
    	vol -= exc;
fprintf(stderr,"VOL: %lf _A^3\n",vol);
    	free( crd ); free( r ); free( r2 );
fprintf( stderr, "TIM: %ld _sec\n", time( NULL ) - t0 );
    	return( Py_BuildValue( "d", vol ) );
    }
    Py_INCREF( Py_None ); return( Py_None );
}


static struct PyMethodDef methods [] = {
    { "molecular",      (PyCFunction)__molecular_grid, METH_VARARGS },
    { 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moddef = {
    PyModuleDef_HEAD_INIT,
    "_volume",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__volume( void ) {
    PyObject    *my_module;
    my_module = PyModule_Create( &moddef );
    return( my_module );
}
#else
void init_volume( void ) {
    PyObject    *my_module;
    my_module = Py_InitModule( "_volume", methods );
}
#endif

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


//pthread_mutex_t	lock;



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


static PyObject* __molecular( PyObject *self, PyObject *args ){
	PyObject	*o_xyz, *o_rad;
	long		i, j, k, l, l3, cpu, siz, rad, rd2, *rng, nit, knd = 0;
	double		mad, tmp, cnt, vol = 0.0, err = 0.0, cte = 4.1887902047863905;
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

	if( PyArg_ParseTuple( args, "lOOd", &cpu, &o_rad, &o_xyz, &dsp ) ) {

		t0 = time( NULL );
fprintf(stderr,"CPU: %ld\n",cpu);
fprintf(stderr,"DSP: %lf\n",dsp);
		d3  = dsp * dsp * dsp;
		siz = (long) PyList_Size( o_rad );
		crd = (double*) malloc( 3 * siz * sizeof( double ) );

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
		mad *= 0.75;
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

		for( l = 0; l < npt[0] * npt[1] * npt[2]; l++ ) if( grd[l] > 0 ) { vol += 1.0; }
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


static struct PyMethodDef methods [] = {
    { "molecular", (PyCFunction)__molecular, METH_VARARGS },
    { 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_volume",
	NULL,
	-1,
	methods,
	NULL,
	NULL,
	NULL,
	NULL
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

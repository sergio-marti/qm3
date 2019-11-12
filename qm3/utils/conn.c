#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

double r_cov[110] = { 0.31, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
        1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76,
        1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22,
        1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75,
        1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39,
        1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01,
        1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87,
        1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32,
        1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06,
        2.00, 1.96, 1.90, 1.87, 1.80, 1.69, 1.60, 1.60, 1.60, 1.60,
        1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60 };



void __ij( long w, long n, long *i, long *j ) {
	int		f = 0;
	double	ii, jj, nn, ww;

	*i = -1; *j = -1;
	nn = (double) n; ww = (double) w;
	if( ww < 0.0 || ww >= nn * ( nn - 1.0 ) / 2.0 ) { return; }
	jj = nn - 1.0;
	while( jj >= .0 && f == 0 ) {
		ii = ( 2.0 * nn - 3.0 - sqrt( 8.0 * jj - 8.0 * ww + 4.0 * nn * ( nn - 3.0 ) + 1.0 ) ) / 2.0;
		if( fmod( ii, 1.0 ) == 0.0 && ii < jj ) { f = 1; *i = (long) ii; *j = (long) jj; }
		else { jj -= 1.0; }
	}
}



typedef struct con_bnd_node { long i,j; struct con_bnd_node *n; } con_bnd;
typedef struct { long siz, _i0, _if; long *num; double *xyz; con_bnd *bnd; } con_arg;



void* __connectivity( void *args ) {
	con_arg		*arg = (con_arg*) args;
	long		w, i, j, i3, j3;
	double		dr, r2;
	con_bnd		*p;

	p = arg->bnd;
	for( w = arg->_i0; w < arg->_if; w++ ) {
		__ij( w, arg->siz, &i, &j );
		if( i == -1 || j == -1 ) { continue; }
		if( arg->num[i] == 1 && arg->num[j] == 1 ) { continue; }
		i3  = 3 * i;
		j3  = 3 * j;
		r2  = ( r_cov[arg->num[i]] + r_cov[arg->num[j]] + 0.1 ) * ( r_cov[arg->num[i]] + r_cov[arg->num[j]] + 0.1 );
		dr  = ( arg->xyz[i3] - arg->xyz[j3] ) * ( arg->xyz[i3] - arg->xyz[j3] ) +
				( arg->xyz[i3+1] - arg->xyz[j3+1] ) * ( arg->xyz[i3+1] - arg->xyz[j3+1] ) +
				( arg->xyz[i3+2] - arg->xyz[j3+2] ) * ( arg->xyz[i3+2] - arg->xyz[j3+2] );
		if( dr <= r2 ) {
			arg->bnd->i++;
			p->n    = (con_bnd*) malloc( sizeof( con_bnd ) );
			p->n->i = i;
			p->n->j = j;
			p->n->n = NULL;
			p       = p->n;
		}
	}
	return( NULL );
}

static PyObject* w_connectivity( PyObject *self, PyObject *args ) {
	PyObject	*out, *object, *anum, *coor;
	double		*xyz;
	long		*num, i, j, n3, n, cpu, *lst, dsp, nit;
	pthread_t	*pid;
	con_arg		*arg;
	con_bnd		*ptr;

	if( PyArg_ParseTuple( args, "lO", &cpu, &object ) ) {

		coor = PyObject_GetAttrString( object, "coor" );
		n3   = PyList_Size( coor );
		n    = n3 / 3;
		xyz  = (double*) malloc( n3 * sizeof( double ) );
		for( i = 0; i < n3; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( coor, i ) );
		Py_DECREF( coor );

		anum = PyObject_GetAttrString( object, "anum" );
		num  = (long*) malloc( n * sizeof( long ) );
		for( i = 0; i < n; i++ ) num[i] = PyLong_AsLong( PyList_GetItem( anum, i ) );
		Py_DECREF( anum );

		nit = n * ( n - 1 ) / 2;
		dsp = (long) ((float)nit / (float)cpu);
		lst = (long*) malloc( (cpu+1) * sizeof( long ) );
		if( dsp == 0 ) { for( i = 0; i < cpu+1; i++ ) lst[i] = 0;  for( i = 0; i < nit + 1; i++ ) lst[i] = i; }
		else { for( i = 0; i < cpu; i++ ) lst[i] = i * dsp; lst[cpu] = nit; }
		pid = (pthread_t*) malloc( cpu * sizeof( pthread_t ) );
		arg = (con_arg*) malloc( cpu * sizeof( con_arg ) );
		for( i = 0; i < cpu; i++ ) {
			arg[i].siz    = n;
			arg[i]._i0    = lst[i];
			arg[i]._if    = lst[i+1];
			arg[i].num    = num;
			arg[i].xyz    = xyz;
			arg[i].bnd    = (con_bnd*) malloc( sizeof( con_bnd ) ); 
			arg[i].bnd->i = 0;
			arg[i].bnd->j = 0;
			arg[i].bnd->n = NULL;
			pthread_create( &pid[i], NULL, __connectivity, (void*) &arg[i] );
		}
		for( i = 0; i < cpu; i++ ) pthread_join( pid[i], NULL );

		n = 0;
		for( i = 0; i < cpu; i++ ) { n += arg[i].bnd->i; }
		out = PyList_New( n );
		j = 0;
		for( i = 0; i < cpu; i++ ) {
			ptr = arg[i].bnd->n;
			while( ptr != NULL ) {
				PyList_SetItem( out, j++, Py_BuildValue( "[l,l]", ptr->i, ptr->j ) );
				ptr = ptr->n;
			}
		}

		free( num ); free( xyz ); free( lst ); free( pid );
		for( i = 0; i < cpu; i++ ) {
			ptr = arg[i].bnd;
			while( ptr != NULL ) {
				ptr = ptr->n;
				free( arg[i].bnd );
				arg[i].bnd = ptr;
			}
		}
		free( arg );
		return( out );
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}




static struct PyMethodDef methods [] = {
	{ "connectivity", (PyCFunction)w_connectivity, METH_VARARGS },
    { 0, 0, 0 }
};



#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_conn",
	NULL,
	-1,
	methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC PyInit__conn( void ) {
	PyObject    *my_module;
	my_module = PyModule_Create( &moddef );
	return( my_module );
}
#else
void init_conn( void ) {
	PyObject	*my_module;
	my_module = Py_InitModule( "_conn", methods );
}
#endif

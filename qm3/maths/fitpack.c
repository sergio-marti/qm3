#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


void   curv1_( long*, double*, double*, double*, double*, long*, double*, double*, double*, long* );
void   curvs_( long*, double*, double*, double*, long*, double*, double*, double*, double*, double*, double*, long* );
double curv2_( double*, long*, double*, double*, double*, double* );
double curvd_( double*, long*, double*, double*, double*, double* );


static PyObject* __init( PyObject *self, PyObject *args ){
    PyObject	*o_x, *o_y, *o_y2;
    double		s, *x, *y, *y2, rn, ep, *ys, *tt, d = 1.0;
    long		i, n, er;
    long		ks = 1, k1 = 3;

    if( PyArg_ParseTuple( args, "OOd", &o_x, &o_y, &s ) ) {
    	n = (long) PyList_Size( o_x );

    	rn = (double) n;
    	ep = sqrt( 2.0 / rn );

    	x  = (double*) malloc( n * sizeof( double ) );
    	y  = (double*) malloc( n * sizeof( double ) );
    	y2 = (double*) malloc( n * sizeof( double ) );
    	tt = (double*) malloc( 9 * n * sizeof( double ) );
    	ys = (double*) malloc( n * sizeof( double ) );

    	for( i = 0; i < n; i++ ) {
    		x[i] = PyFloat_AsDouble( PyList_GetItem( o_x, i ) );
    		y[i] = PyFloat_AsDouble( PyList_GetItem( o_y, i ) );
    	}

//    	curv1_( &n, x, y, &d, &d, &k1, y2, tt, &s, &er );
    	curvs_( &n, x, y, &d, &ks, &rn, &ep, ys, y2, &s, tt, &er );

    	free( x ); free( y ); free( tt );
    	free( ys );

    	if( (int) er == 0 ) {

    		o_y2 = PyList_New( n );
    		for( i = 0; i < n; i++ ) PyList_SetItem( o_y2, i, PyFloat_FromDouble( y2[i] ) );
    		free( y2 );
    		return( o_y2 );

    	} else {

    		free( y2 );
    		Py_INCREF( Py_None ); return( Py_None );

    	}
    }
    Py_INCREF( Py_None ); return( Py_None );
}


static PyObject* __calc( PyObject *self, PyObject *args ){
    PyObject	*o_x, *o_y, *o_y2;
    double		s, *x, *y, *y2, rx, ry, rd;
    long		i, n;

    if( PyArg_ParseTuple( args, "dOOOd", &rx, &o_x, &o_y, &o_y2, &s ) ) {
    	n = (long) PyList_Size( o_x );

    	x  = (double*) malloc( n * sizeof( double ) );
    	y  = (double*) malloc( n * sizeof( double ) );
    	y2 = (double*) malloc( n * sizeof( double ) );

    	for( i = 0; i < n; i++ ) {
    		x[i]  = PyFloat_AsDouble( PyList_GetItem( o_x, i ) );
    		y[i]  = PyFloat_AsDouble( PyList_GetItem( o_y, i ) );
    		y2[i] = PyFloat_AsDouble( PyList_GetItem( o_y2, i ) );
    	}

    	ry = curv2_( &rx, &n, x, y, y2, &s );
    	rd = curvd_( &rx, &n, x, y, y2, &s );

    	free( x ); free( y ); free( y2 );

    	return( Py_BuildValue( "(d,d)", ry, rd ) );
    }
    Py_INCREF( Py_None ); return( Py_None );
}


static struct PyMethodDef methods [] = {
    { "init", (PyCFunction)__init, METH_VARARGS },
    { "calc", (PyCFunction)__calc, METH_VARARGS },
    { 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moddef = {
    PyModuleDef_HEAD_INIT,
    "_fitpack",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__fitpack( void ) {
    PyObject    *my_module;
    my_module = PyModule_Create( &moddef );
    return( my_module );
}
#else
void init_fitpack( void ) {
    PyObject    *my_module;
    my_module = Py_InitModule( "_fitpack", methods );
}
#endif

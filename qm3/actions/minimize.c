#include<Python.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdio.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


void cgp_cgfam_( long*, double*, double*, double*, double*, double*, double*, double*, long*, long*, long* );


static PyObject* w_cgp( PyObject *self, PyObject *args ) {
    long		size, step_number, print_frequency;
    double		gradient_tolerance;
    char		message[256];
    long		iflag, irest, imeth, i, k;
    double		grms, func;
    double		*x, *g, *w, *o, *d;
    PyObject 	*object, *fun_log, *obj_coor, *obj_grad, *obj_func;

    if( PyArg_ParseTuple( args, "OldlllO", 
    	&object, 
    	&step_number, 
    	&gradient_tolerance,
    	&print_frequency,
    	&irest,
    	&imeth,
    	&fun_log ) ) {

//    	size  = PyInt_AsLong( PyObject_GetAttrString( object, "size" ) );
    	size  = PyLong_AsLong( PyObject_GetAttrString( object, "size" ) );
    	x     = (double*) malloc( sizeof(double) * size );
    	g     = (double*) malloc( sizeof(double) * size );
    	d     = (double*) malloc( sizeof(double) * size );
    	o     = (double*) malloc( sizeof(double) * size );
    	w     = (double*) malloc( sizeof(double) * size );

    	PyObject_CallMethod( object, "get_grad", NULL );
    	obj_func = PyObject_GetAttrString( object, "func" );
    	obj_coor = PyObject_GetAttrString( object, "coor" );
    	obj_grad = PyObject_GetAttrString( object, "grad" );
    	func     = PyFloat_AsDouble( obj_func );
    	grms     = .0;
    	for( i = 0; i < size; i++ ) {
    		x[i]   = PyFloat_AsDouble( PyList_GetItem( obj_coor, i ) );
    		g[i]   = PyFloat_AsDouble( PyList_GetItem( obj_grad, i ) );
    		grms  += g[i] * g[i];
    	}
    	Py_DECREF( obj_func );
    	Py_DECREF( obj_coor );
    	Py_DECREF( obj_grad );
    	grms  = sqrt( grms / (double) size );
    	k     = 0;
    	iflag = 0;

    	snprintf( message, 256, "%30.5lf%20.10lf", func, grms );
    	PyObject_CallFunction( fun_log, "s", message );

    	while( k < step_number && grms > gradient_tolerance ) {

    		cgp_cgfam_( &size, x, &func, g, d, o, &gradient_tolerance, w, &iflag, &irest, &imeth );

    		if( iflag == -3 ) {

    			snprintf( message, 256, "\n -- Improper input parameters..." );
    			PyObject_CallFunction( fun_log, "s", message );
    			k = step_number + 1;

    		} else if( iflag == -2 ) {

    			snprintf( message, 256, "\n -- Descent was not obtained..." );
    			PyObject_CallFunction( fun_log, "s", message );
    			k = step_number + 1;

    		} else if( iflag == -1 ) {

    			snprintf( message, 256, "\n -- Line Search failure..." );
    			PyObject_CallFunction( fun_log, "s", message );
    			k = step_number + 1;

    		} else {

    			while( iflag == 2 )
    				cgp_cgfam_( &size, x, &func, g, d, o, &gradient_tolerance, w, &iflag, &irest, &imeth );

    			obj_coor = PyObject_GetAttrString( object, "coor" );
    			for( i = 0; i < size; i++ )
    				PyList_SetItem( obj_coor, i, PyFloat_FromDouble( x[i] ) );
    			Py_DECREF( obj_coor );
    			PyObject_CallMethod( object, "get_grad", NULL );
    			obj_func = PyObject_GetAttrString( object, "func" );
    			func     = PyFloat_AsDouble( obj_func );
    			obj_grad = PyObject_GetAttrString( object, "grad" );
    			grms     = .0;
    			for( i = 0; i < size; i++ ) {
    				g[i] = PyFloat_AsDouble( PyList_GetItem( obj_grad, i ) );
    				grms += g[i] * g[i];
    			}
    			grms = sqrt( grms / (double) size );
    			Py_DECREF( obj_func );
    			Py_DECREF( obj_grad );

    		}

    		k += 1;
    		
    		if( k%print_frequency == 0 && k <= step_number ) {
    			snprintf( message, 256, "%10ld%20.5lf%20.10lf", k, func, grms );
    			PyObject_CallFunction( fun_log, "s", message );
    		}

    		PyObject_CallMethod( object, "current_step", "(l)", k );

    	}

    	if( k%print_frequency != 0  && k <= step_number ) {
    		snprintf( message, 256, "%10ld%20.5lf%20.10lf", k, func, grms );
    		PyObject_CallFunction( fun_log, "s", message );
    	}

    	free( x ); free( g ); free( d ); free( o ); free( w );
    }

    Py_INCREF( Py_None );
    return( Py_None );
}


static struct PyMethodDef methods [] = {
    { "cgp",    (PyCFunction)w_cgp,   METH_VARARGS },
    { 0, 0, 0 }
};


static struct PyModuleDef moddef = {
    PyModuleDef_HEAD_INIT,
    "_minimize",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__minimize( void ) {
    PyObject    *my_module;
    my_module = PyModule_Create( &moddef );
    return( my_module );
}

#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


double pythag( double dx, double dy ) {
	double	x = fabs( dx );
	double	y = fabs( dy );
	if( x > y ) { return( x * sqrt( 1.0 + y * y / ( x * x ) ) ); }
	if( y == 0.0 ) { return( 0.0 ); }
	return( y * sqrt( 1.0 + x * x / ( y * y ) ) );
}


static PyObject* __regular( PyObject *self, PyObject *args ){
	PyObject	*o_x, *o_y, *o_dat, *o_g, *out;
	double		cx, *y, *ox, *oy, *oz, g[2], rz, rw, dd, ww;
	long		i, j, k, l, nx, ny, on;

	if( PyArg_ParseTuple( args, "OOOO", &o_x, &o_y, &o_dat, &o_g ) ) {

		nx = (long) PyList_Size( o_x );
		ny = (long) PyList_Size( o_y );
//fprintf(stderr,"%ld %ld\n",nx,ny);
		y  = (double*) malloc( ny * sizeof( double ) );
		on = (long) PyList_Size( o_dat );
//fprintf(stderr,"%ld\n",on);
		ox = (double*) malloc( on * sizeof( double ) );
		oy = (double*) malloc( on * sizeof( double ) );
		oz = (double*) malloc( on * sizeof( double ) );
		for( i = 0; i < ny; i++ ) { y[i] = PyFloat_AsDouble( PyList_GetItem( o_y, i ) ); }
		for( i = 0; i < on; i++ ) {
			out   = PyList_GetItem( o_dat, i );
			ox[i] = PyFloat_AsDouble( PyList_GetItem( out, 0 ) );
			oy[i] = PyFloat_AsDouble( PyList_GetItem( out, 1 ) );
			oz[i] = PyFloat_AsDouble( PyList_GetItem( out, 2 ) );
//if( i == 0 | i == on -1 ) fprintf(stderr,"%ld %lf %lf %lf\n",i,ox[i],oy[i],oz[i]);
		}
		g[0] = PyFloat_AsDouble( PyTuple_GetItem( o_g, 0 ) );
		g[1] = PyFloat_AsDouble( PyTuple_GetItem( o_g, 1 ) );
//fprintf(stderr,"%lf %lf\n",g[0],g[1]);
		out  = PyList_New( nx * ny );
		l    = 0;
		for( i = 0; i < nx; i++ ) {
			cx = PyFloat_AsDouble( PyList_GetItem( o_x, i ) );
			for( j = 0; j < ny; j++ ) {
				rz = 0.0;
				rw = 0.0;
				for( k = 0; k < on; k++ ) {
					dd = pythag( ( ox[k] - cx ) / g[0], ( oy[k] - y[j] ) / g[1] );
					ww = exp( - dd * dd );
					rz += oz[k] * ww;
					rw += ww;
				}
				PyList_SetItem( out, l++, PyFloat_FromDouble( rz / rw ) );
			}
		}
		free( y ); free( ox ); free( oy ); free( oz );
		return( out );
	}
	Py_INCREF( Py_None ); return( Py_None );
}


static struct PyMethodDef methods [] = {
    { "regular",      (PyCFunction)__regular, METH_VARARGS },
    { 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_grids",
	NULL,
	-1,
	methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC PyInit__grids( void ) {
	PyObject    *my_module;
	my_module = PyModule_Create( &moddef );
	return( my_module );
}
#else
void init_grids( void ) {
	PyObject    *my_module;
	my_module = Py_InitModule( "_grids", methods );
}
#endif

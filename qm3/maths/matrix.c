#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "matrix_core.c"



static PyObject* w_det( PyObject *self, PyObject *args ) {
	PyObject	*lst;
	double		*mat, out;
	long		i, row;

	if( PyArg_ParseTuple( args, "Ol", &lst, &row ) ) {
		if( PyList_Check( lst ) && PyList_Size( lst ) == row*row ) {
			mat = (double*)malloc(row*row*sizeof(double));
			for( i=0; i<row*row; i++ ) mat[i] = PyFloat_AsDouble( PyList_GetItem( lst, i ) );
			out = __det( mat, row );
			free( mat );
			return( Py_BuildValue( "d", out ) );
		} else { Py_INCREF( Py_None ); return( Py_None ); }
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}



static PyObject* w_matmul( PyObject *self, PyObject *args ) {
	PyObject	*a_lst, *b_lst, *out;
	double		*a_mat, *b_mat, *c_mat;
	long		a_row, a_col, b_row, b_col, i;

	if( PyArg_ParseTuple( args, "OllOll", &a_lst, &a_row, &a_col, &b_lst, &b_row, &b_col ) ) {
		if( PyList_Check( a_lst ) && PyList_Check( b_lst ) && a_col == b_row ) {
			a_mat = (double*)malloc((a_row*a_col)*sizeof(double));
			b_mat = (double*)malloc((b_row*b_col)*sizeof(double));
			c_mat = (double*)malloc((a_row*b_col)*sizeof(double));
			for( i=0; i<a_row*a_col; i++ ) a_mat[i] = PyFloat_AsDouble( PyList_GetItem( a_lst, i ) );
			for( i=0; i<b_row*b_col; i++ ) b_mat[i] = PyFloat_AsDouble( PyList_GetItem( b_lst, i ) );
			__mult( a_mat, a_row, a_col, b_mat, b_row, b_col, c_mat );
			out = PyList_New( a_row*b_col );
			for( i=0; i<a_row*b_col; i++ )
				PyList_SetItem( out, i, PyFloat_FromDouble( c_mat[i] ) );
			free( a_mat ); free( b_mat ); free( c_mat );
			return( Py_BuildValue( "O", out ) );
		} else { Py_INCREF( Py_None ); return( Py_None ); }
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}



static PyObject* w_invert( PyObject *self, PyObject *args ) {
	PyObject	*lst, *out;
	double		*mat;
	long		i, row, col; 

	if( PyArg_ParseTuple( args, "Oll", &lst, &row, &col ) ) {
		if( PyList_Check( lst ) && PyList_Size( lst ) == row*col ) {
			mat = (double*)malloc((row*col)*sizeof(double));
			for( i=0; i<row*col; i++ )
				mat[i] = PyFloat_AsDouble( PyList_GetItem( lst, i ) );
			if( row == col ) { __invert( mat, row ); } 
				else { __pseudoinvert( mat, row, col ); }
			out = PyList_New( row*col );
			for( i=0; i<row*col; i++ )
				PyList_SetItem( out, i, PyFloat_FromDouble( mat[i] ) );
			free( mat );
			return( Py_BuildValue( "O", out ) );
		} else { Py_INCREF( Py_None ); return( Py_None ); }
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}



static PyObject* w_diag( PyObject *self, PyObject *args ) {
	PyObject	*lst, *ovec, *oval;
	double		*mat, *vec, *val;
	long		i, row; 

	if( PyArg_ParseTuple( args, "Ol", &lst, &row ) ) {
		if( PyList_Check( lst ) && PyList_Size( lst ) == row*row ) {
			mat = (double*)malloc((row*row)*sizeof(double));
			vec = (double*)malloc((row*row)*sizeof(double));
			val = (double*)malloc((row)*sizeof(double));
			for( i=0; i<row*row; i++ )
				mat[i] = PyFloat_AsDouble( PyList_GetItem( lst, i ) );
			__diag( mat, vec, val, row );
			oval = PyList_New( row );
			for( i=0; i<row; i++ ) PyList_SetItem( oval, i, PyFloat_FromDouble( val[i] ) );
			ovec = PyList_New( row*row );
			for( i=0; i<row*row; i++ ) PyList_SetItem( ovec, i, PyFloat_FromDouble( vec[i] ) );
			free( mat ); free( vec ); free( val );
			return( Py_BuildValue( "(O,O)", oval, ovec ) );
		} else { Py_INCREF( Py_None ); return( Py_None ); }
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}



static PyObject* w_jacobi( PyObject *self, PyObject *args ) {
	PyObject	*lst, *ovec, *oval;
	double		*mat, *vec;
	long		i, row, f; 

	if( PyArg_ParseTuple( args, "Ol", &lst, &row ) ) {
		if( PyList_Check( lst ) && PyList_Size( lst ) == row*row ) {
			mat = (double*)malloc((row*row)*sizeof(double));
			vec = (double*)malloc((row*row)*sizeof(double));
			for( i=0; i<row*row; i++ )
				mat[i] = PyFloat_AsDouble( PyList_GetItem( lst, i ) );
			f = __fjacobi( mat, vec, row );
			oval = PyList_New( row );
			for( i=0; i<row; i++ ) PyList_SetItem( oval, i, PyFloat_FromDouble( mat[i*row+i] ) );
			ovec = PyList_New( row*row );
			for( i=0; i<row*row; i++ ) PyList_SetItem( ovec, i, PyFloat_FromDouble( vec[i] ) );
			free( mat ); free( vec );
			if( f ) return( Py_BuildValue( "(O,O,O)", oval, ovec, Py_True ) );
				else return( Py_BuildValue( "(O,O,O)", oval, ovec, Py_False ) );
		} else { Py_INCREF( Py_None ); return( Py_None ); }
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}



static struct PyMethodDef methods [] = {
	{ "det",     (PyCFunction)w_det,     METH_VARARGS },
	{ "matmul",  (PyCFunction)w_matmul,  METH_VARARGS },
	{ "invert",  (PyCFunction)w_invert,  METH_VARARGS },
	{ "diag",    (PyCFunction)w_diag,    METH_VARARGS },
	{ "jacobi",  (PyCFunction)w_jacobi,  METH_VARARGS },
    { 0, 0, 0 }
};



#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_matrix",
	NULL,
	-1,
	methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC PyInit__matrix( void ) {
	PyObject    *my_module, *d, *t;
	my_module = PyModule_Create( &moddef );
	d = PyModule_GetDict( my_module );
	t = PyUnicode_FromString( MCORE );
	PyDict_SetItemString( d, "version", t );
	Py_DECREF( t );
	return( my_module );
}
#else
void init_matrix( void ) {
	PyObject	*my_module, *d, *t;
	my_module = Py_InitModule( "_matrix", methods );
	d = PyModule_GetDict( my_module );
	t = PyString_FromString( MCORE );
	PyDict_SetItemString( d, "version", t );
	Py_DECREF( t );
}
#endif

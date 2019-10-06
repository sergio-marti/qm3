#include <Python.h>
#include <math.h>
#include <stdio.h>
#include <string.h>


double __a0 = 0.52917721067;
double __fe = 0.00159360164858473; // kcal/mol >> Ha


static PyObject* get_func( PyObject *self, PyObject *args ) {
	PyObject	*oqm = NULL, *omm = NULL;
	long		i, j, i3, j3, nqm, nmm, siz, w_i, w_j;
	double		*xqm = NULL, *xmm = NULL, w_e, w_r, w_f;
	double		dr[3], r2, ss, out = 0.0;
	char		buf[40];
	FILE		*fd;

	if( PyArg_ParseTuple( args, "OO", &oqm, &omm ) ) {

		nqm = PyList_Size( oqm );
		xqm = (double*) malloc( nqm * sizeof( double ) );
		for( i = 0; i < nqm; i++ ) xqm[i] = PyFloat_AsDouble( PyList_GetItem( oqm, i ) );

		nmm = PyList_Size( omm );
		xmm = (double*) malloc( nmm * sizeof( double ) );
		for( i = 0; i < nmm; i++ ) xmm[i] = PyFloat_AsDouble( PyList_GetItem( omm, i ) );

		fd = fopen( "qmmm.nbnd", "rb" );
		bzero( buf, 40 );
		fread( buf, 8, 1, fd );
		memcpy( &siz, &buf, 8 );
		for( i = 0; i < siz; i++ ) {
			bzero( buf, 40 );
			fread( buf, 40, 1, fd );
			memcpy( &w_i,  &buf[0], 8 ); memcpy( &w_j,  &buf[8], 8 );
			memcpy( &w_e, &buf[16], 8 ); memcpy( &w_r, &buf[24], 8 );
			memcpy( &w_f, &buf[32], 8 );
			i3 = w_i * 3;
			j3 = w_j * 3;
			r2 = 0.0;
			for( j = 0; j < 3; j++ ) { dr[j] = xqm[i3+j] - xmm[j3+j]; r2 += dr[j] * dr[j]; }
			ss = w_r / sqrt( r2 );
			ss = ss * ss * ss * ss * ss * ss;
			out += w_e * ss * ( ss - 2.0 ) * w_f;
		}
		fclose( fd );

		free( xqm ); free( xmm );

	}
	return( Py_BuildValue( "d", out * __fe ) );
}


static PyObject* get_grad( PyObject *self, PyObject *args ) {
	PyObject	*oqm = NULL, *omm = NULL, *ogg = NULL;
	long		i, j, i3, j3, nqm, nmm, siz, w_i, w_j;
	double		*xqm = NULL, *xmm = NULL, *xgg = NULL, w_e, w_r, w_f;
	double		dr[3], r2, ss, out = 0.0, df;
	char		buf[40];
	FILE		*fd;
	double 		__fg = __fe * __a0;

	if( PyArg_ParseTuple( args, "OOO", &oqm, &omm, &ogg ) ) {

		nqm = PyList_Size( oqm );
		xqm = (double*) malloc( nqm * sizeof( double ) );
		for( i = 0; i < nqm; i++ ) xqm[i] = PyFloat_AsDouble( PyList_GetItem( oqm, i ) );
		xgg = (double*) malloc( nqm * sizeof( double ) );
		for( i = 0; i < nqm; i++ ) xgg[i] = PyFloat_AsDouble( PyList_GetItem( ogg, i ) );

		nmm = PyList_Size( omm );
		xmm = (double*) malloc( nmm * sizeof( double ) );
		for( i = 0; i < nmm; i++ ) xmm[i] = PyFloat_AsDouble( PyList_GetItem( omm, i ) );

		fd = fopen( "qmmm.nbnd", "rb" );
		bzero( buf, 40 );
		fread( buf, 8, 1, fd );
		memcpy( &siz, &buf, 8 );
		for( i = 0; i < siz; i++ ) {
			bzero( buf, 40 );
			fread( buf, 40, 1, fd );
			memcpy( &w_i,  &buf[0], 8 ); memcpy( &w_j,  &buf[8], 8 );
			memcpy( &w_e, &buf[16], 8 ); memcpy( &w_r, &buf[24], 8 );
			memcpy( &w_f, &buf[32], 8 );
			i3 = w_i * 3;
			j3 = w_j * 3;
			r2 = 0.0;
			for( j = 0; j < 3; j++ ) { dr[j] = xqm[i3+j] - xmm[j3+j]; r2 += dr[j] * dr[j]; }
			ss = w_r / sqrt( r2 );
			ss = ss * ss * ss * ss * ss * ss;
			out += w_e * ss * ( ss - 2.0 ) * w_f;

			df = 12.0 * w_e * ss * ( 1.0 - ss ) / r2 * w_f * __fg;
			for( j = 0; j < 3; j++ ) xgg[i3+j] += df * dr[j];
		}
		fclose( fd );

		for( i = 0; i < nqm; i++ ) PyList_SetItem( ogg, i, PyFloat_FromDouble( xgg[i] ) );

		free( xqm ); free( xmm ); free( xgg );

	}
	return( Py_BuildValue( "d", out * __fe ) );
}


static PyObject* get_hess( PyObject *self, PyObject *args ) {
	PyObject	*oqm = NULL, *omm = NULL, *ogg = NULL, *ohh = NULL;
	long		i, j, i3, j3, nqm, nmm, siz, w_i, w_j, nh, ii;
	double		*xqm = NULL, *xmm = NULL, *xgg = NULL, *xhh = NULL, w_e, w_r, w_f;
	double		dr[3], r2, ss, out = 0.0, df, tt, r4, d2f;
	char		buf[40];
	FILE		*fd;
	double 		__fg = __fe * __a0;
	double		__fh = __fg * __a0;

	if( PyArg_ParseTuple( args, "OOOO", &oqm, &omm, &ogg, &ohh ) ) {

		nqm = PyList_Size( oqm );
		xqm = (double*) malloc( nqm * sizeof( double ) );
		for( i = 0; i < nqm; i++ ) xqm[i] = PyFloat_AsDouble( PyList_GetItem( oqm, i ) );
		xgg = (double*) malloc( nqm * sizeof( double ) );
		for( i = 0; i < nqm; i++ ) xgg[i] = PyFloat_AsDouble( PyList_GetItem( ogg, i ) );
		nh  = PyList_Size( ohh );
		xhh = (double*) malloc( nh * sizeof( double ) );
		for( i = 0; i < nh; i++ ) xhh[i] = PyFloat_AsDouble( PyList_GetItem( ohh, i ) );

		nmm = PyList_Size( omm );
		xmm = (double*) malloc( nmm * sizeof( double ) );
		for( i = 0; i < nmm; i++ ) xmm[i] = PyFloat_AsDouble( PyList_GetItem( omm, i ) );

		fd = fopen( "qmmm.nbnd", "rb" );
		bzero( buf, 40 );
		fread( buf, 8, 1, fd );
		memcpy( &siz, &buf, 8 );
		for( i = 0; i < siz; i++ ) {
			bzero( buf, 40 );
			fread( buf, 40, 1, fd );
			memcpy( &w_i,  &buf[0], 8 ); memcpy( &w_j,  &buf[8], 8 );
			memcpy( &w_e, &buf[16], 8 ); memcpy( &w_r, &buf[24], 8 );
			memcpy( &w_f, &buf[32], 8 );
			i3 = w_i * 3;
			j3 = w_j * 3;
			r2 = 0.0;
			for( j = 0; j < 3; j++ ) { dr[j] = xqm[i3+j] - xmm[j3+j]; r2 += dr[j] * dr[j]; }
			ss = w_r / sqrt( r2 );
			ss = ss * ss * ss * ss * ss * ss;
			out += w_e * ss * ( ss - 2.0 ) * w_f;

			tt = w_e * ss * w_f;
			df = 12.0 * tt * ( 1.0 - ss ) / r2;
			for( j = 0; j < 3; j++ ) xgg[i3+j] += df * dr[j] * __fg;

			r4        = r2 * r2;
			d2f       = 96.0 * tt * ( 1.750 * ss - 1.0 ) / r4;
			ii        = i3 * ( i3 + 1 ) / 2 + i3;
			xhh[ii]   += ( d2f * dr[0] * dr[0] + df ) * __fh;
			ii        += i3 + 1;
			xhh[ii]   += ( d2f * dr[0] * dr[1] ) * __fh;
			xhh[ii+1] += ( d2f * dr[1] * dr[1] + df ) * __fh;
			ii        += i3 + 2;
			xhh[ii]   += ( d2f * dr[0] * dr[2] ) * __fh;
			xhh[ii+1] += ( d2f * dr[1] * dr[2] ) * __fh;
			xhh[ii+2] += ( d2f * dr[2] * dr[2] + df ) * __fh;
		}
		fclose( fd );

		for( i = 0; i < nqm; i++ ) PyList_SetItem( ogg, i, PyFloat_FromDouble( xgg[i] ) );
		for( i = 0;  i < nh; i++ ) PyList_SetItem( ohh, i, PyFloat_FromDouble( xhh[i] ) );

		free( xqm ); free( xmm ); free( xgg ); free( xhh );

	}
	return( Py_BuildValue( "d", out * __fe ) );
}


static struct PyMethodDef methods [] = {
    { "func", (PyCFunction)get_func, METH_VARARGS },
    { "grad", (PyCFunction)get_grad, METH_VARARGS },
    { "hess", (PyCFunction)get_hess, METH_VARARGS },
	{ 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_qmmm",
	NULL,
	-1,
	methods
};

PyMODINIT_FUNC PyInit__qmmm( void ) {
	PyObject    *my_module;
	my_module = PyModule_Create( &moddef );
	return( my_module );
}

#else

void init_qmmm( void ) {
	PyObject    *my_module;
	my_module = Py_InitModule( "_qmmm", methods );
}

#endif

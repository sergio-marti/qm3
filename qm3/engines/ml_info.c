#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double __dist( long i3, long j3, double *xyz ) {
	return( sqrt(
			( xyz[i3]   - xyz[j3]   ) * ( xyz[i3]   - xyz[j3]   ) +
			( xyz[i3+1] - xyz[j3+1] ) * ( xyz[i3+1] - xyz[j3+1] ) +
			( xyz[i3+2] - xyz[j3+2] ) * ( xyz[i3+2] - xyz[j3+2] ) ) );
}


static PyObject* _coul_info( PyObject *self, PyObject *args ) {
    PyObject	*oout, *ocrd, *onum;
	long		siz, dim, i, j, i3, ww;
	double		*num, *xyz, *out;
    if( PyArg_ParseTuple( args, "OO", &onum, &ocrd ) ) {
		siz = PyList_Size( onum );
		num = (double*) malloc( siz * sizeof( double ) );
		xyz = (double*) malloc( 3 * siz * sizeof( double ) );
		for( i = 0; i < siz; i++ ) num[i] = PyFloat_AsDouble( PyList_GetItem( onum, i ) );
		for( i = 0; i < 3 * siz; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( ocrd, i ) );
		dim = siz * ( siz + 1 ) / 2;
		out = (double*) malloc( dim * sizeof( double ) );
		for( i = 0; i < siz; i++ ) {
			i3 = i * 3;
			ww = i * siz - ( ( i - 1 ) * i ) / 2;
			out[ww] = 0.5 * pow( num[i], 2.4 );
			for( j = i + 1; j < siz; j++ ) {
				out[ww+j-i] = num[i] * num[j] / __dist( i3, j * 3, xyz );
			}
		}
		oout = PyList_New( dim );
		for( i = 0; i < dim; i++ ) PyList_SetItem( oout, i, PyFloat_FromDouble( out[i] ) );
		free( num ); free( xyz ); free( out );
		return( oout );
    } else { Py_INCREF( Py_None ); return( Py_None ); }
}


static PyObject* _coul_jaco( PyObject *self, PyObject *args ) {
    PyObject	*oout, *ocrd, *onum;
	long		siz, dim, i, j, i3, j3, row;
	double		*num, *xyz, *out, dr[3], r2, zz;
    if( PyArg_ParseTuple( args, "OO", &onum, &ocrd ) ) {
		siz = PyList_Size( onum );
		num = (double*) malloc( siz * sizeof( double ) );
		xyz = (double*) malloc( 3 * siz * sizeof( double ) );
		for( i = 0; i < siz; i++ ) num[i] = PyFloat_AsDouble( PyList_GetItem( onum, i ) );
		for( i = 0; i < 3 * siz; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( ocrd, i ) );
		dim = siz * ( siz + 1 ) / 2;
		dim *= siz * 3;
		out = (double*) malloc( dim * sizeof( double ) );
		for( i = 0; i < dim; i++ ) out[i] = 0.0;
		row = 0;
		for( i = 0; i < siz; i++ ) {
			i3 = i * 3;
			for( j = i; j < siz; j++ ) {
				if( j != i ) {
					j3 = j * 3;
					dr[0] = xyz[i3]   - xyz[j3];
					dr[1] = xyz[i3+1] - xyz[j3+1];
					dr[2] = xyz[i3+2] - xyz[j3+2];
					r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
					zz = num[i] * num[j] / ( r2 * sqrt( r2 ) );
					out[row+i3]   = - dr[0] * zz;
					out[row+i3+1] = - dr[1] * zz;
					out[row+i3+2] = - dr[2] * zz;
					out[row+j3]   =   dr[0] * zz;
					out[row+j3+1] =   dr[1] * zz;
					out[row+j3+2] =   dr[2] * zz;
				}
				row += siz * 3;
			}
		}
		oout = PyList_New( dim );
		for( i = 0; i < dim; i++ ) PyList_SetItem( oout, i, PyFloat_FromDouble( out[i] ) );
		free( num ); free( xyz ); free( out );
		return( oout );
    } else { Py_INCREF( Py_None ); return( Py_None ); }
}


double __cosang( long i3, long j3, long k3, double *xyz ) {
	double	vij[3], vik[3], mij, mik;
	vij[0] = xyz[j3]   - xyz[i3];
	vij[1] = xyz[j3+1] - xyz[i3+1];
	vij[2] = xyz[j3+2] - xyz[i3+2];
	mij    = vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2];
	vik[0] = xyz[k3]   - xyz[i3];
	vik[1] = xyz[k3+1] - xyz[i3+1];
	vik[2] = xyz[k3+2] - xyz[i3+2];
	mik    = vik[0] * vik[0] + vik[1] * vik[1] + vik[2] * vik[2];
	return( ( vij[0] * vik[0] + vij[1] * vik[1] + vij[2] * vik[2] ) / sqrt( mij * mik ) );
}


double __fcut( double dst, double cut ) {
	if( dst > cut ) { return( 0.0 ); } else { return( 0.5 * ( cos( M_PI * dst / cut ) + 1.0 ) ); }
}


static PyObject* _acsf_info( PyObject *self, PyObject *args ) {
    PyObject	*oout, *ocrd, *oeta2, *oeta5;
	long		siz, dim, i, i3, j, j3, k, k3, l, dd, neta2, neta5;
	double		*xyz, *out, cutx, dse5, pre5, *eta2, *eta5, fij, fik, dij, dik;
    if( PyArg_ParseTuple( args, "dOddOO", &cutx, &oeta2, &dse5, &pre5, &oeta5, &ocrd ) ) {
		siz   = PyList_Size( ocrd ) / 3;
		neta2 = PyList_Size( oeta2 );
		neta5 = PyList_Size( oeta5 );
		dim   = neta2 + neta5;
		xyz   = (double*) malloc( 3 * siz * sizeof( double ) );
		out   = (double*) malloc( siz * dim * sizeof( double ) );
		eta2  = (double*) malloc( neta2 * sizeof( double ) );
		eta5  = (double*) malloc( neta5 * sizeof( double ) );
		for( i = 0; i < siz * 3; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( ocrd, i ) );
		for( i = 0; i < siz * dim; i++ ) out[i] = 0.0;
		for( i = 0; i < neta2; i++ ) eta2[i] = PyFloat_AsDouble( PyList_GetItem( oeta2, i ) );
		for( i = 0; i < neta5; i++ ) eta5[i] = PyFloat_AsDouble( PyList_GetItem( oeta5, i ) );
		for( i = 0; i < siz; i++ ) {
			i3 = i * 3;
			dd = i * dim;
			for( j = 0; j < siz; j++ ) if( j != i ) {
				j3  = j * 3;
				dij = __dist( i3, j3, xyz );
				fij = __fcut( dij, cutx );
				if( fij > 0.0 ) {
					for( l = 0; l < neta2; l++ ) { out[dd+l] += fij * exp( - eta2[l] * dij * dij ); }
					for( k = 0; k < siz; k++ ) if( k != j && k != i ) {
						k3  = k * 3;
						dik = __dist( i3, k3, xyz );
						fik = __fcut( dik, cutx );
						if( fik > 0.0 ) {
							for( l = 0; l < neta5; l++ ) {
								out[dd+neta2+l] += pre5 * fij * fik * pow( 1.0 + __cosang( i3, j3, k3, xyz ), dse5 ) * exp( - eta5[l] * ( dij * dij + dik * dik ) );
							}
						}
					}
				}
			}
		}
		oout = PyList_New( siz * dim );
		for( i = 0; i < siz * dim; i++ ) PyList_SetItem( oout, i, PyFloat_FromDouble( out[i] ) );
		free( eta2 ); free( eta5 ); free( xyz ); free( out );
		return( oout );
    } else { Py_INCREF( Py_None ); return( Py_None ); }
}


static PyObject* _rada_info( PyObject *self, PyObject *args ) {
    PyObject	*oout, *ocrd, *oref;
	long		*ref, siz, i, j, i3, j3, a3, b3, k;
	double		*xyz, rr[3], e1[3], e2[3], e3[3], dd;
    if( PyArg_ParseTuple( args, "OO", &oref, &ocrd ) ) {
		siz = PyList_Size( oref ) / 2;
		xyz = (double*) malloc( 3 * siz * sizeof( double ) );
		ref = (long*) malloc( 2 * siz * sizeof( long ) );
		for( i = 0; i < siz * 2; i++ ) ref[i] = PyLong_AsLong( PyList_GetItem( oref, i ) );
		for( i = 0; i < siz * 3; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( ocrd, i ) );
		k = 0;
		oout = PyList_New( siz * ( siz - 1 ) * 4 );
		for( i = 0; i < siz; i++ ) {
			i3    = i * 3;
			a3    = ref[2*i];
			b3    = ref[2*i+1];

			e1[0] = xyz[a3]   - xyz[i3];
			e1[1] = xyz[a3+1] - xyz[i3+1];
			e1[2] = xyz[a3+2] - xyz[i3+2];
			dd    = sqrt( e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2] );
			e1[0] /= dd; e1[1] /= dd; e1[2] /= dd;

			e2[0] = xyz[b3]   - xyz[i3];
			e2[1] = xyz[b3+1] - xyz[i3+1];
			e2[2] = xyz[b3+2] - xyz[i3+2];
			dd    = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
			e2[0] -= dd * e1[0];
			e2[1] -= dd * e1[1];
			e2[2] -= dd * e1[2];
			dd    = sqrt( e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2] );
			e2[0] /= dd; e2[1] /= dd; e2[2] /= dd;

            e3[0] = e1[1] * e2[2] - e1[2] * e2[1];
            e3[1] = e1[2] * e2[0] - e1[0] * e2[2];
            e3[2] = e1[0] * e2[1] - e1[1] * e2[0];

			for( j = 0; j < siz; j++ ) if( j != i ) {
				j3       = j * 3;
				rr[0]    = xyz[j3]   - xyz[i3];
				rr[1]    = xyz[j3+1] - xyz[i3+1];
				rr[2]    = xyz[j3+2] - xyz[i3+2];
				dd       = sqrt( rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2] );
				PyList_SetItem( oout, k++, PyFloat_FromDouble( 1. / dd ) );
				PyList_SetItem( oout, k++, PyFloat_FromDouble( ( rr[0] * e1[0] + rr[1] * e1[1] + rr[2] * e1[2] ) / dd ) );
				PyList_SetItem( oout, k++, PyFloat_FromDouble( ( rr[0] * e2[0] + rr[1] * e2[1] + rr[2] * e2[2] ) / dd ) );
				PyList_SetItem( oout, k++, PyFloat_FromDouble( ( rr[0] * e3[0] + rr[1] * e3[1] + rr[2] * e3[2] ) / dd ) );
			}
		}
		free( ref ); free( xyz );
		return( oout );
    } else { Py_INCREF( Py_None ); return( Py_None ); }
}


static struct PyMethodDef methods [] = {
	{ "coul_info", (PyCFunction)_coul_info, METH_VARARGS },
	{ "coul_jaco", (PyCFunction)_coul_jaco, METH_VARARGS },
	{ "acsf_info", (PyCFunction)_acsf_info, METH_VARARGS },
	{ "rada_info", (PyCFunction)_rada_info, METH_VARARGS },
	{ 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_ml_info",
	NULL,
	-1,
	methods
};

PyMODINIT_FUNC PyInit__ml_info( void ) {
	PyObject    *my_module;
	my_module = PyModule_Create( &moddef );
	return( my_module );
}

#else

void init_ml_info( void ) {
	PyObject    *my_module;
	my_module = Py_InitModule( "_ml_info", methods );
}

#endif

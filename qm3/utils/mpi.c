#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"


static PyObject* __init( PyObject *self ){
    int			cpu, pid;

    MPI_Init( NULL, NULL );
    MPI_Comm_rank( MPI_COMM_WORLD, &pid );
    MPI_Comm_size( MPI_COMM_WORLD, &cpu );
    return( Py_BuildValue( "(i,i)", pid, cpu ) );
}


static PyObject* __stop( PyObject *self ) {
    MPI_Finalize();
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __barrier( PyObject *self ) {
    MPI_Barrier( MPI_COMM_WORLD );
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __send_r8( PyObject *self, PyObject *args ) {
    PyObject	*o_dat;
    double		*r_dat;
    int			i, who, siz;

    if( PyArg_ParseTuple( args, "iO", &who, &o_dat ) ) {
    	siz = (int) PyList_Size( o_dat );
    	r_dat = (double*) malloc( siz * sizeof( double ) );
    	for( i = 0; i < siz; i++ ) r_dat[i] = PyFloat_AsDouble( PyList_GetItem( o_dat, i ) );
    	MPI_Send( r_dat, siz, MPI_DOUBLE_PRECISION, who, 0, MPI_COMM_WORLD );
    	free( r_dat );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __recv_r8( PyObject *self, PyObject *args ) {
    int			i, siz, who;
    double		*r_dat;
    PyObject	*o_dat;
    MPI_Status  sts;

    if( PyArg_ParseTuple( args, "ii", &who, &siz ) ) {
    	r_dat = (double*) malloc( siz * sizeof( double ) );
    	MPI_Recv( r_dat, siz, MPI_DOUBLE_PRECISION, who, MPI_ANY_TAG, MPI_COMM_WORLD, &sts );
    	o_dat = PyList_New( siz );
    	for( i = 0; i < siz; i++ ) PyList_SetItem( o_dat, i, PyFloat_FromDouble( r_dat[i] ) );
    	free( r_dat );
    	return( o_dat );
    } else { Py_INCREF( Py_None ); return( Py_None ); }
}


static PyObject* __send_i4( PyObject *self, PyObject *args ) {
    PyObject	*o_dat;
    int			*r_dat;
    int			i, who, siz;

    if( PyArg_ParseTuple( args, "iO", &who, &o_dat ) ) {
    	siz = (int) PyList_Size( o_dat );
    	r_dat = (int*) malloc( siz * sizeof( int ) );
    	for( i = 0; i < siz; i++ ) r_dat[i] = (int) PyLong_AsSize_t( PyList_GetItem( o_dat, i ) );
    	MPI_Send( r_dat, siz, MPI_INTEGER, who, 0, MPI_COMM_WORLD );
    	free( r_dat );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __recv_i4( PyObject *self, PyObject *args ) {
    int			i, siz, who;
    int			*r_dat;
    PyObject	*o_dat;
    MPI_Status  sts;

    if( PyArg_ParseTuple( args, "ii", &who, &siz ) ) {
    	r_dat = (int*) malloc( siz * sizeof( int ) );
    	MPI_Recv( r_dat, siz, MPI_INTEGER, who, MPI_ANY_TAG, MPI_COMM_WORLD, &sts );
    	o_dat = PyList_New( siz );
    	for( i = 0; i < siz; i++ ) PyList_SetItem( o_dat, i, PyLong_FromSize_t( r_dat[i] ) );
    	free( r_dat );
    	return( o_dat );
    } else { Py_INCREF( Py_None ); return( Py_None ); }
}


static struct PyMethodDef methods [] = {
    { "init", (PyCFunction)__init, METH_NOARGS },
    { "stop", (PyCFunction)__stop, METH_NOARGS },
    { "barrier", (PyCFunction)__barrier, METH_NOARGS },
    { "send_r8", (PyCFunction)__send_r8, METH_VARARGS },
    { "recv_r8", (PyCFunction)__recv_r8, METH_VARARGS },
    { "send_i4", (PyCFunction)__send_i4, METH_VARARGS },
    { "recv_i4", (PyCFunction)__recv_i4, METH_VARARGS },
    { 0, 0, 0 }
};


static struct PyModuleDef moddef = {
    PyModuleDef_HEAD_INIT,
    "_mpi",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__mpi( void ) {
    PyObject    *my_module;
    my_module = PyModule_Create( &moddef );
    return( my_module );
}

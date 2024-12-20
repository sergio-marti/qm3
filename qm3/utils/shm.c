#define	PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/ipc.h>
#include <sys/shm.h>


// for i in `ipcs -m | grep $USER | awk '{print $2}'`; do ipcrm -m $i; done

// sysctl -w kernel.shmmni=32768


static PyObject* __alloc( PyObject *self, PyObject *args ){
    int			xid;
    Py_ssize_t	i, siz;
    char		*mem;

    if( PyArg_ParseTuple( args, "l", &siz ) ) {
    	xid = shmget( IPC_PRIVATE, siz, IPC_CREAT | 0600 );
    	mem = (char*) shmat( xid, 0, 0 );
    	for( i = 0; i < siz; i++ ) mem[i] = 0;
    	shmdt( mem );
    	return( Py_BuildValue( "i", xid ) );
    } else {
    	Py_INCREF( Py_None );
    	return( Py_None );
    }
}


static PyObject* __clean( PyObject *self, PyObject *args ){
    int		xid;

    if( PyArg_ParseTuple( args, "i", &xid ) ) {
    	shmctl( xid, IPC_RMID, 0 );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __read( PyObject *self, PyObject *args ){
    PyObject	*out;
    int			xid;
	Py_ssize_t	siz;
    char		*mem;

    if( PyArg_ParseTuple( args, "il", &xid, &siz ) ) {
    	mem = (char*) shmat( xid, 0, 0 );
    	out = Py_BuildValue( "y#", (unsigned char*) mem, siz );
    	shmdt( mem );
    	return( out );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __write( PyObject *self, PyObject *args ){
    int			xid;
	Py_ssize_t	i, siz;
    char		*buf, *mem;

    if( PyArg_ParseTuple( args, "iy#", &xid, &buf, &siz ) ) {
    	mem = (char*) shmat( xid, 0, 0 );
    	for( i = 0; i < siz; i++ ) mem[i] = buf[i];
    	shmdt( mem );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __read_r8( PyObject *self, PyObject *args ){
    int			xid;
    Py_ssize_t	i, siz;
    char		*mem;
    PyObject	*out;
	double		tmp;

    if( PyArg_ParseTuple( args, "il", &xid, &siz ) ) {
    	mem = (char*) shmat( xid, 0, 0 );
    	out = PyList_New( siz );
    	for( i = 0; i < siz; i++ ) {
			memcpy( &tmp, &mem[i*8], 8 );
			PyList_SetItem( out, i, PyFloat_FromDouble( tmp ) );
		}
    	shmdt( mem );
    	return( out );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __write_r8( PyObject *self, PyObject *args ){
    int			xid;
    Py_ssize_t	i, siz;
    char		*mem;
    PyObject	*vec;
	double		tmp;

    if( PyArg_ParseTuple( args, "iO", &xid, &vec ) ) {
    	siz = PyList_Size( vec );
    	mem = (char*) shmat( xid, 0, 0 );
    	for( i = 0; i < siz; i++ ) {
			tmp = PyFloat_AsDouble( PyList_GetItem( vec, i ) );
			memcpy( &mem[i*8], &tmp, 8 );
		}
    	shmdt( mem );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static struct PyMethodDef methods [] = {
    { "alloc", (PyCFunction)__alloc, METH_VARARGS },
    { "clean", (PyCFunction)__clean, METH_VARARGS },
    { "read",  (PyCFunction)__read,  METH_VARARGS },
    { "write", (PyCFunction)__write, METH_VARARGS },
    { "read_r8",  (PyCFunction)__read_r8,  METH_VARARGS },
    { "write_r8", (PyCFunction)__write_r8, METH_VARARGS },
    { 0, 0, 0 }
};


static struct PyModuleDef moddef = {
    PyModuleDef_HEAD_INIT,
    "_shm",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__shm( void ) {
    PyObject    *my_module;
    my_module = PyModule_Create( &moddef );
    return( my_module );
}

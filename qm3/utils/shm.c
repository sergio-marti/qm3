#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/ipc.h>
#include <sys/shm.h>


static PyObject* __alloc( PyObject *self, PyObject *args ){
	int		i, siz, xid;
	char	*mem;

	if( PyArg_ParseTuple( args, "i", &siz ) ) {
		xid = shmget( IPC_PRIVATE, siz, IPC_CREAT | 0600 );
		mem = shmat( xid, 0, 0 );
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
	int			xid, siz;
	char		*mem;

	if( PyArg_ParseTuple( args, "ii", &xid, &siz ) ) {
		mem = (char*) shmat( xid, 0, 0 );
#if PY_MAJOR_VERSION >= 3
		out = Py_BuildValue( "y#", (unsigned char*) mem, siz );
#else
		out = Py_BuildValue( "s#", mem, siz );
#endif
		shmdt( mem );
		return( out );
	}
	Py_INCREF( Py_None );
	return( Py_None );
}


static PyObject* __write( PyObject *self, PyObject *args ){
	int		i, xid, siz;
	char	*buf, *mem;

	if( PyArg_ParseTuple( args, "is#", &xid, &buf, &siz ) ) {
		mem = (char*) shmat( xid, 0, 0 );
		for( i = 0; i < siz; i++ ) mem[i] = buf[i];
		shmdt( mem );
	}
	Py_INCREF( Py_None );
	return( Py_None );
}


static PyObject* __read_r8( PyObject *self, PyObject *args ){
	int			i, xid, siz;
	double		*mem;
	PyObject	*out;

	if( PyArg_ParseTuple( args, "ii", &xid, &siz ) ) {
		mem = (double*) shmat( xid, 0, 0 );
		out = PyList_New( siz );
		for( i = 0; i < siz; i++ ) PyList_SetItem( out, i, PyFloat_FromDouble( mem[i] ) );
		shmdt( mem );
		return( out );
	}
	Py_INCREF( Py_None );
	return( Py_None );
}


static PyObject* __write_r8( PyObject *self, PyObject *args ){
	int			i, xid, siz;
	double		*mem;
	PyObject	*vec;

	if( PyArg_ParseTuple( args, "iO", &xid, &vec ) ) {
		siz = (int) PyList_Size( vec );
		mem = (double*) shmat( xid, 0, 0 );
		for( i = 0; i < siz; i++ ) mem[i] = PyFloat_AsDouble( PyList_GetItem( vec, i ) );
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


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_shm",
	NULL,
	-1,
	methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC PyInit__shm( void ) {
	PyObject    *my_module;
	my_module = PyModule_Create( &moddef );
	return( my_module );
}
#else
void init_shm( void ) {
	PyObject    *my_module;
	my_module = Py_InitModule( "_shm", methods );
}
#endif

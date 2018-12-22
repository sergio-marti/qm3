#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


typedef struct {
	PyObject_HEAD
} OBJECT;


static int __init( OBJECT *self, PyObject *args, PyObject *kwds ) {
	return( 0 );
}


static PyObject* __new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
	OBJECT		*self;

	self = (OBJECT*) type->tp_alloc( type, 0 );
	return( (PyObject*) self ) ;
}


static void __dealloc( OBJECT *self ) {
	Py_TYPE( self )->tp_free( (PyObject*) self );
}


static PyObject* METHOD( PyObject *self, PyObject *args ) {
	OBJECT		*obj = NULL;

	obj = (OBJECT*) self;
	Py_INCREF( Py_None );
	return( Py_None );
}


static struct PyMethodDef __methods [] = {
    { "METHOD", (PyCFunction)METHOD, METH_VARARGS },
	{ 0, 0, 0 }
};


static struct PyMemberDef __members [] = {
        { 0, 0, 0, 0 }
};


// --------------------------------------------------------------------------------------


static struct PyMethodDef methods [] = {
	{ 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3

static PyTypeObject TYPE = {
	PyVarObject_HEAD_INIT( NULL, 0 )
	.tp_name = "CLASS",
	.tp_doc = "CLASS DESCRIPTION",
	.tp_basicsize = sizeof( OBJECT ),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	.tp_new = __new,
	.tp_init = (initproc) __init,
	.tp_dealloc = (destructor) __dealloc,
	.tp_members = __members,
	.tp_methods = __methods,
};

static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_MODULE",
	NULL,
	-1,
	methods
};

PyMODINIT_FUNC PyInit__MODULE( void ) {
	PyObject    *my_module;

	my_module = PyModule_Create( &moddef );
	PyType_Ready( &TYPE );
    Py_INCREF( &TYPE );
    PyModule_AddObject( my_module, "CLASS", (PyObject *) &TYPE );
	return( my_module );
}

#else

static PyTypeObject TYPE = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_MODULE",               // tp_name
    sizeof( OBJECT ),        // tp_basicsize
    0,                       // tp_itemsize
    (destructor)__dealloc,   // tp_dealloc
    0,                       // tp_print
    0,                       // tp_getattr
    0,                       // tp_setattr
    0,                       // tp_compare
    0,                       // tp_repr
    0,                       // tp_as_number
    0,                       // tp_as_sequence
    0,                       // tp_as_mapping
    0,                       // tp_hash
    0,                       // tp_call
    0,                       // tp_str
    0,                       // tp_getattro
    0,                       // tp_setattro
    0,                       // tp_as_buffer
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE, // tp_flags
    "CLASS DESCRIPTION",
                             // tp_doc
    0,                       // tp_traverse
    0,                       // tp_clear
    0,                       // tp_richcompare
    0,                       // tp_weaklistoffset
    0,                       // tp_iter
    0,                       // tp_iternext
    __methods,               // tp_methods
    __members,               // tp_members
    0,                       // tp_getset
    0,                       // tp_base
    0,                       // tp_dict
    0,                       // tp_descr_get
    0,                       // tp_descr_set
    0,                       // tp_dictoffset
    (initproc)__init,        // tp_init
    0,                       // tp_alloc
    __new,                   // tp_new
};

void init_MODULE( void ) {
	PyObject    *my_module;

	my_module = Py_InitModule( "_MODULE", methods );
	PyType_Ready( &TYPE );
    Py_INCREF( &TYPE );
    PyModule_AddObject( my_module, "CLASS", (PyObject *) &TYPE );
}

#endif

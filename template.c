#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


typedef struct {
	PyObject_HEAD
} @object@;


static int __init( @object@ *self, PyObject *args, PyObject *kwds ) {
	return( 0 );
}


static PyObject* __new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
	@object@		*self;

	self = (@object@*) type->tp_alloc( type, 0 );
	return( (PyObject*) self ) ;
}


static void __dealloc( @object@ *self ) {
	Py_TYPE( self )->tp_free( (PyObject*) self );
}


static PyObject* @method@( PyObject *self, PyObject *args ) {
	PyObject		*o_;
	@object@		*obj = NULL;

	obj = (@object@*) self;
	if( PyArg_ParseTuple( args, "O", &o_ ) ) {
		return( Py_BuildValue( "O", ... ) );
	}
	Py_INCREF( Py_None );
	return( Py_None );
}


static struct PyMethodDef __methods [] = {
    { "@method@", (PyCFunction)@method@, METH_VARARGS },
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

static PyTypeObject @type@ = {
	PyVarObject_HEAD_INIT( NULL, 0 )
	.tp_name = "@class@",
	.tp_doc = "@class@ description",
	.tp_basicsize = sizeof( @object@ ),
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
	"_@module@",
	NULL,
	-1,
	methods
};

PyMODINIT_FUNC PyInit__@module@( void ) {
	PyObject    *my_module;

	my_module = PyModule_Create( &moddef );
	PyType_Ready( &@type@ );
    Py_INCREF( &@type@ );
    PyModule_AddObject( my_module, "@class@", (PyObject *) &@type@ );
	return( my_module );
}

#else

static PyTypeObject @type@ = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "@class@",                 // tp_name
    sizeof( @object@ ),        // tp_basicsize
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
    "@class@ description",
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

void init_@module@( void ) {
	PyObject    *my_module;

	my_module = Py_InitModule( "_@module@", methods );
	PyType_Ready( &@type@ );
    Py_INCREF( &@type@ );
    PyModule_AddObject( my_module, "@class@", (PyObject *) &@type@ );
}

#endif

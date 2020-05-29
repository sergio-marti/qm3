#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>


typedef void (*Tfunc_initialize)();
typedef void (*Tfunc_update_coor)(double*);
typedef void (*Tfunc_get_func)(double*, double*);
typedef void (*Tfunc_get_grad)(double*, double*, double*);
typedef void (*Tfunc_facade)(double*);


typedef struct {
    PyObject_HEAD
    void				*lib;
    Tfunc_initialize	qm3_initialize;
    Tfunc_update_coor	qm3_update_coor;
    Tfunc_update_coor	qm3_update_chrg;
    Tfunc_get_func		qm3_get_func;
    Tfunc_get_grad		qm3_get_grad;
    Tfunc_facade		qm3_facade;
} oDynamo;


static int __init( oDynamo *self, PyObject *args, PyObject *kwds ) {
    char	*path;

//    path = getenv( "DYNAMO" );
    if( PyArg_ParseTuple( args, "s", &path ) ) {
    	self->lib = dlopen( path, RTLD_NOW | RTLD_GLOBAL );
    	if( self->lib == NULL ) return( -1 );
    	self->qm3_initialize  = (Tfunc_initialize) dlsym( self->lib, "qm3_initialize_" );
    	self->qm3_update_coor = (Tfunc_update_coor) dlsym( self->lib, "qm3_update_coor_" );
    	self->qm3_update_chrg = (Tfunc_update_coor) dlsym( self->lib, "qm3_update_chrg_" );
    	self->qm3_get_func    = (Tfunc_get_func) dlsym( self->lib, "qm3_get_func_" );
    	self->qm3_get_grad    = (Tfunc_get_grad) dlsym( self->lib, "qm3_get_grad_" );
    	self->qm3_facade      = (Tfunc_facade) dlsym( self->lib, "qm3_facade_" );
    	(*(self->qm3_initialize))();
    }
    return( 0 );
}


static PyObject* __new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
    oDynamo		*self;

    self = (oDynamo*) type->tp_alloc( type, 0 );
    self->lib = NULL;
    return( (PyObject*) self ) ;
}


static void __dealloc( oDynamo *self ) {
    if( self->lib ) dlclose( self->lib );
    Py_TYPE( self )->tp_free( (PyObject*) self );
}


static PyObject* qm3_facade( PyObject *self, PyObject *args ) {
    PyObject	*omol, *ocrd;
    long		i, natm;
    double		*coor;
    oDynamo		*obj = NULL;

    obj = (oDynamo*) self;
    if( PyArg_ParseTuple( args, "O", &omol ) ) {
    	natm = PyLong_AsLong( PyObject_GetAttrString( omol, "natm" ) );
    	ocrd = PyObject_GetAttrString( omol, "coor" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	for( i = 0; i < 3 * natm; i++ )
    		coor[i] = PyFloat_AsDouble( PyList_GetItem( ocrd, i ) );
    	Py_DECREF( ocrd );
//    	(*(obj->qm3_update_coor))(coor);
    	(*(obj->qm3_facade))(coor);
    	free( coor );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* qm3_update_coor( PyObject *self, PyObject *args ) {
    PyObject	*omol, *ocrd;
    long		i, natm;
    double		*coor;
    oDynamo		*obj = NULL;

    obj = (oDynamo*) self;
    if( PyArg_ParseTuple( args, "O", &omol ) ) {
    	natm = PyLong_AsLong( PyObject_GetAttrString( omol, "natm" ) );
    	ocrd = PyObject_GetAttrString( omol, "coor" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	for( i = 0; i < 3 * natm; i++ )
    		coor[i] = PyFloat_AsDouble( PyList_GetItem( ocrd, i ) );
    	Py_DECREF( ocrd );
    	(*(obj->qm3_update_coor))(coor);
    	free( coor );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* qm3_update_chrg( PyObject *self, PyObject *args ) {
    PyObject	*omol, *ochg;
    long		i, natm;
    double		*chrg;
    oDynamo		*obj = NULL;

    obj = (oDynamo*) self;
    if( PyArg_ParseTuple( args, "O", &omol ) ) {
    	natm = PyLong_AsLong( PyObject_GetAttrString( omol, "natm" ) );
    	ochg = PyObject_GetAttrString( omol, "chrg" );
    	chrg = (double*) malloc( natm * sizeof( double ) );
    	for( i = 0; i < natm; i++ )
    		chrg[i] = PyFloat_AsDouble( PyList_GetItem( ochg, i ) );
    	Py_DECREF( ochg );
    	(*(obj->qm3_update_chrg))(chrg);
    	free( chrg );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* qm3_get_func( PyObject *self, PyObject *args ) {
    PyObject	*omol, *otmp;
    long		i, natm;
    double		*coor, func, tmp;
    oDynamo		*obj = NULL;

    obj = (oDynamo*) self;
    if( PyArg_ParseTuple( args, "O", &omol ) ) {
    	natm = PyLong_AsLong( PyObject_GetAttrString( omol, "natm" ) );

    	otmp = PyObject_GetAttrString( omol, "coor" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	for( i = 0; i < 3 * natm; i++ )
    		coor[i] = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
    	Py_DECREF( otmp );

//    	(*(obj->qm3_update_coor))(coor);
    	(*(obj->qm3_get_func))(coor, &func);

    	otmp = PyObject_GetAttrString( omol, "func" );
    	tmp = PyFloat_AsDouble( otmp );
    	Py_DECREF( otmp );
    	PyObject_SetAttrString( omol, "func", PyFloat_FromDouble( tmp + func ) );

    	free( coor );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* qm3_get_grad( PyObject *self, PyObject *args ) {
    PyObject	*omol, *otmp;
    long		i, natm;
    double		*coor, *grad, func, tmp;
    oDynamo		*obj = NULL;

    obj = (oDynamo*) self;
    if( PyArg_ParseTuple( args, "O", &omol ) ) {
    	natm = PyLong_AsLong( PyObject_GetAttrString( omol, "natm" ) );

    	otmp = PyObject_GetAttrString( omol, "coor" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	grad = (double*) malloc( 3 * natm * sizeof( double ) );
    	for( i = 0; i < 3 * natm; i++ )
    		coor[i] = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
    	Py_DECREF( otmp );

//    	(*(obj->qm3_update_coor))(coor);
    	(*(obj->qm3_get_grad))(coor, &func, grad);

    	otmp = PyObject_GetAttrString( omol, "func" );
    	tmp = PyFloat_AsDouble( otmp );
    	Py_DECREF( otmp );
    	PyObject_SetAttrString( omol, "func", PyFloat_FromDouble( tmp + func ) );

    	otmp = PyObject_GetAttrString( omol, "grad" );
    	for( i = 0; i < 3 * natm; i++ ) {
    		tmp = PyFloat_AsDouble( PyList_GetItem( otmp, i ) ) + grad[i];
    		PyList_SetItem( otmp, i, PyFloat_FromDouble( tmp ) );
    	}
    	Py_DECREF( otmp );

    	free( coor ); free( grad );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static struct PyMethodDef __methods [] = {
    { "update_coor", (PyCFunction)qm3_update_coor, METH_VARARGS },
    { "update_chrg", (PyCFunction)qm3_update_chrg, METH_VARARGS },
    { "get_func",    (PyCFunction)qm3_get_func, METH_VARARGS },
    { "get_grad",    (PyCFunction)qm3_get_grad, METH_VARARGS },
    { "facade",      (PyCFunction)qm3_facade, METH_VARARGS },
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

static PyTypeObject tDynamo = {
    PyVarObject_HEAD_INIT( NULL, 0 )
    .tp_name = "dynamo",
    .tp_doc = "fDynamo bindings",
    .tp_basicsize = sizeof( oDynamo ),
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
    "_dynamo",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__dynamo( void ) {
    PyObject    *my_module;

    my_module = PyModule_Create( &moddef );
    PyType_Ready( &tDynamo );
    Py_INCREF( &tDynamo );
    PyModule_AddObject( my_module, "dynamo", (PyObject *) &tDynamo );
    return( my_module );
}

#else

static PyTypeObject tDynamo = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "dynamo",               // tp_name
    sizeof( oDynamo ),        // tp_basicsize
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
    "fDynamo bindings",
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

void init_dynamo( void ) {
    PyObject    *my_module;

    my_module = Py_InitModule( "_dynamo", methods );
    PyType_Ready( &tDynamo );
    Py_INCREF( &tDynamo );
    PyModule_AddObject( my_module, "dynamo", (PyObject *) &tDynamo );
}

#endif

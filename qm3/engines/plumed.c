#include <Python.h>
#include "structmember.h"
#include <stdio.h>
#include <math.h>
#include "Plumed.h"


typedef struct {
    PyObject_HEAD
    int		Natoms, box;
    double	*Masses;
    plumed	pmed;
} oPLUMED;


static int __init( oPLUMED *self, PyObject *args, PyObject *kwds ) {
    PyObject	*mol, *tmp, *pbc = Py_False;
    long		i;
    int			RealPrecision;
    double		MDLengthUnits, MDChargeUnits, Temperature = 300.0, Timestep = 0.001, KbT;

    if( PyArg_ParseTuple( args, "O|dO", &mol, &Timestep, &pbc ) ) {

    	if( pbc == Py_True ) { self->box = 1; } else { self->box = 0; }
    	self->Natoms = (int) PyLong_AsLong( PyObject_GetAttrString( mol, "natm" ) );
    	self->Masses = (double*) malloc( self->Natoms * sizeof( double ) );
    	tmp = PyObject_GetAttrString( mol, "mass" );
    	for( i = 0; i < self->Natoms; i++ ) {
    		self->Masses[i] = PyFloat_AsDouble( PyList_GetItem( tmp, i ) );
    	}
    	Py_DECREF( tmp );

    	RealPrecision = 8;
    	MDLengthUnits = 0.1;
    	MDChargeUnits = 1.602176565e-19;

    	self->pmed = plumed_create();
    	// double
    	plumed_cmd( self->pmed, "setRealPrecision", &(RealPrecision) );
    	// Angstrom to nanometers
    	plumed_cmd( self->pmed, "setMDLengthUnits", &(MDLengthUnits) );
    	// atomic to coulomb
    	plumed_cmd( self->pmed, "setMDChargeUnits", &(MDChargeUnits) );
    	// plumed config name
    	plumed_cmd( self->pmed, "setPlumedDat", "plumed.dat" );
    	// number of atmos
    	plumed_cmd( self->pmed, "setNatoms", &(self->Natoms) );
    	// log filename
    	plumed_cmd( self->pmed, "setLogFile", "plumed.log" );
    	// MD time step (in ps)
    	plumed_cmd( self->pmed, "setTimestep", &(Timestep) );
    	// KbT
//    	KbT = Temperature * 0.008314462145468951;
//    	plumed_cmd( self->pmed, "setKbT", &(KbT) );
    	// go!
    	plumed_cmd( self->pmed, "init", NULL );
    }
    return( 0 );
}


static PyObject* __new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
    oPLUMED	*self;
    int		i;

    self = (oPLUMED*) type->tp_alloc( type, 0 );
    self->Masses = NULL;
    return( (PyObject*) self ) ;
}


static void __dealloc( oPLUMED *self ) {
    plumed_finalize( self->pmed );
    free( self->Masses );
    Py_TYPE( self )->tp_free( (PyObject*) self );
}


static PyObject* __calc( PyObject *self, PyObject *args ) {
    PyObject	*mol, *tmp;
    long		i, stp;
    double		*xyz, *frz, ene, acc = 0.0;
    double		box[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double		vir[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    oPLUMED		*obj = NULL;

    obj = (oPLUMED*) self;
    if( PyArg_ParseTuple( args, "Ol", &mol, &stp ) ) {

    	if( obj->box == 1 ) {
    		// box latice vectors in nanometers...
    		tmp = PyObject_GetAttrString( mol, "boxl" );
    		for( i = 0; i < 3; i++ ) box[3*i+i] = PyFloat_AsDouble( PyList_GetItem( tmp, i ) ) * 0.1;
    		Py_DECREF( tmp );
    	}

    	tmp = PyObject_GetAttrString( mol, "coor" );
    	xyz = (double*) malloc( 3 * obj->Natoms * sizeof( double ) );
    	frz = (double*) malloc( 3 * obj->Natoms * sizeof( double ) );
    	for( i = 0; i < 3 * obj->Natoms; i++ ) {
    		xyz[i] = PyFloat_AsDouble( PyList_GetItem( tmp, i ) );
    		frz[i] = 0.0;
    	}
    	Py_DECREF( tmp );

    	tmp = PyObject_GetAttrString( mol, "func" );
    	ene = PyFloat_AsDouble( tmp );
    	Py_DECREF( tmp );

    	plumed_cmd( obj->pmed, "setStep", &(stp) );
    	plumed_cmd( obj->pmed, "setPositions", xyz );
    	plumed_cmd( obj->pmed, "setMasses", obj->Masses );

    	// charge based constraints...
    	//plumed_cmd( obj->pmed, "setCharges", ??? );

    	if( obj->box == 1 ) plumed_cmd( obj->pmed, "setBox", &box[0] );

    	// energy based constraints...
    	//plumed_cmd( obj->pmed, "setEnergy", &ene );

    	plumed_cmd( obj->pmed, "setForces", frz );
    	plumed_cmd( obj->pmed, "setVirial", &vir[0] );
    	plumed_cmd( obj->pmed, "calc", NULL );
    	plumed_cmd( obj->pmed, "getBias", &acc );

    	PyObject_SetAttrString( mol, "func", PyFloat_FromDouble( ene + acc ) );

    	tmp = PyObject_GetAttrString( mol, "grad" );
    	for( i = 0; i < 3*obj->Natoms; i++ ) {
    		acc = PyFloat_AsDouble( PyList_GetItem( tmp, i ) );
    		PyList_SetItem( tmp, i, PyFloat_FromDouble( acc - frz[i] ) );
    	}
    	Py_DECREF( tmp );

    	free( xyz ); free( frz );
    	
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static struct PyMethodDef PLUMED_methods [] = {
    { "get_grad", (PyCFunction)__calc, METH_VARARGS },
    { 0, 0, 0 }
};

static struct PyMethodDef methods [] = {
    { 0, 0, 0 }
};

static struct PyMemberDef members [] = {
    { 0, 0, 0, 0 }
};


// --------------------------------------------------------------------------------------


#if PY_MAJOR_VERSION >= 3

static PyTypeObject tPLUMED = {
    PyVarObject_HEAD_INIT( NULL, 0 )
    .tp_name = "Plumed",
    .tp_doc = "Plumed object",
    .tp_basicsize = sizeof( oPLUMED ),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = __new,
    .tp_init = (initproc) __init,
    .tp_dealloc = (destructor) __dealloc,
    .tp_members = members,
    .tp_methods = PLUMED_methods,
};

static struct PyModuleDef moddef = {
    PyModuleDef_HEAD_INIT,
    "_plumed",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__plumed( void ) {
    PyObject    *my_module;

    my_module = PyModule_Create( &moddef );
    PyType_Ready( &tPLUMED );
    Py_INCREF( &tPLUMED );
    PyModule_AddObject( my_module, "Plumed", (PyObject *) &tPLUMED );
    return( my_module );
}

#else

static PyTypeObject tPLUMED = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Plumed",                  /* tp_name */
    sizeof( oPLUMED ),         /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)__dealloc,     /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_compare */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Plumed object",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    PLUMED_methods,            /* tp_methods */
    members,                   /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)__init,          /* tp_init */
    0,                         /* tp_alloc */
    __new,                     /* tp_new */
};

void init_plumed( void ) {
    PyObject    *my_module;

    my_module = Py_InitModule( "_plumed", methods );
    PyType_Ready( &tPLUMED );
    Py_INCREF( &tPLUMED );
    PyModule_AddObject( my_module, "Plumed", (PyObject *) &tPLUMED );
}

#endif

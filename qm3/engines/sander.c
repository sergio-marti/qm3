#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "sander.h"


#define min(a,b) (((a)<(b))?(a):(b))


typedef struct {
	PyObject_HEAD
    sander_input		mm_opt;
    qmmm_input_options	qm_opt;
    pot_ene				ene;
	double				box[6];
} oSander;


// --------------------------------------------------------------------------------------


static int sander__init( oSander *self, PyObject *args, PyObject *kwds ) {
	PyObject	*o_mol, *o_sel, *pbc, *tmp;
	double		cut, *crd = NULL;
	long		i, k, n, charge;
	char		*prmtop, *method;

	o_sel = NULL;
	method = NULL;
	charge = 0;

	if( PyArg_ParseTuple( args, "OsdO|Osl", &o_mol, &prmtop, &cut, &pbc, &o_sel, &method, &charge ) ) {

//fprintf(stderr,"[%s]\n",prmtop);
//fprintf(stderr,"[%lf]\n",cut);
//fprintf(stderr,"[%d]\n",pbc == Py_True || pbc == Py_False || pbc == Py_None );
//if( PyList_Check( o_sel ) ) fprintf(stderr,"{%ld}\n",PyList_Size( o_sel ) );
//if( method != NULL ) fprintf(stderr,"{%s}\n",method);
//fprintf(stderr,"{%ld}\n",charge);

		// PBC conditions
		if( pbc == Py_False || pbc == Py_None || method == NULL || ( method != NULL && strncmp( method, "EXTERN", 6 ) == 0 ) ) {
			gas_sander_input( &(self->mm_opt), 6 );
//fprintf(stderr," -- gas phase --\n");
		} else {
			pme_sander_input( &(self->mm_opt) );
//fprintf(stderr," -- PME --\n");
		}
		self->mm_opt.cut = cut;
		self->mm_opt.jfastw = 4;

		// MM or QM/MM calculation
		if( o_sel == NULL || o_sel == Py_None || ( PyList_Check( o_sel ) && PyList_Size( o_sel ) == 0 ) ) {
			self->mm_opt.ifqnt = 0;
//fprintf(stderr," -- MM --\n");
		} else {
			self->mm_opt.ifqnt = 1;
//fprintf(stderr," -- QM/MM --\n");
			qm_sander_input( &(self->qm_opt) );

			strncpy( self->qm_opt.qm_theory, method, min( 12, strlen( method ) ) );
//fprintf(stderr,"(%s)\n",self->qm_opt.qm_theory);

			self->qm_opt.spin = 1;
			self->qm_opt.qmcharge = ((int) charge);
			self->qm_opt.qmcut = cut;

			self->qm_opt.itrmax = 1000;
			self->qm_opt.scfconv = 1.0e-10;

			if( strncmp( method, "EXTERN", 6 ) == 0 ) {
				self->qm_opt.adjust_q = 0;
				self->qm_opt.qmmm_int = 1;
				self->qm_opt.qm_ewald = 0;
			} else {
				self->qm_opt.qmmm_int = 5;
				self->qm_opt.qm_ewald = 1;
			}

			n = min( PyList_Size( o_sel ), MAX_QUANTUM_ATOMS );
			for( i = 0; i < n; i++ ) {
				self->qm_opt.iqmatoms[i] = ((int) PyLong_AsLong( PyList_GetItem( o_sel, i ) )) + 1;
			}
		}

		// BOX (from o_mol )
		tmp = PyObject_GetAttrString( o_mol, "boxl" );
		for( k = 0; k < 3; k++ ) self->box[k] = PyFloat_AsDouble( PyList_GetItem( tmp, k ) );
		self->box[3] = 90.0; self->box[4] = 90.0; self->box[5] = 90.0;
		Py_DECREF( tmp );

		// COORDINATES (from o_mol)
		n   = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
		tmp = PyObject_GetAttrString( o_mol, "coor" );
		crd = (double*) malloc( 3 * n * sizeof( double ) );
		for( i = 0; i < n; i++ ) {
			for( k = 0; k < 3; k++ )
				crd[3*i+k] = PyFloat_AsDouble( PyList_GetItem( tmp, 3*i+k ) );
		}
		Py_DECREF( tmp );

		// go!
    	sander_setup( prmtop, crd, self->box, &(self->mm_opt), &(self->qm_opt) );
		free( crd );
	}
	return( 0 );
}


static PyObject* sander__new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
	oSander		*self;

	self = (oSander*) type->tp_alloc( type, 0 );
	return( (PyObject*) self ) ;
}


static void sander__dealloc( oSander *self ) {
	Py_TYPE( self )->tp_free( (PyObject*) self );
	if( is_setup() ) sander_cleanup();
}


static PyObject* sander__get_func( PyObject *self, PyObject *args ) {
	PyObject	*o_mol, *o_coor, *o_tmp;
	long		i, k, natm;
	double		*coor = NULL, *grad = NULL, tmp;
	oSander		*obj = NULL;

	obj = (oSander*) self;
	if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
		natm = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
		o_coor = PyObject_GetAttrString( o_mol, "coor" );
		coor = (double*) malloc( 3 * natm * sizeof( double ) );
		grad = (double*) malloc( 3 * natm * sizeof( double ) );
		for( i = 0; i < natm; i++ ) {
			for( k = 0; k < 3; k++ ) {
				coor[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_coor, 3*i+k ) );
				grad[3*i+k] = 0.0;
			}
		}
		Py_DECREF( o_coor );
		
		set_positions( coor );
		energy_forces( &(obj->ene), grad );

		o_tmp = PyObject_GetAttrString( o_mol, "func" );
		tmp = PyFloat_AsDouble( o_tmp );
		Py_DECREF( o_tmp );
		PyObject_SetAttrString( o_mol, "func", PyFloat_FromDouble( tmp + obj->ene.tot * 4.184 ) );
		free( coor ); free( grad );
	}
	Py_INCREF( Py_None );
	return( Py_None );
}


static PyObject* sander__get_grad( PyObject *self, PyObject *args ) {
	PyObject	*o_mol, *o_coor, *o_tmp;
	long		i, k, natm;
	double		*coor = NULL, *grad = NULL, tmp;
	oSander		*obj = NULL;

	obj = (oSander*) self;
	if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
		natm = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
		o_coor = PyObject_GetAttrString( o_mol, "coor" );
		coor = (double*) malloc( 3 * natm * sizeof( double ) );
		grad = (double*) malloc( 3 * natm * sizeof( double ) );
		for( i = 0; i < natm; i++ ) {
			for( k = 0; k < 3; k++ ) {
				coor[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_coor, 3*i+k ) );
				grad[3*i+k] = 0.0;
			}
		}
		Py_DECREF( o_coor );
		
		set_positions( coor );
		energy_forces( &(obj->ene), grad );

		o_tmp = PyObject_GetAttrString( o_mol, "func" );
		tmp = PyFloat_AsDouble( o_tmp );
		Py_DECREF( o_tmp );
		PyObject_SetAttrString( o_mol, "func", PyFloat_FromDouble( tmp + obj->ene.tot * 4.184 ) );

		o_tmp = PyObject_GetAttrString( o_mol, "grad" );
		for( i = 0; i < 3 * natm; i++ ) {
			tmp = PyFloat_AsDouble( PyList_GetItem( o_tmp, i ) ) - grad[i] * 4.184;
			PyList_SetItem( o_tmp, i, PyFloat_FromDouble( tmp ) );
		}
		Py_DECREF( o_tmp );

		free( coor ); free( grad );
	}
	Py_INCREF( Py_None );
	return( Py_None );
}


static struct PyMethodDef sander_methods [] = {
    { "get_func", (PyCFunction)sander__get_func, METH_VARARGS },
    { "get_grad", (PyCFunction)sander__get_grad, METH_VARARGS },
	{ 0, 0, 0 }
};


static struct PyMemberDef sander_members [] = {
	{ 0, 0, 0, 0 }
};


// --------------------------------------------------------------------------------------


static struct PyMethodDef methods [] = {
	{ 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3

static PyTypeObject tSander = {
	PyVarObject_HEAD_INIT( NULL, 0 )
	.tp_name = "sander",
	.tp_doc = "Sander python wrapping",
	.tp_basicsize = sizeof( oSander ),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	.tp_new = sander__new,
	.tp_init = (initproc) sander__init,
	.tp_dealloc = (destructor) sander__dealloc,
	.tp_members = sander_members,
	.tp_methods = sander_methods,
};

static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_sander",
	NULL,
	-1,
	methods
};

PyMODINIT_FUNC PyInit__sander( void ) {
	PyObject    *my_module;

	my_module = PyModule_Create( &moddef );
	PyType_Ready( &tSander );
    Py_INCREF( &tSander );
    PyModule_AddObject( my_module, "sander", (PyObject *) &tSander );
	return( my_module );
}

#else

static PyTypeObject tSander = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "sander",                     /* tp_name */
    sizeof( oSander ),            /* tp_basicsize */
    0,                            /* tp_itemsize */
    (destructor)sander__dealloc, /* tp_dealloc */
    0,                            /* tp_print */
    0,                            /* tp_getattr */
    0,                            /* tp_setattr */
    0,                            /* tp_compare */
    0,                            /* tp_repr */
    0,                            /* tp_as_number */
    0,                            /* tp_as_sequence */
    0,                            /* tp_as_mapping */
    0,                            /* tp_hash */
    0,                            /* tp_call */
    0,                            /* tp_str */
    0,                            /* tp_getattro */
    0,                            /* tp_setattro */
    0,                            /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,      /* tp_flags */
    "Sander python wrapping",
                                  /* tp_doc */
    0,                            /* tp_traverse */
    0,                            /* tp_clear */
    0,                            /* tp_richcompare */
    0,                            /* tp_weaklistoffset */
    0,                            /* tp_iter */
    0,                            /* tp_iternext */
    sander_methods,               /* tp_methods */
    sander_members,               /* tp_members */
    0,                            /* tp_getset */
    0,                            /* tp_base */
    0,                            /* tp_dict */
    0,                            /* tp_descr_get */
    0,                            /* tp_descr_set */
    0,                            /* tp_dictoffset */
    (initproc)sander__init,       /* tp_init */
    0,                            /* tp_alloc */
    sander__new,                  /* tp_new */
};

void init_sander( void ) {
	PyObject    *my_module;

	my_module = Py_InitModule( "_sander", methods );
	PyType_Ready( &tSander );
    Py_INCREF( &tSander );
    PyModule_AddObject( my_module, "sander", (PyObject *) &tSander );
}

#endif

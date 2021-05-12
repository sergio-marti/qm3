#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <math.h>


typedef struct {
    PyObject_HEAD
    long	natm, ni;
    long	*qm, *mm;
    double	*sc;
    double	*epsi, *rmin;
} oQMLJ;


static int QMLJ__init( oQMLJ *self, PyObject *args, PyObject *kwds ) {
    PyObject	*o_qm = NULL, *o_mm = NULL, *o_ex = NULL, *tmp = NULL;
    PyObject	*o_mol = NULL, *o_epsi = NULL, *o_rmin = NULL;
    long		i, j, k, c, n_qm, n_mm, n_ex, *e_qm = NULL, *e_mm = NULL, i_qm, j_mm;
    double		*e_sc = NULL;

    if( PyArg_ParseTuple( args, "OOOO", &o_mol, &o_qm, &o_mm, &o_ex ) ) {

    	if( PyList_Check( o_qm ) && PyList_Check( o_mm ) && PyList_Check( o_ex ) ) {

    		o_epsi = PyObject_GetAttrString( o_mol, "epsi" );
    		o_rmin = PyObject_GetAttrString( o_mol, "rmin" );
    		self->natm = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
    		self->epsi = (double*) malloc( self->natm * sizeof( double ) );
    		self->rmin = (double*) malloc( self->natm * sizeof( double ) );
    		for( i = 0; i < self->natm; i++ ) {
    			self->epsi[i] = PyFloat_AsDouble( PyList_GetItem( o_epsi, i ) );
    			self->rmin[i] = PyFloat_AsDouble( PyList_GetItem( o_rmin, i ) );
    		}
    		Py_DECREF( o_epsi );
    		Py_DECREF( o_rmin );

    		n_qm = PyList_Size( o_qm );
    		n_mm = PyList_Size( o_mm );
    		n_ex = PyList_Size( o_ex );

    		if( n_ex > 0 ) {
    			e_qm = (long*) malloc( n_ex * sizeof( long ) );
    			e_mm = (long*) malloc( n_ex * sizeof( long ) );
    			e_sc = (double*) malloc( n_ex * sizeof( long ) );
    			for( i = 0; i < n_ex; i++ ) {
    				tmp = PyList_GetItem( o_ex, i );
    				e_qm[i] = PyLong_AsLong( PyList_GetItem( tmp, 0 ) );
    				e_mm[i] = PyLong_AsLong( PyList_GetItem( tmp, 1 ) );
    				e_sc[i] = PyFloat_AsDouble( PyList_GetItem( tmp, 2 ) );
    			}
    		}

    		self->ni = n_qm * n_mm;
    		self->qm = (long*) malloc( self->ni * sizeof( long ) );
    		self->mm = (long*) malloc( self->ni * sizeof( long ) );
    		self->sc = (double*) malloc( self->ni * sizeof( double ) );
    		c = 0;
    		for( i = 0; i < n_qm; i++ ) {
    			i_qm = PyLong_AsLong( PyList_GetItem( o_qm, i ) );
    			for( j = 0; j < n_mm; j++ ) {
    				j_mm = PyLong_AsLong( PyList_GetItem( o_mm, j ) );
    				self->qm[c] = i_qm;
    				self->mm[c] = j_mm;
    				self->sc[c] = 1.0;
    				for( k = 0; k < n_ex; k++ )
    					if( i_qm == e_qm[k] && j_mm == e_mm[k] )
    						self->sc[c] = e_sc[k];
    				c++;
    			}
    		}

    		free( e_qm ); free( e_mm ); free( e_sc );
    	}
    }
    return( 0 );
}


static PyObject* QMLJ__new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
    oQMLJ	*self;

    self = (oQMLJ*) type->tp_alloc( type, 0 );
    self->ni   = 0;
    self->qm   = NULL;
    self->mm   = NULL;
    self->sc   = NULL;
    self->natm = 0;
    self->epsi = NULL;
    self->rmin = NULL;
    return( (PyObject*) self ) ;
}


static void QMLJ__dealloc( oQMLJ *self ) {
    free( self->qm ); free( self->mm ); free( self->sc ); free( self->epsi ); free( self->rmin );
    self->ni   = 0;
    self->qm   = NULL;
    self->mm   = NULL;
    self->sc   = NULL;
    self->natm = 0;
    self->epsi = NULL;
    self->rmin = NULL;
    Py_TYPE( self )->tp_free( (PyObject*) self );
}


static PyObject* QMLJ__get_func( PyObject *self, PyObject *args ) {
    PyObject	*o_mol, *o_coor, *o_ener, *o_boxl;
    long		i, k;
    double		*coor = NULL, *ener = NULL;
    double		boxl[3], dr[3], r2, ss;
    oQMLJ		*obj = NULL;

    obj = (oQMLJ*) self;
    if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
    	o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
    	for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
    	Py_DECREF( o_boxl );
    	o_coor = PyObject_GetAttrString( o_mol, "coor" );
    	coor = (double*) malloc( 3 * obj->natm * sizeof( double ) );
    	ener = (double*) malloc( obj->natm * sizeof( double ) );
    	for( i = 0; i < 3 * obj->natm; i++ )
    		coor[i] = PyFloat_AsDouble( PyList_GetItem( o_coor, i ) );
    	for( i = 0; i < obj->natm; i++ )
			ener[i] = 0.0;
    	Py_DECREF( o_coor );

    	for( i = 0; i < obj->ni; i++ ) {
    		r2 = 0.0;
    		for( k = 0; k < 3; k++ ) {
    			dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
//    			if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
//    			if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
    			dr[k] -= boxl[k] * round( dr[k] / boxl[k] );
    			r2 += dr[k] * dr[k];
    		}
    		ss = ( obj->rmin[obj->qm[i]] + obj->rmin[obj->mm[i]] ) / sqrt( r2 );
    		ss = ss * ss * ss * ss * ss * ss;
    		ener[obj->mm[i]] += obj->epsi[obj->qm[i]] * obj->epsi[obj->mm[i]] * ss * ( ss - 2.0 ) * obj->sc[i];
    	}

    	o_ener = PyObject_GetAttrString( o_mol, "evdw" );
		for( i = 0; i < obj->natm; i++ ) PyList_SetItem( o_ener, i, PyFloat_FromDouble( ener[i] ) );
    	Py_DECREF( o_ener );

    	free( coor );
		free( ener );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static struct PyMethodDef QMLJ_methods [] = {
    { "get_func", (PyCFunction)QMLJ__get_func, METH_VARARGS },
    { 0, 0, 0 }
};


static struct PyMemberDef QMLJ_members [] = {
    { 0, 0, 0, 0 }
};


static struct PyMethodDef methods [] = {
    { 0, 0, 0 }
};


static PyTypeObject TQMLJ = {
    PyVarObject_HEAD_INIT( NULL, 0 )
    .tp_name = "QMLJ",
    .tp_doc = "Truncated Non-Bonded (QM:Lennard-Jones)",
    .tp_basicsize = sizeof( oQMLJ ),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = QMLJ__new,
    .tp_init = (initproc) QMLJ__init,
    .tp_dealloc = (destructor) QMLJ__dealloc,
    .tp_members = QMLJ_members,
    .tp_methods = QMLJ_methods,
};

static struct PyModuleDef moddef = {
    PyModuleDef_HEAD_INIT,
    "_mmdec",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__mmdec( void ) {
    PyObject    *my_module;

    my_module = PyModule_Create( &moddef );
    PyType_Ready( &TQMLJ );
    Py_INCREF( &TQMLJ );
    PyModule_AddObject( my_module, "QMLJ", (PyObject *) &TQMLJ );
    return( my_module );
}

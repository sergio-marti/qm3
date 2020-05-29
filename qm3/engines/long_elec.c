#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <math.h>

/*
técnicamente, si quiero calcular la componente long_elec sobre los átomos QM, 
bastaría con evaluar la electrostática de coulomb con el resto de átomos que
no he incluido en el cálculo QM... la aproximación que se hace es clásica,
ya que estoy teniendo en cuenta las cargas derivadas de la psi_polarizada para
calcular la interacción con lo que me queda de sistema...

Si el sistema es muy grande, se podría hacer uso de sumas de ewald:
1) hay que eliminar la componente de interacción QM/MM a partir de las cargas clásicas: n_qm * n_mm
2) incorporar la interacción directa (escalada por el erfc correspondiente): n_qm * n_mm
3) incorporar la interacción recíproca: n_kvec * ( N + n_qm )
4) incorporar el tinfoil (si se considera necesario, Martin no lo hace por defecto...): N

todo esto aplicaría sobre el vector gradiente de los átomos QM
si se modifica la energía, hay que calcular además la componente de las cargas al cuadrado...
*/
#define    _EC_	1389.35484620709144110151



static PyObject* coulomb( PyObject *self, PyObject *args ) {
    PyObject	*o_mol, *o_coor, *o_chrg, *o_grad, *o_boxl, *o_sqm, *o_smm;
    long		i, j, i3, I3, jj, j3, k, natm, nqm, nmm, *sqm;
    double		*coor = NULL, *chrg = NULL, *grd = NULL;
    double		boxl[3], dr[3], r2, df, tmp;

    if( PyArg_ParseTuple( args, "OOO", &o_mol, &o_sqm, &o_smm ) ) {

    	natm = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
    	o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
    	for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
    	Py_DECREF( o_boxl );
    	o_coor = PyObject_GetAttrString( o_mol, "coor" );
    	o_chrg = PyObject_GetAttrString( o_mol, "chrg" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	chrg = (double*) malloc( natm * sizeof( double ) );
    	for( i = 0; i < natm; i++ ) {
    		for( j = 0; j < 3; j++ )
    			coor[3*i+j] = PyFloat_AsDouble( PyList_GetItem( o_coor, 3*i+j ) );
    		chrg[i] = PyFloat_AsDouble( PyList_GetItem( o_chrg, i ) );
    	}
    	Py_DECREF( o_coor );
    	Py_DECREF( o_chrg );

    	nqm = (long) PyList_Size( o_sqm );
    	sqm = (long*) malloc( nqm * sizeof( long ) );
    	grd = (double*) malloc( 3 * nqm * sizeof( double ) );
    	for( i = 0; i < nqm; i++ ) {
    		sqm[i] = PyLong_AsLong( PyList_GetItem( o_sqm, i ) );
    		I3 = i * 3;
    		grd[I3] = 0.0; grd[I3+1] = 0.0; grd[I3+2] = 0.0;
    	}
    	nmm = (long) PyList_Size( o_smm );
    	for( j = 0; j < nmm; j++ ) {
    		jj = PyLong_AsLong( PyList_GetItem( o_smm, j ) );
    		j3 = jj * 3; 
    		for( i = 0; i < nqm; i++ ) {
    			I3 = i * 3;
    			i3 = sqm[i] * 3;
    			for( k = 0; k < 3; k++ ) {
    				dr[k] = coor[i3+k] - coor[j3+k];
    				dr[k] -= boxl[k] * round( dr[k] / boxl[k] );
    			}
    			r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    			df = - chrg[jj] * chrg[sqm[i]] / ( sqrt( r2 ) * r2 );
    			for( k = 0; k < 3; k++ ) grd[I3+k] += df * dr[k];
    		}
    	}

    	o_grad = PyObject_GetAttrString( o_mol, "grad" );
    	for( i = 0; i < nqm; i++ ) {
    		I3 = i * 3;
    		i3 = sqm[i] * 3;
    		for( k = 0; k < 3; k++ ) {
    			tmp = PyFloat_AsDouble( PyList_GetItem( o_grad, i3+k ) );
    			PyList_SetItem( o_grad, i3+k, PyFloat_FromDouble( tmp + grd[I3+k] * _EC_ ) );
    		}
    	}
    	Py_DECREF( o_grad );

    	free( coor ); free( chrg ); free( grd ); free( sqm );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* coulomb_remain( PyObject *self, PyObject *args ) {
    PyObject	*o_mol, *o_coor, *o_chrg, *o_grad, *o_boxl, *o_sqm, *o_smm;
    long		i, j, i3, I3, j3, k, natm, nqm, nmm, *sqm, *smm;
    double		*coor = NULL, *chrg = NULL, *grd = NULL;
    double		boxl[3], dr[3], r2, df, tmp;

    if( PyArg_ParseTuple( args, "OOO", &o_mol, &o_sqm, &o_smm ) ) {

    	natm = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
    	o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
    	for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
    	Py_DECREF( o_boxl );
    	o_coor = PyObject_GetAttrString( o_mol, "coor" );
    	o_chrg = PyObject_GetAttrString( o_mol, "chrg" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	chrg = (double*) malloc( natm * sizeof( double ) );
    	for( i = 0; i < natm; i++ ) {
    		for( j = 0; j < 3; j++ )
    			coor[3*i+j] = PyFloat_AsDouble( PyList_GetItem( o_coor, 3*i+j ) );
    		chrg[i] = PyFloat_AsDouble( PyList_GetItem( o_chrg, i ) );
    	}
    	Py_DECREF( o_coor );
    	Py_DECREF( o_chrg );

    	nqm = (long) PyList_Size( o_sqm );
    	sqm = (long*) malloc(  nqm * sizeof( long ) );
    	smm = (long*) malloc( natm * sizeof( long ) );
    	grd = (double*) malloc( 3 * nqm * sizeof( double ) );
    	for( i = 0; i < natm; i++ ) smm[i] = 1;
    	for( i = 0; i < nqm; i++ ) {
    		sqm[i] = PyLong_AsLong( PyList_GetItem( o_sqm, i ) );
    		I3 = i * 3;
    		grd[I3] = 0.0; grd[I3+1] = 0.0; grd[I3+2] = 0.0;
    		smm[sqm[i]] = 0;
    	}
    	nmm = (long) PyList_Size( o_smm );
    	for( j = 0; j < nmm; j++ ) smm[(long) PyLong_AsLong( PyList_GetItem( o_smm, j ) )] = 0;
    	for( j = 0; j < natm; j++ ) {
    		if( smm[j] == 1 ) {	
    			j3 = j * 3; 
    			for( i = 0; i < nqm; i++ ) {
    				I3 = i * 3;
    				i3 = sqm[i] * 3;
    				for( k = 0; k < 3; k++ ) {
    					dr[k] = coor[i3+k] - coor[j3+k];
    					dr[k] -= boxl[k] * round( dr[k] / boxl[k] );
    				}
    				r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    				df = - chrg[j] * chrg[sqm[i]] / ( sqrt( r2 ) * r2 );
    				for( k = 0; k < 3; k++ ) grd[I3+k] += df * dr[k];
    			}
    		}
    	}

    	o_grad = PyObject_GetAttrString( o_mol, "grad" );
    	for( i = 0; i < nqm; i++ ) {
    		I3 = i * 3;
    		i3 = sqm[i] * 3;
    		for( k = 0; k < 3; k++ ) {
    			tmp = PyFloat_AsDouble( PyList_GetItem( o_grad, i3+k ) );
    			PyList_SetItem( o_grad, i3+k, PyFloat_FromDouble( tmp + grd[I3+k] * _EC_ ) );
    		}
    	}
    	Py_DECREF( o_grad );

    	free( coor ); free( chrg ); free( grd ); free( sqm ); free( smm );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


// Taken from AmberTools MDGX: SpecialMath.c
double fperfc( double x ) {
    double		i, c, p, q;

    // Most likely, x will be between 0.5 and 4.0
    if( x > 0.5 && x < 4.0 ) {
    	p = ( ( ( ( ( ( ( -0.136864857382717e-6 ) * x + 0.564195517478974 ) * x + 7.21175825088309 ) * x + 43.1622272220567 ) * x +
    		152.989285046940 ) * x + 339.320816734344 ) * x + 451.918953711873 ) * x + 300.459261020162;
    	q = ( ( ( ( ( ( x + 12.7827273196294 ) * x + 77.0001529352295 ) * x + 277.585444743988 ) * x + 638.980264465631 ) * x +
    		931.354094850610 ) * x + 790.950925327898 ) * x + 300.459260956983;
    	return( exp( -x * x ) * p / q );
    }
    // In a minority of cases, x will be less than 0.5
    if( x <= 0.5 ) {
    	c = x * x;
    	p = ( ( ( -0.356098437018154e-1 ) * c + 6.99638348861914 ) * c + 21.9792616182942 ) * c + 242.667955230532;
    	q = ( ( c + 15.0827976304078 ) * c + 91.1649054045149 ) * c + 215.058875869861;
    	return( 1.0 - x * p / q );
    }
    // Otherwise, we have one more condition to evaluate about x
    if( x >= 4.0 && x < 5.65 ) {
    	i = 1.0 / x;
    	c = i * i;
    	p =( ( ( 0.0223192459734185 * c + 0.278661308609648 ) * c + 0.226956593539687 ) * c + 0.0494730910623251 ) * c + 0.00299610707703542;
    	q =( ( ( c + 1.98733201817135 ) * c + 1.05167510706793 ) * c + 0.191308926107830 ) * c + 0.0106209230528468;
    	c = ( - c * p / q + 0.564189583547756 ) * i;
    	return( exp( - x * x ) * c );
    }
    // If x is bigger than 5.65, the answer is easy
    return( 0.0 );
}



typedef struct {
    PyObject_HEAD
    long	nkvec;
    double	beta;
    double	*rkvec, *rkexp;
} oEwald;



static int ewald__init( oEwald *self, PyObject *args, PyObject *kwds ) {
    PyObject	*o_bxl, *o_kmx;
    double		boxl[3], tmp[3], r2, fac;
    long		kmax[3];
    long		i, j, k, l, n;
    int			f;
    
    if( PyArg_ParseTuple( args, "OOd", &o_bxl, &o_kmx, &(self->beta) ) ) {
    	for( i = 0; i < 3; i++ ) {
    		boxl[i] = PyFloat_AsDouble( PyList_GetItem( o_bxl, i ) );
    		kmax[i] = PyLong_AsLong( PyList_GetItem( o_kmx, i ) );
    	}
    	self->nkvec = ( ( 2 * kmax[0] + 1 ) * ( 2 * kmax[1] + 1 ) * ( 2 * kmax[2] + 1 ) - 1 ) / 2;
    	self->rkvec = (double*) malloc( 3 * self->nkvec * sizeof( double ) );
    	self->rkexp = (double*) malloc(     self->nkvec * sizeof( double ) );
    	n = 0;
    	for( k = -kmax[2]; k <= kmax[2]; k++ ) {
    		for( j = -kmax[1]; j <= kmax[1]; j++ ) {
    			for( i = -kmax[0]; i <= kmax[0]; i++ ) {
    				if( k == 0 && j == 0 && i == 0 ) continue;
    				tmp[0] = (double)i;
    				tmp[1] = (double)j;
    				tmp[2] = (double)k;
    				f = 0; l = 0;
    				while( l < n && f == 0 ) {
    					f = ( tmp[0] == - self->rkvec[3*l] && tmp[1] == - self->rkvec[3*l+1] && tmp[2] == - self->rkvec[3*l+2] );
    					l++;
    				}
    				if( f == 0 ) {
    					for( l = 0; l < 3; l++ ) self->rkvec[3*n+l] = tmp[l];
    					n++;
    				}
    			}
    		}
    	}
    	fac = _EC_ * 4.0 * M_PI / ( boxl[0] * boxl[1] * boxl[2] );
    	for( i = 0; i < self->nkvec; i++ ) {
    		r2 = 0.0;
    		for( j = 0; j < 3; j++ ) {
    			self->rkvec[3*i+j] *= 2.0 * M_PI / boxl[j];
    			r2 += self->rkvec[3*i+j] * self->rkvec[3*i+j];
    		}
    		self->rkexp[i] = fac * exp( - r2 / ( 4.0 * self->beta * self->beta ) ) / r2;
    	}
    }
    return( 0 );
}


static PyObject* ewald__new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
    oEwald	*self;

    self = (oEwald*) type->tp_alloc( type, 0 );
    self->beta  = 0.0;
    self->nkvec = 0;
    self->rkvec = NULL;
    self->rkexp = NULL;
    return( (PyObject*) self ) ;
}


static void ewald__dealloc( oEwald *self ) {
    free( self->rkvec );
    free( self->rkexp );
    Py_TYPE( self )->tp_free( (PyObject*) self );
}


static PyObject* ewald__calc( PyObject *self, PyObject *args ) {
    PyObject	*o_mol, *o_coor, *o_chrg, *o_grad, *o_boxl, *o_sqm, *o_smm;
    long		i, j, i3, I3, jj, j3, k, natm, nqm, nmm, *sqm;
    double		*coor = NULL, *chrg = NULL, *cos_ = NULL, *sin_ = NULL, *grd = NULL;
    double		boxl[3], dr[3], r2, r1, df, sqpi = sqrt( M_PI );
    oEwald		*obj = NULL;

    obj = (oEwald*) self;
    if( PyArg_ParseTuple( args, "OOO", &o_mol, &o_sqm, &o_smm ) ) {
    	natm = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
    	o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
    	for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
    	Py_DECREF( o_boxl );
    	o_coor = PyObject_GetAttrString( o_mol, "coor" );
    	o_chrg = PyObject_GetAttrString( o_mol, "chrg" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	chrg = (double*) malloc( natm * sizeof( double ) );
    	for( i = 0; i < natm; i++ ) {
    		i3 = i * 3;
    		for( j = 0; j < 3; j++ )
    			coor[i3+j] = PyFloat_AsDouble( PyList_GetItem( o_coor, i3+j ) );
    		chrg[i] = PyFloat_AsDouble( PyList_GetItem( o_chrg, i ) );
    	}
    	Py_DECREF( o_coor );
    	Py_DECREF( o_chrg );

    	// direct
    	nqm = (long) PyList_Size( o_sqm );
    	sqm = (long*) malloc( nqm * sizeof( long ) );
    	grd = (double*) malloc( 3 * nqm * sizeof( double ) );
    	for( i = 0; i < nqm; i++ ) {
    		sqm[i] = PyLong_AsLong( PyList_GetItem( o_sqm, i ) );
    		I3 = i * 3;
    		grd[I3] = 0.0; grd[I3+1] = 0.0; grd[I3+2] = 0.0;
    	}
    	nmm = (long) PyList_Size( o_smm );
    	for( j = 0; j < nmm; j++ ) {
    		jj = PyLong_AsLong( PyList_GetItem( o_smm, j ) );
    		j3 = jj * 3; 
    		for( i = 0; i < nqm; i++ ) {
    			I3 = i * 3;
    			i3 = sqm[i] * 3;
    			for( k = 0; k < 3; k++ ) {
    				dr[k] = coor[i3+k] - coor[j3+k];
    				dr[k] -= boxl[k] * round( dr[k] / boxl[k] );
    			}
    			r2  = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    			r1  = sqrt( r2 );
    			// df  = _EC_ * chrg[jj] * chrg[sqm[i]] * ( erfc( obj->beta * r1 ) / r1  + 2.0 * obj->beta * exp( - obj->beta * obj->beta * r2 ) / sqpi ) / r2;
    			// unscaled direct gradient removed...
    			df  = _EC_ * chrg[jj] * chrg[sqm[i]] * ( ( erfc( obj->beta * r1 ) - 1.0 ) / r1  + 2.0 * obj->beta * exp( - obj->beta * obj->beta * r2 ) / sqpi ) / r2;
    			for( k = 0; k < 3; k++ ) grd[I3+k] -= df * dr[k];
    		}
    	}

    	// reciprocal
    	cos_ = (double*) malloc( obj->nkvec * sizeof( double ) );
    	sin_ = (double*) malloc( obj->nkvec * sizeof( double ) );
    	for( i = 0; i < obj->nkvec; i++ ) {
    		i3 = i * 3;
    		cos_[i] = 0.0; sin_[i] = 0.0;
    		for( j = 0; j < natm; j++ ) {
    			j3 = 3 * j;
    			r2 = obj->rkvec[i3] * coor[j3] + obj->rkvec[i3+1] * coor[j3+1] + obj->rkvec[i3+2] * coor[j3+2];
    			cos_[i] += chrg[j] * cos( r2 );
    			sin_[i] += chrg[j] * sin( r2 );
    		}
    	}
    	for( i = 0; i < nqm; i++ ) {
    		I3  = i * 3;
    		i3  = sqm[i] * 3;
    		dr[0] = 0.0; dr[1] = 0.0; dr[2] = 0.0;
    		for( j = 0; j < obj->nkvec; j++ ) {
    			j3 = j * 3;
    			r2 = obj->rkvec[j3] * coor[i3] + obj->rkvec[j3+1] * coor[i3+1] + obj->rkvec[j3+2] * coor[i3+2];
    			df = obj->rkexp[j] * ( - cos_[j] * sin( r2 ) + sin_[j] * cos( r2 ) );
    			for( k = 0; k < 3; k++ ) dr[k] += df * obj->rkvec[j3+k];
    		}
    		for( k = 0; k < 3; k++ ) grd[I3+k] += 2.0 * chrg[sqm[i]] * dr[k];
    	}
    	free( cos_ ); free( sin_ );

    	o_grad = PyObject_GetAttrString( o_mol, "grad" );
    	for( i = 0; i < nqm; i++ ) {
    		I3 = i * 3;
    		i3 = sqm[i] * 3;
    		for( k = 0; k < 3; k++ ) {
    			df = PyFloat_AsDouble( PyList_GetItem( o_grad, i3+k ) );
    			PyList_SetItem( o_grad, i3+k, PyFloat_FromDouble( df + grd[I3+k] ) );
    		}
    	}
    	Py_DECREF( o_grad );

    	free( coor ); free( chrg ); free( grd );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static struct PyMethodDef ewald_methods [] = {
    { "calc", (PyCFunction)ewald__calc, METH_VARARGS },
    { 0, 0, 0 }
};


static struct PyMemberDef ewald_members [] = {
    { 0, 0, 0, 0 }
};


static struct PyMethodDef methods [] = {
    { "coulomb", (PyCFunction)coulomb, METH_VARARGS },
    { "coulomb_remain", (PyCFunction)coulomb_remain, METH_VARARGS },
    { 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3

static PyTypeObject Tewald = {
    PyVarObject_HEAD_INIT( NULL, 0 )
    .tp_name = "long_elec",
    .tp_doc = "Long Electrostatic corrections for QM/MM calculations",
    .tp_basicsize = sizeof( oEwald ),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = ewald__new,
    .tp_init = (initproc) ewald__init,
    .tp_dealloc = (destructor) ewald__dealloc,
    .tp_members = ewald_members,
    .tp_methods = ewald_methods,
};

static struct PyModuleDef moddef = {
    PyModuleDef_HEAD_INIT,
    "_long_elec",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__colvar_v( void ) {
    PyObject    *my_module;

    my_module = PyModule_Create( &moddef );
    PyType_Ready( &Tewald );
    Py_INCREF( &Tewald );
    PyModule_AddObject( my_module, "ewald", (PyObject *) &Tewald );
    return( my_module );
}

#else

static PyTypeObject Tewald = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "ewald",                    /* tp_name */
    sizeof( oEwald ),              /* tp_basicsize */
    0,                            /* tp_itemsize */
    (destructor)ewald__dealloc, /* tp_dealloc */
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
    "Long Electrostatic corrections for QM/MM calculations",
                                  /* tp_doc */
    0,                            /* tp_traverse */
    0,                            /* tp_clear */
    0,                            /* tp_richcompare */
    0,                            /* tp_weaklistoffset */
    0,                            /* tp_iter */
    0,                            /* tp_iternext */
    ewald_methods,              /* tp_methods */
    ewald_members,              /* tp_members */
    0,                            /* tp_getset */
    0,                            /* tp_base */
    0,                            /* tp_dict */
    0,                            /* tp_descr_get */
    0,                            /* tp_descr_set */
    0,                            /* tp_dictoffset */
    (initproc)ewald__init,      /* tp_init */
    0,                            /* tp_alloc */
    ewald__new,                 /* tp_new */
};

void init_long_elec( void ) {
    PyObject    *my_module;

    my_module = Py_InitModule( "_long_elec", methods );
    PyType_Ready( &Tewald );
    Py_INCREF( &Tewald );
    PyModule_AddObject( my_module, "ewald", (PyObject *) &Tewald );
}

#endif

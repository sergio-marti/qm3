#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <math.h>


typedef struct {
    PyObject_HEAD
    long	ni, nn;
    long	*qm, *mm;
    double	*sc, *lc;
    double	kmb, ref;
    double	con, cof;	// only in force_switch...
} oVInt;


// --------------------------------------------------------------------------------------


static int coulomb__init( oVInt *self, PyObject *args, PyObject *kwds ) {
    PyObject	*o_qm = NULL, *o_mm = NULL, *o_ex = NULL, *tmp = NULL;
    PyObject	*o_lc = NULL;
    long		i, j, k, c, n_mm, n_ex, *e_qm = NULL, *e_mm = NULL, i_qm, j_mm;
    double		*e_sc = NULL;

    if( PyArg_ParseTuple( args, "ddOOOO", &(self->kmb), &(self->ref), &o_qm, &o_mm, &o_lc, &o_ex ) ) {

    	if( PyList_Check( o_qm ) && PyList_Check( o_mm ) && PyList_Check( o_lc ) && PyList_Check( o_ex ) ) {

    		self->nn = PyList_Size( o_qm );
    		n_mm     = PyList_Size( o_mm );
    		self->ni = self->nn * n_mm;
    		n_ex     = PyList_Size( o_ex );

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
    		
    		self->qm = (long*) malloc( self->ni * sizeof( long ) );
    		self->mm = (long*) malloc( self->ni * sizeof( long ) );
    		self->sc = (double*) malloc( self->ni * sizeof( double ) );
    		self->lc = (double*) malloc( self->nn * sizeof( double ) );
    		c = 0;
    		for( i = 0; i < self->nn; i++ ) {
    			i_qm = PyLong_AsLong( PyList_GetItem( o_qm, i ) );
    			self->lc[i] = PyFloat_AsDouble( PyList_GetItem( o_lc, i ) );
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


static PyObject* coulomb__new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
    oVInt	*self;

    self = (oVInt*) type->tp_alloc( type, 0 );
    self->nn  = 0;
    self->ni  = 0;
    self->kmb = 0.0;
    self->ref = 0.0;
    self->con = 0.0;
    self->cof = 0.0;
    self->qm  = NULL;
    self->mm  = NULL;
    self->lc  = NULL;
    self->sc  = NULL;
    return( (PyObject*) self ) ;
}


static void coulomb__dealloc( oVInt *self ) {
    free( self->qm ); free( self->mm ); free( self->sc ); free( self->lc );
    self->nn  = 0;
    self->ni  = 0;
    self->kmb = 0.0;
    self->ref = 0.0;
    self->con = 0.0;
    self->cof = 0.0;
    self->qm  = NULL;
    self->mm  = NULL;
    self->lc  = NULL;
    self->sc  = NULL;
    Py_TYPE( self )->tp_free( (PyObject*) self );
}


static PyObject* coulomb__get_func( PyObject *self, PyObject *args ) {
    PyObject	*o_mol, *o_coor, *o_chrg, *o_boxl, *o_tmp;
    long		i, k, natm, lst_i, cur_i;
    double		*coor = NULL, *chrg = NULL;
    double		boxl[3], dr[3], r2, *pot = NULL, vv, tmp;
    double		EC = 1389.35484620709144110151;
    oVInt		*obj = NULL;

    obj = (oVInt*) self;
    if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
    	natm = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
    	o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
    	for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
    	Py_DECREF( o_boxl );
    	o_coor = PyObject_GetAttrString( o_mol, "coor" );
    	o_chrg = PyObject_GetAttrString( o_mol, "chrg" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	chrg = (double*) malloc( natm * sizeof( double ) );
    	for( i = 0; i < natm; i++ ) {
    		for( k = 0; k < 3; k++ )
    			coor[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_coor, 3*i+k ) );
    		chrg[i] = PyFloat_AsDouble( PyList_GetItem( o_chrg, i ) );
    	}
    	Py_DECREF( o_coor );
    	Py_DECREF( o_chrg );

    	pot = (double*) malloc( obj->nn * sizeof( double ) );
    	for( i = 0; i < obj->nn; i++ ) pot[i] = 0.0;

    	lst_i = obj->qm[0];
    	cur_i = 0;
    	for( i = 0; i < obj->ni; i++ ) {
    		if( lst_i != obj->qm[i] ) { cur_i++; lst_i = obj->qm[i]; }

    		r2 = 0.0;
    		for( k = 0; k < 3; k++ ) {
    			dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
    			if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
    			if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
    			r2 += dr[k] * dr[k];
    		}
    		pot[cur_i] += chrg[obj->mm[i]] / sqrt( r2 ) * obj->sc[i];
    	}
    	vv = 0.0; for( i = 0; i < obj->nn; i++ ) vv += pot[i] * obj->lc[i]; vv *= EC;

    	o_tmp = PyObject_GetAttrString( o_mol, "func" );
    	tmp = PyFloat_AsDouble( o_tmp );
    	Py_DECREF( o_tmp );
    	PyObject_SetAttrString( o_mol, "func", PyFloat_FromDouble( tmp + 0.5 * obj->kmb * ( vv - obj->ref ) * ( vv - obj->ref ) ) );

    	free( pot ); free( coor ); free( chrg );
    	
    	return( Py_BuildValue( "d", vv ) );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* coulomb__get_grad( PyObject *self, PyObject *args ) {
    PyObject	*o_mol, *o_coor, *o_chrg, *o_grad, *o_boxl, *o_tmp;
    long		i, k, natm, lst_i, cur_i;
    double		*coor = NULL, *chrg = NULL, *grad = NULL;
    double		boxl[3], dr[3], r2, df, *pot = NULL, vv, tmp;
    double		EC = 1389.35484620709144110151;
    oVInt		*obj = NULL;

    obj = (oVInt*) self;
    if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
    	natm = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
    	o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
    	for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
    	Py_DECREF( o_boxl );
    	o_coor = PyObject_GetAttrString( o_mol, "coor" );
    	o_grad = PyObject_GetAttrString( o_mol, "grad" );
    	o_chrg = PyObject_GetAttrString( o_mol, "chrg" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	grad = (double*) malloc( 3 * natm * sizeof( double ) );
    	chrg = (double*) malloc( natm * sizeof( double ) );
    	for( i = 0; i < natm; i++ ) {
    		for( k = 0; k < 3; k++ ) {
    			coor[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_coor, 3*i+k ) );
    			grad[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_grad, 3*i+k ) );
    		}
    		chrg[i] = PyFloat_AsDouble( PyList_GetItem( o_chrg, i ) );
    	}
    	Py_DECREF( o_coor );
    	Py_DECREF( o_chrg );

    	pot = (double*) malloc( obj->nn * sizeof( double ) );
    	for( i = 0; i < obj->nn; i++ ) pot[i] = 0.0;

    	lst_i = obj->qm[0];
    	cur_i = 0;
    	for( i = 0; i < obj->ni; i++ ) {
    		if( lst_i != obj->qm[i] ) { cur_i++; lst_i = obj->qm[i]; }

    		r2 = 0.0;
    		for( k = 0; k < 3; k++ ) {
    			dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
    			if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
    			if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
    			r2 += dr[k] * dr[k];
    		}
    		pot[cur_i] += chrg[obj->mm[i]] / sqrt( r2 ) * obj->sc[i];
    	}
    	vv = 0.0; for( i = 0; i < obj->nn; i++ ) vv += pot[i] * obj->lc[i]; vv *= EC;

    	o_tmp = PyObject_GetAttrString( o_mol, "func" );
    	tmp = PyFloat_AsDouble( o_tmp );
    	Py_DECREF( o_tmp );
    	PyObject_SetAttrString( o_mol, "func", PyFloat_FromDouble( tmp + 0.5 * obj->kmb * ( vv - obj->ref ) * ( vv - obj->ref ) ) );

    	tmp   = obj->kmb * ( vv - obj->ref );
    	lst_i = obj->qm[0];
    	cur_i = 0;
    	for( i = 0; i < obj->ni; i++ ) {
    		if( lst_i != obj->qm[i] ) { lst_i++; cur_i = obj->qm[i]; }

    		r2 = 0.0;
    		for( k = 0; k < 3; k++ ) {
    			dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
    			if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
    			if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
    			r2 += dr[k] * dr[k];
    		}
    		df = tmp * obj->lc[cur_i] * chrg[obj->mm[i]] / ( r2 * sqrt( r2 ) ) * obj->sc[i];
    		for( k = 0; k < 3; k++ ) {
    			grad[3*obj->qm[i]+k] -= df * dr[k];
    			grad[3*obj->mm[i]+k] += df * dr[k];
    		}
    	}

    	for( i = 0; i < 3 * natm; i++ ) PyList_SetItem( o_grad, i, PyFloat_FromDouble( grad[i] ) );
    	Py_DECREF( o_grad );
    	free( pot ); free( coor ); free( chrg ); free( grad );
    	
    	return( Py_BuildValue( "d", vv ) );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static struct PyMethodDef coulomb_methods [] = {
    { "get_func", (PyCFunction)coulomb__get_func, METH_VARARGS },
    { "get_grad", (PyCFunction)coulomb__get_grad, METH_VARARGS },
    { 0, 0, 0 }
};


static struct PyMemberDef coulomb_members [] = {
    { 0, 0, 0, 0 }
};


// --------------------------------------------------------------------------------------


static int fswitch__init( oVInt *self, PyObject *args, PyObject *kwds ) {
    PyObject	*o_qm = NULL, *o_mm = NULL, *o_ex = NULL, *tmp = NULL;
    PyObject	*o_lc = NULL;
    long		i, j, k, c, n_mm, n_ex, *e_qm = NULL, *e_mm = NULL, i_qm, j_mm;
    double		*e_sc = NULL;

// cut_on cut_off
    if( PyArg_ParseTuple( args, "ddddOOOO", &(self->kmb), &(self->ref),
    				&(self->con), &(self->cof),
    				&o_qm, &o_mm, &o_lc, &o_ex ) ) {

    	if( PyList_Check( o_qm ) && PyList_Check( o_mm ) && PyList_Check( o_lc ) && PyList_Check( o_ex ) ) {

    		self->nn = PyList_Size( o_qm );
    		n_mm     = PyList_Size( o_mm );
    		self->ni = self->nn * n_mm;
    		n_ex     = PyList_Size( o_ex );

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
    		
    		self->qm = (long*) malloc( self->ni * sizeof( long ) );
    		self->mm = (long*) malloc( self->ni * sizeof( long ) );
    		self->sc = (double*) malloc( self->ni * sizeof( double ) );
    		self->lc = (double*) malloc( self->nn * sizeof( double ) );
    		c = 0;
    		for( i = 0; i < self->nn; i++ ) {
    			i_qm = PyLong_AsLong( PyList_GetItem( o_qm, i ) );
    			self->lc[i] = PyFloat_AsDouble( PyList_GetItem( o_lc, i ) );
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


static PyObject* fswitch__new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
    oVInt	*self;

    self = (oVInt*) type->tp_alloc( type, 0 );
    self->nn  = 0;
    self->ni  = 0;
    self->kmb = 0.0;
    self->ref = 0.0;
    self->con = 0.0;
    self->cof = 0.0;
    self->qm  = NULL;
    self->mm  = NULL;
    self->lc  = NULL;
    self->sc  = NULL;
    return( (PyObject*) self ) ;
}


static void fswitch__dealloc( oVInt *self ) {
    free( self->qm ); free( self->mm ); free( self->sc ); free( self->lc );
    self->nn  = 0;
    self->ni  = 0;
    self->kmb = 0.0;
    self->ref = 0.0;
    self->con = 0.0;
    self->cof = 0.0;
    self->qm  = NULL;
    self->mm  = NULL;
    self->lc  = NULL;
    self->sc  = NULL;
    Py_TYPE( self )->tp_free( (PyObject*) self );
}


static PyObject* fswitch__get_func( PyObject *self, PyObject *args ) {
    PyObject	*o_mol, *o_coor, *o_chrg, *o_boxl, *o_tmp;
    long		i, k, natm, lst_i, cur_i;
    double		*coor = NULL, *chrg = NULL;
    double		boxl[3], dr[3], r2, *pot = NULL, vv, tmp;
    double		EC = 1389.35484620709144110151;
    oVInt		*obj = NULL;
    double		c2on, c2of, _g, _a, _b, _c, _d, _el1, _el2, r, r3, r5;

    obj = (oVInt*) self;
    if( PyArg_ParseTuple( args, "O", &o_mol ) ) {

    	c2on = obj->con * obj->con;
    	c2of = obj->cof * obj->cof;
    	_g   = pow( c2of - c2on, 3.0 );
    	_a   = c2of * c2of * ( c2of - 3.0 * c2on ) / _g;
    	_b   = 6.0 * c2of * c2on / _g;
    	_c   = - ( c2of + c2on ) / _g;
    	_d   = 0.4 / _g;
    	_el1 = 8.0 * ( c2of * c2on * ( obj->cof - obj->con ) - 0.2 * ( obj->cof * c2of * c2of - obj->con * c2on * c2on ) ) / _g;
    	_el2 = - _a / obj->cof + _b * obj->cof + _c * obj->cof * c2of + _d * obj->cof * c2of * c2of;

    	natm = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
    	o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
    	for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
    	Py_DECREF( o_boxl );
    	o_coor = PyObject_GetAttrString( o_mol, "coor" );
    	o_chrg = PyObject_GetAttrString( o_mol, "chrg" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	chrg = (double*) malloc( natm * sizeof( double ) );
    	for( i = 0; i < natm; i++ ) {
    		for( k = 0; k < 3; k++ )
    			coor[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_coor, 3*i+k ) );
    		chrg[i] = PyFloat_AsDouble( PyList_GetItem( o_chrg, i ) );
    	}
    	Py_DECREF( o_coor );
    	Py_DECREF( o_chrg );

    	pot = (double*) malloc( obj->nn * sizeof( double ) );
    	for( i = 0; i < obj->nn; i++ ) pot[i] = 0.0;

    	lst_i = obj->qm[0];
    	cur_i = 0;
    	for( i = 0; i < obj->ni; i++ ) {
    		if( lst_i != obj->qm[i] ) { cur_i++; lst_i = obj->qm[i]; }
//====================================
    		r2 = 0.0;
    		for( k = 0; k < 3; k++ ) {
    			dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
    			if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
    			if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
    			r2 += dr[k] * dr[k];
    		}
    		if( r2 > c2of ) { continue; }
    		r = sqrt( r2 );
    		if( r2 <= c2on ) {
    			pot[cur_i] += chrg[obj->mm[i]] * ( 1.0 / r + _el1 ) * obj->sc[i];
    		} else {
    			r3 =  r * r2;
    			r5 = r3 * r2;
    			pot[cur_i] += chrg[obj->mm[i]] * ( _a / r - _b * r - _c * r3 - _d * r5 + _el2 ) * obj->sc[i];
    		}
//====================================
    	}
    	vv = 0.0; for( i = 0; i < obj->nn; i++ ) vv += pot[i] * obj->lc[i]; vv *= EC;

    	o_tmp = PyObject_GetAttrString( o_mol, "func" );
    	tmp = PyFloat_AsDouble( o_tmp );
    	Py_DECREF( o_tmp );
    	PyObject_SetAttrString( o_mol, "func", PyFloat_FromDouble( tmp + 0.5 * obj->kmb * ( vv - obj->ref ) * ( vv - obj->ref ) ) );

    	free( pot ); free( coor ); free( chrg );
    	
    	return( Py_BuildValue( "d", vv ) );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* fswitch__get_grad( PyObject *self, PyObject *args ) {
    PyObject	*o_mol, *o_coor, *o_chrg, *o_grad, *o_boxl, *o_tmp;
    long		i, k, natm, lst_i, cur_i;
    double		*coor = NULL, *chrg = NULL, *grad = NULL;
    double		boxl[3], dr[3], r2, df, *pot = NULL, vv, tmp;
    double		EC = 1389.35484620709144110151;
    oVInt		*obj = NULL;
    double		c2on, c2of, _g, _a, _b, _c, _d, _el1, _el2, r, r3, r5;

    obj = (oVInt*) self;
    if( PyArg_ParseTuple( args, "O", &o_mol ) ) {

    	c2on = obj->con * obj->con;
    	c2of = obj->cof * obj->cof;
    	_g   = pow( c2of - c2on, 3.0 );
    	_a   = c2of * c2of * ( c2of - 3.0 * c2on ) / _g;
    	_b   = 6.0 * c2of * c2on / _g;
    	_c   = - ( c2of + c2on ) / _g;
    	_d   = 0.4 / _g;
    	_el1 = 8.0 * ( c2of * c2on * ( obj->cof - obj->con ) - 0.2 * ( obj->cof * c2of * c2of - obj->con * c2on * c2on ) ) / _g;
    	_el2 = - _a / obj->cof + _b * obj->cof + _c * obj->cof * c2of + _d * obj->cof * c2of * c2of;

    	natm = PyLong_AsLong( PyObject_GetAttrString( o_mol, "natm" ) );
    	o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
    	for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
    	Py_DECREF( o_boxl );
    	o_coor = PyObject_GetAttrString( o_mol, "coor" );
    	o_grad = PyObject_GetAttrString( o_mol, "grad" );
    	o_chrg = PyObject_GetAttrString( o_mol, "chrg" );
    	coor = (double*) malloc( 3 * natm * sizeof( double ) );
    	grad = (double*) malloc( 3 * natm * sizeof( double ) );
    	chrg = (double*) malloc( natm * sizeof( double ) );
    	for( i = 0; i < natm; i++ ) {
    		for( k = 0; k < 3; k++ ) {
    			coor[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_coor, 3*i+k ) );
    			grad[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_grad, 3*i+k ) );
    		}
    		chrg[i] = PyFloat_AsDouble( PyList_GetItem( o_chrg, i ) );
    	}
    	Py_DECREF( o_coor );
    	Py_DECREF( o_chrg );

    	pot = (double*) malloc( obj->nn * sizeof( double ) );
    	for( i = 0; i < obj->nn; i++ ) pot[i] = 0.0;

    	lst_i = obj->qm[0];
    	cur_i = 0;
    	for( i = 0; i < obj->ni; i++ ) {
    		if( lst_i != obj->qm[i] ) { cur_i++; lst_i = obj->qm[i]; }
//====================================
    		r2 = 0.0;
    		for( k = 0; k < 3; k++ ) {
    			dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
    			if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
    			if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
    			r2 += dr[k] * dr[k];
    		}
    		if( r2 > c2of ) { continue; }
    		r = sqrt( r2 );
    		if( r2 <= c2on ) {
    			pot[cur_i] += chrg[obj->mm[i]] * ( 1.0 / r + _el1 ) * obj->sc[i];
    		} else {
    			r3 =  r * r2;
    			r5 = r3 * r2;
    			pot[cur_i] += chrg[obj->mm[i]] * ( _a / r - _b * r - _c * r3 - _d * r5 + _el2 ) * obj->sc[i];
    		}
//====================================
    	}
    	vv = 0.0; for( i = 0; i < obj->nn; i++ ) vv += pot[i] * obj->lc[i]; vv *= EC;

    	o_tmp = PyObject_GetAttrString( o_mol, "func" );
    	tmp = PyFloat_AsDouble( o_tmp );
    	Py_DECREF( o_tmp );
    	PyObject_SetAttrString( o_mol, "func", PyFloat_FromDouble( tmp + 0.5 * obj->kmb * ( vv - obj->ref ) * ( vv - obj->ref ) ) );

    	tmp   = obj->kmb * ( vv - obj->ref );
    	lst_i = obj->qm[0];
    	cur_i = 0;
    	for( i = 0; i < obj->ni; i++ ) {
    		if( lst_i != obj->qm[i] ) { lst_i++; cur_i = obj->qm[i]; }
//====================================
    		r2 = 0.0;
    		for( k = 0; k < 3; k++ ) {
    			dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
    			if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
    			if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
    			r2 += dr[k] * dr[k];
    		}
    		if( r2 > c2of ) { continue; }
    		r = sqrt( r2 );
    		if( r2 <= c2on ) {
    			df = tmp * obj->lc[cur_i] * chrg[obj->mm[i]] / ( r2 * sqrt( r2 ) ) * obj->sc[i];
    		} else {
    			r3 =  r * r2;
    			df = tmp * obj->lc[cur_i] * chrg[obj->mm[i]] * ( _a / r3 + _b / r + 3.0 * _c * r + 5.0 * _d * r3 ) * obj->sc[i];
    		}
    		for( k = 0; k < 3; k++ ) {
    			grad[3*obj->qm[i]+k] -= df * dr[k];
    			grad[3*obj->mm[i]+k] += df * dr[k];
    		}
//====================================
    	}

    	for( i = 0; i < 3 * natm; i++ ) PyList_SetItem( o_grad, i, PyFloat_FromDouble( grad[i] ) );
    	Py_DECREF( o_grad );
    	free( pot ); free( coor ); free( chrg ); free( grad );
    	
    	return( Py_BuildValue( "d", vv ) );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static struct PyMethodDef fswitch_methods [] = {
    { "get_func", (PyCFunction)fswitch__get_func, METH_VARARGS },
    { "get_grad", (PyCFunction)fswitch__get_grad, METH_VARARGS },
    { 0, 0, 0 }
};


static struct PyMemberDef fswitch_members [] = {
    { 0, 0, 0, 0 }
};


// --------------------------------------------------------------------------------------


static struct PyMethodDef methods [] = {
    { 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3

static PyTypeObject Tcoulomb = {
    PyVarObject_HEAD_INIT( NULL, 0 )
    .tp_name = "coulomb",
    .tp_doc = "Electrostatic potential collective variable (Coulomb)",
    .tp_basicsize = sizeof( oVInt ),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = coulomb__new,
    .tp_init = (initproc) coulomb__init,
    .tp_dealloc = (destructor) coulomb__dealloc,
    .tp_members = coulomb_members,
    .tp_methods = coulomb_methods,
};

static PyTypeObject Tfswitch = {
    PyVarObject_HEAD_INIT( NULL, 0 )
    .tp_name = "fswitch",
    .tp_doc = "Electrostatic potential collective variable (force switch)",
    .tp_basicsize = sizeof( oVInt ),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = fswitch__new,
    .tp_init = (initproc) fswitch__init,
    .tp_dealloc = (destructor) fswitch__dealloc,
    .tp_members = fswitch_members,
    .tp_methods = fswitch_methods,
};

static struct PyModuleDef moddef = {
    PyModuleDef_HEAD_INIT,
    "_colvar_v",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__colvar_v( void ) {
    PyObject    *my_module;

    my_module = PyModule_Create( &moddef );
    PyType_Ready( &Tcoulomb );
    Py_INCREF( &Tcoulomb );
    PyModule_AddObject( my_module, "coulomb", (PyObject *) &Tcoulomb );
    PyType_Ready( &Tfswitch );
    Py_INCREF( &Tfswitch );
    PyModule_AddObject( my_module, "fswitch", (PyObject *) &Tfswitch );
    return( my_module );
}

#else

static PyTypeObject Tcoulomb = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "coulomb",                    /* tp_name */
    sizeof( oVInt ),              /* tp_basicsize */
    0,                            /* tp_itemsize */
    (destructor)coulomb__dealloc, /* tp_dealloc */
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
    "Electrostatic potential collective variable (Coulomb)",
                                  /* tp_doc */
    0,                            /* tp_traverse */
    0,                            /* tp_clear */
    0,                            /* tp_richcompare */
    0,                            /* tp_weaklistoffset */
    0,                            /* tp_iter */
    0,                            /* tp_iternext */
    coulomb_methods,              /* tp_methods */
    coulomb_members,              /* tp_members */
    0,                            /* tp_getset */
    0,                            /* tp_base */
    0,                            /* tp_dict */
    0,                            /* tp_descr_get */
    0,                            /* tp_descr_set */
    0,                            /* tp_dictoffset */
    (initproc)coulomb__init,      /* tp_init */
    0,                            /* tp_alloc */
    coulomb__new,                 /* tp_new */
};

static PyTypeObject Tfswitch = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "fswitch",                    /* tp_name */
    sizeof( oVInt ),              /* tp_basicsize */
    0,                            /* tp_itemsize */
    (destructor)fswitch__dealloc, /* tp_dealloc */
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
    "Electrostatic potential collective variable (force switch)",
                                  /* tp_doc */
    0,                            /* tp_traverse */
    0,                            /* tp_clear */
    0,                            /* tp_richcompare */
    0,                            /* tp_weaklistoffset */
    0,                            /* tp_iter */
    0,                            /* tp_iternext */
    fswitch_methods,              /* tp_methods */
    fswitch_members,              /* tp_members */
    0,                            /* tp_getset */
    0,                            /* tp_base */
    0,                            /* tp_dict */
    0,                            /* tp_descr_get */
    0,                            /* tp_descr_set */
    0,                            /* tp_dictoffset */
    (initproc)fswitch__init,      /* tp_init */
    0,                            /* tp_alloc */
    fswitch__new,                 /* tp_new */
};

void init_colvar_v( void ) {
    PyObject    *my_module;

    my_module = Py_InitModule( "_colvar_v", methods );
    PyType_Ready( &Tcoulomb );
    Py_INCREF( &Tcoulomb );
    PyModule_AddObject( my_module, "coulomb", (PyObject *) &Tcoulomb );
    PyType_Ready( &Tfswitch );
    Py_INCREF( &Tfswitch );
    PyModule_AddObject( my_module, "fswitch", (PyObject *) &Tfswitch );
}

#endif

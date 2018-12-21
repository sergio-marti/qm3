#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <math.h>


/* --------------------------------------------------------------------------------------
 *
 *	QM-LJ interactions (gradient + hessian)
 */


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


static PyObject* QMLJ__get_grad( PyObject *self, PyObject *args ) {
	PyObject	*o_mol, *o_coor, *o_grad, *o_boxl;
	long		i, k;
	double		*coor = NULL, *grad = NULL;
	double		boxl[3], dr[3], r2, s, ss, df;
	oQMLJ		*obj = NULL;

	obj = (oQMLJ*) self;
	if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
		o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
		for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
		Py_DECREF( o_boxl );
		o_coor = PyObject_GetAttrString( o_mol, "coor" );
		o_grad = PyObject_GetAttrString( o_mol, "grad" );
		coor = (double*) malloc( 3 * obj->natm * sizeof( double ) );
		grad = (double*) malloc( 3 * obj->natm * sizeof( double ) );
		for( i = 0; i < 3 * obj->natm; i++ ) {
			coor[i] = PyFloat_AsDouble( PyList_GetItem( o_coor, i ) );
			grad[i] = PyFloat_AsDouble( PyList_GetItem( o_grad, i ) );
		}
		Py_DECREF( o_coor );

		for( i = 0; i < obj->ni; i++ ) {
			r2 = 0.0;
			for( k = 0; k < 3; k++ ) {
				dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
//				if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
//				if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
				dr[k] -= boxl[k] * round( dr[k] / boxl[k] );
				r2 += dr[k] * dr[k];
			}
			s  = 1.0 / sqrt( r2 );
			ss = ( obj->rmin[obj->qm[i]] + obj->rmin[obj->mm[i]] ) * s;
			ss = ss * ss * ss * ss * ss * ss;
			df = 12.0 * obj->epsi[obj->qm[i]] * obj->epsi[obj->mm[i]] * ss * ( 1.0 - ss ) / r2 * obj->sc[i];
			for( k = 0; k < 3; k++ ) grad[3*obj->qm[i]+k] += df * dr[k];
		}

		for( i = 0; i < 3 * obj->natm; i++ ) PyList_SetItem( o_grad, i, PyFloat_FromDouble( grad[i] ) );
		Py_DECREF( o_grad );
		free( coor ); free( grad );
		
	}
	Py_INCREF( Py_None );
	return( Py_None );
}


static PyObject* QMLJ__get_hess( PyObject *self, PyObject *args ) {
	PyObject	*o_mol, *o_coor, *o_grad, *o_hess, *o_boxl;
	long		i, k, n3, nh, lst_i, cur_i;
	double		*coor = NULL, *grad = NULL, *hess = NULL;
	double		boxl[3], dr[3], r2, s, ss, tt, df, r4, d2f;
	oQMLJ		*obj = NULL;

	obj = (oQMLJ*) self;
	if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
		o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
		for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
		Py_DECREF( o_boxl );
		o_coor = PyObject_GetAttrString( o_mol, "coor" );
		o_grad = PyObject_GetAttrString( o_mol, "grad" );
		o_hess = PyObject_GetAttrString( o_mol, "hess" );
		coor = (double*) malloc( 3 * obj->natm * sizeof( double ) );
		grad = (double*) malloc( 3 * obj->natm * sizeof( double ) );
		nh   = PyList_Size( o_hess );
		n3   = roundl( sqrt( nh ) );
		hess = (double*) malloc( nh * sizeof( double ) );
		for( i = 0; i < 3 * obj->natm; i++ ) {
			coor[i] = PyFloat_AsDouble( PyList_GetItem( o_coor, i ) );
			grad[i] = PyFloat_AsDouble( PyList_GetItem( o_grad, i ) );
		}
		for( k = 0; k < nh; k++ ) hess[k] = PyFloat_AsDouble( PyList_GetItem( o_hess, k ) );
		Py_DECREF( o_coor );

		lst_i = obj->qm[0];
		cur_i = 0;
		for( i = 0; i < obj->ni; i++ ) {
			if( lst_i != obj->qm[i] ) { cur_i++; lst_i = obj->qm[i]; }

			r2 = 0.0;
			for( k = 0; k < 3; k++ ) {
				dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
//				if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
//				if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
				dr[k] -= boxl[k] * round( dr[k] / boxl[k] );
				r2 += dr[k] * dr[k];
			}
			s  = 1.0 / sqrt( r2 );
			ss = ( obj->rmin[obj->qm[i]] + obj->rmin[obj->mm[i]] ) * s;
			ss = ss * ss * ss * ss * ss * ss;
			tt = obj->epsi[obj->qm[i]] * obj->epsi[obj->mm[i]] * ss * obj->sc[i];
			df = 12.0 * tt * ( 1.0 - ss ) / r2;
			for( k = 0; k < 3; k++ ) grad[3*obj->qm[i]+k] += df * dr[k];

			r4         = r2 * r2;
			d2f        = 96.0 * tt * ( 1.750 * ss - 1.0 ) / r4;
			k          = ( cur_i * 3 + 0 ) * n3 + ( cur_i * 3 + 0 );
			hess[k]   += d2f * dr[0] * dr[0] + df;
			hess[k+1] += d2f * dr[0] * dr[1];
			hess[k+2] += d2f * dr[0] * dr[2];
			k          = ( cur_i * 3 + 1 ) * n3 + ( cur_i * 3 + 1 );
			hess[k]   += d2f * dr[1] * dr[1] + df;
			hess[k+1] += d2f * dr[1] * dr[2];
			k          = ( cur_i * 3 + 2 ) * n3 + ( cur_i * 3 + 2 );
			hess[k]   += d2f * dr[2] * dr[2] + df;
		}

		for( i = 0; i < 3 * obj->natm; i++ ) PyList_SetItem( o_grad, i, PyFloat_FromDouble( grad[i] ) );
		Py_DECREF( o_grad );
		for( i = 0; i < n3 - 1; i++ ) for( k = i + 1; k < n3; k++ ) hess[k*n3+i] = hess[i*n3+k];
		for( i = 0; i < nh; i++ ) PyList_SetItem( o_hess, i, PyFloat_FromDouble( hess[i] ) );
		Py_DECREF( o_hess );
		free( coor ); free( grad ); free( hess );
		
	}
	Py_INCREF( Py_None );
	return( Py_None );
}


static struct PyMethodDef QMLJ_methods [] = {
    { "get_grad", (PyCFunction)QMLJ__get_grad, METH_VARARGS },
    { "get_hess", (PyCFunction)QMLJ__get_hess, METH_VARARGS },
	{ 0, 0, 0 }
};


static struct PyMemberDef QMLJ_members [] = {
	{ 0, 0, 0, 0 }
};


/* --------------------------------------------------------------------------------------
 *
 *	QM-LJ interactions (gradient + hessian)
 *	MM-EL interactions (gradient)
 */


typedef struct {
	PyObject_HEAD
	long	ni, natm;
	long	*qm, *mm;
	double	*sc;
	double	*epsi, *rmin;
} oQMLJ_MMEL;


static int QMLJ_MMEL__init( oQMLJ_MMEL *self, PyObject *args, PyObject *kwds ) {
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


static PyObject* QMLJ_MMEL__new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
	oQMLJ_MMEL	*self;

	self = (oQMLJ_MMEL*) type->tp_alloc( type, 0 );
	self->ni   = 0;
	self->qm   = NULL;
	self->mm   = NULL;
	self->sc   = NULL;
	self->natm = 0;
	self->epsi = NULL;
	self->rmin = NULL;
	return( (PyObject*) self ) ;
}


static void QMLJ_MMEL__dealloc( oQMLJ_MMEL *self ) {
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


static PyObject* QMLJ_MMEL__get_grad( PyObject *self, PyObject *args ) {
	PyObject	*o_mol, *o_coor, *o_chrg, *o_grad, *o_boxl;
	long		i, k;
	double		*coor = NULL, *chrg = NULL, *grad = NULL;
	double		boxl[3], dr[3], r2, s, ss, df;
	double		EC = 1389.35484620709144110151;
	oQMLJ_MMEL	*obj = NULL;

	obj = (oQMLJ_MMEL*) self;
	if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
		o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
		for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
		Py_DECREF( o_boxl );
		o_coor = PyObject_GetAttrString( o_mol, "coor" );
		o_grad = PyObject_GetAttrString( o_mol, "grad" );
		o_chrg = PyObject_GetAttrString( o_mol, "chrg" );
		coor = (double*) malloc( 3 * obj->natm * sizeof( double ) );
		grad = (double*) malloc( 3 * obj->natm * sizeof( double ) );
		chrg = (double*) malloc( obj->natm * sizeof( double ) );
		for( i = 0; i < obj->natm; i++ ) {
			for( k = 0; k < 3; k++ ) {
				coor[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_coor, 3*i+k ) );
				grad[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_grad, 3*i+k ) );
			}
			chrg[i] = PyFloat_AsDouble( PyList_GetItem( o_chrg, i ) );
		}
		Py_DECREF( o_coor );
		Py_DECREF( o_chrg );

		for( i = 0; i < obj->ni; i++ ) {
			r2 = 0.0;
			for( k = 0; k < 3; k++ ) {
				dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
//				if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
//				if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
				dr[k] -= boxl[k] * round( dr[k] / boxl[k] );
				r2 += dr[k] * dr[k];
			}
			s  = 1.0 / sqrt( r2 );
			ss = ( obj->rmin[obj->qm[i]] + obj->rmin[obj->mm[i]] ) * s;
			ss = ss * ss * ss * ss * ss * ss;
			df = 12.0 * obj->epsi[obj->qm[i]] * obj->epsi[obj->mm[i]] * ss * ( 1.0 - ss ) / r2 * obj->sc[i];
			for( k = 0; k < 3; k++ ) grad[3*obj->qm[i]+k] += df * dr[k];

			df = EC * chrg[obj->qm[i]] * chrg[obj->mm[i]] / ( r2 * sqrt( r2 ) ) * obj->sc[i];
			for( k = 0; k < 3; k++ ) grad[3*obj->mm[i]+k] += df * dr[k];
		}

		for( i = 0; i < 3 * obj->natm; i++ ) PyList_SetItem( o_grad, i, PyFloat_FromDouble( grad[i] ) );
		Py_DECREF( o_grad );
		free( coor ); free( chrg ); free( grad );
		
	}
	Py_INCREF( Py_None );
	return( Py_None );
}


static PyObject* QMLJ_MMEL__get_hess( PyObject *self, PyObject *args ) {
	PyObject	*o_mol, *o_coor, *o_chrg, *o_grad, *o_hess, *o_boxl;
	long		i, k, n3, nh, lst_i, cur_i;
	double		*coor = NULL, *chrg = NULL, *grad = NULL, *hess = NULL;
	double		boxl[3], dr[3], r2, s, ss, tt, df, r4, d2f;
	double		EC = 1389.35484620709144110151;
	oQMLJ_MMEL	*obj = NULL;

	obj = (oQMLJ_MMEL*) self;
	if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
		o_boxl = PyObject_GetAttrString( o_mol, "boxl" );
		for( k = 0; k < 3; k++ ) boxl[k] = PyFloat_AsDouble( PyList_GetItem( o_boxl, k ) );
		Py_DECREF( o_boxl );
		o_coor = PyObject_GetAttrString( o_mol, "coor" );
		o_grad = PyObject_GetAttrString( o_mol, "grad" );
		o_hess = PyObject_GetAttrString( o_mol, "hess" );
		o_chrg = PyObject_GetAttrString( o_mol, "chrg" );
		coor = (double*) malloc( 3 * obj->natm * sizeof( double ) );
		grad = (double*) malloc( 3 * obj->natm * sizeof( double ) );
		chrg = (double*) malloc( obj->natm * sizeof( double ) );
		nh   = PyList_Size( o_hess );
		n3   = roundl( sqrt( nh ) );
		hess = (double*) malloc( nh * sizeof( double ) );
		for( i = 0; i < obj->natm; i++ ) {
			for( k = 0; k < 3; k++ ) {
				coor[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_coor, 3*i+k ) );
				grad[3*i+k] = PyFloat_AsDouble( PyList_GetItem( o_grad, 3*i+k ) );
			}
			chrg[i] = PyFloat_AsDouble( PyList_GetItem( o_chrg, i ) );
		}
		for( k = 0; k < nh; k++ ) hess[k] = PyFloat_AsDouble( PyList_GetItem( o_hess, k ) );
		Py_DECREF( o_coor );
		Py_DECREF( o_chrg );

		lst_i = obj->qm[0];
		cur_i = 0;
		for( i = 0; i < obj->ni; i++ ) {

			if( lst_i != obj->qm[i] ) { cur_i++; lst_i = obj->qm[i]; }

			r2 = 0.0;
			for( k = 0; k < 3; k++ ) {
				dr[k] = coor[3*obj->qm[i]+k] - coor[3*obj->mm[i]+k];
//				if( dr[k] >    boxl[k] * 0.5 ) { dr[k] -= boxl[k]; }
//				if( dr[k] <= - boxl[k] * 0.5 ) { dr[k] += boxl[k]; }
				dr[k] -= boxl[k] * round( dr[k] / boxl[k] );
				r2 += dr[k] * dr[k];
			}
			s  = 1.0 / sqrt( r2 );
			ss = ( obj->rmin[obj->qm[i]] + obj->rmin[obj->mm[i]] ) * s;
			ss = ss * ss * ss * ss * ss * ss;
			tt = obj->epsi[obj->qm[i]] * obj->epsi[obj->mm[i]] * ss * obj->sc[i];
			df = 12.0 * tt * ( 1.0 - ss ) / r2;
			for( k = 0; k < 3; k++ ) grad[3*obj->qm[i]+k] += df * dr[k];

			r4         = r2 * r2;
			d2f        = 96.0 * tt * ( 1.750 * ss - 1.0 ) / r4;
			k          = ( cur_i * 3 + 0 ) * n3 + ( cur_i * 3 + 0 );
			hess[k]   += d2f * dr[0] * dr[0] + df;
			hess[k+1] += d2f * dr[0] * dr[1];
			hess[k+2] += d2f * dr[0] * dr[2];
			k          = ( cur_i * 3 + 1 ) * n3 + ( cur_i * 3 + 1 );
			hess[k]   += d2f * dr[1] * dr[1] + df;
			hess[k+1] += d2f * dr[1] * dr[2];
			k          = ( cur_i * 3 + 2 ) * n3 + ( cur_i * 3 + 2 );
			hess[k]   += d2f * dr[2] * dr[2] + df;

			df = EC * chrg[obj->qm[i]] * chrg[obj->mm[i]] / ( r2 * sqrt( r2 ) ) * obj->sc[i];
			for( k = 0; k < 3; k++ ) grad[3*obj->mm[i]+k] += df * dr[k];
		}

		for( i = 0; i < 3 * obj->natm; i++ ) PyList_SetItem( o_grad, i, PyFloat_FromDouble( grad[i] ) );
		Py_DECREF( o_grad );
		for( i = 0; i < n3 - 1; i++ ) for( k = i + 1; k < n3; k++ ) hess[k*n3+i] = hess[i*n3+k];
		for( i = 0; i < nh; i++ ) PyList_SetItem( o_hess, i, PyFloat_FromDouble( hess[i] ) );
		Py_DECREF( o_hess );
		free( coor ); free( chrg ); free( grad ); free( hess );
		
	}
	Py_INCREF( Py_None );
	return( Py_None );
}


static struct PyMethodDef QMLJ_MMEL_methods [] = {
    { "get_grad", (PyCFunction)QMLJ_MMEL__get_grad, METH_VARARGS },
    { "get_hess", (PyCFunction)QMLJ_MMEL__get_hess, METH_VARARGS },
	{ 0, 0, 0 }
};


static struct PyMemberDef QMLJ_MMEL_members [] = {
	{ 0, 0, 0, 0 }
};


// --------------------------------------------------------------------------------------


static struct PyMethodDef methods [] = {
	{ 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3

static PyTypeObject TQMLJ = {
	PyVarObject_HEAD_INIT( NULL, 0 )
	.tp_name = "Int_QMLJ",
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

static PyTypeObject TQMLJ_MMEL = {
	PyVarObject_HEAD_INIT( NULL, 0 )
	.tp_name = "Int_QMLJ_MMEL",
	.tp_doc = "Truncated Non-Bonded (QM:Lennard-Jones + MM:Electrostatics)",
	.tp_basicsize = sizeof( oQMLJ_MMEL ),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	.tp_new = QMLJ_MMEL__new,
	.tp_init = (initproc) QMLJ_MMEL__init,
	.tp_dealloc = (destructor) QMLJ_MMEL__dealloc,
	.tp_members = QMLJ_MMEL_members,
	.tp_methods = QMLJ_MMEL_methods,
};


static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_qmmm",
	NULL,
	-1,
	methods
};

PyMODINIT_FUNC PyInit__qmmm( void ) {
	PyObject    *my_module;

	my_module = PyModule_Create( &moddef );
	PyType_Ready( &TQMLJ );
    Py_INCREF( &TQMLJ );
    PyModule_AddObject( my_module, "Int_QMLJ", (PyObject *) &TQMLJ );
	PyType_Ready( &TQMLJ_MMEL );
    Py_INCREF( &TQMLJ_MMEL );
    PyModule_AddObject( my_module, "Int_QMLJ_MMEL", (PyObject *) &TQMLJ_MMEL );
	return( my_module );
}

#else

static PyTypeObject TQMLJ = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Int_QMLJ",                /* tp_name */
    sizeof( oQMLJ ),           /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)QMLJ__dealloc, /* tp_dealloc */
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
    "Truncated Non-Bonded (QM:Leniard-Jones)",
                               /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    QMLJ_methods,              /* tp_methods */
    QMLJ_members,              /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)QMLJ__init,      /* tp_init */
    0,                         /* tp_alloc */
    QMLJ__new,                 /* tp_new */
};

static PyTypeObject TQMLJ_MMEL = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "Int_QMLJ_MMEL",                /* tp_name */
    sizeof( oQMLJ_MMEL ),           /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor)QMLJ_MMEL__dealloc, /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_compare */
    0,                              /* tp_repr */
    0,                              /* tp_as_number */
    0,                              /* tp_as_sequence */
    0,                              /* tp_as_mapping */
    0,                              /* tp_hash */
    0,                              /* tp_call */
    0,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,        /* tp_flags */
	"Truncated Non-Bonded (QM:Leniard-Jones + MM:Electrostatics)",
                                    /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    QMLJ_MMEL_methods,              /* tp_methods */
    QMLJ_MMEL_members,              /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    (initproc)QMLJ_MMEL__init,      /* tp_init */
    0,                              /* tp_alloc */
    QMLJ_MMEL__new,                 /* tp_new */
};

void init_qmmm( void ) {
	PyObject    *my_module;

	my_module = Py_InitModule( "_qmmm", methods );
	PyType_Ready( &TQMLJ );
    Py_INCREF( &TQMLJ );
    PyModule_AddObject( my_module, "Int_QMLJ", (PyObject *) &TQMLJ );
	PyType_Ready( &TQMLJ_MMEL );
    Py_INCREF( &TQMLJ_MMEL );
    PyModule_AddObject( my_module, "Int_QMLJ_MMEL", (PyObject *) &TQMLJ_MMEL );
}

#endif

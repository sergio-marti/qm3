#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))



void __ij( long w, long n, long *i, long *j ) {
	int		f = 0;
	double	ii, jj, nn, ww;

	*i = -1; *j = -1;
	nn = (double) n; ww = (double) w;
	if( ww < 0.0 || ww >= nn * ( nn - 1.0 ) / 2.0 ) { return; }
	jj = nn - 1.0;
	while( jj >= .0 && f == 0 ) {
		ii = ( 2.0 * nn - 3.0 - sqrt( 8.0 * jj - 8.0 * ww + 4.0 * nn * ( nn - 3.0 ) + 1.0 ) ) / 2.0;
		if( fmod( ii, 1.0 ) == 0.0 && ii < jj ) { f = 1; *i = (long) ii; *j = (long) jj; }
		else { jj -= 1.0; }
	}
}



typedef struct con_bnd_node { long i,j; struct con_bnd_node *n; } con_bnd;

typedef struct con_ang_node { long i,j,k; struct con_ang_node *n; } con_ang;
typedef struct { long siz, _i0, _if; long *lst; con_ang *ang; } ang_arg;

typedef struct con_dih_node { long i,j,k,l; struct con_dih_node *n; } con_dih;
typedef struct { long siz, _i0, _if; long *lst; con_dih *dih; } dih_arg;

typedef struct { long siz, _i0, _if; long n_bnd, n_ang, n_dih; long *qms, *bnd, *ang, *dih;
				double cut, *xyz, *box; con_bnd *nbn; } nbn_arg;

typedef struct { long who, _i0, _if, n_lst, *lst, n_dat, *ind; double *xyz, *grd, ene, *dat; } ene_arg;

typedef struct { long who, _i0, _if, n_lst, *lst, n_dat; double *box, *xyz, *grd, ele, vdw,
				*dat, *scl, con, cof, eps; } int_arg;



// ####################################################################################################################


void* __angles( void *args ) {
	ang_arg		*arg = (ang_arg*) args;
	long		w, i, j;
	con_ang		*p;

	p = arg->ang;
	for( w = arg->_i0; w < arg->_if; w++ ) {
		__ij( w, arg->siz, &i, &j );
		if( i == -1 || j == -1 ) { continue; }
		if( arg->lst[2*i] == arg->lst[2*j] ) {
			arg->ang->i++;
			p->n    = (con_ang*) malloc( sizeof( con_ang ) );
			p->n->i = arg->lst[2*i+1];
			p->n->j = arg->lst[2*i];
			p->n->k = arg->lst[2*j+1];
			p->n->n = NULL;
			p       = p->n;
		} else if( arg->lst[2*i] == arg->lst[2*j+1] ) {
			arg->ang->i++;
			p->n    = (con_ang*) malloc( sizeof( con_ang ) );
			p->n->i = arg->lst[2*i+1];
			p->n->j = arg->lst[2*i];
			p->n->k = arg->lst[2*j];
			p->n->n = NULL;
			p       = p->n;
		} else if( arg->lst[2*i+1] == arg->lst[2*j] ) {
			arg->ang->i++;
			p->n    = (con_ang*) malloc( sizeof( con_ang ) );
			p->n->i = arg->lst[2*i];
			p->n->j = arg->lst[2*i+1];
			p->n->k = arg->lst[2*j+1];
			p->n->n = NULL;
			p       = p->n;
		} else if( arg->lst[2*i+1] == arg->lst[2*j+1] ) {
			arg->ang->i++;
			p->n    = (con_ang*) malloc( sizeof( con_ang ) );
			p->n->i = arg->lst[2*i];
			p->n->j = arg->lst[2*i+1];
			p->n->k = arg->lst[2*j];
			p->n->n = NULL;
			p       = p->n;
		}
	}
	return( NULL );
}

static PyObject* w_guess_angles( PyObject *self, PyObject *args ) {
	PyObject	*out, *object, *otmp;
	long		i, j, cpu, *lst, siz, *rng, dsp, nit;
	pthread_t	*pid;
	ang_arg		*arg;
	con_ang		*ptr;

	if( PyArg_ParseTuple( args, "O", &object ) ) {
//		cpu  = PyInt_AsLong( PyObject_GetAttrString( object, "ncpu" ) );
		cpu  = PyLong_AsLong( PyObject_GetAttrString( object, "ncpu" ) );

		otmp = PyObject_GetAttrString( object, "bond" );
		siz  = PyList_Size( otmp );
		lst  = (long*) malloc( 2*siz * sizeof( long ) );
		for( i = 0; i < siz; i++ )
			for( j = 0; j < 2; j++ ) lst[2*i+j] = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), j ) );
		Py_DECREF( otmp );

		nit = siz * ( siz - 1 ) / 2;
		dsp = (long) ((float)nit / (float)cpu);
		rng = (long*) malloc( (cpu+1) * sizeof( long ) );
		if( dsp == 0 ) { for( i = 0; i < cpu+1; i++ ) rng[i] = 0;  for( i = 0; i < nit + 1; i++ ) rng[i] = i; }
		else { for( i = 0; i < cpu; i++ ) rng[i] = i * dsp; rng[cpu] = nit; }
		pid = (pthread_t*) malloc( cpu * sizeof( pthread_t ) );
		arg = (ang_arg*) malloc( cpu * sizeof( ang_arg ) );
		for( i = 0; i < cpu; i++ ) {
			arg[i].siz    = siz;
			arg[i]._i0    = rng[i];
			arg[i]._if    = rng[i+1];
			arg[i].lst    = lst;
			arg[i].ang    = (con_ang*) malloc( sizeof( con_ang ) ); 
			arg[i].ang->i = 0;
			arg[i].ang->j = 0;
			arg[i].ang->k = 0;
			arg[i].ang->n = NULL;
			pthread_create( &pid[i], NULL, __angles, (void*) &arg[i] );
		}
		for( i = 0; i < cpu; i++ ) pthread_join( pid[i], NULL );

// access to the connectivity array to check possible bonds between i & k 
long nat = PyLong_AsLong( PyObject_GetAttrString( object, "natm" ) );
long **con, *nel;
PyObject *olst;
otmp = PyObject_GetAttrString( object, "conn" );
con = (long**) malloc( nat * sizeof( long* ) );
nel = (long*) malloc( nat * sizeof( long ) );
for( i = 0 ; i < nat; i++ ) {
	olst = PyList_GetItem( otmp, i );
	nel[i] = PyList_Size( olst );
	con[i] = (long*) malloc( nel[i] * sizeof( long ) );
	for( j = 0; j < nel[i]; j++ ) con[i][j] = PyLong_AsLong( PyList_GetItem( olst, j ) );
}
Py_DECREF( otmp );
out = PyList_New( 0 );
for( i = 0; i < cpu; i++ ) {
	ptr = arg[i].ang->n;
	while( ptr != NULL ) {
			siz = 0; for( j = 0; j < nel[ptr->i]; j++ ) if( ptr->k == con[ptr->i][j] ) siz++;
			if( siz == 0 ) PyList_Append( out, Py_BuildValue( "[l,l,l]", ptr->i, ptr->j, ptr->k ) );
		ptr = ptr->n;
	}
}
free( nel ); for( i = 0; i < nat; i++ ) free( con[i] ); free( con );

		free( rng ); free( lst ); free( pid );
		for( i = 0; i < cpu; i++ ) {
			ptr = arg[i].ang;
			while( ptr != NULL ) {
				ptr = ptr->n;
				free( arg[i].ang );
				arg[i].ang = ptr;
			}
		}
		free( arg );
		return( out );
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}


// ####################################################################################################################


void* __dihedrals( void *args ) {
	dih_arg		*arg = (dih_arg*) args;
	long		w, i, j;
	con_dih		*p;

	p = arg->dih;
	for( w = arg->_i0; w < arg->_if; w++ ) {
		__ij( w, arg->siz, &i, &j );
		if( i == -1 || j == -1 ) { continue; }
		if( arg->lst[3*i+1] == arg->lst[3*j] && arg->lst[3*i+2] == arg->lst[3*j+1] ) {
			arg->dih->i++;
			p->n    = (con_dih*) malloc( sizeof( con_dih ) );
			p->n->i = arg->lst[3*i];
			p->n->j = arg->lst[3*i+1];
			p->n->k = arg->lst[3*i+2];
			p->n->l = arg->lst[3*j+2];
			p->n->n = NULL;
			p       = p->n;
		} else if( arg->lst[3*i+1] == arg->lst[3*j+2] && arg->lst[3*i+2] == arg->lst[3*j+1] ) {
			arg->dih->i++;
			p->n    = (con_dih*) malloc( sizeof( con_dih ) );
			p->n->i = arg->lst[3*i];
			p->n->j = arg->lst[3*i+1];
			p->n->k = arg->lst[3*i+2];
			p->n->l = arg->lst[3*j];
			p->n->n = NULL;
			p       = p->n;
		} else if( arg->lst[3*i+1] == arg->lst[3*j] && arg->lst[3*i] == arg->lst[3*j+1] ) {
			arg->dih->i++;
			p->n    = (con_dih*) malloc( sizeof( con_dih ) );
			p->n->i = arg->lst[3*i+2];
			p->n->j = arg->lst[3*i+1];
			p->n->k = arg->lst[3*i];
			p->n->l = arg->lst[3*j+2];
			p->n->n = NULL;
			p       = p->n;
		} else if( arg->lst[3*i+1] == arg->lst[3*j+2] && arg->lst[3*i] == arg->lst[3*j+1] ) {
			arg->dih->i++;
			p->n    = (con_dih*) malloc( sizeof( con_dih ) );
			p->n->i = arg->lst[3*i+2];
			p->n->j = arg->lst[3*i+1];
			p->n->k = arg->lst[3*i];
			p->n->l = arg->lst[3*j];
			p->n->n = NULL;
			p       = p->n;
		}
	}
	return( NULL );
}

static PyObject* w_guess_dihedrals( PyObject *self, PyObject *args ) {
	PyObject	*out, *object, *otmp;
	long		i, j, cpu, *lst, siz, *rng, dsp, nit;
	pthread_t	*pid;
	dih_arg		*arg;
	con_dih		*ptr;

	if( PyArg_ParseTuple( args, "O", &object ) ) {
//		cpu  = PyInt_AsLong( PyObject_GetAttrString( object, "ncpu" ) );
		cpu  = PyLong_AsLong( PyObject_GetAttrString( object, "ncpu" ) );

		otmp = PyObject_GetAttrString( object, "angl" );
		siz  = PyList_Size( otmp );
		lst  = (long*) malloc( 3*siz * sizeof( long ) );
		for( i = 0; i < siz; i++ )
			for( j = 0; j < 3; j++ ) lst[3*i+j] = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), j ) );
		Py_DECREF( otmp );

		nit = siz * ( siz - 1 ) / 2;
		dsp = (long) ((float)nit / (float)cpu);
		rng = (long*) malloc( (cpu+1) * sizeof( long ) );
		if( dsp == 0 ) { for( i = 0; i < cpu+1; i++ ) rng[i] = 0;  for( i = 0; i < nit + 1; i++ ) rng[i] = i; }
		else { for( i = 0; i < cpu; i++ ) rng[i] = i * dsp; rng[cpu] = nit; }
		pid = (pthread_t*) malloc( cpu * sizeof( pthread_t ) );
		arg = (dih_arg*) malloc( cpu * sizeof( dih_arg ) );
		for( i = 0; i < cpu; i++ ) {
			arg[i].siz    = siz;
			arg[i]._i0    = rng[i];
			arg[i]._if    = rng[i+1];
			arg[i].lst    = lst;
			arg[i].dih    = (con_dih*) malloc( sizeof( con_dih ) ); 
			arg[i].dih->i = 0;
			arg[i].dih->j = 0;
			arg[i].dih->k = 0;
			arg[i].dih->l = 0;
			arg[i].dih->n = NULL;
			pthread_create( &pid[i], NULL, __dihedrals, (void*) &arg[i] );
		}
		for( i = 0; i < cpu; i++ ) pthread_join( pid[i], NULL );

// access to the connectivity array to check possible bonds between i & k 
long nat = PyLong_AsLong( PyObject_GetAttrString( object, "natm" ) );
long **con, *nel;
PyObject *olst;
otmp = PyObject_GetAttrString( object, "conn" );
con = (long**) malloc( nat * sizeof( long* ) );
nel = (long*) malloc( nat * sizeof( long ) );
for( i = 0 ; i < nat; i++ ) {
	olst = PyList_GetItem( otmp, i );
	nel[i] = PyList_Size( olst );
	con[i] = (long*) malloc( nel[i] * sizeof( long ) );
	for( j = 0; j < nel[i]; j++ ) con[i][j] = PyLong_AsLong( PyList_GetItem( olst, j ) );
}
Py_DECREF( otmp );
out = PyList_New( 0 );
for( i = 0; i < cpu; i++ ) {
	ptr = arg[i].dih->n;
	while( ptr != NULL ) {
			siz = 0; for( j = 0; j < nel[ptr->i]; j++ ) if( ptr->l == con[ptr->i][j] ) siz++;
			if( siz == 0 ) PyList_Append( out, Py_BuildValue( "[l,l,l,l]", ptr->i, ptr->j, ptr->k, ptr->l ) );
		ptr = ptr->n;
	}
}
free( nel ); for( i = 0; i < nat; i++ ) free( con[i] ); free( con );

		free( rng ); free( lst ); free( pid );
		for( i = 0; i < cpu; i++ ) {
			ptr = arg[i].dih;
			while( ptr != NULL ) {
				ptr = ptr->n;
				free( arg[i].dih );
				arg[i].dih = ptr;
			}
		}
		free( arg );
		return( out );
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}


// ####################################################################################################################


void* __update_non_bonded( void *args ) {
	nbn_arg		*arg = (nbn_arg*) args;
	long		w, i, j, k, i3, j3, f;
	double		r2, dr; //i_xyz[3], j_xyz[3];
	con_bnd		*p;

	p = arg->nbn;
	for( w = arg->_i0; w < arg->_if; w++ ) {
		__ij( w, arg->siz, &i, &j );
		if( i == -1 || j == -1 ) { continue; }
		if( arg->qms[i] == 1 && arg->qms[j] == 1 ) { continue; }
		i3  = 3 * i;
		j3  = 3 * j;
/*
		for( k = 0; k < 3; k++ ) {
			i_xyz[k] = arg->xyz[i3+k] - arg->box[k] * round( arg->xyz[i3+k] / arg->box[k] );
			j_xyz[k] = arg->xyz[j3+k] - arg->box[k] * round( arg->xyz[j3+k] / arg->box[k] );
		}
		r2  = ( i_xyz[0] - j_xyz[0] ) * ( i_xyz[0] - j_xyz[0] ) +
				( i_xyz[1] - j_xyz[1] ) * ( i_xyz[1] - j_xyz[1] ) +
				( i_xyz[2] - j_xyz[2] ) * ( i_xyz[2] - j_xyz[2] );
*/
		r2 = 0.0;
		for( k = 0; k < 3; k++ ) {
			dr  = arg->xyz[i3+k] - arg->xyz[j3+k];
			dr -= arg->box[k] * round( dr / arg->box[k] );
			r2 += dr * dr;
		}
		if( r2 <= arg->cut ) {
			f = 0;
			k = 0;
			while( k < arg->n_bnd && f == 0 ) {
				f |= ( ( i == arg->bnd[2*k] && j == arg->bnd[2*k+1] ) || ( i == arg->bnd[2*k+1] && j == arg->bnd[2*k] ) );
				k++;
			}
			k = 0;
			while( k < arg->n_ang && f == 0 ) {
				f |= ( ( i == arg->ang[2*k] && j == arg->ang[2*k+1] ) || ( i == arg->ang[2*k+1] && j == arg->ang[2*k] ) );
				k++;
			}
			k = 0;
			while( k < arg->n_dih && f == 0 ) {
				f |= ( ( i == arg->dih[2*k] && j == arg->dih[2*k+1] ) || ( i == arg->dih[2*k+1] && j == arg->dih[2*k] ) );
				k++;
			}
			if( f == 0 ) {
				arg->nbn->i++;
				p->n    = (con_bnd*) malloc( sizeof( con_bnd ) );
				p->n->i = i;
				p->n->j = j;
				p->n->n = NULL;
				p       = p->n;
			}
		}
	}
	return( NULL );
}

static PyObject* w_update_non_bonded( PyObject *self, PyObject *args ) {
	PyObject	*out, *object, *molecule, *otmp;
	double		*xyz, cut, box[3];
	long		*bnd, *ang, *dih, n_bnd, n_ang, n_dih, i, j, n3, n, cpu, *lst, dsp, nit, *qms;
	pthread_t	*pid;
	nbn_arg		*arg;
	con_bnd		*ptr;

	if( PyArg_ParseTuple( args, "OO", &object, &molecule ) ) {
//		cpu  = PyInt_AsLong( PyObject_GetAttrString( object, "ncpu" ) );
		cpu  = PyLong_AsLong( PyObject_GetAttrString( object, "ncpu" ) );

		otmp = PyObject_GetAttrString( object, "cut_list" );
		cut = PyFloat_AsDouble( otmp );
		if( cut> 0.0 ) { cut *= cut; } else { cut = 1.0e99; }
		Py_DECREF( otmp );

		otmp = PyObject_GetAttrString( molecule, "boxl" );
		for( i = 0; i < 3; i++ ) box[i] = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );

		otmp = PyObject_GetAttrString( molecule, "coor" );
		n3   = PyList_Size( otmp );
		n    = n3 / 3;
		xyz  = (double*) malloc( n3 * sizeof( double ) );
		for( i = 0; i < n3; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );

		otmp = PyObject_GetAttrString( object, "qmat" );
		qms  = (long*) malloc( n * sizeof( long ) );
		for( i = 0; i < n; i++ ) qms[i] = ( Py_True == PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );

		otmp  = PyObject_GetAttrString( object, "bond" );
		n_bnd = PyList_Size( otmp );
		bnd   = (long*) malloc( 2*n_bnd * sizeof( long ) );
		for( i = 0; i < n_bnd; i++ ) {
			bnd[2*i]   = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), 0 ) );
			bnd[2*i+1] = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), 1 ) );
		}
		Py_DECREF( otmp );

		otmp  = PyObject_GetAttrString( object, "angl" );
		n_ang = PyList_Size( otmp );
		ang   = (long*) malloc( 2*n_ang * sizeof( long ) );
		for( i = 0; i < n_ang; i++ ) {
			ang[2*i]   = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), 0 ) );
			ang[2*i+1] = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), 2 ) );
		}
		Py_DECREF( otmp );

		otmp  = PyObject_GetAttrString( object, "dihe" );
		n_dih = PyList_Size( otmp );
		dih   = (long*) malloc( 2*n_dih * sizeof( long ) );
		for( i = 0; i < n_dih; i++ ) {
			dih[2*i]   = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), 0 ) );
			dih[2*i+1] = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), 3 ) );
		}
		Py_DECREF( otmp );

		nit = n * ( n - 1 ) / 2;
		dsp = (long) ((float)nit / (float)cpu);
		lst = (long*) malloc( (cpu+1) * sizeof( long ) );
		if( dsp == 0 ) { for( i = 0; i < cpu+1; i++ ) lst[i] = 0;  for( i = 0; i < nit + 1; i++ ) lst[i] = i; }
		else { for( i = 0; i < cpu; i++ ) lst[i] = i * dsp; lst[cpu] = nit; }
		pid = (pthread_t*) malloc( cpu * sizeof( pthread_t ) );
		arg = (nbn_arg*) malloc( cpu * sizeof( nbn_arg ) );
		for( i = 0; i < cpu; i++ ) {
			arg[i].siz    = n;
			arg[i]._i0    = lst[i];
			arg[i]._if    = lst[i+1];
			arg[i].cut    = cut;
			arg[i].xyz    = xyz;
			arg[i].box    = box;
			arg[i].qms    = qms;
			arg[i].bnd    = bnd;
			arg[i].ang    = ang;
			arg[i].dih    = dih;
			arg[i].n_bnd  = n_bnd;
			arg[i].n_ang  = n_ang;
			arg[i].n_dih  = n_dih;
			arg[i].nbn    = (con_bnd*) malloc( sizeof( con_bnd ) ); 
			arg[i].nbn->i = 0;
			arg[i].nbn->j = 0;
			arg[i].nbn->n = NULL;
			pthread_create( &pid[i], NULL, __update_non_bonded, (void*) &arg[i] );
		}
		for( i = 0; i < cpu; i++ ) pthread_join( pid[i], NULL );

		n = 0;
		for( i = 0; i < cpu; i++ ) { n += arg[i].nbn->i; }
		out = PyList_New( n );
		j = 0;
		for( i = 0; i < cpu; i++ ) {
			ptr = arg[i].nbn->n;
			while( ptr != NULL ) {
				PyList_SetItem( out, j++, Py_BuildValue( "[l,l,d]", ptr->i, ptr->j, 1.0 ) );
				ptr = ptr->n;
			}
		}

		free( qms ); free( bnd ); free( ang ); free( dih ); free( xyz ); free( lst ); free( pid );
		for( i = 0; i < cpu; i++ ) {
			ptr = arg[i].nbn;
			while( ptr != NULL ) {
				ptr = ptr->n;
				free( arg[i].nbn );
				arg[i].nbn = ptr;
			}
		}
		free( arg );
		return( out );
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}


// ####################################################################################################################


void* __energy_bond( void *args ) {
	ene_arg		*arg = (ene_arg*) args;
	long		i, j, ai, aj;
	double		vec[3], val, dif, tmp;

	for( i = arg->_i0; i < arg->_if; i++ ) {
		ai = 3 * arg->lst[2*i];
		aj = 3 * arg->lst[2*i+1];
		for( j = 0; j < 3; j++ ) vec[j] = arg->xyz[ai+j] - arg->xyz[aj+j];
		val = sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
		dif = val - arg->dat[arg->ind[i]*2+1];
		tmp = dif * arg->dat[arg->ind[i]*2];
		arg->ene += tmp * dif;
		if( arg->grd != NULL ) {
			tmp *= 2.0 / val;
			for( j = 0; j < 3; j++ ) {
				arg->grd[arg->who+ai+j] += tmp * vec[j];
				arg->grd[arg->who+aj+j] -= tmp * vec[j];
			}
		}
	}
	return( NULL );
}

static PyObject* w_energy_bond( PyObject *self, PyObject *args ) {
	PyObject	*gradient, *object, *molecule, *otmp;
	double		*xyz, *grd, tmp;
	long		i, j, n3, cpu;
	long		*lst, n_lst, n_dat, *ind;
	double		*dat, out = 0.0;
	long		*rng, dsp, nit;
	pthread_t	*pid;
	ene_arg		*arg;

	if( PyArg_ParseTuple( args, "OOO", &object, &molecule, &gradient ) ) {
//		cpu  = PyInt_AsLong( PyObject_GetAttrString( object, "ncpu" ) );
		cpu  = PyLong_AsLong( PyObject_GetAttrString( object, "ncpu" ) );

		otmp = PyObject_GetAttrString( molecule, "coor" );
		n3   = PyList_Size( otmp );
		xyz  = (double*) malloc( n3 * sizeof( double ) );
		for( i = 0; i < n3; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );

		if( gradient != Py_True ) { grd = NULL; }
		else { grd = (double*) malloc( cpu*n3 * sizeof( double ) ); for( i = 0; i < cpu*n3; i++ ) grd[i] = 0.0; }

		otmp  = PyObject_GetAttrString( object, "bond" );
		n_lst = PyList_Size( otmp );
		lst   = (long*) malloc( 2*n_lst * sizeof( long ) );
		for( i = 0; i < n_lst; i++ )
			for( j = 0; j < 2; j++ ) lst[2*i+j] = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), j ) );
		Py_DECREF( otmp );
		otmp  = PyObject_GetAttrString( object, "bond_data" );
		n_dat = PyList_Size( otmp );
		dat    = (double*) malloc( 2*n_dat * sizeof( double ) );
		for( i = 0; i < n_dat; i++ ) 
			for( j = 0; j < 2; j++ ) dat[2*i+j] = PyFloat_AsDouble( PyList_GetItem( PyList_GetItem( otmp, i ), j ) );
		Py_DECREF( otmp );
		otmp   = PyObject_GetAttrString( object, "bond_indx" );
		ind    = (long*) malloc( n_lst * sizeof( long ) );
		for( i = 0; i < n_lst; i++ ) ind[i] = PyLong_AsLong( PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );
		nit = n_lst;
		dsp = (long) ((float)nit / (float)cpu);
		rng = (long*) malloc( (cpu+1) * sizeof( long ) );
		if( dsp == 0 ) { for( i = 0; i < cpu+1; i++ ) rng[i] = 0;  for( i = 0; i < nit + 1; i++ ) rng[i] = i; }
		else { for( i = 0; i < cpu; i++ ) rng[i] = i * dsp; rng[cpu] = nit; }
		pid = (pthread_t*) malloc( cpu * sizeof( pthread_t ) );
		arg = (ene_arg*) malloc( cpu * sizeof( ene_arg ) );
		for( i = 0; i < cpu; i++ ) {
			arg[i].ene   = 0.0;
			arg[i].who   = i * n3; 
			arg[i]._i0   = rng[i];
			arg[i]._if   = rng[i+1];
			arg[i].xyz   = xyz;
			arg[i].grd   = grd;
			arg[i].lst   = lst;
			arg[i].n_lst = n_lst;
			arg[i].dat   = dat;
			arg[i].n_dat = n_dat;
			arg[i].ind   = ind;
			pthread_create( &pid[i], NULL, __energy_bond, (void*) &arg[i] );
		}
		for( i = 0; i < cpu; i++ ) pthread_join( pid[i], NULL );
		for( i = 0; i < cpu; i++ ) out += arg[i].ene;
		free( xyz ); free( rng ); free( pid ); free ( lst ); free( dat ); free( ind ); free( arg );

		if( grd != NULL ) {
			otmp = PyObject_GetAttrString( molecule, "grad" );
			for( i = 0; i < n3; i++ ) {
				tmp = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
				for( j = 0; j < cpu; j++ ) tmp += grd[i+j*n3];
				PyList_SetItem( otmp, i, PyFloat_FromDouble( tmp ) );
			}
			Py_DECREF( otmp );
		}

		free( grd );
		return( Py_BuildValue( "d", out ) );
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}


// ####################################################################################################################


void* __energy_angle( void *args ) {
	ene_arg		*arg = (ene_arg*) args;
	long		i, j, ai, aj, ak;
	double		dij[3], rij, dkj[3], rkj, val, dif, tmp, fac, dti[3], dtj[3], dtk[3];

	for( i = arg->_i0; i < arg->_if; i++ ) {
		ai = 3 * arg->lst[3*i];
		aj = 3 * arg->lst[3*i+1];
		ak = 3 * arg->lst[3*i+2];
		for( j = 0; j < 3; j++ ) dij[j] = arg->xyz[ai+j] - arg->xyz[aj+j];
		rij = sqrt( dij[0]*dij[0] + dij[1]*dij[1] + dij[2]*dij[2] );
		for( j = 0; j < 3; j++ ) dij[j] /= rij;
		for( j = 0; j < 3; j++ ) dkj[j] = arg->xyz[ak+j] - arg->xyz[aj+j];
		rkj = sqrt( dkj[0]*dkj[0] + dkj[1]*dkj[1] + dkj[2]*dkj[2] );
		for( j = 0; j < 3; j++ ) dkj[j] /= rkj;
		for( fac = 0.0, j = 0; j < 3; j++ ) fac += dij[j] * dkj[j];
		fac = min( fabs( fac ), 1.0 - 1.0e-6 ) * fac / fabs( fac );
		val = acos( fac );
		dif = val - arg->dat[arg->ind[i]*2+1];
		tmp = dif * arg->dat[arg->ind[i]*2];
		arg->ene += tmp * dif;
		if( arg->grd != NULL ) {
			tmp *= -2.0 / sqrt( 1.0 - fac * fac );
			for( j = 0; j < 3; j++ ) {
				dti[j] = ( dkj[j] - fac * dij[j] ) / rij;
				dtk[j] = ( dij[j] - fac * dkj[j] ) / rkj;
				dtj[j] = - ( dti[j] + dtk[j] );
			}
			for( j = 0; j < 3; j++ ) {
				arg->grd[arg->who+ai+j] += tmp * dti[j];
				arg->grd[arg->who+aj+j] += tmp * dtj[j];
				arg->grd[arg->who+ak+j] += tmp * dtk[j];
			}
		}
	}
	return( NULL );
}

static PyObject* w_energy_angle( PyObject *self, PyObject *args ) {
	PyObject	*gradient, *object, *molecule, *otmp;
	double		*xyz, *grd, tmp;
	long		i, j, n3, cpu;
	long		*lst, n_lst, n_dat, *ind;
	double		*dat, out = 0.0;
	long		*rng, dsp, nit;
	pthread_t	*pid;
	ene_arg		*arg;

	if( PyArg_ParseTuple( args, "OOO", &object, &molecule, &gradient ) ) {
//		cpu  = PyInt_AsLong( PyObject_GetAttrString( object, "ncpu" ) );
		cpu  = PyLong_AsLong( PyObject_GetAttrString( object, "ncpu" ) );

		otmp = PyObject_GetAttrString( molecule, "coor" );
		n3   = PyList_Size( otmp );
		xyz  = (double*) malloc( n3 * sizeof( double ) );
		for( i = 0; i < n3; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );

		if( gradient != Py_True ) { grd = NULL; }
		else { grd = (double*) malloc( cpu*n3 * sizeof( double ) ); for( i = 0; i < cpu*n3; i++ ) grd[i] = 0.0; }

		otmp  = PyObject_GetAttrString( object, "angl" );
		n_lst = PyList_Size( otmp );
		lst   = (long*) malloc( 3*n_lst * sizeof( long ) );
		for( i = 0; i < n_lst; i++ )
			for( j = 0; j < 3; j++ ) lst[3*i+j] = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), j ) );
		Py_DECREF( otmp );
		otmp  = PyObject_GetAttrString( object, "angl_data" );
		n_dat = PyList_Size( otmp );
		dat    = (double*) malloc( 2*n_dat * sizeof( double ) );
		for( i = 0; i < n_dat; i++ ) 
			for( j = 0; j < 2; j++ ) dat[2*i+j] = PyFloat_AsDouble( PyList_GetItem( PyList_GetItem( otmp, i ), j ) );
		Py_DECREF( otmp );
		otmp   = PyObject_GetAttrString( object, "angl_indx" );
		ind    = (long*) malloc( n_lst * sizeof( long ) );
		for( i = 0; i < n_lst; i++ ) ind[i] = PyLong_AsLong( PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );
		nit = n_lst;
		dsp = (long) ((float)nit / (float)cpu);
		rng = (long*) malloc( (cpu+1) * sizeof( long ) );
		if( dsp == 0 ) { for( i = 0; i < cpu+1; i++ ) rng[i] = 0;  for( i = 0; i < nit + 1; i++ ) rng[i] = i; }
		else { for( i = 0; i < cpu; i++ ) rng[i] = i * dsp; rng[cpu] = nit; }
		pid = (pthread_t*) malloc( cpu * sizeof( pthread_t ) );
		arg = (ene_arg*) malloc( cpu * sizeof( ene_arg ) );
		for( i = 0; i < cpu; i++ ) {
			arg[i].ene   = 0.0;
			arg[i].who   = i * n3;
			arg[i]._i0   = rng[i];
			arg[i]._if   = rng[i+1];
			arg[i].xyz   = xyz;
			arg[i].grd   = grd;
			arg[i].lst   = lst;
			arg[i].n_lst = n_lst;
			arg[i].dat   = dat;
			arg[i].n_dat = n_dat;
			arg[i].ind   = ind;
			pthread_create( &pid[i], NULL, __energy_angle, (void*) &arg[i] );
		}
		for( i = 0; i < cpu; i++ ) pthread_join( pid[i], NULL );
		for( i = 0; i < cpu; i++ ) out += arg[i].ene;
		free( xyz ); free( rng ); free( pid ); free ( lst ); free( dat ); free( ind ); free( arg );

		if( grd != NULL ) {
			otmp = PyObject_GetAttrString( molecule, "grad" );
			for( i = 0; i < n3; i++ ) {
				tmp = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
				for( j = 0; j < cpu; j++ ) tmp += grd[i+j*n3];
				PyList_SetItem( otmp, i, PyFloat_FromDouble( tmp ) );
			}
			Py_DECREF( otmp );
		}

		free( grd );
		return( Py_BuildValue( "d", out ) );
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}


// ####################################################################################################################


void* __energy_dihedral( void *args ) {
	ene_arg		*arg = (ene_arg*) args;
	long		i, j, ai, aj, ak, al;
	double		rkj, rt2, ru2, rtu, cd, sd, dph;
	double		dji[3], dkj[3], dlk[3], vt[3], vu[3], vtu[3], dki[3], dlj[3], dvt[3], dvu[3];
	double		cs1, cs2, cs3, cs4, cs5, cs6;
	double		sn1, sn2, sn3, sn4, sn5, sn6;

	for( i = arg->_i0; i < arg->_if; i++ ) {
		ai = 3 * arg->lst[4*i];
		aj = 3 * arg->lst[4*i+1];
		ak = 3 * arg->lst[4*i+2];
		al = 3 * arg->lst[4*i+3];
		for( j = 0; j < 3; j++ ) dji[j] = arg->xyz[aj+j] - arg->xyz[ai+j];
		for( j = 0; j < 3; j++ ) dkj[j] = arg->xyz[ak+j] - arg->xyz[aj+j];
		rkj = sqrt( dkj[0]*dkj[0] + dkj[1]*dkj[1] + dkj[2]*dkj[2] );
		for( j = 0; j < 3; j++ ) dlk[j] = arg->xyz[al+j] - arg->xyz[ak+j];
		vt[0] = dji[1] * dkj[2] - dkj[1] * dji[2];
		vt[1] = dji[2] * dkj[0] - dkj[2] * dji[0];
		vt[2] = dji[0] * dkj[1] - dkj[0] * dji[1];
		rt2 = vt[0]*vt[0] + vt[1]*vt[1] + vt[2]*vt[2];
		vu[0] = dkj[1] * dlk[2] - dlk[1] * dkj[2];
		vu[1] = dkj[2] * dlk[0] - dlk[2] * dkj[0];
		vu[2] = dkj[0] * dlk[1] - dlk[0] * dkj[1];
		ru2 = vu[0]*vu[0] + vu[1]*vu[1] + vu[2]*vu[2];
		vtu[0] = vt[1] * vu[2] - vu[1] * vt[2];
		vtu[1] = vt[2] * vu[0] - vu[2] * vt[0];
		vtu[2] = vt[0] * vu[1] - vu[0] * vt[1];
		rtu = sqrt( rt2 * ru2 );
		if( rtu == 0.0 ) { continue; }
		for( cs1 = 0.0, j = 0; j < 3; j++ ) cs1 += vt[j] * vu[j]; cs1 /= rtu;
		for( sn1 = 0.0, j = 0; j < 3; j++ ) sn1 += dkj[j] * vtu[j]; sn1 /= ( rtu * rkj );
		cs2 = cs1 * cs1 - sn1 * sn1;
		sn2 = 2.0 * cs1 * sn1;
		cs3 = cs1 * cs2 - sn1 * sn2;
		sn3 = cs1 * sn2 + sn1 * cs2;
		cs4 = cs1 * cs3 - sn1 * sn3;
		sn4 = cs1 * sn3 + sn1 * cs3;
		cs5 = cs1 * cs4 - sn1 * sn4;
		sn5 = cs1 * sn4 + sn1 * cs4;
		cs6 = cs1 * cs5 - sn1 * sn5;
		sn6 = cs1 * sn5 + sn1 * cs5;
		dph = 0.0;
		if( arg->dat[12*arg->ind[i]] != 0.0 ) { 
			cd        = cos( arg->dat[12*arg->ind[i]+1] );
			sd        = sin( arg->dat[12*arg->ind[i]+1] );
			dph      += arg->dat[12*arg->ind[i]] * ( cs1 * sd - sn1 * cd );
			arg->ene += arg->dat[12*arg->ind[i]] * ( 1.0 + cs1 * cd + sn1 * sd );
		}
		if( arg->dat[12*arg->ind[i]+2] != 0.0 ) { 
			cd        = cos( arg->dat[12*arg->ind[i]+3] );
			sd        = sin( arg->dat[12*arg->ind[i]+3] );
			dph      += arg->dat[12*arg->ind[i]+2] * 2.0 * ( cs2 * sd - sn2 * cd );
			arg->ene += arg->dat[12*arg->ind[i]+2] * ( 1.0 + cs2 * cd + sn2 * sd );
		}
		if( arg->dat[12*arg->ind[i]+4] != 0.0 ) { 
			cd        = cos( arg->dat[12*arg->ind[i]+5] );
			sd        = sin( arg->dat[12*arg->ind[i]+5] );
			dph      += arg->dat[12*arg->ind[i]+4] * 3.0 * ( cs3 * sd - sn3 * cd );
			arg->ene += arg->dat[12*arg->ind[i]+4] * ( 1.0 + cs3 * cd + sn3 * sd );
		}
		if( arg->dat[12*arg->ind[i]+6] != 0.0) { 
			cd        = cos( arg->dat[12*arg->ind[i]+7] );
			sd        = sin( arg->dat[12*arg->ind[i]+7] );
			dph      += arg->dat[12*arg->ind[i]+6] * 4.0 * ( cs4 * sd - sn4 * cd );
			arg->ene += arg->dat[12*arg->ind[i]+6] * ( 1.0 + cs4 * cd + sn4 * sd );
		}
		if( arg->dat[12*arg->ind[i]+8] != 0.0 ) { 
			cd        = cos( arg->dat[12*arg->ind[i]+9] );
			sd        = sin( arg->dat[12*arg->ind[i]+9] );
			dph      += arg->dat[12*arg->ind[i]+8] * 5.0 * ( cs5 * sd - sn5 * cd );
			arg->ene += arg->dat[12*arg->ind[i]+8] * ( 1.0 + cs5 * cd + sn5 * sd );
		}
		if( arg->dat[12*arg->ind[i]+10] != 0.0 ) { 
			cd        = cos( arg->dat[12*arg->ind[i]+11] );
			sd        = sin( arg->dat[12*arg->ind[i]+11] );
			dph      += arg->dat[12*arg->ind[i]+10] * 6.0 * ( cs6 * sd - sn6 * cd );
			arg->ene += arg->dat[12*arg->ind[i]+10] * ( 1.0 + cs6 * cd + sn6 * sd );
		}
		if( arg->grd != NULL ) {
			for( j = 0; j < 3; j++ ) dki[j] = arg->xyz[ak+j] - arg->xyz[ai+j];
			for( j = 0; j < 3; j++ ) dlj[j] = arg->xyz[al+j] - arg->xyz[aj+j];
			dvt[0] = ( vt[1] * dkj[2] - dkj[1] * vt[2] ) / ( rt2 * rkj );
			dvt[1] = ( vt[2] * dkj[0] - dkj[2] * vt[0] ) / ( rt2 * rkj );
			dvt[2] = ( vt[0] * dkj[1] - dkj[0] * vt[1] ) / ( rt2 * rkj );
			dvu[0] = ( vu[1] * dkj[2] - dkj[1] * vu[2] ) / ( ru2 * rkj );
			dvu[1] = ( vu[2] * dkj[0] - dkj[2] * vu[0] ) / ( ru2 * rkj );
			dvu[2] = ( vu[0] * dkj[1] - dkj[0] * vu[1] ) / ( ru2 * rkj );
			arg->grd[arg->who+ai]   += ( dkj[2] * dvt[1] - dkj[1] * dvt[2] ) * dph;
			arg->grd[arg->who+ai+1] += ( dkj[0] * dvt[2] - dkj[2] * dvt[0] ) * dph;
			arg->grd[arg->who+ai+2] += ( dkj[1] * dvt[0] - dkj[0] * dvt[1] ) * dph;
			arg->grd[arg->who+aj]   += ( dki[1] * dvt[2] - dki[2] * dvt[1] - dlk[2] * dvu[1] + dlk[1] * dvu[2] ) * dph;
			arg->grd[arg->who+aj+1] += ( dki[2] * dvt[0] - dki[0] * dvt[2] - dlk[0] * dvu[2] + dlk[2] * dvu[0] ) * dph;
			arg->grd[arg->who+aj+2] += ( dki[0] * dvt[1] - dki[1] * dvt[0] - dlk[1] * dvu[0] + dlk[0] * dvu[1] ) * dph;
			arg->grd[arg->who+ak]   += ( dji[2] * dvt[1] - dji[1] * dvt[2] - dlj[1] * dvu[2] + dlj[2] * dvu[1] ) * dph;
			arg->grd[arg->who+ak+1] += ( dji[0] * dvt[2] - dji[2] * dvt[0] - dlj[2] * dvu[0] + dlj[0] * dvu[2] ) * dph;
			arg->grd[arg->who+ak+2] += ( dji[1] * dvt[0] - dji[0] * dvt[1] - dlj[0] * dvu[1] + dlj[1] * dvu[0] ) * dph;
			arg->grd[arg->who+al]   += ( - dkj[2] * dvu[1] + dkj[1] * dvu[2] ) * dph;
			arg->grd[arg->who+al+1] += ( - dkj[0] * dvu[2] + dkj[2] * dvu[0] ) * dph;
			arg->grd[arg->who+al+2] += ( - dkj[1] * dvu[0] + dkj[0] * dvu[1] ) * dph;
		}
	}
	return( NULL );
}

static PyObject* w_energy_dihedral( PyObject *self, PyObject *args ) {
	PyObject	*gradient, *object, *molecule, *otmp;
	double		*xyz, *grd, tmp;
	long		i, j, n3, cpu;
	long		*lst, n_lst, n_dat, *ind;
	double		*dat, out = 0.0;
	long		*rng, dsp, nit;
	pthread_t	*pid;
	ene_arg		*arg;

	if( PyArg_ParseTuple( args, "OOO", &object, &molecule, &gradient ) ) {
//		cpu  = PyInt_AsLong( PyObject_GetAttrString( object, "ncpu" ) );
		cpu  = PyLong_AsLong( PyObject_GetAttrString( object, "ncpu" ) );

		otmp = PyObject_GetAttrString( molecule, "coor" );
		n3   = PyList_Size( otmp );
		xyz  = (double*) malloc( n3 * sizeof( double ) );
		for( i = 0; i < n3; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );

		if( gradient != Py_True ) { grd = NULL; }
		else { grd = (double*) malloc( cpu*n3 * sizeof( double ) ); for( i = 0; i < cpu*n3; i++ ) grd[i] = 0.0; }

		otmp  = PyObject_GetAttrString( object, "dihe" );
		n_lst = PyList_Size( otmp );
		lst   = (long*) malloc( 4*n_lst * sizeof( long ) );
		for( i = 0; i < n_lst; i++ )
			for( j = 0; j < 4; j++ ) lst[4*i+j] = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), j ) );
		Py_DECREF( otmp );
	
		otmp  = PyObject_GetAttrString( object, "dihe_data" );
		n_dat = PyList_Size( otmp );
		dat    = (double*) malloc( 12*n_dat * sizeof( double ) );
		for( i = 0; i < 12*n_dat; i++ ) dat[i] = 0.0;
		for( i = 0; i < n_dat; i++ )
			for( j = 0; j < 12; j++ ) dat[12*i+j] = PyFloat_AsDouble( PyList_GetItem( PyList_GetItem( otmp, i ), j ) );
		Py_DECREF( otmp );
		otmp   = PyObject_GetAttrString( object, "dihe_indx" );
		ind    = (long*) malloc( n_lst * sizeof( long ) );
		for( i = 0; i < n_lst; i++ ) ind[i] = PyLong_AsLong( PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );
		nit = n_lst;
		dsp = (long) ((float)nit / (float)cpu);
		rng = (long*) malloc( (cpu+1) * sizeof( long ) );
		if( dsp == 0 ) { for( i = 0; i < cpu+1; i++ ) rng[i] = 0;  for( i = 0; i < nit + 1; i++ ) rng[i] = i; }
		else { for( i = 0; i < cpu; i++ ) rng[i] = i * dsp; rng[cpu] = nit; }
		pid = (pthread_t*) malloc( cpu * sizeof( pthread_t ) );
		arg = (ene_arg*) malloc( cpu * sizeof( ene_arg ) );
		for( i = 0; i < cpu; i++ ) {
			arg[i].ene   = 0.0;
			arg[i].who   = i * n3;
			arg[i]._i0   = rng[i];
			arg[i]._if   = rng[i+1];
			arg[i].xyz   = xyz;
			arg[i].grd   = grd;
			arg[i].lst   = lst;
			arg[i].n_lst = n_lst;
			arg[i].dat   = dat;
			arg[i].n_dat = n_dat;
			arg[i].ind   = ind;
			pthread_create( &pid[i], NULL, __energy_dihedral, (void*) &arg[i] );
		}
		for( i = 0; i < cpu; i++ ) pthread_join( pid[i], NULL );
		for( i = 0; i < cpu; i++ ) out += arg[i].ene;
		free( xyz ); free( rng ); free( pid ); free ( lst ); free( dat ); free( ind ); free( arg );
	
		if( grd != NULL ) {
			otmp = PyObject_GetAttrString( molecule, "grad" );
			for( i = 0; i < n3; i++ ) {
				tmp = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
				for( j = 0; j < cpu; j++ ) tmp += grd[i+j*n3];
				PyList_SetItem( otmp, i, PyFloat_FromDouble( tmp ) );
			}
			Py_DECREF( otmp );
		}

		free( grd );
		return( Py_BuildValue( "d", out ) );
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}


// ####################################################################################################################


void* __energy_non_bonded( void *args ) {
	int_arg		*arg = (int_arg*) args;
	double		epsf = 1389.35484620709144110151 / arg->eps;
	long		i, j, ii, jj, ai, aj;
	double		dr[3], r2, eij, sij, qij, r, s6, tmp, df;
	double		c2on, c2of, _g, _a, _b, _c, _d, _el1, _el2, k6, k12, _lj1, _lj2, r3, r5, s, s3, s12;

	if( arg->con > 0.0 && arg->cof > arg->con ) {
		// atom-based force-switched
		c2on = arg->con * arg->con;
		c2of = arg->cof * arg->cof;
		_g   = pow( c2of - c2on, 3.0 );
		_a   = c2of * c2of * ( c2of - 3.0 * c2on ) / _g;
		_b   = 6.0 * c2of * c2on / _g;
		_c   = - ( c2of + c2on ) / _g;
		_d   = 0.4 / _g;
		_el1 = 8.0 * ( c2of * c2on * ( arg->cof - arg->con ) - 0.2 * ( arg->cof * c2of * c2of - arg->con * c2on * c2on ) ) / _g;
		_el2 = - _a / arg->cof + _b * arg->cof + _c * arg->cof * c2of + _d * arg->cof * c2of * c2of;
		k6   = ( arg->cof * c2of ) / ( arg->cof * c2of - arg->con * c2on );
		k12  = pow( c2of, 3.0 ) / ( pow( c2of, 3.0 ) - pow( c2on, 3.0 ) );
		for( i = arg->_i0; i < arg->_if; i++ ) {
			ii = arg->lst[2*i];
			jj = arg->lst[2*i+1];
			ai = 3 * ii;
			aj = 3 * jj;
			for( j = 0; j < 3; j++ ) dr[j] = arg->xyz[ai+j] - arg->xyz[aj+j];
			for( j = 0; j < 3; j++ ) dr[j] -= arg->box[j] * round( dr[j] / arg->box[j] );
			r2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
			if( r2 > c2of ) { continue; }
			eij =   arg->dat[3*ii] * arg->dat[3*jj];
			sij = arg->dat[3*ii+1] + arg->dat[3*jj+1];
			qij = arg->dat[3*ii+2] * arg->dat[3*jj+2] * epsf;
			r   = sqrt( r2 );
			s   = 1.0 / r;
			s3  = pow( sij * s, 3.0 );
			s6  = s3 * s3;
			if( r2 <= c2on ) {
				tmp = qij * s;
				arg->ele += arg->scl[i] * ( tmp + qij * _el1 );
				s12  = s6 * s6;
				_lj1 = pow( sij / arg->cof * sij / arg->con, 3.0 );
				_lj2 = _lj1 * _lj1;
				arg->vdw += arg->scl[i] * eij * ( ( s12 - _lj2 ) - 2.0 * ( s6 - _lj1 ) );
				df   = ( 12.0 * eij * ( s6 - s12 ) - tmp ) / r2;
			} else {
				r3   =  r * r2;
				r5   = r3 * r2;
				arg->ele += arg->scl[i] * qij * ( _a * s - _b * r - _c * r3 - _d * r5 + _el2 );
				_lj1 = pow( sij / arg->cof, 3.0 );
				_lj2 = _lj1 * _lj1;
				arg->vdw += arg->scl[i] * eij * ( k12 * pow( s6 - _lj2, 2.0 ) - 2.0 * k6 * pow( s3 - _lj1, 2.0 ) );
				df   = - qij * ( _a / r3 + _b * s + 3.0 * _c * r + 5.0 * _d * r3 ) ;
				df  -= 12.0 * eij * ( k12 * s6 * ( s6 - _lj2 ) - k6 * s3 * ( s3 - _lj1 ) ) / r2;
			}
			if( arg->grd != NULL ) {
				for( j = 0; j < 3; j++ ) {
					arg->grd[arg->who+ai+j] += arg->scl[i] * df * dr[j];
					arg->grd[arg->who+aj+j] -= arg->scl[i] * df * dr[j];
				}
			}
		}

	} else {
		// atom-based all-atoms
		for( i = arg->_i0; i < arg->_if; i++ ) {
			ii = arg->lst[2*i];
			jj = arg->lst[2*i+1];
			ai = 3 * ii;
			aj = 3 * jj;
			for( j = 0; j < 3; j++ ) dr[j] = arg->xyz[ai+j] - arg->xyz[aj+j];
			for( j = 0; j < 3; j++ ) dr[j] -= arg->box[j] * round( dr[j] / arg->box[j] );
			r2  = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
			eij =   arg->dat[3*ii] * arg->dat[3*jj];
			sij = arg->dat[3*ii+1] + arg->dat[3*jj+1];
			qij = arg->dat[3*ii+2] * arg->dat[3*jj+2] * epsf;
			s   = 1.0 / sqrt( r2 );
			s6  = pow( sij * s, 6.0 );
			tmp = qij * s;
			arg->ele += arg->scl[i] * tmp;
			arg->vdw += arg->scl[i] * eij * s6 * ( s6 - 2.0 );
			if( arg->grd != NULL ) {
				df = arg->scl[i] * ( 12.0 * eij * s6 * ( 1.0 - s6 ) - tmp ) / r2;
				for( j = 0; j < 3; j++ ) {
					arg->grd[arg->who+ai+j] += df * dr[j];
					arg->grd[arg->who+aj+j] -= df * dr[j];
				}
			}
		}

	}
	return( NULL );
}

static PyObject* w_energy_non_bonded( PyObject *self, PyObject *args ) {
	PyObject	*gradient, *object, *molecule, *otmp, *ptmp, *qtmp;
	double		*grd, *xyz, *dat, *scl, tmp;
	long		i, j, n3, cpu;
	long		*lst, n_lst, n_dat;
	double		oel = 0.0, olj = 0.0, con, cof, box[3], epsi;
	long		*rng, dsp, nit;
	pthread_t	*pid;
	int_arg		*arg;

	if( PyArg_ParseTuple( args, "OOOd", &object, &molecule, &gradient, &epsi ) ) {
//		cpu  = PyInt_AsLong( PyObject_GetAttrString( object, "ncpu" ) );
		cpu  = PyLong_AsLong( PyObject_GetAttrString( object, "ncpu" ) );

		otmp = PyObject_GetAttrString( molecule, "boxl" );
		for( i = 0; i < 3; i++ ) box[i] = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );

		otmp = PyObject_GetAttrString( molecule, "coor" );
		n3   = PyList_Size( otmp );
		xyz  = (double*) malloc( n3 * sizeof( double ) );
		for( i = 0; i < n3; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
		Py_DECREF( otmp );

		otmp  = PyObject_GetAttrString( molecule, "epsi" );
		ptmp  = PyObject_GetAttrString( molecule, "rmin" );
		qtmp  = PyObject_GetAttrString( molecule, "chrg" );
		n_dat = PyList_Size( otmp );
		dat   = (double*) malloc( 3*n_dat * sizeof( double ) );
		for( i = 0; i < n_dat; i++ ) {
			dat[3*i]   = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
			dat[3*i+1] = PyFloat_AsDouble( PyList_GetItem( ptmp, i ) );
			dat[3*i+2] = PyFloat_AsDouble( PyList_GetItem( qtmp, i ) );
		}
		Py_DECREF( otmp ); Py_DECREF( ptmp ); Py_DECREF( qtmp );

		if( gradient != Py_True ) { grd = NULL; }
		else { grd = (double*) malloc( cpu*n3 * sizeof( double ) ); for( i = 0; i < cpu*n3; i++ ) grd[i] = 0.0; }

		otmp = PyObject_GetAttrString( object, "cut_on" );
		con = PyFloat_AsDouble( otmp );
		Py_DECREF( otmp );
	
		otmp = PyObject_GetAttrString( object, "cut_off" );
		cof = PyFloat_AsDouble( otmp );
		Py_DECREF( otmp );
	
		// non_bonded
		otmp  = PyObject_GetAttrString( object, "nbnd" );
		n_lst = PyList_Size( otmp );
		lst   = (long*) malloc( 2*n_lst * sizeof( long ) );
		scl   = (double*) malloc( n_lst * sizeof( double ) );
		for( i = 0; i < n_lst; i++ ) {
			lst[2*i]   = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), 0 ) );
			lst[2*i+1] = PyLong_AsLong( PyList_GetItem( PyList_GetItem( otmp, i ), 1 ) );
			scl[i]     = PyFloat_AsDouble( PyList_GetItem( PyList_GetItem( otmp, i ), 2 ) );
		}
		Py_DECREF( otmp );
	
		nit = n_lst;
		dsp = (long) ((float)nit / (float)cpu);
		rng = (long*) malloc( (cpu+1) * sizeof( long ) );
		if( dsp == 0 ) { for( i = 0; i < cpu+1; i++ ) rng[i] = 0;  for( i = 0; i < nit + 1; i++ ) rng[i] = i; }
		else { for( i = 0; i < cpu; i++ ) rng[i] = i * dsp; rng[cpu] = nit; }
		pid = (pthread_t*) malloc( cpu * sizeof( pthread_t ) );
		arg = (int_arg*) malloc( cpu * sizeof( int_arg ) );
		for( i = 0; i < cpu; i++ ) {
			arg[i].ele   = 0.0;
			arg[i].vdw   = 0.0;
			arg[i].con   = con;
			arg[i].cof   = cof;
			arg[i].eps   = epsi;
			arg[i].who   = i * n3;
			arg[i]._i0   = rng[i];
			arg[i]._if   = rng[i+1];
			arg[i].box   = box;
			arg[i].xyz   = xyz;
			arg[i].grd   = grd;
			arg[i].lst   = lst;
			arg[i].scl   = scl;
			arg[i].n_lst = n_lst;
			arg[i].dat   = dat;
			arg[i].n_dat = n_dat;
			pthread_create( &pid[i], NULL, __energy_non_bonded, (void*) &arg[i] );
		}
		for( i = 0; i < cpu; i++ ) pthread_join( pid[i], NULL );
		for( i = 0; i < cpu; i++ ) { oel += arg[i].ele; olj += arg[i].vdw; }
		free( rng ); free( pid ); free ( lst ); free( scl ); free( arg );
	
		if( grd != NULL ) {
			otmp = PyObject_GetAttrString( molecule, "grad" );
			for( i = 0; i < n3; i++ ) {
				tmp = PyFloat_AsDouble( PyList_GetItem( otmp, i ) );
				for( j = 0; j < cpu; j++ ) tmp += grd[i+j*n3];
				PyList_SetItem( otmp, i, PyFloat_FromDouble( tmp ) );
			}
			Py_DECREF( otmp );
		}

		free( grd );
		return( Py_BuildValue( "(d,d)", oel, olj ) );
	} else { Py_INCREF( Py_None ); return( Py_None ); }
}


// ####################################################################################################################


static struct PyMethodDef methods [] = {
	{ "guess_angles",      (PyCFunction)w_guess_angles,      METH_VARARGS },
	{ "guess_dihedrals",   (PyCFunction)w_guess_dihedrals,   METH_VARARGS },
	{ "update_non_bonded", (PyCFunction)w_update_non_bonded, METH_VARARGS },
	{ "ebond",             (PyCFunction)w_energy_bond,       METH_VARARGS },
	{ "eangle",            (PyCFunction)w_energy_angle,      METH_VARARGS },
	{ "edihedral",         (PyCFunction)w_energy_dihedral,   METH_VARARGS },
	{ "enonbonded",        (PyCFunction)w_energy_non_bonded, METH_VARARGS },
    { 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_mol_mech",
	NULL,
	-1,
	methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC PyInit__mol_mech( void ) {
	PyObject    *my_module;
	my_module = PyModule_Create( &moddef );
	return( my_module );
}
#else
void init_mol_mech( void ) {
	PyObject	*my_module;
	my_module = Py_InitModule( "_mol_mech", methods );
}
#endif

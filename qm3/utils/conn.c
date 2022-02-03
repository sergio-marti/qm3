#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

double r_cov[110] = { 0.31, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
        1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76,
        1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22,
        1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75,
        1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39,
        1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01,
        1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87,
        1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32,
        1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06,
        2.00, 1.96, 1.90, 1.87, 1.80, 1.69, 1.60, 1.60, 1.60, 1.60,
        1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60 };



typedef struct one_lst_node { long   i; struct one_lst_node *n; } one_lst;
typedef struct two_lst_node { long i,j; struct two_lst_node *n; } two_lst;
typedef struct { one_lst *idx; long siz; long *num; double *xyz; two_lst *bnd; } con_arg;



void* __connectivity( void *args ) {
    con_arg		*arg = (con_arg*) args;
    long		i, j, i3, j3;
    double		dr, r2;
    two_lst		*pt2;
    one_lst		*pt1;

    pt2 = arg->bnd;
    for( pt1 = arg->idx; pt1 != NULL ; pt1 = pt1->n ) {
        i = pt1->i;
		for( j = i + 1; j < arg->siz; j++ ) {
    		if( arg->num[i] == 1 && arg->num[j] == 1 ) { continue; }
	    	i3  = 3 * i;
	    	j3  = 3 * j;
	    	r2  = ( r_cov[arg->num[i]] + r_cov[arg->num[j]] + 0.1 ) * ( r_cov[arg->num[i]] + r_cov[arg->num[j]] + 0.1 );
	    	dr  = ( arg->xyz[i3] - arg->xyz[j3] ) * ( arg->xyz[i3] - arg->xyz[j3] ) +
	    			( arg->xyz[i3+1] - arg->xyz[j3+1] ) * ( arg->xyz[i3+1] - arg->xyz[j3+1] ) +
	    			( arg->xyz[i3+2] - arg->xyz[j3+2] ) * ( arg->xyz[i3+2] - arg->xyz[j3+2] );
	    	if( dr <= r2 ) {
	    		pt2->n    = (two_lst*) malloc( sizeof( two_lst ) );
	    		pt2->n->i = i;
	    		pt2->n->j = j;
	    		pt2->n->n = NULL;
	    		pt2       = pt2->n;
	    	}
	    }
	}
    return( NULL );
}


static PyObject* w_connectivity( PyObject *self, PyObject *args ) {
    PyObject	*out, *object, *anum, *coor, *otmp;
    double		*xyz, dr, r2;
    long		*num, i, j, k, n3, nat, cpu, i3, j3;
    pthread_t	*pid;
    con_arg		*arg;
    one_lst     *lst, *pt1;
    two_lst		*pt2;

    if( PyArg_ParseTuple( args, "lO", &cpu, &object ) ) {

    	coor = PyObject_GetAttrString( object, "coor" );
    	n3   = PyList_Size( coor );
    	nat  = n3 / 3;
    	xyz  = (double*) malloc( n3 * sizeof( double ) );
    	for( i = 0; i < n3; i++ ) xyz[i] = PyFloat_AsDouble( PyList_GetItem( coor, i ) );
    	Py_DECREF( coor );

    	anum = PyObject_GetAttrString( object, "anum" );
    	num  = (long*) malloc( nat * sizeof( long ) );
    	for( i = 0; i < nat; i++ ) num[i] = PyLong_AsLong( PyList_GetItem( anum, i ) );
    	Py_DECREF( anum );

        if( nat >= cpu * 2 ) {
            // ======================================================================
            // threaded version
            lst = (one_lst*) malloc( cpu * sizeof( one_lst ) );
            for( i = 0; i < cpu; i++ ) lst[i].n = NULL;
            i = 0;
            j = nat - 2;
            k = 0;
            while( i < j ) {
                pt1 = &lst[k%cpu];
                while( pt1->n != NULL ) pt1 = pt1->n;
                pt1->n = (one_lst*) malloc( sizeof( one_lst ) );
                pt1->n->i = i;
                pt1->n->n = NULL;
                pt1 = pt1->n;
                pt1->n = (one_lst*) malloc( sizeof( one_lst ) );
                pt1->n->i = j;
                pt1->n->n = NULL;
                i++; j--; k++;
            }
            if( i == j ) {
                pt1 = &lst[k%cpu];
                while( pt1->n != NULL ) pt1 = pt1->n;
                pt1->n = (one_lst*) malloc( sizeof( one_lst ) );
                pt1->n->i = i;
                pt1->n->n = NULL;
            }

//j = 0; for( i = 0; i < cpu; i++ ) { printf( "[%2ld]: ", i ); k = 0; pt1 = lst[i].n; while( pt1 != NULL ) { k += nat - 1 - pt1->i; pt1 = pt1->n; } printf( "%ld\n", k ); j += k; } printf( ">> %ld vs %ld\n", j, nat * ( nat - 1 ) / 2 );

            pid = (pthread_t*) malloc( cpu * sizeof( pthread_t ) );
            arg = (con_arg*) malloc( cpu * sizeof( con_arg ) );
            for( i = 0; i < cpu; i++ ) {
                arg[i].siz    = nat;
                arg[i].idx    = lst[i].n;
	    		arg[i].num    = num;
	    		arg[i].xyz    = xyz;
	    		arg[i].bnd    = (two_lst*) malloc( sizeof( two_lst ) ); 
	    		arg[i].bnd->n = NULL;
                pthread_create( &pid[i], NULL, __connectivity, (void*) &arg[i] );
            }
            for( i = 0; i < cpu; i++ ) pthread_join( pid[i], NULL );
    
            out = PyList_New( 0 );
            for( i = 0; i < cpu; i++ ) {
                pt2 = arg[i].bnd->n;
                while( pt2 != NULL ) {
                    otmp = Py_BuildValue( "[l,l]", pt2->i, pt2->j );
                    PyList_Append( out, otmp );
    				Py_DECREF( otmp );
                    pt2 = pt2->n;
                }
            }

            for( i = 0; i < cpu; i++ ) {
                pt1 = lst[i].n;
                while( pt1 != NULL ) {
                    pt1 = pt1->n;
                    free( lst[i].n );
                    lst[i].n = pt1;
                }
            }
            free( lst );

            for( i = 0; i < cpu; i++ ) {
                pt2 = arg[i].bnd;
                while( pt2 != NULL ) {
                    pt2 = pt2->n;
                    free( arg[i].bnd );
                    arg[i].bnd = pt2;
                }
            }
            free( arg ); free( pid );

        } else {
            // ======================================================================
            // serial version
            out = PyList_New( 0 );
            for( i = 0; i < nat - 1; i++ ) {
                for( j = i + 1; j < nat; j++ ) {
                    if( num[i] == 1 && num[j] == 1 ) { continue; }
	    			i3  = 3 * i;
			    	j3  = 3 * j;
			    	r2  = ( r_cov[num[i]] + r_cov[num[j]] + 0.1 ) * ( r_cov[num[i]] + r_cov[num[j]] + 0.1 );
			    	dr  = ( xyz[i3] - xyz[j3] ) * ( xyz[i3] - xyz[j3] ) +
			    			( xyz[i3+1] - xyz[j3+1] ) * ( xyz[i3+1] - xyz[j3+1] ) +
			    			( xyz[i3+2] - xyz[j3+2] ) * ( xyz[i3+2] - xyz[j3+2] );
			    	if( dr <= r2 ) {
                    	otmp = Py_BuildValue( "[l,l]", i, j );
	                    PyList_Append( out, otmp );
	    				Py_DECREF( otmp );
                    }
                }
            }
        }

    	free( num ); free( xyz );
    	return( out );
    } else { Py_INCREF( Py_None ); return( Py_None ); }
}




static struct PyMethodDef methods [] = {
    { "connectivity", (PyCFunction)w_connectivity, METH_VARARGS },
    { 0, 0, 0 }
};



static struct PyModuleDef moddef = {
    PyModuleDef_HEAD_INIT,
    "_conn",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__conn( void ) {
    PyObject    *my_module;
    my_module = PyModule_Create( &moddef );
    return( my_module );
}

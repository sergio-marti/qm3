#include<Python.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdio.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


double get_func( long nx, double dx, double* grdx, long ny, double dy, double* grdy, double* coor ) {
	double	out = 0.0, tmp_i, tmp_j;
	long	i, j;

	for( i = 0; i < nx; i++ ) {
		for( j = 0; j < ny; j++ ) {

			if( i == 0 ) {
				tmp_i = grdx[i*ny+j] - 1.0 / dx * ( coor[(i+1)*ny+j] - coor[i*ny+j] );
			} else if( i == nx - 1 ) {
				tmp_i = grdx[i*ny+j] - 1.0 / dx * ( coor[i*ny+j] - coor[(i-1)*ny+j] );
			} else {
				tmp_i = grdx[i*ny+j] - 0.5 / dx * ( coor[(i+1)*ny+j] - coor[(i-1)*ny+j] );
			}

			if( j == 0 ) {
				tmp_j = grdy[i*ny+j] - 1.0 / dy * ( coor[i*ny+j+1] - coor[i*ny+j] );
			} else if( j == ny - 1 ) {
				tmp_j = grdy[i*ny+j] - 1.0 / dy * ( coor[i*ny+j] - coor[i*ny+j-1] );
			} else {
				tmp_j = grdy[i*ny+j] - 0.5 / dy * ( coor[i*ny+j+1] - coor[i*ny+j-1] );
			}

			out += tmp_i * tmp_i + tmp_j * tmp_j;
		}
	}
	return( out );
}


static PyObject* w_lsfe_2d( PyObject *self, PyObject *args ) {
	long		nx, ny;
	double		dx, dy, delt, c_bak, f_bak, f_for, func;
	long		size, step_number, print_frequency;
	double		tolerance, step_size;
	PyObject 	*object, *output;
	double		*coor, *grad, *step, *velo, *grdx, *grdy;
	double		norm, grms, tmp, alph, ssiz, vsiz, vfac;
	long		it, k, nstp;

	if( PyArg_ParseTuple( args, "ldldOdd", &nx, &dx, &ny, &dy, &object, &step_size, &delt ) ) {

		size = nx * ny;
		coor = (double*) malloc( sizeof(double) * size );
		grad = (double*) malloc( sizeof(double) * size );
		step = (double*) malloc( sizeof(double) * size );
		velo = (double*) malloc( sizeof(double) * size );

		grdx = (double*) malloc( sizeof(double) * size );
		grdy = (double*) malloc( sizeof(double) * size );

		for( k = 0; k < size; k++ ) {
			coor[k] = 0.0;
			step[k] = 0.0;
			velo[k] = 0.0;
			grdx[k] = PyFloat_AsDouble( PyList_GetItem( PyList_GetItem( object, k ), 0 ) );
			grdy[k] = PyFloat_AsDouble( PyList_GetItem( PyList_GetItem( object, k ), 1 ) );
		}
		
		if( delt < 0.0 ) delt = min( 1.0e-4, min( dx, dy ) / 100.0 );
		printf( "---- Least-Squares Finite Elements Integration / 2D (FIRE)\n" );
		printf( "nx = %20ld\nny = %20ld\n", nx, ny );
		printf( "dx = %20.10lf\ndy = %20.10lf\ndq = %20.10lf\n", dx, dy, delt );
		norm = 0.0;
		grms = 0.0;
		func = get_func( nx, dx, grdx, ny, dy, grdy, coor );
		for( k = 0; k < size; k++ ) {
			c_bak   = coor[k];
			coor[k] = c_bak + delt;
			f_for   = get_func( nx, dx, grdx, ny, dy, grdy, coor );
			coor[k] = c_bak - delt;
			f_bak   = get_func( nx, dx, grdx, ny, dy, grdy, coor );
			coor[k] = c_bak;
			grad[k] = 0.5 * ( f_for - f_bak ) / delt;
			norm   += grad[k] * grad[k];
			grms    = max( grms, fabs( grad[k] ) );
		}
		norm = sqrt( norm );

		step_number = 100000;
		print_frequency = 100;
		tolerance = 0.01;

		nstp = 0;
		ssiz = step_size;
		alph = 0.1;

		it = 0;
		printf( "ss = %20.10lf\nfg = %20.10lf\n", step_size, tolerance );
		printf( "----------------------------------------------------------------------\n" );
		printf( "%10s%20s%20s%20s\n", "Iter", "Function", "Gradient", "Step Size" );
		printf( "%10ld%20.5lf%20.10lf%20.10lf\n", it, func, grms, ssiz );
		while( it < step_number && ( func > tolerance || grms > tolerance ) ) {
			vsiz = 0.0; vfac = 0.0;
			for( k = 0; k < size; k++ ) { vsiz += velo[k] * velo[k]; vfac -= velo[k] * grad[k]; }
			vsiz = sqrt( vsiz );
			if( vfac > 0.0 ) {
				for( k = 0; k < size; k++ ) velo[k] = ( 1.0 - alph ) * velo[k] - alph * grad[k] / norm * vsiz;
				if( nstp > 5 ) {
					ssiz = min( ssiz * 1.1, step_size );
					alph *= 0.99;
				}
				nstp++;
			} else {
				for( k = 0; k < size; k++ ) velo[k] = 0.0;
				alph = 0.1;
				ssiz *= 0.5;
				nstp = 0;
			}
			for( k = 0; k < size; k++ ) {
				velo[k] -= ssiz * grad[k];
				step[k] = ssiz * velo[k];
			}
			tmp = 0.0; for( k = 0; k < size; k++ ) tmp += step[k] * step[k]; tmp = sqrt( tmp );
			if( tmp > ssiz ) { for( k = 0; k < size; k++ ) step[k] *= ssiz / tmp; }
			for( k = 0; k < size; k++ ) { coor[k] += step[k]; }

			norm = 0.0;
			grms = 0.0;
			func = get_func( nx, dx, grdx, ny, dy, grdy, coor );
			for( k = 0; k < size; k++ ) {
				c_bak = coor[k];
				coor[k] = c_bak + delt;
				f_for   = get_func( nx, dx, grdx, ny, dy, grdy, coor );
				coor[k] = c_bak - delt;
				f_bak   = get_func( nx, dx, grdx, ny, dy, grdy, coor );
				coor[k] = c_bak;
				grad[k] = 0.5 * ( f_for - f_bak ) / delt;
				norm   += grad[k] * grad[k];
				grms    = max( grms, fabs( grad[k] ) );
			}
			norm = sqrt( norm );

			it++;
			if( it%print_frequency == 0 ) { printf( "%10ld%20.5lf%20.10lf%20.10lf\n", it, func, grms, ssiz ); }

		}

		if( it%print_frequency != 0 ) { printf( "%10ld%20.5lf%20.10lf%20.10lf\n", it, func, grms, ssiz ); }
		printf( "----------------------------------------------------------------------\n" );

		output = PyList_New( size );
		for( k = 0; k < size; k++ ) PyList_SetItem( output, k, PyFloat_FromDouble( coor[k] ) );

		free( coor ); free( grad ); free( grdx ); free( grdy ); free( velo );  free( step );

		return( output );

	} else {

		Py_INCREF( Py_None ); return( Py_None );
	}
}



static struct PyMethodDef methods [] = {
	{ "least_squares_finite_elements_2d",   (PyCFunction)w_lsfe_2d,  METH_VARARGS },
    { 0, 0, 0 }
};



#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moddef = {
	PyModuleDef_HEAD_INIT,
	"_ode",
	NULL,
	-1,
	methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC PyInit__ode( void ) {
	PyObject    *my_module;
	my_module = PyModule_Create( &moddef );
	return( my_module );
}
#else
void init_ode( void ) {
	PyObject	*my_module;
	my_module = Py_InitModule( "_ode", methods );
}
#endif

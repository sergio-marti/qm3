#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


void __swap( long flg, char* txt ) {
    char    s;
    if( flg == 1 ) { s = txt[0]; txt[0] = txt[3]; txt[3] = s; s = txt[1]; txt[1] = txt[2]; txt[2] = s; }
}


typedef struct {
	PyObject_HEAD
	FILE	*fdes;
	long	natm, free, curr, head, crys, *sele, fsiz, swap;
	char	*buff, fnam[2048];
} oDCD;


static int __init( oDCD *self, PyObject *args, PyObject *kwds ) {
	return( 0 );
}


static PyObject* __new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
	oDCD		*self;

	self = (oDCD*) type->tp_alloc( type, 0 );
	self->fdes = NULL;
	self->natm = 0;
	self->free = 0;
	self->curr = 0;
	self->head = 0;
	self->crys = 0;
	self->sele = NULL;
	self->fsiz = 0;
	self->buff = NULL;
	self->swap = 0;
    bzero( self->fnam, 2048 );
	return( (PyObject*) self ) ;
}


static void __dealloc( oDCD *self ) {
	free( self->sele );
	free( self->buff );
	self->fdes = NULL;
	self->natm = 0;
	self->free = 0;
	self->curr = 0;
	self->head = 0;
	self->crys = 0;
	self->sele = NULL;
	self->fsiz = 0;
	self->buff = NULL;
	self->swap = 0;
    bzero( self->fnam, 2048 );
	Py_TYPE( self )->tp_free( (PyObject*) self );
}


static PyObject* __read( PyObject *self, PyObject *args ) {
	PyObject	*o_flg;
	char		*fname, buf[2048];
	oDCD		*obj = NULL;
	int			blk;
	long		i, tmp, fix, nfr = 0;
	struct stat	inf;

	obj = (oDCD*) self;
	o_flg = Py_True;
	if( PyArg_ParseTuple( args, "s|O", &fname, &o_flg ) ) {
		obj->fdes = fopen( fname, "rb" );
		if( obj->fdes == NULL ) { Py_INCREF( Py_None ); return( Py_None ); }
		stat( fname, &inf );
		obj->fsiz = (long) inf.st_size;
		if( o_flg == Py_True ) {
			printf( "* [%s] ", fname );
			for( i = 0; i < 55 - (long)strlen( fname ); i++ ) printf( "-" );
			printf( "\n" );
		}
        fread( buf, 1, 4, obj->fdes ); memcpy( &blk, &buf[0], 4 );
		if( blk == 84 ) { 
			obj->swap = 0;
			if( o_flg == Py_True ) { printf( "+ \tSame Endian\n" ); }
		} else {
			obj->swap = 1;
			if( o_flg == Py_True ) { printf( "+ \tSwapping Endian\n" ); }
		}
		fread( buf, 1, 4, obj->fdes );
		fread( buf, 1, 4, obj->fdes ); __swap( obj->swap, buf ); memcpy( &blk, &buf[0], 4 );
		nfr = (long) blk;
		if( o_flg == Py_True ) { printf( "+ \tNFrames: %ld\n", nfr ); }
		fread( buf, 1, 28, obj->fdes );
		fread( buf, 1, 4, obj->fdes ); __swap( obj->swap, buf ); memcpy( &blk, &buf[0], 4 );
		fix = (long) blk;
		fread( buf, 1, 4, obj->fdes );
		fread( buf, 1, 4, obj->fdes ); __swap( obj->swap, buf ); memcpy( &blk, &buf[0], 4 );
		obj->crys = (long) blk;
		fread( buf, 1, 40, obj->fdes );
		fread( buf, 1, 4, obj->fdes ); __swap( obj->swap, buf ); memcpy( &blk, &buf[0], 4 );
		tmp = 8 + (long) blk;
		fread( buf, 1, tmp, obj->fdes );
		fread( buf, 1, 4, obj->fdes ); __swap( obj->swap, buf ); memcpy( &blk, &buf[0], 4 );
		obj->natm = (long) blk;
		if( o_flg == Py_True ) { printf( "+ \tAtoms: %ld\n", obj->natm ); }
		fread( buf, 1, 4, obj->fdes );
		obj->free = obj->natm - (long)fix;
		if( fix > 0 ) {
			obj->sele = (long*) malloc( obj->free * sizeof( long ) );
			if( o_flg == Py_True ) {
				printf( "+ \tFixed Atoms: %ld\n", (long)fix );
				printf( "+ \tFree  Atoms: %ld\n", obj->free );
			}
			fread( buf, 1, 4, obj->fdes );
			for( i = 0; i < obj->free; i++ ) {
				fread( buf, 1, 4, obj->fdes ); __swap( obj->swap, buf ); memcpy( &blk, &buf[0], 4 );
				obj->sele[i] = (long) blk - 1;

			}
			fread( buf, 1, 4, obj->fdes );
		}
		obj->curr = 0;
		obj->head = ftell( obj->fdes );
		if( o_flg == Py_True )
			{ printf( "------------------------------------------------------------\n" ); }
		obj->buff = (char*) malloc( 4 * obj->natm * sizeof( char ) );
    	bzero( obj->fnam, 2048 );
	}

	return(  Py_BuildValue( "l", nfr ) );
}


static PyObject* __next( PyObject *self, PyObject *args ) {
	PyObject	*o_mol, *o_crd;
    float       blk;
    long        i, j, k;
    oDCD        *obj = NULL;

    obj = (oDCD*) self;
	if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
		if( obj->fdes == NULL ) {
            printf( "* No DCD open...\n" );
            Py_INCREF( Py_False ); return( Py_False );
		}
        if( fread( obj->buff, 1, 4, obj->fdes ) < 4 ) {
            printf( "* End-Of-File reached...\n" );
            Py_INCREF( Py_False ); return( Py_False );
        }
        if( obj->crys == 1 ) {
            if( fread( obj->buff, 1, 56, obj->fdes ) < 56 ) {
            	printf( "* End-Of-File reached...\n" );
                Py_INCREF( Py_False ); return( Py_False );
            }
        }
        if( obj->natm != obj->free && obj->curr > 0 ) {
    		k = 4 * obj->free;
    		if( obj->fsiz - ftell( obj->fdes ) - 20 - 3 * k < 0 ) {
            	printf( "* End-Of-File reached...\n" );
                Py_INCREF( Py_False ); return( Py_False );
    		}
			o_crd = PyObject_GetAttrString( o_mol, "coor" );
    		fread( obj->buff, 1, k, obj->fdes );
    		for( i = 0, j = 0; i < obj->free; i++, j += 4 ) {
				__swap( obj->swap, &(obj->buff[j]) ); memcpy( &blk, &(obj->buff[j]), 4 );
            	PyList_SetItem( o_crd, 3 * obj->sele[i], PyFloat_FromDouble( (double) blk ) );
    		}
    		fread( obj->buff, 1, 8, obj->fdes );
    		fread( obj->buff, 1, k, obj->fdes );
    		for( i = 0, j = 0; i < obj->free; i++, j += 4 ) {
				__swap( obj->swap, &(obj->buff[j]) ); memcpy( &blk, &(obj->buff[j]), 4 );
            	PyList_SetItem( o_crd, 3 * obj->sele[i] + 1, PyFloat_FromDouble( (double) blk ) );
    		}
    		fread( obj->buff, 1, 8, obj->fdes );
    		fread( obj->buff, 1, k, obj->fdes );
    		for( i = 0, j = 0; i < obj->free; i++, j += 4 ) {
				__swap( obj->swap, &(obj->buff[j]) ); memcpy( &blk, &(obj->buff[j]), 4 );
            	PyList_SetItem( o_crd, 3 * obj->sele[i] + 2, PyFloat_FromDouble( (double) blk ) );
    		}
    		fread( obj->buff, 1, 4, obj->fdes );
			Py_DECREF( o_crd );
		} else {
    		k = 4 * obj->natm;
    		if( obj->fsiz - ftell( obj->fdes ) - 20 - 3 * k < 0 ) {
            	printf( "* End-Of-File reached...\n" );
                Py_INCREF( Py_False ); return( Py_False );
    		}
			o_crd = PyObject_GetAttrString( o_mol, "coor" );
    		fread( obj->buff, 1, k, obj->fdes );
    		for( i = 0, j = 0; i < obj->natm; i++, j += 4 ) {
				__swap( obj->swap, &(obj->buff[j]) ); memcpy( &blk, &(obj->buff[j]), 4 );
            	PyList_SetItem( o_crd, 3 * i, PyFloat_FromDouble( (double) blk ) );
    		}
    		fread( obj->buff, 1, 8, obj->fdes );
    		fread( obj->buff, 1, k, obj->fdes );
    		for( i = 0, j = 0; i < obj->natm; i++, j += 4 ) {
				__swap( obj->swap, &(obj->buff[j]) ); memcpy( &blk, &(obj->buff[j]), 4 );
            	PyList_SetItem( o_crd, 3 * i + 1, PyFloat_FromDouble( (double) blk ) );
    		}
    		fread( obj->buff, 1, 8, obj->fdes );
    		fread( obj->buff, 1, k, obj->fdes );
    		for( i = 0, j = 0; i < obj->natm; i++, j += 4 ) {
				__swap( obj->swap, &(obj->buff[j]) ); memcpy( &blk, &(obj->buff[j]), 4 );
            	PyList_SetItem( o_crd, 3 * i + 2, PyFloat_FromDouble( (double) blk ) );
    		}
    		fread( obj->buff, 1, 4, obj->fdes );
			Py_DECREF( o_crd );
		}
    	obj->curr += 1;
	}

    Py_INCREF( Py_True );
    return( Py_True );
}


static PyObject* __close( PyObject *self, PyObject *args ) {
    FILE    *fd;
    char    buf[4];
    oDCD    *obj = NULL;

    obj = (oDCD*) self;
    if( obj->fdes != NULL ) { fclose( obj->fdes ); obj->fdes = NULL; }
    if( obj->fnam[0] != 0 ) {
        fd = fopen( obj->fnam, "rb+" );
        memcpy( &buf[0], &(obj->curr), 4 ); __swap( obj->swap, buf );
        fseek( fd, 8, SEEK_SET );
        fwrite( buf, 1, 4, fd );
        fseek( fd, 20, SEEK_SET );
        fwrite( buf, 1, 4, fd );
		fclose( fd );
    	bzero( obj->fnam, 2048 );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __goto( PyObject *self, PyObject *args ) {
    long        num = 0, dsp;
    oDCD        *obj = NULL;

    obj = (oDCD*) self;
	if( PyArg_ParseTuple( args, "l", &num ) ) {
		if( obj->fdes != NULL && num >= 0 ) {
			if( obj->natm != obj->free && num > 1 ) {
                dsp = 4 + obj->crys * 56 + 3 * 4 * obj->free + 20;
			} else {
                dsp = 4 + obj->crys * 56 + 3 * 4 * obj->natm + 20;
			}
			obj->curr = num;
        	fseek( obj->fdes, obj->head + dsp * num, SEEK_SET );
		}
	}
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __write( PyObject *self, PyObject *args ) {
    PyObject    *o_sele = NULL;
    char        *fname;
    long        natom;
    char        buf[256];
    int         blk;
    long        i, fix;
    oDCD        *obj = NULL;

    obj = (oDCD*) self;
    if( PyArg_ParseTuple( args, "sl|O", &fname, &natom, &o_sele ) ) {
        obj->natm = natom;
        obj->curr = 0;
        snprintf( obj->fnam, 2048, "%s", fname );
fprintf(stderr,"[%s] %ld\n", obj->fnam, obj->natm );

        obj->fdes  = fopen( obj->fnam, "wb" );
        blk = 84; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        memcpy( &buf[0], "CORD", 4 );         fwrite( buf, 1, 4, obj->fdes );
        blk = -1; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk =  1; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk =  1; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk = -1; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk =  0; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk =  0; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk =  0; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );

        obj->free = natom;
        obj->sele = NULL;
        fix = 0;
        if( o_sele != NULL && PyList_Check( o_sele ) ) {
            i = (long) PyList_Size( o_sele );
            if( i > 0 && i < natom ) {
                obj->free = i;
                obj->sele = (long*) malloc( obj->free * sizeof( long ) );
                for( i = 0; i < obj->free; i++ )
                    obj->sele[i] = PyLong_AsLong( PyList_GetItem( o_sele, i ) );
                fix = natom - obj->free;
            }
        }

        blk = (int)( 3 * obj->free ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk = (int) fix; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk = 0; memcpy( &buf[0], &blk, 4 ); for( i = 0; i < 11; i++ ) fwrite( buf, 1, 4, obj->fdes );
        blk = 84; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes ); fwrite( buf, 1, 4, obj->fdes );
        blk = 1; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        memcpy( &buf[0], "* Created with QMCube (qm3)                                                     ", 80 );
        fwrite( buf, 1, 80, obj->fdes );
        blk = 84; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk = 4; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk = (int) natom; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
        blk = 4; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
		if( obj->sele != NULL ) {
        	blk = (int)( 4 * obj->free ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
			for( i = 0; i < obj->free; i++ ) {
				blk = (int)( obj->sele[i] + 1 ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
			}
        	blk = (int)( 4 * obj->free ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fdes );
		}

		obj->buff = (char*) malloc( 4 * obj->natm * sizeof( char ) );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* __append( PyObject *self, PyObject *args ) {
	PyObject	*o_mol, *o_crd;
    char        buf[4];
    int         cnt;
    float       blk;
    long        i, j;
    oDCD        *obj = NULL;

    obj = (oDCD*) self;
	if( PyArg_ParseTuple( args, "O", &o_mol ) ) {
	    if( obj->fdes == NULL ) {
	        printf( "* No DCD open...\n" );
	        Py_INCREF( Py_False ); return( Py_False );
	    }
		o_crd = PyObject_GetAttrString( o_mol, "coor" );
        if( obj->natm != obj->free && obj->curr > 0 ) {
            cnt = (int)( 4 * obj->free );
            memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fdes );
            for( i = 0, j = 0; i < obj->free; i++, j += 4 ) {
            	blk = (float) PyFloat_AsDouble( PyList_GetItem( o_crd, 3 * obj->sele[i] ) );
				memcpy( &(obj->buff[j]), &blk, 4 );
            }
			fwrite( obj->buff, 1, cnt, obj->fdes );
            fwrite( buf, 1, 4, obj->fdes ); fwrite( buf, 1, 4, obj->fdes );
            for( i = 0, j = 0; i < obj->free; i++, j += 4 ) {
            	blk = (float) PyFloat_AsDouble( PyList_GetItem( o_crd, 3 * obj->sele[i] + 1 ) );
				memcpy( &(obj->buff[j]), &blk, 4 );
            }
			fwrite( obj->buff, 1, cnt, obj->fdes );
            fwrite( buf, 1, 4, obj->fdes ); fwrite( buf, 1, 4, obj->fdes );
            for( i = 0, j = 0; i < obj->free; i++, j += 4 ) {
            	blk = (float) PyFloat_AsDouble( PyList_GetItem( o_crd, 3 * obj->sele[i] + 2 ) );
				memcpy( &(obj->buff[j]), &blk, 4 );
            }
			fwrite( obj->buff, 1, cnt, obj->fdes );
            fwrite( buf, 1, 4, obj->fdes );
        } else {
            cnt = (int)( 4 * obj->natm );
            memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fdes );
            for( i = 0, j = 0; i < obj->natm; i++, j += 4 ) {
            	blk = (float) PyFloat_AsDouble( PyList_GetItem( o_crd, 3 * i ) );
				memcpy( &(obj->buff[j]), &blk, 4 );
            }
			fwrite( obj->buff, 1, cnt, obj->fdes );
            fwrite( buf, 1, 4, obj->fdes ); fwrite( buf, 1, 4, obj->fdes );
            for( i = 0, j = 0; i < obj->natm; i++, j += 4 ) {
            	blk = (float) PyFloat_AsDouble( PyList_GetItem( o_crd, 3 * i + 1 ) );
				memcpy( &(obj->buff[j]), &blk, 4 );
            }
			fwrite( obj->buff, 1, cnt, obj->fdes );
            fwrite( buf, 1, 4, obj->fdes ); fwrite( buf, 1, 4, obj->fdes );
            for( i = 0, j = 0; i < obj->natm; i++, j += 4 ) {
            	blk = (float) PyFloat_AsDouble( PyList_GetItem( o_crd, 3 * i + 2 ) );
				memcpy( &(obj->buff[j]), &blk, 4 );
            }
			fwrite( obj->buff, 1, cnt, obj->fdes );
            fwrite( buf, 1, 4, obj->fdes );
        }
		Py_DECREF( o_crd );
		obj->curr++;
	}
    Py_INCREF( Py_True );
    return( Py_True );
}


static struct PyMethodDef __methods [] = {
    { "open_read", (PyCFunction)__read, METH_VARARGS },
    { "open_write", (PyCFunction)__write, METH_VARARGS },
    { "close", (PyCFunction)__close, METH_VARARGS },
    { "next", (PyCFunction)__next, METH_VARARGS },
    { "goto", (PyCFunction)__goto, METH_VARARGS },
    { "append", (PyCFunction)__append, METH_VARARGS },
	{ 0, 0, 0 }
};


static struct PyMemberDef __members [] = {
        { "curr", T_LONG, offsetof( oDCD, curr ), 0 },
        { "natm", T_LONG, offsetof( oDCD, natm ), 0 },
        { 0, 0, 0, 0 }
};


// --------------------------------------------------------------------------------------


static struct PyMethodDef methods [] = {
	{ 0, 0, 0 }
};


static PyTypeObject tDCD = {
	PyVarObject_HEAD_INIT( NULL, 0 )
	.tp_name = "dcd",
	.tp_doc = "(fast) DCD IO support",
	.tp_basicsize = sizeof( oDCD ),
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
	"_dcd",
	NULL,
	-1,
	methods
};

PyMODINIT_FUNC PyInit__dcd( void ) {
	PyObject    *my_module;

	my_module = PyModule_Create( &moddef );
	PyType_Ready( &tDCD );
    Py_INCREF( &tDCD );
    PyModule_AddObject( my_module, "dcd", (PyObject *) &tDCD );
	return( my_module );
}

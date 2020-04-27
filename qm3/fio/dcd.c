#include <Python.h>
#include <structmember.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))


void __swap( long flg, char* txt ) {
    char    s;
    if( flg == 1 ) { s = txt[0]; txt[0] = txt[3]; txt[3] = s; s = txt[1]; txt[1] = txt[2]; txt[2] = s; }
}


typedef struct {
    PyObject_HEAD
    PyObject    *x, *y, *z;
    long        natm, nfre, curr, fcry, fswp;
    long        *sele;
    FILE        *fd;
    char        fnam[2048];
} oDCD;


static int __init( oDCD *self, PyObject *args, PyObject *kwds ) {
    bzero( self->fnam, 2048 );
    self->natm = 0;
    self->nfre = 0;
    self->curr = 0;
    self->fcry = 0;
    self->fswp = 0;
    self->x    = NULL;
    self->y    = NULL;
    self->z    = NULL;
    self->sele = NULL;
    self->fd   = NULL;
    return( 0 );
}


static PyObject* __new( PyTypeObject *type, PyObject *args, PyObject *kwds ) {
    oDCD        *self;

    self = (oDCD*) type->tp_alloc( type, 0 );
    return( (PyObject*) self ) ;
}


static void __dealloc( oDCD *self ) {
    bzero( self->fnam, 2048 );
    self->natm = 0;
    self->nfre = 0;
    self->curr = 0;
    self->fcry = 0;
    self->fswp = 0;
    self->fd   = NULL;
    free( self->sele );
    Py_TYPE( self )->tp_free( self->x );
    Py_TYPE( self )->tp_free( self->y );
    Py_TYPE( self )->tp_free( self->z );
    Py_TYPE( self )->tp_free( (PyObject*) self );
}


static PyObject* dcd_read( PyObject *self, PyObject *args ) {
    PyObject    *o_flg = Py_True;
    char        *fname;
    long        i, nfr, fix, tmp;
    int         blk;
    char        buf[2048];
    oDCD        *obj = NULL;

    obj = (oDCD*) self;
    if( PyArg_ParseTuple( args, "s|O", &fname, &o_flg ) ) {
        if( o_flg == Py_True ) {
            printf( "* [%s] ", fname );
            for( i = 0; i < 55 - strlen( fname ); i++ ) printf( "-" );
            printf( "\n" );
        }
        if( ( obj->fd = fopen( fname, "rb" ) ) == NULL ) {
            Py_INCREF( Py_None );
            return( Py_None );
        }
        fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
        if( blk != 84 ) {
            obj->fswp = 1;
            if( o_flg == Py_True ) printf( "+ \tSwapping Endian\n" );
        } else { if( o_flg == Py_True ) printf( "+ \tSame Endian\n" ); }
        fread( buf, 1, 4, obj->fd );
        fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
        nfr = (long) blk;
        if( o_flg == Py_True ) printf( "+ \tNFrames: %ld\n", nfr );
        fread( buf, 1, 28, obj->fd );
        fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
        fix = (long) blk;
        fread( buf, 1, 4, obj->fd );
        fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
        obj->fcry = (long) blk;
        fread( buf, 1, 40, obj->fd );
        fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
        tmp = (long) blk;
        fread( buf, tmp + 8, 1, obj->fd );
        fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
        obj->natm = (long) blk;
        if( o_flg == Py_True ) printf( "+ \tAtoms: %ld\n", obj->natm );
        fread( buf, 1, 4, obj->fd );
        obj->nfre = obj->natm - fix;
        if( fix > 0 ) {
            obj->sele = (long*) malloc( obj->nfre * sizeof( long ) );
            if( o_flg == Py_True ) printf( "+ \tFixed Atoms: %ld\n+ \tFree  Atoms: %ld", fix, obj->nfre );
            fread( buf, 1, 4, obj->fd );
            for( i = 0; i < obj->nfre; i++ ) {
                fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
                obj->sele[i] = (long) blk - 1;
            }
            fread( buf, 1, 4, obj->fd );
        }
        obj->curr = 0;
        obj->x = PyList_New( 0 );
        obj->y = PyList_New( 0 );
        obj->z = PyList_New( 0 );
        for( i = 0; i < obj->natm; i++ ) { 
            PyList_Append( obj->x, PyFloat_FromDouble( 9999.0 ) );
            PyList_Append( obj->y, PyFloat_FromDouble( 9999.0 ) );
            PyList_Append( obj->z, PyFloat_FromDouble( 9999.0 ) );
        }
        if( o_flg == Py_True ) printf( "------------------------------------------------------------\n" );
        return( Py_BuildValue( "l", nfr ) );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* dcd_next( PyObject *self, PyObject *args ) {
    float       blk;
    char        buf[80];
    long        i;
    oDCD        *obj = NULL;

    obj = (oDCD*) self;

    if( fread( buf, 1, 4, obj->fd ) < 4 ) {
        Py_INCREF( Py_False );
        return( Py_False );
    }
    if( obj->fcry == 1 ) {
        if( fread( buf, 1, 56, obj->fd ) < 56 ) {
            Py_INCREF( Py_False );
            return( Py_False );
        }
    }
    if( obj->natm != obj->nfre && obj->curr > 0 ) {
        for( i = 0; i < obj->natm; i++ ) {
            fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
            PyList_SetItem( obj->x, obj->sele[i], PyFloat_FromDouble( (double) blk ) );
        }
        fread( buf, 1, 8, obj->fd );
        for( i = 0; i < obj->natm; i++ ) {
            fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
            PyList_SetItem( obj->y, obj->sele[i], PyFloat_FromDouble( (double) blk ) );
        }
        fread( buf, 1, 8, obj->fd );
        for( i = 0; i < obj->natm; i++ ) {
            fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
            PyList_SetItem( obj->z, obj->sele[i], PyFloat_FromDouble( (double) blk ) );
        }
        fread( buf, 1, 4, obj->fd );
    } else {
        for( i = 0; i < obj->natm; i++ ) {
            fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
            PyList_SetItem( obj->x, i, PyFloat_FromDouble( (double) blk ) );
        }
        fread( buf, 1, 8, obj->fd );
        for( i = 0; i < obj->natm; i++ ) {
            fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
            PyList_SetItem( obj->y, i, PyFloat_FromDouble( (double) blk ) );
        }
        fread( buf, 1, 8, obj->fd );
        for( i = 0; i < obj->natm; i++ ) {
            fread( buf, 1, 4, obj->fd ); __swap( obj->fswp, buf ); memcpy( &blk, &buf[0], 4 );
            PyList_SetItem( obj->z, i, PyFloat_FromDouble( (double) blk ) );
        }
        fread( buf, 1, 4, obj->fd );
    }

    obj->curr += 1;

    Py_INCREF( Py_True );
    return( Py_True );
}


static PyObject* dcd_close( PyObject *self, PyObject *args ) {
    FILE    *fd;
    char    buf[4];
    oDCD    *obj = NULL;

    obj = (oDCD*) self;
    if( obj->fd != NULL ) { fclose( obj->fd ); obj->fd = NULL; }
    if( obj->fnam[0] != 0 ) {
        fd = fopen( obj->fnam, "rb+" );
        fseek( fd, 8, SEEK_SET );
        memcpy( &buf[0], &(obj->curr), 4 );
        fwrite( buf, 1, 4, fd );
        fseek( fd, 20, SEEK_SET );
        fwrite( buf, 1, 4, fd );
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* dcd_write( PyObject *self, PyObject *args ) {
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

        obj->fd  = fopen( obj->fnam, "wb" );
        blk = 84; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        memcpy( &buf[0], "CORD", 4 ); fwrite( buf, 1, 4, obj->fd );
        blk = -1; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk =  1; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk =  1; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk = -1; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk =  0; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk =  0; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk =  0; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );

        obj->nfre = natom;
        obj->sele = NULL;
        fix = 0;
        if( PyList_Check( o_sele ) ) {
            i = (long) PyList_Size( o_sele );
            if( i > 0 && i < natom ) {
                obj->nfre = i;
                obj->sele = (long*) malloc( obj->nfre * sizeof( long ) );
                for( i = 0; i < obj->nfre; i++ )
                    obj->sele[i] = PyLong_AsLong( PyList_GetItem( o_sele, i ) );
                fix = natom - obj->nfre;
            }
        }

        blk = (int)( 3 * obj->nfre ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk = (int) fix; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk = 0; memcpy( &buf[0], &blk, 4 ); for( i = 0; i < 11; i++ ) fwrite( buf, 1, 4, obj->fd );
        blk = 84; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd ); fwrite( buf, 1, 4, obj->fd );
        blk = 1; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        memcpy( &buf[0], "* Created with QMCube (qm3)                                                     ", 80 );
        fwrite( buf, 1, 80, obj->fd );
        blk = 84; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk = 4; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk = (int) natom; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk = 4; memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        blk = (int)( 4 * obj->nfre ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        for( i = 0; i < obj->nfre; i++ ) {
            blk = (int)( obj->sele[i] + 1 ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        }
        blk = (int)( 4 * obj->nfre ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );

        obj->x = PyList_New( 0 );
        obj->y = PyList_New( 0 );
        obj->z = PyList_New( 0 );
        for( i = 0; i < obj->natm; i++ ) { 
            PyList_Append( obj->x, PyFloat_FromDouble( 9999.0 ) );
            PyList_Append( obj->y, PyFloat_FromDouble( 9999.0 ) );
            PyList_Append( obj->z, PyFloat_FromDouble( 9999.0 ) );
        }
    }
    Py_INCREF( Py_None );
    return( Py_None );
}


static PyObject* dcd_append( PyObject *self, PyObject *args ) {
    char        buf[4];
    int         cnt;
    float       blk;
    long        i;
    oDCD        *obj = NULL;

    obj = (oDCD*) self;
    if( ftell( obj->fd ) == -1 ) {
        printf( "* No DCD open...\n" );
        Py_INCREF( Py_False );
        return( Py_False );
    }

    if( obj->natm != obj->nfre && obj->curr > 0 ) {
        cnt = (int)( 4 * obj->nfre );
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
        for( i = 0; i < obj->nfre; i++ ) {
            blk = (float) PyFloat_AsDouble( PyList_GetItem( obj->x, obj->sele[i] ) ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        }
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
        for( i = 0; i < obj->nfre; i++ ) {
            blk = (float) PyFloat_AsDouble( PyList_GetItem( obj->y, obj->sele[i] ) ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        }
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
        for( i = 0; i < obj->nfre; i++ ) {
            blk = (float) PyFloat_AsDouble( PyList_GetItem( obj->z, obj->sele[i] ) ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        }
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
    } else {
        cnt = (int)( 4 * obj->natm );
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
        for( i = 0; i < obj->natm; i++ ) {
            blk = (float) PyFloat_AsDouble( PyList_GetItem( obj->x, i ) ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        }
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
        for( i = 0; i < obj->natm; i++ ) {
            blk = (float) PyFloat_AsDouble( PyList_GetItem( obj->y, i ) ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        }
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
        for( i = 0; i < obj->natm; i++ ) {
            blk = (float) PyFloat_AsDouble( PyList_GetItem( obj->z, i ) ); memcpy( &buf[0], &blk, 4 ); fwrite( buf, 1, 4, obj->fd );
        }
        memcpy( &buf[0], &cnt, 4 ); fwrite( buf, 1, 4, obj->fd );
    }

    Py_INCREF( Py_True );
    return( Py_True );
}


static struct PyMethodDef __methods [] = {
    { "read", (PyCFunction)dcd_read, METH_VARARGS },
    { "next", (PyCFunction)dcd_next, METH_VARARGS },
    { "close", (PyCFunction)dcd_close, METH_VARARGS },
    { "write", (PyCFunction)dcd_write, METH_VARARGS },
    { "append", (PyCFunction)dcd_append, METH_VARARGS },
    { 0, 0, 0 }
};


static struct PyMemberDef __members [] = {
        { "X", T_OBJECT, offsetof( oDCD, x ), 0 },
        { "Y", T_OBJECT, offsetof( oDCD, y ), 0 },
        { "Z", T_OBJECT, offsetof( oDCD, z ), 0 },
        { "N", T_LONG, offsetof( oDCD, natm ), 0 },
        { 0, 0, 0, 0 }
};


// --------------------------------------------------------------------------------------


static struct PyMethodDef methods [] = {
    { 0, 0, 0 }
};


#if PY_MAJOR_VERSION >= 3

static PyTypeObject Tdcd = {
    PyVarObject_HEAD_INIT( NULL, 0 )
    .tp_name = "dcd",
    .tp_doc = "DCD",
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
    PyType_Ready( &Tdcd );
    Py_INCREF( &Tdcd );
    PyModule_AddObject( my_module, "dcd", (PyObject *) &Tdcd );
    return( my_module );
}

#else

static PyTypeObject Tdcd = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "dcd",                 // tp_name
    sizeof( oDCD ),        // tp_basicsize
    0,                       // tp_itemsize
    (destructor)__dealloc,   // tp_dealloc
    0,                       // tp_print
    0,                       // tp_getattr
    0,                       // tp_setattr
    0,                       // tp_compare
    0,                       // tp_repr
    0,                       // tp_as_number
    0,                       // tp_as_sequence
    0,                       // tp_as_mapping
    0,                       // tp_hash
    0,                       // tp_call
    0,                       // tp_str
    0,                       // tp_getattro
    0,                       // tp_setattro
    0,                       // tp_as_buffer
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE, // tp_flags
    "DCD",
                             // tp_doc
    0,                       // tp_traverse
    0,                       // tp_clear
    0,                       // tp_richcompare
    0,                       // tp_weaklistoffset
    0,                       // tp_iter
    0,                       // tp_iternext
    __methods,               // tp_methods
    __members,               // tp_members
    0,                       // tp_getset
    0,                       // tp_base
    0,                       // tp_dict
    0,                       // tp_descr_get
    0,                       // tp_descr_set
    0,                       // tp_dictoffset
    (initproc)__init,        // tp_init
    0,                       // tp_alloc
    __new,                   // tp_new
};

void init_dcd( void ) {
    PyObject    *my_module;

    my_module = Py_InitModule( "_dcd", methods );
    PyType_Ready( &Tdcd );
    Py_INCREF( &Tdcd );
    PyModule_AddObject( my_module, "dcd", (PyObject *) &Tdcd );
}

#endif

# -*- coding: iso-8859-1 -*-
import math
import re
import qm3.fio
try:
    from qm3.engines._mmint import *
except:
    pass


def p_sander( data ):
    out = {}
    if( type( data ) == list ):
        lst = data[:]
    else:
        lst = [ data ]
    for fname in lst:
        f = qm3.fio.open_r( fname )
        l = f.readline()
        while( l != "" ):
            if( l[0:4] == "MOD4" or l[0:4] == "NONB" ):
                t = f.readline().split()
                while( len( t ) >= 3 ):
                    out[t[0].upper()] = [ math.sqrt( float( t[2] ) * qm3.constants.K2J ), float( t[1] ) ]
                    t = f.readline().split()
            l = f.readline()
        qm3.fio.close( f, fname )
    return( out )


def p_charmm( data ):
    out = {}
    if( type( data ) == list ):
        lst = data[:]
    else:
        lst = [ data ]
    p = re.compile( "^[\ ]*([^\ ^!]+)[\ ]+[0-9\.]+[\ ]+-([0-9\.]+)[\ ]+([0-9\.]+)" ) 
    for fname in lst:
        f = qm3.fio.open_r( fname )
        l = f.readline()
        while( l != "" ):
            if( p.match( l ) ):
                t = p.findall( l )[0]
                if( len( t ) == 3 ):
                    out[t[0].upper()] = [ math.sqrt( float( t[1] ) * qm3.constants.K2J ), float( t[2] ) ]
            l = f.readline()
        qm3.fio.close( f, fname )
    return( out )


def p_dynamo( data ):
    out = {}
    if( type( data ) == list ):
        lst = data[:]
    else:
        lst = [ data ]
    for fname in lst:
        f = qm3.fio.open_r( fname )
#        l = f.readline()
#        while( l != "" ):
#            if( l[0:4] == "MOD4" ):
#                t = f.readline().split()
#                while( len( t ) >= 3 ):
#                    out[t[0].upper()] = [ math.sqrt( float( t[2] ) * qm3.constants.K2J ), float( t[1] ) ]
#                    t = f.readline().split()
#            l = f.readline()
        qm3.fio.close( f, fname )
    return( out )


def p_default( data ):
    out = {}
    pat = re.compile( "([^\ ^\n]+)[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)" )
    if( type( data ) == list ):
        lst = data[:]
    else:
        lst = [ data ]
    for fname in lst:
        f = qm3.fio.open_r( fname )
        for i,j,k in pat.findall( f.read() ):
            out[i.upper()] = [ math.sqrt( float( j ) * qm3.constants.K2J ), float( k ) ]
        qm3.fio.close( f, fname )
    return( out )


def non_bonded( mol, data, parser = p_default ):
    out = True
    if( mol.type == [] ):
        print( "- Molecule types undefined!" )
        return( False )
    mol.epsi = []
    mol.rmin = []
    nbd = parser( data )
    for i in range( mol.natm ):
        t = mol.type[i].upper()
        if( t in nbd ):
            mol.epsi.append( nbd[t][0] )
            mol.rmin.append( nbd[t][1] )
        else:
            mol.epsi.append( None )
            mol.rmin.append( None )
            print( "- Unknown atom type [%s] at index: %d"%( t, i+1 ) )
            out = False
    return( out )



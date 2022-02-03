# -*- coding: iso-8859-1 -*-
import sys


def open_r( fname ):
    if( fname != None ):
        if( type( fname ) == str ):
            return( open( fname, "rt" ) )
        else:
            return( fname )
    else:
        return( sys.stdin )



def open_w( fname, append = False ):
    if( fname != None ):
        if( type( fname ) == str ):
            if( append ):
                return( open( fname, "at" ) )
            else:
                return( open( fname, "wt" ) )
        else:
            return( fname )
    else:
        return( sys.stdout )



def close( fd, fname ):
    if( fname != None ):
        if( type( fname ) == str ):
            fd.close()



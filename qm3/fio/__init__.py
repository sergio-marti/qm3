# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange



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



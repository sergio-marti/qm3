# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange


import	random as py_random
import	os
import	struct

import	numpy



def random():
#	return( py_random.random() )
#	return( float( struct.unpack( "Q", os.urandom( 8 ) )[0] ) / 18446744073709551616.0 )
	return( numpy.random.random() )


def sample( lst, num ):
#	return( py_random.sample( lst, num ) )
	return( list( numpy.random.choice( lst, num, replace = False ) ) )


def randint( a, b ):
#	return( py_random.randint( a, b ) )
	return( numpy.random.random_integers( a, b ) )


# --- Box–Muller:
# two_pi = 2.0 * 3.14159265358979323846;
# u1 = rand() * ( 1.0 / RAND_MAX );
# u2 = rand() * ( 1.0 / RAND_MAX );
# z0 = sqrt( -2.0 * log( u1 ) ) * cos( two_pi * u2 );
# z1 = sqrt( -2.0 * log( u1 ) ) * sin( two_pi * u2 );
# return( z0 * sigma + mu );
# ----------------------------------------------
def gauss( mean, stdv ):
#	return( py_random.gauss( mean, stdv ) )
	return( numpy.random.normal( mean, stdv ) )


py_random.seed( os.urandom( 128 ) )
numpy.random.seed()

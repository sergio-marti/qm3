# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import	sys
if( sys.version_info[0] == 2 ):
	range = xrange


import	random as py_random
import	os
import	struct
import	math

import	numpy



def random():
#	return( float( struct.unpack( "Q", os.urandom( 8 ) )[0] ) / 18446744073709551616.0 )
#	return( py_random.random() )
	return( numpy.random.random() )


def sample( lst, num ):
#	return( py_random.sample( lst, num ) )
	return( list( numpy.random.choice( lst, num, replace = False ) ) )


def randint( a, b ):
#	return( py_random.randint( a, b ) )
	return( numpy.random.random_integers( a, b ) )


def gauss( mean, stdv ):
#	return( mean + stdv * math.sqrt( -2.0 * math.log( random() ) ) * math.cos( 2.0 * math.pi * random() ) )
#	return( py_random.gauss( mean, stdv ) )
	return( numpy.random.normal( mean, stdv ) )


py_random.seed( os.urandom( 128 ) )
numpy.random.seed()

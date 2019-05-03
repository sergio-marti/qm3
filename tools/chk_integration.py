#!/usr/bin/env python

from __future__ import print_function, division

import	sys
import	qm3.actions.string
import	matplotlib.pyplot as plt
from	matplotlib.backends.backend_pdf import PdfPages


pdf = PdfPages( "dF_cmp.pdf" )

f = open( "string.dat", "rt" )
t = len( f.readline().strip().split() )
n = len( f.readlines() ) + 1
f.close()

crd = int( sys.argv[1] )
win = t // crd
dsp = 1000

for i in range( int( round( n / float( dsp ), 0 ) ) ):
	qm3.actions.string.string_integrate( crd, win, i * dsp, -1 )
	f = open( "string.mfep", "rt" )
	ene = [ float( l.strip() ) for l in f ]
	f.close()
	plt.clf()
	plt.grid( True )
	plt.title( "%d:"%( i * dsp ), loc = "left" )
	plt.plot( [ ene[i] / 4.184 for i in xrange( win ) ], "-o" )
	pdf.savefig()
	plt.close()

pdf.close()

#!/usr/bin/env python3
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

his = []
for i in range( int( round( n / float( dsp ), 0 ) ) ):
	qm3.actions.string.string_integrate( crd, win, i * dsp, -1 )
	f = open( "string.mfep", "rt" )
	ene = [ float( l.strip() ) / 4.184 for l in f ]
	his.append( ene[:] )
	f.close()
	plt.clf()
	plt.grid( True )
	plt.title( "%d:"%( i * dsp ), loc = "left" )
	plt.plot( ene, "-o" )
	pdf.savefig()
	plt.close()

plt.clf()
plt.grid( True )
for i in range( len( his ) ):
	plt.plot( his[i], "-", label = "%d:"%( i * dsp ) )
plt.legend( loc = "lower left" )
pdf.savefig()
plt.close()

pdf.close()

#!/usr/bin/env python3
import  sys
import  os
import  qm3.actions.string
import  matplotlib.pyplot as plt
from    matplotlib.backends.backend_pdf import PdfPages


pdf = PdfPages( "plot.pdf" )


#
# show convergence...
#
f = open( "convergence.dat", "rt" )
y = [ float( l.strip() ) for l in f ]
f.close()

plt.clf()
plt.grid( True )
plt.plot( y )
pdf.savefig()
if( len( sys.argv ) == 1 ):
    plt.show()
plt.close()


#
# set the equilibration time and integrate...
#
crd = 5
win = 60
if( len( sys.argv ) == 2 ):
    rng = [ int( i ) for i in sys.argv[1].split( ":" ) ] 
else:
    sys.stdout.write( "Range [0:-1]: " )
    sys.stdout.flush()
    rng = [ int( i ) for i in sys.stdin.readline().strip().split( ":" ) ] 
qm3.actions.string.string_integrate( crd, win, rng[0], rng[1] )


#
# save and plot the resulting FEP
#
f = open( "string.mfep", "rt" )
ene = [ float( l.strip() ) for l in f ]
f.close()

plt.clf()
plt.grid( True )
plt.plot( [ ene[i] / 4.184 for i in range( win ) ], "-o" )
pdf.savefig()
plt.close()

f = open( "string.colvars", "rt" )
cvs = []
for l in f:
    cvs += [ float( i ) for i in l.strip().split() ]
f.close()

plt.clf()
plt.grid( True )
plt.ylim( 0.5, 4.5 )
for i in range( crd ):
    plt.plot( [ cvs[j*crd+i] for j in range( win ) ], "-o" )
pdf.savefig()
plt.close()

f = open( "string.avevars", "rt" )
cvs = []
for l in f:
    cvs += [ float( i ) for i in l.strip().split() ]
f.close()

plt.clf()
plt.grid( True )
plt.ylim( 0.5, 4.5 )
for i in range( crd ):
    plt.plot( [ cvs[j*crd+i] for j in range( win ) ], "-o" )
pdf.savefig()
plt.close()


##
## plot each coordinate from string.dat (don't show!)
##
f = open( "string.dat", "rt" )
dat = []
for l in f:
    dat += [ float( i ) for i in l.strip().split() ]
f.close()

npt = len( dat ) // ( crd * win )
for i in range( crd ):
    plt.clf()
    plt.grid( True )
    for j in range( win ):
#        plt.plot( [ dat[i+crd*j+k*crd*win] for k in range( npt ) ] )
        plt.plot( [ dat[i+crd*j+k*crd*win] for k in range( rng[0], rng[1] ) ] )
    pdf.savefig()
    plt.close()

pdf.close()

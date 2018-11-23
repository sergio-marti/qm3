import qm3.actions.string
import matplotlib.pyplot as plt
from   matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages( "plot.pdf" )

f = open( "convergence.dat", "rt" )
y = [ float( l.strip() ) for l in f ]
f.close()

plt.clf()
plt.grid( True )
plt.plot( y )
pdf.savefig()
plt.show()
plt.close()

crd = 2
win = 24
#rng = [ int( i ) for i in raw_input( "Range [0:-1]: " ).strip().split( ":" ) ] 
rng = [ 15000, -1 ]
qm3.actions.string.string_integrate( crd, win, rng[0], rng[1] )

f = open( "string.mfep", "rt" )
ene = [ float( l.strip() ) for l in f ]
f.close()

plt.clf()
plt.grid( True )
f, = plt.plot( [ ene[i] / 4.184 for i in range( win ) ], "-o" )
f.set_label( "mfep" )
plt.legend()
pdf.savefig()
plt.close()

f = open( "string.colvars", "rt" )
cvs = []
for l in f:
    cvs += [ float( i ) for i in l.strip().split() ]
f.close()

plt.clf()
plt.grid( True )
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
for i in range( crd ):
    plt.plot( [ cvs[j*crd+i] for j in range( win ) ], "-o" )
pdf.savefig()
plt.close()

pdf.close()

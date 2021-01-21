py_exe=python3

$py_exe << EOD
import  sys
import  os
import  qm3.actions.string
import  matplotlib.pyplot as plt
from    matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages( "04.string.pdf" )
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
win = 64
qm3.actions.string.string_integrate( crd, win, 500, 1500 )
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
plt.ylim( 1.4, 3.6 )
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
plt.ylim( 1.4, 3.6 )
for i in range( crd ):
    plt.plot( [ cvs[j*crd+i] for j in range( win ) ], "-o" )
pdf.savefig()
plt.close()

pdf.close()
EOD


for ff in string.metrics string.avevars \
        string.colvars string.dFdz \
        string.dzds string.mfep; do
    mv -vf $ff 04.$ff
done

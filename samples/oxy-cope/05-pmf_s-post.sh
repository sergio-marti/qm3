py_exe=python3

$py_exe << EOD
import glob
import math
import matplotlib.pyplot as pyplot
import qm3.utils.free_energy

lst = glob.glob( "05.dat_s.??" )
pmf = qm3.utils.free_energy.umbint( lst )
pmf.setup()
pmf.integrate()
f = open( "05.pmf_s", "wt" )
for i in range( len( pmf.crd ) ):
    pmf.pmf[i] /= 4.184
    f.write( "%20.10lf%20.10lf%20.10lf\n"%( pmf.crd[i], pmf.pmf[i], pmf.rms[i] / 4.184 ) )
f.close()

pyplot.clf()
pyplot.grid( True )
pyplot.plot( pmf.crd, pmf.pmf, 'o-' )
pyplot.savefig( "05.pmf_s.pdf" )
pyplot.show()

qm3.utils.free_energy.plot_data( lst )
EOD

mv -vf overlap_coor.pdf 05.pmf_s.overlap_coor.pdf
mv -vf overlap_gaus.pdf 05.pmf_s.overlap_gaus.pdf

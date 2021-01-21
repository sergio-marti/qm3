py_exe=python3

$py_exe << EOD
from matplotlib import pyplot
x = []
for i in range( 24 ):
    f = open( "03.neb.%02d.pdb"%( i ), "rt" )
    x.append( float( f.readline().split()[1] ) / 4.184 )
    f.close()
t = min( x )
x = [ i - t for i in x ]

pyplot.clf()
pyplot.grid( True )
pyplot.plot( x, '-o' )
pyplot.savefig( "03.neb.pdf" )
pyplot.show()
EOD

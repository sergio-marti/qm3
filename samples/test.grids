import qm3.maths.grids
try:
    import cStringIO as io
except:
    import io
import  random

random.seed()
nx = 40
ny = 30
f = io.StringIO()
for i in [ random.uniform( -10, 10 ) for k in range( nx ) ]:
    for j in [ random.uniform( -10, 10 ) for k in range( ny ) ]:
        f.write( "%18.6lf%18.6lf%18.6lf\n"%( i, j, i*i+j*j ) )
f.seek( 0 )
obj = qm3.maths.grids.grid()
obj.regular( f, points = ( nx, ny ), gauss = ( 1.0, 1.0 ) )
f.close()
obj.plot()
print( obj.calc( 4, 4 ) )

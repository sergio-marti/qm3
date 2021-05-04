#!/usr/bin/env python3
import  sys
import  math
import  matplotlib.pyplot as plt

f = open( "04.colvar", "rt" )
vpot, sigm = [ float( i ) for i in f.readline().split() ]
valu = [ float( l.strip() ) for l in f ]
f.close()

plt.clf()
plt.grid( True )
plt.plot( valu, '-' )
plt.savefig( "04.mdyn_coor.pdf" )
plt.show()

# search for the 3rd recrossing
t = input( "range [0:x] = " )
valu = valu[0:int(t)]

rnge = [ -1.5, 1.5 ]
delx = 0.01
nptx = int( round( ( rnge[1] - rnge[0] ) / delx, 0 ) ) + 1
x    = [ rnge[0] + delx * i for i in range( nptx ) ]
y    = []
for i in range( nptx ):
    t = 0.0
    for v in valu:
        t += math.exp( - 0.5 * math.pow( ( x[i] - v ) / sigm, 2.0 ) )
    y.append( - vpot * t / 4.184 )

t = min(y[0:nptx//2])
y = [ i - t for i in y ]

plt.clf()
plt.grid( True )
plt.plot( x, y, '-' )
plt.savefig( "04.mdyn_fes.pdf" )
plt.show()

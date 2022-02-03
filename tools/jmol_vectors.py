#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-
import  sys
import  qm3.mol
import  pickle
import  math
import  os


m = qm3.mol.molecule()
m.xyz_read( "xyz" )
m.guess_atomic_numbers()
s = 3 * m.natm


if( not os.path.isfile( "valvec.pk" ) ):
    f = open( "hessian.dump", "rb" )
    hess = pickle.load( f )
    f.close()
    val, vec = qm3.utils.hessian_frequencies( m.mass, m.coor, hess )
    val = [ i * 0.880 for i in val ]
    f = open( "valvec.pk", "wb" )
    pickle.dump( val, f )
    pickle.dump( vec, f )
    f.close()
else:
    f = open( "valvec.pk", "rb" )
    val = pickle.load( f )
    vec = pickle.load( f )
    f.close()


dsp = 10.0
wid = 0.1
col = "#33cc33"
for w in sys.argv[1:]:
    who = int( w )
    f = open( "jmol.%d"%( who ), "wt" )
#    f.write( "load \"nmode.%03d\"\n"%( who ) )
    f.write( "load \"xyz\"\n" )
    f.write( "cpk 10%\nwireframe 15%\ncolor background white\n" )
    for i in range( m.natm ):
        i3 = 3 * i
        v0 = []
        vf = []
        for k in [0, 1, 2]:
            v0.append( m.coor[i3+k] - dsp * vec[who+s*(i3+k)] )
            vf.append( m.coor[i3+k] + dsp * vec[who+s*(i3+k)] )
        if( math.sqrt( sum( [ (ii-jj)*(ii-jj) for ii,jj in zip( vf, v0 ) ] ) ) >= 1.0 ):
            f.write( "draw id \"a%03d\" arrow { %.3lf %.3lf %.3lf } { %.3lf %.3lf %.3lf } width %.2lf color \"%s\"\n"%( i, v0[0], v0[1], v0[2], vf[0], vf[1], vf[2], wid, col ) )
    f.close()

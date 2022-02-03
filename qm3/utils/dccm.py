# -*- coding: iso-8859-1 -*-
import math
import qm3.fio
import qm3.utils


def dynamic_cross_correlation_map( mol, dcd, fname ):
    ref = []
    cas = []
    sel = []
    for i in range( mol.natm ):
        if( mol.labl[i] == "CA" ):
            ref += mol.coor[3*i:3*i+3][:]
            sel.append( i )
            cas.append( sel.index( i ) )
        elif( mol.labl[i] == "C" or mol.labl[i] == "N" ):
            ref += mol.coor[3*i:3*i+3][:]
            sel.append( i )
    nca = len( cas )
    mas = [ 1.0 for i in sel ]
    dri = []
    nfr = 0
    while( mol.dcd_read( dcd ) ):
        cur = []
        for i in sel:
            cur += mol.coor[3*i:3*i+3][:]
        qm3.utils.superimpose_kabsch( mas, ref, cur )
        for i in cas:
            i3 = 3 * i
            for j in [0, 1, 2]:
                dri.append( cur[i3+j] - ref[i3+j] )
        nfr += 1
    dri2 = []
    for i in range( nca ):
        t = 0.0
        for j in range( nfr ):
            for k in [0, 1, 2]:
                t += dri[j*nca*3+i*3+k] * dri[j*nca*3+i*3+k]
        dri2.append( math.sqrt( t / float( nfr ) ) )
    f = qm3.fio.open_w( fname )
    for i in range( nca ):
        for j in range( nca ):
            t = 0.0
            for k in range( nfr ):
                for l in [0, 1, 2]:
                    t += dri[k*nca*3+i*3+l] * dri[k*nca*3+j*3+l]
            f.write( "%10d%10d%20.10lf\n"%( i+1, j+1, t / ( float( nfr ) * dri2[i] * dri2[j] ) ) )
        f.write( "\n" )
    qm3.fio.close( f, fname )

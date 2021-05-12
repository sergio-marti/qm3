#!/usr/bin/env python3
import  sys
import  pickle
import  struct
import  collections
import  qm3.mol
import  qm3.fio.dcd
import  qm3.engines.mopac
import  _mmdec


def to_center( mol ):
    box = [ mol.coor[0:3][:], mol.coor[0:3][:] ]
    cen = [ .0, .0, .0 ]
    for i in range( mol.natm ):
        i3 = i * 3
        for j in [0, 1, 2 ]:
            box[0][j] = min( box[0][j], mol.coor[i3+j] )
            box[1][j] = max( box[1][j], mol.coor[i3+j] )
            cen[j] += mol.coor[i3+j]
    cen = [ i / mol.natm for i in cen ]
    mol.boxl = [ i-j for i,j in zip( box[1], box[0] ) ]
    for i in range( mol.natm ):
        i3 = i * 3
        for j in [0, 1, 2]:
            mol.coor[i3+j] -= cen[j]


def decompose( mol, res, eqm, emm, fds ):
    mol.func = 0
    eqm.get_func( mol, 1000 )
    print( "    full_QM: %20.10lf"%( mol.func ) )
    chg = mol.chrg[:]
    mol.chrg = [ 0.0 for i in range( mol.natm ) ]
    mol.func = 0
    eqm.get_func( mol, -1 )
    vac = mol.func
    print( "     vac_QM: %20.10lf"%( mol.func ) )
    mol.evdw = [ 0.0 for i in range( mol.natm ) ]
    emm.get_func( mol )
    print( "    full_MM: %20.10lf"%( sum( mol.evdw ) ) )
    _qm = []
    _mm = []
    for k in res.keys():
        mol.chrg = [ 0.0 for i in range( mol.natm ) ]
        for j in res[k]:
            mol.chrg[j] = chg[j]
        mol.func = 0
        eqm.get_func( mol, -1 )
        _qm.append( mol.func )
        _mm.append( sum( [ mol.evdw[j] for j in res[k] ] ) )
        print( "             %20.10lf%20.10lf  %s"%( mol.func - vac, _mm[-1], k ) )
    mol.chrg = chg[:]
    fds[0].write( "%20.10lf\n"%( vac ) )
    fds[0].flush()
    fds[1].write( "".join( [ "%20.10lf"%( i ) for i in _qm ] ) + "\n" )
    fds[1].flush()
    fds[2].write( "".join( [ "%20.10lf"%( i ) for i in _mm ] ) + "\n" )
    fds[2].flush()



f = open( "molec.pk", "rb" )
mol = pickle.load( f )
f.close()

f = open( "data", "rb" )
f.read( 4 )
mol.epsi = list( struct.unpack( "%dd"%( mol.natm ), f.read( mol.natm * 8 ) ) )
f.read( 8 )
mol.rmin = list( struct.unpack( "%dd"%( mol.natm ), f.read( mol.natm * 8 ) ) )
f.close()

mol.coor = [ 0.0 for i in range( 3 * mol.natm ) ]
mol.xyz_read( "xyz", replace = True )
to_center( mol )

f = open( "sele_QM.pk", "rb" )
sqm = pickle.load( f )
f.close()
f = open( "sele_LA.pk", "rb" )
sla = pickle.load( f )
f.close()
f = open( "sele_EX.pk", "rb" )
sex = pickle.load( f )
f.close()
f = open( "sele_MM.pk", "rb" )
smm = pickle.load( f )
f.close()
qqm = int( round( sum( [ mol.chrg[i] for i in sqm ] ), 0 ) )

eqm = qm3.engines.mopac.dl_mopac( mol, "PM3", qqm, 1, sqm, smm, sla, 16.0, 18.5 )
emm = _mmdec.QMLJ( mol, sqm, smm, sex )


# trm = dcd_frames // ncpus  OR  dcd_frames for serial calculations (no arguments)
if( len( sys.argv ) > 1 ):
    cpu = int( sys.argv[1] )
    trm = 500
    fds = [ open( w, "wt" ) for w in [ "__.vac.%d"%( cpu ), "__.ele.%d"%( cpu ), "__.vdw.%d"%( cpu ) ] ]
else:
    cpu = 0
    trm = 1000
    fds = [ open( w, "wt" ) for w in [ "__.vac", "__.ele", "__.vdw" ] ]


res = collections.OrderedDict()
for i in smm:
    if( mol.resn[i] in [ "WAT" ] ):
        key = "BLK"
    elif( mol.resn[i] in [ "Na+" ] ):
        key = "CIO"
    else:
        key = "%s:%s:%d"%( mol.segn[i], mol.resn[i], mol.resi[i] )
    if( not key in [ "A:CYS:244" ] ):
        if( key in res ):
            res[key].append( i )
        else:
            res[key] = [ i ]

#to_center( mol )
#decompose( mol, res, eqm, emm, fds )

dcd = qm3.fio.dcd.dcd()
dcd.open_read( "dcd" )
dcd.goto( cpu * trm )
i = 0
while( dcd.next( mol ) and i < trm ):
    print( "-- frame: ", i )
    to_center( mol )
    decompose( mol, res, eqm, emm, fds )
    i += 1
dcd.close()

for w in fds:
    w.close()

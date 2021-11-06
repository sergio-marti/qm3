#!/bin/bash

prt=enzyme.pdb
lig=ligand.mol2

source ~/Devel/amber/rc
source ~/Devel/openbabel/rc
vina=~/Devel/docking/AutoDock-Vina/build/mac/release/vina

rm -f borra.*
cp $prt borra.pdb
cat > inp << EOD
source leaprc.gaff
source oldff/leaprc.ff03
source leaprc.water.tip3p
prt = loadpdb borra.pdb
savepdb prt borra.pdb
saveamberparm prt borra.prmtop borra.inpcrd
quit
EOD
tleap -f inp
obabel -imol2 $lig -opdbqt | egrep -v "^USER|^TER" > lig.x
obabel -ient borra.pdb -opdbqt | egrep "^ATOM|^HETATM" > borra.pdbqt

python3 << EOD
import  re
import  numpy as np
__frmt = re.compile( "[aAiIeEdD]([0-9]+)" )
# -- parse coordinates and atom labels --
coor = []
indx = {}
cont = 0
f = open( "borra.pdb", "rt" )
for l in f:
    if( l[0:4] == "ATOM" ):
        t = l.split()
        coor += [ float( t[5] ), float( t[6] ), float( t[7] ) ]
        k = "%s-%s-%s"%( t[3], t[4], t[2] )
        if( not k in indx ):
            indx["%s-%s-%s"%( t[3], t[4], t[2] )] = cont
        else:
            print( k )
        cont += 1
f.close()
natm = len( indx )
coor = np.array( coor, dtype = float ).reshape( ( natm, 3 ) )
# -- parse charges and masses --
f = open( "borra.prmtop", "rt" )
l = f.readline()
while( l != "" ):
    if( l[0:12] == "%FLAG CHARGE" ):
        chrg = []
        dsp = int( __frmt.findall( f.readline() )[0] )
        while( len( chrg ) < natm ):
            l = f.readline()
            chrg += [ float( l[i:i+dsp] ) / 18.2223 for i in range( 0, len( l ) - 1, dsp ) ]
    elif( l[0:10] == "%FLAG MASS" ):
        mass = []
        dsp = int( __frmt.findall( f.readline() )[0] )
        while( len( mass ) < natm ):
            l = f.readline()
            mass += [ float( l[i:i+dsp] ) for i in range( 0, len( l ) - 1, dsp ) ]
    l = f.readline()
f.close()
mass = np.array( mass, dtype = float ).reshape( ( natm, 1 ) )
# -- to principal axes --
mc = np.sum( mass * coor, axis = 0 ) / np.sum( mass )
coor -= mc
xx = 0.0; xy = 0.0; xz = 0.0; yy = 0.0; yz = 0.0; zz = 0.0
for i in range( natm ):
    xx += mass[i] * coor[i,0] * coor[i,0]
    xy += mass[i] * coor[i,0] * coor[i,1]
    xz += mass[i] * coor[i,0] * coor[i,2]
    yy += mass[i] * coor[i,1] * coor[i,1]
    yz += mass[i] * coor[i,1] * coor[i,2]
    zz += mass[i] * coor[i,2] * coor[i,2]
rr = np.array( [ yy + zz, -xy, -xz, -xy, xx + zz, -yz, -xz, -yz, xx + yy ] ).reshape( ( 3, 3 ) )
val, vec = np.linalg.eigh( rr )
if( np.linalg.det( vec ) < 0.0 ):
    vec[:,0] *= -1.0
for i in range( natm ):
    coor[i,:] = np.dot( coor[i,:], vec )
# -- flush PDBQT --
f = open( "borra.pdbqt", "rt" )
g = open( "prt.x", "wt" )
for l in f:
    t = l.split()
    w = indx["%s-%s-%s"%( t[3], t[4], t[2] )]
    g.write( l[0:30] + "%8.3lf%8.3lf%8.3lf  0.00  0.00    %6.3lf"%( coor[w,0], coor[w,1], coor[w,2], chrg[w] ) + l[76:] )
g.close()
f.close()
f = open( "borra.pdb", "rt" )
g = open( "prt.pdb", "wt" )
w = 0
for l in f:
    if( l[0:4] == "ATOM" ):
        g.write( l[0:30] + "%8.3lf%8.3lf%8.3lf  0.00  0.00\n"%( coor[w,0], coor[w,1], coor[w,2] ) )
        w += 1
g.close()
f.close()
# -- setup as many centers as needed --
def d_indx( w ):
    t = w // 2
    o = list( range( - t, t + 1 ) )
    if( w % 2 == 0 ):
        del o[o.index( 0 )]
    return( o )
cent = ( np.max( coor, axis = 0 ) - np.min( coor, axis = 0 ) ) // 30 + 1
c = 0
for i in d_indx( int( cent[0] ) ):
    for j in d_indx( int( cent[1] ) ):
        for k in d_indx( int( cent[2] ) ):
            f = open( "inp_%04d"%( c ), "wt" )
            f.write( """receptor = prt.x
ligand = lig.x
center_x = %.1lf
center_y = %.1lf
center_z = %.1lf
size_x = 30.
size_y = 30.
size_z = 30.
energy_range = 10
num_modes = 10
out = out_%04d
"""%( 15.0 * i, 15.0 * j, 15.0 * k, c ) )
            c += 1
EOD

echo "load prt.x, format=pdbqt" > view.pml
rm -f view.vmd
for ff in inp_????; do
	gg=`echo $ff | cut -c5-`
	$vina --cpu 4 --config $ff | tee log_$gg
	cat >> view.vmd << EOD
mol new out_$gg type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation Licorice 0.100000 12.000000 12.000000
mol color Name
mol selection {all}
mol material Opaque
mol addrep top
EOD
	echo "load out_$gg, format=pdbqt" >> view.pml
done

cat >> view.vmd << EOD
mol new prt.x type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation Surf 1.400000 0.000000
mol color Name
mol selection {all}
mol material Opaque
mol addrep top
EOD

cat >> view.pml << EOD
util.cba( 144, "all", _self=cmd )
util.cba(  33, "prt", _self=cmd )
cmd.hide( "cartoon", "prt" )
cmd.show( "surface", "prt" )
set all_states, on
cmd.zoom( "all", animate=-1 )
clip atoms, 5, all
EOD

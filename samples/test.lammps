import qm3.engines._lammps
import os
import time
import qm3.maths.matrix

f = open( "data", "wt" )
f.write( """LAMMPS data file

6 atoms
2 bonds
1 angles
0 dihedrals
0 impropers

2 atom types
1 bond types
1 angle types
0 dihedral types
0 improper types

-50.0 50.0 xlo xhi
-50.0 50.0 ylo yhi
-50.0 50.0 zlo zhi

Masses

1 15.9994
2 1.008

Pair Coeffs

1 0.102 3.188
2 0.000 0.000

Bond Coeffs

1 450 0.9572

Angle Coeffs

1 55.0 104.52

Atoms

1 1 1  0.00000       0.12109      0.06944     -0.22458
2 1 2  0.00000      -0.52694      0.16499     -0.94583
3 1 2  0.00000       0.70159     -0.63565     -0.54677
4 2 1 -0.83400      -0.45114      1.12675      2.21102 
5 2 2  0.41700      -0.29157      0.59483      1.39876 
6 2 2  0.41700       0.05804      1.92714      2.01036 

Bonds

1 1 4 5
2 1 4 6

Angles

1 1 5 4 6
""" )
f.close()

mol = qm3.engines._lammps.lammps_read( "data" )

f = open( "lammps.inp", "wt" )
f.write( """units           real
atom_style      full
boundary        p p p
pair_style      lj/cut/coul/long 10.0
kspace_style    pppm 1e-4
pair_modify     mix arithmetic
bond_style      harmonic
angle_style     harmonic
neighbor        2.0 bin
neigh_modify    delay 10

read_data       data
read_dump       lammps.xyzq 0 x y z q box no format native

group           qmatm id 1-3
neigh_modify    exclude group qmatm qmatm

reset_timestep  0
timestep        1.
thermo_style    multi
thermo          1
run             0
print           $(pe) file lammps.ener screen no
write_dump      all custom lammps.force id fx fy fz modify sort id format line "%d %.10lf %.10lf %.10lf"
""" )
f.close()
eng = qm3.engines._lammps.lammps()
eng.exe = "mpirun -n 4 ./bin/lmp_mpi-16Mar18 -in lammps.inp -sc none -log lammps.log"
mol.func = 0
mol.grad = [ 0 for i in range( 3 * mol.natm ) ]
eng.get_grad( mol )
print( mol.func )
qm3.maths.matrix.mprint( mol.grad, mol.natm, 3 )

print( 80*"-" )
os.unlink( "lammps.ener" )
os.unlink( "lammps.force" )
os.unlink( "lammps.log" )
try:
    os.mkfifo( "lammps.pipe" )
except:
    pass
f = open( "lammps.inp", "wt" )
f.write( """units           real
atom_style      full
boundary        f f f
pair_style      lj/charmm/coul/charmm 4.5 6.0
pair_modify     mix arithmetic
bond_style      harmonic
angle_style     harmonic
neighbor        2.0 bin
neigh_modify    delay 10

read_data       data

group           qmatm id 1-3
neigh_modify    exclude group qmatm qmatm

reset_timestep  0
timestep        1.
thermo_style    multi
thermo          1
""" )
f.close()
os.system( "./bin/lmp_mpi-16Mar18 -in lammps.pipe -sc none -log lammps.log &" )
time.sleep( 10 )
eng = qm3.engines._lammps.lammps_pipe( "lammps.inp" )
mol.func = 0
mol.grad = [ 0 for i in range( 3 * mol.natm ) ]
eng.get_grad( mol )
print( mol.func )
qm3.maths.matrix.mprint( mol.grad, mol.natm, 3 )
eng.stop()

print( 80*"-" )
eng = qm3.engines._lammps.py_lammps( "lammps.inp" )
mol.func = 0
mol.grad = [ 0 for i in range( 3 * mol.natm ) ]
eng.get_grad( mol )
print( mol.func )
qm3.maths.matrix.mprint( mol.grad, mol.natm, 3 )
eng.stop()

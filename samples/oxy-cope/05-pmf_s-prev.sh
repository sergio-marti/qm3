py_exe=python3

cat > 05.pmf_s.cnf << EOD
2 64
dist         1         2
dist         4         5
EOD

cp -vf 04.string.avevars 05.pmf_s.str
cp -vf 04.string.metrics 05.pmf_s.met

$py_exe << EOD | tee 05.pmf_s.rng
import  qm3.mol
import  qm3.engines.colvar_s
mol = qm3.mol.molecule( "04.string.21.pdb" )
qm3.engines.colvar_s.colvar_s( 0.0, 0.0, "05.pmf_s.cnf", "05.pmf_s.str", "05.pmf_s.met", mol )
EOD

$py_exe << EOD
import  os
import  math
import  glob
import  qm3.mol
import  qm3.engines.colvar_s

win = 64
dat = []
mol = qm3.mol.molecule( "04.string.00.pdb" )
eng = qm3.engines.colvar_s.colvar_s( 0.0, 0.0, "05.pmf_s.cnf", "05.pmf_s.str", "05.pmf_s.met", mol )
mol.func = 0.0
dat.append( eng.get_func( mol ) )
for i in range( 1, win ):
    mol.pdb_read( "04.string.%02d.pdb"%( i ) )
    dat.append( eng.get_func( mol ) )
for i in range( win ):
    rx = 0.110703125 * i
    w  = sorted( [ ( math.fabs( dat[j] - rx ), j ) for j in range( win ) ] )[0][1] 
    os.symlink( "04.string.%02d.pdb"%( w ), "05.pmf_s.seed.%02d.pdb"%( i ) )
EOD

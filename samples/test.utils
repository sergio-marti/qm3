import qm3.utils
import qm3.elements
import qm3.maths.matrix

smb = [ "O", "H", "H" ]
mas = [ qm3.elements.mass[qm3.elements.rsymbol[i]] for i in smb ]
crd = [ -0.000050, -0.000050, -0.053402, -0.540040, 0.539989, -0.638619, 0.539989, -0.540040, -0.638619 ]
hes = qm3.maths.matrix.from_upper_diagonal_columns( [ 3209.69558, -3210.83014, 3209.69558, -0.17702, -0.17702,
    4210.58088, -1604.71796, 1605.39887, -1352.55138, 1741.36948, 1605.43127, -1604.97762, 1352.7284,
     -1741.34898, 1741.62915, -1739.13819, 1739.31522, -2105.29043, 1545.8494, -1546.01718, 1995.47763,
     -1604.97762, 1605.43127, 1352.7284, -136.65153, 135.9177, 193.28879, 1741.62915, 1605.39887,
     -1604.71796, -1352.55137, 135.95011, -136.65153, -193.29803, -1741.34898, 1741.36948, 1739.31521,
     -1739.13819, -2105.29043, -193.29803, 193.28879, 109.81281, -1546.01718, 1545.8494, 1995.47763 ], 9 )

print( "Distances:\n", qm3.utils.distance( crd[0:3], crd[3:6] ), qm3.utils.distance( crd[0:3], crd[6:9] ) )

print( "\nAngle:\n", qm3.utils.angle( crd[3:6], crd[0:3], crd[6:9] ) )

tmp = crd[:]
print( "Geometrical-center:\n", qm3.utils.center( [ 1., 1., 1. ], tmp ) )

tmp = crd[:]
print( "\nMass-center:\n", qm3.utils.center( mas, tmp ) )
print( "\nMass-centered coordinates:\n", tmp )

tmp = crd[:]
print( "\nMass-center & rotation-qm3.maths.matrix:\n", qm3.utils.moments_of_inertia( mas, tmp ) )
print( "\nMoments-of-inertia coordinates:\n", tmp )

qm3.utils.superimpose_quaternion( mas, crd, tmp )
print( "\nSuperimposed (quaternion) coordinates:\n", tmp )

qm3.utils.moments_of_inertia( mas, tmp )
qm3.utils.superimpose_kabsch( mas, crd, tmp )
print( "\nSuperimposed (kabsch) coordinates:\n", tmp )

qm3.utils.rotate( tmp, tmp[0:3], [ tmp[3] - tmp[0], tmp[4] - tmp[1], tmp[5] - tmp[2] ], 90. )
print( "\nO->H1 90_deg rotation:\n", tmp )

frq, mds = qm3.utils.hessian_frequencies( mas, crd, hes, True )
print( "\nFrequencies (cm^-1):" )
print( frq )

print( "\nGibbs energy (298K / 1atm): ", qm3.utils.gibbs_rrho( mas, crd, frq ), "_kJ/mol" )

r, F = qm3.utils.force_constants( mas, frq, mds )
print( "\nReduced masses (g/mol):" )
print( r )
print( "\nForce constants (mDyne/A):" )
print( F )

#!/usr/bin/env python3
import  sys
import  re


#
# This script builds a "personalized" mol_mech.prm from a psf and a parameters files
#


if( len( sys.argv ) == 1 ):
    print( "%s psf par"%( sys.argv[0] ) )
    sys.exit( 1 )



kind = []
chrg = []
f = open( sys.argv[1], "rt" )
f.readline()
f.readline()
for i in range( int( f.readline().strip().split()[0] ) + 1 ):
    f.readline()
natm = int( f.readline().strip().split()[0] )
for i in range( natm ):
    t = f.readline().strip().split()
    kind.append( t[5] )
    chrg.append( float( t[6] ) )
f.readline()
n = int( f.readline().strip().split()[0] )
m_bond = []
while( len( m_bond ) < n ):
    t = [ int( i ) - 1 for i in f.readline().strip().split() ]
    for i in range( len( t ) // 2 ):
        m_bond.append( [ t[2*i], t[2*i+1] ] )
f.readline()
n = int( f.readline().strip().split()[0] )
m_angl = []
while( len( m_angl ) < n ):
    t = [ int( i ) - 1 for i in f.readline().strip().split() ]
    for i in range( len( t ) // 3 ):
        m_angl.append( [ t[3*i], t[3*i+1], t[3*i+2] ] )
f.readline()
n = int( f.readline().strip().split()[0] )
m_dihe = []
while( len( m_dihe ) < n ):
    t = [ int( i ) - 1 for i in f.readline().strip().split() ]
    for i in range( len( t ) // 4 ):
        m_dihe.append( [ t[4*i], t[4*i+1], t[4*i+2], t[4*i+3] ] )
f.readline()
n = int( f.readline().strip().split()[0] )
m_impr = []
while( len( m_impr ) < n ):
    t = [ int( i ) - 1 for i in f.readline().strip().split() ]
    for i in range( len( t ) // 4 ):
        m_impr.append( [ t[4*i], t[4*i+1], t[4*i+2], t[4*i+3] ] )
f.close()


p_bnd = re.compile( "^([A-Z0-9\']+)[\ ]+([A-Z0-9\']+)[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)" )
p_ang = re.compile( "^([A-Z0-9\']+)[\ ]+([A-Z0-9\']+)[\ ]+([A-Z0-9\']+)[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)" )
p_die = re.compile( "^([A-Z0-9\']+)[\ ]+([A-Z0-9\']+)[\ ]+([A-Z0-9\']+)[\ ]+([A-Z0-9\']+)[\ ]+([0-9\-\.]+)[\ ]+([1-6])[\ ]+([0-9\.]+)" )
p_imp = re.compile( "^([A-Z0-9\']+)[\ ]+([A-Z0-9\']+)[\ ]+([A-Z0-9\']+)[\ ]+([A-Z0-9\']+)[\ ]+([0-9\.]+)[\ ]+0[\ ]+([0-9\.]+)" )
p_nbn = re.compile( "^([A-Z0-9\']+)[\ ]+[0\.]+[\ ]+-([0-9\.]+)[\ ]+([0-9\.]+)" )

bond = {}
angl = {}
dihe = {}
impr = {}
nbnd = {}
f = open( sys.argv[2], "rt" )
for l in f:
    if( p_bnd.match( l ) ):
        ai,aj,k,r = p_bnd.findall( l )[0]
        bond["%s-%s"%( ai, aj )] = [ float( k ), float( r ) ]
    elif( p_ang.match( l ) ):
        ai,aj,ak,k,r = p_ang.findall( l )[0]
        angl["%s-%s-%s"%( ai, aj, ak )] = [ float( k ), float( r ) ]
    elif( p_die.match( l ) ):
        ai,aj,ak,al,k,n,d = p_die.findall( l )[0]
        t = "%s-%s-%s-%s"%( ai, aj, ak, al )
        if( t in dihe ):
            dihe[t].append( [ float( k ), int( n ), float( d ) ] )
        else:
            dihe[t] = [ [ float( k ), int( n ), float( d ) ] ]
    elif( p_imp.match( l ) ):
        ai,aj,ak,al,k,r = p_imp.findall( l )[0]
        impr["%s-%s-%s-%s"%( ai, aj, ak, al )] = [ float( k ), float( r ) ]
    elif( p_nbn.match( l ) ):
        ai,k,r = p_nbn.findall( l )[0]
        nbnd[ai] = [ float( k ), float( r ) ]
f.close()

ok_bnd = {}
for i,j in m_bond:
    ij = "%s-%s"%( kind[i], kind[j] )
    ji = "%s-%s"%( kind[j], kind[i] )
    if( ij in bond ):
        ok_bnd[ij] = bond[ij]
    elif( ji in bond ):
        ok_bnd[ji] = bond[ji]

ok_ang = {}
for i,j,k in m_angl:
    ijk = "%s-%s-%s"%( kind[i], kind[j], kind[k] )
    kji = "%s-%s-%s"%( kind[k], kind[j], kind[i] )
    if( ijk in angl ):
        ok_ang[ijk] = angl[ijk]
    elif( kji in angl ):
        ok_ang[kji] = angl[kji]

for k in dihe:
    dihe[k].sort( key = lambda var: var[1] )

ok_die = {}
for i,j,k,l in m_dihe:
    ijkl = "%s-%s-%s-%s"%( kind[i], kind[j], kind[k], kind[l] )
    lkji = "%s-%s-%s-%s"%( kind[l], kind[k], kind[j], kind[i] )
    xjkx = "X-%s-%s-X"%( kind[j], kind[k] )
    xkjx = "X-%s-%s-X"%( kind[k], kind[j] )
    if( ijkl in dihe ):
        ok_die[ijkl] = dihe[ijkl]
    elif( lkji in dihe ):
        ok_die[lkji] = dihe[lkji]
    elif( xjkx in dihe ):
        ok_die[xjkx] = dihe[xjkx]
    elif( xkjx in dihe ):
        ok_die[xkjx] = dihe[xkjx]

ok_nbn = { i: nbnd[i] for i in kind }

f = open( "mol_mech.prm", "wt" )
for t in ok_nbn:
    f.write( "%-10s%12.6lf%12.6lf\n"%( t, ok_nbn[t][0], ok_nbn[t][1] ) )
f.write( "\n" )
for t in ok_bnd:
    i,j = t.split( "-" )
    f.write( "%-10s%-10s%12.2lf%12.3lf\n"%( i, j, ok_bnd[t][0], ok_bnd[t][1] ) )
f.write( "\n" )
for t in ok_ang:
    i,j,k = t.split( "-" )
    f.write( "%-10s%-10s%-10s%12.2lf%12.2lf\n"%( i, j, k, ok_ang[t][0], ok_ang[t][1] ) )
f.write( "\n" )
for t in ok_die:
    i,j,k,l = t.split( "-" )
    if( i == "X" ):
        i = "*"
    if( l == "X" ):
        l = "*"
    f.write( "%-10s%-10s%-10s%-10s"%( i, j, k ,l ) +
        "".join( [ "%12.6lf%4d%8.2lf"%( dihe[t][i][0], dihe[t][i][1], dihe[t][i][2] ) for i in range( len( dihe[t] ) ) ] ) + "\n" )
f.close()

f = open( "type_chrg_impr", "wt" )
for i in range( natm ):
    f.write( "%-10s%12.6lf\n"%( kind[i], chrg[i] ) )
f.write( "\n" )
for i,j,k,l in m_impr:
    ijkl = "%s-%s-%s-%s"%( kind[i], kind[j], kind[k], kind[l] )
    lkji = "%s-%s-%s-%s"%( kind[l], kind[k], kind[j], kind[i] )
    ixxl = "%s-X-X%s"%( kind[i], kind[l] )
    if( ijkl in impr ):
        f.write( "[ %d, %d, %d, %d, %.4lf, %.2lf ]\n"%( i, j, k, l, impr[ijkl][0], impr[ijkl][1] ) )
    elif( lkji in impr ):
        f.write( "[ %d, %d, %d, %d, %.4lf, %.2lf ]\n"%( i, j, k, l, impr[lkji][0], impr[lkji][1] ) )
    elif( xjkx in impr ):
        f.write( "[ %d, %d, %d, %d, %.4lf, %.2lf ]\n"%( i, j, k, l, impr[ixxl][0], impr[ixxl][1] ) )
f.close()

# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import qm3.elements
import qm3.fio



__vina = []


def psf_read( mol, fname = None ):
    """
 ------------------------------------------------------------------------
 VINA types
 ------------------------------------------------------------------------
const atom_kind atom_kind_data[] = {
    { "C",    2.00000,    0.15000,   -0.00143,   33.51030,   0.77}, //  0    Non H-bonding Aliphatic Carbon
    { "A",    2.00000,    0.15000,   -0.00052,   33.51030,   0.77}, //  1    Non H-bonding Aromatic Carbon
 ------------------------------------------------------------------------
    { "N",    1.75000,    0.16000,   -0.00162,   22.44930,   0.75}, //  2    Non H-bonding Nitrogen
    {"NA",    1.75000,    0.16000,   -0.00162,   22.44930,   0.75}, //  9    Acceptor 1 H-bond Nitrogen
 ------------------------------------------------------------------------
    { "O",    1.60000,    0.20000,   -0.00251,   17.15730,   0.73}, //  3    Acceptor S Spherical Oxygen
    {"OA",    1.60000,    0.20000,   -0.00251,   17.15730,   0.73}, // 10    Acceptor 2 H-bonds Oxygen
 ------------------------------------------------------------------------
    { "H",    1.00000,    0.02000,    0.00051,    0.00000,   0.37}, //  6    Non H-bonding Hydrogen
    {"HD",    1.00000,    0.02000,    0.00051,    0.00000,   0.37}, // 12    Donor 1 H-bond Hydrogen
 ------------------------------------------------------------------------
    { "P",    2.10000,    0.20000,   -0.00110,   38.79240,   1.06}, //  4
 ------------------------------------------------------------------------
    { "S",    2.00000,    0.20000,   -0.00214,   33.51030,   1.02}, //  5    Non H-bonding Sulphur
    {"SA",    2.00000,    0.20000,   -0.00214,   33.51030,   1.02}, // 11    Acceptor 2 H-bonds Sulphur
 ------------------------------------------------------------------------
    { "F",    1.54500,    0.08000,   -0.00110,   15.44800,   0.71}, //  7
    {"Cl",    2.04500,    0.27600,   -0.00110,   35.82350,   0.99}, // 18
    {"Br",    2.16500,    0.38900,   -0.00110,   42.56610,   1.14}  // 19
    { "I",    2.36000,    0.55000,   -0.00110,   55.05850,   1.33}, //  8
    {"Mg",    0.65000,    0.87500,   -0.00110,    1.56000,   1.30}, // 13
    {"Mn",    0.65000,    0.87500,   -0.00110,    2.14000,   1.39}, // 14
    {"Zn",    0.74000,    0.55000,   -0.00110,    1.70000,   1.31}, // 15
    {"Ca",    0.99000,    0.55000,   -0.00110,    2.77000,   1.74}, // 16
    {"Fe",    0.65000,    0.01000,   -0.00110,    1.84000,   1.25}, // 17
};
 ------------------------------------------------------------------------
    """
    global    __vina
    __vina = []
    conn = []
    f = qm3.fio.open_r( fname )
    if( f.readline().split()[0] == "PSF" ):
        f.readline()
        for i in range( int( f.readline().split()[0] ) + 1 ):
            f.readline()
        if( mol.natm == int( f.readline().split()[0] ) ):
            mol.chrg = []
            mol.anum = []
            for i in range( mol.natm ):
                t = f.readline().split()
                if( mol.segn[i] == t[1] and mol.resi[i] == int( t[2] ) and mol.resn[i] == t[3] and mol.labl[i] == t[4]  ):
                    mol.chrg.append( float( t[6] ) )
                    t[7] = float( t[7] )
                    mol.anum.append( sorted( [ ( math.fabs( qm3.elements.mass[j] - t[7] ), j ) for j in iter( qm3.elements.mass ) if j > 0 ] )[0][1] )
                    __vina.append( qm3.elements.symbol[mol.anum[i]] )
                    conn.append( [] )
                else:
                    print( "- Wrong data (%d): %s/%s %d/%s %s/%s %s/%s"%( i+1, mol.segn[i], t[1], mol.resi[i], t[2], mol.resn[i], t[3], mol.labl[i], t[4] ) )
            # parse connectivity...
            f.readline()
            n = int( f.readline().strip().split()[0] )
            i = 0
            while( i < n ):
                t = [ int( j ) - 1 for j in f.readline().strip().split() ]
                for j in range( 0, len( t ), 2 ):
                    conn[t[j]].append( t[j+1] )
                    conn[t[j+1]].append( t[j] )
                i += len( t ) // 2
        else:
            print( "- Invalid number of atoms in PSF!" )
    qm3.fio.close( f, fname )
    for i in range( mol.natm ):
        if( mol.anum[i] == 1 ):
            if( 6 in [ mol.anum[j] for j in conn[i] ] ):
                mol.chrg[conn[i][0]] += mol.chrg[i]
            else:
                __vina[i] = "HD"
        if( mol.anum[i] == 6 ):
            if( len( conn[i] ) == 3 ):
                if( not 8 in [ mol.anum[j] for j in conn[i] ] ):
                    __vina[i] = "A"
        if( mol.anum[i] == 7 and len( conn[i] ) < 4 ):
            __vina[i] = "NA"
        if( mol.anum[i] == 8 and len( conn[i] ) > 1 ):
            __vina[i] = "OA"
        if( mol.anum[i] == 16 ):
            __vina[i] = "SA"


def receptor( mol, fname = None, sele = [] ):
    global    __vina
    f = qm3.fio.open_w( fname )
    f.write( "REMARK --------------------------------------------------------------------------\n" )
    f.write( "REMARK                            x       y       z     vdW  Elec       q    Type\n" )
    f.write( "REMARK                         _______ _______ _______ _____ _____    ______ ____\n" )
    if( len( sele ) > 0 ):
        t_sel = sorted( sele )
    else:
        t_sel = range( mol.natm )
    k = 0
    for i in t_sel:
        i3 = i * 3
        if( __vina[i] != "H" ):
            f.write( "ATOM  %5d  %-2s  %-3s A%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%10.3lf %-4s\n"%( 
                ( k % 99999 ) + 1, qm3.elements.symbol[mol.anum[i]],
                mol.resn[i][0:min(4,len(mol.resn[i]))], mol.resi[i],
                mol.coor[i3], mol.coor[i3+1], mol.coor[i3+2], 0.0, 0.0, mol.chrg[i], __vina[i] ) )
            k += 1
    qm3.fio.close( f, fname )



def ligand( mol, fname = None, sele = [] ):
    """
    Ligands can be treated as flexible in AutoDock, and we use the idea of a "torsion tree" to represent the
    rigid and rotatable pieces. There is always one "root", and zero or more "branches". Branches can be nested.
    Every branch defines one rotatable bond. The torsion tree is represented in the PDBQT with the following
    records, and the placement of these records is important, and usually means reordering the ATOM/HETATM records:

    The ROOT record precedes the rigid part of the molecule, from which zero or more rotatable bonds may emanate,
    and the ENDROOT record follows the last atom in the rigid "root".
    The ROOT/ENDROOT block of atoms should be given first in the PDBQT file.

    Sets of atoms that are moved by rotatable bonds are enclosed by BRANCH and ENDBRANCH records.
    Both BRANCH records and ENDBRANCH records should give two integers separated by spaces, which are the serial numbers
    of the first and second atoms involved in the rotatable bond.
    It is possible to nest BRANCH/ENDBRANCH blocks.

    The last line of the PDBQT file contains a TORSDOF record, which is followed by an integer. This is the number of
    torsional degrees of freedom in the ligand (the number of rotatable bonds in the ligand).

REMARK                            x       y       z     vdW  Elec       q    Type
REMARK                         _______ _______ _______ _____ _____    ______ ____
ROOT
ATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00    +0.000 C 
ATOM      2  C   UNL     1       0.866   0.866   0.866  0.00  0.00    +0.000 C 
ENDROOT
BRANCH   2   3
ATOM      3  C   UNL     1       1.732   0.000   1.732  0.00  0.00    +0.000 C 
ATOM      4  C   UNL     1       2.598   0.866   2.598  0.00  0.00    +0.000 C 
ENDBRANCH   2   3
TORSDOF 1

    """
    f = qm3.fio.open_w( fname )
    f.write( "REMARK --------------------------------------------------------------------------\n" )
    f.write( "REMARK                            x       y       z     vdW  Elec       q    Type\n" )
    f.write( "REMARK                         _______ _______ _______ _____ _____    ______ ____\n" )
    f.write( "ROOT\n" )
    if( len( sele ) > 0 ):
        t_sel = sorted( sele )
    else:
        t_sel = range( mol.natm )
    k = 0
    for i in t_sel:
        i3 = i * 3
        f.write( "ATOM  %5d  %-2s  %-3s A%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%10.3lf %-4s\n"%( 
            ( k % 99999 ) + 1, qm3.elements.symbol[mol.anum[i]], "LIG", 1,
            mol.coor[i3], mol.coor[i3+1], mol.coor[i3+2], 0.0, 0.0, mol.chrg[i], __vina[i] ) )
        k += 1
    f.write( "ENDROOT\n" )
    f.write( "TORSDOF 0\n" )
    qm3.fio.close( f, fname )


def config( mol, center, fname = "vina.conf" ):
    f = qm3.fio.open_w( fname )
    f.write( """receptor = _RECEPTOR_
ligand = _LIGAND_
center_x = %8.3lf
center_y = %8.3lf
center_z = %8.3lf
size_x = 30.
size_y = 30.
size_z = 30.
energy_range = 5
num_modes = 10"""%( center[0], center[1], center[2] ) )
    qm3.fio.close( f, fname )
    if( type( fname ) == str ):
        print( "\nedit '%s' and run with:\n\t\tvina --cpu N --config %s\n\n"%( fname, fname ) )


def parse_dock( mol, fname, sele, model = 0 ):
    f = qm3.fio.open_r( fname )
    l = f.readline()
    while( l != "" ):
        t = l.strip().split()
        if( t[0] == "MODEL" and int( t[1] ) == model ):
            i = 0
            while( i < len( sele ) ):
                if( t[0] == "ATOM" or t[0] == "HETATM" ):
                    mol.coor[3*sele[i]]   = float( t[6] )
                    mol.coor[3*sele[i]+1] = float( t[7] )
                    mol.coor[3*sele[i]+2] = float( t[8] )
                    i += 1
                t = f.readline().strip().split()
        l = f.readline()
    qm3.fio.close( f, fname )


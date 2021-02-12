# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import os
import pickle
import inspect
import qm3.elements
import qm3.maths.matrix
import qm3.mol
import qm3.utils
import qm3.constants
import qm3.engines.mmres
import qm3.fio

try:
    import qm3.engines._mol_mech
    mol_mech_so = True
    print( ">> mol_mech: binary threaded version" )
except:
    mol_mech_so = False
    print( ">> mol_mech: Python version" )



###################################################################################################
# SFF: Simple Force Field
#
class simple_force_field( object ):
    def __init__( self, mol, ncpu = os.sysconf( 'SC_NPROCESSORS_ONLN' ) ):
        self.ncpu     = ncpu
#        self.cut_on   = 10.0
#        self.cut_off  = 12.0
#        self.cut_list = 14.0
        self.cut_on   = -1
        self.cut_off  = -1
        self.cut_list = -1
        self.natm     = mol.natm
        self.nbnd     = []
        self.qmat     = [ False for i in range( self.natm ) ]
        self.free     = [ True  for i in range( self.natm ) ]
        self.path     = os.path.abspath( os.path.dirname( inspect.getfile( self.__class__ ) ) ) + os.sep


    def topology( self, mol, bond = [], angl = [], dihe = [], impr = [], qtyp = True, qchg = True ):
        """
        impr = [ [ central_i, j, k, l, kmb (kal/mol.rad^2), ref (deg) ], ... ]
        """
        # guess (or not) bonds
        if( len( bond ) > 0 ):
            self.bond = bond[:]
        else:
            self.bond = qm3.utils.connectivity( mol, self.ncpu )
        # build connectivity
        self.conn = [ [] for i in range( self.natm ) ]
        for i,j in self.bond:
            self.conn[i].append( j )
            self.conn[j].append( i )
        # guess (or not) angles
        if( len( angl ) > 0 ):
            self.angl = angl[:]
        else:
            self.x__angles()
        # guess (or not) dihedrals
        if( len( dihe ) > 0 ):
            self.dihe = dihe[:]
        else:
            self.x__dihedrals()
        # copy impropers (if any...)
        self.impr = impr[:]
        # guess (or not) atom types
        if( qtyp ):
            self.x__types( mol )
        # guess (or not) partial charges
        if( qchg and len( mol.chrg ) > 0 ):
            self.x__charges( mol )

    
    # SYBYL atom types (kinda)
    # http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
    # uses FORMAL CHARGES present in MOL.CHRG
    def x__types( self, mol ):
        def __any( lst_a, lst_b ):
            return( len( set( lst_a ).intersection( set( lst_b ) ) ) > 0 )
        mol.type = []
        # default to atomic symbol...
        for i in range( mol.natm ):
            mol.type.append( qm3.elements.symbol[mol.anum[i]] )
            nb = len( self.conn[i] )
            if( mol.anum[i] == 1 ):
                if( mol.anum[self.conn[i][0]] == 6 ):
                    mol.type[i] = "Hn"
            elif( mol.anum[i] == 14 ):
                mol.type[i] = "Si.3"
            elif( mol.anum[i] == 15 ):
                mol.type[i] = "P.3"
            elif( mol.anum[i] in [ 6, 7, 8, 16 ] ):
                mol.type[i] += "_%d"%( len( self.conn[i] ) )
        # 2nd pass
        for i in range( mol.natm ):
            if( mol.type[i] in [ "C_2", "C_1" ] ):
                mol.type[i] = "C.1"
            elif( mol.type[i] == "C_4" ):
                mol.type[i] = "C.3"

            elif( mol.type[i] == "C_3" ):
                if( __any( [ "C_3", "C.ar" ], [ mol.type[j] for j in self.conn[i] ] ) ):
                    mol.type[i] = "C.ar"
                else:
                    if( __any( [ "O_1", "O.co2", "O.2" ], [ mol.type[j] for j in self.conn[i] ] ) ):
                        mol.type[i] = "C.co"
                    else:
                        mol.type[i] = "C.2"
            elif( mol.type[i] == "O_1" ):
                if( mol.type[self.conn[i][0]] in [ "C_3", "C.co" ] and mol.chrg[i] == -0.5 ):
                    mol.type[i] = "O.co2"
                elif( mol.chrg[i] == -1.0 ):
                    mol.type[i] = "O.x"
                else:
                    mol.type[i] = "O.2"
            elif( mol.type[i] == "O_2" ):
                if( 1 in [ mol.anum[j] for j in self.conn[i] ] ):
                    mol.type[i] = "O.h"
                else:
                    mol.type[i] = "O.3"
            elif( mol.type[i] == "N_4" ):
                mol.type[i] = "N.4"
            elif( mol.type[i] == "N_3" ):
                if( mol.chrg[i] == 1.0 ):
                    mol.type[i] = "N.pl"
                else:
                    mol.type[i] = "N.3"
            elif( mol.type[i] == "N_2" ):
                mol.type[i] = "N.2"
            elif( mol.type[i] == "N_1" ):
                mol.type[i] = "N.1"
            elif( mol.type[i] == "S_2" ):
                if( 1 in [ mol.anum[j] for j in self.conn[i] ] ):
                    mol.type[i] = "S.h"
                else:
                    mol.type[i] = "S.3"
            elif( mol.type[i] == "S_1" ):
                if( mol.chrg[i] == -1.0 ):
                    mol.type[i] = "S.x"
                else:
                    mol.type[i] = "S.2"
            elif( mol.type[i] == "S_3" and __any( [ "O_1", "O.2" ], [ mol.type[j] for j in self.conn[i] ] ) ):
                mol.type[i] = "S.o"
            elif( mol.type[i] == "S_4" and __any( [ "O_1", "O.2" ], [ mol.type[j] for j in self.conn[i] ] ) ):
                mol.type[i] = "S.o2"


    def x__angles( self ):
        if( mol_mech_so ):
            self.angl = qm3.engines._mol_mech.guess_angles( self )
        else:
            self.angl = []
            for i in range( len( self.bond ) - 1 ):
                for j in range( i + 1, len( self.bond ) ):
                    if( self.bond[i][0] == self.bond[j][0] ):
                        self.angl.append( [ self.bond[i][1], self.bond[i][0], self.bond[j][1] ] )
                    elif( self.bond[i][0] == self.bond[j][1] ):
                        self.angl.append( [ self.bond[i][1], self.bond[i][0], self.bond[j][0] ] )
                    elif( self.bond[i][1] == self.bond[j][0] ):
                        self.angl.append( [ self.bond[i][0], self.bond[i][1], self.bond[j][1] ] )
                    elif( self.bond[i][1] == self.bond[j][1] ):
                        self.angl.append( [ self.bond[i][0], self.bond[i][1], self.bond[j][0] ] )
            for i in range( len( self.angl )-1, -1, -1 ):
                if( self.angl[i][0] in self.conn[self.angl[i][2]] ):
                    del self.angl[i]


    def x__dihedrals( self ):
        if( mol_mech_so ):
            self.dihe = qm3.engines._mol_mech.guess_dihedrals( self )
        else:
            self.dihe = []
            for i in range( len( self.angl ) - 1 ):
                for j in range( i + 1, len( self.angl ) ):
                    if( self.angl[i][1] == self.angl[j][0] and self.angl[i][2] == self.angl[j][1] ):
                        self.dihe.append( [ self.angl[i][0], self.angl[i][1], self.angl[i][2], self.angl[j][2] ] )
                    elif( self.angl[i][1] == self.angl[j][2] and self.angl[i][2] == self.angl[j][1] ):
                        self.dihe.append( [ self.angl[i][0], self.angl[i][1], self.angl[i][2], self.angl[j][0] ] )
                    elif( self.angl[i][1] == self.angl[j][0] and self.angl[i][0] == self.angl[j][1] ):
                        self.dihe.append( [ self.angl[i][2], self.angl[i][1], self.angl[i][0], self.angl[j][2] ] )
                    elif( self.angl[i][1] == self.angl[j][2] and self.angl[i][0] == self.angl[j][1] ):
                        self.dihe.append( [ self.angl[i][2], self.angl[i][1], self.angl[i][0], self.angl[j][0] ] )
            for i in range( len( self.dihe )-1, -1, -1 ):
                if( self.dihe[i][0] in self.conn[self.dihe[i][3]] ):
                    del self.dihe[i]


    # uses FORMAL CHARGES present in MOL.CHRG
    def x__charges( self, mol ):
        # Electronegativity Equalization Method (B3LYP_6-311G_NPA.par) [10.1186/s13321-015-0107-1]
        f = open( self.path + "mol_mech.eem", "rt" )
        kap = float( f.readline().strip() )
        prm = {}
        for l in f:
            t = l.strip().split()
            prm[t[0]] = [ float( t[1] ), float( t[2] ) ]
        f.close()
        mat = []
        vec = []
        for i in range( mol.natm ):
            for j in range( mol.natm ):
                if( j == i ):
                    mat.append( prm[mol.type[i]][1] )
                else:
                    mat.append( kap / qm3.utils.distance( mol.coor[3*i:3*i+3], mol.coor[3*j:3*j+3] ) )
            mat.append( -1 )
            vec.append( - prm[mol.type[i]][0] )
        mat += [ 1 ] * mol.natm + [ 0 ]
        vec.append( sum( mol.chrg ) )
        mol.chrg = qm3.maths.matrix.solve( mat, vec )[0:mol.natm]


    def psf_read( self, mol, fname ):
        self.impr = []
        impr = []
        fd = qm3.fio.open_r( fname )
        fd.readline()
        fd.readline()
        for i in range( int( fd.readline().strip().split()[0] ) + 1 ):
            fd.readline()
        if( mol.natm == int( fd.readline().split()[0] ) ):
            mol.type = []
            mol.chrg = []
            mol.mass = []
            for i in range( mol.natm ):
                t = fd.readline().strip().split()
                mol.type.append( t[5] )
                mol.chrg.append( float( t[6] ) )
                mol.mass.append( float( t[7] ) )
            fd.readline()
            self.bond = []
            n = int( fd.readline().strip().split()[0] )
            while( len( self.bond ) < n ):
                t = [ int( i ) - 1 for i in fd.readline().strip().split() ]
                for i in range( len( t ) // 2 ):
                    self.bond.append( [ t[2*i], t[2*i+1] ] )
            self.conn = [ [] for i in range( self.natm ) ]
            for i,j in self.bond:
                self.conn[i].append( j )
                self.conn[j].append( i )
            fd.readline()
            self.angl = []
            n = int( fd.readline().strip().split()[0] )
            while( len( self.angl ) < n ):
                t = [ int( i ) - 1 for i in fd.readline().strip().split() ]
                for i in range( len( t ) // 3 ):
                    self.angl.append( [ t[3*i], t[3*i+1], t[3*i+2] ] )
            fd.readline()
            self.dihe = []
            n = int( fd.readline().strip().split()[0] )
            while( len( self.dihe ) < n ):
                t = [ int( i ) - 1 for i in fd.readline().strip().split() ]
                for i in range( len( t ) // 4 ):
                    self.dihe.append( [ t[4*i], t[4*i+1], t[4*i+2], t[4*i+3] ] )
            fd.readline()
            n = int( fd.readline().strip().split()[0] )
            while( len( impr ) < n ):
                t = [ int( i ) - 1 for i in fd.readline().strip().split() ]
                for i in range( len( t ) // 4 ):
                    impr.append( [ t[4*i], t[4*i+1], t[4*i+2], t[4*i+3] ] )
        else:
            print( "- Invalid number of atoms in PSF!" )
        qm3.fio.close( fd, fname )
        return( impr )


    def parameters( self, mol, path = None ):
        if( path != None ):
            self.path = path + os.sep
        out = True
        self.bond_data = []
        self.bond_indx = []
        self.angl_data = []
        self.angl_indx = []
        self.dihe_data = []
        self.dihe_indx = []
        self.impr_data = []
        self.impr_indx = []
        f = open( self.path + "mol_mech.prm", "rt" )
        tmp_typ = {}
        cnt_bnd = 0
        tmp_bnd = {}
        cnt_ang = 0
        tmp_ang = {}
        cnt_dih = 0
        tmp_dih = {}
        cnt_imp = 0
        tmp_imp = {}
        for l in f:
            t = l.strip().split()
            if( len( t ) > 0 and t[0][0] != "#" ):
                if( len( t ) == 3 ):
                    tmp_typ[t[0]] = [ math.sqrt( float( t[1] ) * qm3.constants.K2J ), float( t[2] ) ]
                elif( len( t ) == 4 ):
                    self.bond_data.append( [ float( t[2] ) * qm3.constants.K2J, float( t[3] ) ] )
                    tmp_bnd["%s:%s"%( t[0], t[1] )] = cnt_bnd
                    cnt_bnd += 1
                elif( len( t ) == 5 ):
                    self.angl_data.append( [ float( t[3] ) * qm3.constants.K2J, float( t[4] ) / qm3.constants.R2D ] )
                    tmp_ang["%s:%s:%s"%( t[0], t[1], t[2] )] = cnt_ang
                    cnt_ang += 1
                elif( len( t ) >= 7 ):
                    tmp = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
                    for i in range( 4, len( t ), 3 ):
                        n = int( t[i+1] ) - 1
                        if( n >= 0 and n < 6 ):
                            tmp[2*n]   = float( t[i]   ) * qm3.constants.K2J
                            tmp[2*n+1] = float( t[i+2] ) / qm3.constants.R2D
                    self.dihe_data.append( tmp[:] )
                    tmp_dih["%s:%s:%s:%s"%( t[0], t[1], t[2], t[3] )] = cnt_dih
                    cnt_dih += 1
        f.close()
        mol.epsi = []
        mol.rmin = []
        for i in range( mol.natm ):
            if( mol.type[i] in tmp_typ ):
                mol.epsi.append( tmp_typ[mol.type[i]][0] )
                mol.rmin.append( tmp_typ[mol.type[i]][1] )
            else:
                mol.epsi.append( None )
                mol.rmin.append( None )
                print( "- missing atom type [%s]: %d"%( mol.type[i], i+1 ) )
                out = False
        for i,j in self.bond:
            td = "%s:%s"%( mol.type[i], mol.type[j] )
            ti = "%s:%s"%( mol.type[j], mol.type[i] )
            if( td in tmp_bnd ):
                self.bond_indx.append( tmp_bnd[td] )
            elif( ti in tmp_bnd ):
                self.bond_indx.append( tmp_bnd[ti] )
            else:
                self.bond_indx.append( None )
                print( "- missing parameter [bond]: ", td )
                out = False
        for i,j,k in self.angl:
            td = "%s:%s:%s"%( mol.type[i], mol.type[j], mol.type[k] )
            ti = "%s:%s:%s"%( mol.type[k], mol.type[j], mol.type[i] )
            ts = "*:%s:*"%( mol.type[j] )
            if( td in tmp_ang ):
                self.angl_indx.append( tmp_ang[td] )
            elif( ti in tmp_ang ):
                self.angl_indx.append( tmp_ang[ti] )
            elif( ts in tmp_ang ):
                self.angl_indx.append( tmp_ang[ts] )
            else:
                self.angl_indx.append( None )
                print( "- missing parameter [angl]: ", td )
                out = False
        for i,j,k,l in self.dihe:
            td = "%s:%s:%s:%s"%( mol.type[i], mol.type[j], mol.type[k], mol.type[l] )
            ti = "%s:%s:%s:%s"%( mol.type[l], mol.type[k], mol.type[j], mol.type[i] )
            ts = "*:%s:%s:*"%( mol.type[j], mol.type[k] )
            tz = "*:%s:%s:*"%( mol.type[k], mol.type[j] )
            if( td in tmp_dih ):
                self.dihe_indx.append( tmp_dih[td] )
            elif( ti in tmp_dih ):
                self.dihe_indx.append( tmp_dih[ti] )
            elif( ts in tmp_dih ):
                self.dihe_indx.append( tmp_dih[ts] )
            elif( tz in tmp_dih ):
                self.dihe_indx.append( tmp_dih[tz] )
            else:
                self.dihe_indx.append( None )
                print( "- missing parameter [dihe]: ", td )
                out = False
        return( out )


    def qm_atoms( self, sele ):
        for i in sele:
            self.qmat[i] = True
        # delete QM-QM bonds
        for i in range( len( self.bond ) -1, -1, -1 ):
            if( self.qmat[self.bond[i][0]] and self.qmat[self.bond[i][1]] ):
                del self.bond[i]
                del self.bond_indx[i]
        # delete QM-QM-MM angles
        for i in range( len( self.angl ) -1, -1, -1 ):
            if( sum( [ self.qmat[j] for j in self.angl[i] ] ) >= 2 ):
                del self.angl[i]
                del self.angl_indx[i]
        # delete QM-QM-QM-MM dihedrals
        for i in range( len( self.dihe ) -1, -1, -1 ):
            if( sum( [ self.qmat[j] for j in self.dihe[i] ] ) >= 3 ):
                del self.dihe[i]
                del self.dihe_indx[i]
        # delete QM-QM-QM-MM impropers
        for i in range( len( self.impr ) -1, -1, -1 ):
            if( sum( [ self.qmat[j] for j in self.impr[i][0:4] ] ) >= 3 ):
                del self.impr[i]
        # delete QM-QM non_bonded
        for i in range( len( self.nbnd ) -1, -1, -1 ):
            if( self.qmat[self.nbnd[i][0]] and self.qmat[self.nbnd[i][1]] ):
                del self.nbnd[i]
        

    def x__ebond( self, mol, gradient = False ):
        if( self.bond == [] ):
            return( 0.0 )
        if( mol_mech_so ):
            out = qm3.engines._mol_mech.ebond( self, mol, gradient )
        else:
            bak = mol.func
            out = 0.0
            for i in range( len( self.bond ) ):
                if( self.free[self.bond[i][0]] or self.free[self.bond[i][1]] ):
                    mol.func = 0.0
                    qm3.engines.mmres.mm_bond( mol, self.bond_data[self.bond_indx[i]][0],
                        self.bond_data[self.bond_indx[i]][1],
                        self.bond[i][0], self.bond[i][1],
                        ffac = 1.0, grad = gradient, gfac = [ 1.0, 1.0 ] )
                    out += mol.func
            mol.func = bak
        return( out )


    def x__eangle( self, mol, gradient = False ):
        if( self.angl == [] ):
            return( 0.0 )
        if( mol_mech_so ):
            out = qm3.engines._mol_mech.eangle( self, mol, gradient )
        else:
            bak = mol.func
            out = 0.0
            for i in range( len( self.angl ) ):
                if( self.free[self.angl[i][0]] or self.free[self.angl[i][1]] or self.free[self.angl[i][2]] ):
                    mol.func = 0.0
                    qm3.engines.mmres.mm_angle( mol, self.angl_data[self.angl_indx[i]][0],
                        self.angl_data[self.angl_indx[i]][1],
                        self.angl[i][0], self.angl[i][1], self.angl[i][2],
                        ffac = 1.0, grad = gradient, gfac = [ 1.0, 1.0, 1.0 ] )
                    out += mol.func
            mol.func = bak
        return( out )


    def x__edihedral( self, mol, gradient = False ):
        if( self.dihe == [] ):
            return( 0.0 )
        if( mol_mech_so ):
            out = qm3.engines._mol_mech.edihedral( self, mol, gradient )
        else:
            out = 0.0
            bak = mol.func
            for i in range( len( self.dihe ) ):
                if( self.free[self.dihe[i][0]] or self.free[self.dihe[i][1]] or self.free[self.dihe[i][2]]
                                        or self.free[self.dihe[i][3]] ):
                    mol.func = 0.0
                    qm3.engines.mmres.mm_dihedral( mol, self.dihe_data[self.dihe_indx[i]],
                        self.dihe[i][0], self.dihe[i][1], self.dihe[i][2], self.dihe[i][3],
                        ffac = 1.0, grad = gradient, gfac = [ 1.0, 1.0, 1.0, 1.0 ] )
                    out += mol.func
            mol.func = bak
        return( out )


    def x__eimproper( self, mol, gradient = False ):
        if( self.impr == [] ):
            return( 0.0 )
        out = 0.0
        bak = mol.func
        for i in range( len( self.impr ) ):
            if( self.free[self.impr[i][0]] or self.free[self.impr[i][1]] or self.free[self.impr[i][2]]
                                    or self.free[self.impr[i][3]] ):
                mol.func = 0.0
                qm3.engines.mmres.mm_improper( mol, self.impr[i][4] * qm3.constants.K2J, self.impr[i][5],
                    self.impr[i][0], self.impr[i][1], self.impr[i][2], self.impr[i][3],
                    ffac = 1.0, grad = gradient, gfac = [ 1.0, 1.0, 1.0, 1.0 ] )
                out += mol.func
        mol.func = bak
        return( out )


    def update_non_bonded( self, mol ):
        if( mol_mech_so ):
            self.nbnd = qm3.engines._mol_mech.update_non_bonded( self, mol )
        else:
            self.nbnd = []
            if( self.cut_list > 0.0 ):
                c2l = self.cut_list * self.cut_list
            else:
                c2l = 1.0e99
            for i in range( mol.natm - 1 ):
                i3 = 3 * i
                i_crd = [ mol.coor[i3+k] - mol.boxl[k] * round( mol.coor[i3+k] / mol.boxl[k], 0 ) for k in [ 0, 1, 2 ] ]
                for j in range( i+1, mol.natm ):
                    if( not ( self.qmat[i] and self.qmat[j] ) and ( self.free[i] or self.free[j] ) ):
                        j_crd = [ mol.coor[j*3+k] - mol.boxl[k] * round( mol.coor[j*3+k] / mol.boxl[k], 0 ) for k in [ 0, 1, 2 ] ]
                        if( qm3.utils.distanceSQ( i_crd, j_crd ) <= c2l ):
                            f = False
                            n = len( self.bond )
                            k = 0
                            while( k < n and not f ):
                                f |= ( ( i == self.bond[k][0] and j == self.bond[k][1] ) or
                                        ( i == self.bond[k][1] and j == self.bond[k][0] )  )
                                k += 1
                            n = len( self.angl )
                            k = 0
                            while( k < n and not f ):
                                f |= ( ( i == self.angl[k][0] and j == self.angl[k][2] ) or
                                        ( i == self.angl[k][2] and j == self.angl[k][0] )  )
                                k += 1
                            n = len( self.dihe )
                            k = 0
                            while( k < n and not f ):
                                f |= ( ( i == self.dihe[k][0] and j == self.dihe[k][3] ) or
                                        ( i == self.dihe[k][3] and j == self.dihe[k][0] )  )
                                k += 1
                            if( not f ):
                                self.nbnd.append( [ i, j, 1.0 ] )
            for i,j,k,l in self.dihe:
                if( not ( self.qmat[i] and self.qmat[l] ) ):
                    self.nbnd.append( [ i, l, 0.5 ] )
                        

    def x__enonbonded( self, mol, gradient = False, epsilon = 1.0 ):
        if( not self.nbnd ):
            self.update_non_bonded( mol )
        if( mol_mech_so ):
            oel, olj = qm3.engines._mol_mech.enonbonded( self, mol, gradient, epsilon )
        else:
            epsf = 1389.35484620709144110151 / epsilon
            oel  = 0.0
            olj  = 0.0
            if( self.cut_on > 0.0 and self.cut_off > self.cut_on ):
                c2on =  self.cut_on * self.cut_on
                c2of = self.cut_off * self.cut_off
                _g   = math.pow( c2of - c2on, 3.0 )
                _a   = c2of * c2of * ( c2of - 3.0 * c2on ) / _g
                _b   = 6.0 * c2of * c2on / _g
                _c   = - ( c2of + c2on ) / _g
                _d   = 0.4 / _g
                _el1 = 8.0 * ( c2of * c2on * ( self.cut_off - self.cut_on ) - 0.2 * ( self.cut_off * c2of * c2of - self.cut_on * c2on * c2on ) ) / _g
                _el2 = - _a / self.cut_off + _b * self.cut_off + _c * self.cut_off * c2of + _d * self.cut_off * c2of * c2of
                k6   = ( self.cut_off * c2of ) / ( self.cut_off * c2of - self.cut_on * c2on )
                k12  = math.pow( c2of, 3.0 ) / ( math.pow( c2of, 3.0 ) - math.pow( c2on, 3.0 ) )
                for i,j,f in self.nbnd:
                    ai = i * 3
                    aj = j * 3
                    dr = [ ii-jj for ii,jj in zip( mol.coor[ai:ai+3], mol.coor[aj:aj+3] ) ]
                    dr = [ dr[k] - mol.boxl[k] * round( dr[k] / mol.boxl[k], 0 ) for k in [ 0, 1, 2] ]
                    r2 = sum( [ ii*ii for ii in dr ] )
                    if( r2 > c2of ):
                        continue
                    eij = mol.epsi[i] * mol.epsi[j]
                    sij = mol.rmin[i] + mol.rmin[j]
                    qij = mol.chrg[i] * mol.chrg[j] * epsf * ( not self.qmat[i] ) * ( not self.qmat[j] )
                    r   = math.sqrt( r2 )
                    s   = 1.0 / r
                    s3  = math.pow( sij * s, 3.0 )
                    s6  = s3 * s3
                    if( r2 <= c2on ):
                        tmp = qij * s
                        oel += f * ( tmp + qij * _el1 )
                        s12  = s6 * s6
                        _lj1 = math.pow( sij / self.cut_off * sij / self.cut_on, 3.0 )
                        _lj2 = _lj1 * _lj1
                        olj += f * eij * ( ( s12 - _lj2 ) - 2.0 * ( s6 - _lj1 ) )
                        df   = ( 12.0 * eij * ( s6 - s12 ) - tmp ) / r2
                    else:
                        r3   = r * r2
                        r5   = r3 * r2
                        oel += f * qij * ( _a * s - _b * r - _c * r3 - _d * r5 + _el2 )
                        _lj1 = math.pow( sij / self.cut_off, 3.0 )
                        _lj2 = _lj1 * _lj1
                        olj += f * eij * ( k12 * math.pow( s6 - _lj2, 2.0 ) - 2.0 * k6 * math.pow( s3 - _lj1, 2.0 ) )
                        df   = - qij * ( _a / r3 + _b * s + 3.0 * _c * r + 5.0 * _d * r3 ) 
                        df  -= 12.0 * eij * ( k12 * s6 * ( s6 - _lj2 ) - k6 * s3 * ( s3 - _lj1 ) ) / r2
                    if( gradient ):
                        for j in [0, 1, 2]:
                            mol.grad[ai+j] += f * df * dr[j]
                            mol.grad[aj+j] -= f * df * dr[j]
            else:
                for i,j,f in self.nbnd:
                    ai  = i * 3
                    aj  = j * 3
                    dr  = [ ii-jj for ii,jj in zip( mol.coor[ai:ai+3], mol.coor[aj:aj+3] ) ]
                    dr = [ dr[k] - mol.boxl[k] * round( dr[k] / mol.boxl[k], 0 ) for k in [ 0, 1, 2] ]
                    r2  = sum( [ ii*ii for ii in dr ] )
                    eij = mol.epsi[i] * mol.epsi[j]
                    sij = mol.rmin[i] + mol.rmin[j]
                    qij = mol.chrg[i] * mol.chrg[j] * epsf * ( not self.qmat[i] ) * ( not self.qmat[j] )
                    r   = 1.0 / math.sqrt( r2 )
                    s6  = math.pow( sij * r, 6.0 )
                    tmp = qij * r
                    oel += f * tmp
                    olj += f * eij * s6 * ( s6 - 2.0 )
                    if( gradient ):
                        df = f * ( 12.0 * eij * s6 * ( 1.0 - s6 ) - tmp ) / r2;
                        for j in [0, 1, 2]:
                            mol.grad[ai+j] += df * dr[j]
                            mol.grad[aj+j] -= df * dr[j]
        return( oel, olj )



    def get_func( self, mol, qprint = False ):
        e_bond = self.x__ebond( mol, gradient = False )
        e_angl = self.x__eangle( mol, gradient = False )
        e_dihe = self.x__edihedral( mol, gradient = False )
        e_impr = self.x__eimproper( mol, gradient = False )
        e_elec, e_vdwl = self.x__enonbonded( mol, gradient = False )
        mol.func += e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl
        if( qprint ):
            print( "ETot:", e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl, "_kJ/mol" )
            print( "   Bond:%18.4lf   Angl:%18.4lf   Dihe:%18.4lf"%( e_bond, e_angl, e_dihe ) )
            print( "   Impr:%18.4lf   Elec:%18.4lf   VdWl:%18.4lf"%( e_impr, e_elec, e_vdwl ) )


    def get_grad( self, mol, qprint = False ):
        e_bond = self.x__ebond( mol, gradient = True )
        e_angl = self.x__eangle( mol, gradient = True )
        e_dihe = self.x__edihedral( mol, gradient = True )
        e_impr = self.x__eimproper( mol, gradient = True )
        e_elec, e_vdwl = self.x__enonbonded( mol, gradient = True )
        mol.func += e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl
        if( qprint ):
            print( "ETot:", e_bond + e_angl + e_dihe + e_impr + e_elec + e_vdwl, "_kJ/mol" )
            print( "   Bond:%18.4lf   Angl:%18.4lf   Dihe:%18.4lf"%( e_bond, e_angl, e_dihe ) )
            print( "   Impr:%18.4lf   Elec:%18.4lf   VdWl:%18.4lf"%( e_impr, e_elec, e_vdwl ) )


    # use before calling "qm_atoms" method: store pure MM system
    def system_write( self, fname, mol = None ):
        f = open( fname, "wb" )
        pickle.dump( self.natm, f )
        pickle.dump( self.bond, f )
        pickle.dump( self.conn, f )
        pickle.dump( self.angl, f )
        pickle.dump( self.dihe, f )
        pickle.dump( self.impr, f )
        pickle.dump( self.bond_data, f )
        pickle.dump( self.bond_indx, f )
        pickle.dump( self.angl_data, f )
        pickle.dump( self.angl_indx, f )
        pickle.dump( self.dihe_data, f )
        pickle.dump( self.dihe_indx, f )
        pickle.dump( self.impr_data, f )
        pickle.dump( self.impr_indx, f )
        if( mol != None ):
            if( mol.natm == self.natm ):
                pickle.dump( mol.anum, f )
                pickle.dump( mol.type, f )
                pickle.dump( mol.chrg, f )
                pickle.dump( mol.epsi, f )
                pickle.dump( mol.rmin, f )
        f.close()


    # use before calling "qm_atoms" method: store pure MM system
    def system_read( self, fname, mol = None ):
        f = open( fname, "rb" )
        self.natm = pickle.load( f )
        self.bond = pickle.load( f )
        self.conn = pickle.load( f )
        self.angl = pickle.load( f )
        self.dihe = pickle.load( f )
        self.impr = pickle.load( f )
        self.bond_data = pickle.load( f )
        self.bond_indx = pickle.load( f )
        self.angl_data = pickle.load( f )
        self.angl_indx = pickle.load( f )
        self.dihe_data = pickle.load( f )
        self.dihe_indx = pickle.load( f )
        self.impr_data = pickle.load( f )
        self.impr_indx = pickle.load( f )
        if( mol != None ):
            if( mol.natm == self.natm ):
                mol.anum = pickle.load( f )
                mol.type = pickle.load( f )
                mol.chrg = pickle.load( f )
                mol.epsi = pickle.load( f )
                mol.rmin = pickle.load( f )
        f.close()

# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import os
import math
import qm3.constants
import qm3.fio
import qm3.utils
try:
    import cPickle as pickle
except:
    import pickle



class qmbase( object ):

    def __init__( self, mol, inp, sele, nbnd = [], link = [] ):
        self._cx = qm3.constants.A0
        self._ce = qm3.constants.H2J
        self._cg = self._ce / qm3.constants.A0
        self._ch = self._cg / qm3.constants.A0

        self.sel = sorted( sele )
        self.lnk = link[:]
#        t = [ j for i,j in self.lnk ]
#        self.nbn = sorted( [ i for i in nbnd if not i in t ] )
        self.nbn = sorted( set( nbnd ).difference( set( sele + sum( link, [] ) ) ) )

        self.exe = ""
        self.vla = []
#        self.smb = [ i.title() for i in mol.guess_symbols( sele ) ]
        self.smb = mol.guess_symbols( self.sel )

        if( not mol.chrg ):
            mol.chrg = [ 0.0 for i in range( mol.natm ) ]

        f = qm3.fio.open_r( inp )
        self.inp = f.read()
        qm3.fio.close( f, inp )


    def mk_input( self, mol, run ):
        pass


    def parse_log( self, mol, run ):
        pass


    def get_func( self, mol ):
        self.mk_input( mol, "ener" )
        os.system( self.exe )
        self.parse_log( mol, "ener" )


    def get_grad( self, mol ):
        self.mk_input( mol, "grad" )
        os.system( self.exe )
        self.parse_log( mol, "grad" )


    def get_hess( self, mol ):
        self.mk_input( mol, "hess" )
        os.system( self.exe )
        self.parse_log( mol, "hess" )



# ----------------------------------------------------------------------------------
# Link-Atom stuff
#
def LA_coordinates( qm_i, mm_j, mol, dst = 1.1 ):
    v = [ ii-jj for ii,jj in zip( mol.coor[3*mm_j:3*mm_j+3], mol.coor[3*qm_i:3*qm_i+3] ) ]
    m = math.sqrt( sum( [ k*k for k in v ] ) )
    v = [ k / m * dst for k in v ]
    return( [ mol.coor[3*qm_i+k] - mol.boxl[k] * round( ( mol.coor[3*qm_i] + v[k] ) / mol.boxl[k], 0 ) + v[k] for k in [0, 1, 2] ], [ -v[0], -v[1], -v[2] ] )


def LA_gradient( lnk, grd ):
    for qm_i,mm_j,vec in lnk:
        m = math.sqrt( sum( [ k*k for k in vec ] ) )
        t = sum( [ vec[k] * ( grd[3*qm_i+k] - grd[3*mm_j+k] ) for k in [0, 1, 2] ] ) * 0.5 / m
        grd[3*qm_i:3*qm_i+3] = [ grd[3*qm_i+k] - t * vec[k] / m for k in [0, 1, 2] ]



# ----------------------------------------------------------------------------------
# Exclusions
#
def exclusions( sele_QM, molec, bonds = [] ):
    if( len( bonds ) > 0 ):
        bond = bonds
        natm = 0
        for i,j in bonds:
            natm = max( natm, max( i, j ) )
        natm += 1
    else:
        bond = qm3.utils.connectivity( molec )
        natm = molec.natm
    atmm = []
    conn = []
    for i in range( natm ):
        atmm.append( True )
        conn.append( [] )
    for i in sele_QM:
        atmm[i] = False
    for i,j in bond:
        conn[i].append( j )
        conn[j].append( i )
    latm = []
    excl = []
    nx12 = 0
    nx13 = 0
    nx14 = 0
    buf12 = ""
    buf13 = ""
    buf14 = ""
    for i in sele_QM:
        for j in conn[i]:
            if( j != i and atmm[j] ):
                latm.append( [ i, j ] )
                excl.append( [ i, j, 0.0 ] )
                nx12 += 1
                if( molec.type != [] ):
                    buf12 += "        # %s - %s || %s - %s\n"%( molec.labl[i], molec.labl[j], molec.type[i], molec.type[j] )
                else:
                    buf12 += "        # %s - %s\n"%( molec.labl[i], molec.labl[j] )
                buf12 += "        self.exc.append( qm3.engines.mmres.distance( kumb_kJ/mol.A^2, xref_A, [ %d, %d ] ) )\n"%( i, j )
                buf12 += "        self.exc[-1].ffac = 0.0\n"
                buf12 += "        self.exc[-1].gfac = [ 1.0, 0.0 ]\n"
#                buf12 += "        self.exc[-1].hind = [ self.sele.index( %d ), -1 ]\n"%( i )
            for k in conn[j]:
                if( k != i and atmm[k] ):
                    excl.append( [ i, k, 0.0 ] )
                    nx13 += 1
                    if( atmm[j] and atmm[k] ):
                        if( molec.type != [] ):
                            buf13 += "        # %s - %s - %s || %s - %s - %s\n"%( molec.labl[i], molec.labl[j],
                                molec.labl[k], molec.type[i], molec.type[j], molec.type[k] )
                        else:
                            buf13 += "        # %s - %s - %s\n"%( molec.labl[i], molec.labl[j], molec.labl[k] )
                        buf13 += "        self.exc.append( qm3.engines.mmres.angle( kumb_kJ/mol.rad^2, xref_deg, [ %d, %d, %d ] ) )\n"%( i, j, k )
                        buf13 += "        self.exc[-1].ffac = 0.0\n"
                        buf13 += "        self.exc[-1].gfac = [ 1.0, 0.0, 0.0 ]\n"
#                        buf13 += "        self.exc[-1].hind = [ self.sele.index( %d ), -1, -1 ]\n"%( i )
                for l in conn[k]:
                    if( k != i and l != j and l != i and atmm[l] ):
                        excl.append( [ i, l, 0.5 ] )
                        nx14 += 1
                        if( atmm[k] and atmm[l] ):
                            if( molec.type != [] ):
                                buf14 += "        # %s - %s - %s - %s || %s - %s - %s - %s\n"%( molec.labl[i], molec.labl[j], molec.labl[k],
                                        molec.labl[l], molec.type[i], molec.type[j], molec.type[k], molec.type[l] )
                            else:
                                buf14 += "        # %s - %s - %s - %s\n"%( molec.labl[i], molec.labl[j], molec.labl[k], molec.labl[l] )
                            buf14 += "        self.exc.append( qm3.engines.mmres.dihedral( { n: [ frc_kJ/mol, dsp_deg ] }, [ %d, %d, %d, %d ] ) )\n"%( i, j, k, l )
                            buf14 += "        self.exc[-1].ffac = 0.0\n"
                            if( atmm[j] ):
                                buf14 += "        self.exc[-1].gfac = [ 1.0, 0.0, 0.0, 0.0 ]\n"
#                                buf14 += "        self.exc[-1].hind = [ self.sele.index( %d ), -1, -1, -1 ]\n"%( i )
                            else:
                                buf14 += "        self.exc[-1].gfac = [ 1.0, 1.0, 0.0, 0.0 ]\n"
#                                buf14 += "        self.exc[-1].hind = [ self.sele.index( %d ), self.sele.index( %d ), -1, -1 ]\n"%( i, j )
    fd = open( "exclusions.src", "wt" )
    fd.write( "        self.exc = []\n" )
    fd.write( "        #------------------------------------------------------------------\n" )
    fd.write( buf12 )
    fd.write( "        #------------------------------------------------------------------\n" )
    fd.write( buf13 )
    fd.write( "        #------------------------------------------------------------------\n" )
    fd.write( buf14 )
    fd.write( """


    #------------------------------------------------------------------
    def exc_hess( self, sel, mol, dsp = 1.e-4 ):
        _func = mol.func
        _grad = mol.grad[:]
        _hess = mol.hess[:]
        zero  = [ 0.0 for i in range( 3 * mol.natm ) ]
        hess  = []
        for i in sel:
            i3 = i * 3
            for j in [0, 1, 2]:
                bak = mol.coor[i3+j]
                mol.coor[i3+j] = bak + dsp
                mol.grad = zero[:]
                for itm in self.exc:
                    itm.get_grad( mol )
                gp = mol.grad[:]
                mol.coor[i3+j] = bak - dsp
                mol.grad = zero[:]
                for itm in self.exc:
                    itm.get_grad( mol )
                mol.coor[i3+j] = bak
                hess.append( [] )
                for k in sel:
                    k3 = k * 3
                    for l in [0, 1, 2]:
                        hess[-1].append( ( gp[k3+l] - mol.grad[k3+l] ) / ( 2. * dsp ) )
        mol.func = _func
        mol.grad = _grad[:]
        for itm in self.exc:
            itm.get_grad( mol )
        mol.hess = []
        k = 0
        for i in range( self.size ):
            for j in range( self.size ):
                mol.hess.append( _hess[k] + 0.5 * ( hess[i][j] + hess[j][i] ) )
                k += 1
""" )
    fd.close()
    f = open( "sele_LA.pk", "wb" )
    pickle.dump( latm, f )
    f.close()
    excl.sort()
    f = open( "sele_EX.pk", "wb" )
    pickle.dump( excl, f )
    f.close()
    print( "\n>> %d exclusions generated (1-2:%d, 1-3:%d, 1-4:%d)"%( nx12 + nx13 + nx14, nx12, nx13, nx14 ) )
    return( latm )



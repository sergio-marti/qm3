# -*- coding: iso-8859-1 -*-
import sys
import os
import math
import qm3.constants
import qm3.engines

import  numpy
import  pyscf.gto
import  pyscf.dft
import  pyscf.qmmm
import  pyscf.grad


options = { "basis": "def2-svp",
            "conv_tol": 1.e-9,
            "charge": 0,
            "spin": 0,
            "method": "b3lypg",
            "memory": 4096, # MB
            "grid": 3,
            "max_cyc": 50 }


class run_native( object ):

    def __init__( self, mol, opts, sele, nbnd = [], link = [] ):
        self._ce = qm3.constants.H2J
        self._cx = 1.0 / qm3.constants.A0
        self._cg = self._ce * self._cx
        self.sel = sele[:]
        self.lnk = link[:]
        t = [ j for i,j in self.lnk ]
        self.nbn = [ i for i in nbnd if not i in t ]

        self.smb = mol.guess_symbols( sele )
        if( not mol.chrg ):
            mol.chrg = [ 0.0 for i in range( mol.natm ) ]

        aQM = pyscf.gto.Mole()
        aQM.unit = "Angstrom"
        aQM.symmetry = False
        aQM.basis = opts.pop( "basis" )
        aQM.spin = opts.pop( "spin" )
        aQM.charge = opts.pop( "charge" )
        aQM.atom = ""
        j = 0
        for i in sele:
            i3 = i * 3
            aQM.atom += "%-2s%20.10lf%20.10lf%20.10lf\n"%( self.smb[j], 
                mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) )
            j += 1
        self.vla = []
        if( len( self.lnk ) > 0 ):
            k = len( self.sel )
            for i,j in self.lnk:
                c, v = qm3.engines.LA_coordinates( i, j, mol )
                self.aQM.atom += "%-2s%20.10lf%20.10lf%20.10lf\n"%( "H", c[0], c[1], c[2] )
                self.vla.append( ( self.sel.index( i ), k, v[:] ) )
                k += 1
        aQM.build()
        if( aQM.spin == 0 ):
            self.dft = pyscf.dft.RKS( aQM )
        else:
            self.dft = pyscf.dft.UKS( aQM )
        self.dft.verbose = False
        self.dft.direct_scf = True
        self.dft.conv_tol = opts.pop( "conv_tol" )
        self.dft.max_cycle = opts.pop( "max_cyc" )
        self.dft.grids.level = opts.pop( "grid" )
        self.dft.xc = opts.pop( "method" )
        self.dft.max_memory = int( opts.pop( "memory" ) )
        if( len( self.nbn ) > 0 ):
            crd = []
            chg = []
            for i in self.nbn:
                i3 = i * 3
                chg.append( mol.chrg[i] )
                crd.append( [
                    mol.coor[i3]   - mol.boxl[0] * round( mol.coor[i3]   / mol.boxl[0], 0 ), 
                    mol.coor[i3+1] - mol.boxl[1] * round( mol.coor[i3+1] / mol.boxl[1], 0 ),
                    mol.coor[i3+2] - mol.boxl[2] * round( mol.coor[i3+2] / mol.boxl[2], 0 ) ] )
            self.scf = pyscf.qmmm.mm_charge( self.dft, numpy.array( crd ) * self._cx, numpy.array( chg ), unit = "Bohr" )
        else:
            self.scf = self.dft


    def update_coor( self, mol ):
        crd = []
        for i in self.sel:
            i3 = i * 3
            crd.append( [ ( mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) ) for j in [0, 1, 2] ] )
        self.vla = []
        if( len( self.lnk ) > 0 ):
            k = len( self.sel )
            for i,j in self.lnk:
                c, v = qm3.engines.LA_coordinates( i, j, mol )
                crd.append( c[:] )
                self.vla.append( ( self.sel.index( i ), k, v[:] ) )
                k += 1
        self.scf.mol.set_geom_( numpy.array( crd ) )
        if( len( self.nbn ) > 0 ):
            crd = []
            for i in self.nbn:
                i3 = i * 3
                crd.append( [ ( mol.coor[i3+j] - mol.boxl[j] * round( mol.coor[i3+j] / mol.boxl[j], 0 ) ) for j in [0, 1, 2] ] )
            self.scf.mm_mol.set_geom_( numpy.array( crd ) * self._cx )


    def get_func( self, mol ):
        self.update_coor( mol )
        mol.func += self.scf.kernel() * self._ce
        chg = self.scf.mulliken_pop( verbose = 0 )[1].tolist()
        for i in range( len( self.sel ) ):
            mol.chrg[self.sel[i]] = chg[i]


    def get_grad( self, mol ):
        self.update_coor( mol )
        mol.func += self.scf.kernel() * self._ce
        chg = self.scf.mulliken_pop( verbose = 0 )[1].tolist()
        for i in range( len( self.sel ) ):
            mol.chrg[self.sel[i]] = chg[i]
        g = self.scf.Gradients().run( grid_response = True ).de.flatten().tolist()
        qm3.engines.LA_gradient( self.vla, g )
        for i in range( len( self.sel ) ):
            i3 = 3 * i
            I3 = 3 * self.sel[i]
            for j in [0, 1, 2]:
                mol.grad[I3+j] += g[i3+j] * self._cg
        if( len( self.nbn ) > 0 ):
            den = self.scf.make_rdm1()
            dr  = self.scf.mol.atom_coords()[:,None,:] - self.scf.mm_mol.atom_coords()
            r   = numpy.linalg.norm( dr, axis = 2 )
            g   = numpy.einsum( "r,R,rRx,rR->Rx", self.scf.mol.atom_charges(), self.scf.mm_mol.atom_charges(), dr, r ** -3 )
            if( len( den.shape ) == 3 ):
                for i,q in enumerate( self.scf.mm_mol.atom_charges() ):
                    with self.scf.mol.with_rinv_origin( self.scf.mm_mol.atom_coord( i ) ):
                        v = self.scf.mol.intor( "int1e_iprinv" )
                    g[i] += ( numpy.einsum( "ij,xji->x", den[0], v ) + numpy.einsum( "ij,xij->x", den[0], v.conj() ) ) * -q
                    g[i] += ( numpy.einsum( "ij,xji->x", den[1], v ) + numpy.einsum( "ij,xij->x", den[1], v.conj() ) ) * -q
            else:
                for i,q in enumerate( self.scf.mm_mol.atom_charges() ):
                    with self.scf.mol.with_rinv_origin( self.scf.mm_mol.atom_coord( i ) ):
                        v = self.scf.mol.intor( "int1e_iprinv" )
                    g[i] += ( numpy.einsum( "ij,xji->x", den, v ) + numpy.einsum( "ij,xij->x", den, v.conj() ) ) * -q
            g = g.flatten().tolist()
            for i in range( len( self.nbn ) ):
                i3 = 3 * i
                I3 = 3 * self.nbn[i]
                for j in [0, 1, 2]:
                    mol.grad[I3+j] += g[i3+j] * self._cg

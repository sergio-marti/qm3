# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import    sys
if( sys.version_info[0] == 2 ):
    range = xrange
import    math
import    qm3.utils
import    qm3.constants
import    qm3.maths.matrix
import    qm3.fio
import    qm3.elements
import    re
import    collections


MXLAT = 1.0e300


class molecule( object ):

    def __init__( self, fname = None ):
        self.initialize()
        if( fname ): 
            self.pdb_read( fname )


    def initialize( self ):
        self.natm = 0
        # -- basic properties
        self.labl = []
        self.coor = []
        # -- sequence properties
        self.segn = []
        self.resi = []
        self.resn = []
        # -- atomic properties: empty if not available...
        self.anum = []
        self.chrg = []
        self.mass = []
        self.epsi = []
        self.rmin = []
        self.type = []
        # -- faster sequence stuff
        self.res_lim = []
        self.seg_lim = []
        self.indx = collections.OrderedDict()

        self.boxl = [ MXLAT, MXLAT, MXLAT ]


    def atom_number( self, segn, resi, labl ):
        # segname
        i_sgn = -1
        i = 0
        j = len( self.seg_lim )
        while( i < j and i_sgn == -1 ):
            if( self.segn[self.res_lim[self.seg_lim[i]]] == segn ):
                i_sgn = i
            i += 1
        if( i_sgn == -1 ):
            return( -1 )
        # resid
        i_rsi = -1
        i = self.seg_lim[i_sgn]
        j = self.seg_lim[i_sgn+1]
        while( i < j and i_rsi == -1 ):
            if( self.resi[self.res_lim[i]] == resi ):
                i_rsi = i
            i += 1
        if( i_rsi == -1 ):
            return( -1 )
        # label
        i_lbl = -1
        i = self.res_lim[i_rsi]
        j = self.res_lim[i_rsi+1]
        while( i < j and i_lbl == -1 ):
            if( self.labl[i] == labl ):
                i_lbl = i
            i += 1
        return( i_lbl )


    def norm_resid( self, segn = None ):
        for i in range( len( self.seg_lim ) - 1 ):
            if( self.segn[self.res_lim[self.seg_lim[i]]] == segn or segn == None ):
                l = 0
                for j in range( self.seg_lim[i], self.seg_lim[i+1] ):
                    l += 1
                    for k in range( self.res_lim[j], self.res_lim[j+1] ):
                        self.resi[k] = l
        return( self )


    def sph_sel( self, sele, radius ):
        crd = []
        cen = [ 0.0, 0.0, 0.0 ]
        for i in sele:
            i3 = i * 3
            crd += self.coor[i3:i3+3][:]
            for j in [ 0, 1, 2]:
                cen[j] += self.coor[i3+j]
        cen = [ i / float( len( sele ) ) for i in cen ]
        dsp = 0.0
        siz = len( crd ) // 3
        for k in range( siz ):
            dsp = max( dsp, qm3.utils.distanceSQ( cen, crd[3*k:3*k+3], self.boxl ) )
        cut = math.pow( radius + math.sqrt( dsp ) + 0.1, 2.0 )
        rad = radius * radius
        res = []
        for k0 in range( len( self.res_lim ) - 1 ):
            k1 = self.res_lim[k0]
            kn = self.res_lim[k0+1]
            kf = False
            while( k1 < kn and not kf ):
                kf |= qm3.utils.distanceSQ( cen, self.coor[3*k1:3*k1+3], self.boxl ) <= cut
                k1 += 1
            if( kf ):
                k1 = self.res_lim[k0]
                kf = False
                while( k1 < kn and not kf ):
                    i1 = 0
                    while( i1 < siz and not kf ):
                        kf |= qm3.utils.distanceSQ( self.coor[3*k1:3*k1+3], crd[3*i1:3*i1+3], self.boxl ) <= rad
                        i1 += 1
                    k1 += 1
                if( kf ):
                    res.append( k0 )
        out = []
        for i in res:
            for j in range( self.res_lim[i], self.res_lim[i+1] ):
                out.append( j )
        return( out )


    def rotate( self, center, axis, theta, sele = None ):
        if( sele ):
            t_sel = sele[:]
        else:
            t_sel = range( self.natm )
        t_crd = []
        for i in t_sel:
            i3 = i * 3
            t_crd += self.coor[i3:i3+3][:]
        qm3.utils.rotate( t_crd, center, axis, theta )
        for i in range( len( t_sel ) ):
            i3 = i * 3
            self.coor[3*t_sel[i]:3*t_sel[i]+3] = t_crd[i3:i3+3][:]


    def prune( self, sele ):
        out = molecule()
        out.natm = len( sele )
        for i in sele:
            i3 = i * 3
            out.segn.append( self.segn[i] )
            out.resn.append( self.resn[i] )
            out.resi.append( self.resi[i] )
            out.labl.append( self.labl[i] )
            out.coor += self.coor[i3:i3+3]
            if( self.anum ):
                out.anum.append( self.anum[i] )
            if( self.type ):
                out.type.append( self.type[i] )
            if( self.mass ):
                out.mass.append( self.mass[i] )
            if( self.chrg ):
                out.chrg.append( self.chrg[i] )
            if( self.epsi ):
                out.epsi.append( self.epsi[i] )
            if( self.rmin ):
                out.rmin.append( self.rmin[i] )
        out.settle()
        return( out )


    def append( self, molec ):
        self.natm += molec.natm
        self.segn += molec.segn[:]
        self.resn += molec.resn[:]
        self.resi += molec.resi[:]
        self.labl += molec.labl[:]
        self.coor += molec.coor[:]
        self.type += molec.type[:]
        self.anum += molec.anum[:]
        self.mass += molec.mass[:]
        self.chrg += molec.chrg[:]
        self.epsi += molec.epsi[:]
        self.rmin += molec.rmin[:]
        self.settle()


    def settle( self ):
        self.res_lim = []
        self.seg_lim = []
        l_seg = None
        l_rsn = None
        l_rsi = None
        for i in range( self.natm ):
            if( l_seg != self.segn[i] or l_rsn != self.resn[i] or l_rsi != self.resi[i] ):
                self.res_lim.append( i )
                l_seg = self.segn[i]
                l_rsi = self.resi[i]
                l_rsn = self.resn[i]
        l_seg = self.segn[0]
        self.seg_lim.append( 0 )
        for i in range( len( self.res_lim ) ):
            if( l_seg != self.segn[self.res_lim[i]] ):
                self.seg_lim.append( i )
                l_seg = self.segn[self.res_lim[i]]
        self.res_lim.append( self.natm )
        self.seg_lim.append( len( self.res_lim ) - 1 )
        ## -- atom indexing... (memory consuming, but faster atom numbers and selections...)
        self.indx = collections.OrderedDict()
        for i in self.seg_lim[:-1]:
            self.indx[self.segn[self.res_lim[i]]] = collections.OrderedDict()
        for i in range( len( self.res_lim ) - 1 ):
            self.indx[self.segn[self.res_lim[i]]][self.resi[self.res_lim[i]]] = collections.OrderedDict()
            for j in range( self.res_lim[i], self.res_lim[i+1] ):
                self.indx[self.segn[self.res_lim[i]]][self.resi[self.res_lim[i]]][self.labl[j]] = j


    # -- this only works if the system propagates on the main axes: orthorhombic box along X,Y,Z
    def guess_boxl( self ):
        m_box = self.coor[0:3][:]
        M_box = self.coor[0:3][:]
        for i in range( self.natm ):
            i3 = i * 3
            for j in [0, 1, 2]:
                m_box[j] = min( m_box[j], self.coor[i3+j] )
                M_box[j] = max( M_box[j], self.coor[i3+j] )
        self.boxl = [ M_box[j] - m_box[j] for j in [0, 1, 2] ]


    def guess_atomic_numbers( self ):
        if( self.mass  ):
            self.anum = []
            for i in range( self.natm ):
                self.anum.append( sorted( [ ( math.fabs( qm3.elements.mass[j] - self.mass[i] ), j ) for j in iter( qm3.elements.mass ) if j > 0 ] )[0][1] )
        else:
            self.mass = []
            self.anum = []
            for i in range( self.natm ):
                self.anum.append( qm3.elements.rsymbol[str.join( "", [ j for j in self.labl[i] if j.isalpha() ] ).title()] )
                self.mass.append( qm3.elements.mass[self.anum[i]] )


    def fill_masses( self ):
        if( self.anum ):
            self.mass = [ qm3.elements.mass[i] for i in self.anum ]
        else:
            self.guess_atomic_numbers()


    def guess_symbols( self, sele = None ):
        if( sele ):
            t_sel = sele[:]
        else:
            t_sel = range( self.natm )
        smb = []
        if( self.anum ):
            for i in t_sel:
                smb.append( qm3.elements.symbol[self.anum[i]] )
        elif( self.mass  ):
            for i in t_sel:
                smb.append( qm3.elements.symbol[ sorted( [ ( math.fabs( qm3.elements.mass[j] - self.mass[i] ), j ) for j in iter( qm3.elements.mass ) if j > 0 ] )[0][1] ] )
        else:
            for i in t_sel:
                smb.append( "".join( [ j for j in self.labl[i] if j.isalpha() ] ).title() )
        return( smb )


###################################################################################################
# PDB
#
    def pdb_read( self, fname = None, append = False ):
        """
          1         2         3         4         5         6         7
.123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789.1234
ATOM     31 HG22 ILE     3      -8.509  29.691  -4.228  0.00  0.00      A
ATOM    771  H1  TIP3  253       6.588   9.359  -8.787  1.00  0.00      B    H
ATOM  00000  OH2 HOH  3314     -11.039 -22.605  29.142  0.00  0.00      B4
ATOM  05681  SOD SOD   146     -37.840 -41.531   8.396  0.00  0.00      ION1

HETATM   84  OH2 HOH A  28      52.537   5.370  44.344                          
HETATM   85  H1  HOH A  29       8.127  45.914  57.300                          
HETATM   86  H2  HOH A  29       9.503  46.512  57.945                          
        """
        if( not append ):
            self.initialize()
        f = qm3.fio.open_r( fname )
        for l in f:
            if( l[0:4] == "ATOM" or l[0:4] == "HETA" ):
                self.labl.append( l[12:17].strip() )
                if( l[21] == " " ):
                    self.resn.append( l[17:22].strip() )
                    self.resi.append(  int( l[22:26] ) )
                    self.segn.append( l[72:].strip().split()[0] )
                else:
                    self.resn.append( l[17:21].strip() )
                    self.resi.append(  int( l[22:26] ) )
                    self.segn.append( l[21] )
                self.coor.append( float( l[30:38] ) )
                self.coor.append( float( l[38:46] ) )
                self.coor.append( float( l[46:54] ) )
                self.natm += 1
        qm3.fio.close( f, fname )
        self.settle()


    def pdb_write( self, fname = None, sele = None ):
        f = qm3.fio.open_w( fname )
        if( sele ):
            t_sel = sele[:]
        else:
            t_sel = range( self.natm )
        j = 0
        for i in t_sel:
            i3 = i * 3
            f.write( "ATOM  %5d %-5s%-5s%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf      %-4s\n"%( ( j % 99999 ) + 1, 
                " " * ( len( self.labl[i] ) < 4 ) + self.labl[i],
                self.resn[i], self.resi[i], self.coor[i3], self.coor[i3+1], self.coor[i3+2], 
                0.0, 0.0, self.segn[i] ) )
            j += 1
        f.write( "END\n" )
        qm3.fio.close( f, fname )


###################################################################################################
# DCD
#
    def dcd_read( self, dcd, advance = True ):
        f = False
        if( self.natm == dcd.N ):
            if( advance ):
                f = dcd.next()
            for i in range( dcd.N ):
                i3 = i * 3
                self.coor[i3]   = dcd.X[i]
                self.coor[i3+1] = dcd.Y[i]
                self.coor[i3+2] = dcd.Z[i]
        else:
            print( "- Wrong atom number (%d) vs DCD (%d)"%( self.natm, dcd.N ) )
        return( f )


    def dcd_write( self, dcd ):
        if( self.natm == dcd.N ):
            for i in range( dcd.N ):
                i3 = i * 3
                dcd.X[i] = self.coor[i3]
                dcd.Y[i] = self.coor[i3+1]
                dcd.Z[i] = self.coor[i3+2]
            dcd.append()
        else:
            print( "- Wrong atom number (%d) vs DCD (%d)"%( self.natm, dcd.N ) )


###################################################################################################
# XYZ
#
    def xyz_read( self, fname = None, append = False, replace = False ):
        if( not append and not replace ):
            self.initialize()
        f = qm3.fio.open_r( fname )
        n = int( f.readline().strip() )
        f.readline()
        if( replace and self.natm == n ):
            for i in range( n ):
                i3 = i * 3
                t = f.readline().split()
                for j in [0, 1, 2]:
                    self.coor[i3+j] = float( t[j+1] )
        else:
            for i in range( n ):
                t = f.readline().split()
                self.labl.append( t[0] )
                self.resi.append( 1 )
                self.resn.append( "XXX" )
                self.segn.append( "X" )
                self.coor += [ float( j ) for j in t[1:4] ]
                self.natm += 1
        qm3.fio.close( f, fname )
        self.settle()


    def xyz_write( self, fname = None, sele = None, formt = "%20.10lf" ):
        fmt = "%-4s" + 3 * formt + "\n"
        f = qm3.fio.open_w( fname )
        if( sele ):
            t_sel = sele[:]
        else:
            t_sel = list( range( self.natm ) )
        smb = self.guess_symbols( t_sel )
        siz = len( t_sel )
        f.write( "%d\n\n"%( siz ) )
        for i in range( siz ):
            f.write( fmt%( smb[i], self.coor[3*t_sel[i]], self.coor[3*t_sel[i]+1], self.coor[3*t_sel[i]+2] ) )
        qm3.fio.close( f, fname )


###################################################################################################
# ZMAT
#
    def zmat_read( self, fname = None, orig = [ 0.0, 0.0, 0.0 ], axis = [ 1.0, 0.0, 0.0 ] ):
        self.initialize()
        f = qm3.fio.open_r( fname )
        # 1st atom
        self.labl.append( f.readline().strip().upper() )
        self.coor += [ 0.0, 0.0, 0.0 ]
        self.natm += 1
        # 2nd atom
        t = f.readline().strip().split()
        self.labl.append( t[0].upper() )
        t_d = float( t[2] )
        self.coor += [ float( t[2] ), 0.0, 0.0 ]
        self.natm += 1
        # 3rd atom
        t = f.readline().strip().split()
        self.labl.append( t[0].upper() )
        t_d = float( t[2] )
        t_a = float( t[4] ) / qm3.constants.R2D
        if( t[1] == "1" ):
            self.coor += [ self.coor[3*0+0] + t_d * math.cos( t_a ), t_d * math.sin( t_a ), 0.0 ]
        else:
            self.coor += [ self.coor[3*1+0] - t_d * math.cos( t_a ), t_d * math.sin( t_a ), 0.0 ]
        self.natm += 1
        # 4th and so on...
        for l in f:
            t = l.strip().split()
            self.labl.append( t[0].upper() )
            n_a = int( t[1] ) - 1
            t_d = float( t[2] )
            n_b = int( t[3] ) - 1
            t_a = float( t[4] ) / qm3.constants.R2D
            n_c = int( t[5] ) - 1
            t_t = float( t[6] ) / qm3.constants.R2D
            # ------------------------------------------------------------------------
            cosa = math.cos( t_a )
            pa = self.coor[3*n_a:3*n_a+3][:]
            pb = self.coor[3*n_b:3*n_b+3][:]
            vb = [ i-j for i,j in zip( pb, pa ) ]
            r = 1. / math.sqrt( vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2] )
            if( math.fabs( cosa ) >= 0.9999999991 ):
                r *= ( cosa * t_d )
                self.coor += [ pa[i] + vb[i] + r for i in [0, 1, 2] ]
            else:
                pc = self.coor[3*n_c:3*n_c+3][:]
                va = [ i-j for i,j in zip( pc, pa ) ]
                xyb = math.sqrt( vb[0] * vb[0] + vb[1] * vb[1] )
                flg = 0
                if( xyb <= 0.10 ):
                    xpa = va[2]
                    va[2] = - va[0]
                    va[0] = xpa
                    xpb = vb[2]
                    vb[2] = - vb[0]
                    vb[0] = xpb
                    xyb = math.sqrt( vb[0] * vb[0] + vb[1] * vb[1] )
                    flg = 1
                costh = vb[0] / xyb
                sinth = vb[1] / xyb
                xpa = va[0] * costh + va[1] * sinth
                ypa = va[1] * costh - va[0] * sinth
                sinph = vb[2] * r
                cosph = math.sqrt( math.fabs( 1.0 - sinph * sinph ) )
                xqa = xpa * cosph + va[2] * sinph
                zqa = va[2] * cosph - xpa * sinph
                yza = math.sqrt( ypa * ypa + zqa * zqa )
                coskh = ypa / yza
                sinkh = zqa / yza
                if( yza < 1.0e-10 ):
                    coskh = 1.0
                    sinkh = 0.0
                sina =  math.sin( t_a )
                sind = -math.sin( t_t )
                cosd =  math.cos( t_t )
                vd = [ t_d * cosa, t_d * sina * cosd, t_d * sina * sind ]
                ypd = vd[1] * coskh - vd[2] * sinkh
                zpd = vd[2] * coskh + vd[1] * sinkh
                xpd = vd[0] * cosph - zpd * sinph
                zqd = zpd * cosph + vd[0] * sinph
                xqd = xpd * costh - ypd * sinth
                yqd = ypd * costh + xpd * sinth
                if( flg == 1 ):
                    xrd = -zqd
                    zqd = xqd
                    xqd = xrd
                self.coor += [ xqd + pa[0], yqd + pa[1], zqd + pa[2] ]
            # ------------------------------------------------------------------------
            self.natm += 1
        qm3.fio.close( f, fname )
        # ---------------------------------------------------------------------
        nax = qm3.maths.matrix.norm( axis )
        if( qm3.utils.distanceSQ( nax, [1.0, 0.0, 0.0] ) > 0.0 ):
            qm3.utils.superimpose_vector( nax, [1.0, 0.0, 0.0], self.coor )
            for i in range( self.natm ):
                i3 = i * 3
                for j in [0, 1, 2]:
                    self.coor[i3+j] += orig[j]
        # ---------------------------------------------------------------------
        self.resi = [ 1 for i in range( self.natm ) ]
        self.resn = [ "XXX" for i in range( self.natm ) ]
        self.segn = [ "X" for i in range( self.natm ) ]
        self.settle()


    # conn: [ [ idx_0, na_0, nb_0, nc_0 ], ..., [ idx_N, na_N, nb_N, nc_N ] ]        (C indexes)
    def zmat_write( self, conn, fname = None ):
        temp = [ conn[i][0] for i in range( len( conn ) ) ]
        rcon = []
        for i in range( len( conn ) ):
            rcon.append( [ conn[i][0] ] + [ temp.index( conn[i][j] ) + 1 for j in range( 1, len( conn[i] ) ) ] )
        f = qm3.fio.open_w( fname )
        # -------------
        f.write( "%-2s\n"%( self.labl[conn[0][0]][0] ) )
        f.write( "%-2s%6d%12.3lf\n"%( self.labl[conn[1][0]][0], 
            rcon[1][1], qm3.utils.distance( self.coor[3*conn[1][0]:3*conn[1][0]+3], self.coor[3*conn[1][1]:3*conn[1][1]+3] ) ) )
        f.write( "%-2s%6d%12.3lf%6d%12.3lf\n"%( self.labl[conn[2][0]][0], 
            rcon[2][1], qm3.utils.distance( self.coor[3*conn[2][0]:3*conn[2][0]+3], self.coor[3*conn[2][1]:3*conn[2][1]+3] ),
            rcon[2][2], qm3.utils.angle( self.coor[3*conn[2][0]:3*conn[2][0]+3], self.coor[3*conn[2][1]:3*conn[2][1]+3], self.coor[3*conn[2][2]:3*conn[2][2]+3] ) ) )
        for i in range( 3, self.natm ):
            f.write( "%-2s%6d%12.3lf%6d%12.3lf%6d%12.3lf\n"%( self.labl[conn[i][0]][0], 
                rcon[i][1], qm3.utils.distance( self.coor[3*conn[i][0]:3*conn[i][0]+3], self.coor[3*conn[i][1]:3*conn[i][1]+3] ),
                rcon[i][2], qm3.utils.angle( self.coor[3*conn[i][0]:3*conn[i][0]+3], self.coor[3*conn[i][1]:3*conn[i][1]+3], self.coor[3*conn[i][2]:3*conn[i][2]+3] ),
                rcon[i][3], qm3.utils.dihedral( self.coor[3*conn[i][0]:3*conn[i][0]+3], self.coor[3*conn[i][1]:3*conn[i][1]+3], self.coor[3*conn[i][2]:3*conn[i][2]+3], self.coor[3*conn[i][3]:3*conn[i][3]+3] ) ) )
        # -------------
        qm3.fio.close( f, fname )


###################################################################################################
# TOPOLOGY and PARAMETERS (x-plor PSF)
#
    def psf_read( self, fname ):
        out = True
        self.type = []
        self.chrg = []
        self.mass = []
        f = qm3.fio.open_r( fname )
        if( f.readline().split()[0] == "PSF" ):
            f.readline()
            for i in range( int( f.readline().split()[0] ) + 1 ):
                f.readline()
            if( self.natm == int( f.readline().split()[0] ) ):
                for i in range( self.natm ):
                    t = f.readline().split()
                    if( self.segn[i] == t[1] and self.resi[i] == int( t[2] ) and self.resn[i] == t[3] and self.labl[i] == t[4]  ):
                        self.type.append( t[5] )
                        self.chrg.append( float( t[6] ) )
                        self.mass.append( float( t[7] ) )
                    else:
                        self.type.append( None )
                        self.chrg.append( None )
                        self.mass.append( None )
                        print( "- Wrong PSF data (%d): %s/%s %d/%s %s/%s %s/%s"%( i+1, self.segn[i], t[1], self.resi[i], t[2], self.resn[i], t[3], self.labl[i], t[4] ) )
                        out = False
            else:
                print( "- Invalid number of atoms in PSF!" )
                out = False
        qm3.fio.close( f, fname )
        return( out )
    
    
    def nbnd_read( self, fname ):
        out = True
        if( self.type == [] ):
            print( "- PSF data undefined (load one!)" )
            return( False )
        self.epsi = []
        self.rmin = []
        nbd = {}
        f = qm3.fio.open_r( fname )
        for i,j,k in re.compile( "([A-Z0-9]+)[\ ]+([0-9\.]+)[\ ]+([0-9\.]+)" ).findall( f.read() ):
            nbd[i] = [ math.sqrt( float( j ) * qm3.constants.K2J ), float( k ) ]
        f.close()
        for i in range( self.natm ):
            if( self.type[i] in nbd ):
                self.epsi.append( nbd[self.type[i]][0] )
                self.rmin.append( nbd[self.type[i]][1] )
            else:
                self.epsi.append( None )
                self.rmin.append( None )
                print( "- Atom index %d misses Non-Bonded data..."%( i+1 ) )
                out = False
        qm3.fio.close( f, fname )
        return( out )

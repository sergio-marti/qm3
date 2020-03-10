#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function, division
import    sys
if( sys.version_info[0] == 2 ):
    range = xrange
import  re
try:
    import cStringIO as io
except:
    import io
try:
    import cPickle as pickle
except:
    import pickle

import  qm3.mol
import  qm3.elements
import  qm3.utils
import  qm3.engines


KEY = re.compile( "^>[Qq][Mm]3\:.+$" )
SP0 = re.compile( "^([0-9]+)$" )
SP1 = re.compile( "^([0-9]+)-([0-9]+)$" )
SP2 = re.compile( "^([A-Za-z0-9]+):([0-9]+)$" )
SP3 = re.compile( "^([A-Za-z0-9]+):([0-9]+)-([0-9]+)$" )
SP4 = re.compile( "^([A-Za-z0-9]+)/([0-9]+)/(.+)$" )
SP5 = re.compile( "^([A-Za-z0-9]+):([0-9]+)@([0-9\.]+)$" )


engines = {
    "gaussian": ( "qm3.engines.gaussian.gaussian", "/bin/bash r.gauss" ),
    "demon":    ( "qm3.engines.demon.demon", "/bin/bash r.demon" ),
    "orca":     ( "qm3.engines.orca.orca", "/bin/bash r.orca" ),
    "nwchem":   ( "qm3.engines.nwchem.nwchem", "/bin/bash r.nwchem" ),
    "sqm":      ( "qm3.engines.sqm.dl_sqm", "" ),
    "dftb":     ( "qm3.engines.dftb.dl_dftb", "" ),
    # a fake input must be provided for fDynamo ---------------------------------
    "fdynamo":  ( "qm3.engines.dynamo.py_dynamo", "./dynamo.so" ),
    # ---------------------------------------------------------------------------
    "charmm":   ( "qm3.engines.charmm.charmm_shm", "/bin/bash r.charmm &" ),
    "namd":     ( "qm3.engines.namd.namd_shm", "/bin/bash r.namd &" )
}
    

def __parse( fname ):
    dat = {}
    key = None
    buf = ""
    f = open( "input", "rt" )
    for l in f:
        if( KEY.match( l ) ):
            if( key != None and buf != None ):
                dat[key] = buf
            key = l.strip()[5:].upper()
            buf = ""
        else:
            buf += l
    f.close()
    if( key != None and buf != None ):
        dat[key] = buf
    return( dat )


def __sele( mol, buf ):
    sel  = [ False for i in range( mol.natm ) ]
    _not = False
    for itm in buf.split():
        # -- all
        if( itm == "*" ):
            for i in range( mol.natm ):
                sel[i] = True
        # -- negate ALL selection (at the end...)
        elif( itm == "not" ):
            _not = True
        # -- atom number (NOT C-indexing)
        elif( SP0.match( itm ) ):
            sel[int(itm)-1] = True
        # -- range of atom numbers (NOT C-indexing)
        elif( SP1.match( itm ) ):
            rng = [ int( i ) for i in SP1.findall( itm )[0] ]
            for i in range( rng[0], rng[1] + 1 ):
                sel[i-1] = True
        # -- residue number (by chain)
        elif( SP2.match( itm ) ):
            rng = SP2.findall( itm )[0]
            for i in list( mol.indx[rng[0]][int(rng[1])].values() ):
                sel[i] = True
        # -- range of residue numbers (by chain)
        elif( SP3.match( itm ) ):
            rng = SP3.findall( itm )[0]
            for i in range( int( rng[1] ), int( rng[2] ) + 1 ):
                for j in list( mol.indx[rng[0]][i].values() ):
                    sel[j] = True
        # -- chain / residue_number / atom_label
        elif( SP4.match( itm ) ):
            rng = SP4.findall( itm )[0]
            sel[mol.indx[rng[0]][int(rng[1])][rng[2]]] = True
        # -- radial selection by residue around chain / residue_number
        elif( SP5.match( itm ) ):
            rng = SP5.findall( itm )[0]
            for i in mol.sph_sel( list( mol.indx[rng[0]][int(rng[1])].values() ), float( rng[2] ) ):
                sel[i] = True
        # -- backbone
        elif( itm == "backbone" ):
            for i in range( mol.natm ):
                if( mol.labl[i] in [ "C", "O", "CA", "N" ] ):
                    sel[i] = True
    if( _not ):
        for i in range( mol.natm ):
            sel[i] = not sel[i]
    return( [ i for i in range( mol.natm ) if sel[i] ] )


def __mkinp( fname ):
    # -- parse input file
    dat = __parse( fname )
    if( dat == {} ):
        return
    # -- try to build a molecule...
    mol = None; xyz = None; zmt = None; pdb = None; psf = None; bnd = None
    if( "ZMT" in dat ):
        try:
            zmt = dat["ZMT"].strip()
            g = io.StringIO( zmt )
            mol = qm3.mol.molecule()
            mol.zmat_read( g )
            g.close()
            mol.guess_atomic_numbers()
            mol.boxl = []
        except:
            mol = None
    if( "XYZ" in dat ):
        try:
            xyz = dat["XYZ"].strip()
            xyz = "%d\n\n%s"%( len( xyz.split( "\n" ) ), xyz )
            g = io.StringIO( xyz )
            mol = qm3.mol.molecule()
            mol.xyz_read( g )
            g.close()
            mol.guess_atomic_numbers()
            mol.boxl = []
        except:
            mol = None
    if( "PDB" in dat ):
        tmp = dat["PDB"].strip().split()
        if( len( tmp ) >= 2 ):
            try:
                pdb = tmp[0]
                psf = tmp[1]
                mol = qm3.mol.molecule( pdb )
                mol.psf_read( psf )
                g = open( psf, "rt" )
                l = g.readline()
                while( l != "" and l.find( "!NBOND:" ) < 0 ):
                    l = g.readline()
                bnd = []
                n = int( l.split()[0] )
                i = 0
                while( i < n ):
                    t = [ int( j ) - 1 for j in g.readline().split() ]
                    m = len( t ) // 2
                    i += m
                    for j in range( m ):
                        bnd.append( [ t[2*j], t[2*j+1] ] )
                g.close()
                mol.boxl = []
                if( len( tmp ) == 3 ):
                    x = float( tmp[2] )
                    mol.boxl = [ x, x, x ]
                elif( len( tmp ) == 5 ):
                    mol.boxl = [ float( tmp[2] ), float( tmp[3] ), float( tmp[4] ) ]
            except:
                mol = None
    if( mol == None ):
        return
    # -- active selection: empty [] for ALL atoms
    sel = []
    if( "SELE" in dat ):
        sel = __sele( mol, dat["SELE"] )
        g = open( "sele.pk", "wb" )
        pickle.dump( sel, g )
        g.close()
    # -- search for MM (engine & input template)
    emm = None; imm = None
    if( "MM" in dat and "MMINP" in dat ):
        try:
            emm = engines[dat["MM"].split()[0]]
            imm = dat["MMINP"]
        except:
            emm = None
            imm = None
    # -- search for QM (engine, input, sele & envi): build exclusions
    eqm = None; iqm = None; sqm = []; smm = []
    if( "QM" in dat and "QMINP" in dat ):
        try:
            eqm = engines[dat["QM"].split()[0]]
            iqm = dat["QMINP"]
        except:
            eqm = None
            iqm = None
        if( eqm != None and iqm != None ):
            if( "QMSEL" in dat ):
                sqm = __sele( mol, dat["QMSEL"] )
            else:
                sqm = list( range( mol.natm ) )
            if( sqm != [] ):
                    qm3.engines.exclusions( sqm, mol, bnd )
            if( "QMENV" in dat ):
                smm = __sele( mol, dat["QMENV"] )
                if( smm != [] and sqm != [] ):
                    _tmp = []
                    g = open( "sele_LA.pk", "rb" )
                    for i,j in pickle.load( g ):
                        _tmp.append( j )
                    g.close()
                    smm = list( set( smm ).difference( set( sqm + _tmp ) ) )
        g = open( "sele_QM.pk", "wb" )
        pickle.dump( sqm, g )
        g.close()
        g = open( "sele_MM.pk", "wb" )
        pickle.dump( smm, g )
        g.close()
    if( eqm == None and emm == None ):
        return
    # -- let's go!
    fd = open( "run.py", "wt" )
    fd.write( """from __future__ import print_function, division
import sys
if( sys.version_info[0] == 2 ):
    range = xrange
import math
import os
import time
try:
    import cPickle as pickle
except:
    import pickle
try:
    import cStringIO as io
except:
    import io

import qm3.mol
import qm3.problem
import qm3.utils
import qm3.fio.dcd
""" )
    # -- import: engines
    if( eqm != None ):
        fd.write( "import %s\n"%( ".".join( eqm[0].split( "." )[0:-1] ) ) )
    if( emm != None ):
        fd.write( "import %s\n"%( ".".join( emm[0].split( "." )[0:-1] ) ) )
    if( eqm != None and emm != None ):
        fd.write( "import qm3.engines._mmint\n" )
    if( "REST" in dat ):
        fd.write( "import qm3.engines.mmres\n" )
    # -- import: actions
    if( "TODO" in dat ):
        tmp = dat["TODO"].upper()
        if( tmp.find( "FIRE" ) >= 0 or tmp.find( "BAKER" ) >= 0 ):
            fd.write( "import qm3.actions.minimize\n" )
        if( tmp.find( "PATH" ) >= 0 ):
            fd.write( "import qm3.actions.paths\n" )
        if( tmp.find( "DYNA" ) >= 0 ):
            fd.write( "import qm3.actions.dynamics\n" )
    # -- problem definition
    fd.write( """

class my_problem( qm3.problem.template ):
    def __init__( self ):
        qm3.problem.template.__init__( self )

""" )
    # -- molecule definition
    if( pdb != None and psf != None ):
        fd.write( """        self.mol = qm3.mol.molecule( "%s" )
        self.mol.psf_read( "%s" )
        self.mol.nbnd_read( "_NONBONDED_" )
"""%( pdb, psf ) )
    elif( xyz != None ):
        fd.write( """        f = io.StringIO( \"\"\"%s\"\"\" )
        self.mol = qm3.mol.molecule()
        self.mol.xyz_read( f )
        f.close()
        self.mol.guess_atomic_numbers()
        self.mol.fill_masses()
"""%( xyz ) )
    else:
        fd.write( """        f = io.StringIO( \"\"\"%s\"\"\" )
        self.mol = qm3.mol.molecule()
        self.mol.zmat_read( f )
        f.close()
        self.mol.guess_atomic_numbers()
        self.mol.fill_masses()
"""%( zmt ) )
    if( mol.boxl != [] ):
        fd.write( "        self.mol.boxl = [ %.4lf, %.4lf, %.4lf ]\n"%( mol.boxl[0], mol.boxl[1], mol.boxl[2] ) )
    # -- definition of the active system
    if( sel == [] ):
        fd.write( """
        self.sel  = []
        self.size = 3 * self.mol.natm
        self.coor = self.mol.coor
        self.mass = self.mol.mass
""" )
    else:
        fd.write( """
        f = open( "sele.pk", "rb" )
        self.sel = pickle.load( f )
        f.close()
        self.size = 3 * len( self.sel )
        self.coor = []
        self.mass = []
        for i in self.sel:
            i3 = i * 3
            self.coor += self.mol.coor[i3:i3+3][:]
            self.mass.append( self.mol.mass[i] )
""" )
    # -- generic
    fd.write( """
        self.func = 0.0
        self.grad = []
        self.hess = []
""" )
    # -- MM engine
    if( emm != None ):
        if( emm[0] == "qm3.engines.dynamo.py_dynamo" ):
            fd.write( "\n        self.emm = qm3.engines.dynamo.py_dynamo( \"%s\" )\n"%( emm[1] ) )
        elif( emm[0] == "qm3.engines.charmm.charmm_shm" ):
            fd.write( """
        os.system( "%s" )
        self.emm = qm3.engines.charmm.charmm_shm( "i.charmm" )
"""%( emm[1] ) )
        elif( emm[0] == "qm3.engines.namd.namd_shm" ):
            fd.write( """
        os.system( "%s" )
        while( not os.path.isfile( "namd.shmid" ) ):
            time.sleep( 1 )
        time.sleep( 2 )
        self.emm = qm3.engines.namd.namd_shm()
"""%( emm[1] ) )
    # -- QM engine
    if( eqm != None ):
        fd.write( """
        f = open( "sele_QM.pk", "rb" )
        sqm = pickle.load( f )
        f.close()
        f = open( "sele_LA.pk", "rb" )
        sla = pickle.load( f )
        f.close()
        f = open( "sele_MM.pk", "rb" )
        smm = pickle.load( f )
        f.close()
        f = io.StringIO( \"\"\"%s\"\"\" )
        self.eqm = %s( self.mol, f, sqm, smm, sla )
        f.close()
"""%( iqm, eqm[0] ) )
        if( eqm[1] != "" ):
            fd.write( "        self.eqm.exe = \"%s\"\n"%( eqm[1] ) )
    # -- QM(LJ) fixing
    if( eqm != None and emm != None ):
        fd.write( """
        f = open( "sele_EX.pk", "rb" )
        exc = pickle.load( f )
        f.close()
        self.fix = qm3.engines._mmint.QMLJ( self.mol, sqm, smm, exc )
""" )
    # -- restraints
    fd.write( "\n        self.umb = []\n" )
    if( "REST" in dat ):
        for itm in dat["REST"].strip().upper().split( "\n" ):
            tmp = itm.split()
            if( tmp[0][0:4] == "DIST" and len( tmp ) == 5 ):
                fd.write( "        self.umb.append( qm3.engines.mmres.distance( %.1lf, %.3lf, [ %d, %d ] ) )\n"%(
                    float( tmp[1] ), float( tmp[2] ), int( tmp[3] ) - 1, int( tmp[4] ) - 1 ) )
            elif( tmp[0][0:4] == "ANGL" and len( tmp ) == 6 ):
                fd.write( "        self.umb.append( qm3.engines.mmres.angle( %.1lf, %.1lf, [ %d, %d, %d ] ) )\n"%(
                    float( tmp[1] ), float( tmp[2] ), int( tmp[3] ) - 1, int( tmp[4] ) - 1, int( tmp[5] ) - 1 ) )
            elif( tmp[0][0:4] == "MULT" and len( tmp ) == 7 ):
                fd.write( "        self.umb.append( qm3.engines.mmres.multiple_distance( %.1lf, %.3lf, [ %d, %d, %d, %d ], [ 1.0, -1.0 ] ) )\n"%(
                    float( tmp[1] ), float( tmp[2] ), int( tmp[3] ) - 1, int( tmp[4] ) - 1, int( tmp[5] ) - 1, int( tmp[6] ) - 1 ) )
    # -- QM/MM MM-exclusions
    try:
        fd.write( "\n" )
        g = open( "exclusions.src", "rt" )
        fd.write( g.read() )
        g.close()
    except:
        fd.write( "\n        self.exc = []\n" )
    # -- zero QM charges on the MM engine
    if( sqm != None and emm != None ):
        fd.write( """
        for i in sqm:
            self.mol.chrg[i] = 0.0
        self.emm.update_chrg( self.mol )
""" )
    # -- DCD stuff
    fd.write( """
        self.dcd = qm3.fio.dcd.dcd()
""" )
    if( sel == [] ):
        fd.write( "        self.dcd.write( \"dcd\", self.mol.natm )\n" )
    else:
        fd.write( "        self.dcd.write( \"dcd\", self.mol.natm, self.sel )\n" )
    fd.write( """


    def current_step( self, stp ):
        if( stp % 10 == 0 ):
            self.mol.dcd_write( self.dcd )
""" )
    # -- update_coor method
    fd.write( """

    def update_coor( self ):
        if( self.sel != [] ):
            for i in range( len( self.sel ) ):
                ii = 3 * self.sel[i]
                jj = 3 * i
                for j in [ 0, 1, 2 ]:
                    self.mol.coor[ii+j] = self.coor[jj+j]
""" )
    # -- get_func method
    fd.write( """

    def get_func( self ):
        self.update_coor()
        self.mol.func = 0.0
""" )
    if( emm != None ):
        fd.write( "        self.emm.get_func( self.mol )\n" )
    if( eqm != None ):
        fd.write( "        self.eqm.get_func( self.mol )\n" )
    fd.write( "        for itm in self.umb:\n            itm.get_func( self.mol )\n" )
    fd.write( "        self.func = self.mol.func\n" )
    # -- get_grad method
    fd.write( """

    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( 3 * self.mol.natm ) ]
""" )
    if( emm != None ):
        fd.write( "        self.emm.get_grad( self.mol )\n" )
    if( eqm != None ):
        fd.write( "        self.eqm.get_grad( self.mol )\n" )
    if( emm != None and eqm != None ):
        fd.write( "        self.fix.get_grad( self.mol )\n" )
    fd.write( "        for itm in self.umb:\n            itm.get_grad( self.mol )\n" )
    fd.write( "        for itm in self.exc:\n            itm.get_grad( self.mol )\n" )
    fd.write( "        self.func = self.mol.func\n" )
    if( sel == [] ):
        fd.write( "        self.grad = self.mol.grad\n" )
    else:
        fd.write( """        self.grad = []
        for i in self.sel:
            i3 = i * 3
            self.grad += self.mol.grad[i3:i3+3][:]
        """ )
    # -- get_hess method (numerical)
    fd.write( """

    def get_hess( self ):
        if( os.access( "update.dump", os.R_OK ) ):
            self.get_grad()
            self.hess = [ 0.0 for i in range( self.size * self.size ) ]
            qm3.utils.manage_hessian( self.coor, self.grad, self.hess, should_update = True )
        else:
            self.update_coor()
            self.num_hess( central = True )
            qm3.utils.manage_hessian( self.coor, self.grad, self.hess, should_update = False )
""" )
    # -- setup MM input
    if( emm != None ):
        w = ""
        if( emm[0] == "qm3.engines.charmm.charmm_shm" ):
            w = "i.charmm"
        elif( emm[0] == "qm3.engines.namd.namd_shm" ):
            w = "i.namd"
        if( w != "" ):
            fd.write( """\n\n\n\n
f = open( "%s", "wt" )
f.write( \"\"\"%s
\"\"\" )
f.close()
"""%( w, imm ) )
    # -- instance object
    fd.write( """\n\n\n\n
obj = my_problem()
obj.get_grad()
print( obj.func )
import qm3.maths.matrix
qm3.maths.matrix.mprint( obj.grad, obj.mol.natm, 3 )
""" )
    # -- actions
    if( "TODO" in dat ):
        for itm in dat["TODO"].strip().upper().split( "\n" ):
            tmp = itm.split()
            #fire  step size nprt gtol
            if( tmp[0][0:4] == "FIRE" and len( tmp ) == 5 ):
                fd.write( """
qm3.actions.minimize.fire( obj, step_number = %d, step_size = %.4lf, print_frequency = %d, gradient_tolerance = %.2lf )
obj.mol.pdb_write( "last.pdb" )
"""%( int( tmp[1] ), float( tmp[2] ), int( tmp[3] ), float( tmp[4] ) ) )
            #baker step size nprt gtol [-1|mode]
            if( tmp[0][0:4] == "BAKE" and len( tmp ) >= 5 ):
                i = -1
                if( len( tmp ) == 6 ):
                    i = int( tmp[5] )
                fd.write( """
qm3.actions.minimize.baker( obj, step_number = %d, step_size = %.4lf, print_frequency = %d, gradient_tolerance = %.2lf, follow_mode = %d, allow_overlap = False )
obj.mol.pdb_write( "last.pdb" )
"""%( int( tmp[1] ), float( tmp[2] ), int( tmp[3] ), float( tmp[4] ), i ) )
            #path  step size nprt gtol
            if( tmp[0][0:4] == "PATH" and len( tmp ) == 5 ):
                fd.write( """
qm3.actions.paths.page_mciver( obj, step_number = %d, step_size = %.4lf, print_frequency = %d, gradient_tolerance = %.2lf, avoid_recrossing = False )
obj.mol.pdb_write( "last.pdb" )
"""%( int( tmp[1] ), float( tmp[2] ), int( tmp[3] ), float( tmp[4] ) ) )
            #dyna  step size nprt gamm temp
            if( tmp[0][0:4] == "DYNA" and len( tmp ) == 6 ):
                t = float( tmp[5] )
                fd.write( """
qm3.actions.dynamics.assign_velocities( obj, temperature = %.1lf )
qm3.actions.dynamics.langevin_verlet( obj, step_number = %d, step_size = %.4lf, print_frequency = %d, gamma_factor = %.1lf, temperature = %.1lf )
obj.mol.pdb_write( "last.pdb" )
"""%( t, int( tmp[1] ), float( tmp[2] ), int( tmp[3] ), float( tmp[4] ), t ) )
            #freq  temp [n_atm:mass ...]
            if( tmp[0][0:4] == "FREQ" and len( tmp ) >= 2 ):
                if( len( tmp ) > 2 ):
                    fd.write( "\n" )
                    for itm in tmp[2:]:
                        _t = itm.split( ":" )
                        _e = int( _t[0] ) - 1
                        _m = float( _t[1] )
                        if( _e >= 0 and _e < mol.natm ):
                            fd.write( "obj.mass[%d] = %.4lf\n"%( _e, _m ) )
                fd.write( """
obj.get_hess()
frq, mds = qm3.utils.hessian_frequencies( obj.mass, obj.coor, obj.hess, True )
print( "Frequencies (cm^-1):" )
print( frq )
zpe, gib = qm3.utils.gibbs_rrho( obj.mass, obj.coor, frq, temp = %.2lf )
print( "ZPE:", zpe )
print( "Gibbs (1atm):", gib )
print( "Total:", zpe + gib, "_kJ/mol" )
if( obj.sel == [] ):
    s = obj.mol.guess_symbols()
else:
    s = obj.mol.guess_symbols( obj.sel )
i = 0
while( i < obj.size and math.fabs( frq[i] ) < 10. ):
    i += 1
if( frq[i] < .0 ):
    qm3.utils.normal_mode_view( obj.coor, frq, mds, s, i, afac = 4 )
else:
    qm3.utils.normal_mode_view( obj.coor, frq, mds, s, i )
"""%( float( tmp[1] ) ) )
    # -- stop engines (if needed...)
    fd.write( "\nobj.dcd.close()\n" )
    if( emm != None and emm in [ "qm3.engines.charmm.charmm_shm", "qm3.engines.namd.namd_shm" ] ):
        fd.write( "\nobj.emm.stop()\n" )
    # -- all done
    fd.close()
        





"""
>QM3:PDB
    pdb file path
    psf file path
    [X_box [Y_box Z_box]]

>QM3:XYZ
    cartesian coordinates (for QM)

>QM3:ZMT
    z-matrix coordinates (for QM)

>QM3:QM
   QM program
        "gaussian" == qm3.engines.gaussian.gaussian

>QM3:QMsel
    QM atoms selection

>QM3:QMenv
    charges around QM atoms

>QM3:QMinp
    QM input template

>QM3:MM
    MM program
        "namd" == qm3.enginges.namd.namd_shm

>QM3:MMinp
    MM input template

>QM3:Sele
    active atoms selection

>QM3:Rest
    distance kumb ref atm_1 atm_2
    angle kumb ref atm_1 atm_2 amt_3
    multiple_distance kumb ref atm_1 atm_2 atm_3 atm_4  [1,2 - 3,4]

>QM3:Todo
    fire  step size nprt gtol
    baker step size nprt gtol [-1|mode]
    path  step size nprt gtol
    dyna  step size nprt gamm temp
    freq  temp [n_atm:mass ...]


QM3:QM | QM3:MM     compulsory one of them
QM3:XXinp           compulsory if QM3:XX defined
QM3:Sele            all the atoms by default
QM3:PDB             box size optional (X: cubic | X Y Z: orthorhombic)
                        - overrides: QM3:XYZ and QM3:ZMT
                        - parsing: QM3:QMsel && QM3:QMenv
"""

__mkinp( sys.argv[1] )

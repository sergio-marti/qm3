#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-
import  sys
import  bottle
import  os
import  re
import  io
import  pickle
import  qm3.mol
import  qm3.elements
import  qm3.utils
import  qm3.engines
import  qm3.fio.xplor
import  qm3.engines.sander


QM_engines = {
	"SQM": "qm3.engines.sqm.run_dynlib",
	"xTB": "qm3.engines.xtb.run_dynlib",
	"MOPAC": "qm3.engines.mopac.run_dynlib",
	"PySCF": "qm3.engines.pyscf.run_native",
	"Gaussian": "qm3.engines.gaussian.run_single",
	"NWChem": "qm3.engines.nwchem.run_single",
	"Orca": "qm3.engines.orca.run_single" }


MM_engines = {
	"OpenMM": "qm3.engines.openmm.run_native",
	"NAMD": "qm3.engines.namd.run_shmem" }



cwd = os.path.dirname( sys.argv[0] )
pdb = "None"
knd = None
top = None
if( len( sys.argv ) >= 2 ):
	pdb = sys.argv[1]
if( len( sys.argv ) >= 3 ):
    tmp = sys.argv[2].lower().split( ":" )
    if( tmp[0] in [ "charmm", "amber" ] ):
        knd = tmp[0]
        top = tmp[1]
print( pdb, knd, top )



@bottle.route( "/" )
def __index():
    qm = ""
    for ee in [ j for i,j in sorted( [ k.lower(), k ] for k in QM_engines ) ]:
        qm += "<option>%s</option>"%( ee )
    mm = ""
    for ee in [ j for i,j in sorted( [ k.lower(), k ] for k in MM_engines ) ]:
        mm += "<option>%s</option>"%( ee )
    pdb_view = "None"
    if( pdb != "None" ):
        pdb_view = "<a href=\"/loadpdb\" target=\"_blank\">%s</a>"%( os.path.basename( pdb ) )
    qm_sel = ""
    qm_cor = "none"
    qm_chg = "block"
    qm_env = "block"
    if( pdb == "None" ):
        qm_sel = "*"
        qm_cor = "block"
        qm_chg = "none"
        qm_env = "none"
    return( bottle.template( "index.html", 
        pdb = pdb_view,
        qm_eng = qm,
        qm_sel = qm_sel,
        qm_chg = qm_chg,
        qm_cor = qm_cor,
        qm_env = qm_env,
        mm_eng = mm ) )


@bottle.route( "/jsmol/<sub:path>" )
def __jsmol( sub ):
    if( sub[-3:] == "png" ):
        return( bottle.static_file( sub, root = os.path.join( cwd, "jsmol" ), mimetype = "image/png" ) )
    else:
        return( bottle.static_file( sub, root = os.path.join( cwd, "jsmol" ) ) )


@bottle.route( "/editor" )
def __editor():
    return( bottle.template( "editor.html" ) )


@bottle.route( "/loadpdb" )
def __loadpdb():
    return( bottle.template( "loadpdb.html" ) )


@bottle.route( "/flushpdb" )
def __flushpdb():
    if( os.path.isfile( pdb ) ):
        job = os.path.basename( pdb )
        fld = os.path.dirname( pdb )
        return( bottle.static_file( job, root = fld, mimetype = "application/octet-stream" ) )


SP0 = re.compile( "^([0-9]+)$" )
SP1 = re.compile( "^([0-9]+)-([0-9]+)$" )
SP2 = re.compile( "^([A-Za-z0-9]+):([0-9]+)$" )
SP3 = re.compile( "^([A-Za-z0-9]+):([0-9]+)-([0-9]+)$" )
SP4 = re.compile( "^([A-Za-z0-9]+)/([0-9]+)/(.+)$" )
SP5 = re.compile( "^([A-Za-z0-9]+):([0-9]+)@([0-9\.]+)$" )
def __selection( mol, buf ):
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


@bottle.route( "/mkinput", method = "POST" )
def __mkinput():
    job = "run_qm3"
    # -- PDB coordinates
    if( pdb != "None" ):
        if( knd == "charmm" and top != None ):
            mol = qm3.mol.molecule( pdb )
            qm3.fio.xplor.psf_read( mol, top )
            mol.guess_atomic_numbers()
        elif( knd == "amber" and top != None ):
            mol = qm3.mol.molecule( pdb )
            qm3.engines.sander.topology_read( mol, top )
        else:
            mol = qm3.mol.molecule( pdb )
    else:
        mol = None
    # -- QM cartesian coordinates
    act = bottle.request.forms.get( "Todo" ).strip()
    xyz = bottle.request.forms.get( "QMcor" ).strip()
    if( mol == None and xyz != "" ):
        try:
            xyz = "%d\n\n%s"%( len( xyz.split( "\n" ) ), xyz.replace( "\r", "" ) )
            mol = qm3.mol.molecule()
            f = io.StringIO( xyz )
            f.seek( 0 )
            mol.xyz_read( f )
            f.close()
            mol.guess_atomic_numbers()
        except:
            xyz = ""
            mol = None
    # -- (pickle) selections: Selec
    sel = []
    _sel = bottle.request.forms.get( "Selec" ).strip()
    if( _sel != "" ):
        if( mol != None ):
            sel = __selection( mol, _sel )
            if( len( sel ) == mol.natm ):
                sel = []
    f = open( "sele.pk", "wb" )
    pickle.dump( sel, f )
    f.close()
    # -- (pickle) selections: QMsel
    eqm = bottle.request.forms.get( "QMeng" ).strip()
    f = open( "sele_LA.pk", "wb" )
    pickle.dump( [], f )
    f.close()
    f = open( "sele_EX.pk", "wb" )
    pickle.dump( [], f )
    f.close()
    sqm = []
    _sqm = bottle.request.forms.get( "QMsel" ).strip()
    if( _sqm != "" and eqm != "" ):
        if( mol != None ):
            sqm = __selection( mol, _sqm )
            # -- LAs and exclusions
            if( sqm != [] ):
                qm3.engines.exclusions( sqm, molec = mol )
    f = open( "sele_QM.pk", "wb" )
    pickle.dump( sqm, f )
    f.close()
    # -- (pickle) selections: QMenv
    smm = []
    _smm = bottle.request.forms.get( "QMenv" ).strip()
    if( _smm != "" and eqm != "" ):
        if( mol != None ):
            smm = __selection( mol, _smm )
    if( smm != [] and sqm != [] ):
        _tmp = []
        f = open( "sele_LA.pk", "rb" )
        for i,j in pickle.load( f ):
            _tmp.append( j )
        f.close()
        smm = list( set( smm ).difference( set( sqm + _tmp ) ) )
    f = open( "sele_MM.pk", "wb" )
    pickle.dump( smm, f )
    f.close()
    if( mol != None and sel != [] ):
        tmp = set( list( range( mol.natm ) ) ).difference( set( sel ) )
        if( sqm != [] ):
            tmp = tmp.union( set( sqm ) )
        f = open( "fixed.pk", "wb" )
        pickle.dump( list( tmp ), f )
        f.close()
    # -- let's go!
    f = open( job, "wt" )
    f.write( """#!/usr/bin/env python3
import  sys
import  math
import  os
import  time
import  pickle
import  io
import  qm3.mol
import  qm3.problem
import  qm3.utils
""" )
    # -- imports
    if( eqm != "--" ):
        f.write( "import  %s\n"%( ".".join( QM_engines[eqm].split( "." )[0:-1] ) ) )
    emm = bottle.request.forms.get( "MMeng" ).strip()
    if( emm != "--" ):
        f.write( "import  %s\n"%( ".".join( MM_engines[emm].split( "." )[0:-1] ) ) )
        if( emm == "OpenMM" ):
            f.write( "import  openmm\nimport  openmm.app\n" )
    # -- RR_[0-9]+_[t/{i,j,k,l}/r/u]
    umb = re.compile( "RR_[0-9]+_t" ).findall( " ".join( list( bottle.request.forms.keys() ) ) )
    if( len( umb ) > 0 ):
        f.write( "import  qm3.engines.mmres\n" )
    # -- kind of forcefield
    if( knd == "charmm" ):
        f.write( "import  qm3.fio.xplor\n" )
    else:
        f.write( "import  qm3.engines.sander\n" )
    # -- use qm3.engines.mmint for NAMD
    if( emm == "NAMD" ):
        f.write( "import  qm3.engines.mmint\n" )
    # -- actions
    if( act == "min_fire" ):
        f.write( "import  qm3.actions.minimize\n" )
    elif( act == "dyn_lang" ):
        f.write( "import  qm3.actions.dynamics\n" )
    # -- molecule definition (based on pdb or cartesian input: QMcor)
    f.write( """

class my_problem( qm3.problem.template ):
    def __init__( self ):
        qm3.problem.template.__init__( self )

""" )
    if( xyz != "" ):
        f.write( """        f = io.StringIO( \"\"\"%s\"\"\" )
        f.seek( 0 )
        self.mole = qm3.mol.molecule()
        self.mole.xyz_read( f )
        f.close()
        self.mole.guess_atomic_numbers()
        self.mole.fill_masses()
"""%( xyz ) )
    else:
        if( mol != None and pdb != "None" ):
            f.write( "        self.mole = qm3.mol.molecule( \"%s\" )\n"%( pdb ) )
        else:
            f.write( "        self.mole = qm3.mol.molecule( _PDB_ )\n" )
        # -- get MM info from PSF or PRMTOP
        if( knd == "charmm" ):
            f.write( """
        qm3.fio.xplor.psf_read( self.mole, _PSF_ )
        self.mole.guess_atomic_numbers()
""" )
        else:
            f.write( """
        qm3.engines.sander.topology_read( self.mole, _PRMTOP_ )
""" )
    # -- boxl: box_X, box_Y, box_Z
    tmp = bottle.request.forms.get( "box_X" ).strip()
    try:
        box = [ float( tmp ), float( tmp ), float( tmp ) ]
    except:
        box = []
    try:
        tmp = float( bottle.request.forms.get( "box_Y" ).strip() )
        box[1] = tmp
    except:
        pass
    try:
        tmp = float( bottle.request.forms.get( "box_Z" ).strip() )
        box[2] = tmp
    except:
        pass
    if( box != [] ):
        f.write( "        self.mole.boxl = [ %.2lf, %.2lf, %.2lf ]\n"%( box[0], box[1], box[2] ) )
    # -- definition of the active system (full or sele.pk based)
    if( sel == [] ):
        f.write( """
        self.size = 3 * self.mole.natm
        self.coor = self.mole.coor
        self.mass = self.mole.mass
        self.sele = []
""" )
    else:
        f.write( """
        f = open( "sele.pk", "rb" )
        self.sele = pickle.load( f )
        f.close()
        self.size = 3 * len( self.sele )
        self.coor = []
        self.mass = []
        for i in self.sele:
            i3 = i * 3
            self.coor += self.mole.coor[i3:i3+3][:]
            self.mass.append( self.mole.mass[i] )
""" )
    # -- generic
    f.write( """
        self.func = 0.0
        self.grad = []
        self.hess = []
""" )
    # -- engines: MMeng 
    if( emm == "NAMD" ):
        if( MM_engines[emm] == "qm3.engines.namd.run_shmem" ):
            f.write( """
        os.system( "bash ./r.namd &" )
        while( not os.path.isfile( "namd.shmid" ) ):
            time.sleep( 1 )
        time.sleep( 4 )
        self.emm = %s()
"""%( MM_engines[emm] ) )
    if( eqm != "--" ):
        f.write( """
        f = open( "sele_QM.pk", "rb" )
        sqm = pickle.load( f )
        f.close()
        f = open( "sele_LA.pk", "rb" )
        sla = pickle.load( f )
        f.close()
        f = open( "sele_MM.pk", "rb" )
        smm = pickle.load( f )
        f.close()
""" )
    if( emm == "OpenMM" ):
        if( knd == "charmm" ):
            f.write( """
        _top = openmm.app.charmmpsffile.CharmmPsfFile( _PSF_ )
        _sys = _top.createSystem(
            openmm.app.charmmparameterset.CharmmParameterSet( _PRM_ ),
            nonbondedMethod = openmm.app.CutoffPeriodic,
            nonbondedCutoff = _CUTOFF_ * openmm.unit.angstrom,
            rigidWater = False,
            implicitSolvent = None,
            switchDistance = _SWITCH_ * openmm.unit.angstrom )
""" )
        else:
            f.write( """
        _top = openmm.app.amberprmtopfile.AmberPrmtopFile( _PRMTOP_ )
        _sys = _top.createSystem(
            nonbondedMethod = openmm.app.CutoffPeriodic,
            nonbondedCutoff = _CUTOFF_ * openmm.unit.angstrom,
            rigidWater = False,
            implicitSolvent = None,
            switchDistance = _SWITCH_ * openmm.unit.angstrom )
""" )
        f.write( """
        _top.topology.setUnitCellDimensions(
            openmm.Vec3( self.mole.boxl[0], self.mole.boxl[1], self.mole.boxl[2] ) * openmm.unit.angstrom )

""" )
        if( sqm != [] ):
            f.write( "        self.emm = qm3.engines.openmm.run_native( _sys, _top.topology, qm_atom = sqm, platform = _PLATFORM_ )\n\n" )
        else:
            f.write( "        self.emm = qm3.engines.openmm.run_native( _sys, _top.topology, platform = _PLATFORM_ )\n\n" )
    # -- engines: QMeng
    if( eqm != "--" ):
        if( QM_engines[eqm] in [ "qm3.engines.gaussian.run_single", "qm3.engines.orca.run_single", "qm3.engines.nwchem.run_single" ] ):
            _w = QM_engines[eqm].split( "." )[-2]
            f.write( """
        self.eqm = %s( self.mole, "i.%s", sqm, smm, sla )
        self.eqm.exe = "bash ./r.%s"
"""%( QM_engines[eqm], _w, _w ) )
        elif( QM_engines[eqm] == "qm3.engines.sqm.run_dynlib" ):
            f.write( "        self.eqm = %s( self.mole, \"i.sqm\", sqm, smm, sla )\n"%( QM_engines[eqm] ) )
        elif( QM_engines[eqm] == "qm3.engines.xtb.run_dynlib" ):
            f.write( "        self.eqm = %s( self.mole, _CHARGE_, _UNPAIRED_, sqm, smm, sla )\n"%( QM_engines[eqm] ) )
        elif( QM_engines[eqm] == "qm3.engines.mopac.run_dynlib" ):
            f.write( "        self.eqm = %s( self.mole, _METHOD_, _CHARGE_, _MULT_, sqm, smm, sla )\n"%( QM_engines[eqm] ) )
        elif( QM_engines[eqm] == "qm3.engines.pyscf.run_native" ):
            f.write( "        opt = qm3.engines.pyscf.options\n" )
            f.write( "        self.eqm = %s( self.mole, opt, sqm, smm, sla )\n"%( QM_engines[eqm] ) )
    # -- QM(LJ) fixing
    if( emm not in [ "--", "OpenMM" ] and sqm != [] and smm != [] ):
        if( knd == "charmm" ):
            f.write( "\n        qm3.engines.mmint.non_bonded( self.mole, [ _PARM_ ], qm3.engines.mmint.p_charmm )" )
        else:
            f.write( "\n        qm3.engines.mmint.non_bonded( self.mole, [ _PARM_ ], qm3.engines.mmint.p_sander )" )
        f.write( """
        f = open( "sele_EX.pk", "rb" )
        exc = pickle.load( f )
        f.close()
        self.fix = qm3.engines.mmint.QMLJ( self.mole, sqm, smm, exc )
""" )
    # exclusions [empty by default: suggest the user to use "exclusions.src" after loading missing parameters]
    f.write( "\n        # -- parse the generated file \"exclusions.src\" to address QM/MM boundaries properly\n" )
    f.write( "        self.exc = []\n" )
    # -- engines: restraints
    if( len( umb ) > 0 ):
        f.write( "\n        self.umb = []\n" )
        for itm in umb:
            try:
                _t = bottle.request.forms.get( itm ).strip()
                _u = float( bottle.request.forms.get( itm[0:-1] + "u" ).strip() )
                _r = float( bottle.request.forms.get( itm[0:-1] + "r" ).strip() )
                _i = int( bottle.request.forms.get( itm[0:-1] + "i" ).strip() ) - 1 
                _j = int( bottle.request.forms.get( itm[0:-1] + "j" ).strip() ) - 1
            except:
                _t = "__WRONG__"
            if( _t == "dst" ):
                f.write( "        self.umb.append( qm3.engines.restraints.distance( %.1lf, %.3lf, [ %d, %d ] ) )\n"%( _u, _r, _i, _j ) )
            elif( _t == "ang" ):
                try:
                    _k = int( bottle.request.forms.get( itm[0:-1] + "k" ).strip() ) - 1
                    f.write( "        self.umb.append( qm3.engines.restraints.angle( %.1lf, %.2lf, [ %d, %d, %d ] ) )\n"%( _u, _r, _i, _j, _k ) )
                except:
                    pass
            elif( _t == "mul" ):
                try:
                    _k = int( bottle.request.forms.get( itm[0:-1] + "k" ).strip() ) - 1
                    _l = int( bottle.request.forms.get( itm[0:-1] + "l" ).strip() ) - 1
                    f.write( "        self.umb.append( qm3.engines.restraints.multiple_distance( %.1lf, %.3lf, [ %d, %d, %d, %d ], [ 1.0, -1.0 ] ) )\n"%( _u, _r, _i, _j, _k, _l ) )
                except:
                    pass
    # -- zero QM charges on the MM engine (if any...)
    if( sqm != [] and emm != "--" ):
        f.write( """
        for i in sqm:
            self.mole.chrg[i] = 0.0
        self.emm.update_chrg( self.mole )
""" )
    # -- update_coor method
    f.write( """

    def update_coor( self ):
        if( self.sele != [] ):
            for i in range( len( self.sele ) ):
                ii = 3 * self.sele[i]
                jj = 3 * i
                for j in [ 0, 1, 2 ]:
                    self.mole.coor[ii+j] = self.coor[jj+j]
""" )
    # -- get_func method
    f.write( """

    def get_func( self ):
        self.update_coor()
        self.mole.func = 0.0
""" )
    if( emm != "--" ):
        f.write( "        self.emm.get_func( self.mole )\n" )
    if( sqm != [] ):
        f.write( "        self.eqm.get_func( self.mole )\n" )
    if( emm not in [ "--", "OpenMM" ] and sqm != [] and smm != [] ):
        f.write( "        self.fix.get_func( self.mole )\n" )
    f.write( "        for exc in self.exc:\n            exc.get_func( self.mole )\n" )
    if( len( umb ) > 0 ):
        f.write( "        for umb in self.umb:\n            umb.get_func( self.mole )\n" )
    f.write( "        self.func = self.mole.func\n" )
    # -- get_grad method
    f.write( """

    def get_grad( self ):
        self.update_coor()
        self.mole.func = 0.0
        self.mole.grad = [ 0.0 for i in range( 3 * self.mole.natm ) ]
""" )
    if( emm != "--" ):
        f.write( "        self.emm.get_grad( self.mole )\n" )
    if( sqm != [] ):
        f.write( "        self.eqm.get_grad( self.mole )\n" )
    if( emm not in [ "--", "OpenMM" ] and sqm != [] and smm != [] ):
        f.write( "        self.fix.get_grad( self.mole )\n" )
    f.write( "        for exc in self.exc:\n            exc.get_grad( self.mole )\n" )
    if( len( umb ) > 0 ):
        f.write( "        for umb in self.umb:\n            umb.get_grad( self.mole )\n" )
    f.write( "        self.func = self.mole.func\n" )
    if( sel == [] ):
        f.write( "        self.grad = self.mole.grad\n" )
    else:
        f.write( """        self.grad = []
            for i in self.sele:
                i3 = i * 3
                self.grad += self.mole.grad[i3:i3+3][:]
        """ )
    # -- get_hess method (numerical)
    f.write( """

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
    # -- object generation
    f.write( """




obj = my_problem()
""" )
    # -- actions
    if( act == "min_fire" ):
        try:
            _stp =   int( bottle.request.forms.get( "fire_stp" ).strip() )
            _siz = float( bottle.request.forms.get( "fire_siz" ).strip() )
            _prt =   int( bottle.request.forms.get( "fire_prt" ).strip() )
            _tol = float( bottle.request.forms.get( "fire_tol" ).strip() )
        except:
            _stp = 1000
            _siz = 0.1
            _prt = 1
            _tol = 1.0
        f.write( """
qm3.actions.minimize.fire( obj, step_number = %d, step_size = %.4lf, print_frequency = %d, gradient_tolerance = %.2lf )
obj.mole.pdb_write( "last.pdb" )
"""%( _stp, _siz, _prt, _tol ) )
    elif( act == "dyn_lang" ):
        try:
            _stp =   int( bottle.request.forms.get( "lang_stp" ).strip() )
            _siz = float( bottle.request.forms.get( "lang_siz" ).strip() )
            _prt =   int( bottle.request.forms.get( "lang_prt" ).strip() )
            _gam = float( bottle.request.forms.get( "lang_gam" ).strip() )
            _tmp = float( bottle.request.forms.get( "lang_tmp" ).strip() )
        except:
            _stp = 1000
            _siz = 0.001
            _prt = 1
            _gam = 50.0
            _tmp = 300.0
        f.write( """
qm3.actions.dynamics.assign_velocities( obj, temperature = %.1lf )
qm3.actions.dynamics.langevin_verlet( obj, step_number = %d, step_size = %.4lf, print_frequency = %d, gamma_factor = %.1lf, temperature = %.1lf )
obj.mole.pdb_write( "last.pdb" )
"""%( _tmp, _stp, _siz, _prt, _gam, _tmp ) )
    elif( act == "nor_mode" and emm == "--" and eqm != "--" ):
        try:
            _tmp = float( bottle.request.forms.get( "mode_tmp" ).strip() )
            _pre = float( bottle.request.forms.get( "mode_pre" ).strip() )
        except:
            _tmp = 300.0
            _pre = 1.0
        # -- user masses
        for itm in bottle.request.forms.get( "mode_mas" ).strip().split():
            try:
                _t = itm.split( ":" )
                _e = int( _t[0] ) - 1
                _m = float( _t[1] )
                if( _e >= 0 and _e < mol.natm ):
                    f.write( "obj.mass[%d] = %.4lf\n"%( _e, _m ) )
            except:
                pass
        f.write( """
obj.get_hess()
frq, mds = qm3.utils.hessian_frequencies( obj.mass, obj.coor, obj.hess, True )
print( "Frequencies (cm^-1):" )
print( frq )
zpe, gib = qm3.utils.gibbs_rrho( obj.mass, obj.coor, frq, temp = %.2lf, press = %.2lf )
print( "ZPE:", zpe )
print( "Gibbs (298K / 1atm):", gib )
print( "Total:", zpe + gib, "_kJ/mol" )

i = 0
while( i < obj.size and math.fabs( frq[i] ) < 10. ):
    i += 1
if( obj.sele == [] ):
    s = obj.mole.guess_symbols()
else:
    s = obj.mole.guess_symbols( obj.sele )
if( frq[i] < .0 ):
    qm3.utils.normal_mode_view( obj.coor, frq, mds, s, i, afac = 4 )
else:
    qm3.utils.normal_mode_view( obj.coor, frq, mds, s, i )
"""%( _tmp, _pre ) )
    # -- stop engines (if needed...)
    if( emm != "--" and MM_engines[emm] in [ "qm3.engines.charmm.run_shmem", "qm3.engines.namd.run_shmem" ] ):
        f.write( "\nobj.emm.stop()\n" )
    f.close()
    return( "<h1>Check the script generated:</h1><h1>%s</h1>"%( job ) )


bottle.TEMPLATE_PATH.append( os.path.join( cwd, "views" ) )
bottle.run( host = "localhost", port = 8080 )

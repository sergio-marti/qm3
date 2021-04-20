import  qm3.mol
import  qm3.elements
import  qm3.problem
import  qm3.engines.mol_mech
import  qm3.engines.xtb
import  qm3.actions.minimize
import  qm3.actions.neb
import  qm3.utils.msi
import  os
import  sys
import  pickle
import  math
import  qm3.utils



class parmsi_neb( qm3.problem.template ):
    def __init__( self, problem, sele, kumb, conn, ncpu, guess ):
        qm3.problem.template.__init__( self )

        self.node = len( guess ) - 2
        self.conn = conn
        self.wami = conn.node
        self.ncpu = ncpu
        self.chnk = self.node // ncpu

        self.kumb = kumb
        self.sele = sele[:]
        self.crd0 = guess[ 0][:]
        self.crdf = guess[-1][:]
        self.prob = problem

        self.dime = len( self.crd0 )
        self.size = self.dime * self.node
        self.coor = []
        for k in range( 1, self.node + 1 ):
            self.coor += guess[k][:]
        self.func = 0
        self.grad = []

        # end-points potential
        self.prob.coor = self.crd0[:]
        self.prob.get_func()
        self.pot0 = self.prob.func
        self.prob.coor = self.crdf[:]
        self.prob.get_func()
        self.potf = self.prob.func
        if( self.wami == 0 ):
            print( self.pot0, self.potf )


    def get_grad( self ):
        # ----------------------------------------------------------------------
        def __calc_tau( potm, poti, potp, crdm, crdi, crdp ):
            dcM = [ ii-jj for ii,jj in zip( crdp, crdi ) ]
            dcm = [ ii-jj for ii,jj in zip( crdi, crdm ) ]
            dpM = max( math.fabs( potp - poti ), math.fabs( potm - poti ) )
            dpm = min( math.fabs( potp - poti ), math.fabs( potm - poti ) )
            if( potp > poti and poti > potm ):
                tau = dcM[:]
            elif( potp < poti and poti < potm ):
                tau = dcm[:]
            else:
                if( potp > potm ):
                    tau = [ dpM * ii + dpm * jj for ii,jj in zip( dcM, dcm ) ]
                else:
                    tau = [ dpm * ii + dpM * jj for ii,jj in zip( dcM, dcm ) ]
            tmp = math.sqrt( sum( [ ii * ii for ii in tau ] ) )
            tau = [ ii / tmp for ii in tau ]
            gum = [ self.kumb * ( ii - jj ) for ii,jj in zip( dcm, dcM ) ]
            tmp = sum( [ ii * jj for ii,jj in zip( gum, tau ) ] )
            gum = [ tmp * ii for ii in tau ]
            return( tau, gum )
        # ----------------------------------------------------------------------
        vpot = []
        gpot = []
        # get potential energy and gradients from the chain elements
        for who in range( self.wami * self.chnk, ( self.wami + 1 ) * self.chnk ):
            self.prob.coor = self.coor[self.dime*who:self.dime*who+self.dime][:]
            self.prob.get_grad()
            vpot.append( self.prob.func )
            gpot += self.prob.grad[:]
            try:
                self.prob.neb_data( who )
            except:
                pass
        # sync vpot from and to nodes
        self.conn.barrier()
        if( self.wami == 0 ):
            for i in range( 1, self.ncpu ):
                vpot += self.conn.recv( i, self.chnk )
            for i in range( 1, self.ncpu ):
                self.conn.send( i, vpot )
        else:
            self.conn.send( 0, vpot )
            vpot = self.conn.recv( 0, self.node )
        # sync gpot from and to nodes
        self.conn.barrier()
        if( self.wami == 0 ):
            for i in range( 1, self.ncpu ):
                gpot += self.conn.recv( i, self.dime * self.chnk )
            for i in range( 1, self.ncpu ):
                self.conn.send( i, gpot )
        else:
            self.conn.send( 0, gpot )
            gpot = self.conn.recv( 0, self.size )
        self.func = sum( vpot )
        self.grad = []
        # calculate neb components    
        for who in range( self.node ):
            if( who == 0 ):
                # first node
                tau, gum = __calc_tau( self.pot0, vpot[who], vpot[who+1],
                            self.crd0,
                            self.coor[self.dime*who:self.dime*who+self.dime],
                            self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime] )
            elif( who == self.node - 1 ):
                # last node
                tau, gum = __calc_tau( vpot[who-1], vpot[who], self.potf,
                            self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime],
                            self.coor[self.dime*who:self.dime*who+self.dime],
                            self.crdf )
            else:
                # intermediates
                tau, gum = __calc_tau( vpot[who-1], vpot[who], vpot[who+1],
                            self.coor[self.dime*(who-1):self.dime*(who-1)+self.dime],
                            self.coor[self.dime*who:self.dime*who+self.dime],
                            self.coor[self.dime*(who+1):self.dime*(who+1)+self.dime] )
            # common to all nodes
            gsp = sum( [ ii * jj for ii,jj in zip( tau, gpot[self.dime*who:self.dime*who+self.dime] ) ] )
            grd = [ ii - gsp * jj + kk for ii,jj,kk in zip( gpot[self.dime*who:self.dime*who+self.dime], tau, gum ) ]
            qm3.utils.project_RT_modes( self.prob.mass, self.coor[self.dime*who:self.dime*who+self.dime], grd, [] )
            self.grad += grd[:]



class envi_problem( qm3.problem.template ):
    def __init__( self, molec ):
        qm3.problem.template.__init__( self )

        self.mol = molec

        self.emm = qm3.engines.mol_mech.simple_force_field( self.mol )
        self.emm.cut_on   = 12.0
        self.emm.cut_off  = 14.0
        self.emm.cut_list = 18.0
        self.emm.system_read( "../02.mm_data.pk", self.mol )

        f = open( "../02.sele_MM.pk", "rb" )
        self.sele = pickle.load( f )
        f.close()

        f = open( "../01.sele_QM.pk", "rb" )
        for i in pickle.load( f ):
            self.emm.free[i] = False
        f.close()
        f = open( "../02.fixed.pk", "rb" )
        for i in pickle.load( f ):
            self.emm.free[i] = False
        f.close()

        self.size = len( self.sele ) * 3
        self.coor = []
        for i in self.sele:
            i3 = i * 3
            self.coor += self.mol.coor[i3:i3+3]

        self.log = open( "log", "wt" )


    def feed_log( self, txt ):
        self.log.write( txt + "\n" )


    def update_coor( self ):
        for i in range( len( self.sele ) ):
            i3 = i * 3
            I3 = self.sele[i] * 3
            for j in [0, 1, 2]:
                self.mol.coor[I3+j] = self.coor[i3+j]


    def get_func( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.emm.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.mol.natm * 3 ) ]
        self.emm.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = []
        for i in self.sele:
            i3 = i * 3
            self.grad += self.mol.grad[i3:i3+3]



class core_problem( qm3.problem.template ):
    def __init__( self ):
        qm3.problem.template.__init__( self )

        self.mol = qm3.mol.molecule( "../03.pes.02.15.pdb" )
        self.mol.boxl = [ 40., 40., 40. ]

        self.envi = envi_problem( self.mol )

        f = open( "../01.sele_QM.pk", "rb" )
        self.sele = pickle.load( f )
        f.close()
        f = open( "../02.sele_MM.pk", "rb" )
        smm = pickle.load( f )
        f.close()

        self.eqm = qm3.engines.xtb.dl_xtb( self.mol, 0, 0, self.sele, smm )

        self.size = len( self.sele ) * 3
        self.coor = []
        self.mass = []
        for i in self.sele:
            i3 = i * 3
            self.coor += self.mol.coor[i3:i3+3]
            self.mass.append( qm3.elements.mass[self.mol.anum[i]] )

        self.mm_flg = False
        self.get_func()
        self.mm_flg = True


    def neb_data( self, who ):
        f = open( "../03.neb.%02d.pdb"%( who ), "wt" )
        f.write( "REMARK: %.6lf\n"%( self.func ) )
        self.mol.pdb_write( f )
        f.close()


    def update_coor( self ):
        for i in range( len( self.sele ) ):
            i3 = i * 3
            I3 = self.sele[i] * 3
            for j in [0, 1, 2]:
                self.mol.coor[I3+j] = self.coor[i3+j]
        if( self.mm_flg ):
            self.envi.feed_log( str( self.mol.chrg[0:17] ) )
            qm3.actions.minimize.fire( self.envi, print_frequency = 10, gradient_tolerance = 0.2,
                step_size = 0.1, step_number = 200, log_function = self.envi.feed_log )


    def get_func( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.eqm.get_func( self.mol )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.mol.natm * 3 ) ]
        self.eqm.get_grad( self.mol )
        self.func = self.mol.func
        self.grad = []
        for i in self.sele:
            i3 = i * 3
            self.grad += self.mol.grad[i3:i3+3]





ncpu = 40
node = int( sys.argv[1] )

os.mkdir( "scratch_%02d"%( node ) )
os.chdir( "scratch_%02d"%( node ) )

conn = qm3.utils.msi.client( node, "../qm3_msi" )

if( node == 0 ):
    def gen_log( txt ):
        sys.stdout.write( txt + "\n" )
        sys.stdout.flush()
else:
    def gen_log( txt ):
        pass

f = open( "../01.sele_QM.pk", "rb" )
sel = pickle.load( f )
f.close()

crd = []
mol = qm3.mol.molecule( "../03.pes.02.15.pdb" )
tmp = []
for i in sel:
    tmp += mol.coor[3*i:3*i+3][:]
crd.append( tmp[:] )

mol.pdb_read( "../03.pes.16.01.pdb" )
tmp = []
for i in sel:
    tmp += mol.coor[3*i:3*i+3][:]
crd.append( tmp[:] )

crd = qm3.actions.neb.distribute( 24, crd )

obj = parmsi_neb( core_problem(), sel, 200., conn, ncpu, crd )
qm3.actions.minimize.fire( obj, step_number = 1000, print_frequency = 1, step_size = 0.1,
        gradient_tolerance = 0.1, log_function = gen_log, fire2 = True, exit_uphill = True )

conn.barrier()

#
# mpirun -n 8 python 03-scan.py
#
import  qm3.mol
import  qm3.problem
import  qm3.engines.mol_mech
import  qm3.engines.xtb
import  qm3.engines.mmres
import  qm3.actions.minimize
import  qm3.utils._mpi
import  os
import  time
import  pickle
import  random


# ---------------------------------------------------
# >> problem config

class my_problem( qm3.problem.template ):
    def __init__( self ):
        qm3.problem.template.__init__( self )

        self.mol = qm3.mol.molecule( "../02.pdb" )
        self.mol.boxl = [ 40., 40., 40. ]

        self.emm = qm3.engines.mol_mech.simple_force_field( self.mol, 8 )
        self.emm.cut_on   = 12.0
        self.emm.cut_off  = 14.0
        self.emm.cut_list = 18.0
        self.emm.system_read( "../02.mm_data.pk", self.mol )

        f = open( "../01.sele_QM.pk", "rb" )
        sqm = pickle.load( f )
        f.close()
        f = open( "../02.sele_MM.pk", "rb" )
        smm = pickle.load( f )
        f.close()
        f = open( "../02.sele.pk", "rb" )
        self.sele = pickle.load( f )
        f.close()
        f = open( "../02.fixed.pk", "rb" )
        for i in pickle.load( f ):
            self.emm.free[i] = False
        f.close()

        self.emm.qm_atoms( sqm )

        self.eqm = qm3.engines.xtb.dl_xtb( self.mol, 0, 0, sqm, smm )

        self.umb = []
        self.dat = []

        self.size = len( self.sele ) * 3
        self.coor = []
        for i in self.sele:
            i3 = i * 3
            self.coor += self.mol.coor[i3:i3+3]


    def grid_point( self, ci, cj, fname = None ):
        self.umb = []
        self.umb.append( qm3.engines.mmres.distance( 5000.0, rx1[0] + ci * dx1,
            [ self.mol.indx["A"][1]["C2"], self.mol.indx["A"][1]["C3"] ] ) )
        self.umb.append( qm3.engines.mmres.distance( 5000.0, rx2[0] + cj * dx2,
            [ self.mol.indx["A"][1]["C5"], self.mol.indx["A"][1]["C6"] ] ) )
        if( fname != None ):
            tmp = qm3.mol.molecule( fname )
            self.mol.coor = tmp.coor[:]
            self.coor = []
            for i in self.sele:
                i3 = i * 3
                self.coor += self.mol.coor[i3:i3+3]


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
        self.eqm.get_func( self.mol )
        self.dat = []
        for itm in self.umb:
            self.dat.append( itm.get_func( self.mol ) )
        self.func = self.mol.func


    def get_grad( self ):
        self.update_coor()
        self.mol.func = 0.0
        self.mol.grad = [ 0.0 for i in range( self.mol.natm * 3 ) ]
        self.emm.get_grad( self.mol )
        self.eqm.get_grad( self.mol )
        self.dat = []
        for itm in self.umb:
            self.dat.append( itm.get_grad( self.mol ) )
        self.func = self.mol.func
        self.grad = []
        for i in self.sele:
            i3 = i * 3
            self.grad += self.mol.grad[i3:i3+3]



# ---------------------------------------------------
# >> pes config

rx1 = [ 1.4, 3.5 ] 
dx1 = 0.1

rx2 = [ 1.4, 3.5 ]
dx2 = 0.1




# ---------------------------------------------------
# scan!

random.seed()

nx1 = int( ( rx1[1] - rx1[0] ) / dx1 ) + 1
nx2 = int( ( rx2[1] - rx2[0] ) / dx2 ) + 1

pid, cpu = qm3.utils._mpi.init()

os.mkdir( "scratch_%02d"%( pid ) )
os.chdir( "scratch_%02d"%( pid ) )

efd = open( "../03.ene.%02d"%( pid ), "wt" )
lfd = open( "log", "wt" )

def feed_log( txt ):
    lfd.write( txt + "\n" )

obj = my_problem()

if( pid == 0 ):
    print( "rnge:", nx1, nx2 )
    print( "ncpu:", cpu )

    a1 = obj.mol.indx["A"][1]["C2"] * 3
    a2 = obj.mol.indx["A"][1]["C3"] * 3
    ci = int( round( ( qm3.utils.distance( obj.mol.coor[a1:a1+3], obj.mol.coor[a2:a2+3] ) - rx1[0] ) / dx1, 0 ) )
    ci = min( max( 0, ci ), nx1 - 1 )
    a1 = obj.mol.indx["A"][1]["C5"] * 3
    a2 = obj.mol.indx["A"][1]["C6"] * 3
    cj = int( round( ( qm3.utils.distance( obj.mol.coor[a1:a1+3], obj.mol.coor[a2:a2+3] ) - rx2[0] ) / dx2, 0 ) )
    cj = min( max( 0, cj ), nx2 - 1 )
    
    feed_log( "%d %d %d"%( pid, ci, cj ) )
    obj.grid_point( ci, cj )
    qm3.actions.minimize.fire( obj, print_frequency = 100,
        gradient_tolerance = 0.1, step_number = 5000, fire2 = True, log_function = feed_log )
    obj.mol.pdb_write( "../03.pes.%02d.%02d.pdb"%( ci, cj ) )
    efd.write( "%20.10lf%20.10lf%20.10lf\n"%( obj.dat[0], obj.dat[1], obj.func ) )
    efd.flush()
    
    flg = []
    par = []
    idx = []
    n   = nx1 * nx2
    for i in range( nx1 ):
        for j in range( nx2 ):
            flg.append( False )
            par.append( (i,j) )
            idx.append( i*nx2+j )
    flg[ci*nx2+cj] = True
    mov = [ (-1,0), (+1,0), (0,-1), (0,+1) ]
    for i in range( max( n // 100, 10 ) ):
        random.shuffle( mov )
        random.shuffle( idx )
    lst = []
    for j in range( 1, n ):
        random.shuffle( idx )
        q = True
        i = 0
        while( q and i < n ):
            if( not flg[idx[i]] ):
                random.shuffle( mov )
                k = 0
                while( q and k < 4 ):
                    I = par[idx[i]][0] + mov[k][0]
                    J = par[idx[i]][1] + mov[k][1]
                    w = I * nx2 + J
                    if( I >= 0 and I < nx1 and J >= 0 and J < nx2 and flg[w] ):
                        q = False
                        lst.append( ( par[idx[i]][0], par[idx[i]][1], "../03.pes.%02d.%02d.pdb"%( I, J ) ) )
                        flg[idx[i]] = True
                    k += 1
            i += 1
    
    tsk = [ [] for i in range( cpu ) ]
    for i in range( n - 1 ):
        tsk[i%cpu].append( lst[i] )

    f = open( "../03.tasks.pk", "wb" )
    pickle.dump( tsk, f )
    f.close()

    qm3.utils._mpi.barrier()

else:
    qm3.utils._mpi.barrier()

    f = open( "../03.tasks.pk", "rb" )
    tsk = pickle.load( f )
    f.close()


for i,j,s in tsk[pid]:
    feed_log( "%d %d %d %s"%( pid, i, j, s ) )
    while( not os.path.isfile( s ) ):
        if( pid == 0 ):
            feed_log( ">> %d %d waiting for: %s"%( i, j, s ) )
        time.sleep( 1 )
    time.sleep( 2 )
    obj.grid_point( i, j, s )
    qm3.actions.minimize.fire( obj, print_frequency = 100, gradient_tolerance = 0.1,
        step_number = 1000, fire2 = True, log_function = feed_log )
    obj.mol.pdb_write( "../03.pes.%02d.%02d.pdb"%( i, j ) )
    efd.write( "%20.10lf%20.10lf%20.10lf\n"%( obj.dat[0], obj.dat[1], obj.func ) )
    efd.flush()

efd.close()
lfd.close()

qm3.utils._mpi.barrier()
qm3.utils._mpi.stop()

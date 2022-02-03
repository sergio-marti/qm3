# -*- coding: iso-8859-1 -*-
import qm3.constants
import qm3.fio
import re
import os
import time



#    def xyz_read( self, fname = None ):
#        self.conn = []
#        self.atyp = []
#        f = qm3.fio.open_r( fname )
#        self.natm = int( f.readline().split()[0] )
#        for i in range( self.natm ):
#            t = f.readline().split()
#            self.labl.append( t[1] )
#            self.coor += [ float( t[2] ), float( t[3] ), float( t[4] ) ]
#            self.atyp.append( int( t[5] ) )
#            self.conn.append( t[6:][:] )
#        qm3.fio.close( f, fname )
#    
#
#    def xyz_write( self, fname = None ):
#        f = qm3.fio.open_w( fname )
#        f.write( "%8d\n"%( self.natm ) )
#        for i in range( self.natm ):
#            i3 = i * 3
#            f.write( "%8d  %3s%14.6lf%14.6lf%14.6lf%6d%s\n"%( i+1, self.labl[i],
#                self.coor[i3], self.coor[i3+1], self.coor[i3+2], self.atyp[i], 
#                "".join( [ "%8s"%( j ) for j in self.conn[i] ] ) ) )
#        qm3.fio.close( f, fname )



class tinker( object ):

    def __init__( self, ini, cpu = os.sysconf( 'SC_NPROCESSORS_ONLN' ) ):
        f = open( "tinker.key", "wt" )
        f.write( ini )
        f.close()
        self.exe = "mpirun -n %d /Users/smarti/Devel/tinker/hp-1.0/bin/testgrad tinker Y N N > tinker.log"%( cpu )


    def update_coor( self, mol ):
        mol.xyz_write( "tinker.xyz" )


    def get_func( self, mol ):
        self.update_coor( mol )
        os.system( self.exe )
        f = open( "tinker.log", "rt" )
        l = f.readline()
        while( l != "" ):
            L = l.strip()
            if( L[0:22] == "Total Potential Energy" ):
                mol.func += float( l.split()[4] ) * qm3.constants.K2J
                l = ""
            else:
                l = f.readline()
        f.close()


    def get_grad( self, mol ):
        self.update_coor( mol )
        self.get_func( mol )
        f = open( "tinker.log", "rt" )
        p = re.compile( "Anlyt[\ ]*([0-9]+)[\ ]*([0-9\.\-]+)[\ ]*([0-9\.\-]+)[\ ]*([0-9\.\-]+)[\ ]*[0-9\.]+" )
        l = f.readline()
        while( l != "" ):
            L = l.strip()
            if( p.match( L ) ):
                t = L.split()
                i = 3 * ( int( t[1] ) - 1 )
                for j in [0, 1, 2]:
                    mol.grad[i+j] += float( t[2+j] ) * qm3.constants.K2J
            l = f.readline()
        f.close()








TINKER_INP = """parameters    params/amoebabio09
verbose
a-axis        62.23
vdw-cutoff    9.0 
ewald
ewald-cutoff  7.0
pme-grid      64 64 64
"""

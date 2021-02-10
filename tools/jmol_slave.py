#!/usr/bin/env python3
import  socket
import  json
import  time
import  subprocess


class jmol:
    def __init__( self ):
        f = open( "jmol", "wt" )
        f.write( "sync -6969\nset platformSpeed 2\n" )
        f.close()
        self.proc = subprocess.Popen( ["/usr/bin/java", "-Xmx1024m",
            "-jar", "/Applications/JMol.app/Contents/MacOS/Jmol.jar", "-s", "jmol" ],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL )
        self.sckt = socket.socket( socket.AF_INET, socket.SOCK_STREAM )
        time.sleep( 5 )
        self.sckt.connect( ( "localhost", 6969 ) )

    def mol_view( self, mol, sel ):
        mol.pdb_write( "jmol", sele = sel )
        self.sckt.send( bytes( json.dumps( {"type" : "command", "command" : "load jmol\nlabel %a" }  ) + "\r\n", "ascii" ) )
        input( "-- hit enter --" )

    def stop( self ):
        self.sckt.close()
        self.proc.terminate()

import qm3.utils.pes_samples
import qm3.actions.paths


class my_obj( qm3.utils.pes_samples.muller_brown ):
    def __init__( self ):
        qm3.utils.pes_samples.muller_brown.__init__( self )
        self.fd = open( "path.log", "wt" )

    def current_step( self, istep ):
        self.fd.write( "%20.10lf%20.10lf%20.10lf\n"%( self.coor[0], self.coor[1], self.func ) )
        

obj = my_obj()

obj.coor = [-0.822001668377818, 0.6243127881478763]
qm3.actions.paths.page_mciver( obj, step_number = 1000, step_size =  0.0028, project_RT = False )
print( obj.coor )
obj.fd.write( "\n" )

obj.coor = [-0.822001668377818, 0.6243127881478763]
qm3.actions.paths.page_mciver( obj, step_number = 1000, step_size = -0.0028, project_RT = False )
print( obj.coor )
obj.fd.close()

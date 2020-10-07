import qm3.actions.dynamics
import qm3.problem
import math


class sample( qm3.problem.template ):
    def __init__( self ):
        self.size = 6
        self.func = 0.0
        self.mass = [ 1.0, 1.0 ]
        self.coor = [ 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 ]
        self.grad = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
        self.velo = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
        self.hess = None
        self.ref = 1.0
        self.kmb = 10.0

    def get_func( self ):
        dst = math.sqrt( sum( [ math.pow( self.coor[i] - self.coor[3+i], 2.0 ) for i in [ 0, 1, 2 ] ] ) )
        self.func = 0.5 * self.kmb * ( dst - self.ref ) * ( dst - self.ref ) 
    
    def get_grad( self ):
        dr  = [ self.coor[i] - self.coor[3+i] for i in [ 0, 1, 2 ] ]
        dst = math.sqrt( sum( [ dr[i] * dr[i] for i in [ 0, 1, 2 ] ] ) )
        self.func = 0.5 * self.kmb * ( dst - self.ref ) * ( dst - self.ref ) 
        tmp = self.kmb * ( dst - self.ref )
        self.grad = [ tmp * dr[0], tmp * dr[1], tmp * dr[2], - tmp * dr[0], - tmp * dr[1], - tmp * dr[2] ]


obj = sample()

obj.coor = [ 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 ]
qm3.actions.dynamics.assign_velocities( obj, 10.0, project = False )
job = qm3.actions.dynamics.langevin_verlet( obj, step_size = 0.001, temperature = 10.0, gamma_factor = 100.0,
    print_frequency = 100, project = False, step_number = 1000 )

obj.coor = [ 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 ]
qm3.actions.dynamics.assign_velocities( obj, 10.0, project = False )
job = qm3.actions.dynamics.velocity_verlet( obj, step_size = 0.001, temperature = 10.0, scale_frequency = 100,
    print_frequency = 50, project = False )
for i in range( 1000 ):
    job.integrate()
job.stats()

obj.coor = [ 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 ]
qm3.actions.dynamics.assign_velocities( obj, 10.0, project = False )
job = qm3.actions.dynamics.velocity_verlet( obj, step_size = 0.001, temperature = 10.0, scale_frequency = 0,
    temperature_coupling = 0.1, print_frequency = 50, project = False, step_number = 1000 )

obj.coor = [ 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 ]
qm3.actions.dynamics.assign_velocities( obj, 10.0, project = False )
qm3.actions.dynamics.velocity_verlet( obj, step_size = 0.001, temperature = 10.0, scale_frequency = 0,
    temperature_coupling = 0.0, print_frequency = 100, project = False, step_number = 1000 )

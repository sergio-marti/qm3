import qm3.actions.minimize
import qm3.utils.pes_samples

obj = qm3.utils.pes_samples.muller_brown()

obj.coor = [ -0.5440, 0.5094 ]
qm3.actions.minimize.baker( obj, follow_mode = 0 )
print( obj.coor )

obj.coor = [ -0.5440, 0.5094 ]
qm3.actions.minimize.rfo( obj, follow_mode = 0 )
print( obj.coor )

obj.coor = [ -1.1067, 0.8016 ]
qm3.actions.minimize.baker( obj )
print( obj.coor )

obj.coor = [ -1.1067, 0.8016 ]
qm3.actions.minimize.fire( obj )
print( obj.coor )

obj.coor = [ -1.1067, 0.8016 ]
qm3.actions.minimize.l_bfgs( obj )
print( obj.coor )

obj.coor = [ -1.1067, 0.8016 ]
qm3.actions.minimize.conjugate_gradient_plus( obj )
print( obj.coor )

obj.coor = [ -1.1067, 0.8016 ]
qm3.actions.minimize.steepest_descent( obj, step_number = 1000, print_frequency = 100,
    step_size = 0.001, gradient_tolerance = 0.1 )
print( obj.coor )

obj.coor = [ -1.1067, 0.8016 ]
qm3.actions.minimize.adam( obj, step_number = 1000, print_frequency = 10,
    step_size = 0.1, gradient_tolerance = 0.1 )
print( obj.coor )

obj.coor = [ -1.1067, 0.8016 ]
qm3.actions.minimize.downhill( obj )
print( obj.coor )

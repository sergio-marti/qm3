import qm3.problem

class my_problem( qm3.problem.template ):
    def __init__( self ):
        self.size = 2
        self.coor = [ .0, .0 ]
        self.func = .0
        self.grad = []
        self.hess = []

    def get_func( self ):
        self.func = self.coor[0] * self.coor[0] + self.coor[1] * self.coor[1]

    def get_grad( self ):
        self.get_func()
        self.grad = [ 2.0 * self.coor[0], 2.0 * self.coor[1] ]

    def get_hess( self ):
        self.get_grad()
        self.hess = [ 2.0, 0.0, 0.0, 2.0 ]


obj = my_problem()
obj.coor = [ 1.0, 1.0 ]
obj.get_hess()
print( obj.func, obj.grad, obj.hess )
obj.num_hess()
print( obj.hess )

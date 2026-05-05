from resolution_cube import ResolutionCube
from khov_algebra import multiply, comultiply

class KhovanovComplex:
    
    def __init__(self, dtcode):
        self.cube = ResolutionCube(dtcode)
        self.n = self.cube.n
        self.state_circles = {}
        self.state_basis = {}

    def cache_state_circles(self, state):
        if state not in self.state_circles:
            self.state_circles[state] = self.cube.get_circles(state)
        return self.state_circles[state]
    
    def basis_for_state(self, state):
        if state in self.state_basis:
            return self.state_basis[state]
        
        # we have k circles in this state
        k = len(self.cache_state_circles(state)) 
        basis = []

        for labels in range(pow(2, k)):
            labels_tuple = tuple((labels>>i) & 1 for i in reversed(range(k)))
            basis.append((state, labels_tuple))
        
        self.state_basis[state] = basis
        return basis
    
    def chain_basis(self, degree):
        basis = []
        for state in range(1<<self.n):
            if state.bit_count() == degree:
                basis.extend(self.basis_for_state(state))
        return basis
    
    


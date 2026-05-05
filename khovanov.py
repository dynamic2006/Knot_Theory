from resolution_cube import ResolutionCube
from khov_algebra import multiply, comultiply

class KhovanovComplex:
    
    def __init__(self, dtcode):
        self.cube = ResolutionCube(dtcode)
        self.n = self.cube.n
        self.state_circles = {}
        self.state_basis = {}

    def cache_state_circles(self, state):
        """
        Returns list of circles for given state
        """
        if state not in self.state_circles:
            self.state_circles[state] = self.cube.get_circles(state)
        return self.state_circles[state]
    
    def get_basis_for_state(self, state):
        """
        Returns basis: 
        list of basis vectors (each circle can be 1 or X assigned)
        """
        if state in self.state_basis:
            return self.state_basis[state]
        
        # we have k circles in this state
        k = len(self.cache_state_circles(state)) 
        basis = []

        for labels in range(pow(2, k)):
            # recall convention that 0=>1 and 1=>X
            # we generate all possible labeling assignments to circles
            # basis is the list of all of these for a given state
            labels_tuple = tuple((labels>>i) & 1 for i in reversed(range(k)))
            basis.append((state, labels_tuple))
        
        self.state_basis[state] = basis
        return basis
    
    def get_chain_basis_for_degree(self, degree):
        """
        Combines all parts of the cube that are together
        ie, have same degree/height, giving full chain basis
        """
        basis = []
        for state in range(1<<self.n):
            if state.bit_count() == degree:
                basis.extend(self.get_basis_for_state(state))
        return basis
    
    def get_state_from_resolving_crossing(self, state, crossing):
        return state | (1<<crossing)
    
    def get_edge_kind(self, state, crossing):
        
        src_state = state
        dest_state = self.get_state_from_resolving_crossing(state, crossing)
        src_circles = self.cache_state_circles(src_state)
        dest_circles = self.cache_state_circles(dest_state)

        if len(dest_circles) == len(src_circles)-1:
            return "merge"
        elif len(dest_circles) == len(src_circles)+1:
            return "split"

    def get_vertex_to_circle_map(self, state):
        circles = self.cache_state_circles(state)
        vertex_to_circle = {}
        for idx, circle in enumerate(circles):
            for vertex in circle:
                vertex_to_circle[vertex] = idx
        return vertex_to_circle

    def get_edge_map_info(self, state, crossing):
        """
        Returns details for what circles are involved in crossing resolution.
        merge edge: ("merge", src_circ_a, src_circ_b, dest_circ)
        split edge: ("merge", src_circ, dest_circ_a, dest_circ_b)
        """
        dest_state = self.get_state_from_resolving_crossing(state, crossing)
        kind = self.get_edge_kind(state, crossing)
        src_vertex_to_circ_map = self.get_vertex_to_circle_map(state)
        dest_vertex_to_circ_map = self.get_vertex_to_circle_map(dest_state)

        # the vertices of our concern -- the ones at this crossing!
        vertices = [
            self.cube.get_vertex_id(crossing, 0),
            self.cube.get_vertex_id(crossing, 1),
            self.cube.get_vertex_id(crossing, 2),
            self.cube.get_vertex_id(crossing, 3)
        ]

        src_circs = sorted({src_vertex_to_circ_map[v] for v in vertices})
        dest_circs = sorted({dest_vertex_to_circ_map[v] for v in vertices})

        if kind == "merge":
            return ("merge", src_circs[0], src_circs[1], dest_circs[0])
        else:
            return ("split", src_circs[0], dest_circs[0], dest_circs[1])
    
    def get_src_to_dest_circs_map(self, state, crossing):
        dest_state = self.get_state_from_resolving_crossing(state, crossing)
        src_circs = self.cache_state_circles(state)
        vertex_to_circ_map = self.get_vertex_to_circle_map(dest_state)

        src_to_dest_circs_map = {}

        for i, circ in enumerate(src_circs):
            dest_circs = {vertex_to_circ_map[v] for v in circ}
            src_to_dest_circs_map[i] = dest_circs

        return src_to_dest_circs_map


    def apply_edge_map(self, state, crossing, labels):
        """
        Apply edge map to one basis vector:
        - at starting state
        - resolve crossing
        Returns list of (coefficient, dest_labels)
        """
        dest_state = self.get_state_from_resolving_crossing(state, crossing)
        info = self.get_edge_map_info(state, crossing)

        kind = info[0]

        if kind == "merge":
            _, a, b, c = info
            outputs = multiply(labels[a], )
    

K = KhovanovComplex("BCA")

for state in range(1 << K.n):
    for crossing in range(K.n):
        if ((state >> crossing) & 1) == 0:
            bits = format(state, f"0{K.n}b")
            print(bits, crossing, K.get_edge_map_info(state, crossing))
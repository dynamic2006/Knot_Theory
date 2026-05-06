from resolution_cube import ResolutionCube
from khov_algebra import multiply, comultiply

# being able to index on row and col makes life easier
# so that we can just do cols[col][row] = coeff
def sparse_entries_to_cols(M, num_cols):
    """
    Convert: M[col] = [(row, coeff), ...]
    into: cols[col] = {row: coeff}
    """
    cols = []
    for col in range(num_cols):
        col_dict = {}
        for row, coeff in M.get(col, []):
            col_dict[row] = col_dict.get(row, 0) + coeff
            if col_dict[row] == 0:
                del col_dict[row] # need to delete so pivoting works
        cols.append(col_dict)
    return cols

# compute how many columns are linearly indep
def sparse_rank_q(cols):
    """
    Sparse Column Reduction
    Computes rank over Q
    cols: list of dicts {row : coeff}
    """
    pivots = {} # carrying out gaussian elimination
    rank = 0

    for col in cols:
        col = dict(col) # copy so we don't modify og
        while col:
            pivot_row = max(col.keys())
            # if we find a new indep col
            if pivot_row not in pivots:
                pivots[pivot_row] = col
                rank+=1
                break
            
            # otherwise reduce
            pivot_col = pivots[pivot_row]
            a = col[pivot_row]
            b = pivot_col[pivot_row]

            # b*col - a*pivot_col ==> col
            new_col = {}
            for r, v in col.items():
                new_col[r] = new_col.get(r, 0) + b*v
            for r, v in pivot_col.items():
                new_col[r] = new_col.get(r,0) - a*v
            col = {r : v for r,v in new_col.items() if v != 0}
    
    return rank

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
        split edge: ("split", src_circ, dest_circ_a, dest_circ_b)
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

    # whole point of this is to make cube circle labeling cohesive globally
    # so later when we build the sparse matrix
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
        src_to_dest_circs_map = self.get_src_to_dest_circs_map(state, crossing)
        dest_circ_count = len(self.cache_state_circles(dest_state))

        dest_results = []

        if kind == "merge":
            _, a, b, c = info
            results = multiply(labels[a], labels[b])
            for coeff, merged_label in results:
                dest_labels = [None]*dest_circ_count
                # Handle unchanged circles
                for src_circ, dest_circs in src_to_dest_circs_map.items():
                    if src_circ == a or src_circ == b:
                        continue # this was a CHANGED circle
                    dest_circ = next(iter(dest_circs))
                    dest_labels[dest_circ] = labels[src_circ] # inherits label from src_circ
                dest_labels[c] = merged_label # Handle changed circles
                dest_results.append((coeff, tuple(dest_labels)))
            
        if kind == "split":
            _, a, b, c = info
            results = comultiply(labels[a])
            for coeff, split_label in results:
                dest_labels = [None]*dest_circ_count
                # Handle unchanged circles
                for src_circ, dest_circs in src_to_dest_circs_map.items():
                    if src_circ == a:
                        continue # this was a CHANGED circle
                    dest_circ = next(iter(dest_circs))
                    dest_labels[dest_circ] = labels[src_circ]
                dest_labels[b] = split_label[0]
                dest_labels[c] = split_label[1]
                dest_results.append((coeff, tuple(dest_labels)))

        return dest_results

    # see Bar-Natan p.343  ==> we need this so every square face anticommutes
    # d_j d_i = - d_i d_j
    def cube_sign(self, state, crossing):
        ones_before = (state & ((1<<crossing)-1)).bit_count()
        return -1 if ones_before%2 else 1

    def differential_sparse(self, degree):
        """
        Builds sparse matrix for degree: C_degree -> C_{degree+1}
        Returns : M, num_rows, num_cols
        M[col] = list of (row, coeff)
        """
        src_basis = self.get_chain_basis_for_degree(degree)
        dest_basis = self.get_chain_basis_for_degree(degree+1)

        dest_basis_to_rownum_map = {}
        for row, basis_vector in enumerate(dest_basis):
            dest_basis_to_rownum_map[basis_vector] = row

        M = {}

        for col, (state, labels) in enumerate(src_basis):
            entries = []
            for crossing in range(self.n):
                # if valid edge in cube
                if((state>>crossing) & 1) == 0:
                    dest_state = self.get_state_from_resolving_crossing(state, crossing)
                    dest_results = self.apply_edge_map(state, crossing, labels)

                    sign = self.cube_sign(state, crossing)
                    for coeff, dest_labels in dest_results:
                        dest_basis_vector = (dest_state, dest_labels)
                        row = dest_basis_to_rownum_map[dest_basis_vector]
                        entries.append((row, sign*coeff))
            if entries:
                M[col] = entries
        
        return M, len(dest_basis), len(src_basis)

    def rank_differential_q(self, degree):
        M, num_rows, num_cols = self.differential_sparse(degree)
        cols = sparse_entries_to_cols(M, num_cols)
        return sparse_rank_q(cols)

    def get_free_rank(self, degree):
        """
        Free rank of H_degree
        [IGNORING TORSION]
        """
        C_dim = len(self.get_chain_basis_for_degree(degree))
        
        # check if we've gone off the right edge
        if degree < self.n:
            rank_d_i = self.rank_differential_q(degree)
        else:
            rank_d_i = 0
        
        # check if going left goes off the left edge
        if degree > 0:
            rank_d_prev = self.rank_differential_q(degree-1)
        else:
            rank_d_prev = 0
        
        return C_dim - rank_d_i - rank_d_prev
    
    def print_free_rank(self):
        for degree in range(self.n+1):
            print(f"rank H_{degree} = {self.get_free_rank(degree)}")
    

K = KhovanovComplex("bCdA")
# K = KhovanovComplex("bfjihgaedc") # stress test

# for state in range(1 << K.n):
#     for crossing in range(K.n):
#         if ((state >> crossing) & 1) == 0:
#             bits = format(state, f"0{K.n}b")
#             print(bits, crossing, K.get_edge_map_info(state, crossing))


# for state in range(1 << K.n):
#     for crossing in range(K.n):
#         if ((state >> crossing) & 1) == 0:
#             for basis in K.get_basis_for_state(state):
#                 print(
#                     basis,
#                     "->",
#                     K.apply_edge_map(state, crossing, basis[1])
#                 )


# for degree in range(K.n):
#     M, rows, cols = K.differential_sparse(degree)

#     print(f"d_{degree}: {rows} x {cols}")
#     print(M)

for degree in range(K.n):
    print(f"rank d_{degree} =", K.rank_differential_q(degree))

K.print_free_rank()
import regina
from cube import Cube
from kalgebra import multiply, comultiply
from utils import crossingport_to_vertexidx
from fractions import Fraction

class KhovanovComplex:

    def __init__(self, code):
        
        self.cube = Cube(code)
        self.knot = self.cube.get_knot()
        self.crossings = self.cube.get_crossings()
        
        self.state_circles = {} # {key: state, value: [circles]}
        self.state_basis_vectors = {} # {key: state, value: [labeling tuples]}
        self.chain_basis = {} # {key: (i,j), value: [(state, labeling tuple)]}
        self.quantum_gradings_in_degree = {} # {key: i, value: {j}}
        self.vertexidx_to_statecircle = {} 
        # vertexidx_to_statecircle[state][vertexidx] = circleidx
        # so vertexidx belongs to circle: state_circles[state][circleidx]

        self.positive_crossings = []
        self.negative_crossings = []

        self.setup()

    # ==== SETUP / HELPER FUNCTIONS ====

    def setup(self):
        """
        Setup for the cube that we build our chain complex on
        """
        # initialize postive and negative crossings
        for crossingidx in range(self.crossings):
            if self.knot.crossing(crossingidx).sign() == 1:
                self.positive_crossings.append(crossingidx)
            else:
                self.negative_crossings.append(crossingidx)

        self.fill_state_circles()
        self.fill_state_basis()
        self.fill_chain_basis()
    
    def homological_grading(self, state):
        return state.bit_count() - len(self.negative_crossings)
    
    def quantum_grading(self, state, labeling):
        vminus = sum(labeling) # number of X's (X has degree -1)
        vplus = len(labeling) - vminus # number of 1's is just the other circs
        return state.bit_count() + vplus - vminus + len(self.positive_crossings) - 2*len(self.negative_crossings)

    def fill_state_circles(self):
        for state in range(2**self.crossings):
            self.vertexidx_to_statecircle[state] = {}
            state_circles = self.cube.get_circles(state)
            self.state_circles[state] = state_circles
            for circleidx, circle in enumerate(state_circles):
                for vertexidx in circle:
                    self.vertexidx_to_statecircle[state][vertexidx] = circleidx

    def fill_state_basis(self):
        for state in range(2**self.crossings):
            circle_count = len(self.state_circles[state])
            # can label each circle with 1 (rep'd by 0) or X (rep'd by 1)
            labelings = []
            for labels in range(2**circle_count):
                labels_tuple = tuple((labels>>i) & 1 for i in range(circle_count))
                labelings.append(labels_tuple)
            self.state_basis_vectors[state] = labelings

    def fill_chain_basis(self):
        """
        Returns list of all basis vectors with state info (labeling tuples) from homological grading i and quantum grading j.
        """
        for state in range(2**self.crossings):
            for labeling in self.state_basis_vectors[state]:
                i = self.homological_grading(state)
                j = self.quantum_grading(state, labeling)
                if i not in self.quantum_gradings_in_degree:
                    self.quantum_gradings_in_degree[i] = set()
                self.quantum_gradings_in_degree[i].add(j)
                if (i, j) not in self.chain_basis:
                    self.chain_basis[(i,j)] = []
                self.chain_basis[(i,j)].append((state, labeling))
                
    # ==== EDGE LOGIC ====

    def get_edge_map_info(self, state, crossingidx):
        """
        Edge is srcstate -> deststate s.t. crossing at crossingidx changes 0->1.
        
        Returns:
            ("merge", src_circ_a, src_circ_b, dest_circ)
        or
            ("split", src_circ, dest_circ_a, dest_circ_b)
        """

        srcstate = state
        deststate = srcstate | (1<<crossingidx)

        critical_vertices = [
            crossingport_to_vertexidx(crossingidx, port)
            for port in range(4)
        ]

        src_circs = {
            self.vertexidx_to_statecircle[srcstate][v]
            for v in critical_vertices
        }
        dest_circs = {
            self.vertexidx_to_statecircle[deststate][v]
            for v in critical_vertices
        }

        if len(src_circs) == 2 and len(dest_circs) == 1:
            # two circles merged into one
            a, b = sorted(src_circs)
            dest = next(iter(dest_circs))
            return ("merge", a, b, dest)
        elif len(src_circs) == 1 and len(dest_circs) == 2:
            # one circle split into two
            a, b = sorted(dest_circs)
            src = next(iter(src_circs))
            return ("split", src, a, b)
                
    def apply_edge_map(self, state, labeling, crossingidx):
        """
        Applies edge map from srcstate -> deststate s.t. 
        crossing at crossingidx changes 0->1.

        Returns:
            (deststate, [(coeff, dest_labeling)])
        """

        srcstate = state
        deststate = srcstate | (1<<crossingidx)
        info = self.get_edge_map_info(state, crossingidx)

        kind = info[0]
        results = []

        if kind == "merge":

            _, src_a, src_b, dest = info
            dest_labeling = [None] * len(self.state_circles[deststate])

            # handle unchanged circles
            for srcCircidx, srcCirc in enumerate(self.state_circles[srcstate]):
                if srcCircidx == src_a or srcCircidx == src_b:
                    continue
                representative_vertexidx = srcCirc[0]
                destCircidx = self.vertexidx_to_statecircle[deststate][representative_vertexidx]
                dest_labeling[destCircidx] = labeling[srcCircidx]

            label_a = labeling[src_a]
            label_b = labeling[src_b]

            for res_label in multiply(label_a, label_b):
                dest_labeling[dest] = res_label
                results.append(tuple(dest_labeling))

        elif kind == "split":

            _, src, dest_a, dest_b = info
            dest_labeling = [None] * len(self.state_circles[deststate])

            # handle unchanged circles
            for srcCircidx, srcCirc in enumerate(self.state_circles[srcstate]):
                if srcCircidx == src:
                    continue
                representative_vertexidx = srcCirc[0]
                destCircidx = self.vertexidx_to_statecircle[deststate][representative_vertexidx]
                dest_labeling[destCircidx] = labeling[srcCircidx]

            label = labeling[src]

            for res_label_a, res_label_b in comultiply(label):
                dest_labeling[dest_a] = res_label_a
                dest_labeling[dest_b] = res_label_b
                results.append(tuple(dest_labeling))

        return results

    def edge_sign(self, state, crossingidx):
        """
        To make anticommutative cube faces: (-1)^{#1-bits before crossingidx}
        """
        count = 0
        for k in range(crossingidx):
            count += (state >> k) & 1
        return -1 if count%2 else 1
    
    # ==== SPARSE MATRIX ====

    def build_differential_sparse_matrix(self, i, j):
        """
        Builds d: C^{i,j} -> C^{i+1,j} over Q.

        Returns:
            M, cSrc, cDest
            M[col] = {key : row, value : coeff}
            cSrc = C^{i,j} chain basis
            cDest = C^{i+1,j} chain basis
        """
        
        cSrc = self.chain_basis.get((i,j), [])
        cDest = self.chain_basis.get((i+1,j), [])

        # literally cuz iirc map lookup is considerably faster than list
        cDestMap = {
            basis_vector : idx 
            for idx, basis_vector in enumerate(cDest)
        }

        M = []

        for state, labeling in cSrc:
            col = {}
            for crossingidx in range(self.crossings):
                if ((state>>crossingidx) & 1) == 1:
                    continue
                coeff = Fraction(self.edge_sign(state, crossingidx))
                deststate = state | (1<<crossingidx)
                for dest_labeling in self.apply_edge_map(state, labeling, crossingidx):
                    basis_vector = (deststate, dest_labeling)
                    if basis_vector not in cDestMap:
                        continue
                    row = cDestMap[basis_vector]
                    col[row] = col.get(row, Fraction(0)) + coeff
                    if col[row] == 0:
                        del col[row] # sparsify!
            M.append(col)

        return M, cSrc, cDest
    
    # ==== LINALG | GAUSSIAN ELIMINATION ====

    def get_differential_rank(self, M):
        """
        M is a sparse matrix, where entries are fractions
        M[col][row] = coeff

        Proceed via gaussian elimination over Q.
        """

        pivots = {}

        for col in M:
            col = dict(col) # copy construct so we don't mess with og
            while col:
                pivot_row = max(col.keys())
                pivot_coeff = col[pivot_row]
                if pivot_row not in pivots:
                    # found a new independent
                    inv = Fraction(1,1) / pivot_coeff
                    for r in list(col.keys()):
                        col[r] *= inv
                        if col[r] == 0:
                            del col[r] # sparsify!
                    pivots[pivot_row] = col
                    break
                pivot_col = pivots[pivot_row]
                factor = col[pivot_row]
                for r, val in pivot_col.items():
                    col[r] = col.get(r, Fraction(0)) - factor * val
                    if col[r] == 0:
                        del col[r] # sparsify!
        
        return len(pivots)
    
    def get_differential_generators(self, M):
        """
        M is a sparse matrix, where entries are fractions
        M[col][row] = coeff

        Proceed via gaussian elimination over Q.
        
        Returns:
            rank, pivots, pivot_generators, kernel_generators
        """

        pivots = {}
        
        pivot_generators = {}
        # pivot_generators[pivot_row] = {key: colidx of M, val: coeff}
        kernel_generators = []
        # kernel_generators = [{key: colidx of M, val: coeff}]
        # combining all entries in each list element gives a kernel element

        for colidx, col in enumerate(M):
            col = dict(col) # copy construct so we don't mess with og
            gen = {colidx : Fraction(1)} # before any elim, col is just og basis vector at M[colidx]
            while col:
                pivot_row = max(col.keys())
                pivot_coeff = col[pivot_row]
                if pivot_row not in pivots:
                    # found a new independent
                    inv = Fraction(1,1) / pivot_coeff
                    for r in list(col.keys()):
                        col[r] *= inv
                        if col[r] == 0:
                            del col[r] # sparsify!
                    for k in list(gen.keys()):
                        gen[k] *= inv
                        if gen[k] == 0:
                            del gen[k] # sparsify!
                    pivots[pivot_row] = col
                    pivot_generators[pivot_row] = gen
                    break
                
                pivot_col = pivots[pivot_row]
                pivot_gen = pivot_generators[pivot_row]
                factor = col[pivot_row]
                
                for r, val in pivot_col.items():
                    col[r] = col.get(r, Fraction(0)) - factor * val
                    if col[r] == 0:
                        del col[r] # sparsify!
                for k, val in pivot_gen.items():
                    gen[k] = gen.get(k, Fraction(0)) - factor * val
                    if gen[k] == 0:
                        del gen[k] # sparsify!
            
            # if col reduced to zero, then this combo is in kernel
            if not col:
                kernel_generators.append(gen)
        
        return len(pivots), pivots, pivot_generators, kernel_generators
    
    def quotient_cycles_by_boundaries(self, cycles, boundaries):
        """
        cycles: list of sparse vectors in C^{i,j}
        boundaries: list of sparse vectors in C^{i,j}

        Returns homology representatives (generators): cycles mod boundaries.
        """

        pivots = {}
        survivors = []

        # First put all boundaries into the pivot table
        for b in boundaries:
            v = dict(b) # copy so we don't mess with og
            while v:
                pivot = max(v.keys())
                coeff = v[pivot]
                if pivot not in pivots:
                    inv = Fraction(1,1) / coeff
                    for k in list(v.keys()):
                        v[k] *= inv
                        if v[k] == 0:
                            del v[k] # sparsify!
                    pivots[pivot] = v
                    break
                for k, val in pivots[pivot]:
                    v[k] = v.get(k, Fraction(0)) + (-coeff)*val
                    if v[k] == 0:
                        del v[k] # sparsify!
        
        # Next reduce cycles modulo those boundaries
        for c in cycles:
            v = dict(c) # copy so we don't mess with og
            while v:
                pivot = max(v.keys())
                coeff = v[pivot]
                if pivot not in pivots:
                    # this cycle is NOT a boundary, so it SURVIVES
                    inv = Fraction(1,1) / coeff
                    for k in list(v.keys()):
                        v[k] *= inv
                        if v[k] == 0:
                            del v[k] # sparsify!
                    pivots[pivot] = v
                    survivors.append(v)
                    break
                for k, val in pivots[pivot].items():
                    v[k] = v.get(k, Fraction(0)) + (-coeff)*val
                    if v[k] == 0:
                        del v[k] # sparsify!
        
        return survivors

    # ==== HOMOLOGY CALCULATION ====
    
    def homology_rank(self, i, j):
        Cij_dim = len(self.chain_basis.get((i, j), []))
        d_cur, _, _ = self.build_differential_sparse_matrix(i, j)
        d_prev, _, _ = self.build_differential_sparse_matrix(i - 1, j)
        rank_d_cur = self.get_differential_rank(d_cur)
        rank_d_prev = self.get_differential_rank(d_prev)

        return Cij_dim - rank_d_cur - rank_d_prev
    
    def homology_generators(self, i, j):
        
        d_cur, cCur, cNext = self.build_differential_sparse_matrix(i,j)
        d_prev, cPrev, _ = self.build_differential_sparse_matrix(i-1,j)

        _, _, _, cycles = self.get_differential_generators(d_cur)
        _, boundary_cols, _, _ = self.get_differential_generators(d_prev)

        boundaries = list(boundary_cols.values())

        survivors = self.quotient_cycles_by_boundaries(cycles, boundaries)
        return survivors, cCur
    
    def print_homology(self, i, j):
        r = self.homology_rank(i,j)
        if r != 0:
            print(f"rank Kh^({i},{j}) = {r}")

    def print_homology_generators(self, i, j):
        survivors, cCur = self.homology_generators(i,j)
        if len(survivors) == 0:
            return
        print(f"Kh^({i},{j}) has {len(survivors)} generators")
        for genidx, gen in enumerate(survivors):
            print(f"generator {genidx}:")
            for basis_idx, coeff in sorted(gen.items()):
                state, labeling = cCur[basis_idx]
                bits = format(state, f"0{self.crossings}b")
                prettylabeling = [None]*len(labeling)
                for i, val in enumerate(labeling):
                    if val == 0:
                        prettylabeling[i] = '1'
                    else:
                        prettylabeling[i] = 'X'
                print(f"  {coeff} * state={bits}, labeling={prettylabeling}")
        print()

    def print_full_homology(self):
        for i, js in sorted(self.quantum_gradings_in_degree.items()):
            for j in sorted(js):
                self.print_homology(i, j)
    
    def print_full_homology_generators(self):
        for i, js in sorted(self.quantum_gradings_in_degree.items()):
            for j in sorted(js):
                self.print_homology_generators(i, j)

# Usage
# K = KhovanovComplex("eabcdbadcvbZa")
# K.print_full_homology()
# print("===========================")
# K.print_full_homology_generators()
import regina
from utils import *

class DSU:

    def __init__(self, n):
        self.n = n
        self.parents = []
        self.sizes = []
        for i in range(n):
            self.parents.append(i)
            self.sizes.append(1)
    
    def find(self, x):
        if self.parents[x] != x:
            self.parents[x] = self.find(self.parents[x])
        return self.parents[x]
    
    def unite(self, x, y):
        x_root = self.find(x)
        y_root = self.find(y)
        if x_root == y_root:
            return False
        if self.sizes[x_root] < self.sizes[y_root]:
            x_root, y_root = y_root, x_root
        self.sizes[x_root] += self.sizes[y_root]
        self.parents[y_root] = x_root
        return True
    
    def connected(self, x, y):
        return self.find(x) == self.find(y)
    
    def count_components(self):
        roots = {self.find(i) for i in range(self.n)}
        return len(roots)
    
    def get_components(self):
        comps = {}
        for i in range(self.n):
            root = self.find(i)
            if root not in comps:
                comps[root] = []
            comps[root].append(i)
        return [comps[key] for key in comps]

class Cube:

    def __init__(self, code):
        self.knot = get_knot(code)
        self.crossings = self.knot.size()
        self.n = 4*self.crossings

    # ==== helper functions ====
    
    def get_knot(self):
        return self.knot
    
    # ==== the good stuff ====
    
    def build_fresh_graph(self):
        """
        Builds and returns a fresh graph for the knot.
        Graph has no resolutions. Note this means the crossing endpoints are all disjoint, so this is NOT a picture of the knot.

        Returns edge list for graph.
        """
        edges = []
        for crossingidx in range(self.crossings):
            
            crossing = self.knot.crossing(crossingidx)
            
            # handle the upper strand
            # connect the outgoing (next) upper strand at this crossing
            # to wherever it is incoming (prev) at its next crossing
            upperstrand = crossing.upper()
            strand_to_crossing_and_type
            uppersrc = strand_endpoint_to_vertexidx(upperstrand, "next")
            upperdest = strand_endpoint_to_vertexidx(upperstrand.next(), "prev")
            edges.append((uppersrc, upperdest))
            
            # handle the lower strand
            # connect the outgoing (next) lower strand at this crossing
            # to wherever it is incoming (prev) at its next crossing
            lowerstrand = crossing.lower()
            lowersrc = strand_endpoint_to_vertexidx(lowerstrand, "next")
            lowerdest = strand_endpoint_to_vertexidx(lowerstrand.next(), "prev")
            edges.append((lowersrc, lowerdest))

        return edges
    
    def smooth_crossing(self, crossingidx, restype):
        """
        crossingidx: crossing idx of the crossing to resolve
        restype: resolution type. 0/1 value for 0 or 1 resolution respectively.

        Returns:
            list of edges to add to graph to simulate desired smoothing
        """
        
        p0 = crossingport_to_vertexidx(crossingidx, 0)
        p1 = crossingport_to_vertexidx(crossingidx, 1)
        p2 = crossingport_to_vertexidx(crossingidx, 2)
        p3 = crossingport_to_vertexidx(crossingidx, 3)

        if self.knot.crossing(crossingidx).sign() == 1:
            # positive crossing -- oriented res is zero res
            if restype == 0:
                return [(p0, p3), (p1, p2)]
            else:
                return [(p0, p2), (p1, p3)]
        else:
            # negative crossing -- oriented res is one res
            if restype == 1:
                return [(p0, p3), (p1, p2)]
            else:
                return [(p0, p2), (p1, p3)]
            
    def build_state(self, state):
        dsu = DSU(self.n)
        edges = self.build_fresh_graph()
        for a, b in edges:
            dsu.unite(a, b)
        for crossingidx in range(self.crossings):
            restype = (state>>crossingidx) & 1
            res_edges = self.smooth_crossing(crossingidx, restype)
            for a, b in res_edges:
                dsu.unite(a, b)
        return dsu
    
    def count_circles(self, state):
        dsu = self.build_state(state)
        return dsu.count_components()
    
    def get_circles(self, state):
        dsu = self.build_state(state)
        return dsu.get_components()

    def print_state(self, state):
        bits = format(state, f"0{self.crossings}b")
        circles = self.get_circles(state)
        print(bits, len(circles), circles)

    def print_cube(self):
        states = [i for i in range(2**self.crossings)]
        sorted_states = sorted(states, key=lambda x: x.bit_count())
        for state in sorted_states:
            self.print_state(state)
            print()

K = Cube("eabcdbadcvbZa")
K.print_cube()
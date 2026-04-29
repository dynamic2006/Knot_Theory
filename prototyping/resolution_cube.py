import regina

class DSU:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra != rb:
            self.parent[rb] = ra

    def count_components(self):
        return len({self.find(i) for i in range(len(self.parent))})


class ResolutionCube:
    def __init__(self, dtcode: str):
        self.knot = regina.Link.fromDT(dtcode)
        self.n = self.knot.size()
        self.num_vertices = 4 * self.n

    def port_id(self, crossing_index, port):
        """
        crossing i has four ports:

            0 = lower.prev side
            1 = lower.next side
            2 = upper.prev side
            3 = upper.next side
        """
        return 4 * crossing_index + port

    def strand_to_crossing_and_type(self, s):
        """
        Return (crossing_index, kind), where kind is '^' or '_'.
        """
        text = str(s)
        kind = text[0]
        idx = int(text[1:])
        return idx, kind

    def endpoint_port_id(self, s, direction):
        """
        Convert a Regina strand endpoint into the corresponding 4-port vertex.

        direction is either:
            "prev" or "next"

        For crossing i:
            _i.prev  -> port 0
            _i.next  -> port 1
            ^i.prev  -> port 2
            ^i.next  -> port 3
        """
        idx, kind = self.strand_to_crossing_and_type(s)

        if kind == "_" and direction == "prev":
            return self.port_id(idx, 0)
        if kind == "_" and direction == "next":
            return self.port_id(idx, 1)
        if kind == "^" and direction == "prev":
            return self.port_id(idx, 2)
        if kind == "^" and direction == "next":
            return self.port_id(idx, 3)

    def build_graph(self):
        """
        Add in the original edges to build the graph
        For each local crossing, consider each strand (over/under). Connect:
        - current strand's next-side port  -- current strand's prev-side port at next crossing
        - note: we don't actually connect the prev and next sides of each strand
        why? because DSU can't break edges! so just leave this part for DSU
        """
        edges = []

        for i in range(self.n):
            cross = self.knot.crossing(i)

            for strand in [cross.upper(), cross.lower()]:
                start = self.endpoint_port_id(strand, "next")
                end = self.endpoint_port_id(strand.next(), "prev")
                edges.append((start, end))

        return edges

    def smoothing_edges(self, crossing_index, bit):
        """
        Internal smoothing edges at one crossing.

        Ports:
            0 = lower.prev
            1 = lower.next
            2 = upper.prev
            3 = upper.next

        0-resolution:
            lower.prev -- upper.next
            upper.prev -- lower.next

        1-resolution:
            lower.prev -- upper.prev
            lower.next -- upper.next
        """
        p0 = self.port_id(crossing_index, 0)
        p1 = self.port_id(crossing_index, 1)
        p2 = self.port_id(crossing_index, 2)
        p3 = self.port_id(crossing_index, 3)

        if bit == 0:
            return [(p0, p3), (p2, p1)]
        else:
            return [(p0, p2), (p1, p3)]

    def circle_count(self, state):
        dsu = DSU(self.num_vertices)

        for a, b in self.build_graph():
            dsu.union(a, b)

        for i in range(self.n):
            bit = (state >> i) & 1

            for a, b in self.smoothing_edges(i, bit):
                dsu.union(a, b)

        return dsu.count_components()

    def print_cube(self):
        for state in range(2 ** self.n):
            bits = format(state, f"0{self.n}b")
            print(bits, self.circle_count(state))


K = ResolutionCube("BCA")
K.print_cube()
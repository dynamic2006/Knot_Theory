import regina
from config import *

def get_knot(code):
    if CODETYPE == "sig":
        return regina.Link.fromSig(code)
    elif CODETYPE == "dt":
        return regina.Link.fromDT(code)
    else:
        raise ValueError("Unknown CODETYPE")
    
def crossingport_to_vertexidx(crossingidx, port):
    """
    Crossing i has four ports:
        0 = lower.prev side
        1 = lower.next side
        2 = upper.prev side
        3 = upper.next side
    """
    return 4*crossingidx + port

def strand_to_crossing_and_type(strand):
    """
    Parses strand object
    Input:
        strand: regina strand object. Ex: ^0
    Returns:
        (crossingidx, kind). Ex: (0, ^)
    """
    text = str(strand)
    kind = text[0]
    crossingidx = int(text[1:])
    return crossingidx, kind

def strand_endpoint_to_vertexidx(strand, direction):
    """
    Convert regina strand endpoint => vertex idx

    direction: "prev" or "next"

    For crossing i:
        _i.prev  -> port 0
        _i.next  -> port 1
        ^i.prev  -> port 2
        ^i.next  -> port 3
    """
    crossingidx, kind = strand_to_crossing_and_type(strand)

    if kind == "_" and direction == "prev":
        return crossingport_to_vertexidx(crossingidx, 0)
    if kind == "_" and direction == "next":
        return crossingport_to_vertexidx(crossingidx, 1)
    if kind == "^" and direction == "prev":
        return crossingport_to_vertexidx(crossingidx, 2)
    if kind == "^" and direction == "next":
        return crossingport_to_vertexidx(crossingidx, 3)
def multiply(a, b):
    if a == b:
        return [0]
    else:
        return [1]

def comultiply(a):
    if a == 0:
        return [(0,1), 
                (1,0)]
    else:
        return [(1,1),
                (0,0)]
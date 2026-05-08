def multiply(a, b):
    if a == 0 and b == 0:
        return [0]
    elif a == 0 and b == 1:
        return [1]
    elif a == 1 and b == 0:
        return [1]
    else:
        return []

def comultiply(a):
    if a == 0:
        return [(1,0), 
                (0,1)]
    else:
        return [(1,1)]
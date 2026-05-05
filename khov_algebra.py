def convert(label: str):
    if label == '1':
        return 0
    elif label == 'X':
        return 1

def multiply(a, b):
    if a == '1' and b == '1':
        return [(1, convert('1'))]
    elif a == '1' and b == 'X':
        return [(1, convert('X'))]
    elif a == 'X' and b == '1':
        return [(1, convert('X'))]
    else:
        return []

def comultiply(a):
    if a == '1':
        return [(1, (convert('X'), convert('1'))), 
                (1, (convert('1'), convert('X')))]
    else:
        return [(1, (convert('X'), convert('X')))]
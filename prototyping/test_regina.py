import regina

K = regina.Link.fromDT("bca")
print(K)
print(K.size())

c = K.crossing(0)
print(dir(c))

for i in range(K.size()):
    c = K.crossing(i)
    print(i,c)

    u = c.upper()
    l = c.lower()

    print(u)

    print("INFORMATION")
    print("upper next:", u.next())
    print("upper prev:", u.prev())
    print("lower next:", l.next())
    print("lower prev:", l.prev())
    print("==========")
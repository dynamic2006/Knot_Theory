import regina
from cube import ResolutionCube

knot_sigs = ["gabcadebcfdefvvpd"]

# with open("../docs/knots.txt", "r") as f:
#     for line in f:
#         parts = line.split()

#         if len(parts) >= 2:
#             knot_sigs.append(parts[1])

# print(knot_sigs)

best_sig = ""
for sig in knot_sigs:
    K_zero = regina.Link.fromSig(sig)
    K_oriented = regina.Link.fromSig(sig)

    # zero resolve the K_zero knot
    zero_cube = ResolutionCube(sig)
    zero_cube_count = zero_cube.count_circles(0)
    # orient resolve the K_oriented knot
    
    oriented_state = 0
    for i in range(K_oriented.size()-1, -1, -1):
        print(i, K_oriented.crossing(i))
        oriented_state = oriented_state << 1
        if K_oriented.crossing(i).sign() != 1:
            oriented_state+=1
    
    oriented_cube = ResolutionCube(sig)
    oriented_cube_count = oriented_cube.count_circles(oriented_state)

    if(zero_cube_count != oriented_cube_count):
        bits = format(oriented_state, f"0{6}b")
        print(bits)
        print("ZERO RES: " + str(zero_cube_count))
        print("ORIENTED RES: " + str(oriented_cube_count))
        if(len(sig) < len(best_sig) or best_sig==""):
            best_sig = sig

print(best_sig)
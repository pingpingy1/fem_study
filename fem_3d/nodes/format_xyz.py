filename = "../IGZO.xyz"

coords = []

with open(filename, "r") as f:
    N = int(f.readline().strip())
    f.readline()
    for _ in range(N):
        coords.append(tuple(map(float, f.readline().strip().split()[1:])))

with open("nodes.node", "w") as f:
    f.write(f"{len(coords)} 3 0 0\n")
    for i, p in enumerate(coords):
        f.write(f"{i + 1} {p[0]} {p[1]} {p[2]}\n")

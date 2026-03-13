N = 49
a = 1.0 / (N + 1)
nodes = []
boundary = []

for i in range(1, N + 1):
    for j in range(1, N + 1):
        nodes.append((a * i, a * j))

for i in range(N + 2):
    boundary.append((0.0, a * i))
    boundary.append((1.0, a * i))
    if i != 0 and i != N + 1:
        boundary.append((a * i, 0.0))
        boundary.append((a * i, 1.0))

with open("nodes.txt", "w") as f:
    f.write(f"{len(nodes)} {len(boundary)}\n")
    for i, p in enumerate(nodes):
        f.write(f"{p[0]} {p[1]}\n")
    for i, p in enumerate(boundary):
        f.write(f"{p[0]} {p[1]}\n")


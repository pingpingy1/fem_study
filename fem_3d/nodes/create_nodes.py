N = 19
a = 1.0 / (N + 1)
nodes = []
boundary = [
    (0.0, 0.0, 0.0),
    (0.0, 0.0, 1.0),
    (0.0, 1.0, 0.0),
    (0.0, 1.0, 1.0),
    (1.0, 0.0, 0.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 0.0),
    (1.0, 1.0, 1.0),
]

for i in range(1, N + 1):
    for j in range(1, N + 1):
        for k in range(1, N + 1):
            nodes.append((a * i, a * j, a * k))

with open("nodes.node", "w") as f:
    f.write(f"{len(boundary) + len(nodes)} 3 0 0\n")
    for i, p in enumerate(boundary):
        f.write(f"{i+1} {p[0]} {p[1]} {p[2]}\n")
    for i, p in enumerate(nodes):
        f.write(f"{i+len(boundary)+1} {p[0]} {p[1]} {p[2]}\n")

N = 24
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

with open("nodes.poly", "w") as f:
    f.write(f"{len(boundary) + len(nodes)} 3 0 1\n")
    for i, p in enumerate(boundary):
        f.write(f"{i+1} {p[0]} {p[1]} {p[2]} 1\n")
    for i, p in enumerate(nodes):
        f.write(f"{i+len(boundary)+1} {p[0]} {p[1]} {p[2]} 0\n")
    f.write("""
6 0
1
4 1 2 4 3
1
4 1 2 6 5
1
4 1 3 7 5
1
4 2 4 8 6
1
4 3 4 8 7
1
4 5 6 8 7

0""")

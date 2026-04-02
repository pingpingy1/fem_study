from math import ceil
import os

TEST = False

filename = "../test.xyz" if TEST else "../a-IGZO.xyz"

coords = []

L = 1 if TEST else 10.98
W = 1 if TEST else 10.98
H = 1 if TEST else 10.98

with open(filename, "r") as f:
    N = int(f.readline().strip())
    f.readline()
    for _ in range(N):
        point = tuple(map(float, f.readline().strip().split()[1:]))
        coords.append((
            point[0] % L,
            point[1] % W,
            point[2] % H,
        ))

"""
L = int(ceil(max(x for (x, _, _) in coords)))
W = int(ceil(max(y for (_, y, _) in coords)))
H = int(ceil(max(z for (_, _, z) in coords)))
"""

print("Dimensions: " + str(L) + " * " + str(W) + " * " + str(H))

# ==================================
# Triangulate Interfaces
# ==================================
"""
To enforce periodic (Bloch) boundary conditions,
we algorithmically triangulate the unit cell boundary
with a simple rectangular grid.
"""

SURFACE_RES = 0.5

def generate_boundary_mesh(L, W, H, Nx, Ny, Nz):
    nodes = []
    node_id = {}
    def add_node(x, y, z):
        key = (x, y, z)
        if key not in node_id:
            node_id[key] = len(nodes)
            nodes.append(key)
        return node_id[key]

    faces = []

    def triangulate_face(origin, du, dv, Nu, Nv):
        """
        origin: starting poiint
        du, dv: direction vectors
        Nu, Nv: number of divisions
        """
        face = []
        for i in range(int(Nu)):
            for j in range(int(Nv)):
                p00 = (origin[0] +  i      * du[0] +  j      * dv[0],
                       origin[1] +  i      * du[1] +  j      * dv[1],
                       origin[2] +  i      * du[2] +  j      * dv[2])
                p10 = (origin[0] + (i + 1) * du[0] +  j      * dv[0],
                       origin[1] + (i + 1) * du[1] +  j      * dv[1],
                       origin[2] + (i + 1) * du[2] +  j      * dv[2])
                p01 = (origin[0] +  i      * du[0] + (j + 1) * dv[0],
                       origin[1] +  i      * du[1] + (j + 1) * dv[1],
                       origin[2] +  i      * du[2] + (j + 1) * dv[2])
                p11 = (origin[0] + (i + 1) * du[0] + (j + 1) * dv[0],
                       origin[1] + (i + 1) * du[1] + (j + 1) * dv[1],
                       origin[2] + (i + 1) * du[2] + (j + 1) * dv[2])

                n00 = add_node(*p00)
                n01 = add_node(*p01)
                n10 = add_node(*p10)
                n11 = add_node(*p11)

                if (i + j) % 2 == 0:
                    face.append((n00, n10, n11))
                    face.append((n00, n11, n01))
                else:
                    face.append((n00, n10, n01))
                    face.append((n10, n11, n01))
        return face

    dx = L / Nx
    dy = W / Ny
    dz = H / Nz

    # x = 0
    faces.append(triangulate_face(
        origin=(0, 0, 0),
        du=(0, dy, 0),
        dv=(0, 0, dz),
        Nu=Ny,
        Nv=Nz,
    ))

    # x = L
    faces.append(triangulate_face(
        origin=(L, 0, 0),
        du=(0, dy, 0),
        dv=(0, 0, dz),
        Nu=Ny,
        Nv=Nz,
    ))

    # y = 0
    faces.append(triangulate_face(
        origin=(0, 0, 0),
        du=(dx, 0, 0),
        dv=(0, 0, dz),
        Nu=Nx,
        Nv=Nz,
    ))

    # y = W
    faces.append(triangulate_face(
        origin=(0, W, 0),
        du=(dx, 0, 0),
        dv=(0, 0, dz),
        Nu=Nx,
        Nv=Nz,
    ))

    # z = 0
    faces.append(triangulate_face(
        origin=(0, 0, 0),
        du=(dx, 0, 0),
        dv=(0, dy, 0),
        Nu=Nx,
        Nv=Ny,
    ))

    # z = H
    faces.append(triangulate_face(
        origin=(0, 0, H),
        du=(dx, 0, 0),
        dv=(0, dy, 0),
        Nu=Nx,
        Nv=Ny,
    ))

    return nodes, faces


surface_nodes, surface_mesh = generate_boundary_mesh(
    L, W, H, ceil(L / SURFACE_RES), ceil(W / SURFACE_RES), ceil(H / SURFACE_RES)
)


# =================================
# Build Wire Mesh
# =================================

with open("nodes.poly", "w") as f:
    f.write("{} 3 0 1\n".format(len(coords) + len(surface_nodes)))

    for i, (x, y, z) in enumerate(surface_nodes):
        f.write("{} {} {} {} 1\n".format(i + 1, x, y, z))
    for i, (x, y, z) in enumerate(coords):
        f.write("{} {} {} {} 0\n".format(i + 1 + len(surface_nodes), x, y, z))

    f.write("\n6 0\n")

    for i in range(6):
        f.write("{}\n".format(len(surface_mesh[i])))
        for j, (p1, p2, p3) in enumerate(surface_mesh[i]):
            f.write("3 {} {} {}\n".format(p1 + 1, p2 + 1, p3 + 1))

    f.write("\n0")

os.system("./tetgen -pqYkV nodes.poly")

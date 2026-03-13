import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.linalg import eigh

x = np.loadtxt("nodes.txt", skiprows=1)
H = np.loadtxt("H.txt")
M = np.loadtxt("M.txt")

E, psi = eigh(H, M)
print(E.shape)
print(psi.shape)

# ===================================

def load_triangles(filename):
    tri_list = []

    with open(filename, "r") as f:
        N_tri = int(f.readline())

        for line in f:
            parts = line.strip().split()
            assert len(parts) == 3
            indices = map(lambda s: int(s) - 1, parts)
            tri_list.append(list(indices))
    
    return np.array(tri_list)

tri_indices = load_triangles("triangles.txt")
triangulation = tri.Triangulation(
    x[:,0],
    x[:,1],
    triangles=tri_indices,
)

fig, axes = plt.subplots(1, 5, figsize=(20, 4))
fig.suptitle(f"First 5 Eigenfunctions")
cmap = "RdBu_r"

psi_full = np.vstack([psi, np.zeros((len(x) - len(E), len(E)))])

for i in range(5):
    ax = axes[i]
    z = psi_full[:,i]

    limit = np.max(np.abs(z))
    levels = np.linspace(-limit, limit, 50)

    cf = ax.tricontourf(triangulation, z, levels=levels, cmap=cmap)

    ax.set_title(f"E[{i}] = {E[i]:.4f}")
    ax.set_aspect("equal")
    ax.axis("off")
    plt.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig("solution.png")

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

x = np.loadtxt("nodes.txt", skiprows=1)
H = np.loadtxt("H.txt")
M = np.loadtxt("M.txt")

E, psi = eigh(H, M)
# print(E.shape)
# print(psi.shape)

n_plot = min(5, len(E))
plt.figure(figsize=(6,4))

for i in range(n_plot):
    plt.plot(
        x,
        psi[:,i],
        marker="o",
        label=f"n={i+1}, E={E[i]:.4f}",
    )

plt.xlabel("x")
plt.ylabel("psi(x)")
plt.title("Eigenstates")
plt.legend()
plt.grid(True)

plt.savefig("solution.png")

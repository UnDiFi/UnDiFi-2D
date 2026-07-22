# Compares su2_capturing's converged solution against this directory's
# own iteration-501 final state (na00501.1.node/.ele, sh99_final.dat --
# copied here from the scratch run directory the actual 501-iteration
# run used; see README.md). Run from this directory.
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

GAM = 1.4

HERE = os.path.dirname(os.path.abspath(__file__))
BASE_CAP = os.path.join(HERE, "..", "su2_capturing")
BASE_FIT = HERE
FIT_NODE = os.path.join(HERE, "na00501.1.node")
FIT_ELE = os.path.join(HERE, "na00501.1.ele")
FIT_SHOCK = os.path.join(HERE, "sh99_final.dat")


def mach_from_prim(rho, u, v, p):
    a = np.sqrt(GAM * p / rho)
    return np.sqrt(u * u + v * v) / a


def read_su2_mesh(path):
    with open(path) as f:
        lines = f.read().split("\n")
    idx = 0
    assert lines[idx].startswith("NDIME")
    idx += 1
    nelem = int(lines[idx].split("=")[1])
    idx += 1
    tris = []
    for _ in range(nelem):
        parts = lines[idx].split()
        idx += 1
        tris.append((int(parts[1]), int(parts[2]), int(parts[3])))
    assert lines[idx].startswith("NPOIN")
    npoin = int(lines[idx].split("=")[1])
    idx += 1
    xy = np.zeros((npoin, 2))
    for _ in range(npoin):
        parts = lines[idx].split()
        idx += 1
        pid = int(parts[2])
        xy[pid, 0] = float(parts[0])
        xy[pid, 1] = float(parts[1])
    return xy, np.array(tris, dtype=int)


def load_capturing():
    xy, tris = read_su2_mesh(f"{BASE_CAP}/na00.1.su2")
    data = np.loadtxt(f"{BASE_CAP}/na00.1_restart_converged.csv", delimiter=",", skiprows=1)
    pid = data[:, 0].astype(int)
    rho = data[:, 3]
    u = data[:, 4] / rho
    v = data[:, 5] / rho
    rhoE = data[:, 6]
    p = (GAM - 1) * (rhoE - 0.5 * rho * (u * u + v * v))
    order = np.argsort(pid)
    M = np.zeros_like(rho)
    M[pid[order]] = mach_from_prim(rho[order], u[order], v[order], p[order])
    return xy, tris, M


def load_fitting():
    with open(f"{BASE_FIT}/na00501.1.node") as f:
        header = f.readline().split()
        npoin = int(header[0])
        xy = np.zeros((npoin + 1, 2))
        M = np.zeros(npoin + 1)
        for _ in range(npoin):
            parts = f.readline().split()
            i = int(parts[0])
            x, y = float(parts[1]), float(parts[2])
            z1, z2, z3, z4 = (float(parts[3]), float(parts[4]),
                              float(parts[5]), float(parts[6]))
            rho = z1 * z1
            H = z2 / z1
            u = z3 / z1
            v = z4 / z1
            q2 = u * u + v * v
            p = rho * (GAM - 1) / GAM * (H - 0.5 * q2)
            xy[i] = (x, y)
            M[i] = mach_from_prim(rho, u, v, p)
    with open(f"{BASE_FIT}/na00501.1.ele") as f:
        nelem = int(f.readline().split()[0])
        tris = np.zeros((nelem, 3), dtype=int)
        for k in range(nelem):
            parts = f.readline().split()
            tris[k] = (int(parts[1]), int(parts[2]), int(parts[3]))
    return xy, tris, M


def load_shock(path):
    with open(path) as f:
        lines = f.read().split("\n")
    npts = int(lines[1].split()[0])
    pts = np.array([list(map(float, lines[2 + k].split()))[:2] for k in range(npts)])
    return pts


xy_c, tri_c, M_c = load_capturing()
xy_f, tri_f, M_f = load_fitting()
shock = load_shock(FIT_SHOCK)

vmax = 3.0
levels = np.linspace(0, vmax, 41)

fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharex=True, sharey=True)

for ax, xy, tri, M, title in [
    (axes[0], xy_c, tri_c, M_c, "Shock-capturing (SU2, converged)"),
    (axes[1], xy_f, tri_f, M_f, "Shock-fitting (SU2, iteration 501)"),
]:
    triang = mtri.Triangulation(xy[:, 0], xy[:, 1], tri)
    cf = ax.tricontourf(triang, M, levels=levels, cmap="turbo", extend="max")
    ax.tricontour(triang, M, levels=[1.0], colors="k", linewidths=0.8,
                   linestyles="--")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_aspect("equal")
    ax.set_xlim(-0.3, 2.2)
    ax.set_ylim(-2.2, 2.2)
    wall = plt.Circle((0, 0), 1.0, fill=True, color="0.15", zorder=5)
    ax.add_patch(wall)

axes[1].plot(shock[:, 0], shock[:, 1], "-", color="white", lw=1.6, zorder=6,
             label="fitted shock front")
axes[1].legend(loc="upper left", fontsize=8, framealpha=0.85)
axes[0].set_ylabel("y")

cbar = fig.colorbar(cf, ax=axes, orientation="vertical", fraction=0.035, pad=0.02)
cbar.set_label("Mach number (capped at 3; freestream is M=20)")

fig.suptitle("CircularCylinder-1, M$_\\infty$=20 -- SU2 shock-capturing vs shock-fitting\n"
              "dashed line: sonic (M=1)")
out = os.path.join(HERE, "capturing_vs_fitting.png")
fig.savefig(out, dpi=170, bbox_inches="tight")
print("wrote", out)
print("Mach range capturing:", M_c.min(), M_c.max())
print("Mach range fitting:", M_f.min(), M_f.max())

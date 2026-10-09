"""Plots of the 2D template mesh (overall C-domain and close-ups)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.collections import PolyCollection

d = np.load("template2d.npz")
nodes = d["nodes"]
views = [("Domain (C-mesh)", (-6.5, 9.6, -6.5, 6.5)),
         ("Airfoil", (-0.05, 0.42, -0.12, 0.17)),
         ("Leading edge", (-0.012, 0.02, -0.012, 0.018)),
         ("Trailing edge", (0.285, 0.315, -0.012, 0.01))]
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
for ax, (title, (x0, x1, y0, y1)) in zip(axs.flat, views):
    for cells, col in ((d["quads"], "tab:blue"), (d["ftri"], "tab:orange")):
        p = nodes[cells]
        keep = ((p[:, :, 0].max(1) > x0) & (p[:, :, 0].min(1) < x1)
                & (p[:, :, 1].max(1) > y0) & (p[:, :, 1].min(1) < y1))
        ax.add_collection(PolyCollection(p[keep], facecolor="none",
                                         edgecolor=col, lw=0.25))
    ax.fill(d["foil"][:, 0], d["foil"][:, 1], color="0.6")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
    ax.set_title(title); ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")
fig.suptitle("2D section mesh: quads (blue, 5 cm wall layer + wake strip), triangles (orange)")
fig.tight_layout()
os.makedirs("results", exist_ok=True)
fig.savefig("results/mesh_2d_section.png", dpi=150)

"""Figure for Section 'The split branch on the tetraquadric': the five charge coordinates of the fixed datum
with the possibly present singlets as arrows, and the eleven-datum verdict."""
import json, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

R = json.load(open("yukawa_charges_results.json"))
fig = plt.figure(figsize=(11, 4.6))
ax = fig.add_subplot(1, 2, 1)
ax.set_aspect("equal"); ax.axis("off")
pos = {i: (math.cos(2 * math.pi * i / 5 + math.pi / 2), math.sin(2 * math.pi * i / 5 + math.pi / 2)) for i in range(5)}
bcl = R[0]
sources = {p for (p, q) in bcl["singlets"]}
for i in range(5):
    x, y = pos[i]
    col = "#c0392b" if i not in sources else "#2c3e50"
    ax.add_patch(plt.Circle((x, y), 0.17, color=col, zorder=3))
    ax.text(x, y, f"$e_{i+1}$", color="white", ha="center", va="center", fontsize=13, zorder=4)
for (p, q) in bcl["singlets"]:
    (x1, y1), (x2, y2) = pos[p], pos[q]
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=16, shrinkA=14, shrinkB=14,
                                 connectionstyle="arc3,rad=0.15", color="#7f8c8d", lw=1.4, zorder=2))
lab = {1: "10, 5bar, 5", 4: "10, 5bar, 5", 3: "5bar", 0: "5bar", 2: "no matter"}
for i, t in lab.items():
    x, y = pos[i]; r = 1.42
    ax.text(r * x, r * y, t.replace("5bar", r"$\overline{\mathbf{5}}$").replace("10", r"$\mathbf{10}$").replace(", 5", r", $\mathbf{5}$"),
            ha="center", va="center", fontsize=10)
ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.6, 1.65)
ax.set_title("Fixed datum: singlet charges $e_p-e_q$ as arrows $p\\to q$\nno arrow leaves $e_1$ or $e_3$", fontsize=10)

ax2 = fig.add_subplot(1, 2, 2); ax2.axis("off")
names = ["(fixed)"] + [str(i) for i in range(1, 11)]
cols = ["down", "up"]
for j, c in enumerate(cols):
    ax2.text(j + 1.0, 11.4, c, ha="center", va="center", fontsize=11, weight="bold")
for i, o in enumerate(R):
    y = 10.5 - i
    ax2.text(0.0, y, names[i], ha="center", va="center", fontsize=10)
    for j, key in enumerate(["down_allowed", "up_allowed"]):
        ok = o[key] > 0
        ax2.add_patch(plt.Rectangle((j + 0.6, y - 0.4), 0.8, 0.8, color="#27ae60" if ok else "#c0392b", alpha=0.85))
        ax2.text(j + 1.0, y, "allowed" if ok else "forbidden", ha="center", va="center", fontsize=8, color="white")
ax2.set_xlim(-0.5, 2.8); ax2.set_ylim(-0.3, 12)
ax2.set_title("Eleven split data: perturbative selection rules\n(all possibly present singlets inserted)", fontsize=10)
fig.tight_layout()
import os, sys
out = sys.argv[1] if len(sys.argv) > 1 else "."
fig.savefig(os.path.join(out, "fig_charges.png"), dpi=200)
print("ok")

"""Regenerate the figures in docs/figures."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from corrugated_buckling import Sheet, trapezoidal, equivalent_plate, analytical, frame_fe, finite_strip
from corrugated_buckling.dissertation import SHEET, LENGTH, WIDTH, PROFILES, ANSYS, REPORTED

OUT = pathlib.Path(__file__).resolve().parents[1] / "docs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
BLUE, RED, GREY, SAND = "#2A6FB0", "#C4472D", "#456072", "#D9A441"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 150})

# 1. Buckling modes across the corrugations
fig, axes = plt.subplots(3, 2, figsize=(10, 5.2), sharex=True)
for row, (name, p) in enumerate(PROFILES.items()):
    for col, ends in enumerate(["pinned", "fixed"]):
        N, xz, mode = frame_fe.buckle(SHEET, p, LENGTH, ends)
        if mode[np.abs(mode[:, 1]).argmax(), 1] < 0: mode = -mode
        ax = axes[row, col]
        ax.plot(xz[:, 0], xz[:, 1], color=GREY, lw=1, alpha=0.5)
        ax.plot(xz[:, 0] + 0.08 * mode[:, 0], xz[:, 1] + 0.08 * mode[:, 1], color=RED, lw=1.6)
        ax.set_aspect("equal"); ax.set_ylim(-0.07, 0.15); ax.grid(False); ax.set_yticks([])
        ax.set_title(f"{name}, {ends} ends: {N/1e3:.1f} kN/m", fontsize=10, loc="left")
axes[2, 0].set_xlabel("x across the corrugations [m]"); axes[2, 1].set_xlabel("x across the corrugations [m]")
fig.tight_layout(); fig.savefig(OUT / "modes_across.png"); plt.close(fig)

# 2. Reported, corrected, numerical and Ansys values
fig, ax = plt.subplots(figsize=(10, 4.2))
keys = list(ANSYS); x = np.arange(len(keys)); w = 0.2
rep = [REPORTED[k] for k in keys]; cor = []; num = []
for name, d, ends in keys:
    p = PROFILES[name]; st = equivalent_plate(SHEET, p)
    if d == "x":
        cor.append(analytical.wide_strut(st.D11, LENGTH, ends)); num.append(frame_fe.buckle(SHEET, p, LENGTH, ends)[0])
    else:
        cor.append(analytical.wide_strut(st.D22, LENGTH, ends)); num.append(finite_strip.panel_load(SHEET, p, WIDTH, LENGTH, ends, "free")[0])
for i, (vals, lab, colr) in enumerate([(rep, "dissertation, analytical", GREY), (cor, "corrected analytical", BLUE),
                                       (num, "frame FE / finite strip", RED), ([ANSYS[k] for k in keys], "dissertation, Ansys", SAND)]):
    ax.bar(x + (i - 1.5) * w, np.array(vals) / 1e3, w, label=lab, color=colr)
ax.set_yscale("log"); ax.set_ylabel("critical load [kN/m]")
ax.set_xticks(x); ax.set_xticklabels([f"{n[:4]}\n{d}, {e}" for n, d, e in keys], fontsize=8)
ax.legend(ncol=4, fontsize=8, frameon=False, loc="upper left"); ax.set_ylim(10, 1e6)
fig.tight_layout(); fig.savefig(OUT / "comparison.png"); plt.close(fig)

# 3. Finite strip signature curves along the corrugations
hw = np.geomspace(0.01, 4, 60)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
for ax, t in zip(axes, [0.00635, 0.0007]):
    sheet = Sheet(t=t)
    for name, colr in [("trapezoidal", BLUE), ("round", RED)]:
        p = PROFILES[name]
        ax.loglog(hw, finite_strip.signature_curve(sheet, p, WIDTH, hw) / 1e3, color=colr, label=name)
        ax.axhline(analytical.squash_load(sheet, p) / 1e3, color=colr, ls="--", lw=1)
    ax.axvline(LENGTH, color=GREY, lw=0.8, ls=":")
    ax.set_title(f"t = {t*1e3:g} mm (dashed: yield of the section)", fontsize=10, loc="left")
    ax.set_xlabel("buckle half-wavelength [m]")
axes[0].set_ylabel("critical load [kN/m]"); axes[0].legend(frameon=False)
fig.tight_layout(); fig.savefig(OUT / "signature_curves.png"); plt.close(fig)

# 4. Effect of corrugation depth
f = np.linspace(0.001, 0.0245, 40); Nx, Ny, mass = [], [], []
for fi in f:
    p = trapezoidal(fi, 0.0508, np.pi / 4); st = equivalent_plate(SHEET, p)
    Nx.append(analytical.wide_strut(st.D11, LENGTH)); Ny.append(analytical.wide_strut(st.D22, LENGTH))
    mass.append(analytical.mass_per_area(SHEET, p))
fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.semilogy(f * 2e3, np.array(Nx) / 1e3, color=BLUE, label="across the corrugations")
ax.semilogy(f * 2e3, np.array(Ny) / 1e3, color=RED, label="along the corrugations")
ax.axhline(analytical.wide_strut(SHEET.D11, LENGTH) / 1e3, color=GREY, ls=":", label="flat plate")
ax.set_xlabel("corrugation depth 2f [mm]"); ax.set_ylabel("Euler load, pinned ends [kN/m]"); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(OUT / "depth_study.png"); plt.close(fig)
print("figures written to", OUT)

"""Recompute every case of the dissertation and compare with what it reported."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from corrugated_buckling import equivalent_plate, analytical, frame_fe, finite_strip
from corrugated_buckling.dissertation import SHEET, LENGTH, WIDTH, PROFILES, ANSYS, REPORTED

print("Equivalent plate stiffness (A in N/m, D in N m)")
for name, p in PROFILES.items():
    st = equivalent_plate(SHEET, p)
    print(f"  {name:12s} " + "  ".join(f"{k}={v:.4g}" for k, v in st.as_dict().items()))

hdr = f"{'panel':12s} {'dir':3s} {'ends':7s} {'reported':>11s} {'corrected':>11s} {'numerical':>11s} {'Ansys 2023':>11s} {'squash':>10s}"
print("\nCritical line loads in kN/m\n" + hdr + "\n" + "-" * len(hdr))
for (name, d, ends), ansys in ANSYS.items():
    p = PROFILES[name]; st = equivalent_plate(SHEET, p)
    if d == "x":
        corrected = analytical.wide_strut(st.D11, LENGTH, ends)
        numerical = frame_fe.buckle(SHEET, p, LENGTH, ends)[0]
    else:
        corrected = analytical.wide_strut(st.D22, LENGTH, ends)
        numerical = finite_strip.panel_load(SHEET, p, WIDTH, LENGTH, ends, "free")[0]
    squash = analytical.squash_load(SHEET, p, d)
    print(f"{name:12s} {d:3s} {ends:7s} {REPORTED[(name, d, ends)]/1e3:11.1f} {corrected/1e3:11.1f} "
          f"{numerical/1e3:11.1f} {ansys/1e3:11.1f} {squash/1e3:10.1f}")
print("\nnumerical = frame finite elements (x) or finite strips with free unloaded edges (y)")

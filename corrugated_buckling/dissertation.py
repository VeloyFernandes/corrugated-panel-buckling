"""Reference case and the Ansys eigenvalue results reported in the 2023 dissertation."""
import numpy as np
from .homogenise import Sheet
from .profiles import flat, trapezoidal, round_

SHEET = Sheet(E=2.1e11, nu=0.3, t=0.00635)
LENGTH = WIDTH = 1.016                       # panel side [m]
PROFILES = {
    "flat": flat(),
    "trapezoidal": trapezoidal(f=0.0127, c=0.0508, theta=np.pi / 4),
    "round": round_(R=0.0254, L=0.0254),
}

# Ansys Workbench eigenvalue buckling loads [N/m]; keys are (panel, direction, ends)
ANSYS = {
    ("flat", "x", "pinned"): 47075, ("flat", "x", "fixed"): 188250,
    ("flat", "y", "pinned"): 47075, ("flat", "y", "fixed"): 188250,
    ("trapezoidal", "x", "pinned"): 38843, ("trapezoidal", "x", "fixed"): 155760,
    ("trapezoidal", "y", "pinned"): 154670, ("trapezoidal", "y", "fixed"): 154670,
    ("round", "x", "pinned"): 19642, ("round", "x", "fixed"): 78398,
    ("round", "y", "pinned"): 567070, ("round", "y", "fixed"): 567070,
}

# Analytical values printed in the dissertation [N/m]
REPORTED = {
    ("flat", "x", "pinned"): 47079.3179, ("flat", "x", "fixed"): 188317.2715,
    ("flat", "y", "pinned"): 47079.3179, ("flat", "y", "fixed"): 188317.2715,
    ("trapezoidal", "x", "pinned"): 39001.7839, ("trapezoidal", "x", "fixed"): 156007.1357,
    ("trapezoidal", "y", "pinned"): 1702730.0189, ("trapezoidal", "y", "fixed"): 6810920.0758,
    ("round", "x", "pinned"): 18313.1263, ("round", "x", "fixed"): 73252.5053,
    ("round", "y", "pinned"): 105188.3902, ("round", "y", "fixed"): 420753.5608,
}

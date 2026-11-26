import numpy as np
import pytest
from corrugated_buckling import Sheet, flat, trapezoidal, round_, equivalent_plate, analytical, frame_fe, finite_strip

S = Sheet(); A = 1.016
TRAP = trapezoidal(0.0127, 0.0508, np.pi / 4); RND = round_(0.0254, 0.0254)


def test_profile_integrals_match_polyline():
    for p in (TRAP, round_(0.0254, 0.0254, n_arc=200)):
        d = np.diff(p.points, axis=0); ds = np.hypot(*d.T); zm = (p.points[:-1, 1] + p.points[1:, 1]) / 2
        assert np.isclose(ds.sum(), 2 * p.l, rtol=1e-3)
        assert np.isclose((d[:, 0] ** 2 / ds).sum(), p.I1, rtol=1e-3)
        assert np.isclose((zm ** 2 * ds).sum() + (d[:, 1] ** 2 * ds).sum() / 12, p.I2, rtol=1e-3)


def test_dissertation_stiffness_values():
    st = equivalent_plate(S, TRAP)
    assert st.D11 == pytest.approx(4079.2, rel=1e-4)
    assert st.D22 == pytest.approx(1.7809e5, rel=1e-4)
    assert equivalent_plate(S, RND).D11 == pytest.approx(1915.4, rel=1e-4)


def test_flat_limit_recovers_isotropic_plate():
    st = equivalent_plate(S, flat())
    assert st.D11 == pytest.approx(S.D11) and st.D22 == pytest.approx(S.D11)
    assert st.A11 == pytest.approx(S.A11) and st.A66 == pytest.approx(S.A66)


@pytest.mark.parametrize("ends", ["pinned", "fixed", "pinned-fixed", "fixed-free"])
@pytest.mark.parametrize("p", [flat(), TRAP, RND])
def test_frame_fe_matches_homogenised_strut(p, ends):
    st = equivalent_plate(S, p)
    assert frame_fe.buckle(S, p, A, ends)[0] == pytest.approx(analytical.wide_strut(st.D11, A, ends), rel=5e-3)


def test_finite_strip_square_plate_k4():
    N = finite_strip.buckle(S, flat(), A, A, "supported")[0]
    assert N == pytest.approx(4 * np.pi ** 2 * S.D11 / A ** 2, rel=1e-3)


def test_finite_strip_global_mode_along_corrugations():
    st = equivalent_plate(S, TRAP)
    N = finite_strip.buckle(S, TRAP, A, A, "supported")[0]
    assert N == pytest.approx(analytical.wide_strut(st.D22, A), rel=0.03)

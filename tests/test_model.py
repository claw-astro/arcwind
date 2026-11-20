import numpy as np
import pytest

from diskwind import ConicalWindModel, ViewGeometry
from diskwind.geometry import cone_cut_heights
from diskwind.laws import VelocityLaw, available_laws, get_law, register_law

G = ViewGeometry(inclination_deg=49.0, s=-1, distance_pc=137.5, vsys=1.9)


def test_cut_heights_lie_on_the_cut():
    phi = np.linspace(0, 2 * np.pi, 50)
    y0 = 1.0 * G.distance_pc
    z, den = cone_cut_heights(y0, 20.0, 0.4, phi, G)
    R = 20.0 + 0.4 * z
    Y = z * G.sini + G.s * R * np.sin(phi) * G.cosi
    assert np.allclose(Y[den > 0], y0)


def test_edge_on_curve_is_ellipse():
    g = ViewGeometry(90.0, 1, 100.0, 0.0)
    m = ConicalWindModel("constant", radial="proportional", R0_au=50, tan_theta=0.0,
                         Vz=5, Vr=3, j=100)
    c = m.curve(1.0, g)
    ok = c.valid
    # v = Vr sin(phi) + (j/R) cos(phi): ellipse in (X, v)
    assert np.allclose(c.v[ok], 3 * np.sin(c.phi[ok]) + 2 * np.cos(c.phi[ok]))


def test_broken_power_law_values():
    law = get_law("broken-power-law", delta=0.1)
    p = dict(Vz_b=10.0, z_b=300.0, q1=0.2, q2=1.0)
    assert np.isclose(law([300.0], p)[0], 10.0)
    a, b, c, d = law([3.0, 30.0, 3000.0, 30000.0], p)
    assert np.isclose(np.log10(b / a), 0.2, atol=0.01)      # slope q1 well below z_b
    assert np.isclose(np.log10(d / c), 1.0, atol=0.01)      # slope q2 well above z_b


def test_radial_modes_and_names():
    assert "Vr" not in ConicalWindModel("power-law", radial="cone").names
    assert "Vr" in ConicalWindModel("power-law", radial="proportional").names
    m = ConicalWindModel("broken-power-law", radial="independent")
    assert {"Vr_b", "r1", "r2"} <= set(m.names) and m.names.count("z_b") == 1
    m.set(tan_theta=0.5)
    Vz, Vr = ConicalWindModel("power-law", radial="cone", tan_theta=0.5, Vz=4).velocities([100])
    assert np.isclose(Vr[0], 2.0)


def test_fix_free_bounds_and_draws():
    m = ConicalWindModel("power-law")
    m.fix(xslope=0.0).set_bounds(Vz=(0, 20))
    assert "xslope" not in m.free_names and m.bounds["Vz"] == (0, 20)
    m.free("xslope")
    assert "xslope" in m.free_names
    draws = np.tile(m.vector(), (5, 1))
    Vz, Vr = m.velocities([50, 100, 200], draws)
    assert Vz.shape == (5, 3)
    with pytest.raises(KeyError):
        m.set(not_a_param=1)


def test_custom_law_registration():
    @register_law
    class Linear(VelocityLaw):
        name = "linear-test"
        param_names = ("V0", "dVdz")
        amplitude = "V0"
        radial_map = {"V0": "Vr0"}

        def evaluate(self, z, p):
            return p["V0"] + p["dVdz"] * z

        def default_bounds(self, z_cuts):
            return {"V0": (0, 50), "dVdz": (-1, 1)}

        def guess(self, v_ref, z_cuts):
            return {"V0": v_ref, "dVdz": 0.0}

    assert "linear-test" in available_laws()
    m = ConicalWindModel("linear-test", V0=5, dVdz=0.01)
    assert np.isclose(m.velocities([100])[0][0], 6.0)
    c = m.curve(1.0, G)
    assert c.valid.any()

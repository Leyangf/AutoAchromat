"""
test_thermal.py – Unit tests for the thermal analysis module.

All numerical reference values are derived from real SCHOTT AGF data
for N-BK7 and N-SF5.  See plan section 1.2 for worked calculations.
"""
from __future__ import annotations

import math
import tempfile
import os

import pytest

from autoachromat.glass_reader import Glass, read_agf
from autoachromat.thermal import (
    ThermalMetrics,
    dn_dT,
    thermo_optical_coeff,
    system_thermal_power_derivative,
    required_housing_cte,
    thermal_defocus,
    compute_thermal_metrics,
)

# ---------------------------------------------------------------------------
# Real-data fixtures (SCHOTT AGF values, verified in plan §1.2)
# ---------------------------------------------------------------------------

LAM_D = 0.58756   # d-line [µm]

# N-BK7: nd=1.5168, ν=64.17, CTE=7.1 ppm/K
# TD: D0=1.86e-6, E0=4.34e-7, λ_tk=0.170
BK7 = Glass(
    name="N-BK7",
    catalog="SCHOTT",
    formula_id=2,
    cd=[
        1.03961212, 6.00069867e-3, 0.231792344,
        2.00179144e-2, 1.01046945, 103.560653,
    ],
    cte_m40_20=7.1,
    td_D0=1.86e-6,
    td_D1=1.31e-8,
    td_D2=-1.37e-11,
    td_E0=4.34e-7,
    td_E1=6.27e-10,
    td_ltk=0.170,
    td_T_ref=20.0,
)

# N-SF5: nd=1.67271, ν=32.25, CTE=7.94 ppm/K
# TD: D0=-2.51e-7, E0=7.85e-7, λ_tk=0.278
SF5 = Glass(
    name="N-SF5",
    catalog="SCHOTT",
    formula_id=2,
    cd=[
        1.52481889, 1.1254756e-2, 0.187085527,
        5.88995392e-2, 1.42729015, 129.141675,
    ],
    cte_m40_20=7.94,
    td_D0=-2.51e-7,
    td_D1=1.07e-8,
    td_D2=-2.40e-11,
    td_E0=7.85e-7,
    td_E1=1.15e-9,
    td_ltk=0.278,
    td_T_ref=20.0,
)

N_BK7 = 1.5168   # n_d for N-BK7
N_SF5  = 1.67271  # n_d for N-SF5

# Achromatic normalised powers for BK7+SF5 doublet
# phi1 = ν1/(ν1-ν2), phi2 = -ν2/(ν1-ν2)
_NU1, _NU2 = 64.17, 32.25
_DELTA_NU = _NU1 - _NU2
PHI1 = _NU1 / _DELTA_NU   # ≈ +2.010
PHI2 = -_NU2 / _DELTA_NU  # ≈ -1.010


# ---------------------------------------------------------------------------
# TestDnDT
# ---------------------------------------------------------------------------

class TestDnDT:
    def test_bk7_value(self):
        """N-BK7 @ d-line: dn/dT ≈ 1.39e-6 /K (plan §1.2)."""
        result = dn_dT(BK7, LAM_D, N_BK7)
        assert result is not None
        assert abs(result - 1.39e-6) < 0.05e-6   # within 5% of reference

    def test_sf5_value(self):
        """N-SF5 @ d-line: dn/dT ≈ 1.44e-6 /K."""
        result = dn_dT(SF5, LAM_D, N_SF5)
        assert result is not None
        assert abs(result - 1.44e-6) < 0.05e-6

    def test_missing_D0_returns_none(self):
        g = Glass(name="X", catalog="TEST", td_E0=4.34e-7, td_ltk=0.170)
        assert dn_dT(g, LAM_D, 1.5) is None

    def test_missing_E0_returns_none(self):
        g = Glass(name="X", catalog="TEST", td_D0=1.86e-6, td_ltk=0.170)
        assert dn_dT(g, LAM_D, 1.5) is None

    def test_missing_ltk_returns_none(self):
        g = Glass(name="X", catalog="TEST", td_D0=1.86e-6, td_E0=4.34e-7)
        assert dn_dT(g, LAM_D, 1.5) is None

    def test_wavelength_dependence(self):
        """dn/dT should increase as wavelength approaches λ_tk."""
        r_far = dn_dT(BK7, 0.8, N_BK7)
        r_near = dn_dT(BK7, 0.25, N_BK7)
        assert r_far is not None and r_near is not None
        # Closer to λ_tk=0.170 → larger E₀/(λ²-λ_tk²) term
        assert abs(r_near) > abs(r_far)


# ---------------------------------------------------------------------------
# TestThermoOpticalCoeff
# ---------------------------------------------------------------------------

class TestThermoOpticalCoeff:
    def test_bk7_negative(self):
        """N-BK7: V should be negative (≈ −4.4 ppm/K)."""
        V = thermo_optical_coeff(BK7, LAM_D, N_BK7)
        assert V is not None
        assert V < 0
        assert abs(V * 1e6 - (-4.4)) < 0.5   # within 0.5 ppm/K

    def test_sf5_negative(self):
        """N-SF5: V should be negative (≈ −5.8 ppm/K)."""
        V = thermo_optical_coeff(SF5, LAM_D, N_SF5)
        assert V is not None
        assert V < 0
        assert abs(V * 1e6 - (-5.8)) < 0.5

    def test_missing_cte_returns_none(self):
        g = Glass(
            name="X", catalog="TEST",
            td_D0=1.86e-6, td_E0=4.34e-7, td_ltk=0.170,
        )  # no cte_m40_20
        assert thermo_optical_coeff(g, LAM_D, 1.5) is None

    def test_missing_td_returns_none(self):
        g = Glass(name="X", catalog="TEST", cte_m40_20=7.1)  # no TD
        assert thermo_optical_coeff(g, LAM_D, 1.5) is None

    def test_identity(self):
        """V = dn/dT/(n-1) - alpha should be self-consistent."""
        V = thermo_optical_coeff(BK7, LAM_D, N_BK7)
        dndt = dn_dT(BK7, LAM_D, N_BK7)
        assert V is not None and dndt is not None
        alpha = BK7.cte_m40_20 * 1e-6
        expected = dndt / (N_BK7 - 1.0) - alpha
        assert abs(V - expected) < 1e-15


# ---------------------------------------------------------------------------
# TestSystemThermalPowerDerivative
# ---------------------------------------------------------------------------

class TestSystemThermalPowerDerivative:
    def test_bk7_sf5_doublet_negative(self):
        """BK7+SF5 achromat: dΦ/dT_norm should be negative."""
        V1 = thermo_optical_coeff(BK7, LAM_D, N_BK7)
        V2 = thermo_optical_coeff(SF5, LAM_D, N_SF5)
        dphi = system_thermal_power_derivative(V1, PHI1, V2, PHI2)
        assert dphi is not None
        assert dphi < 0

    def test_linear_combination(self):
        """dΦ/dT_norm = V1*phi1 + V2*phi2 by definition."""
        V1 = thermo_optical_coeff(BK7, LAM_D, N_BK7)
        V2 = thermo_optical_coeff(SF5, LAM_D, N_SF5)
        dphi = system_thermal_power_derivative(V1, PHI1, V2, PHI2)
        assert dphi is not None
        assert abs(dphi - (V1 * PHI1 + V2 * PHI2)) < 1e-20

    def test_none_V1_propagates(self):
        assert system_thermal_power_derivative(None, PHI1, -5e-6, PHI2) is None

    def test_none_V2_propagates(self):
        assert system_thermal_power_derivative(-4.4e-6, PHI1, None, PHI2) is None

    def test_both_none_propagates(self):
        assert system_thermal_power_derivative(None, PHI1, None, PHI2) is None


# ---------------------------------------------------------------------------
# TestRequiredHousingCte
# ---------------------------------------------------------------------------

class TestRequiredHousingCte:
    def test_positive_for_typical_doublet(self):
        """α_h,req = −dΦ/dT_norm > 0 for BK7+SF5."""
        V1 = thermo_optical_coeff(BK7, LAM_D, N_BK7)
        V2 = thermo_optical_coeff(SF5, LAM_D, N_SF5)
        dphi = system_thermal_power_derivative(V1, PHI1, V2, PHI2)
        alpha_req = required_housing_cte(dphi)
        assert alpha_req is not None
        assert alpha_req > 0

    def test_realistic_range(self):
        """BK7+SF5: α_h,req should be in (1, 10) ppm/K."""
        V1 = thermo_optical_coeff(BK7, LAM_D, N_BK7)
        V2 = thermo_optical_coeff(SF5, LAM_D, N_SF5)
        dphi = system_thermal_power_derivative(V1, PHI1, V2, PHI2)
        alpha_req = required_housing_cte(dphi)
        assert alpha_req is not None
        alpha_ppm = alpha_req * 1e6
        assert 1.0 < alpha_ppm < 10.0

    def test_negation_of_dphi(self):
        dphi = -3.0e-6
        assert abs(required_housing_cte(dphi) - 3.0e-6) < 1e-20

    def test_none_propagates(self):
        assert required_housing_cte(None) is None


# ---------------------------------------------------------------------------
# TestThermalDefocus
# ---------------------------------------------------------------------------

class TestThermalDefocus:
    def test_zero_when_matched(self):
        """δf' = 0 when α_h_actual = α_h_required."""
        alpha_req = 5.0e-6
        df = thermal_defocus(200.0, alpha_req, alpha_req, 20.0)
        assert df is not None
        assert abs(df) < 1e-15

    def test_aluminum_negative(self):
        """Aluminium (23 ppm/K) > α_h,req → δf' < 0."""
        alpha_req = 3.0e-6
        df = thermal_defocus(200.0, alpha_req, 23.0e-6, 20.0)
        assert df is not None
        assert df < 0

    def test_linear_in_delta_T(self):
        alpha_req = 3.0e-6
        df1 = thermal_defocus(200.0, alpha_req, 10.0e-6, 10.0)
        df2 = thermal_defocus(200.0, alpha_req, 10.0e-6, 20.0)
        assert df1 is not None and df2 is not None
        assert abs(df2 - 2 * df1) < 1e-12

    def test_formula_explicit(self):
        """δf' = f' × (α_req − α_actual) × ΔT."""
        f = 150.0
        alpha_req = 4e-6
        alpha_act = 11e-6
        dT = 30.0
        expected = f * (alpha_req - alpha_act) * dT
        result = thermal_defocus(f, alpha_req, alpha_act, dT)
        assert result is not None
        assert abs(result - expected) < 1e-12

    def test_none_propagates(self):
        assert thermal_defocus(200.0, None, 10.0e-6, 20.0) is None


# ---------------------------------------------------------------------------
# TestComputeThermalMetrics
# ---------------------------------------------------------------------------

class TestComputeThermalMetrics:
    def test_full_data(self):
        m = compute_thermal_metrics(BK7, SF5, N_BK7, N_SF5, PHI1, PHI2, LAM_D)
        assert m.thermal_data_available is True
        assert m.V1 is not None and m.V1 < 0
        assert m.V2 is not None and m.V2 < 0
        assert m.dphi_dT_norm is not None and m.dphi_dT_norm < 0
        assert m.alpha_housing_required is not None and m.alpha_housing_required > 0

    def test_missing_td(self):
        g_no_td = Glass(name="X", catalog="TEST", cte_m40_20=7.1)
        m = compute_thermal_metrics(g_no_td, SF5, 1.5, N_SF5, PHI1, PHI2, LAM_D)
        assert m.thermal_data_available is False
        assert m.V1 is None

    def test_missing_cte(self):
        g_no_cte = Glass(
            name="X", catalog="TEST",
            td_D0=1.86e-6, td_E0=4.34e-7, td_ltk=0.170,
        )
        m = compute_thermal_metrics(g_no_cte, SF5, 1.5, N_SF5, PHI1, PHI2, LAM_D)
        assert m.thermal_data_available is False
        assert m.V1 is None

    def test_never_raises(self):
        """compute_thermal_metrics must not raise for any input combination."""
        bad = Glass(name="BAD", catalog="TEST")
        m = compute_thermal_metrics(bad, bad, 0.0, 0.0, 0.0, 0.0, LAM_D)
        assert isinstance(m, ThermalMetrics)

    def test_alpha_req_consistent(self):
        """alpha_housing_required == −dphi_dT_norm."""
        m = compute_thermal_metrics(BK7, SF5, N_BK7, N_SF5, PHI1, PHI2, LAM_D)
        assert m.dphi_dT_norm is not None and m.alpha_housing_required is not None
        assert abs(m.alpha_housing_required + m.dphi_dT_norm) < 1e-20


# ---------------------------------------------------------------------------
# TestAgfThermalParsing
# ---------------------------------------------------------------------------

class TestAgfThermalParsing:
    """Verify that glass_reader correctly parses ED CTE and TD constants."""

    _AGF_TEMPLATE = """\
CC Test catalog
NM N-BK7 2 517642.251 1.5168 64.17 0 1
ED 7.100000 8.300000 2.511000 -0.000900 0
CD 1.039612120E+00 6.000698670E-03 2.317923440E-01 2.001791440E-02 1.010469450E+00 1.035606530E+02
TD 1.860000E-06 1.310000E-08 -1.370000E-11 4.340000E-07 6.270000E-10 1.700000E-01 2.000000E+01
LD 3.00000E-01 2.50000E+00
NM N-SF5 2 673323.286 1.67271 32.25 0 1
ED 7.940000 9.210000 2.858000 0.008800 0
CD 1.524818890E+00 1.125475600E-02 1.870855270E-01 5.889953920E-02 1.427290150E+00 1.291416750E+02
TD -2.510000E-07 1.070000E-08 -2.400000E-11 7.850000E-07 1.150000E-09 2.780000E-01 2.000000E+01
LD 3.70000E-01 2.50000E+00
"""

    _AGF_NO_TD = """\
CC No-TD catalog
NM GLASS-A 2 500000.000 1.5 55.0 0 1
ED 8.0 9.0 2.5 0.001 0
CD 1.0 0.01 0.2 0.02 1.0 100.0
LD 0.35 2.5
"""

    def _write_agf(self, content: str) -> str:
        fd, path = tempfile.mkstemp(suffix=".agf")
        with os.fdopen(fd, "w") as f:
            f.write(content)
        return path

    def test_ed_cte_parsed(self):
        path = self._write_agf(self._AGF_TEMPLATE)
        try:
            _, glasses = read_agf(path)
        finally:
            os.unlink(path)
        bk7 = next(g for g in glasses if g.name == "N-BK7")
        assert bk7.cte_m40_20 == pytest.approx(7.1)
        assert bk7.cte_20_300 == pytest.approx(8.3)

    def test_td_constants_parsed(self):
        path = self._write_agf(self._AGF_TEMPLATE)
        try:
            _, glasses = read_agf(path)
        finally:
            os.unlink(path)
        bk7 = next(g for g in glasses if g.name == "N-BK7")
        assert bk7.td_D0 == pytest.approx(1.86e-6, rel=1e-4)
        assert bk7.td_E0 == pytest.approx(4.34e-7, rel=1e-4)
        assert bk7.td_ltk == pytest.approx(0.170, rel=1e-4)
        assert bk7.td_T_ref == pytest.approx(20.0)

    def test_sf5_td_negative_D0(self):
        path = self._write_agf(self._AGF_TEMPLATE)
        try:
            _, glasses = read_agf(path)
        finally:
            os.unlink(path)
        sf5 = next(g for g in glasses if g.name == "N-SF5")
        assert sf5.td_D0 == pytest.approx(-2.51e-7, rel=1e-4)

    def test_no_td_line_gives_none(self):
        path = self._write_agf(self._AGF_NO_TD)
        try:
            _, glasses = read_agf(path)
        finally:
            os.unlink(path)
        g = glasses[0]
        assert g.td_D0 is None
        assert g.td_E0 is None
        assert g.td_ltk is None
        # CTE should still be parsed from ED line
        assert g.cte_m40_20 == pytest.approx(8.0)

    def test_dpgf_still_parsed(self):
        """Existing dpgf parsing must be preserved after ED modification."""
        path = self._write_agf(self._AGF_TEMPLATE)
        try:
            _, glasses = read_agf(path)
        finally:
            os.unlink(path)
        bk7 = next(g for g in glasses if g.name == "N-BK7")
        assert bk7.dpgf == pytest.approx(-0.0009, rel=1e-3)

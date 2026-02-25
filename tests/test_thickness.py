"""
Unit tests for the thickness / diameter helpers in thickening.py.

Covers:
- Table 10-2 outside‐diameter lookup (retaining ring)
- Table 10-3 min edge / centre thickness lookup
- Signed sag with strict rejection
- element_thickness for positive, negative, and invalid lenses
"""

from __future__ import annotations

import math
import pytest

from autoachromat.thickening import (
    _sag,
    _lens_power,
    _piecewise_lookup,
    lookup_delta,
    lookup_t_edge_min,
    lookup_t_center_min,
    outside_diameter,
    element_thickness,
    correct_radii_for_thickness,
    _reconcile_cemented_radius,
    _SAG_MARGIN,
    _MIN_CENTER_ABSOLUTE,
)


# ===================================================================
# Table 10-2: outside diameter (retaining ring)
# ===================================================================


class TestLookupDelta:
    """lookup_delta(D, mount='retaining_ring')"""

    def test_D_less_than_6_returns_none(self):
        assert lookup_delta(5.0) is None
        assert lookup_delta(6.0) is None

    def test_D_8p5_returns_1(self):
        """Prompt example: D=8.5 → Δ=1.0"""
        assert lookup_delta(8.5) == 1.0

    @pytest.mark.parametrize(
        "D,expected",
        [
            (7.0, 1.0),
            (10.0, 1.0),
            (15.0, 1.5),
            (18.0, 1.5),
            (25.0, 2.0),
            (30.0, 2.0),
            (40.0, 2.5),
            (50.0, 2.5),
            (60.0, 3.0),
            (80.0, 3.0),
            (100.0, 3.5),
            (120.0, 3.5),
            (150.0, 4.5),
        ],
    )
    def test_piecewise_values(self, D, expected):
        assert lookup_delta(D) == expected

    def test_unsupported_mount_raises(self):
        with pytest.raises(ValueError, match="Unsupported mount"):
            lookup_delta(50.0, mount="glue")


class TestOutsideDiameter:
    def test_D_8p5(self):
        """Prompt example: D=8.5, retaining ring → φ=9.5"""
        assert outside_diameter(8.5) == pytest.approx(9.5)

    def test_D_50(self):
        assert outside_diameter(50.0) == pytest.approx(52.5)

    def test_D_too_small_returns_none(self):
        assert outside_diameter(5.0) is None


# ===================================================================
# Table 10-3: min edge / centre thickness
# ===================================================================


class TestLookupThicknessMinima:
    @pytest.mark.parametrize(
        "D,expected",
        [
            (5.0, 0.4),
            (8.0, 0.6),
            (12.0, 0.8),
            (25.0, 1.2),
            (50.0, 1.8),
            (70.0, 2.4),
            (100.0, 3.0),
            (200.0, 4.0),
        ],
    )
    def test_edge_min(self, D, expected):
        assert lookup_t_edge_min(D) == expected

    @pytest.mark.parametrize(
        "D,expected",
        [
            (5.0, 0.6),
            (8.0, 0.8),
            (12.0, 1.0),
            (25.0, 1.5),
            (50.0, 2.2),
            (70.0, 3.5),
            (100.0, 5.0),
            (200.0, 8.0),
        ],
    )
    def test_center_min(self, D, expected):
        assert lookup_t_center_min(D) == expected


# ===================================================================
# Signed sag
# ===================================================================


class TestSag:
    def test_flat_surface(self):
        """Infinite / zero radius → sag = 0."""
        assert _sag(float("inf"), 10.0) == 0.0
        assert _sag(0.0, 10.0) == 0.0

    def test_convex_positive_R(self):
        """Positive R → positive sag."""
        # R=100, a=25: sag = 100 - sqrt(100^2 - 25^2) = 100 - sqrt(9375)
        expected = 100.0 - math.sqrt(100.0**2 - 25.0**2)
        assert _sag(100.0, 25.0) == pytest.approx(expected)
        assert _sag(100.0, 25.0) > 0.0

    def test_concave_negative_R(self):
        """Negative R → negative sag."""
        sag = _sag(-100.0, 25.0)
        assert sag < 0.0
        assert sag == pytest.approx(-_sag(100.0, 25.0))

    def test_reject_when_R_too_small(self):
        """|R| < a + margin → ValueError."""
        with pytest.raises(ValueError, match="cannot exist"):
            _sag(20.0, 25.0)  # |R| = 20 < a = 25

    def test_reject_borderline(self):
        """|R| == a (exactly) → still rejected (within margin)."""
        with pytest.raises(ValueError, match="cannot exist"):
            _sag(25.0, 25.0)


# ===================================================================
# Lens power helper
# ===================================================================


class TestLensPower:
    def test_biconvex_positive(self):
        """R1 > 0, R2 < 0 → Φ > 0."""
        phi = _lens_power(100.0, -100.0, 1.5)
        assert phi > 0.0

    def test_biconcave_negative(self):
        """R1 < 0, R2 > 0 → Φ < 0."""
        phi = _lens_power(-100.0, 100.0, 1.5)
        assert phi < 0.0

    def test_plano_convex(self):
        """R1 > 0, R2 = ∞ → Φ > 0."""
        phi = _lens_power(100.0, float("inf"), 1.5)
        assert phi > 0.0

    def test_flat_is_zero(self):
        phi = _lens_power(float("inf"), float("inf"), 1.5)
        assert phi == pytest.approx(0.0)


# ===================================================================
# element_thickness – positive lens
# ===================================================================


class TestElementThicknessPositive:
    """Biconvex positive lens: centre thick > edge thin."""

    def test_basic_biconvex(self):
        """D=50, biconvex R1=100 R2=-100 n=1.5 → positive lens."""
        result = element_thickness(100.0, -100.0, 50.0, 1.5)
        assert result is not None
        t_center, t_edge = result
        assert t_center > 0
        assert t_edge > 0
        # Edge thickness should satisfy min for D=50 (1.8 mm)
        assert t_edge >= lookup_t_edge_min(50.0) - 1e-9

    def test_center_increases_with_aperture(self):
        """Larger aperture → larger centre thickness for same radii."""
        r1 = element_thickness(100.0, -100.0, 30.0, 1.5)
        r2 = element_thickness(100.0, -100.0, 50.0, 1.5)
        assert r1 is not None and r2 is not None
        assert r2[0] > r1[0]  # t_center grows with D

    def test_edge_equals_threshold(self):
        """For positive lens, t_edge should equal t_edge_min (exact match
        when the edge constraint dominates)."""
        # Use a strongly curved biconvex so edge constraint dominates
        result = element_thickness(60.0, -60.0, 50.0, 1.5)
        assert result is not None
        t_center, t_edge = result
        te_min = lookup_t_edge_min(50.0)
        # t_edge must be at least te_min
        assert t_edge >= te_min - 1e-9
        # If edge constraint dominated (no center-min override),
        # t_edge should be very close to te_min
        # (This only holds when t_from_edge > _MIN_CENTER_ABSOLUTE)
        if t_center > _MIN_CENTER_ABSOLUTE + 0.01:
            assert t_edge == pytest.approx(te_min, abs=1e-9)


# ===================================================================
# element_thickness – negative lens
# ===================================================================


class TestElementThicknessNegative:
    """Negative meniscus: centre thin, edge thick."""

    def test_basic_biconcave(self):
        """D=50, biconcave R1=-100 R2=100 n=1.5 → negative lens."""
        result = element_thickness(-100.0, 100.0, 50.0, 1.5)
        assert result is not None
        t_center, t_edge = result
        assert t_center > 0
        assert t_edge > 0
        # Centre thickness should be floored by lookup
        tc_min = lookup_t_center_min(50.0)
        assert t_center >= tc_min - 1e-9

    def test_center_floored_by_lookup(self):
        """For negative lens the centre-thickness table value dominates."""
        result = element_thickness(-200.0, 200.0, 50.0, 1.5)
        assert result is not None
        t_center, t_edge = result
        tc_min = lookup_t_center_min(50.0)  # 2.2 mm
        assert t_center >= tc_min - 1e-9
        # Edge is naturally larger for negative lens
        assert t_edge >= t_center


# ===================================================================
# element_thickness – rejection cases
# ===================================================================


class TestElementThicknessReject:
    def test_R_too_small_rejects(self):
        """|R| < D/2 → None."""
        # R=20 with D=50 → semi-aperture=25 > 20
        assert element_thickness(20.0, -100.0, 50.0, 1.5) is None

    def test_R_borderline_rejects(self):
        """|R| == a → still rejected."""
        assert element_thickness(25.0, -100.0, 50.0, 1.5) is None

    def test_both_surfaces_invalid(self):
        assert element_thickness(10.0, -10.0, 50.0, 1.5) is None


# ===================================================================
# Radius correction  –  keep Φ_thick = Φ_thin (Strategy A)
# ===================================================================


class TestCorrectRadiiBasic:
    """Basic properties of correct_radii_for_thickness."""

    def test_biconvex_returns_valid(self):
        """Biconvex lens with reasonable thickness → correction succeeds."""
        result = correct_radii_for_thickness(100.0, -100.0, 5.0, 1.5)
        assert result is not None
        R1p, R2p = result
        assert math.isfinite(R1p)
        assert math.isfinite(R2p)

    def test_flat_surface_no_change(self):
        """One flat surface → thickness term vanishes → no correction."""
        R1, R2 = 100.0, float("inf")
        result = correct_radii_for_thickness(R1, R2, 5.0, 1.5)
        assert result is not None
        assert result[0] == pytest.approx(R1)
        assert result[1] == R2  # inf stays inf

    def test_zero_thickness_no_change(self):
        """t=0 → B_prod=0 → no correction needed."""
        result = correct_radii_for_thickness(100.0, -100.0, 0.0, 1.5)
        assert result is not None
        assert result[0] == pytest.approx(100.0)
        assert result[1] == pytest.approx(-100.0)


class TestCorrectRadiiPowerPreservation:
    """Verify that Φ_thick(corrected) ≈ Φ_thin(original)."""

    @staticmethod
    def _thick_power(R1, R2, t, n):
        """Lensmaker's equation for thick lens."""

        def c(R):
            if not math.isfinite(R) or abs(R) < 1e-12:
                return 0.0
            return 1.0 / R

        c1, c2 = c(R1), c(R2)
        return (n - 1) * (c1 - c2 + (n - 1) * t * c1 * c2 / n)

    @staticmethod
    def _thin_power(R1, R2, n):
        def c(R):
            if not math.isfinite(R) or abs(R) < 1e-12:
                return 0.0
            return 1.0 / R

        return (n - 1) * (c(R1) - c(R2))

    @pytest.mark.parametrize(
        "R1,R2,t,n",
        [
            (100.0, -100.0, 5.0, 1.5),  # symmetric biconvex
            (60.0, -120.0, 8.0, 1.7),  # asymmetric positive
            (-80.0, -200.0, 4.0, 1.6),  # positive meniscus
            (200.0, 80.0, 6.0, 1.5),  # negative meniscus
            (-100.0, 100.0, 3.0, 1.5),  # symmetric biconcave
        ],
    )
    def test_power_preserved(self, R1, R2, t, n):
        """After correction, thick-lens power should match thin-lens power."""
        phi_thin = self._thin_power(R1, R2, n)
        result = correct_radii_for_thickness(R1, R2, t, n)
        assert result is not None
        R1p, R2p = result
        phi_thick = self._thick_power(R1p, R2p, t, n)
        assert phi_thick == pytest.approx(phi_thin, rel=1e-8)

    @pytest.mark.parametrize(
        "R1,R2,t,n",
        [
            (100.0, -100.0, 5.0, 1.5),
            (60.0, -120.0, 8.0, 1.7),
        ],
    )
    def test_shape_factor_preserved(self, R1, R2, t, n):
        """Shape factor B = (c1+c2)/(c1-c2) should be invariant."""

        def sf(Ra, Rb):
            ca = 1.0 / Ra if abs(Ra) > 1e-12 else 0.0
            cb = 1.0 / Rb if abs(Rb) > 1e-12 else 0.0
            denom = ca - cb
            if abs(denom) < 1e-15:
                return 0.0
            return (ca + cb) / denom

        B_orig = sf(R1, R2)
        result = correct_radii_for_thickness(R1, R2, t, n)
        assert result is not None
        B_corr = sf(result[0], result[1])
        assert B_corr == pytest.approx(B_orig, rel=1e-8)


class TestCorrectRadiiScale:
    """Verify the scale factor is close to 1 for typical lenses."""

    def test_scale_near_unity(self):
        """For D=50, t~5 mm, the correction should be small."""
        R1, R2, t, n = 100.0, -100.0, 5.0, 1.5
        result = correct_radii_for_thickness(R1, R2, t, n)
        assert result is not None
        # scale s = R1 / R1' (since R1' = R1/s)
        s = R1 / result[0]
        assert 0.95 < s < 1.05

    def test_larger_thickness_larger_correction(self):
        """Thicker lens → larger correction (s further from 1)."""
        R1, R2, n = 100.0, -100.0, 1.5
        res_thin = correct_radii_for_thickness(R1, R2, 2.0, n)
        res_thick = correct_radii_for_thickness(R1, R2, 15.0, n)
        assert res_thin is not None and res_thick is not None
        s_thin = R1 / res_thin[0]
        s_thick = R1 / res_thick[0]
        # thicker lens needs more correction
        assert abs(s_thick - 1.0) > abs(s_thin - 1.0)


class TestCorrectRadiiEdgeCases:
    """Edge cases and rejection."""

    def test_both_flat_no_change(self):
        result = correct_radii_for_thickness(float("inf"), float("inf"), 5.0, 1.5)
        assert result is not None
        assert result == (float("inf"), float("inf"))

    def test_zero_power_no_change(self):
        """c1 == c2 (same curvature) → zero power → s=1."""
        result = correct_radii_for_thickness(100.0, 100.0, 5.0, 1.5)
        assert result is not None
        assert result[0] == pytest.approx(100.0)
        assert result[1] == pytest.approx(100.0)


class TestReconcileCementedRadius:
    def test_average(self):
        assert _reconcile_cemented_radius(50.0, 52.0) == pytest.approx(51.0)

    def test_identical(self):
        assert _reconcile_cemented_radius(100.0, 100.0) == pytest.approx(100.0)

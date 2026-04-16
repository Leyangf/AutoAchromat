"""Unit tests for seidel_refine.py — thickness-aware Seidel refinement."""

from __future__ import annotations

import pytest

from autoachromat.models import Inputs, Candidate, ThickPrescription, ElementRx
from autoachromat.seidel_refine import (
    trace_marginal_cemented,
    trace_marginal_spaced,
    _P_from_Q_cemented,
    refine_seidel,
)
from autoachromat.cemented import _solve_Q_roots, radii_from_Q
from autoachromat.thickening import thicken_cemented


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def inp_cemented():
    return Inputs.with_defaults(
        fprime=200.0, D=50.0, system_type="cemented", P0=0.0, W0=0.0,
    )


@pytest.fixture
def inp_spaced():
    return Inputs.with_defaults(
        fprime=200.0, D=50.0, system_type="spaced", P0=0.0, W0=0.0, air_gap=10.0,
    )


def _make_cemented_candidate(inp: Inputs, n1=1.5168, n2=1.62004, nu1=64.17, nu2=36.37):
    """Build a realistic cemented candidate from Seidel synthesis."""
    from autoachromat.optics import achromat_power
    phi1, phi2 = achromat_power(nu1, nu2, inp.C0)

    Q_roots = _solve_Q_roots(n1, n2, phi1, phi2, inp.P0)
    assert len(Q_roots) > 0, "No Q roots for test glass pair"
    Q = Q_roots[0]

    R1, R2, R3 = radii_from_Q(inp, n1, n2, phi1, Q)

    return Candidate(
        system_type="cemented",
        glass1="N-BK7", glass2="N-F2",
        catalog1="SCHOTT", catalog2="SCHOTT",
        n1=n1, n2=n2, nu1=nu1, nu2=nu2,
        phi1=phi1, phi2=phi2,
        Q=Q,
        radii=[R1, R2, R3],
    )


# ---------------------------------------------------------------------------
# Trace tests
# ---------------------------------------------------------------------------


class TestTraceMarginal:
    """Verify paraxial ray trace properties."""

    def test_thin_limit_h_unity(self):
        """With near-zero thickness, all ray heights should be ~1."""
        trace = trace_marginal_cemented(
            R1=100.0, R2=-200.0, R3=-300.0,
            t1=1e-6, t2=1e-6,
            n1=1.5, n2=1.6,
        )
        for k in range(3):
            assert abs(trace.h[k] - 1.0) < 1e-4

    def test_h_deviates_with_thickness(self):
        """With real thickness, h2 and h3 should differ from 1."""
        trace = trace_marginal_cemented(
            R1=100.0, R2=-200.0, R3=-300.0,
            t1=8.0, t2=4.0,
            n1=1.5, n2=1.6,
        )
        assert trace.h[0] == 1.0
        assert trace.h[1] != 1.0
        assert trace.h[2] != 1.0

    def test_spaced_4_surfaces(self):
        """Spaced trace returns 4 surface states."""
        trace = trace_marginal_spaced(
            R1=100.0, R2=-200.0, R3=200.0, R4=-100.0,
            t1=6.0, t2=4.0, air_gap=5.0,
            n1=1.5, n2=1.6,
        )
        assert len(trace.h) == 4
        assert len(trace.u_before) == 4
        assert len(trace.u_after) == 4
        assert trace.h[0] == 1.0


# ---------------------------------------------------------------------------
# Seidel P tests
# ---------------------------------------------------------------------------


class TestThickSeidelP:
    """Verify thick-lens Seidel SA computation."""

    def test_thin_limit_matches_quadratic(self, inp_cemented):
        """With t→0, thick P should match the thin-lens quadratic result."""
        from autoachromat.optics import achromat_power
        n1, n2 = 1.5168, 1.62004
        nu1, nu2 = 64.17, 36.37
        phi1, phi2 = achromat_power(nu1, nu2, inp_cemented.C0)

        Q_roots = _solve_Q_roots(n1, n2, phi1, phi2, inp_cemented.P0)
        Q = Q_roots[0]

        # Thin-lens P should be ~P0 by construction
        P_thick = _P_from_Q_cemented(
            Q, phi1, n1, n2,
            t1=1e-6, t2=1e-6,
            f=inp_cemented.fprime,
        )
        assert P_thick is not None
        assert abs(P_thick - inp_cemented.P0) < 1e-4

    def test_thick_P_differs_from_thin(self):
        """With real thickness on a fast system, thick P deviates from P0."""
        # Use f/2 (D=50, f=100) to magnify the thick-lens effect
        inp_fast = Inputs.with_defaults(
            fprime=100.0, D=50.0, system_type="cemented", P0=0.0, W0=0.0,
        )
        cand = _make_cemented_candidate(inp_fast)
        rx = thicken_cemented(cand, inp_fast)
        assert rx is not None

        e1, e2 = rx.elements[0], rx.elements[1]
        assert cand.Q is not None
        P_thick = _P_from_Q_cemented(
            cand.Q, cand.phi1, cand.n1, cand.n2,
            t1=e1.t_center, t2=e2.t_center,
            f=inp_fast.fprime,
        )
        assert P_thick is not None
        # With real thickness on a fast system, deviation should be measurable
        assert abs(P_thick - inp_fast.P0) > 1e-8


# ---------------------------------------------------------------------------
# Newton solver tests
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Integration: refine_seidel
# ---------------------------------------------------------------------------


class TestRefineSeidel:
    """End-to-end refinement tests."""

    def test_cemented_refinement_changes_radii(self):
        """Refinement should produce a different prescription on a fast system."""
        # Use f/1.6 (D=50, f=80) to get a meaningful thick-lens correction
        inp_fast = Inputs.with_defaults(
            fprime=80.0, D=50.0, system_type="cemented", P0=0.0, W0=0.0,
        )
        cand = _make_cemented_candidate(inp_fast)
        rx_orig = thicken_cemented(cand, inp_fast)
        assert rx_orig is not None

        rx_refined = refine_seidel(rx_orig, cand, inp_fast)
        assert rx_refined is not None
        # The refined prescription should be a new object (not the same rx)
        assert rx_refined is not rx_orig

    def test_thin_limit_no_change(self, inp_cemented):
        """With very thin elements, refinement should barely change anything."""
        from autoachromat.optics import achromat_power
        n1, n2 = 1.5168, 1.62004
        nu1, nu2 = 64.17, 36.37
        phi1, phi2 = achromat_power(nu1, nu2, inp_cemented.C0)

        Q_roots = _solve_Q_roots(n1, n2, phi1, phi2, 0.0)
        Q = Q_roots[0]
        R1, R2, R3 = radii_from_Q(inp_cemented, n1, n2, phi1, Q)

        # Build a minimal-thickness prescription manually
        rx_thin = ThickPrescription(
            system_type="cemented",
            elements=[
                ElementRx(R_front=R1, R_back=R2, t_center=0.01, t_edge=0.01,
                          nd=n1, vd=nu1),
                ElementRx(R_front=R2, R_back=R3, t_center=0.01, t_edge=0.01,
                          nd=n2, vd=nu2),
            ],
            air_gap=None,
            back_focus_guess=200.0,
            D=50.0,
            wavelengths=(0.48613, 0.58756, 0.65627),
        )
        cand = Candidate(
            system_type="cemented",
            glass1="N-BK7", glass2="N-F2",
            catalog1="SCHOTT", catalog2="SCHOTT",
            n1=n1, n2=n2, nu1=nu1, nu2=nu2,
            phi1=phi1, phi2=phi2,
            Q=Q, radii=[R1, R2, R3],
        )

        rx_result = refine_seidel(rx_thin, cand, inp_cemented)
        # Should return original (Q change < threshold)
        assert rx_result is rx_thin

    def test_returns_original_on_none_Q(self, inp_cemented):
        """Candidate with Q=None should return rx unchanged."""
        cand = Candidate(
            system_type="cemented",
            glass1="A", glass2="B",
            catalog1="C", catalog2="C",
            n1=1.5, n2=1.6, nu1=60.0, nu2=30.0,
            phi1=0.5, phi2=0.5,
            Q=None,
            radii=[100.0, -200.0, -300.0],
        )
        rx = ThickPrescription(
            system_type="cemented",
            elements=[
                ElementRx(R_front=100.0, R_back=-200.0, t_center=5.0, t_edge=2.0,
                          nd=1.5, vd=60.0),
                ElementRx(R_front=-200.0, R_back=-300.0, t_center=3.0, t_edge=2.0,
                          nd=1.6, vd=30.0),
            ],
            air_gap=None,
            back_focus_guess=200.0,
            D=50.0,
            wavelengths=(0.48613, 0.58756, 0.65627),
        )
        assert refine_seidel(rx, cand, inp_cemented) is rx

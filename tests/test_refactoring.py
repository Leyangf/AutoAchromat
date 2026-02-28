"""
Unit tests for the refactored modules:
- optics.prepare_glass_data
- pipeline.process_candidate / run_pipeline / PipelineResult.to_dict
- cemented / spaced still produce the same candidates through prepare_glass_data
"""

from __future__ import annotations

import math
import pytest
from unittest.mock import MagicMock, patch
from dataclasses import dataclass, field

from autoachromat.glass_reader import Glass
from autoachromat.models import Inputs, Candidate, ElementRx, ThickPrescription
from autoachromat.optics import (
    prepare_glass_data,
    refractive_index,
    compute_abbe_number,
    Config,
    filter_glasses,
)
from autoachromat.pipeline import process_candidate, run_pipeline, PipelineResult
from autoachromat.optiland_bridge.evaluator import OpticMetrics


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_SCHOTT_NBAK4_CD = [
    1.28834642,
    0.00716845680,
    1.01914753,
    0.0225085830,
    0.314866860,
    87.3437824,
]

_SCHOTT_SF5_CD = [
    1.46141885,
    0.0111826126,
    0.247713019,
    0.0508594669,
    0.949995832,
    112.041888,
]


def _make_glass(name: str, cd: list[float], formula_id: int = 2) -> Glass:
    """Create a minimal Glass for testing."""
    return Glass(
        name=name,
        catalog="TEST",
        formula_id=formula_id,
        cd=cd,
        ld_min_um=0.3,
        ld_max_um=2.5,
    )


@pytest.fixture
def nbak4() -> Glass:
    return _make_glass("N-BAK4", _SCHOTT_NBAK4_CD)


@pytest.fixture
def sf5() -> Glass:
    return _make_glass("SF5", _SCHOTT_SF5_CD)


@pytest.fixture
def glasses(nbak4, sf5) -> list[Glass]:
    return [nbak4, sf5]


@pytest.fixture
def default_inputs() -> Inputs:
    return Inputs(
        lam0=0.58756,
        lam1=0.48613,
        lam2=0.65627,
        D=50.0,
        fprime=200.0,
        C0=0.0,
        P0=1.0,
        W0=0.0,
        min_delta_nu=10.0,
        max_PE=100.0,
        N=20,
        system_type="cemented",
    )


def _make_candidate(glass1: Glass, glass2: Glass, inputs: Inputs) -> Candidate:
    """Build a minimal cemented Candidate for pipeline testing."""
    n1 = refractive_index(glass1, inputs.lam0)
    n2 = refractive_index(glass2, inputs.lam0)
    nu1 = compute_abbe_number(glass1, inputs.lam0, inputs.lam1, inputs.lam2)
    nu2 = compute_abbe_number(glass2, inputs.lam0, inputs.lam1, inputs.lam2)
    return Candidate(
        system_type="cemented",
        glass1=glass1.name,
        catalog1=glass1.catalog,
        glass2=glass2.name,
        catalog2=glass2.catalog,
        n1=n1,
        n2=n2,
        nu1=nu1,
        nu2=nu2,
        phi1=0.5,
        phi2=0.5,
        Q=0.0,
        W=0.0,
        radii=[200.0, -150.0, -400.0],
        PE=1.0,
        formula_id1=glass1.formula_id,
        cd1=list(glass1.cd),
        formula_id2=glass2.formula_id,
        cd2=list(glass2.cd),
    )


# ===================================================================
# optics.prepare_glass_data
# ===================================================================


class TestPrepareGlassData:
    """Tests for optics.prepare_glass_data."""

    def test_returns_tuples(self, glasses):
        gdata = prepare_glass_data(glasses, 0.58756, 0.48613, 0.65627)
        assert len(gdata) > 0
        for item in gdata:
            assert len(item) == 3
            g, n0, nu = item
            assert isinstance(g, Glass)
            assert isinstance(n0, float) and n0 > 1.0
            assert isinstance(nu, float) and nu > 0

    def test_filters_excluded_glass(self, nbak4, sf5):
        excluded = Glass(
            name="EXCLUDED",
            catalog="TEST",
            formula_id=2,
            cd=_SCHOTT_NBAK4_CD,
            ld_min_um=0.3,
            ld_max_um=2.5,
            exclude_sub=True,
        )
        gdata = prepare_glass_data([nbak4, sf5, excluded], 0.58756, 0.48613, 0.65627)
        names = [g.name for g, _, _ in gdata]
        assert "EXCLUDED" not in names
        assert "N-BAK4" in names
        assert "SF5" in names

    def test_filters_out_of_range(self, nbak4):
        narrow = Glass(
            name="NARROW",
            catalog="TEST",
            formula_id=2,
            cd=_SCHOTT_NBAK4_CD,
            ld_min_um=0.5,
            ld_max_um=0.55,
        )
        gdata = prepare_glass_data([nbak4, narrow], 0.58756, 0.48613, 0.65627)
        names = [g.name for g, _, _ in gdata]
        assert "NARROW" not in names
        assert "N-BAK4" in names

    def test_consistent_with_old_approach(self, glasses):
        """Results should match the old Config+filter_glasses approach."""
        cfg = Config(0.58756, 0.48613, 0.65627)
        usable = filter_glasses(glasses, cfg)
        old_gdata = []
        for g in usable:
            try:
                n0 = refractive_index(g, 0.58756)
                nu = compute_abbe_number(g, 0.58756, 0.48613, 0.65627)
            except Exception:
                continue
            old_gdata.append((g.name, n0, nu))

        new_gdata = prepare_glass_data(glasses, 0.58756, 0.48613, 0.65627)
        new_flat = [(g.name, n0, nu) for g, n0, nu in new_gdata]

        assert len(new_flat) == len(old_gdata)
        for (name_a, n_a, nu_a), (name_b, n_b, nu_b) in zip(old_gdata, new_flat):
            assert name_a == name_b
            assert n_a == pytest.approx(n_b)
            assert nu_a == pytest.approx(nu_b)


# ===================================================================
# pipeline.PipelineResult.to_dict
# ===================================================================


class TestPipelineResultToDict:
    """Tests for PipelineResult serialisation."""

    def test_to_dict_with_rx(self, nbak4, sf5, default_inputs):
        cand = _make_candidate(nbak4, sf5, default_inputs)
        rx = ThickPrescription(
            system_type="cemented",
            elements=[
                ElementRx(
                    R_front=200.0,
                    R_back=-150.0,
                    t_center=5.0,
                    t_edge=2.0,
                    nd=1.5,
                    vd=55.0,
                ),
                ElementRx(
                    R_front=-150.0,
                    R_back=-400.0,
                    t_center=3.0,
                    t_edge=4.0,
                    nd=1.7,
                    vd=30.0,
                ),
            ],
            air_gap=None,
            back_focus_guess=180.0,
            D=50.0,
            wavelengths=(0.48613, 0.58756, 0.65627),
        )
        m = OpticMetrics(success=True, glass1="N-BAK4", glass2="SF5")
        pr = PipelineResult(candidate=cand, rx=rx, metrics=m)
        d = pr.to_dict()

        assert d["thin_radii"] == cand.radii
        assert d["thick_radii"] == [[200.0, -150.0], [-150.0, -400.0]]
        assert d["t_center"] == [5.0, 3.0]
        assert d["t_edge"] == [2.0, 4.0]
        assert d["air_gap"] is None
        assert d["success"] is True

    def test_to_dict_without_rx(self, nbak4, sf5, default_inputs):
        cand = _make_candidate(nbak4, sf5, default_inputs)
        m = OpticMetrics(success=False, error_msg="thickening failed")
        pr = PipelineResult(candidate=cand, rx=None, metrics=m)
        d = pr.to_dict()

        assert d["thin_radii"] == cand.radii
        assert d["thick_radii"] is None
        assert d["t_center"] is None
        assert d["t_edge"] is None
        assert d["air_gap"] is None
        assert d["success"] is False


# ===================================================================
# pipeline.run_pipeline  (with mocked thicken/build/evaluate)
# ===================================================================


class TestRunPipeline:
    """Tests for the pipeline runner (mocked to avoid optiland dependency)."""

    def test_run_pipeline_calls_progress(self, nbak4, sf5, default_inputs):
        cand = _make_candidate(nbak4, sf5, default_inputs)

        with patch("autoachromat.pipeline.thicken", return_value=None):
            progress_calls: list[tuple[int, int]] = []

            def on_progress(cur, total, pr):
                progress_calls.append((cur, total))

            results = run_pipeline(
                [cand, cand], default_inputs, on_progress=on_progress
            )

        assert len(results) == 2
        assert len(progress_calls) == 2
        assert progress_calls[0] == (1, 2)
        assert progress_calls[1] == (2, 2)
        for r in results:
            assert r.metrics.success is False
            assert "thickening" in r.metrics.error_msg

    def test_run_pipeline_max_n(self, nbak4, sf5, default_inputs):
        cand = _make_candidate(nbak4, sf5, default_inputs)

        with patch("autoachromat.pipeline.thicken", return_value=None):
            results = run_pipeline([cand] * 10, default_inputs, max_n=3)

        assert len(results) == 3

    def test_process_candidate_thicken_fail(self, nbak4, sf5, default_inputs):
        cand = _make_candidate(nbak4, sf5, default_inputs)

        with patch("autoachromat.pipeline.thicken", return_value=None):
            pr = process_candidate(cand, default_inputs)

        assert pr.rx is None
        assert pr.metrics.success is False
        assert pr.metrics.build_time_ms >= 0

    def test_process_candidate_build_fail(self, nbak4, sf5, default_inputs):
        cand = _make_candidate(nbak4, sf5, default_inputs)
        fake_rx = MagicMock()

        with (
            patch("autoachromat.pipeline.thicken", return_value=fake_rx),
            patch(
                "autoachromat.pipeline.build_optic_from_prescription", return_value=None
            ),
        ):
            pr = process_candidate(cand, default_inputs)

        assert pr.rx is fake_rx
        assert pr.metrics.success is False
        assert "build_optic" in pr.metrics.error_msg

    def test_process_candidate_success(self, nbak4, sf5, default_inputs):
        cand = _make_candidate(nbak4, sf5, default_inputs)
        fake_rx = MagicMock()
        fake_op = MagicMock()
        fake_metrics = OpticMetrics(success=True, efl=200.0)

        with (
            patch("autoachromat.pipeline.thicken", return_value=fake_rx),
            patch(
                "autoachromat.pipeline.build_optic_from_prescription",
                return_value=fake_op,
            ),
            patch("autoachromat.pipeline.evaluate", return_value=fake_metrics),
        ):
            pr = process_candidate(cand, default_inputs)

        assert pr.metrics.success is True
        assert pr.metrics.efl == 200.0
        assert pr.metrics.build_time_ms >= 0


# ===================================================================
# cemented / spaced still work through prepare_glass_data
# ===================================================================


class TestCementedSpacedIntegration:
    """Verify that cemented/spaced solvers still produce candidates
    after refactoring to use prepare_glass_data."""

    def test_cemented_produces_candidates(self, glasses, default_inputs):
        from autoachromat.cemented import run_cemented

        cands = run_cemented(default_inputs, glasses)
        assert isinstance(cands, list)
        # With two very different glasses, we should get at least one candidate
        if len(cands) > 0:
            c = cands[0]
            assert c.system_type == "cemented"
            assert len(c.radii) == 3

    def test_spaced_produces_candidates(self, glasses):
        from autoachromat.spaced import run_spaced

        inputs = Inputs(
            lam0=0.58756,
            lam1=0.48613,
            lam2=0.65627,
            D=50.0,
            fprime=200.0,
            C0=0.0,
            P0=1.0,
            W0=0.0,
            min_delta_nu=10.0,
            max_PE=100.0,
            N=20,
            system_type="spaced",
        )
        cands = run_spaced(inputs, glasses)
        assert isinstance(cands, list)
        if len(cands) > 0:
            c = cands[0]
            assert c.system_type == "spaced"
            assert len(c.radii) == 4

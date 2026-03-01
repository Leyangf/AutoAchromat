"""
Contract tests – verify serialization key stability, failure semantics
and public API invariants that downstream consumers may rely on.
"""

from __future__ import annotations

import pytest
from unittest.mock import patch, MagicMock

from autoachromat.models import Inputs, Candidate, ElementRx, ThickPrescription
from autoachromat.pipeline import (
    PipelineResult,
    run_pipeline,
    run_design,
    DesignResult,
)
from autoachromat.optiland_bridge.evaluator import OpticMetrics


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def default_inputs() -> Inputs:
    return Inputs.with_defaults()


@pytest.fixture
def sample_candidate() -> Candidate:
    return Candidate(
        system_type="cemented",
        glass1="N-BAK4",
        catalog1="SCHOTT",
        glass2="SF5",
        catalog2="SCHOTT",
        n1=1.5688,
        n2=1.6727,
        nu1=55.98,
        nu2=32.21,
        phi1=0.5,
        phi2=0.5,
        Q=0.0,
        W=0.0,
        radii=[200.0, -150.0, -400.0],
        PE=1.0,
    )


@pytest.fixture
def sample_rx() -> ThickPrescription:
    return ThickPrescription(
        system_type="cemented",
        elements=[
            ElementRx(
                R_front=200.0,
                R_back=-150.0,
                t_center=5.0,
                t_edge=2.0,
                nd=1.5688,
                vd=55.98,
            ),
            ElementRx(
                R_front=-150.0,
                R_back=-400.0,
                t_center=3.0,
                t_edge=4.0,
                nd=1.6727,
                vd=32.21,
            ),
        ],
        air_gap=None,
        back_focus_guess=180.0,
        D=50.0,
        wavelengths=(0.48613, 0.58756, 0.65627),
    )


# ===================================================================
# Serialization key stability  (to_dict output contract)
# ===================================================================


class TestToDictKeyStability:
    """to_dict() must always produce certain keys that consumers rely on."""

    REQUIRED_KEYS = {
        "success",
        "error_msg",
        "thin_radii",
        "thick_radii",
        "t_center",
        "t_edge",
        "air_gap",
        "efl",
        "fno",
        "bfd",
        "rms_spot_radius",
        "geo_spot_radius",
        "SA",
        "CC",
        "LchC",
        "TchC",
        "PE",
        "glass1",
        "glass2",
        "build_time_ms",
        "eval_time_ms",
        # Thermal analysis keys
        "thermal_available",
        "dn_dT_1",
        "dn_dT_2",
        "V1_ppm_K",
        "V2_ppm_K",
        "dphi_dT_norm",
        "alpha_housing_required_ppm_K",
    }

    def test_to_dict_has_required_keys_with_rx(self, sample_candidate, sample_rx):
        m = OpticMetrics(success=True)
        pr = PipelineResult(candidate=sample_candidate, rx=sample_rx, metrics=m)
        d = pr.to_dict()
        missing = self.REQUIRED_KEYS - set(d.keys())
        assert not missing, f"Missing required keys: {missing}"

    def test_to_dict_has_required_keys_without_rx(self, sample_candidate):
        m = OpticMetrics(success=False, error_msg="test failure")
        pr = PipelineResult(candidate=sample_candidate, rx=None, metrics=m)
        d = pr.to_dict()
        missing = self.REQUIRED_KEYS - set(d.keys())
        assert not missing, f"Missing required keys: {missing}"

    def test_to_dict_thick_radii_structure(self, sample_candidate, sample_rx):
        """thick_radii must be a list of [R_front, R_back] pairs."""
        m = OpticMetrics(success=True)
        pr = PipelineResult(candidate=sample_candidate, rx=sample_rx, metrics=m)
        d = pr.to_dict()
        assert isinstance(d["thick_radii"], list)
        for pair in d["thick_radii"]:
            assert isinstance(pair, list) and len(pair) == 2

    def test_to_dict_none_fields_without_rx(self, sample_candidate):
        """When rx is None, prescription fields must be None (not missing)."""
        m = OpticMetrics(success=False)
        pr = PipelineResult(candidate=sample_candidate, rx=None, metrics=m)
        d = pr.to_dict()
        assert d["thick_radii"] is None
        assert d["t_center"] is None
        assert d["t_edge"] is None
        assert d["air_gap"] is None


# ===================================================================
# Failure semantics
# ===================================================================


class TestFailureSemantics:
    """Pipeline must never raise; failures are captured in metrics."""

    def test_process_candidate_never_raises(self, sample_candidate, default_inputs):
        """process_candidate must return PipelineResult even on failure."""
        from autoachromat.pipeline import process_candidate

        with patch("autoachromat.pipeline.thicken", side_effect=RuntimeError("boom")):
            pr = process_candidate(sample_candidate, default_inputs)
        # Should not raise — failure is in metrics
        assert isinstance(pr, PipelineResult)

    def test_run_pipeline_returns_same_count(self, sample_candidate, default_inputs):
        """run_pipeline must return exactly len(candidates) results."""
        with patch("autoachromat.pipeline.thicken", return_value=None):
            results = run_pipeline(
                [sample_candidate, sample_candidate, sample_candidate],
                default_inputs,
            )
        assert len(results) == 3

    def test_failed_metrics_have_error_msg(self, sample_candidate, default_inputs):
        """On failure, error_msg must be non-empty."""
        with patch("autoachromat.pipeline.thicken", return_value=None):
            pr = run_pipeline([sample_candidate], default_inputs)[0]
        assert not pr.metrics.success
        assert pr.metrics.error_msg  # non-empty string


# ===================================================================
# Inputs.with_defaults contract
# ===================================================================


class TestInputsWithDefaults:
    """Inputs.with_defaults() must produce a valid frozen Inputs."""

    def test_returns_inputs(self):
        inp = Inputs.with_defaults()
        assert isinstance(inp, Inputs)

    def test_overrides_applied(self):
        inp = Inputs.with_defaults(fprime=100.0, system_type="spaced")
        assert inp.fprime == 100.0
        assert inp.system_type == "spaced"
        # Other defaults intact
        assert inp.lam0 == 0.58756
        assert inp.D == 50.0

    def test_frozen(self):
        inp = Inputs.with_defaults()
        with pytest.raises(AttributeError):
            inp.fprime = 999.0  # type: ignore[misc]

    def test_unknown_field_raises(self):
        with pytest.raises(TypeError):
            Inputs.with_defaults(nonexistent_field=42)


# ===================================================================
# DesignResult structure
# ===================================================================


class TestDesignResult:
    """run_design must return a well-formed DesignResult."""

    def test_run_design_returns_design_result(self, default_inputs):
        with (
            patch("autoachromat.pipeline.load_catalog", return_value=[]),
            patch("autoachromat.pipeline.run_cemented", return_value=[]),
        ):
            dr = run_design(default_inputs, ["dummy.agf"])
        assert isinstance(dr, DesignResult)
        assert isinstance(dr.candidates, list)
        assert isinstance(dr.results, list)
        assert isinstance(dr.n_glasses, int)
        assert isinstance(dr.synth_time_ms, float)

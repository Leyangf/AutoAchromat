"""
gui.py – Tkinter GUI for AutoAchromat doublet designer.

Pipeline:  config → thin-lens synthesis → thickening → ray tracing → results

Usage:
    python -m autoachromat.gui
"""

from __future__ import annotations

import json
import math
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Any, Optional

# ---------------------------------------------------------------------------
# Matplotlib – set non-interactive backend BEFORE any project import that
# pulls in optiland (which itself imports matplotlib.pyplot at module level).
# FigureCanvasTkAgg renders into a tkinter widget regardless of backend.
# ---------------------------------------------------------------------------
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# ---------------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------------
from .models import Candidate, Inputs, ThickPrescription
from .pipeline import run_design, run_stage_b, PipelineResult
from .optiland_bridge.evaluator import OpticMetrics
from .optiland_bridge.builder import build_optic_from_prescription


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _fmt(v: Optional[float], fmt: str = ".4f") -> str:
    """Format nullable float for display."""
    if v is None:
        return "—"
    try:
        if not math.isfinite(v):
            return "—"
    except (TypeError, ValueError):
        return "—"
    return f"{v:{fmt}}"


def _interp_transmittance(
    trans_data: list[tuple[float, float, float]], lam_um: float,
) -> float | None:
    """Linearly interpolate internal transmittance at a given wavelength.

    trans_data is [(wavelength_um, transmittance, thickness_mm), ...].
    Returns transmittance (0-1) or None if out of range.
    """
    if not trans_data:
        return None
    # Sort by wavelength
    pts = sorted(trans_data, key=lambda t: t[0])
    wls = [p[0] for p in pts]
    trs = [p[1] for p in pts]

    if lam_um < wls[0] or lam_um > wls[-1]:
        return None
    if len(wls) == 1:
        return trs[0] if abs(lam_um - wls[0]) < 1e-6 else None

    for i in range(len(wls) - 1):
        if wls[i] <= lam_um <= wls[i + 1]:
            t = (lam_um - wls[i]) / (wls[i + 1] - wls[i])
            return trs[i] + t * (trs[i + 1] - trs[i])
    return None


# ---------------------------------------------------------------------------
# Result row  (one per candidate)
# ---------------------------------------------------------------------------


class ResultRow:
    """Bundles a PipelineResult with an index for display."""

    __slots__ = ("result", "idx", "stage", "src_idx")

    def __init__(self, idx: int, result: PipelineResult, stage: str = "A", src_idx: int = -1):
        self.idx = idx
        self.result = result
        self.stage = stage
        self.src_idx = src_idx  # Stage A row index that this was optimised from

    @property
    def cand(self) -> Candidate:
        return self.result.candidate

    @property
    def rx(self) -> Optional[ThickPrescription]:
        return self.result.rx

    @property
    def metrics(self) -> OpticMetrics:
        return self.result.metrics


# ---------------------------------------------------------------------------
# Main application window
# ---------------------------------------------------------------------------

_PAD: dict[str, Any] = dict(padx=6, pady=3)


class AutoAchromatGUI(tk.Tk):
    """Main GUI window for AutoAchromat."""

    # (name, lam1_nm, lam0_nm, lam2_nm)
    _WL_PRESETS = [
        ("Visible d/F/C", 486.13, 587.56, 656.27),
        ("Visible e/F'/C'", 479.99, 546.07, 643.85),
        ("NIR (0.8-1.0)", 800, 900, 1000),
        ("SWIR (1.0-1.7)", 1000, 1350, 1700),
        ("MWIR (3-5)", 3000, 4000, 5000),
        ("LWIR (8-12)", 8000, 10000, 12000),
    ]

    def __init__(self) -> None:
        super().__init__()
        self.title("AutoAchromat – Doublet Designer")
        self.geometry("1280x820")
        self.minsize(960, 640)

        # State – auto-load SCHOTT catalog if present
        _default_schott = (
            Path(__file__).resolve().parent.parent.parent
            / "data"
            / "catalogs"
            / "SCHOTT.AGF"
        )
        self._catalog_paths: list[str] = (
            [str(_default_schott)] if _default_schott.is_file() else []
        )
        self._results: list[ResultRow] = []
        self._running = False
        self._sort_col: str = "rms"
        self._sort_reverse: bool = False

        self._build_ui()
        self._refresh_catalog_display()

    # ===================================================================
    # UI construction
    # ===================================================================

    def _build_ui(self) -> None:
        top = ttk.Frame(self)
        top.pack(fill=tk.X, padx=8, pady=4)

        self._build_global_frame(top)

        # Stage A (left) and Stage B (right) — equal width
        stage_row = ttk.Frame(top)
        stage_row.pack(fill=tk.X, **_PAD)
        stage_row.columnconfigure(0, weight=1)
        stage_row.columnconfigure(1, weight=1)
        self._build_stage_a_frame(stage_row)
        self._build_stage_b_frame(stage_row)

        self._build_status_bar(top)

        # PanedWindow: ③ Results | ④ Details | ⑤ Drawing  (equal weight)
        self._main_pane = ttk.PanedWindow(self, orient=tk.VERTICAL)
        self._main_pane.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        self._build_result_table(self._main_pane)
        self._build_detail_panel(self._main_pane)
        self._build_drawing_panel(self._main_pane)

        # After first render, force equal sash positions
        self.after(50, self._equalize_pane_heights)

    # ------------------------------------------------------------------
    # Global parameters (row 0 of top area)
    # ------------------------------------------------------------------

    def _build_global_frame(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text=" ① Global ")
        frame.pack(fill=tk.X, **_PAD)

        # Row 0: geometry + wavelengths — use grid for consistent spacing
        for c in (1, 3, 5, 7):
            frame.columnconfigure(c, weight=0)
        frame.columnconfigure(9, weight=1)

        r = 0
        ttk.Label(frame, text="f' [mm]:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        self._var_fprime = tk.DoubleVar(value=200.0)
        ttk.Entry(frame, textvariable=self._var_fprime, width=7).grid(
            row=r, column=1, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="D [mm]:").grid(row=r, column=2, sticky=tk.E, **_PAD)
        d_frame = ttk.Frame(frame)
        d_frame.grid(row=r, column=3, sticky=tk.W, **_PAD)
        self._var_D = tk.DoubleVar(value=50.0)
        ttk.Entry(d_frame, textvariable=self._var_D, width=7).pack(side=tk.LEFT)
        self._var_fno_display = tk.StringVar(value="F/4.0")
        ttk.Label(d_frame, textvariable=self._var_fno_display, width=5).pack(
            side=tk.LEFT, padx=(4, 0),
        )
        ttk.Label(frame, text="Field [°]:").grid(row=r, column=4, sticky=tk.E, **_PAD)
        self._var_half_field = tk.DoubleVar(value=1.0)
        ttk.Entry(frame, textvariable=self._var_half_field, width=5).grid(
            row=r, column=5, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="Air gap [mm]:").grid(row=r, column=6, sticky=tk.E, **_PAD)
        self._var_air_gap = tk.DoubleVar(value=1.0)
        self._ent_air_gap = ttk.Entry(frame, textvariable=self._var_air_gap, width=5)
        self._ent_air_gap.grid(row=r, column=7, sticky=tk.W, **_PAD)
        self._ent_air_gap.configure(state=tk.DISABLED)

        ttk.Label(frame, text="λ [nm]:").grid(row=r, column=8, sticky=tk.E, **_PAD)
        wl_frame = ttk.Frame(frame)
        wl_frame.grid(row=r, column=9, sticky=tk.W, **_PAD)
        self._var_wl_preset = tk.StringVar(value="Visible d/F/C")
        wl_combo = ttk.Combobox(
            wl_frame, textvariable=self._var_wl_preset, width=14,
            state="readonly",
            values=[name for name, *_ in self._WL_PRESETS],
        )
        wl_combo.pack(side=tk.LEFT, padx=(0, 4))
        wl_combo.bind("<<ComboboxSelected>>", lambda _: self._on_wl_preset())
        self._var_lam1_nm = tk.DoubleVar(value=486.13)
        self._var_lam0_nm = tk.DoubleVar(value=587.56)
        self._var_lam2_nm = tk.DoubleVar(value=656.27)
        ttk.Label(wl_frame, text="λ₁").pack(side=tk.LEFT)
        ttk.Entry(wl_frame, textvariable=self._var_lam1_nm, width=7).pack(side=tk.LEFT, padx=1)
        ttk.Label(wl_frame, text="λ₀").pack(side=tk.LEFT, padx=(3, 0))
        ttk.Entry(wl_frame, textvariable=self._var_lam0_nm, width=7).pack(side=tk.LEFT, padx=1)
        ttk.Label(wl_frame, text="λ₂").pack(side=tk.LEFT, padx=(3, 0))
        ttk.Entry(wl_frame, textvariable=self._var_lam2_nm, width=7).pack(side=tk.LEFT, padx=1)

        # Row 1: type, top-N, catalogs
        r = 1
        ttk.Label(frame, text="Type:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        type_frame = ttk.Frame(frame)
        type_frame.grid(row=r, column=1, columnspan=3, sticky=tk.W, **_PAD)
        self._var_type = tk.StringVar(value="cemented")
        ttk.Radiobutton(
            type_frame, text="Cemented", variable=self._var_type,
            value="cemented", command=self._on_type_changed,
        ).pack(side=tk.LEFT, padx=2)
        ttk.Radiobutton(
            type_frame, text="Spaced", variable=self._var_type,
            value="spaced", command=self._on_type_changed,
        ).pack(side=tk.LEFT, padx=2)

        ttk.Label(frame, text="Top-N:").grid(row=r, column=4, sticky=tk.E, **_PAD)
        self._var_N = tk.IntVar(value=20)
        ttk.Entry(frame, textvariable=self._var_N, width=5).grid(
            row=r, column=5, sticky=tk.W, **_PAD
        )

        ttk.Label(frame, text="Catalogs:").grid(row=r, column=6, sticky=tk.E, **_PAD)
        cat_frame = ttk.Frame(frame)
        cat_frame.grid(row=r, column=7, columnspan=3, sticky=tk.W, **_PAD)
        self._var_cat_display = tk.StringVar(value="(none)")
        ttk.Label(
            cat_frame, textvariable=self._var_cat_display, width=30, relief="sunken",
        ).pack(side=tk.LEFT)
        ttk.Button(cat_frame, text="+", width=3, command=self._add_catalog).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Button(cat_frame, text="−", width=3, command=self._clear_catalogs).pack(
            side=tk.LEFT, padx=2
        )

        self._var_fprime.trace_add("write", lambda *_: self._update_fno_display())
        self._var_D.trace_add("write", lambda *_: self._update_fno_display())
        self._var_D.trace_add("write", lambda *_: self._auto_fill_thickness())

    # ------------------------------------------------------------------
    # Stage A (left) and Stage B (right)
    # ------------------------------------------------------------------

    def _build_stage_a_frame(self, parent: ttk.Frame) -> None:
        outer = ttk.LabelFrame(parent, text=" ② Stage A: Thin-Lens Synthesis ")
        outer.grid(row=0, column=0, sticky="nsew", **_PAD)

        # Main row: manufacturing + run (fixed layout, never changes)
        main_row = ttk.Frame(outer)
        main_row.pack(fill=tk.X, **_PAD)

        ttk.Label(main_row, text="te min:").pack(side=tk.LEFT, **_PAD)
        self._var_te_min = tk.DoubleVar(value=0.0)
        ttk.Entry(main_row, textvariable=self._var_te_min, width=6).pack(
            side=tk.LEFT, **_PAD
        )
        ttk.Label(main_row, text="tc min:").pack(side=tk.LEFT, **_PAD)
        self._var_tc_min = tk.DoubleVar(value=0.0)
        ttk.Entry(main_row, textvariable=self._var_tc_min, width=6).pack(
            side=tk.LEFT, **_PAD
        )
        ttk.Button(
            main_row, text="Auto", width=5, command=self._auto_fill_thickness,
        ).pack(side=tk.LEFT, **_PAD)

        self._btn_run = ttk.Button(
            main_row, text=" ▶ Run Stage A ", command=self._on_run
        )
        self._btn_run.pack(side=tk.RIGHT, padx=6, pady=3)

        # Advanced toggle
        self._adv_visible = tk.BooleanVar(value=False)
        self._btn_adv = ttk.Checkbutton(
            outer, text="▶ Advanced", variable=self._adv_visible,
            command=self._toggle_advanced,
        )
        self._btn_adv.pack(anchor=tk.W, padx=6)

        # Advanced fields (hidden by default, packed below toggle)
        self._adv_frame = ttk.Frame(outer)

        self._var_C0 = tk.DoubleVar(value=0.0)
        self._var_P0 = tk.DoubleVar(value=0.0)
        self._var_W0 = tk.DoubleVar(value=0.0)
        self._var_min_dnu = tk.DoubleVar(value=10.0)
        self._var_max_PE = tk.DoubleVar(value=100.0)

        ttk.Label(self._adv_frame, text="SA target (0=corr):").pack(side=tk.LEFT, **_PAD)
        ttk.Entry(self._adv_frame, textvariable=self._var_P0, width=6).pack(side=tk.LEFT, **_PAD)
        ttk.Label(self._adv_frame, text="Coma:").pack(side=tk.LEFT, **_PAD)
        ttk.Entry(self._adv_frame, textvariable=self._var_W0, width=6).pack(side=tk.LEFT, **_PAD)
        ttk.Label(self._adv_frame, text="LCA:").pack(side=tk.LEFT, **_PAD)
        ttk.Entry(self._adv_frame, textvariable=self._var_C0, width=6).pack(side=tk.LEFT, **_PAD)
        ttk.Label(self._adv_frame, text="Δν min:").pack(side=tk.LEFT, **_PAD)
        ttk.Entry(self._adv_frame, textvariable=self._var_min_dnu, width=5).pack(side=tk.LEFT, **_PAD)
        ttk.Label(self._adv_frame, text="PE max:").pack(side=tk.LEFT, **_PAD)
        ttk.Entry(self._adv_frame, textvariable=self._var_max_PE, width=5).pack(side=tk.LEFT, **_PAD)

        # Auto-fill manufacturing thickness from D on startup
        self._auto_fill_thickness()

    def _build_stage_b_frame(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text=" ② Stage B: Optimisation (0 = auto) ")
        frame.grid(row=0, column=1, sticky="nsew", **_PAD)

        r = 0
        ttk.Label(frame, text="t₁:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        rng1 = ttk.Frame(frame)
        rng1.grid(row=r, column=1, sticky=tk.W, **_PAD)
        self._var_t1_min = tk.DoubleVar(value=0.0)
        ttk.Entry(rng1, textvariable=self._var_t1_min, width=5).pack(side=tk.LEFT)
        ttk.Label(rng1, text="–").pack(side=tk.LEFT, padx=2)
        self._var_t1_max = tk.DoubleVar(value=0.0)
        ttk.Entry(rng1, textvariable=self._var_t1_max, width=5).pack(side=tk.LEFT)

        ttk.Label(frame, text="t₂:").grid(row=r, column=2, sticky=tk.E, **_PAD)
        rng2 = ttk.Frame(frame)
        rng2.grid(row=r, column=3, sticky=tk.W, **_PAD)
        self._var_t2_min = tk.DoubleVar(value=0.0)
        ttk.Entry(rng2, textvariable=self._var_t2_min, width=5).pack(side=tk.LEFT)
        ttk.Label(rng2, text="–").pack(side=tk.LEFT, padx=2)
        self._var_t2_max = tk.DoubleVar(value=0.0)
        ttk.Entry(rng2, textvariable=self._var_t2_max, width=5).pack(side=tk.LEFT)

        ttk.Label(frame, text="Air gap:").grid(row=r, column=4, sticky=tk.E, **_PAD)
        rng_gap = ttk.Frame(frame)
        rng_gap.grid(row=r, column=5, sticky=tk.W, **_PAD)
        self._var_gap_min = tk.DoubleVar(value=0.0)
        self._ent_gap_min = ttk.Entry(rng_gap, textvariable=self._var_gap_min, width=5)
        self._ent_gap_min.pack(side=tk.LEFT)
        ttk.Label(rng_gap, text="–").pack(side=tk.LEFT, padx=2)
        self._var_gap_max = tk.DoubleVar(value=0.0)
        self._ent_gap_max = ttk.Entry(rng_gap, textvariable=self._var_gap_max, width=5)
        self._ent_gap_max.pack(side=tk.LEFT)
        self._ent_gap_min.configure(state=tk.DISABLED)
        self._ent_gap_max.configure(state=tk.DISABLED)

        r = 1
        ttk.Label(frame, text="W_efl:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        self._var_w_efl = tk.DoubleVar(value=10.0)
        ttk.Entry(frame, textvariable=self._var_w_efl, width=5).grid(
            row=r, column=1, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="W_rms:").grid(row=r, column=2, sticky=tk.E, **_PAD)
        self._var_w_rms = tk.DoubleVar(value=4.0)
        ttk.Entry(frame, textvariable=self._var_w_rms, width=5).grid(
            row=r, column=3, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="W_field:").grid(row=r, column=4, sticky=tk.E, **_PAD)
        self._var_w_field = tk.DoubleVar(value=1.0)
        ttk.Entry(frame, textvariable=self._var_w_field, width=5).grid(
            row=r, column=5, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="Iter:").grid(row=r, column=6, sticky=tk.E, **_PAD)
        self._var_maxiter = tk.IntVar(value=200)
        ttk.Entry(frame, textvariable=self._var_maxiter, width=5).grid(
            row=r, column=7, sticky=tk.W, **_PAD
        )
        self._btn_stage_b = ttk.Button(
            frame, text=" ⚡ Optimize ",
            command=self._on_stage_b, state=tk.DISABLED,
        )
        self._btn_stage_b.grid(row=r, column=8, sticky=tk.E, padx=6, pady=3)

    # ------------------------------------------------------------------
    # Status bar (below Stage A/B)
    # ------------------------------------------------------------------

    def _build_status_bar(self, parent: ttk.Frame) -> None:
        bar = ttk.Frame(parent)
        bar.pack(fill=tk.X, **_PAD)

        self._btn_export_json = ttk.Button(
            bar, text="Export JSON", command=self._on_export_json, state=tk.DISABLED
        )
        self._btn_export_json.pack(side=tk.LEFT, **_PAD)

        self._btn_export_csv = ttk.Button(
            bar, text="Export CSV", command=self._on_export_csv, state=tk.DISABLED
        )
        self._btn_export_csv.pack(side=tk.LEFT, **_PAD)

        self._btn_export_zmx = ttk.Button(
            bar, text="Export .zmx", command=self._on_export_zmx, state=tk.DISABLED
        )
        self._btn_export_zmx.pack(side=tk.LEFT, **_PAD)

        ttk.Separator(bar, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=8, pady=2
        )

        self._var_filter_shape = tk.BooleanVar(value=False)
        self._chk_filter_shape = ttk.Checkbutton(
            bar, text="Hide extreme shapes",
            variable=self._var_filter_shape,
            command=self._on_filter_shapes_toggle,
            state=tk.DISABLED,
        )
        self._chk_filter_shape.pack(side=tk.LEFT, **_PAD)

        ttk.Separator(bar, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=8, pady=2
        )

        self._progress = ttk.Progressbar(bar, length=200, mode="determinate")
        self._progress.pack(side=tk.LEFT, **_PAD)

        self._var_status = tk.StringVar(value="Ready")
        ttk.Label(bar, textvariable=self._var_status).pack(
            side=tk.LEFT, fill=tk.X, expand=True, **_PAD
        )

    def _toggle_advanced(self) -> None:
        """Show/hide the advanced Seidel parameters row."""
        if self._adv_visible.get():
            self._adv_frame.pack(fill=tk.X, padx=6, pady=(0, 3))
            self._btn_adv.configure(text="▼ Advanced")
        else:
            self._adv_frame.pack_forget()
            self._btn_adv.configure(text="▶ Advanced")

    def _on_wl_preset(self) -> None:
        """Fill wavelength entries from the selected preset."""
        name = self._var_wl_preset.get()
        for preset_name, lam1, lam0, lam2 in self._WL_PRESETS:
            if preset_name == name:
                self._var_lam1_nm.set(lam1)
                self._var_lam0_nm.set(lam0)
                self._var_lam2_nm.set(lam2)
                break

    def _update_fno_display(self) -> None:
        """Refresh the F/# label and window title from current f' and D."""
        try:
            fprime = self._var_fprime.get()
            D = self._var_D.get()
            if D > 0:
                fno = fprime / D
                self._var_fno_display.set(f"F/{fno:.1f}")
                stype = self._var_type.get()
                self.title(
                    f"AutoAchromat — f/{fno:.1f} {stype} "
                    f"D={D:.0f}mm f'={fprime:.0f}mm"
                )
            else:
                self._var_fno_display.set("")
        except tk.TclError:
            self._var_fno_display.set("")

    def _on_type_changed(self) -> None:
        """Enable/disable spaced-only fields based on system type."""
        is_spaced = self._var_type.get() == "spaced"
        state = tk.NORMAL if is_spaced else tk.DISABLED
        self._ent_air_gap.configure(state=state)
        self._ent_gap_min.configure(state=state)
        self._ent_gap_max.configure(state=state)
        self._update_filter_shape_state()
        self._update_fno_display()  # refresh title with new type

    def _update_filter_shape_state(self) -> None:
        """Enable shape filter checkbox only for spaced type with results."""
        if self._var_type.get() == "spaced" and self._results:
            self._chk_filter_shape.configure(state=tk.NORMAL)
        else:
            self._var_filter_shape.set(False)
            self._chk_filter_shape.configure(state=tk.DISABLED)

    def _auto_fill_thickness(self) -> None:
        """Fill te_min / tc_min from Table 10-3 based on current D."""
        from .thickening import lookup_t_edge_min, lookup_t_center_min

        try:
            D = self._var_D.get()
        except tk.TclError:
            return
        if D <= 0:
            return
        self._var_te_min.set(round(lookup_t_edge_min(D), 2))
        self._var_tc_min.set(round(lookup_t_center_min(D), 2))

    # ------------------------------------------------------------------
    # ③ Results table
    # ------------------------------------------------------------------

    _COLUMNS = [
        ("#", "idx", 40, tk.CENTER),
        ("From", "src", 35, tk.CENTER),
        ("Glass 1", "glass1", 130, tk.W),
        ("Glass 2", "glass2", 130, tk.W),
        ("Type", "type", 55, tk.CENTER),
        ("EFL", "efl", 70, tk.E),
        ("FNO", "fno", 55, tk.E),
        ("t₁", "t1", 50, tk.E),
        ("t₂", "t2", 50, tk.E),
        ("gap", "gap", 45, tk.E),
        ("RMS [µm]", "rms", 70, tk.E),
        ("GEO [µm]", "geo", 70, tk.E),
        ("SA", "sa", 70, tk.E),
        ("LchC", "lchc", 70, tk.E),
        ("2nd Sp [µm]", "ss", 75, tk.E),
        ("TchC", "tchc", 70, tk.E),
        ("PE", "pe", 60, tk.E),
        ("α_h [ppm/K]", "alpha_h", 85, tk.E),
        ("Cost1", "cost1", 55, tk.E),
        ("Cost2", "cost2", 55, tk.E),
    ]

    def _build_result_table(self, parent: ttk.PanedWindow) -> None:
        frame = ttk.LabelFrame(self, text=" ③ Results ")
        parent.add(frame, weight=1)

        col_ids = [c[1] for c in self._COLUMNS]
        self._tree = ttk.Treeview(
            frame, columns=col_ids, show="headings", selectmode="extended"
        )

        for display, cid, width, anchor in self._COLUMNS:
            self._tree.heading(
                cid, text=display, command=lambda c=cid: self._on_sort(c)
            )
            self._tree.column(cid, width=width, anchor=anchor, minwidth=40)  # type: ignore[arg-type]

        # Scrollbars
        vsb = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=self._tree.yview)
        hsb = ttk.Scrollbar(frame, orient=tk.HORIZONTAL, command=self._tree.xview)
        self._tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self._tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        self._tree.bind("<<TreeviewSelect>>", self._on_row_select)

    # ------------------------------------------------------------------
    # ④ Detail panel
    # ------------------------------------------------------------------

    def _build_detail_panel(self, parent: ttk.PanedWindow) -> None:
        frame = ttk.LabelFrame(self, text=" ④ Selected Design Details ")
        parent.add(frame, weight=1)

        # Use grid with uniform columns to force proportional widths
        frame.columnconfigure(0, weight=1, uniform="detail")
        frame.columnconfigure(1, weight=1, uniform="detail")
        frame.columnconfigure(2, weight=1, uniform="detail")
        frame.columnconfigure(3, weight=1, uniform="detail")
        frame.rowconfigure(0, weight=1)

        # Prescription & First-Order
        left = ttk.LabelFrame(frame, text="Prescription & First-Order")
        left.grid(row=0, column=0, sticky="nsew", padx=4, pady=2)
        self._detail_rx = tk.Text(
            left, height=10, width=1, font=("Consolas", 9), state=tk.DISABLED
        )
        self._detail_rx.pack(fill=tk.BOTH, expand=True, padx=2, pady=2)

        # Aberrations
        mid = ttk.LabelFrame(frame, text="Aberrations")
        mid.grid(row=0, column=1, sticky="nsew", padx=4, pady=2)
        self._detail_ab = tk.Text(
            mid, height=10, width=1, font=("Consolas", 9), state=tk.DISABLED
        )
        self._detail_ab.pack(fill=tk.BOTH, expand=True, padx=2, pady=2)

        # Glass & Thermal
        right = ttk.LabelFrame(frame, text="Glass & Thermal")
        right.grid(row=0, column=2, sticky="nsew", padx=4, pady=2)
        self._detail_cmp = tk.Text(
            right, height=10, width=1, font=("Consolas", 9), state=tk.DISABLED
        )
        self._detail_cmp.pack(fill=tk.BOTH, expand=True, padx=2, pady=2)

        # n(λ) dispersion curve
        self._nlam_frame = ttk.LabelFrame(frame, text="Dispersion")
        self._nlam_frame.grid(row=0, column=3, sticky="nsew", padx=4, pady=2)
        self._nlam_fig = None
        self._nlam_canvas: Optional[FigureCanvasTkAgg] = None

    def _build_drawing_panel(self, parent: ttk.PanedWindow) -> None:
        outer = ttk.Frame(self)
        parent.add(outer, weight=1)
        outer.columnconfigure(0, weight=3)
        outer.columnconfigure(1, weight=2)
        outer.columnconfigure(2, weight=2)
        outer.rowconfigure(0, weight=1)

        self._draw_frame = ttk.LabelFrame(outer, text=" ⑤ 2D Optical Layout ")
        self._draw_frame.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        self._mpl_fig = None
        self._mpl_canvas: Optional[FigureCanvasTkAgg] = None

        self._cfs_frame = ttk.LabelFrame(outer, text=" Chromatic Focal Shift ")
        self._cfs_frame.grid(row=0, column=1, sticky="nsew", padx=2, pady=2)
        self._cfs_fig = None
        self._cfs_canvas: Optional[FigureCanvasTkAgg] = None

        self._seidel_frame = ttk.LabelFrame(outer, text=" Seidel per Surface ")
        self._seidel_frame.grid(row=0, column=2, sticky="nsew", padx=2, pady=2)
        self._seidel_fig = None
        self._seidel_canvas: Optional[FigureCanvasTkAgg] = None

    def _equalize_pane_heights(self) -> None:
        """Set equal initial heights for the three main panes."""
        total = self._main_pane.winfo_height()
        if total <= 1:
            # Window not yet mapped; retry after a short delay
            self.after(50, self._equalize_pane_heights)
            return
        third = total // 3
        self._main_pane.sashpos(0, third)
        self._main_pane.sashpos(1, 2 * third)

    # ===================================================================
    # Catalog management
    # ===================================================================

    def _add_catalog(self) -> None:
        paths = filedialog.askopenfilenames(
            title="Select AGF Catalog File(s)",
            filetypes=[("AGF Files", "*.agf *.AGF"), ("All Files", "*.*")],
        )
        for p in paths:
            if p and p not in self._catalog_paths:
                self._catalog_paths.append(p)
        self._refresh_catalog_display()

    def _clear_catalogs(self) -> None:
        self._catalog_paths.clear()
        self._refresh_catalog_display()

    def _refresh_catalog_display(self) -> None:
        if not self._catalog_paths:
            self._var_cat_display.set("(none)")
        else:
            names = [Path(p).stem for p in self._catalog_paths]
            self._var_cat_display.set(", ".join(names))

    # ===================================================================
    # Build Inputs from GUI fields
    # ===================================================================

    def _make_inputs(self) -> Inputs:
        return Inputs(
            lam0=self._var_lam0_nm.get() / 1000.0,
            lam1=self._var_lam1_nm.get() / 1000.0,
            lam2=self._var_lam2_nm.get() / 1000.0,
            D=self._var_D.get(),
            fprime=self._var_fprime.get(),
            C0=self._var_C0.get(),
            P0=self._var_P0.get(),
            W0=self._var_W0.get(),
            min_delta_nu=self._var_min_dnu.get(),
            max_PE=self._var_max_PE.get(),
            N=self._var_N.get(),
            system_type=self._var_type.get(),  # type: ignore[arg-type]
            air_gap=self._var_air_gap.get(),
            half_field_angle=self._var_half_field.get(),
            te_min=self._var_te_min.get(),
            tc_min=self._var_tc_min.get(),
            t1_min=self._var_t1_min.get(),
            t1_max=self._var_t1_max.get(),
            t2_min=self._var_t2_min.get(),
            t2_max=self._var_t2_max.get(),
            gap_min=self._var_gap_min.get(),
            gap_max=self._var_gap_max.get(),
            w_efl=self._var_w_efl.get(),
            w_rms=self._var_w_rms.get(),
            w_field=self._var_w_field.get(),
            maxiter=self._var_maxiter.get(),
        )

    # ===================================================================
    # Pipeline execution (background thread)
    # ===================================================================

    def _on_run(self) -> None:
        if self._running:
            return
        if not self._catalog_paths:
            messagebox.showwarning(
                "No Catalogs", "Please add at least one AGF catalog file."
            )
            return

        # Build Inputs on the main thread (tkinter vars are not thread-safe)
        try:
            inputs = self._make_inputs()
        except Exception as e:
            messagebox.showerror("Input Error", str(e))
            return

        self._running = True
        self._btn_run.configure(state=tk.DISABLED)
        self._btn_export_json.configure(state=tk.DISABLED)
        self._btn_export_csv.configure(state=tk.DISABLED)
        self._btn_export_zmx.configure(state=tk.DISABLED)
        self._btn_stage_b.configure(state=tk.DISABLED)
        self._update_filter_shape_state()
        self._var_filter_shape.set(False)
        self._results.clear()
        self._clear_tree()
        self._clear_details()
        self._var_status.set("Running design pipeline...")

        thread = threading.Thread(
            target=self._run_pipeline, args=(inputs,), daemon=True,
        )
        thread.start()

    def _run_pipeline(self, inputs: Inputs) -> None:
        """Execute full pipeline in a background thread."""
        try:

            self._set_status("Running design pipeline...")

            def _on_progress(current: int, total: int, pr: PipelineResult) -> None:
                row = ResultRow(current - 1, pr)
                self._results.append(row)
                self._insert_tree_row(row)
                self._set_progress(current, total)
                self._set_status(
                    f"Processing {current}/{total}: "
                    f"{pr.candidate.glass1} + {pr.candidate.glass2} ..."
                )

            dr = run_design(
                inputs,
                self._catalog_paths,
                on_progress=_on_progress,
            )

            n_ok = sum(1 for r in self._results if r.metrics.success)
            self._set_status(
                f"Done — {n_ok}/{len(dr.results)} succeeded  |  "
                f"Synthesis: {len(dr.candidates)} candidates  |  "
                f"{dr.synth_time_ms:.0f} ms"
            )

        except Exception as e:
            self._set_status(f"ERROR: {e}")
            self.after(
                0, lambda err=e: messagebox.showerror("Pipeline Error", str(err))
            )

        finally:
            self.after(0, self._pipeline_done)

    def _pipeline_done(self) -> None:
        self._running = False
        self._btn_run.configure(state=tk.NORMAL)
        if self._results:
            self._btn_export_json.configure(state=tk.NORMAL)
            self._btn_export_csv.configure(state=tk.NORMAL)
            self._btn_export_zmx.configure(state=tk.NORMAL)
            self._btn_stage_b.configure(state=tk.NORMAL)
            self._update_filter_shape_state()

    # ------------------------------------------------------------------
    # Shape filter
    # ------------------------------------------------------------------

    @staticmethod
    def _is_extreme_shape(row: "ResultRow") -> bool:
        """Return True if any element has a problematic shape for manufacturing.

        Criteria (any triggers rejection)
        ----------------------------------
        For spaced doublets:
        - Meniscus (both radii same sign) with convex side facing the
          air gap.  Elem 1: both R < 0; Elem 2: both R > 0.
        For all types:
        - Near-hemispherical surface  |R| < semi_aperture * 1.2
        """
        rx = row.rx
        if rx is None:
            return False

        semi_ap = rx.D / 2.0
        r_min = semi_ap * 1.2

        # Near-hemispherical check
        for elem in rx.elements:
            for R in (elem.R_front, elem.R_back):
                if abs(R) < r_min:
                    return True

        # Meniscus with convex facing the air gap (spaced only)
        if rx.system_type == "spaced" and len(rx.elements) == 2:
            e1, e2 = rx.elements
            if e1.R_front < 0 and e1.R_back < 0:
                return True
            if e2.R_front > 0 and e2.R_back > 0:
                return True

        return False

    def _on_filter_shapes_toggle(self) -> None:
        """Show/hide results with extreme shapes (non-destructive)."""
        self._refresh_tree()
        self._clear_details()
        if self._var_filter_shape.get():
            hidden = sum(1 for r in self._results if self._is_extreme_shape(r))
            shown = len(self._results) - hidden
            self._var_status.set(
                f"Shape filter ON: showing {shown}, hiding {hidden}"
            )
        else:
            self._var_status.set(
                f"Shape filter OFF: showing all {len(self._results)} results"
            )

    # ------------------------------------------------------------------
    # Stage B: thick-lens optimisation
    # ------------------------------------------------------------------

    def _on_stage_b(self) -> None:
        if self._running:
            return

        # Get selected rows from the treeview
        sel_iids = self._tree.selection()
        if not sel_iids:
            messagebox.showwarning(
                "No Selection",
                "Please select one or more rows in the results table to optimise.",
            )
            return

        # Resolve selected indices to Stage-A ResultRows with successful metrics
        selected_rows: list[ResultRow] = []
        for iid in sel_iids:
            try:
                idx = int(iid)
            except ValueError:
                continue
            row = next((r for r in self._results if r.idx == idx), None)
            if row is not None and row.stage == "A" and row.metrics.success:
                selected_rows.append(row)

        if not selected_rows:
            messagebox.showwarning(
                "No Valid Selection",
                "Please select Stage-A rows with successful evaluation.",
            )
            return

        # Build Inputs on the main thread (tkinter vars are not thread-safe)
        try:
            inputs = self._make_inputs()
        except Exception as e:
            messagebox.showerror("Input Error", str(e))
            return

        self._running = True
        self._btn_run.configure(state=tk.DISABLED)
        self._btn_stage_b.configure(state=tk.DISABLED)
        self._update_filter_shape_state()

        self._stage_b_selected = selected_rows
        thread = threading.Thread(
            target=self._run_stage_b_bg, args=(inputs,), daemon=True,
        )
        thread.start()

    def _run_stage_b_bg(self, inputs: Inputs) -> None:
        """Execute Stage B optimisation in a background thread."""
        try:
            selected_results = [r.result for r in self._stage_b_selected]
            n_sel = len(selected_results)

            # Build lookup: (glass1, glass2) → Stage A display index
            _src_map: dict[tuple[str, str], int] = {}
            for r in self._stage_b_selected:
                key = (r.cand.glass1, r.cand.glass2)
                if key not in _src_map:
                    _src_map[key] = r.idx

            def _on_progress(current: int, total: int, pr: PipelineResult) -> None:
                idx = len(self._results)
                src = _src_map.get((pr.candidate.glass1, pr.candidate.glass2), -1)
                row = ResultRow(idx, pr, stage="B", src_idx=src)
                self._results.append(row)
                self._insert_tree_row(row)
                self._set_progress(current, total)
                self._set_status(
                    f"Stage B {current}/{total}: "
                    f"{pr.candidate.glass1} + {pr.candidate.glass2} ..."
                )

            self._set_status(f"Stage B: optimising {n_sel} selected candidates...")
            results_b = run_stage_b(
                selected_results,
                inputs,
                top_n=n_sel,
                on_progress=_on_progress,
            )

            n_ok = sum(1 for r in results_b if r.metrics.success)
            self._set_status(
                f"Stage B done — {n_ok}/{len(results_b)} optimised successfully"
            )

        except Exception as e:
            self._set_status(f"Stage B ERROR: {e}")
            self.after(
                0, lambda err=e: messagebox.showerror("Stage B Error", str(err))
            )

        finally:
            self.after(0, self._stage_b_done)

    def _stage_b_done(self) -> None:
        self._running = False
        self._btn_run.configure(state=tk.NORMAL)
        self._btn_stage_b.configure(state=tk.NORMAL)
        if self._results:
            self._update_filter_shape_state()

    # ===================================================================
    # Thread-safe UI updates
    # ===================================================================

    def _set_status(self, msg: str) -> None:
        self.after(0, lambda: self._var_status.set(msg))

    def _set_progress(self, current: int, total: int) -> None:
        def _update():
            self._progress["maximum"] = total
            self._progress["value"] = current

        self.after(0, _update)

    # ===================================================================
    # Treeview helpers
    # ===================================================================

    def _clear_tree(self) -> None:
        for item in self._tree.get_children():
            self._tree.delete(item)

    def _row_values(self, row: ResultRow) -> tuple:
        m = row.metrics
        rx = row.rx
        t1 = _fmt(rx.elements[0].t_center, ".2f") if rx else "—"
        t2 = _fmt(rx.elements[1].t_center, ".2f") if rx else "—"
        gap = _fmt(rx.air_gap, ".2f") if rx and rx.air_gap is not None else "—"
        tp = "Cem" if row.cand.system_type == "cemented" else "Sp"
        if row.stage == "B":
            tp += "★"

        th = row.cand.thermal
        alpha_h = (
            _fmt(th.alpha_housing_required * 1e6, ".1f")
            if th is not None and th.alpha_housing_required is not None
            else "—"
        )

        src = str(row.src_idx + 1) if row.src_idx >= 0 else ""

        return (
            row.idx + 1,
            src,
            f"{row.cand.catalog1}:{row.cand.glass1}",
            f"{row.cand.catalog2}:{row.cand.glass2}",
            tp,
            _fmt(m.efl, ".1f"),
            _fmt(m.fno, ".3f"),
            t1,
            t2,
            gap,
            _fmt(m.rms_spot_radius, ".2f"),
            _fmt(m.geo_spot_radius, ".2f"),
            _fmt(m.SA, ".3f"),
            _fmt(m.LchC, ".3f"),
            _fmt(m.secondary_spectrum * 1000, ".1f") if m.secondary_spectrum is not None else "—",
            _fmt(m.TchC, ".3f"),
            _fmt(m.PE, ".2f"),
            alpha_h,
            _fmt(row.cand.cost1, ".2f"),
            _fmt(row.cand.cost2, ".2f"),
        )

    def _insert_tree_row(self, row: ResultRow) -> None:
        values = self._row_values(row)
        tag = "ok" if row.metrics.success else "fail"

        def _insert():
            iid = self._tree.insert(
                "", tk.END, iid=str(row.idx), values=values, tags=(tag,)
            )
            self._tree.tag_configure("fail", foreground="#999999")
            self._tree.see(iid)

        self.after(0, _insert)

    def _refresh_tree(self) -> None:
        """Re-populate tree with current sort order, respecting shape filter."""
        self._clear_tree()
        filtering = self._var_filter_shape.get()
        for row in self._results:
            if filtering and self._is_extreme_shape(row):
                continue
            values = self._row_values(row)
            tag = "ok" if row.metrics.success else "fail"
            self._tree.insert("", tk.END, iid=str(row.idx), values=values, tags=(tag,))
        self._tree.tag_configure("fail", foreground="#999999")

    # ===================================================================
    # Sorting
    # ===================================================================

    _SORT_KEY_MAP = {
        "idx": lambda r: r.idx,
        "src": lambda r: r.src_idx if r.src_idx >= 0 else 1e9,
        "glass1": lambda r: f"{r.cand.catalog1}:{r.cand.glass1}",
        "glass2": lambda r: f"{r.cand.catalog2}:{r.cand.glass2}",
        "type": lambda r: r.cand.system_type,
        "efl": lambda r: r.metrics.efl if r.metrics.efl is not None else 1e9,
        "fno": lambda r: r.metrics.fno if r.metrics.fno is not None else 1e9,
        "t1": lambda r: r.rx.elements[0].t_center if r.rx else 1e9,
        "t2": lambda r: r.rx.elements[1].t_center if r.rx else 1e9,
        "gap": lambda r: r.rx.air_gap if r.rx and r.rx.air_gap is not None else 1e9,
        "rms": lambda r: (
            r.metrics.rms_spot_radius if r.metrics.rms_spot_radius is not None else 1e9
        ),
        "geo": lambda r: (
            r.metrics.geo_spot_radius if r.metrics.geo_spot_radius is not None else 1e9
        ),
        "sa": lambda r: abs(r.metrics.SA) if r.metrics.SA is not None else 1e9,
        "lchc": lambda r: abs(r.metrics.LchC) if r.metrics.LchC is not None else 1e9,
        "ss": lambda r: (
            abs(r.metrics.secondary_spectrum)
            if r.metrics.secondary_spectrum is not None
            else 1e9
        ),
        "tchc": lambda r: abs(r.metrics.TchC) if r.metrics.TchC is not None else 1e9,
        "pe": lambda r: abs(r.metrics.PE) if r.metrics.PE is not None else 1e9,
        "alpha_h": lambda r: (
            r.cand.thermal.alpha_housing_required
            if r.cand.thermal is not None
            and r.cand.thermal.alpha_housing_required is not None
            else 1e9
        ),
        "cost1": lambda r: r.cand.cost1 if r.cand.cost1 is not None else 1e9,
        "cost2": lambda r: r.cand.cost2 if r.cand.cost2 is not None else 1e9,
    }

    def _on_sort(self, col: str) -> None:
        if self._sort_col == col:
            self._sort_reverse = not self._sort_reverse
        else:
            self._sort_col = col
            self._sort_reverse = False

        key_fn = self._SORT_KEY_MAP.get(col, lambda r: r.idx)
        self._results.sort(key=key_fn, reverse=self._sort_reverse)
        self._refresh_tree()

    # ===================================================================
    # Row selection → detail panel
    # ===================================================================

    def _on_row_select(self, _event: tk.Event) -> None:
        sel = self._tree.selection()
        if not sel:
            return

        # Update Stage B button text with selection count
        n_valid = 0
        for iid in sel:
            try:
                idx = int(iid)
            except ValueError:
                continue
            row = next((r for r in self._results if r.idx == idx), None)
            if row is not None and row.stage == "A" and row.metrics.success:
                n_valid += 1
        if n_valid > 0:
            self._btn_stage_b.configure(text=f" ⚡ Optimize ({n_valid}) ")
        else:
            self._btn_stage_b.configure(text=" ⚡ Optimize ")

        # Fill detail for first selected row
        try:
            idx = int(sel[0])
        except (ValueError, IndexError):
            return

        row = next((r for r in self._results if r.idx == idx), None)
        if row is None:
            return

        self._fill_detail(row)

    def _clear_details(self) -> None:
        for widget in (self._detail_rx, self._detail_ab, self._detail_cmp):
            widget.configure(state=tk.NORMAL)
            widget.delete("1.0", tk.END)
            widget.configure(state=tk.DISABLED)
        self._clear_drawing()
        self._clear_cfs()
        self._clear_seidel()
        self._clear_nlam()

    def _fill_detail(self, row: ResultRow) -> None:
        self._fill_prescription(row)
        self._fill_aberrations(row)
        self._fill_comparison(row)
        self._fill_drawing(row)
        self._fill_cfs(row)
        self._fill_seidel(row)
        self._fill_nlam(row)

    def _fill_prescription(self, row: ResultRow) -> None:
        """Fill the Prescription text widget with surface table."""
        w = self._detail_rx
        w.configure(state=tk.NORMAL)
        w.delete("1.0", tk.END)

        rx = row.rx
        if rx is None:
            w.insert(tk.END, "  (thickening failed)\n")
            w.configure(state=tk.DISABLED)
            return

        e1, e2 = rx.elements
        lines = []
        lines.append(f"{'Surf':>4}  {'R [mm]':>10}  {'t [mm]':>8}  {'Material':<12}")
        lines.append("─" * 44)

        # Object
        lines.append(f" OBJ  {'∞':>10}  {'∞':>8}  Air")

        if rx.system_type == "cemented":
            lines.append(
                f" STO  {e1.R_front:10.3f}  {e1.t_center:8.3f}  n={e1.nd:.4f} ν={e1.vd:.1f}"
            )
            lines.append(
                f"   3  {e1.R_back:10.3f}  {e2.t_center:8.3f}  n={e2.nd:.4f} ν={e2.vd:.1f}"
            )
            _bfd = (
                row.metrics.bfd if row.metrics.bfd is not None else rx.back_focus_guess
            )
            lines.append(f"   4  {e2.R_back:10.3f}  {_bfd:8.3f}  Air")
        else:
            lines.append(
                f" STO  {e1.R_front:10.3f}  {e1.t_center:8.3f}  n={e1.nd:.4f} ν={e1.vd:.1f}"
            )
            lines.append(f"   3  {e1.R_back:10.3f}  {rx.air_gap or 0:8.3f}  Air")
            lines.append(
                f"   4  {e2.R_front:10.3f}  {e2.t_center:8.3f}  n={e2.nd:.4f} ν={e2.vd:.1f}"
            )
            _bfd = (
                row.metrics.bfd if row.metrics.bfd is not None else rx.back_focus_guess
            )
            lines.append(f"   5  {e2.R_back:10.3f}  {_bfd:8.3f}  Air")

        lines.append(f" IMG  {'∞':>10}  {'—':>8}")

        # First-order data
        m = row.metrics
        lines.append("")
        lines.append("  First-order")
        lines.append("  ───────────────────────")
        lines.append(f"  EFL (optiland): {_fmt(m.efl, '.3f')} mm")
        lines.append(f"  FNO: {_fmt(m.fno, '.4f')}")
        lines.append(f"  BFD: {_fmt(m.bfd, '.3f')} mm")
        if rx.actual_efl is not None:
            dev_pct = (rx.efl_deviation or 0.0) * 100.0
            lines.append(f"  EFL (ABCD):     {rx.actual_efl:.3f} mm")
            lines.append(f"  EFL deviation:  {dev_pct:+.3f} %")

        # Thin-lens synthesis parameters
        cand = row.cand
        lines.append("")
        lines.append("  Synthesis Parameters")
        lines.append("  ───────────────────────")
        lines.append(f"  φ₁ = {cand.phi1:.6f}")
        lines.append(f"  φ₂ = {cand.phi2:.6f}")
        if cand.Q is not None:
            lines.append(f"  Q  = {cand.Q:.6f}")
        if cand.Q1 is not None:
            lines.append(f"  Q₁ = {cand.Q1:.6f}")
            lines.append(f"  Q₂ = {cand.Q2:.6f}")
        if cand.W is not None:
            lines.append(f"  W  = {cand.W:.6f}")
        lines.append(f"  PE = {_fmt(cand.PE, '.4f')}")

        w.insert(tk.END, "\n".join(lines))
        w.configure(state=tk.DISABLED)

    def _fill_aberrations(self, row: ResultRow) -> None:
        """Fill the Aberrations text widget."""
        w = self._detail_ab
        w.configure(state=tk.NORMAL)
        w.delete("1.0", tk.END)

        m = row.metrics
        if not m.success:
            w.insert(tk.END, f"  Error: {m.error_msg}\n")
            w.configure(state=tk.DISABLED)
            return

        lines = []
        lines.append("  Seidel Coefficients")
        lines.append("  ───────────────────────")
        lines.append(f"  SA  (spherical):  {_fmt(m.SA, '.4f')}")
        lines.append(f"  CC  (coma):       {_fmt(m.CC, '.4f')}")
        lines.append(f"  AC  (astigmatism):{_fmt(m.AC, '.4f')}")
        lines.append(f"  PC  (Petzval):    {_fmt(m.PC, '.4f')}")
        lines.append(f"  DC  (distortion): {_fmt(m.DC, '.4f')}")
        lines.append("")
        lines.append("  Chromatic")
        lines.append("  ───────────────────────")
        lines.append(f"  LchC (longit.):   {_fmt(m.LchC, '.4f')}")
        lines.append(f"  TchC (transv.):   {_fmt(m.TchC, '.4f')}")
        ss_um = m.secondary_spectrum * 1000.0 if m.secondary_spectrum is not None else None
        lines.append(f"  2nd spectrum:     {_fmt(ss_um, '.1f')} µm")
        lines.append("")
        lines.append("  Spot Diagram")
        lines.append("  ───────────────────────")
        lines.append(f"  RMS spot: {_fmt(m.rms_spot_radius, '.2f')} µm")
        lines.append(f"  GEO spot: {_fmt(m.geo_spot_radius, '.2f')} µm")
        if m.fno is not None:
            try:
                lam0_um = self._make_inputs().lam0
            except Exception:
                lam0_um = 0.58756
            airy_um = 1.22 * lam0_um * m.fno
            lines.append(f"  Airy rad: {airy_um:.2f} µm")
            if m.rms_spot_radius is not None and airy_um > 0:
                ratio = m.rms_spot_radius / airy_um
                status = "diffr. limited" if ratio < 1.0 else f"{ratio:.1f}× Airy"
                lines.append(f"  RMS/Airy: {status}")

        w.insert(tk.END, "\n".join(lines))
        w.configure(state=tk.DISABLED)

    def _fill_comparison(self, row: ResultRow) -> None:
        """Fill the Glass & Geometry panel."""
        w = self._detail_cmp
        w.configure(state=tk.NORMAL)
        w.delete("1.0", tk.END)

        cand = row.cand
        lines = []

        # Glass properties
        try:
            lam0_nm = self._var_lam0_nm.get()
        except tk.TclError:
            lam0_nm = 587.56
        lines.append(f"  Glass Properties (at λ₀={lam0_nm:.1f}nm)")
        lines.append("  ─────────────────────────")

        # Compute partial dispersion for each glass
        P_vals = [None, None]
        try:
            from .optics import refractive_index
            from .glass_reader import Glass
            inputs = self._make_inputs()
            for i, (fid, cd) in enumerate([
                (cand.formula_id1, cand.cd1),
                (cand.formula_id2, cand.cd2),
            ]):
                if fid is not None and cd:
                    cat = [cand.catalog1, cand.catalog2][i]
                    name = [cand.glass1, cand.glass2][i]
                    g = Glass(name=name, catalog=cat, formula_id=fid, cd=cd)
                    n0 = refractive_index(g, inputs.lam0)
                    n1 = refractive_index(g, inputs.lam1)
                    n2 = refractive_index(g, inputs.lam2)
                    dn = n1 - n2
                    if abs(dn) > 1e-10:
                        P_vals[i] = (n0 - n1) / dn
        except Exception:
            pass

        p1_str = f"P={P_vals[0]:.4f}" if P_vals[0] is not None else ""
        p2_str = f"P={P_vals[1]:.4f}" if P_vals[1] is not None else ""

        lines.append(f"  {cand.glass1}:")
        lines.append(f"    n₀={cand.n1:.6f}  ν={cand.nu1:.2f}  {p1_str}")
        lines.append(f"  {cand.glass2}:")
        lines.append(f"    n₀={cand.n2:.6f}  ν={cand.nu2:.2f}  {p2_str}")
        lines.append(f"  Δν={abs(cand.nu1 - cand.nu2):.2f}")
        if P_vals[0] is not None and P_vals[1] is not None:
            lines.append(f"  ΔP={abs(P_vals[0] - P_vals[1]):.4f}")

        # Internal transmittance at design wavelengths
        lines.append("")
        lines.append("  Transmittance (internal)")
        lines.append("  ─────────────────────────")
        try:
            inputs = self._make_inputs()
            for name, trans_data in [
                (cand.glass1, cand.trans1),
                (cand.glass2, cand.trans2),
            ]:
                if not trans_data:
                    lines.append(f"  {name}: (no IT data)")
                    continue
                lines.append(f"  {name}:")
                for lam in [inputs.lam1, inputs.lam0, inputs.lam2]:
                    t_val = _interp_transmittance(trans_data, lam)
                    lam_nm = lam * 1000
                    t_str = f"{t_val:.4f}" if t_val is not None else "—"
                    lines.append(f"    {lam_nm:>7.1f} nm : {t_str}")
        except Exception:
            lines.append("  (error)")

        # Thermal analysis
        th = cand.thermal
        lines.append("")
        lines.append("  Thermal Analysis")
        lines.append("  ─────────────────────────")
        if th is None or not th.thermal_data_available:
            lines.append("  (no TD/ED data in catalog)")
        else:
            lines.append(f"  dn/dT₁:  {_fmt(th.dn_dT_1, '.3e')} /K")
            lines.append(f"  dn/dT₂:  {_fmt(th.dn_dT_2, '.3e')} /K")
            v1_ppm = th.V1 * 1e6 if th.V1 is not None else None
            v2_ppm = th.V2 * 1e6 if th.V2 is not None else None
            lines.append(f"  V₁:      {_fmt(v1_ppm, '.2f')} ppm/K")
            lines.append(f"  V₂:      {_fmt(v2_ppm, '.2f')} ppm/K")
            dphi_ppm = th.dphi_dT_norm * 1e6 if th.dphi_dT_norm is not None else None
            lines.append(f"  dΦ/dT:   {_fmt(dphi_ppm, '.2f')} ppm/K")
            alpha_ppm = (
                th.alpha_housing_required * 1e6
                if th.alpha_housing_required is not None
                else None
            )
            lines.append(f"  α_h req: {_fmt(alpha_ppm, '.2f')} ppm/K")

        # Cost
        if cand.cost1 is not None or cand.cost2 is not None:
            lines.append("")
            lines.append("  Cost")
            lines.append("  ─────────────────────────")
            lines.append(f"  Glass 1: {_fmt(cand.cost1, '.2f')}")
            lines.append(f"  Glass 2: {_fmt(cand.cost2, '.2f')}")

        w.insert(tk.END, "\n".join(lines))
        w.configure(state=tk.DISABLED)

    # ===================================================================
    # 2D optical layout drawing
    # ===================================================================

    def _clear_drawing(self) -> None:
        """Release current matplotlib figure and clear the drawing frame."""
        if self._mpl_fig is not None:
            plt.close(self._mpl_fig)
            self._mpl_fig = None
        self._mpl_canvas = None
        for w in self._draw_frame.winfo_children():
            w.destroy()

    def _fill_drawing(self, row: ResultRow) -> None:
        """Generate and embed the 2D optical layout for the selected design."""
        self._clear_drawing()

        if row.rx is None or not row.metrics.success:
            ttk.Label(
                self._draw_frame,
                text="(no layout — evaluation failed)",
                foreground="#999999",
            ).pack(expand=True)
            return

        try:
            op = build_optic_from_prescription(row.rx)
            if op is None:
                ttk.Label(
                    self._draw_frame,
                    text="(builder returned None)",
                    foreground="#999999",
                ).pack(expand=True)
                return

            fig, _ = op.draw(num_rays=5, figsize=(12, 3))
            fig.tight_layout()
            self._mpl_fig = fig

            canvas = FigureCanvasTkAgg(fig, master=self._draw_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            self._mpl_canvas = canvas

        except Exception as exc:
            ttk.Label(
                self._draw_frame,
                text=f"Drawing error: {type(exc).__name__}: {exc}",
                foreground="red",
                wraplength=300,
            ).pack(expand=True)

    # ===================================================================
    # Chromatic focal shift plot
    # ===================================================================

    def _clear_cfs(self) -> None:
        """Release chromatic focal shift figure."""
        if self._cfs_fig is not None:
            plt.close(self._cfs_fig)
            self._cfs_fig = None
        self._cfs_canvas = None
        for w in self._cfs_frame.winfo_children():
            w.destroy()

    def _fill_cfs(self, row: ResultRow) -> None:
        """Plot chromatic focal shift curve for the selected design."""
        from .optiland_bridge.evaluator import chromatic_focal_shift

        self._clear_cfs()

        if row.rx is None or not row.metrics.success:
            ttk.Label(
                self._cfs_frame,
                text="(no data)",
                foreground="#999999",
            ).pack(expand=True)
            return

        inputs = self._make_inputs()
        result = chromatic_focal_shift(row.rx, row.cand, inputs)
        if result is None:
            ttk.Label(
                self._cfs_frame,
                text="(no dispersion data)",
                foreground="#999999",
            ).pack(expand=True)
            return

        wavelengths_um, delta_bfd_mm = result
        wavelengths_nm = [w * 1000.0 for w in wavelengths_um]
        delta_bfd = [d * 1000.0 for d in delta_bfd_mm]  # mm → µm

        try:
            fig, ax = plt.subplots(figsize=(5, 3))
            ax.plot(wavelengths_nm, delta_bfd, "b-", linewidth=1.5)
            ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")

            # Mark the three design wavelengths
            for lam_var, marker, label in [
                (inputs.lam1, "v", r"$\lambda_1$"),
                (inputs.lam0, "o", r"$\lambda_0$"),
                (inputs.lam2, "^", r"$\lambda_2$"),
            ]:
                bfd_val = None
                for wl, db in zip(wavelengths_um, delta_bfd):
                    if abs(wl - lam_var) < 1e-6:
                        bfd_val = db
                        break
                if bfd_val is None:
                    from .optiland_bridge.evaluator import _bfd_at_wavelength
                    bfd0 = _bfd_at_wavelength(row.rx, row.cand, inputs.lam0)
                    bfd_v = _bfd_at_wavelength(row.rx, row.cand, lam_var)
                    if bfd0 is not None and bfd_v is not None:
                        bfd_val = (bfd_v - bfd0) * 1000.0  # mm → µm
                if bfd_val is not None:
                    ax.plot(lam_var * 1000, bfd_val, marker, color="red", markersize=6)

            ax.set_xlabel("Wavelength [nm]")
            ax.set_ylabel("Focal shift [µm]")
            ax.set_title(
                f"{row.cand.glass1} + {row.cand.glass2}",
                fontsize=9,
            )
            fig.tight_layout()

            self._cfs_fig = fig
            canvas = FigureCanvasTkAgg(fig, master=self._cfs_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            self._cfs_canvas = canvas

        except Exception as exc:
            ttk.Label(
                self._cfs_frame,
                text=f"Plot error: {type(exc).__name__}: {exc}",
                foreground="red",
                wraplength=200,
            ).pack(expand=True)

    # ===================================================================
    # n(λ) dispersion curve
    # ===================================================================

    def _clear_nlam(self) -> None:
        if self._nlam_fig is not None:
            plt.close(self._nlam_fig)
            self._nlam_fig = None
        self._nlam_canvas = None
        for w in self._nlam_frame.winfo_children():
            w.destroy()

    def _fill_nlam(self, row: ResultRow) -> None:
        """Plot n(λ) for both glasses of the selected design."""
        from .optics import refractive_index
        from .glass_reader import Glass

        self._clear_nlam()

        cand = row.cand
        if not cand.cd1 or not cand.cd2:
            ttk.Label(
                self._nlam_frame, text="(no dispersion data)", foreground="#999999",
            ).pack(expand=True)
            return
        if cand.formula_id1 is None or cand.formula_id2 is None:
            ttk.Label(
                self._nlam_frame, text="(no formula)", foreground="#999999",
            ).pack(expand=True)
            return

        try:
            inputs = self._make_inputs()
            lam_min = min(inputs.lam1, inputs.lam2)
            lam_max = max(inputs.lam1, inputs.lam2)
            margin = (lam_max - lam_min) * 0.1
            lam_lo = lam_min - margin
            lam_hi = lam_max + margin

            g1 = Glass(name=cand.glass1, catalog=cand.catalog1,
                       formula_id=cand.formula_id1, cd=cand.cd1)
            g2 = Glass(name=cand.glass2, catalog=cand.catalog2,
                       formula_id=cand.formula_id2, cd=cand.cd2)

            n_pts = 30
            wls = [lam_lo + (lam_hi - lam_lo) * i / (n_pts - 1) for i in range(n_pts)]
            wls_nm = [w * 1000 for w in wls]

            # Compute n(λ) - n(λ₀) for both glasses
            try:
                n1_ref = refractive_index(g1, inputs.lam0)
                n2_ref = refractive_index(g2, inputs.lam0)
            except Exception:
                n1_ref = n2_ref = 0.0

            dn1_vals, dn2_vals = [], []
            for w in wls:
                try:
                    dn1_vals.append(refractive_index(g1, w) - n1_ref)
                except Exception:
                    dn1_vals.append(float("nan"))
                try:
                    dn2_vals.append(refractive_index(g2, w) - n2_ref)
                except Exception:
                    dn2_vals.append(float("nan"))

            fig, ax = plt.subplots(figsize=(4, 3))
            ax.plot(wls_nm, dn1_vals, "b-", linewidth=1.2,
                    label=cand.glass1)
            ax.plot(wls_nm, dn2_vals, "r-", linewidth=1.2,
                    label=cand.glass2)
            ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")

            for lam_var, marker in [
                (inputs.lam1, "v"), (inputs.lam0, "o"), (inputs.lam2, "^"),
            ]:
                lnm = lam_var * 1000
                try:
                    ax.plot(lnm, refractive_index(g1, lam_var) - n1_ref,
                            marker, color="blue", markersize=4)
                    ax.plot(lnm, refractive_index(g2, lam_var) - n2_ref,
                            marker, color="red", markersize=4)
                except Exception:
                    pass

            ax.set_xlabel("Wavelength [nm]", fontsize=7)
            ax.set_ylabel("n(λ) − n(λ₀)", fontsize=7)
            ax.legend(fontsize=6)
            ax.tick_params(labelsize=6)
            fig.tight_layout()

            self._nlam_fig = fig
            canvas = FigureCanvasTkAgg(fig, master=self._nlam_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            self._nlam_canvas = canvas

        except Exception as exc:
            ttk.Label(
                self._nlam_frame,
                text=f"Plot error: {type(exc).__name__}: {exc}",
                foreground="red", wraplength=200,
            ).pack(expand=True)

    # ===================================================================
    # Seidel per-surface bar chart
    # ===================================================================

    def _clear_seidel(self) -> None:
        if self._seidel_fig is not None:
            plt.close(self._seidel_fig)
            self._seidel_fig = None
        self._seidel_canvas = None
        for w in self._seidel_frame.winfo_children():
            w.destroy()

    def _fill_seidel(self, row: ResultRow) -> None:
        """Plot per-surface Seidel aberration bar chart."""
        self._clear_seidel()

        m = row.metrics
        if not m.success or not m.SA_per_surf:
            ttk.Label(
                self._seidel_frame, text="(no data)", foreground="#999999",
            ).pack(expand=True)
            return

        try:
            import numpy as _np

            n_surf = len(m.SA_per_surf)
            colors = plt.cm.tab10.colors[:n_surf]  # type: ignore[attr-defined]
            surf_labels = [f"S{s + 1}" for s in range(n_surf)]

            fig, (ax_sa, ax_rest) = plt.subplots(
                2, 1, figsize=(5, 3), height_ratios=[1, 1],
            )

            # Top: SA + LchC (large values, grouped by aberration)
            top_labels = ["SA", "LchC"]
            top_data = [m.SA_per_surf, m.LchC_per_surf]
            x_top = _np.arange(len(top_labels))
            width_top = 0.8 / n_surf

            for s in range(n_surf):
                vals = [top_data[a][s] if s < len(top_data[a]) else 0
                        for a in range(len(top_labels))]
                offset = (s - (n_surf - 1) / 2) * width_top
                ax_sa.bar(
                    x_top + offset, vals, width_top,
                    label=surf_labels[s], color=colors[s],
                )

            ax_sa.set_xticks(x_top)
            ax_sa.set_xticklabels(top_labels)
            ax_sa.axhline(0, color="gray", linewidth=0.5)

            # Total sum annotations
            for a_idx, (label, total) in enumerate([
                ("SA", m.SA), ("LchC", m.LchC),
            ]):
                if total is not None:
                    ax_sa.annotate(
                        f"Σ={total:.3f}", xy=(a_idx, total),
                        fontsize=6, ha="center", va="bottom" if total >= 0 else "top",
                        color="black",
                    )
                    ax_sa.plot(
                        [a_idx - 0.4, a_idx + 0.4], [total, total],
                        "k--", linewidth=0.8,
                    )

            ax_sa.set_ylabel("Coeff.", fontsize=8)
            ax_sa.legend(fontsize=7, ncol=n_surf)
            ax_sa.tick_params(labelsize=7)
            ax_sa.set_title(
                f"{row.cand.glass1} + {row.cand.glass2}", fontsize=9,
            )

            # Bottom: CC, AC, PC, DC, TchC (small values, shared scale)
            rest_labels = ["CC", "AC", "PC", "DC", "TchC"]
            rest_data = [
                m.CC_per_surf, m.AC_per_surf,
                m.PC_per_surf, m.DC_per_surf,
                m.TchC_per_surf,
            ]
            x = _np.arange(len(rest_labels))
            width = 0.8 / n_surf

            for s in range(n_surf):
                vals = [rest_data[a][s] if s < len(rest_data[a]) else 0
                        for a in range(len(rest_labels))]
                offset = (s - (n_surf - 1) / 2) * width
                ax_rest.bar(
                    x + offset, vals, width,
                    label=surf_labels[s], color=colors[s],
                )

            ax_rest.set_xticks(x)
            ax_rest.set_xticklabels(rest_labels)
            ax_rest.axhline(0, color="gray", linewidth=0.5)
            ax_rest.set_ylabel("Coeff.", fontsize=8)
            ax_rest.tick_params(labelsize=7)
            fig.align_ylabels([ax_sa, ax_rest])
            fig.tight_layout()

            self._seidel_fig = fig
            canvas = FigureCanvasTkAgg(fig, master=self._seidel_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            self._seidel_canvas = canvas

        except Exception as exc:
            ttk.Label(
                self._seidel_frame,
                text=f"Plot error: {type(exc).__name__}: {exc}",
                foreground="red", wraplength=200,
            ).pack(expand=True)

    # ===================================================================
    # Export
    # ===================================================================

    def _export_data(self) -> list[dict]:
        """Build JSON-serialisable list from results."""
        return [row.result.to_dict() for row in self._results]

    def _on_export_json(self) -> None:
        path = filedialog.asksaveasfilename(
            title="Export Results as JSON",
            defaultextension=".json",
            filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")],
        )
        if not path:
            return
        data = self._export_data()
        Path(path).write_text(
            json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        self._var_status.set(f"Exported {len(data)} results to {Path(path).name}")

    def _on_export_csv(self) -> None:
        path = filedialog.asksaveasfilename(
            title="Export Results as CSV",
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")],
        )
        if not path:
            return

        header = [c[0] for c in self._COLUMNS]
        lines = [",".join(header)]
        for row in self._results:
            vals = self._row_values(row)
            # Escape commas in values
            escaped = []
            for v in vals:
                s = str(v)
                if "," in s:
                    s = f'"{s}"'
                escaped.append(s)
            lines.append(",".join(escaped))

        Path(path).write_text("\n".join(lines), encoding="utf-8")
        self._var_status.set(
            f"Exported {len(self._results)} results to {Path(path).name}"
        )

    def _on_export_zmx(self) -> None:
        """Export the selected design as a Zemax .zmx file."""
        sel = self._tree.selection()
        if not sel:
            messagebox.showwarning(
                "No Selection", "Please select a row to export."
            )
            return
        try:
            idx = int(sel[0])
        except (ValueError, IndexError):
            return
        row = next((r for r in self._results if r.idx == idx), None)
        if row is None or row.rx is None:
            messagebox.showwarning(
                "No Prescription", "Selected design has no thick-lens prescription."
            )
            return

        path = filedialog.asksaveasfilename(
            title="Export as Zemax .zmx",
            defaultextension=".zmx",
            initialfile=f"{row.cand.glass1}_{row.cand.glass2}.zmx",
            filetypes=[("Zemax Files", "*.zmx"), ("All Files", "*.*")],
        )
        if not path:
            return

        from .zemax_export import export_zmx

        try:
            inputs = self._make_inputs()
            export_zmx(row.rx, row.cand, inputs, path)
            self._var_status.set(f"Exported .zmx to {Path(path).name}")
        except Exception as e:
            messagebox.showerror("Export Error", str(e))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    app = AutoAchromatGUI()
    app.mainloop()


if __name__ == "__main__":
    main()

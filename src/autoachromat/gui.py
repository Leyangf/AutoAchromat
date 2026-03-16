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


# ---------------------------------------------------------------------------
# Result row  (one per candidate)
# ---------------------------------------------------------------------------


class ResultRow:
    """Bundles a PipelineResult with an index for display."""

    __slots__ = ("result", "idx", "stage")

    def __init__(self, idx: int, result: PipelineResult, stage: str = "A"):
        self.idx = idx
        self.result = result
        self.stage = stage

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
        frame = ttk.LabelFrame(parent, text=" Global ")
        frame.pack(fill=tk.X, **_PAD)

        # Use uniform column weights so entries align
        for c in (1, 3, 5, 7, 9):
            frame.columnconfigure(c, weight=1, uniform="val")

        r = 0
        ttk.Label(frame, text="f' [mm]:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        self._var_fprime = tk.DoubleVar(value=200.0)
        ttk.Entry(frame, textvariable=self._var_fprime, width=9).grid(
            row=r, column=1, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="D [mm]:").grid(row=r, column=2, sticky=tk.E, **_PAD)
        self._var_D = tk.DoubleVar(value=50.0)
        ttk.Entry(frame, textvariable=self._var_D, width=9).grid(
            row=r, column=3, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="λ₀ [µm]:").grid(row=r, column=4, sticky=tk.E, **_PAD)
        self._var_lam0 = tk.DoubleVar(value=0.58756)
        ttk.Entry(frame, textvariable=self._var_lam0, width=9).grid(
            row=r, column=5, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="λ₁ [µm]:").grid(row=r, column=6, sticky=tk.E, **_PAD)
        self._var_lam1 = tk.DoubleVar(value=0.48613)
        ttk.Entry(frame, textvariable=self._var_lam1, width=9).grid(
            row=r, column=7, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="λ₂ [µm]:").grid(row=r, column=8, sticky=tk.E, **_PAD)
        self._var_lam2 = tk.DoubleVar(value=0.65627)
        ttk.Entry(frame, textvariable=self._var_lam2, width=9).grid(
            row=r, column=9, sticky=tk.W, **_PAD
        )

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
        ttk.Entry(frame, textvariable=self._var_N, width=9).grid(
            row=r, column=5, sticky=tk.W, **_PAD
        )

        ttk.Label(frame, text="Catalogs:").grid(row=r, column=6, sticky=tk.E, **_PAD)
        cat_frame = ttk.Frame(frame)
        cat_frame.grid(row=r, column=7, columnspan=3, sticky=tk.W, **_PAD)
        self._var_cat_display = tk.StringVar(value="(none)")
        ttk.Label(
            cat_frame, textvariable=self._var_cat_display, width=24, relief="sunken"
        ).pack(side=tk.LEFT)
        ttk.Button(cat_frame, text="+", width=3, command=self._add_catalog).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Button(cat_frame, text="−", width=3, command=self._clear_catalogs).pack(
            side=tk.LEFT, padx=2
        )
        self._btn_load = ttk.Button(
            cat_frame, text="Load", width=5, command=self._on_load_config
        )
        self._btn_load.pack(side=tk.LEFT, padx=(8, 2))

    # ------------------------------------------------------------------
    # Stage A (left) and Stage B (right)
    # ------------------------------------------------------------------

    def _build_stage_a_frame(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text=" Stage A: Thin-Lens Synthesis ")
        frame.grid(row=0, column=0, sticky="nsew", **_PAD)

        r = 0
        ttk.Label(frame, text="C₀:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        self._var_C0 = tk.DoubleVar(value=0.0)
        ttk.Entry(frame, textvariable=self._var_C0, width=8).grid(
            row=r, column=1, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="P₀:").grid(row=r, column=2, sticky=tk.E, **_PAD)
        self._var_P0 = tk.DoubleVar(value=0.0)
        ttk.Entry(frame, textvariable=self._var_P0, width=8).grid(
            row=r, column=3, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="W₀:").grid(row=r, column=4, sticky=tk.E, **_PAD)
        self._var_W0 = tk.DoubleVar(value=0.0)
        ttk.Entry(frame, textvariable=self._var_W0, width=8).grid(
            row=r, column=5, sticky=tk.W, **_PAD
        )

        r = 1
        ttk.Label(frame, text="Δν min:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        self._var_min_dnu = tk.DoubleVar(value=10.0)
        ttk.Entry(frame, textvariable=self._var_min_dnu, width=8).grid(
            row=r, column=1, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="PE max:").grid(row=r, column=2, sticky=tk.E, **_PAD)
        self._var_max_PE = tk.DoubleVar(value=100.0)
        ttk.Entry(frame, textvariable=self._var_max_PE, width=8).grid(
            row=r, column=3, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="Air gap [mm]:").grid(row=r, column=4, sticky=tk.E, **_PAD)
        self._var_air_gap = tk.DoubleVar(value=1.0)
        self._ent_air_gap = ttk.Entry(frame, textvariable=self._var_air_gap, width=8)
        self._ent_air_gap.grid(row=r, column=5, sticky=tk.W, **_PAD)
        self._ent_air_gap.configure(state=tk.DISABLED)
        self._btn_run = ttk.Button(
            frame, text=" ▶ Run Stage A ", command=self._on_run
        )
        self._btn_run.grid(row=r, column=6, sticky=tk.E, padx=6, pady=3)

    def _build_stage_b_frame(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text=" Stage B: Thick-Lens Optimisation (0 = auto) ")
        frame.grid(row=0, column=1, sticky="nsew", **_PAD)

        r = 0
        ttk.Label(frame, text="Field [°]:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        self._var_half_field = tk.DoubleVar(value=1.0)
        ttk.Entry(frame, textvariable=self._var_half_field, width=6).grid(
            row=r, column=1, sticky=tk.W, **_PAD
        )
        ttk.Label(frame, text="Air gap [mm]:").grid(row=r, column=2, sticky=tk.E, **_PAD)
        rng_gap = ttk.Frame(frame)
        rng_gap.grid(row=r, column=3, sticky=tk.W, **_PAD)
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
        ttk.Label(frame, text="t₁ [mm]:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        rng1 = ttk.Frame(frame)
        rng1.grid(row=r, column=1, sticky=tk.W, **_PAD)
        self._var_t1_min = tk.DoubleVar(value=0.0)
        ttk.Entry(rng1, textvariable=self._var_t1_min, width=5).pack(side=tk.LEFT)
        ttk.Label(rng1, text="–").pack(side=tk.LEFT, padx=2)
        self._var_t1_max = tk.DoubleVar(value=0.0)
        ttk.Entry(rng1, textvariable=self._var_t1_max, width=5).pack(side=tk.LEFT)

        ttk.Label(frame, text="t₂ [mm]:").grid(row=r, column=2, sticky=tk.E, **_PAD)
        rng2 = ttk.Frame(frame)
        rng2.grid(row=r, column=3, sticky=tk.W, **_PAD)
        self._var_t2_min = tk.DoubleVar(value=0.0)
        ttk.Entry(rng2, textvariable=self._var_t2_min, width=5).pack(side=tk.LEFT)
        ttk.Label(rng2, text="–").pack(side=tk.LEFT, padx=2)
        self._var_t2_max = tk.DoubleVar(value=0.0)
        ttk.Entry(rng2, textvariable=self._var_t2_max, width=5).pack(side=tk.LEFT)

        self._btn_stage_b = ttk.Button(
            frame, text=" ⚡ Optimize Selected ",
            command=self._on_stage_b, state=tk.DISABLED,
        )
        self._btn_stage_b.grid(row=r, column=4, sticky=tk.E, padx=6, pady=3)

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

        ttk.Separator(bar, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=8, pady=2
        )

        self._progress = ttk.Progressbar(bar, length=200, mode="determinate")
        self._progress.pack(side=tk.LEFT, **_PAD)

        self._var_status = tk.StringVar(value="Ready")
        ttk.Label(bar, textvariable=self._var_status).pack(
            side=tk.LEFT, fill=tk.X, expand=True, **_PAD
        )

    def _on_type_changed(self) -> None:
        """Enable/disable spaced-only fields based on system type."""
        is_spaced = self._var_type.get() == "spaced"
        state = tk.NORMAL if is_spaced else tk.DISABLED
        self._ent_air_gap.configure(state=state)
        self._ent_gap_min.configure(state=state)
        self._ent_gap_max.configure(state=state)

    # ------------------------------------------------------------------
    # ③ Results table
    # ------------------------------------------------------------------

    _COLUMNS = [
        ("#", "idx", 40, tk.CENTER),
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

        # Left: Prescription table
        left = ttk.LabelFrame(frame, text="Prescription")
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=4, pady=2)
        self._detail_rx = tk.Text(
            left, height=10, width=45, font=("Consolas", 9), state=tk.DISABLED
        )
        self._detail_rx.pack(fill=tk.BOTH, expand=True, padx=2, pady=2)

        # Middle: Aberrations
        mid = ttk.LabelFrame(frame, text="Aberrations")
        mid.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=4, pady=2)
        self._detail_ab = tk.Text(
            mid, height=10, width=35, font=("Consolas", 9), state=tk.DISABLED
        )
        self._detail_ab.pack(fill=tk.BOTH, expand=True, padx=2, pady=2)

        # Right: Thickness / Radii comparison
        right = ttk.LabelFrame(frame, text="Thin → Thick Comparison")
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=4, pady=2)
        self._detail_cmp = tk.Text(
            right, height=10, width=40, font=("Consolas", 9), state=tk.DISABLED
        )
        self._detail_cmp.pack(fill=tk.BOTH, expand=True, padx=2, pady=2)

    def _build_drawing_panel(self, parent: ttk.PanedWindow) -> None:
        self._draw_frame = ttk.LabelFrame(self, text=" ⑤ 2D Optical Layout ")
        parent.add(self._draw_frame, weight=1)
        self._mpl_fig = None
        self._mpl_canvas: Optional[FigureCanvasTkAgg] = None

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
    # Load config
    # ===================================================================

    def _on_load_config(self) -> None:
        path = filedialog.askopenfilename(
            title="Select config JSON",
            filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")],
        )
        if not path:
            return
        try:
            d = json.loads(Path(path).read_text(encoding="utf-8"))
        except Exception as e:
            messagebox.showerror("Config Error", str(e))
            return

        # Fill fields
        self._var_fprime.set(d.get("fprime", 200.0))
        self._var_D.set(d.get("D", 50.0))
        self._var_lam0.set(d.get("lam0", 0.58756))
        self._var_lam1.set(d.get("lam1", 0.48613))
        self._var_lam2.set(d.get("lam2", 0.65627))
        self._var_C0.set(d.get("C0", 0.0))
        self._var_P0.set(d.get("P0", 0.0))
        self._var_W0.set(d.get("W0", 0.0))
        self._var_min_dnu.set(d.get("min_delta_nu", 10.0))
        self._var_max_PE.set(d.get("max_PE", 100.0))
        self._var_N.set(d.get("N", 20))
        self._var_type.set(d.get("system_type", "cemented"))
        self._var_air_gap.set(d.get("air_gap", 1.0))
        self._var_half_field.set(d.get("half_field_angle", 1.0))
        self._var_t1_min.set(d.get("t1_min", 0.0))
        self._var_t1_max.set(d.get("t1_max", 0.0))
        self._var_t2_min.set(d.get("t2_min", 0.0))
        self._var_t2_max.set(d.get("t2_max", 0.0))
        self._var_gap_min.set(d.get("gap_min", 0.0))
        self._var_gap_max.set(d.get("gap_max", 0.0))
        self._on_type_changed()

        # Catalog paths (resolve relative to config dir)
        config_dir = Path(path).resolve().parent
        self._catalog_paths.clear()
        for cp in d.get("catalogs", []):
            pp = Path(cp)
            if not pp.is_absolute():
                pp = config_dir / pp
            self._catalog_paths.append(str(pp.resolve()))
        self._refresh_catalog_display()

        self._var_status.set(f"Loaded config: {Path(path).name}")

    # ===================================================================
    # Build Inputs from GUI fields
    # ===================================================================

    def _make_inputs(self) -> Inputs:
        return Inputs(
            lam0=self._var_lam0.get(),
            lam1=self._var_lam1.get(),
            lam2=self._var_lam2.get(),
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
            t1_min=self._var_t1_min.get(),
            t1_max=self._var_t1_max.get(),
            t2_min=self._var_t2_min.get(),
            t2_max=self._var_t2_max.get(),
            gap_min=self._var_gap_min.get(),
            gap_max=self._var_gap_max.get(),
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

        self._running = True
        self._btn_run.configure(state=tk.DISABLED)
        self._btn_export_json.configure(state=tk.DISABLED)
        self._btn_export_csv.configure(state=tk.DISABLED)
        self._results.clear()
        self._clear_tree()
        self._clear_details()

        thread = threading.Thread(target=self._run_pipeline, daemon=True)
        thread.start()

    def _run_pipeline(self) -> None:
        """Execute full pipeline in a background thread."""
        try:
            inputs = self._make_inputs()

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
            self._btn_stage_b.configure(state=tk.NORMAL)

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

        self._running = True
        self._btn_run.configure(state=tk.DISABLED)
        self._btn_stage_b.configure(state=tk.DISABLED)

        self._stage_b_selected = selected_rows
        thread = threading.Thread(target=self._run_stage_b_bg, daemon=True)
        thread.start()

    def _run_stage_b_bg(self) -> None:
        """Execute Stage B optimisation in a background thread."""
        try:
            inputs = self._make_inputs()

            selected_results = [r.result for r in self._stage_b_selected]
            n_sel = len(selected_results)

            def _on_progress(current: int, total: int, pr: PipelineResult) -> None:
                idx = len(self._results)
                row = ResultRow(idx, pr, stage="B")
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

        return (
            row.idx + 1,
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
        """Re-populate tree with current sort order."""
        self._clear_tree()
        for row in self._results:
            values = self._row_values(row)
            tag = "ok" if row.metrics.success else "fail"
            self._tree.insert("", tk.END, iid=str(row.idx), values=values, tags=(tag,))
        self._tree.tag_configure("fail", foreground="#999999")

    # ===================================================================
    # Sorting
    # ===================================================================

    _SORT_KEY_MAP = {
        "idx": lambda r: r.idx,
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

    def _fill_detail(self, row: ResultRow) -> None:
        self._fill_prescription(row)
        self._fill_aberrations(row)
        self._fill_comparison(row)
        self._fill_drawing(row)

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
        lines.append("")
        lines.append("  Spot Diagram")
        lines.append("  ───────────────────────")
        lines.append(f"  RMS spot: {_fmt(m.rms_spot_radius, '.2f')} µm")
        lines.append(f"  GEO spot: {_fmt(m.geo_spot_radius, '.2f')} µm")
        lines.append("")
        lines.append("  First-order")
        lines.append("  ───────────────────────")
        lines.append(f"  EFL (optiland): {_fmt(m.efl, '.3f')} mm")
        lines.append(f"  FNO: {_fmt(m.fno, '.4f')}")
        lines.append(f"  BFD: {_fmt(m.bfd, '.3f')} mm")
        rx = row.rx
        if rx is not None and rx.actual_efl is not None:
            dev_pct = (rx.efl_deviation or 0.0) * 100.0
            lines.append("")
            lines.append("  Thick-lens EFL (ABCD)")
            lines.append("  ───────────────────────")
            lines.append(f"  EFL (ABCD):     {rx.actual_efl:.3f} mm")
            lines.append(f"  EFL (optiland): {_fmt(m.efl, '.3f')} mm")
            lines.append(f"  Deviation:      {dev_pct:+.3f} %")

        w.insert(tk.END, "\n".join(lines))
        w.configure(state=tk.DISABLED)

    def _fill_comparison(self, row: ResultRow) -> None:
        """Fill the thin → thick comparison text widget."""
        w = self._detail_cmp
        w.configure(state=tk.NORMAL)
        w.delete("1.0", tk.END)

        cand = row.cand
        rx = row.rx

        lines = []
        lines.append("  Thin-lens Radii")
        lines.append("  ─────────────────────────")
        for i, r in enumerate(cand.radii):
            lines.append(f"  R{i + 1} = {r:12.4f} mm")

        if rx is not None:
            e1, e2 = rx.elements
            lines.append("")
            lines.append("  Thick-lens Radii (unchanged)")
            lines.append("  ─────────────────────────")
            if rx.system_type == "cemented":
                lines.append(f"  R1 = {e1.R_front:12.4f} mm")
                lines.append(f"  R2 = {e1.R_back:12.4f} mm")
                lines.append(f"  R3 = {e2.R_back:12.4f} mm")
            else:
                lines.append(f"  R1 = {e1.R_front:12.4f} mm")
                lines.append(f"  R2 = {e1.R_back:12.4f} mm")
                lines.append(f"  R3 = {e2.R_front:12.4f} mm")
                lines.append(f"  R4 = {e2.R_back:12.4f} mm")

            lines.append("")
            lines.append("  Thickness")
            lines.append("  ─────────────────────────")
            lines.append(f"  Elem 1: tc={e1.t_center:.3f}  te={e1.t_edge:.3f}")
            lines.append(f"  Elem 2: tc={e2.t_center:.3f}  te={e2.t_edge:.3f}")
            if rx.air_gap is not None:
                lines.append(f"  Air gap: {rx.air_gap:.3f} mm")
            lines.append(f"  Outside dia φ: {rx.D:.1f} mm")
        else:
            lines.append("\n  (thickening failed)")

        # Thermal analysis
        th = cand.thermal
        lines.append("")
        lines.append("  Thermal Analysis")
        lines.append("  ─────────────────────────")
        if th is None or not th.thermal_data_available:
            lines.append("  (no TD/ED data in catalog)")
        else:
            lines.append(f"  dn/dT₁:  {_fmt(th.dn_dT_1, '.3e')} /K  ({cand.glass1})")
            lines.append(f"  dn/dT₂:  {_fmt(th.dn_dT_2, '.3e')} /K  ({cand.glass2})")
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


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    app = AutoAchromatGUI()
    app.mainloop()


if __name__ == "__main__":
    main()

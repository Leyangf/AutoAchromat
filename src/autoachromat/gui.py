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
import time
import tkinter as tk
from dataclasses import asdict
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Optional

# ---------------------------------------------------------------------------
# Project imports
# ---------------------------------------------------------------------------
from .glass_reader import load_catalog
from .models import Candidate, Inputs, ThickPrescription
from .cemented import run_cemented
from .spaced import run_spaced
from .thickening import thicken
from .optiland_bridge.builder import build_optic_from_prescription
from .optiland_bridge.evaluator import evaluate, OpticMetrics


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
    """Bundles Candidate + ThickPrescription + OpticMetrics for one row."""

    __slots__ = ("cand", "rx", "metrics", "idx")

    def __init__(
        self,
        idx: int,
        cand: Candidate,
        rx: Optional[ThickPrescription],
        metrics: OpticMetrics,
    ):
        self.idx = idx
        self.cand = cand
        self.rx = rx
        self.metrics = metrics


# ---------------------------------------------------------------------------
# Main application window
# ---------------------------------------------------------------------------

_PAD = dict(padx=6, pady=3)


class AutoAchromatGUI(tk.Tk):
    """Main GUI window for AutoAchromat."""

    def __init__(self) -> None:
        super().__init__()
        self.title("AutoAchromat – Doublet Designer")
        self.geometry("1280x820")
        self.minsize(960, 640)

        # State
        self._catalog_paths: list[str] = []
        self._results: list[ResultRow] = []
        self._running = False
        self._sort_col: str = "rms"
        self._sort_reverse: bool = False

        self._build_ui()

    # ===================================================================
    # UI construction
    # ===================================================================

    def _build_ui(self) -> None:
        # Top frame: Input params + Controls
        top = ttk.Frame(self)
        top.pack(fill=tk.X, padx=8, pady=4)

        self._build_input_frame(top)
        self._build_control_frame(top)

        # PanedWindow: Results table (top) + Detail panel (bottom)
        pane = ttk.PanedWindow(self, orient=tk.VERTICAL)
        pane.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        self._build_result_table(pane)
        self._build_detail_panel(pane)

    # ------------------------------------------------------------------
    # ① Input parameters
    # ------------------------------------------------------------------

    def _build_input_frame(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text=" ① Input Parameters ")
        frame.pack(fill=tk.X, **_PAD)

        # --- row 0: Optical ---
        r = 0
        ttk.Label(frame, text="f' [mm]:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        self._var_fprime = tk.DoubleVar(value=200.0)
        ttk.Entry(frame, textvariable=self._var_fprime, width=10).grid(
            row=r, column=1, **_PAD
        )

        ttk.Label(frame, text="D [mm]:").grid(row=r, column=2, sticky=tk.E, **_PAD)
        self._var_D = tk.DoubleVar(value=50.0)
        ttk.Entry(frame, textvariable=self._var_D, width=10).grid(
            row=r, column=3, **_PAD
        )

        # Separator
        ttk.Separator(frame, orient=tk.VERTICAL).grid(
            row=0, column=4, rowspan=3, sticky="ns", padx=10
        )

        # --- Wavelengths ---
        ttk.Label(frame, text="λ₀ [µm]:").grid(row=r, column=5, sticky=tk.E, **_PAD)
        self._var_lam0 = tk.DoubleVar(value=0.58756)
        ttk.Entry(frame, textvariable=self._var_lam0, width=10).grid(
            row=r, column=6, **_PAD
        )

        ttk.Label(frame, text="λ₁ [µm]:").grid(row=r, column=7, sticky=tk.E, **_PAD)
        self._var_lam1 = tk.DoubleVar(value=0.48613)
        ttk.Entry(frame, textvariable=self._var_lam1, width=10).grid(
            row=r, column=8, **_PAD
        )

        ttk.Label(frame, text="λ₂ [µm]:").grid(row=r, column=9, sticky=tk.E, **_PAD)
        self._var_lam2 = tk.DoubleVar(value=0.65627)
        ttk.Entry(frame, textvariable=self._var_lam2, width=10).grid(
            row=r, column=10, **_PAD
        )

        # --- row 1: Targets ---
        r = 1
        ttk.Label(frame, text="C₀:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        self._var_C0 = tk.DoubleVar(value=0.0)
        ttk.Entry(frame, textvariable=self._var_C0, width=10).grid(
            row=r, column=1, **_PAD
        )

        ttk.Label(frame, text="P₀:").grid(row=r, column=2, sticky=tk.E, **_PAD)
        self._var_P0 = tk.DoubleVar(value=1.0)
        ttk.Entry(frame, textvariable=self._var_P0, width=10).grid(
            row=r, column=3, **_PAD
        )

        ttk.Label(frame, text="W₀:").grid(row=r, column=5, sticky=tk.E, **_PAD)
        self._var_W0 = tk.DoubleVar(value=0.0)
        ttk.Entry(frame, textvariable=self._var_W0, width=10).grid(
            row=r, column=6, **_PAD
        )

        ttk.Label(frame, text="Δν min:").grid(row=r, column=7, sticky=tk.E, **_PAD)
        self._var_min_dnu = tk.DoubleVar(value=10.0)
        ttk.Entry(frame, textvariable=self._var_min_dnu, width=10).grid(
            row=r, column=8, **_PAD
        )

        ttk.Label(frame, text="PE max:").grid(row=r, column=9, sticky=tk.E, **_PAD)
        self._var_max_PE = tk.DoubleVar(value=100.0)
        ttk.Entry(frame, textvariable=self._var_max_PE, width=10).grid(
            row=r, column=10, **_PAD
        )

        # --- row 2: Type + Top-N ---
        r = 2
        ttk.Label(frame, text="Type:").grid(row=r, column=0, sticky=tk.E, **_PAD)
        self._var_type = tk.StringVar(value="cemented")
        type_frame = ttk.Frame(frame)
        type_frame.grid(row=r, column=1, columnspan=3, sticky=tk.W, **_PAD)
        ttk.Radiobutton(
            type_frame, text="Cemented", variable=self._var_type, value="cemented"
        ).pack(side=tk.LEFT, padx=4)
        ttk.Radiobutton(
            type_frame, text="Spaced", variable=self._var_type, value="spaced"
        ).pack(side=tk.LEFT, padx=4)

        ttk.Label(frame, text="Top-N:").grid(row=r, column=5, sticky=tk.E, **_PAD)
        self._var_N = tk.IntVar(value=20)
        ttk.Entry(frame, textvariable=self._var_N, width=10).grid(
            row=r, column=6, **_PAD
        )

        # --- Catalog list ---
        ttk.Label(frame, text="Catalogs:").grid(row=r, column=7, sticky=tk.E, **_PAD)
        cat_frame = ttk.Frame(frame)
        cat_frame.grid(row=r, column=8, columnspan=3, sticky=tk.W, **_PAD)

        self._var_cat_display = tk.StringVar(value="(none)")
        ttk.Label(
            cat_frame, textvariable=self._var_cat_display, width=30, relief="sunken"
        ).pack(side=tk.LEFT)
        ttk.Button(cat_frame, text="+", width=3, command=self._add_catalog).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Button(cat_frame, text="−", width=3, command=self._clear_catalogs).pack(
            side=tk.LEFT, padx=2
        )

    # ------------------------------------------------------------------
    # ② Control bar
    # ------------------------------------------------------------------

    def _build_control_frame(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text=" ② Control ")
        frame.pack(fill=tk.X, **_PAD)

        self._btn_run = ttk.Button(
            frame, text=" ▶  Run Pipeline ", command=self._on_run
        )
        self._btn_run.pack(side=tk.LEFT, **_PAD)

        self._btn_load = ttk.Button(
            frame, text="Load Config", command=self._on_load_config
        )
        self._btn_load.pack(side=tk.LEFT, **_PAD)

        ttk.Separator(frame, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=8, pady=2
        )

        self._btn_export_json = ttk.Button(
            frame, text="Export JSON", command=self._on_export_json, state=tk.DISABLED
        )
        self._btn_export_json.pack(side=tk.LEFT, **_PAD)

        self._btn_export_csv = ttk.Button(
            frame, text="Export CSV", command=self._on_export_csv, state=tk.DISABLED
        )
        self._btn_export_csv.pack(side=tk.LEFT, **_PAD)

        ttk.Separator(frame, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=8, pady=2
        )

        self._progress = ttk.Progressbar(frame, length=200, mode="determinate")
        self._progress.pack(side=tk.LEFT, **_PAD)

        self._var_status = tk.StringVar(value="Ready")
        ttk.Label(frame, textvariable=self._var_status, width=60).pack(
            side=tk.LEFT, **_PAD
        )

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
        ("ms", "time", 55, tk.E),
    ]

    def _build_result_table(self, parent: ttk.PanedWindow) -> None:
        frame = ttk.LabelFrame(self, text=" ③ Results ")
        parent.add(frame, weight=3)

        col_ids = [c[1] for c in self._COLUMNS]
        self._tree = ttk.Treeview(
            frame, columns=col_ids, show="headings", selectmode="browse"
        )

        for display, cid, width, anchor in self._COLUMNS:
            self._tree.heading(
                cid, text=display, command=lambda c=cid: self._on_sort(c)
            )
            self._tree.column(cid, width=width, anchor=anchor, minwidth=40)

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
        parent.add(frame, weight=2)

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
        self._var_P0.set(d.get("P0", 1.0))
        self._var_W0.set(d.get("W0", 0.0))
        self._var_min_dnu.set(d.get("min_delta_nu", 10.0))
        self._var_max_PE.set(d.get("max_PE", 100.0))
        self._var_N.set(d.get("N", 20))
        self._var_type.set(d.get("system_type", "cemented"))

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
            system_type=self._var_type.get(),
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

            # Step 1: Load catalogs
            self._set_status("Loading catalogs...")
            glasses = load_catalog(self._catalog_paths)
            n_glasses = (
                sum(len(v) for v in glasses.values())
                if isinstance(glasses, dict)
                else len(glasses)
            )
            self._set_status(
                f"Loaded {n_glasses} glasses from {len(self._catalog_paths)} catalog(s)"
            )

            # Step 2: Thin-lens synthesis
            self._set_status("Running thin-lens synthesis...")
            t0 = time.perf_counter()
            if inputs.system_type == "cemented":
                cands = run_cemented(inputs, glasses)
            elif inputs.system_type == "spaced":
                cands = run_spaced(inputs, glasses)
            else:
                raise ValueError(f"Unknown system_type={inputs.system_type}")
            synth_ms = (time.perf_counter() - t0) * 1000.0
            self._set_status(f"Synthesis: {len(cands)} candidates in {synth_ms:.0f} ms")

            # Step 3-5: Thicken → Build → Evaluate  (top-N)
            top_n = inputs.N
            todo = cands[:top_n]
            self._set_progress(0, len(todo))

            for i, cand in enumerate(todo):
                self._set_status(
                    f"Processing {i + 1}/{len(todo)}: {cand.glass1} + {cand.glass2} ..."
                )

                # Thicken
                t_start = time.perf_counter()
                rx = thicken(cand, inputs)

                if rx is None:
                    m = OpticMetrics(
                        glass1=cand.glass1,
                        glass2=cand.glass2,
                        catalog1=cand.catalog1,
                        catalog2=cand.catalog2,
                        system_type=cand.system_type,
                        PE=cand.PE,
                        success=False,
                        error_msg="thickening failed",
                        build_time_ms=(time.perf_counter() - t_start) * 1000.0,
                    )
                    row = ResultRow(i, cand, None, m)
                    self._results.append(row)
                    self._insert_tree_row(row)
                    self._set_progress(i + 1, len(todo))
                    continue

                # Build
                op = build_optic_from_prescription(rx)
                build_ms = (time.perf_counter() - t_start) * 1000.0

                if op is None:
                    m = OpticMetrics(
                        glass1=cand.glass1,
                        glass2=cand.glass2,
                        catalog1=cand.catalog1,
                        catalog2=cand.catalog2,
                        system_type=cand.system_type,
                        PE=cand.PE,
                        success=False,
                        error_msg="build failed",
                        build_time_ms=build_ms,
                    )
                    row = ResultRow(i, cand, rx, m)
                    self._results.append(row)
                    self._insert_tree_row(row)
                    self._set_progress(i + 1, len(todo))
                    continue

                # Evaluate
                m = evaluate(op, cand, inputs)
                m.build_time_ms = build_ms

                row = ResultRow(i, cand, rx, m)
                self._results.append(row)
                self._insert_tree_row(row)
                self._set_progress(i + 1, len(todo))

            n_ok = sum(1 for r in self._results if r.metrics.success)
            self._set_status(
                f"Done — {n_ok}/{len(todo)} succeeded  |  "
                f"Synthesis: {len(cands)} candidates  |  {synth_ms:.0f} ms"
            )

        except Exception as e:
            self._set_status(f"ERROR: {e}")
            self.after(0, lambda: messagebox.showerror("Pipeline Error", str(e)))

        finally:
            self.after(0, self._pipeline_done)

    def _pipeline_done(self) -> None:
        self._running = False
        self._btn_run.configure(state=tk.NORMAL)
        if self._results:
            self._btn_export_json.configure(state=tk.NORMAL)
            self._btn_export_csv.configure(state=tk.NORMAL)

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
        total_ms = m.build_time_ms + m.eval_time_ms
        tp = "Cem" if row.cand.system_type == "cemented" else "Sp"

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
            f"{total_ms:.0f}",
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
        "time": lambda r: r.metrics.build_time_ms + r.metrics.eval_time_ms,
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

    def _fill_detail(self, row: ResultRow) -> None:
        self._fill_prescription(row)
        self._fill_aberrations(row)
        self._fill_comparison(row)

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
            lines.append(f"   4  {e2.R_back:10.3f}  {rx.back_focus_guess:8.1f}  Air")
        else:
            lines.append(
                f" STO  {e1.R_front:10.3f}  {e1.t_center:8.3f}  n={e1.nd:.4f} ν={e1.vd:.1f}"
            )
            lines.append(f"   3  {e1.R_back:10.3f}  {rx.air_gap or 0:8.3f}  Air")
            lines.append(
                f"   4  {e2.R_front:10.3f}  {e2.t_center:8.3f}  n={e2.nd:.4f} ν={e2.vd:.1f}"
            )
            lines.append(f"   5  {e2.R_back:10.3f}  {rx.back_focus_guess:8.1f}  Air")

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
        lines.append(f"  EFL: {_fmt(m.efl, '.3f')} mm")
        lines.append(f"  FNO: {_fmt(m.fno, '.4f')}")

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
            lines.append("  Thick-lens Radii (corrected)")
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

        w.insert(tk.END, "\n".join(lines))
        w.configure(state=tk.DISABLED)

    # ===================================================================
    # Export
    # ===================================================================

    def _export_data(self) -> list[dict]:
        """Build JSON-serialisable list from results."""
        out = []
        for row in self._results:
            d = asdict(row.metrics)
            d["thin_radii"] = row.cand.radii
            if row.rx is not None:
                elems = row.rx.elements
                d["thick_radii"] = [
                    [elems[0].R_front, elems[0].R_back],
                    [elems[1].R_front, elems[1].R_back],
                ]
                d["t_center"] = [elems[0].t_center, elems[1].t_center]
                d["t_edge"] = [elems[0].t_edge, elems[1].t_edge]
                d["air_gap"] = row.rx.air_gap
                d["back_focus_guess"] = row.rx.back_focus_guess
            else:
                d["thick_radii"] = None
                d["t_center"] = None
                d["t_edge"] = None
                d["air_gap"] = None
                d["back_focus_guess"] = None
            out.append(d)
        return out

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

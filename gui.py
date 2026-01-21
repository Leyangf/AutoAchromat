import sys
from argparse import Namespace
from pathlib import Path

import matplotlib.pyplot as plt
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGroupBox,
    QFormLayout,
    QLineEdit,
    QComboBox,
    QSpinBox,
    QDoubleSpinBox,
    QPushButton,
    QFileDialog,
    QTableWidget,
    QTableWidgetItem,
    QListWidget,
    QLabel,
    QMessageBox,
)

# Add src to path for direct script execution
sys.path.insert(0, str(Path(__file__).parent))

from src import cli, io
from src.synthesis import run, default_cfg
from src.raytrace import OPTILAND_AVAILABLE
from src.raytrace.optiland_interface import compute_realray_for_candidate


class MainWindow(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("AutoAchromat GUI")
        self.resize(1200, 700)

        self.config_path: str | None = None
        self.rows: list[dict] = []  # Store rows for real-ray calculation
        self.cfg = None  # Store config for real-ray calculation

        self._build_ui()
        self._load_defaults()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)

        # Input area
        input_group = QGroupBox("Input Parameters")
        input_layout = QHBoxLayout()

        # Left column: catalogs and system
        left_form = QFormLayout()
        self.agf_list = QListWidget()
        self.agf_list.setMinimumHeight(80)

        agf_buttons = QHBoxLayout()
        self.btn_add_agf = QPushButton("Add AGF")
        self.btn_remove_agf = QPushButton("Remove")
        agf_buttons.addWidget(self.btn_add_agf)
        agf_buttons.addWidget(self.btn_remove_agf)

        left_form.addRow("AGF catalogs", self.agf_list)
        left_form.addRow("", agf_buttons)

        self.system_type = QComboBox()
        self.system_type.addItems(["cemented", "air_spaced"])
        self.N_keep = QSpinBox()
        self.N_keep.setRange(1, 10000)
        self.display_rows = QSpinBox()
        self.display_rows.setRange(1, 10000)
        self.out_dir = QLineEdit()
        self.btn_out_dir = QPushButton("Choose output")

        out_dir_layout = QHBoxLayout()
        out_dir_layout.addWidget(self.out_dir)
        out_dir_layout.addWidget(self.btn_out_dir)

        left_form.addRow("System type", self.system_type)
        left_form.addRow("N_keep", self.N_keep)
        left_form.addRow("Display rows", self.display_rows)
        left_form.addRow("Output dir", out_dir_layout)

        # Middle column: optical parameters
        mid_form = QFormLayout()
        self.focal = self._make_double_box(0.0, 1e6, 3)
        self.aperture = self._make_double_box(0.0, 1e6, 3)
        self.P0 = self._make_double_box(-1e6, 1e6, 6)
        self.W0 = self._make_double_box(-1e6, 1e6, 6)
        self.C0 = self._make_double_box(-1e6, 1e6, 6)

        mid_form.addRow("Focal length f [mm]", self.focal)
        mid_form.addRow("Aperture D [mm]", self.aperture)
        mid_form.addRow("P0", self.P0)
        mid_form.addRow("W0", self.W0)
        mid_form.addRow("C0", self.C0)

        # Right column: wavelengths and filters
        right_form = QFormLayout()
        self.lam0 = self._make_double_box(0.1, 2.0, 6)
        self.lam1 = self._make_double_box(0.1, 2.0, 6)
        self.lam2 = self._make_double_box(0.1, 2.0, 6)
        self.min_delta_nu = self._make_double_box(0.0, 1e6, 4)
        self.max_PE = QLineEdit()
        self.d_air = self._make_double_box(0.0, 1e6, 3)

        right_form.addRow("lam0 [um]", self.lam0)
        right_form.addRow("lam1 [um]", self.lam1)
        right_form.addRow("lam2 [um]", self.lam2)
        right_form.addRow("min_delta_nu", self.min_delta_nu)
        right_form.addRow("max_PE", self.max_PE)
        right_form.addRow("d_air [mm]", self.d_air)

        # Real-ray aberration parameters
        self.crown_thickness = self._make_double_box(0.1, 100.0, 3)
        self.flint_thickness = self._make_double_box(0.1, 100.0, 3)
        self.field_deg = self._make_double_box(0.0, 45.0, 3)
        right_form.addRow("Crown thickness [mm]", self.crown_thickness)
        right_form.addRow("Flint thickness [mm]", self.flint_thickness)
        right_form.addRow("Field angle [deg]", self.field_deg)

        # Action buttons
        action_layout = QVBoxLayout()
        self.btn_load_config = QPushButton("Load config")
        self.btn_run = QPushButton("Run")
        self.btn_realray = QPushButton("Compute All Real-Ray")
        self.btn_realray.setEnabled(False)  # Disabled until results are available
        self.btn_realray.setToolTip(
            "Compute real-ray aberrations for all candidates in the table"
        )
        self.btn_draw_2d = QPushButton("Show 2D Drawing")
        self.btn_draw_2d.setEnabled(False)  # Disabled until results are available
        self.btn_draw_2d.setToolTip("Draw 2D layout of the selected optical system")
        action_layout.addWidget(self.btn_load_config)
        action_layout.addWidget(self.btn_run)
        action_layout.addWidget(self.btn_realray)
        action_layout.addWidget(self.btn_draw_2d)
        action_layout.addStretch()

        input_layout.addLayout(left_form)
        input_layout.addLayout(mid_form)
        input_layout.addLayout(right_form)
        input_layout.addLayout(action_layout)
        input_group.setLayout(input_layout)

        # Table area
        table_group = QGroupBox("Result of calculation")
        table_layout = QVBoxLayout()
        self.table = QTableWidget(0, 0)
        self.table.setAlternatingRowColors(True)
        table_layout.addWidget(self.table)
        table_group.setLayout(table_layout)

        # Status
        self.status = QLabel("Ready")
        self.status.setAlignment(Qt.AlignmentFlag.AlignLeft)

        main_layout.addWidget(input_group)
        main_layout.addWidget(table_group)
        main_layout.addWidget(self.status)

        # Wiring
        self.btn_add_agf.clicked.connect(self._on_add_agf)
        self.btn_remove_agf.clicked.connect(self._on_remove_agf)
        self.btn_out_dir.clicked.connect(self._on_choose_out_dir)
        self.btn_load_config.clicked.connect(self._on_load_config)
        self.btn_run.clicked.connect(self._on_run)
        self.btn_realray.clicked.connect(self._on_compute_all_realray)
        self.btn_draw_2d.clicked.connect(self._on_draw_2d)
        self.system_type.currentTextChanged.connect(self._on_system_changed)

    def _make_double_box(
        self, min_val: float, max_val: float, decimals: int
    ) -> QDoubleSpinBox:
        box = QDoubleSpinBox()
        box.setRange(min_val, max_val)
        box.setDecimals(decimals)
        box.setSingleStep(0.1)
        return box

    def _load_defaults(self) -> None:
        cfg = default_cfg()
        self._set_agf_paths(cfg.agf_paths)
        self.system_type.setCurrentText(cfg.system_type)
        self.N_keep.setValue(cfg.N_keep)
        self.display_rows.setValue(cfg.N_keep)
        self.out_dir.setText(cfg.out_dir)

        self.focal.setValue(cfg.f)
        self.aperture.setValue(cfg.D)
        self.P0.setValue(cfg.P0)
        self.W0.setValue(cfg.W0)
        self.C0.setValue(cfg.C0)

        self.lam0.setValue(cfg.lam0)
        self.lam1.setValue(cfg.lam1)
        self.lam2.setValue(cfg.lam2)
        self.min_delta_nu.setValue(cfg.min_delta_nu)
        self.max_PE.setText(str(cfg.max_PE))
        self.d_air.setValue(cfg.d_air)
        self.crown_thickness.setValue(cfg.crown_lens_thickness_mm)
        self.flint_thickness.setValue(cfg.flint_lens_thickness_mm)
        self.field_deg.setValue(cfg.field_for_aberration_deg)

        self._on_system_changed(cfg.system_type)

    def _set_agf_paths(self, paths: list[str]) -> None:
        self.agf_list.clear()
        for p in paths:
            self.agf_list.addItem(p)

    def _get_agf_paths(self) -> list[str]:
        return [self.agf_list.item(i).text() for i in range(self.agf_list.count())]

    def _on_add_agf(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select AGF catalogs",
            str(Path.cwd()),
            "AGF files (*.AGF *.agf);;All files (*)",
        )
        if not paths:
            return
        existing = set(self._get_agf_paths())
        for p in paths:
            if p not in existing:
                self.agf_list.addItem(p)

    def _on_remove_agf(self) -> None:
        for item in self.agf_list.selectedItems():
            row = self.agf_list.row(item)
            self.agf_list.takeItem(row)

    def _on_choose_out_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select output directory", str(Path.cwd())
        )
        if path:
            self.out_dir.setText(path)

    def _on_system_changed(self, system_type: str) -> None:
        is_air = system_type == "air_spaced"
        self.d_air.setEnabled(is_air)

    def _on_load_config(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select config JSON",
            str(Path.cwd()),
            "JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        try:
            config = cli.load_config_file(path)
        except Exception as exc:
            QMessageBox.critical(self, "Config error", str(exc))
            return

        self.config_path = path
        self._apply_config_to_ui(config)
        self.status.setText(f"Loaded config: {path}")

    def _apply_config_to_ui(self, config: dict) -> None:
        if "agf_paths" in config and config["agf_paths"]:
            self._set_agf_paths(list(config["agf_paths"]))
        if "system_type" in config and config["system_type"]:
            self.system_type.setCurrentText(str(config["system_type"]))
        if "N_keep" in config and config["N_keep"] is not None:
            self.N_keep.setValue(int(config["N_keep"]))
        if "display" in config and config["display"] is not None:
            self.display_rows.setValue(int(config["display"]))
        if "f" in config and config["f"] is not None:
            self.focal.setValue(float(config["f"]))
        if "D" in config and config["D"] is not None:
            self.aperture.setValue(float(config["D"]))
        if "P0" in config and config["P0"] is not None:
            self.P0.setValue(float(config["P0"]))
        if "W0" in config and config["W0"] is not None:
            self.W0.setValue(float(config["W0"]))
        if "C0" in config and config["C0"] is not None:
            self.C0.setValue(float(config["C0"]))
        if "lam0" in config and config["lam0"] is not None:
            self.lam0.setValue(float(config["lam0"]))
        if "lam1" in config and config["lam1"] is not None:
            self.lam1.setValue(float(config["lam1"]))
        if "lam2" in config and config["lam2"] is not None:
            self.lam2.setValue(float(config["lam2"]))
        if "min_delta_nu" in config and config["min_delta_nu"] is not None:
            self.min_delta_nu.setValue(float(config["min_delta_nu"]))
        if "max_PE" in config and config["max_PE"] is not None:
            self.max_PE.setText(str(config["max_PE"]))
        if "d_air" in config and config["d_air"] is not None:
            self.d_air.setValue(float(config["d_air"]))
        if (
            "crown_lens_thickness_mm" in config
            and config["crown_lens_thickness_mm"] is not None
        ):
            self.crown_thickness.setValue(float(config["crown_lens_thickness_mm"]))
        if (
            "flint_lens_thickness_mm" in config
            and config["flint_lens_thickness_mm"] is not None
        ):
            self.flint_thickness.setValue(float(config["flint_lens_thickness_mm"]))
        if (
            "field_for_aberration_deg" in config
            and config["field_for_aberration_deg"] is not None
        ):
            self.field_deg.setValue(float(config["field_for_aberration_deg"]))
        if "out_dir" in config and config["out_dir"]:
            self.out_dir.setText(str(config["out_dir"]))

    def _collect_args(self) -> Namespace:
        max_pe_text = self.max_PE.text().strip()
        if not max_pe_text:
            max_pe_value = None
        else:
            try:
                max_pe_value = float(max_pe_text)
            except ValueError as exc:
                raise ValueError("max_PE must be a number") from exc

        agf_paths = self._get_agf_paths()
        if not agf_paths:
            agf_paths = None

        out_dir = self.out_dir.text().strip() or None

        return Namespace(
            config=self.config_path,
            agf_paths=agf_paths,
            system=self.system_type.currentText(),
            N=self.N_keep.value(),
            f=self.focal.value(),
            D=self.aperture.value(),
            P0=self.P0.value(),
            W0=self.W0.value(),
            C0=self.C0.value(),
            lam0=self.lam0.value(),
            lam1=self.lam1.value(),
            lam2=self.lam2.value(),
            min_delta_nu=self.min_delta_nu.value(),
            max_PE=max_pe_value,
            d_air=self.d_air.value(),
            crown_lens_thickness_mm=self.crown_thickness.value(),
            flint_lens_thickness_mm=self.flint_thickness.value(),
            field_for_aberration_deg=self.field_deg.value(),
            out_dir=out_dir,
            display=self.display_rows.value(),
        )

    def _on_run(self) -> None:
        try:
            args = self._collect_args()
        except ValueError as exc:
            QMessageBox.critical(self, "Input error", str(exc))
            return

        result = cli.build_config(args)
        if result[0] is None:
            QMessageBox.critical(self, "Config error", result[1])
            return

        cfg, display_count = result

        self.status.setText("Running synthesis pipeline...")
        QApplication.processEvents()

        results, stats = run(cfg)
        rows = io.extract_rows(results)

        io.ensure_out_dir(cfg.out_dir)
        io.save_csv(rows, f"{cfg.out_dir}/results_table.csv")
        io.save_report_md(rows, stats, cfg, display_count, f"{cfg.out_dir}/report.md")
        io.save_resolved_config(
            cfg, display_count, f"{cfg.out_dir}/resolved_config.json"
        )

        self._populate_table(rows, display_count)

        # Store for real-ray calculation
        self.rows = rows
        self.cfg = cfg
        self.btn_realray.setEnabled(len(rows) > 0)
        self.btn_draw_2d.setEnabled(len(rows) > 0)

        best_pe = stats.get("best_PE")
        status_msg = f"Done. Candidates: {len(rows)}"
        if best_pe is not None:
            status_msg += f" | best_PE: {io.fmt_float(best_pe, 6)}"
        status_msg += f" | output: {cfg.out_dir}"
        self.status.setText(status_msg)

        if not rows:
            QMessageBox.information(
                self, "No results", "No candidates found matching the criteria."
            )

    def _populate_table(self, rows: list[dict], display_count: int) -> None:
        columns = io.COLUMNS
        headers = [c[1] for c in columns]

        display_rows = rows[:display_count] if display_count else rows

        self.table.setColumnCount(len(columns))
        self.table.setRowCount(len(display_rows))
        self.table.setHorizontalHeaderLabels(headers)

        for r, row in enumerate(display_rows):
            for c, (key, _, _, is_float) in enumerate(columns):
                val = row.get(key)
                if val is None:
                    text = ""
                elif is_float:
                    text = io.fmt_float(val, 6)
                else:
                    text = str(val)
                item = QTableWidgetItem(text)
                if is_float:
                    item.setTextAlignment(
                        Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                    )
                self.table.setItem(r, c, item)

        self.table.resizeColumnsToContents()

    def _on_compute_all_realray(self) -> None:
        """Compute real-ray aberrations for all candidates in the table."""
        # Check Optiland availability
        if not OPTILAND_AVAILABLE:
            QMessageBox.warning(
                self,
                "Optiland not available",
                "Optiland is not installed. Please install it with:\n\n"
                "  pip install optiland\n\n"
                "Real-ray aberration calculation requires Optiland.",
            )
            return

        if not self.rows:
            QMessageBox.information(
                self, "No candidates", "No candidates available. Run synthesis first."
            )
            return

        if self.cfg is None:
            QMessageBox.warning(
                self,
                "No configuration",
                "No configuration loaded. Run synthesis first.",
            )
            return

        # Update config with current GUI values for thickness and field
        self.cfg.crown_lens_thickness_mm = self.crown_thickness.value()
        self.cfg.flint_lens_thickness_mm = self.flint_thickness.value()
        self.cfg.field_for_aberration_deg = self.field_deg.value()

        total_rows = len(self.rows)
        errors = []
        self.btn_realray.setEnabled(False)
        self.btn_draw_2d.setEnabled(False)
        self.btn_run.setEnabled(False)

        try:
            for row_index, row in enumerate(self.rows):
                candidate_rank = row.get("rank", row_index + 1)

                self.status.setText(
                    f"Computing real-ray aberrations: {row_index + 1}/{total_rows} "
                    f"(candidate #{candidate_rank})..."
                )
                QApplication.processEvents()

                try:
                    # Compute real-ray aberrations
                    realray_results = compute_realray_for_candidate(row, self.cfg)

                    # Update the row data
                    self.rows[row_index] = io.update_row_with_realray(
                        row, realray_results
                    )

                    # Update CSV row
                    csv_path = f"{self.cfg.out_dir}/results_table.csv"
                    io.update_csv_row(csv_path, row_index, realray_results)

                    # Update table display
                    self._update_table_row(row_index, realray_results)

                    # Track errors
                    if realray_results.get("_realray_error"):
                        errors.append(
                            f"#{candidate_rank}: {realray_results['_realray_error']}"
                        )

                except Exception as exc:
                    errors.append(f"#{candidate_rank}: {exc}")

            # Show completion message
            if errors:
                error_msg = "\n".join(errors[:10])  # Show first 10 errors
                if len(errors) > 10:
                    error_msg += f"\n... and {len(errors) - 10} more errors"
                QMessageBox.warning(
                    self,
                    "Computation completed with errors",
                    f"Real-ray aberrations computed for {total_rows} candidates.\n\n"
                    f"Errors encountered:\n{error_msg}",
                )
            else:
                QMessageBox.information(
                    self,
                    "Computation complete",
                    f"Real-ray aberrations computed for all {total_rows} candidates.\n\n"
                    f"Results saved to: {self.cfg.out_dir}",
                )

            self.status.setText(
                f"Real-ray aberrations computed for {total_rows} candidates | "
                f"output: {self.cfg.out_dir}"
            )

        except Exception as exc:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to compute real-ray aberrations:\n\n{exc}",
            )
            self.status.setText("Real-ray calculation failed")

        finally:
            self.btn_realray.setEnabled(True)
            self.btn_draw_2d.setEnabled(True)
            self.btn_run.setEnabled(True)

    def _update_table_row(self, row_index: int, realray_results: dict) -> None:
        """Update a single table row with real-ray results."""
        columns = io.COLUMNS

        for c, (key, _, _, is_float) in enumerate(columns):
            if key in realray_results:
                val = realray_results.get(key)
                if val is None:
                    text = ""
                elif is_float:
                    text = io.fmt_float(val, 6)
                else:
                    text = str(val)
                item = QTableWidgetItem(text)
                if is_float:
                    item.setTextAlignment(
                        Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                    )
                self.table.setItem(row_index, c, item)

        self.table.resizeColumnsToContents()

    def _on_draw_2d(self) -> None:
        """Draw 2D layout of the selected optical system."""
        # Check prerequisites
        if not hasattr(self, "rows") or not self.rows:
            QMessageBox.warning(self, "No Results", "No results available to draw.")
            return

        if self.cfg is None:
            QMessageBox.warning(self, "No Config", "No configuration available.")
            return

        # Get selected row
        selected_items = self.table.selectedItems()
        if not selected_items:
            QMessageBox.warning(
                self, "No Selection", "Please select a row from the table first."
            )
            return

        selected_row_index = selected_items[0].row()

        if selected_row_index < 0 or selected_row_index >= len(self.rows):
            QMessageBox.warning(self, "Invalid Selection", "Invalid row selected.")
            return

        # Get the candidate row data
        candidate_row = self.rows[selected_row_index]

        # Import and call the draw function
        try:
            from src.raytrace.optiland_interface import draw_optical_system

            # Update config with current GUI values for thickness
            self.cfg.crown_lens_thickness_mm = self.crown_thickness.value()
            self.cfg.flint_lens_thickness_mm = self.flint_thickness.value()
            self.cfg.field_for_aberration_deg = self.field_deg.value()

            self.status.setText("Drawing 2D layout...")
            QApplication.processEvents()

            fig, error = draw_optical_system(candidate_row, self.cfg, num_rays=5)

            if error:
                QMessageBox.critical(self, "Draw Error", error)
                self.status.setText("Drawing failed")
                return

            # Show the figure
            plt.show()
            self.status.setText("2D layout displayed")

        except ImportError as exc:
            QMessageBox.critical(
                self,
                "Import Error",
                f"Failed to import draw function:\n\n{exc}\n\n"
                "Make sure Optiland is installed: pip install optiland",
            )
            self.status.setText("Drawing failed - Optiland not available")

        except Exception as exc:
            QMessageBox.critical(
                self, "Error", f"Failed to draw optical system:\n\n{exc}"
            )
            self.status.setText("Drawing failed")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())

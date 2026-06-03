"""
Family Selection Widget
=======================
Optional first subtab of the "Analysis in R" tab.

Runs FamilySelection.R to compare Tweedie, Gamma(log), and Negative Binomial
families on the user's data and recommends the best fit based on AIC, BIC,
dispersion, deviance explained, and zero-proportion matching.

If the user skips this step entirely, Tweedie is used as the default for the
BAM analysis (same as the original behaviour).
"""

import glob as _glob
import os
import shutil
import signal
import subprocess
import sys
import threading

import pandas as pd

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QComboBox, QGroupBox, QHBoxLayout, QLabel, QMessageBox,
    QPushButton, QSplitter, QTableWidget, QTableWidgetItem,
    QTextEdit, QVBoxLayout, QWidget,
)

from paths import resource_path, external_data_dir

_APP_R_LIB            = os.path.join(external_data_dir(), "R", "library")
_BUNDLED_RSCRIPT      = os.path.join(
    external_data_dir(), "R", "runtime", "bin",
    "Rscript.exe" if sys.platform == "win32" else "Rscript"
)
_FAMILY_SELECTION_R   = resource_path("R", "scripts", "FamilySelection.R")

# Human-readable names for display
_FAMILY_DISPLAY = {
    "Tweedie":     "Tweedie (compound Poisson-Gamma)",
    "Gamma":       "Gamma (log link)",
    "NegBinomial": "Negative Binomial",
}

# Interpretation text shown after the test
_FAMILY_INTERPRETATION = {
    "Tweedie": (
        "Tweedie (compound Poisson-Gamma) provided the best fit. This is consistent with "
        "movement data that includes genuine zero-activity periods — the data has meaningful "
        "zero mass that Gamma and Negative Binomial cannot capture as well. "
        "The Tweedie assumption is justified for this dataset."
    ),
    "Gamma": (
        "Gamma (log link) provided the best fit. The data has minimal or no zero-activity "
        "periods and the simpler Gamma distribution describes the variance structure better "
        "than Tweedie or Negative Binomial. A small shift was applied to ensure all values "
        "are strictly positive."
    ),
    "NegBinomial": (
        "Negative Binomial provided the best fit. This suggests pixel differences are best "
        "described as overdispersed count data — the discrete, integer nature of pixel counts "
        "is better captured by the Negative Binomial than the continuous Tweedie distribution."
    ),
}


class FamilySelectionWidget(QWidget):
    """
    Optional subtab that runs FamilySelection.R and recommends a distribution
    family for the BAM analysis.  Emits family_changed whenever the active
    family is updated (either from the test result or a manual override).
    """

    # Emits (family_name: str, shift_val: float) on any family change
    family_changed = Signal(str, float)

    _log_line_ready = Signal(str)
    _run_finished   = Signal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._df             = None
        self._output_dir     = None
        self._process        = None
        self._variable_names = ["Genotype", "Drug"]
        self._variable_refs  = {}
        self._ref_condition  = ""
        self._active_family  = "Tweedie"
        self._active_shift   = 0.0

        self._build_ui()
        self._log_line_ready.connect(self._append_log)
        self._run_finished.connect(self._on_run_finished)

    # ── Public API ────────────────────────────────────────────────────────────

    def set_variable_names(self, names: list):
        self._variable_names = list(names)

    def set_references(self, variable_refs: dict, ref_condition: str):
        self._variable_refs = dict(variable_refs)
        self._ref_condition = ref_condition

    def set_output_dir(self, path: str):
        self._output_dir = path
        self._update_run_enabled()

    def load_data(self, df: pd.DataFrame):
        self._df = df
        self._update_run_enabled()

    def get_active_family(self) -> tuple:
        """Returns (family_name, shift_val) for use by the BAM widget."""
        return self._active_family, self._active_shift

    def _update_run_enabled(self):
        has_data = self._df is not None and not self._df.empty
        has_dir  = bool(self._output_dir)
        self._run_button.setEnabled(has_data and has_dir)

    # ── Build UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # ── Info banner ───────────────────────────────────────────────────────
        info = QLabel(
            "<b>Optional: Distribution Family Selection</b><br>"
            "This test fits a simplified model with three candidate distributions and "
            "recommends the best fit based on AIC, BIC, dispersion, deviance explained, "
            "and zero-proportion matching (for families that support zeros).<br><br>"
            "<i>If you skip this step, <b>Tweedie</b> is used by default.</i>"
        )
        info.setWordWrap(True)
        info.setStyleSheet(
            "background: #2a3a2a; border: 1px solid #4a7a4a; "
            "border-radius: 4px; padding: 10px;"
        )
        layout.addWidget(info)

        # ── Run controls ──────────────────────────────────────────────────────
        run_group = QGroupBox("1. Run Family Selection Test")
        rg_layout = QHBoxLayout(run_group)

        self._run_button = QPushButton("Run Family Selection")
        self._run_button.setEnabled(False)
        self._run_button.setMinimumHeight(36)
        self._run_button.setStyleSheet("QPushButton { font-weight: bold; }")
        self._run_button.clicked.connect(self._on_run)

        self._stop_button = QPushButton("Stop")
        self._stop_button.setEnabled(False)
        self._stop_button.clicked.connect(self._on_stop)

        self._status_label = QLabel("Awaiting data and output directory…")

        rg_layout.addWidget(self._run_button)
        rg_layout.addWidget(self._stop_button)
        rg_layout.addWidget(self._status_label, stretch=1)
        layout.addWidget(run_group)

        # ── Splitter: log on top, results on bottom ───────────────────────────
        splitter = QSplitter(Qt.Vertical)

        log_group = QGroupBox("2. R Output Log")
        lg_layout = QVBoxLayout(log_group)
        self._log_text = QTextEdit()
        self._log_text.setReadOnly(True)
        self._log_text.setFont(QFont("Courier New", 9))
        lg_layout.addWidget(self._log_text)
        splitter.addWidget(log_group)

        results_group = QGroupBox("3. Results")
        res_layout    = QVBoxLayout(results_group)

        # Comparison table
        self._table = QTableWidget(0, 7)
        self._table.setHorizontalHeaderLabels([
            "Family", "AIC", "BIC", "Dispersion",
            "Dev. Explained (%)", "Zero Match", "Shift Applied",
        ])
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.setEditTriggers(QTableWidget.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectRows)
        self._table.setAlternatingRowColors(True)
        res_layout.addWidget(self._table)

        # Recommendation text
        res_layout.addWidget(QLabel("<b>Recommendation:</b>"))
        self._interp_text = QLabel(
            "Run the test above to see the recommended family and interpretation."
        )
        self._interp_text.setWordWrap(True)
        self._interp_text.setStyleSheet(
            "background: #1e2a3a; border: 1px solid #3a5a7a; "
            "border-radius: 4px; padding: 10px;"
        )
        res_layout.addWidget(self._interp_text)

        # Family selector (auto-set from test, can be overridden)
        sel_row = QHBoxLayout()
        sel_row.addWidget(QLabel("<b>Family to use for BAM analysis:</b>"))
        self._family_combo = QComboBox()
        self._family_combo.addItem("Tweedie — compound Poisson-Gamma (default)", ("Tweedie",     0.0))
        self._family_combo.addItem("Gamma — log link",                             ("Gamma",       0.0))
        self._family_combo.addItem("Negative Binomial",                            ("NegBinomial", 0.0))
        self._family_combo.currentIndexChanged.connect(self._on_family_combo_changed)
        sel_row.addWidget(self._family_combo)
        sel_row.addStretch()
        res_layout.addLayout(sel_row)

        splitter.addWidget(results_group)
        splitter.setSizes([150, 450])
        layout.addWidget(splitter, stretch=1)

    # ── Slots ─────────────────────────────────────────────────────────────────

    def _on_family_combo_changed(self, index: int):
        family_name, shift_val = self._family_combo.itemData(index)
        self._active_family = family_name
        self._active_shift  = shift_val
        self.family_changed.emit(family_name, shift_val)

    def _on_run(self):
        if self._df is None or self._output_dir is None:
            QMessageBox.warning(
                self, "Not Ready",
                "Please load data and set an output directory in the Data Loading tab first."
            )
            return

        rscript = self._find_rscript()
        if rscript is None:
            QMessageBox.critical(
                self, "Rscript Not Found",
                "Could not locate Rscript.\n\n"
                "Install R from https://cran.r-project.org, then run "
                "the package setup script:\n"
                "  Rscript App/R/install_packages.R"
            )
            return

        if not os.path.isfile(_FAMILY_SELECTION_R):
            QMessageBox.critical(
                self, "R Script Not Found",
                f"FamilySelection.R not found at:\n{_FAMILY_SELECTION_R}"
            )
            return

        # Export CSV if not already present
        csv_path = os.path.join(self._output_dir, "finomena_pre-processed_data.csv")
        if not os.path.isfile(csv_path):
            try:
                df = self._df.copy()
                df["pixel_diff"] = df["pxl_diff"] if "pxl_diff" in df.columns else 0
                df["animal_id"]  = df["plate"].astype(str) + "_" + df["location"].astype(str)
                export_cols = ["time_sec", "location", "loc_coord", "pixel_diff",
                               "Condition", "Phase", "Group", "animal_id", "plate"]
                for var in self._variable_names:
                    if var in df.columns and var not in export_cols:
                        export_cols.append(var)
                missing = [c for c in export_cols if c not in df.columns]
                if missing:
                    QMessageBox.critical(
                        self, "Missing Columns",
                        f"Cannot export CSV — missing columns:\n{missing}"
                    )
                    return
                df[export_cols].to_csv(csv_path, index=False)
            except Exception as e:
                QMessageBox.critical(self, "Export Error", str(e))
                return

        self._run_button.setEnabled(False)
        self._stop_button.setEnabled(True)
        self._status_label.setText("Running…")
        self._log_text.clear()
        self._append_log(f"Script: {_FAMILY_SELECTION_R}")
        self._append_log(f"Input:  {csv_path}")
        self._append_log(f"Output: {self._output_dir}")
        self._append_log("-" * 60)

        var_names_str  = ",".join(self._variable_names)
        ref_values_str = ",".join(
            self._variable_refs.get(v, "") for v in self._variable_names
        )

        cmd = [
            rscript, _FAMILY_SELECTION_R,
            csv_path, self._output_dir,
            var_names_str, ref_values_str, self._ref_condition,
        ]

        kwargs = {}
        if sys.platform == "win32":
            kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW
        else:
            # Put the child in its own process group so killpg can take down
            # the whole R subprocess tree (parallel workers etc.) on Stop or
            # app close.
            kwargs["start_new_session"] = True

        self._process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",   # R emits UTF-8; Windows would otherwise default to cp1252
            errors="replace",   # don't crash if R ever emits a stray non-UTF-8 byte
            bufsize=1,
            env=self._r_env(rscript),
            **kwargs
        )

        threading.Thread(target=self._stream_output, daemon=True).start()

    def request_termination(self) -> bool:
        """
        Force-kill the R process AND its descendants (silent, no UI updates).
        See BamWidget.request_termination for rationale. Used by both the Stop
        button and the MainWindow closeEvent so quitting the app tears down R
        cleanly.
        """
        proc = self._process
        if proc is None or proc.poll() is not None:
            return False
        if sys.platform == "win32":
            try:
                subprocess.run(
                    ["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                    capture_output=True, timeout=5,
                    creationflags=subprocess.CREATE_NO_WINDOW,
                )
            except Exception:
                try:
                    proc.terminate()
                except Exception:
                    pass
        else:
            try:
                pgid = os.getpgid(proc.pid)
                os.killpg(pgid, signal.SIGTERM)
                try:
                    proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    os.killpg(pgid, signal.SIGKILL)
            except Exception:
                try:
                    proc.terminate()
                except Exception:
                    pass
        return True

    def _on_stop(self):
        if self.request_termination():
            self._append_log("--- Stopped by user ---")
        self._run_button.setEnabled(True)
        self._stop_button.setEnabled(False)
        self._status_label.setText("Stopped.")

    def _stream_output(self):
        for line in iter(self._process.stdout.readline, ""):
            self._log_line_ready.emit(line.rstrip())
        self._process.wait()
        self._run_finished.emit(self._process.returncode == 0)

    def _append_log(self, line: str):
        self._log_text.append(line)
        self._log_text.verticalScrollBar().setValue(
            self._log_text.verticalScrollBar().maximum()
        )

    def _on_run_finished(self, success: bool):
        self._run_button.setEnabled(True)
        self._stop_button.setEnabled(False)
        if success:
            self._status_label.setText("Complete ✓")
            self._load_results()
        else:
            self._status_label.setText("Error — see log")
            lines  = self._log_text.toPlainText().splitlines()
            detail = "\n".join(lines[-20:])
            QMessageBox.warning(
                self, "Family Selection Failed",
                "The R script exited with an error.\n\nLast log output:\n\n" + detail
            )

    def _load_results(self):
        results_csv = os.path.join(self._output_dir, "family_selection_results.csv")
        winner_csv  = os.path.join(self._output_dir, "family_selection_winner.csv")

        if not os.path.isfile(results_csv) or not os.path.isfile(winner_csv):
            self._append_log("Result CSVs not found — check the log for errors.")
            return

        try:
            results_df = pd.read_csv(results_csv)
            winner_df  = pd.read_csv(winner_csv)
        except Exception as e:
            self._append_log(f"Could not read results: {e}")
            return

        # Populate comparison table
        cols = ["Family", "AIC", "BIC", "Dispersion", "Dev_Explained", "Zero_Match", "Shift_Applied"]
        self._table.setRowCount(len(results_df))
        winner_name = winner_df["family"].iloc[0]
        for row_idx, row in results_df.iterrows():
            is_winner = row.get("Family") == winner_name
            for col_idx, col in enumerate(cols):
                val  = row.get(col, "")
                item = QTableWidgetItem(str(val))
                item.setTextAlignment(Qt.AlignCenter)
                if is_winner:
                    item.setBackground(Qt.darkGreen)
                self._table.setItem(row_idx, col_idx, item)
        self._table.resizeColumnsToContents()

        # Winner details
        winner_shift = float(winner_df["shift"].iloc[0])
        zero_prop    = float(winner_df["zero_prop_obs"].iloc[0])

        interp = _FAMILY_INTERPRETATION.get(
            winner_name,
            f"{winner_name} was selected as the best-fitting family (lowest AIC)."
        )
        self._interp_text.setText(
            f"<b>Recommended: {_FAMILY_DISPLAY.get(winner_name, winner_name)}</b><br><br>"
            f"{interp}<br><br>"
            f"<i>Observed zero proportion in data: {zero_prop:.1f}%"
            + (f" &nbsp;|&nbsp; Shift applied: {winner_shift:.4f}" if winner_shift > 0 else "")
            + "</i>"
        )

        # Update Gamma's stored shift value in the combo to match what the test found
        gamma_row = results_df[results_df["Family"] == "Gamma"]
        gamma_shift = float(gamma_row["Shift_Applied"].iloc[0]) if not gamma_row.empty else winner_shift
        for i in range(self._family_combo.count()):
            name, _ = self._family_combo.itemData(i)
            if name == "Gamma":
                self._family_combo.setItemData(i, ("Gamma", gamma_shift))
                break

        # Set combo to winner without triggering the signal twice
        family_to_idx = {"Tweedie": 0, "Gamma": 1, "NegBinomial": 2}
        self._family_combo.blockSignals(True)
        self._family_combo.setCurrentIndex(family_to_idx.get(winner_name, 0))
        self._family_combo.blockSignals(False)

        # Directly update active state and emit once
        self._active_family = winner_name
        self._active_shift  = winner_shift
        self.family_changed.emit(winner_name, winner_shift)

    # ── R helpers (mirrors BamWidget) ─────────────────────────────────────────

    @staticmethod
    def _r_env(rscript_path: str) -> dict:
        env      = os.environ.copy()
        existing = env.get("R_LIBS", "")
        env["R_LIBS"] = _APP_R_LIB + (os.pathsep + existing if existing else "")
        if rscript_path == _BUNDLED_RSCRIPT:
            env["R_HOME"] = os.path.dirname(os.path.dirname(rscript_path))  # runtime/bin/../
        return env

    @staticmethod
    def _find_rscript() -> str:
        if os.path.isfile(_BUNDLED_RSCRIPT):
            return _BUNDLED_RSCRIPT
        found = shutil.which("Rscript")
        if found:
            return found
        if sys.platform == "win32":
            for pattern in [
                r"C:\Program Files\R\R-*\bin\Rscript.exe",
                r"C:\Program Files (x86)\R\R-*\bin\Rscript.exe",
            ]:
                matches = _glob.glob(pattern)
                if matches:
                    return sorted(matches)[-1]
        elif sys.platform == "darwin":
            for candidate in (
                "/Library/Frameworks/R.framework/Resources/bin/Rscript",
                "/opt/homebrew/bin/Rscript",
                "/usr/local/bin/Rscript",
            ):
                if os.path.isfile(candidate):
                    return candidate
        else:
            for candidate in ("/usr/bin/Rscript", "/usr/local/bin/Rscript"):
                if os.path.isfile(candidate):
                    return candidate
        return None

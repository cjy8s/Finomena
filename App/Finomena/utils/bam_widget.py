"""
BAM Analysis Widget
===================
Receives the preprocessed full_df from the Python pipeline,
exports it to a CSV compatible with the adapted TweedieAR1 BAM.R script,
runs Rscript as a subprocess with real-time log streaming. Interactive plots
of the results live in the Visualizations tab (visualizations_widget.py),
which reads master_results.csv that R writes alongside other outputs.

R package isolation
-------------------
All required R packages are installed into App/R/library/ (the "app library"),
which is a plain directory inside the project that any user can write to — no
admin privileges required.  When Rscript is launched, R_LIBS is set so R finds
packages there before anywhere else.  This keeps the project fully self-contained
regardless of which system R is installed.
"""

import glob
import json
import os
import shutil
import signal
import subprocess
import sys
import threading

import pandas as pd

from PySide6.QtCore import Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QGroupBox, QHBoxLayout, QLabel, QMessageBox,
    QPushButton, QTextEdit, QVBoxLayout, QWidget
)

# ── Path resolution (works in dev and when frozen by PyInstaller) ─────────────
from paths import resource_path, external_data_dir

# R library lives alongside the executable (too large to bundle inside .exe)
_APP_R_LIB = os.path.join(external_data_dir(), "R", "library")

# Bundled Rscript: lives at R/runtime/bin/ in both dev (App/) and frozen (dist/Finomena/) modes.
# Populated by setup_dev.bat/.sh for development, and by the build scripts for distribution.
_BUNDLED_RSCRIPT = os.path.join(
    external_data_dir(), "R", "runtime", "bin",
    "Rscript.exe" if sys.platform == "win32" else "Rscript"
)

# R script is small — bundled inside the app via resource_path
_DEFAULT_R_SCRIPT = resource_path("R", "scripts", "TweedieAR1 BAM.R")

# Summary CSV that R writes — read at the end of a run so the run log shows
# headline numbers per phase group.
_EXPECTED_CSV = "summary_statistics_by_group.csv"


class BamWidget(QWidget):
    """
    Tab widget for running the TweedieAR1 BAM R analysis.
    """

    # Public signals
    analysis_complete = Signal(str)   # emits output_dir path on success

    # Analysis thread signals
    _log_line_ready = Signal(str)
    _run_finished   = Signal(bool)    # True = success

    def __init__(self, parent=None):
        super().__init__(parent)
        self._df                  = None
        self._csv_path            = None
        self._output_dir          = None
        self._process             = None
        self._variable_names      = ["Genotype", "Drug"]
        self._variable_refs       = {}         # {var_name: ref_value}
        self._ref_condition       = ""
        self._roles               = {}         # {condition_name: role_str}
        self._family_name         = "Tweedie"  # distribution family (set by FamilySelectionWidget)
        self._shift_val           = 0.0        # additive shift applied before fitting
        self._contrast_widget     = None       # reference to ContrastSelectionWidget
        self._correction_widget   = None       # reference to CorrectionWidget
        # Correction strategy comes in via correction.json sidecar (Correction
        # tab writes it). emmeans calls now always run with adjust = "none".
        self._user_stopped        = False      # set by Stop button so _on_run_finished skips the error dialog

        self._build_ui()

        self._log_line_ready.connect(self._append_log)
        self._run_finished.connect(self._on_run_finished)

    # ── Public API ────────────────────────────────────────────────────────────

    def set_variable_names(self, names: list):
        """Stores the descriptor variable names (e.g. ['Genotype', 'Drug', 'Antidote'])."""
        self._variable_names = list(names)

    def set_references(self, variable_refs: dict, ref_condition: str):
        """Stores the reference/control levels to be passed to the R script."""
        self._variable_refs = dict(variable_refs)
        self._ref_condition = ref_condition

    def set_contrast_selection_widget(self, widget):
        """Reference to ContrastSelectionWidget; its kept pairs are read at run time."""
        self._contrast_widget = widget

    def set_correction_widget(self, widget):
        """Reference to CorrectionWidget; its spec is serialized to
        correction.json in the output dir right before the R subprocess
        launches."""
        self._correction_widget = widget

    def on_contrasts_changed(self):
        """
        Called when the user modifies the contrast selection. Marks any prior
        BAM results as stale: resets the status, and shows a warning prompting
        the user to re-run.
        """
        self._r_status_label.setText("")
        self._stale_warning_label.show()

    def reset_stale_warning(self):
        """Hide the 'contrasts changed' warning.

        Called by the app after loading an experiment config so the warning
        does not flash up just from the cascade of signals fired while
        restoring conditions/roles/contrast-selection state.
        """
        self._stale_warning_label.hide()

    def set_roles(self, roles: dict):
        """Stores the {condition_name: role_str} mapping."""
        self._roles = dict(roles)

    def set_family(self, family_name: str, shift_val: float):
        """Receives the chosen distribution family from FamilySelectionWidget."""
        self._family_name = family_name
        self._shift_val   = shift_val
        labels = {
            "Tweedie":     "Tweedie (compound Poisson-Gamma)",
            "Gamma":       "Gamma (log link)",
            "NegBinomial": "Negative Binomial",
        }
        display = labels.get(family_name, family_name)
        shift_note = f" — shift: {shift_val:.4f}" if shift_val != 0 else ""
        self._family_label.setText(f"Family: <b>{display}</b>{shift_note}")

    def set_output_dir(self, path: str):
        """Sets the output directory for R results (called from MainWindow)."""
        self._output_dir = path
        self._update_run_enabled()

    def load_data(self, df: pd.DataFrame):
        """Receives the phase-assigned full_df from the pipeline."""
        self._df = df
        self._update_run_enabled()

    def _update_run_enabled(self):
        """Enables the Run button when both data and output directory are available."""
        has_data = self._df is not None and not self._df.empty
        has_dir = bool(self._output_dir)
        self._run_button.setEnabled(has_data and has_dir)

    # ── Build UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # Correction method choice now lives in the dedicated Correction tab.
        # The BAM widget only runs the model and reports it back; correction
        # gets applied post-hoc against raw emmeans p-values. The within-
        # emmeans Contrast Correction stage (Dunnett/Šidák) is gone entirely.
        info_label = QLabel(
            "<b>Multiple-testing correction:</b> configured in the "
            "<i>Correction</i> tab. This stage produces raw emmeans p-values "
            "that the chosen correction (flat BH/Holm or tree-based "
            "TreeBH/graphicalMCP/Holm-gatekeeping) operates on."
        )
        info_label.setWordWrap(True)
        info_label.setStyleSheet("padding: 6px; color: #aaaaaa;")
        layout.addWidget(info_label)

        # ── 2: Run BAM Analysis ────────────────────────────────────────────────
        row_23 = QHBoxLayout()

        run_group = QGroupBox("2. Run BAM Analysis")
        rg_layout = QVBoxLayout(run_group)

        self._family_label = QLabel("Family: <b>Tweedie (compound Poisson-Gamma)</b> — default")
        self._family_label.setStyleSheet("font-style: italic; color: #aaaaaa;")
        rg_layout.addWidget(self._family_label)

        # Warning shown when contrast selection changes after a BAM run
        self._stale_warning_label = QLabel(
            "⚠ Contrasts changed since the last run — re-run the BAM analysis to refresh results."
        )
        self._stale_warning_label.setStyleSheet(
            "color: #c0392b; font-weight: bold; padding: 4px;"
        )
        self._stale_warning_label.setWordWrap(True)
        self._stale_warning_label.hide()
        rg_layout.addWidget(self._stale_warning_label)

        run_row = QHBoxLayout()
        self._run_button = QPushButton("Run R Analysis")
        self._run_button.setEnabled(False)
        self._run_button.setMinimumHeight(36)
        self._run_button.setStyleSheet("QPushButton { font-weight: bold; }")
        self._run_button.clicked.connect(self._on_run)
        self._stop_button = QPushButton("Stop")
        self._stop_button.setEnabled(False)
        self._stop_button.clicked.connect(self._on_stop)
        self._r_status_label = QLabel("")
        run_row.addWidget(self._run_button)
        run_row.addWidget(self._stop_button)
        run_row.addWidget(self._r_status_label, stretch=1)
        rg_layout.addLayout(run_row)

        row_23.addWidget(run_group, stretch=1)
        layout.addLayout(row_23)

        # ── 3. R Output Log (no figure viewer — see Visualizations tab) ───────
        log_group = QGroupBox("3. R Output Log")
        lg_layout = QVBoxLayout(log_group)
        self._log_text = QTextEdit()
        self._log_text.setReadOnly(True)
        mono_font = QFont("Courier New", 9)
        self._log_text.setFont(mono_font)
        lg_layout.addWidget(self._log_text)
        layout.addWidget(log_group, stretch=1)

    # ── Run R ─────────────────────────────────────────────────────────────────

    def _export_csv(self) -> str:
        """Exports the pre-processed data CSV to the output directory. Returns path or None."""
        if self._df is None or self._df.empty or not self._output_dir:
            return None
        try:
            df = self._df.copy()
            df['pixel_diff'] = df['pxl_diff'] if 'pxl_diff' in df.columns else 0
            df['animal_id'] = df['plate'].astype(str) + "_" + df['location'].astype(str)
            export_cols = ['time_sec', 'location', 'loc_coord', 'pixel_diff',
                           'Condition', 'Phase', 'Group', 'animal_id', 'plate']
            # Also include any individual variable columns present in the dataframe
            for var in self._variable_names:
                if var in df.columns and var not in export_cols:
                    export_cols.append(var)
            missing = [c for c in export_cols if c not in df.columns]
            if missing:
                QMessageBox.critical(
                    self, "Missing Columns",
                    f"Cannot export — missing columns:\n{missing}")
                return None
            path = os.path.join(self._output_dir, "finomena_pre-processed_data.csv")
            df[export_cols].to_csv(path, index=False)
            self._csv_path = path
            return path
        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))
            return None

    def _on_run(self):
        if self._df is None or self._output_dir is None:
            QMessageBox.warning(
                self, "Not Ready",
                "Please load data and set an output directory in the Data Loading tab first."
            )
            return

        # Step 1: Export CSV to output directory
        csv_path = self._export_csv()
        if not csv_path:
            return

        output_dir = self._output_dir

        # Step 1b: Write contrast-selection sidecar JSON. The R script looks for
        # this file in the output directory and, if present:
        #   - filters Full_Interaction comparisons to the listed pairs
        #     (kept_pairs non-empty), or skips the family entirely (empty)
        #   - if posterior_equivalence.enabled, also runs Bayesian equivalence
        #     testing on those same pairs and emits posterior_equivalence_*.png
        #     and posterior_equivalence_summary.csv
        # If the sidecar is missing, R falls back to full pairwise (legacy).
        if self._contrast_widget is not None:
            try:
                sidecar_data = {
                    "kept_pairs":           self._contrast_widget.get_kept_pairs(),
                    "posterior_equivalence": self._contrast_widget.get_posterior_settings(),
                }
                sidecar = os.path.join(output_dir, "contrast_selection.json")
                with open(sidecar, "w", encoding="utf-8") as f:
                    json.dump(sidecar_data, f, indent=2)
            except Exception as e:
                QMessageBox.warning(
                    self, "Contrast Sidecar",
                    f"Could not write contrast_selection.json:\n{e}\n\n"
                    "R will fall back to full pairwise comparisons."
                )

        # Step 1c: Write correction.json sidecar from the Correction tab.
        # Without this the R script falls back to flat BH per family at q=0.05.
        if self._correction_widget is not None:
            try:
                self._correction_widget.set_output_dir(output_dir)
                path = self._correction_widget.write_sidecar()
                if not path:
                    QMessageBox.warning(
                        self, "Correction Sidecar",
                        "Could not write correction.json — falling back to "
                        "flat BH at q=0.05."
                    )
            except Exception as e:
                QMessageBox.warning(
                    self, "Correction Sidecar",
                    f"Could not write correction.json:\n{e}\n\n"
                    "R will fall back to flat BH at q=0.05."
                )

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

        if not os.path.isfile(_DEFAULT_R_SCRIPT):
            QMessageBox.critical(
                self, "R Script Not Found",
                f"The BAM R script was not found at:\n{_DEFAULT_R_SCRIPT}"
            )
            return

        self._run_button.setEnabled(False)
        self._stop_button.setEnabled(True)
        self._user_stopped = False
        self._r_status_label.setText("Running…")
        self._stale_warning_label.hide()
        self._log_text.clear()
        self._append_log(f"Rscript: {rscript}")
        self._append_log(f"Script:  {_DEFAULT_R_SCRIPT}")
        self._append_log("-" * 60)
        # The R script echoes the rest of the args (input, output, vars, refs,
        # roles, family, etc.) on startup — no need to duplicate them here.
        var_names_str = ",".join(self._variable_names)
        ref_values_str = ",".join(
            self._variable_refs.get(v, "") for v in self._variable_names
        )
        # Build roles string: "WT+DMSO=RC,KO+DMSO=EXP,..."
        _ROLE_ABBREV = {
            "Reference Control": "RC",
            "Experimental":      "EXP",
            "Positive Control":  "PC",
            "Potential Rescue":  "RES",
        }
        roles_str = ",".join(
            f"{cond}={_ROLE_ABBREV.get(role, role)}"
            for cond, role in self._roles.items()
            if role  # skip empty/unassigned
        )

        # Correction method/contrast-adjust removed from the CLI: the R
        # script always uses adjust = "none" on emmeans, and the correction
        # strategy is read from correction.json (written by the Correction
        # tab). Positional args kept for back-compat where the R script just
        # ignores them now.
        cmd = [
            rscript, _DEFAULT_R_SCRIPT,
            self._csv_path, output_dir,
            var_names_str, ref_values_str, self._ref_condition,
            "BH", "none",            # legacy args, ignored by the R script
            roles_str,
            self._family_name, str(self._shift_val),
        ]

        kwargs = {}
        if sys.platform == "win32":
            kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW
        else:
            # Put the child in its own process group so we can SIGTERM/SIGKILL
            # the whole tree (R parallel workers etc.) on Stop.
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

        t = threading.Thread(target=self._stream_output, daemon=True)
        t.start()

    def request_termination(self) -> bool:
        """
        Force-kill the R process AND its descendants (silent, no UI updates).

        Plain Popen.terminate() is unreliable here:
          • Windows: it only ends the immediate Rscript.exe; any worker
            processes spawned by R (parallel:: or BLAS threads using their own
            processes) keep running until they finish naturally.
          • Linux/macOS: same problem unless we kill the whole process group.
        We use taskkill /F /T on Windows and SIGKILL on the process group on
        POSIX to guarantee the whole tree is gone.

        Safe to call at any time. Returns True if a process was killed,
        False if there was nothing running. Used by the Stop button AND by
        the MainWindow closeEvent so quitting the app tears down R cleanly.
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
        """Stop button: tear down R and update the UI."""
        self._user_stopped = True
        if not self.request_termination():
            self._stop_button.setEnabled(False)
            return
        self._append_log("--- Stopping R analysis (killing process tree) ---")
        self._run_button.setEnabled(True)
        self._stop_button.setEnabled(False)
        self._r_status_label.setText("Stopped.")

    def _stream_output(self):
        """Background thread: streams stdout from Rscript line by line."""
        for line in iter(self._process.stdout.readline, ""):
            self._log_line_ready.emit(line.rstrip())
        self._process.wait()
        self._run_finished.emit(self._process.returncode == 0)

    def _append_log(self, line: str):
        """Slot: appends a line to the log and auto-scrolls."""
        self._log_text.append(line)
        self._log_text.verticalScrollBar().setValue(
            self._log_text.verticalScrollBar().maximum()
        )

    def _on_run_finished(self, success: bool):
        self._run_button.setEnabled(True)
        self._stop_button.setEnabled(False)
        if self._user_stopped:
            # User clicked Stop — _on_stop already set the status label and
            # logged the action. Don't pop a "failed" dialog for an intentional cancel.
            return
        if success:
            self._r_status_label.setText("Complete ✓")
            self._load_summary()
            self.analysis_complete.emit(self._output_dir)
        else:
            self._r_status_label.setText("Error — see log")
            lines  = self._log_text.toPlainText().splitlines()
            detail = "\n".join(lines[-20:])
            QMessageBox.warning(
                self, "R Analysis Failed",
                "The R script exited with an error.\n\nLast log output:\n\n" + detail
            )

    # ── Load summary CSV into the log ─────────────────────────────────────────

    def _load_summary(self):
        """Reads the per-group summary CSV and echoes it into the run log.

        Interactive figures live in the Visualizations tab — this method only
        surfaces headline counts/effects from summary_statistics_by_group.csv
        so the user gets immediate feedback inline with the R log.
        """
        csv_path = os.path.join(self._output_dir, _EXPECTED_CSV)
        if os.path.isfile(csv_path):
            try:
                summary = pd.read_csv(csv_path)
                self._append_log("\n--- Summary Statistics ---")
                self._append_log(summary.to_string())
            except Exception as exc:
                self._append_log(f"Could not read summary CSV: {exc}")
        else:
            self._append_log(f"[Not found] {_EXPECTED_CSV}")

        master_path = os.path.join(self._output_dir, "master_results.csv")
        if os.path.isfile(master_path):
            self._append_log(
                "\nmaster_results.csv written — open the Visualizations tab "
                "to explore the contrast results."
            )

    # ── Helpers ───────────────────────────────────────────────────────────────

    @staticmethod
    def _r_env(rscript_path: str) -> dict:
        """
        Returns an environment dict for subprocess that prepends the app's
        R library to R_LIBS so packages installed there are found first,
        regardless of any system-wide R configuration.
        """
        env = os.environ.copy()
        existing = env.get("R_LIBS", "")
        env["R_LIBS"] = (
            _APP_R_LIB + (os.pathsep + existing if existing else "")
        )
        # If using the bundled R, set R_HOME so R finds its own base packages in runtime/library/
        if rscript_path == _BUNDLED_RSCRIPT:
            env["R_HOME"] = os.path.dirname(os.path.dirname(rscript_path))  # runtime/bin/../
        return env

    @staticmethod
    def _find_rscript() -> str:
        """Returns the path to Rscript, preferring the bundled copy."""
        # 1. Bundled R inside the project / alongside the exe
        if os.path.isfile(_BUNDLED_RSCRIPT):
            return _BUNDLED_RSCRIPT
        # 2. System PATH
        found = shutil.which("Rscript")
        if found:
            return found
        # 3. Platform-specific common install locations
        if sys.platform == "win32":
            for pattern in [
                r"C:\Program Files\R\R-*\bin\Rscript.exe",
                r"C:\Program Files (x86)\R\R-*\bin\Rscript.exe",
            ]:
                matches = glob.glob(pattern)
                if matches:
                    return sorted(matches)[-1]
        elif sys.platform == "darwin":
            # CRAN .pkg installer, Homebrew Apple Silicon, Homebrew Intel
            for candidate in (
                "/Library/Frameworks/R.framework/Resources/bin/Rscript",
                "/opt/homebrew/bin/Rscript",
                "/usr/local/bin/Rscript",
            ):
                if os.path.isfile(candidate):
                    return candidate
        else:
            # Linux
            for candidate in ("/usr/bin/Rscript", "/usr/local/bin/Rscript"):
                if os.path.isfile(candidate):
                    return candidate
        return None

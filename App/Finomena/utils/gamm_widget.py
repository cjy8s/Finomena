"""
GAMM Analysis Widget
====================
Receives the preprocessed full_df from the Python pipeline,
exports it to a CSV compatible with the adapted TweedieAR1 GAMM.R script,
runs Rscript as a subprocess with real-time log streaming,
and displays the resulting PNG figures inline.

R package isolation
-------------------
All required R packages are installed into App/R/library/ (the "app library"),
which is a plain directory inside the project that any user can write to — no
admin privileges required.  When Rscript is launched, R_LIBS is set so R finds
packages there before anywhere else.  This keeps the project fully self-contained
regardless of which system R is installed.
"""

import glob
import os
import shutil
import subprocess
import sys
import threading

import pandas as pd

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import (
    QComboBox, QGroupBox, QHBoxLayout, QLabel, QMessageBox,
    QPushButton, QSizePolicy, QSplitter, QTextEdit, QVBoxLayout, QWidget
)

from figure_viewer import FigureViewerWidget

# ── Path resolution (works in dev and when frozen by PyInstaller) ─────────────
from paths import resource_path, external_data_dir

# R library lives alongside the executable (too large to bundle inside .exe)
_APP_R_LIB = os.path.join(external_data_dir(), "R", "library")

# If a bundled Rscript lives inside the project, prefer it over system R.
_BUNDLED_RSCRIPT = os.path.join(
    external_data_dir(), "R", "bin",
    "Rscript.exe" if sys.platform == "win32" else "Rscript"
)

# R script is small — bundled inside the app via resource_path
_DEFAULT_R_SCRIPT = resource_path("R", "scripts", "TweedieAR1 GAMM.R")

# Expected output files from the R script (static ones; per-variable plots are discovered dynamically)
_EXPECTED_PNGS_STATIC = [
    "forest_plot_interactions.png",
    "heatmap_diverging_v1.png",
    "rescue_assessment_context_faceted.png",
]
_EXPECTED_CSV = "summary_statistics_by_group.csv"


class GammWidget(QWidget):
    """
    Tab widget for running the TweedieAR1 GAMM R analysis.
    """

    # Public signal
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
        self._global_correction   = "BH"       # default global multiple-testing method
        self._contrast_correction = "dunnett"  # default within-contrast adjustment

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

    def set_roles(self, roles: dict):
        """Stores the {condition_name: role_str} mapping."""
        self._roles = dict(roles)

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

        # ── 1. Statistical Correction Methods ────────────────────────────────
        corr_group = QGroupBox("1. Statistical Correction Methods")
        corr_outer = QHBoxLayout(corr_group)

        # ── Left: Global correction ───────────────────────────────────────────
        global_col = QVBoxLayout()
        global_col.addWidget(QLabel("<b>Global Multiple-Testing Correction</b><br>"
                                    "<small>Applied across all group × contrast results.</small>"))

        self._global_combo = QComboBox()
        self._global_combo.addItem("Benjamini-Hochberg (BH)  —  controls FDR", "BH")
        self._global_combo.addItem("Holm  —  controls FWER (step-down Bonferroni)", "holm")
        self._global_combo.currentIndexChanged.connect(self._on_global_correction_changed)
        global_col.addWidget(self._global_combo)

        self._global_desc = QLabel()
        self._global_desc.setWordWrap(True)
        self._global_desc.setStyleSheet("font-style: italic;")
        global_col.addWidget(self._global_desc)
        global_col.addStretch()
        corr_outer.addLayout(global_col, stretch=1)

        # ── Right: Contrast (within-emmeans) correction ───────────────────────
        contrast_col = QVBoxLayout()
        contrast_col.addWidget(QLabel("<b>Contrast Correction (Treatment vs. Reference)</b><br>"
                                      "<small>Applied inside each emmeans contrast call.</small>"))

        self._contrast_combo = QComboBox()
        self._contrast_combo.addItem("Dunnett's Test", "dunnett")
        self._contrast_combo.addItem("Dunnett-Šidák", "sidak")
        self._contrast_combo.currentIndexChanged.connect(self._on_contrast_correction_changed)
        contrast_col.addWidget(self._contrast_combo)

        self._contrast_desc = QLabel()
        self._contrast_desc.setWordWrap(True)
        self._contrast_desc.setStyleSheet("font-style: italic;")
        contrast_col.addWidget(self._contrast_desc)
        contrast_col.addStretch()
        corr_outer.addLayout(contrast_col, stretch=1)

        layout.addWidget(corr_group)

        # Populate description labels with initial text
        self._on_global_correction_changed(0)
        self._on_contrast_correction_changed(0)

        # ── 2: Run GAMM Analysis ───────────────────────────────────────────────
        row_23 = QHBoxLayout()

        run_group = QGroupBox("2. Run GAMM Analysis")
        rg_layout = QHBoxLayout(run_group)
        self._run_button = QPushButton("Run R Analysis")
        self._run_button.setEnabled(False)
        self._run_button.setMinimumHeight(36)
        self._run_button.setStyleSheet("QPushButton { font-weight: bold; }")
        self._run_button.clicked.connect(self._on_run)
        self._stop_button = QPushButton("Stop")
        self._stop_button.setEnabled(False)
        self._stop_button.clicked.connect(self._on_stop)
        self._r_status_label = QLabel("")
        rg_layout.addWidget(self._run_button)
        rg_layout.addWidget(self._stop_button)
        rg_layout.addWidget(self._r_status_label, stretch=1)

        row_23.addWidget(run_group, stretch=1)
        layout.addLayout(row_23)

        # ── 4 + 5: Resizable splitter between log and figures ─────────────────
        splitter = QSplitter(Qt.Vertical)

        log_group = QGroupBox("3. R Output Log")
        lg_layout = QVBoxLayout(log_group)
        self._log_text = QTextEdit()
        self._log_text.setReadOnly(True)
        mono_font = QFont("Courier New", 9)
        self._log_text.setFont(mono_font)
        lg_layout.addWidget(self._log_text)
        splitter.addWidget(log_group)

        results_group = QGroupBox("4. Results")
        res_layout = QVBoxLayout(results_group)
        self._figure_viewer = FigureViewerWidget(title="GAMM Figures")
        res_layout.addWidget(self._figure_viewer)
        splitter.addWidget(results_group)

        # Start with log ~25% and figures ~75%
        splitter.setSizes([150, 450])
        layout.addWidget(splitter, stretch=1)

    # ── Correction method slots ───────────────────────────────────────────────

    _GLOBAL_DESCRIPTIONS = {
        "BH": (
            "Benjamini-Hochberg (False Discovery Rate): among all comparisons you "
            "call significant, at most 5% are expected to be false positives on "
            "average. The chance of making even one false positive is typically "
            "above 5%, but you gain more power to detect real effects. Best when "
            "exploring many comparisons and accepting that a small fraction of "
            "hits may be wrong — follow-up experiments will validate them."
        ),
        "holm": (
            "Holm (step-down Bonferroni, Family-Wise Error Rate): the probability "
            "of making even one false positive across all your comparisons is kept "
            "below 5%. More conservative than BH — you will miss more real effects, "
            "but nearly every result you call significant will be a true positive. "
            "Best when a single false positive has serious consequences, such as "
            "claiming a drug works when it does not."
        ),
    }

    _CONTRAST_DESCRIPTIONS = {
        "dunnett": (
            "Dunnett's Test: designed specifically for the many-to-one design — "
            "comparing multiple treatment groups against a single reference/control. "
            "The probability of making even one false positive across all "
            "treatment-vs-reference comparisons is kept below 5%. It exploits the "
            "shared correlation structure (all comparisons involve the same "
            "reference group), making it more powerful than Bonferroni or Holm "
            "for this design. Gold standard for pharmacology dose-response "
            "experiments comparing each drug condition to a vehicle control."
        ),
        "sidak": (
            "Dunnett-Šidák: serves the same purpose as Dunnett's — all contrasts "
            "compare treatments against the reference group only, and the "
            "probability of even one false positive is kept below 5%. Uses the "
            "Šidák inequality instead of the exact multivariate t-distribution, "
            "making it very slightly more conservative (misses a marginally "
            "larger fraction of true effects). The practical difference is "
            "negligible for most datasets. Appropriate as a recognized alternative "
            "to Dunnett's when you want an auditable, widely cited correction that "
            "still anchors every comparison to the reference group."
        ),
    }

    def _on_global_correction_changed(self, index: int):
        key = self._global_combo.itemData(index)
        self._global_correction = key
        self._global_desc.setText(self._GLOBAL_DESCRIPTIONS.get(key, ""))

    def _on_contrast_correction_changed(self, index: int):
        key = self._contrast_combo.itemData(index)
        self._contrast_correction = key
        self._contrast_desc.setText(self._CONTRAST_DESCRIPTIONS.get(key, ""))

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

        rscript = self._find_rscript()
        if rscript is None:
            QMessageBox.critical(
                self, "Rscript Not Found",
                "Could not locate Rscript.\n\n"
                "Install R from https://cran.r-project.org, restart the app, "
                "then click 'Install R Packages'."
            )
            return

        if not os.path.isfile(_DEFAULT_R_SCRIPT):
            QMessageBox.critical(
                self, "R Script Not Found",
                f"The GAMM R script was not found at:\n{_DEFAULT_R_SCRIPT}"
            )
            return

        self._run_button.setEnabled(False)
        self._stop_button.setEnabled(True)
        self._r_status_label.setText("Running…")
        self._log_text.clear()
        self._append_log(f"Rscript: {rscript}")
        self._append_log(f"Script:  {_DEFAULT_R_SCRIPT}")
        self._append_log(f"Input:   {self._csv_path}")
        self._append_log(f"Output:  {output_dir}")
        var_names_str = ",".join(self._variable_names)
        ref_values_str = ",".join(
            self._variable_refs.get(v, "") for v in self._variable_names
        )
        self._append_log(f"Vars:    {var_names_str}")
        self._append_log(f"Refs:    {ref_values_str}  Condition={self._ref_condition!r}")
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
        self._append_log(
            f"Stats:   Global={self._global_correction!r}  "
            f"Contrast={self._contrast_correction!r}"
        )
        self._append_log(f"Roles:   {roles_str}")
        self._append_log(f"R_LIBS:  {_APP_R_LIB}")
        self._append_log("-" * 60)

        cmd = [
            rscript, _DEFAULT_R_SCRIPT,
            self._csv_path, output_dir,
            var_names_str, ref_values_str, self._ref_condition,
            self._global_correction, self._contrast_correction,
            roles_str,
        ]

        kwargs = {}
        if sys.platform == "win32":
            kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW

        self._process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=self._r_env(rscript),
            **kwargs
        )

        t = threading.Thread(target=self._stream_output, daemon=True)
        t.start()

    def _on_stop(self):
        if self._process:
            self._process.terminate()
            self._append_log("--- Analysis stopped by user ---")
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
        if success:
            self._r_status_label.setText("Complete ✓")
            self._load_output_figures()
            self.analysis_complete.emit(self._output_dir)
        else:
            self._r_status_label.setText("Error — see log")
            lines  = self._log_text.toPlainText().splitlines()
            detail = "\n".join(lines[-20:])
            QMessageBox.warning(
                self, "R Analysis Failed",
                "The R script exited with an error.\n\nLast log output:\n\n" + detail
            )

    # ── Load output figures ───────────────────────────────────────────────────

    @staticmethod
    def _make_r_name(name: str) -> str:
        """Mimics R's make.names(): replaces non-alphanumeric chars with dots."""
        import re
        s = re.sub(r'[^A-Za-z0-9.]', '.', name)
        if s and s[0].isdigit():
            s = 'X' + s
        return s

    def _load_output_figures(self):
        """Loads R's PNG outputs directly as QPixmap."""
        # Build the full list: per-variable forest plots + static plots
        # R uses make.names() + tolower() for filenames, so we must match
        expected_pngs = []
        for v in self._variable_names:
            r_name = self._make_r_name(v).lower()
            expected_pngs.append(f"forest_plot_{r_name}_effect.png")
        expected_pngs.extend(_EXPECTED_PNGS_STATIC)

        figures = []
        for png_name in expected_pngs:
            path = os.path.join(self._output_dir, png_name)
            if not os.path.isfile(path):
                self._append_log(f"[Not found] {png_name}")
                continue
            try:
                pixmap = QPixmap(path)
                if pixmap.isNull():
                    self._append_log(f"[Error]     {png_name}: failed to load image")
                    continue
                figures.append({
                    'pixmap':   pixmap,
                    'title':    png_name,
                    'filepath': path,
                })
                self._append_log(f"[Loaded]    {png_name}")
            except Exception as exc:
                self._append_log(f"[Error]     {png_name}: {exc}")

        if figures:
            self._figure_viewer.load_figures(figures)
        else:
            self._append_log("No figures were rendered inline.")

        csv_path = os.path.join(self._output_dir, _EXPECTED_CSV)
        if os.path.isfile(csv_path):
            try:
                summary = pd.read_csv(csv_path)
                self._append_log("\n--- Summary Statistics ---")
                self._append_log(summary.to_string())
            except Exception as exc:
                self._append_log(f"Could not read summary CSV: {exc}")

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
        # If using a bundled R, also set R_HOME so base packages resolve correctly
        bundled_home = os.path.dirname(os.path.dirname(rscript_path))  # bin/../
        if os.path.isfile(_BUNDLED_RSCRIPT) and rscript_path == _BUNDLED_RSCRIPT:
            env["R_HOME"] = bundled_home
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

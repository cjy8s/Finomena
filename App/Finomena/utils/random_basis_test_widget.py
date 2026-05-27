"""
Architecture × Method Ablation Widget (diagnostic)
==================================================
Sweeps three architectures (per-phase / global / HGAM) × random-effect basis
× by-variable × AR.start mode, computes a kitchen-sink set of fit/adequacy/
identifiability/random-eats-fixed metrics per fit, then optionally runs
5-fold animal-stratified CV on the Pareto front of each architecture for a
final cross-architecture ranking by holdout deviance.

Held constant across all variants: select=TRUE, discrete=TRUE, rho=0.2,
family (from FamilySelectionWidget), AR.start logic (animal vs phase resets
varied as an axis for B/C; fixed to animal-only for A), m = mgcv default.

Self-contained: deleting this file + removing the import/tab/closeEvent
lines in app.py removes the diagnostic without touching anything else. The
R script's ablation env vars default to empty (= normal run), so the rest
of the app keeps working unchanged.
"""

import json
import math
import os
import signal
import subprocess
import sys
import threading

import pandas as pd

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QFont
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QGroupBox, QHBoxLayout, QHeaderView, QLabel,
    QMessageBox, QPushButton, QScrollArea, QSizePolicy, QSplitter,
    QTableWidget, QTableWidgetItem, QTabWidget, QTextEdit, QVBoxLayout, QWidget,
)

# matplotlib via QtAgg backend — same pattern as visualizations_widget.py
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Reuse BAM widget's path resolution + Rscript discovery + env builder
from bam_widget import (
    BamWidget,
    _APP_R_LIB,
    _BUNDLED_RSCRIPT,
    _DEFAULT_R_SCRIPT,
)

# ── Axis options ──────────────────────────────────────────────────────────
_RE_BASES         = ("none", "re", "re_slope", "fs", "sz")
_BY_VARIABLES     = ("cond_combo", "group_factor", "interaction")
_AR_START_MODES   = ("animal", "phase")
_ARCHITECTURES    = ("A", "B", "C")

_ARCH_LABELS = {
    "A": "Arch A — Per-Phase BAMs (one model per phase Group)",
    "B": "Arch B — Fully Global (smooths on time_sec, no phase resets)",
    "C": "Arch C — HGAM (smooths on time_in_group, phase-local)",
}

_RE_LABELS = {
    "none":     "none — no animal/plate random effect (AR(1) only baseline)",
    "re":       "re — pure random intercept",
    "re_slope": "re_slope — intercept + random slope",
    "fs":       "fs — unconstrained factor smooth",
    "sz":       "sz — sum-to-zero constrained smooth",
}

_BY_LABELS = {
    "cond_combo":   "cond_combo — by=Condition_Combo (Pedersen-canonical single factor)",
    "group_factor": "group_factor — by=Group_factor (phase-only trajectory)",
    "interaction":  "interaction — by=interaction(C, G) (per-cell)",
}

_AR_LABELS = {
    "animal": "animal — restart AR(1) only at animal boundaries",
    "phase":  "phase — restart AR(1) at animal AND phase boundaries",
}

_ROLE_ABBREV = {
    "Reference Control": "RC",
    "Experimental":      "EXP",
    "Positive Control":  "PC",
    "Potential Rescue":  "RES",
}

# Fixed rho during ablation — keeps all variants on the same AR(1) footing
_ABLATION_RHO = "0.2"

# CV folds
_N_CV_FOLDS = 5

# Adequacy thresholds (loose — only catastrophic failures get flagged)
_ADEQUACY = {
    "k_index_pass_rate_min":     0.30,   # < 30% of smooths passing → fail
    "max_overall_concurvity":    0.97,   # > 0.97 → fail
    "max_random_vs_fixed_conc":  0.95,   # > 0.95 → fail
    "max_ar1_residual_lag1":     0.70,   # > 0.7 → fail
    "phi_relative_multiplier":   100.0,  # > 100× min(phi across variants) → fail
}

# Pareto direction: True = minimize, False = maximize
_PARETO_METRICS = {
    "AIC":                  True,   # minimize
    "BIC":                  True,   # minimize
    "dev_explained":        False,  # maximize
    "r_sq_adj":             False,  # maximize
    "max_random_vs_fixed":  True,   # minimize (random-eats-fixed)
}


# ─── Helper: Pareto front identification ─────────────────────────────────
def _pareto_front(records, metrics):
    """records: list of dicts. metrics: dict {name: True/False (minimize)}.
    Returns indices of non-dominated records. NaN values are treated as
    'worst' so a candidate with NaN can't dominate one with a finite value
    on that metric."""
    valid = []
    for i, r in enumerate(records):
        vals = []
        ok = True
        for m, _ in metrics.items():
            v = r.get(m)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                ok = False
                break
            vals.append(v)
        if ok:
            valid.append((i, vals))
    if not valid:
        return []
    metric_dirs = [metrics[m] for m in metrics]  # True = minimize

    def dominates(a, b):
        strictly_better = False
        for av, bv, mn in zip(a, b, metric_dirs):
            if mn:  # minimize: a beats b if av <= bv
                if av > bv:
                    return False
                if av < bv:
                    strictly_better = True
            else:   # maximize
                if av < bv:
                    return False
                if av > bv:
                    strictly_better = True
        return strictly_better

    out = []
    for ia, va in valid:
        dominated = False
        for ib, vb in valid:
            if ib == ia:
                continue
            if dominates(vb, va):
                dominated = True
                break
        if not dominated:
            out.append(ia)
    return out


def _fmt(v, places=4):
    """Format a numeric value for table display, returning '—' for NA/None."""
    if v is None:
        return "—"
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return "—"
    if isinstance(v, (int,)):
        return str(v)
    try:
        return f"{float(v):.{places}g}"
    except (TypeError, ValueError):
        return str(v)


def _safe_get(d, *keys, default=None):
    """Nested dict get — returns default on KeyError/TypeError/None."""
    cur = d
    for k in keys:
        if cur is None or not isinstance(cur, dict):
            return default
        cur = cur.get(k)
    return cur if cur is not None else default


class AblationTestWidget(QWidget):
    """Architecture × method ablation runner with Pareto + CV selection."""

    _log_line_ready  = Signal(str)
    _stage1_done     = Signal(str, bool)   # (run_id, success)
    _stage3_done     = Signal(str, bool)   # (run_id, success)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._bam_widget        = None
        self._contrast_widget   = None
        self._correction_widget = None

        # Run state
        self._process            = None
        self._user_stopped       = False
        self._current_run_id     = None
        self._current_stage      = None     # "stage1" or "stage3"

        # Stage 1
        self._stage1_queue       = []       # list of dicts (run plans)
        self._stage1_total       = 0
        self._stage1_metrics     = {}       # run_id -> metrics dict
        self._stage1_master      = {}       # run_id -> master_results.csv DataFrame
        self._stage1_fail_reason = {}       # run_id -> str

        # Pareto / Stage 3
        self._pareto_candidates  = {}       # arch -> list of run_id (within-arch Pareto)
        self._pareto_union       = []       # list of run_id (across all arch)
        self._stage3_queue       = []
        self._stage3_total       = 0
        self._stage3_results     = {}       # run_id -> list of (fold, deviance)

        self._build_ui()
        self._log_line_ready.connect(self._append_log)
        self._stage1_done.connect(self._on_stage1_done)
        self._stage3_done.connect(self._on_stage3_done)

    # ── Public wiring (called from MainWindow) ────────────────────────────
    def set_bam_widget(self, w):        self._bam_widget = w
    def set_contrast_widget(self, w):   self._contrast_widget = w
    def set_correction_widget(self, w): self._correction_widget = w

    # ── UI construction ───────────────────────────────────────────────────

    def _build_ui(self):
        main = QVBoxLayout(self)

        intro = QLabel(
            "<b>Architecture × Method Ablation (diagnostic)</b><br>"
            "Sweeps three architectures × random-effect basis × by-variable × "
            "AR.start mode, computes the kitchen-sink metrics per variant, "
            "and runs k-fold CV on the Pareto front for a fair cross-arch "
            "ranking. Outputs under <code>&lt;output_dir&gt;/ablation/</code>."
        )
        intro.setWordWrap(True)
        intro.setStyleSheet("padding: 4px; color: #aaaaaa;")
        main.addWidget(intro)

        # Top: axes + run controls (collapsible via splitter)
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(self._build_axes_panel())
        splitter.addWidget(self._build_results_tabs())
        splitter.addWidget(self._build_log_panel())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 3)
        splitter.setStretchFactor(2, 1)
        main.addWidget(splitter, stretch=1)

    def _build_axes_panel(self) -> QWidget:
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMaximumHeight(290)
        inner = QWidget()
        scroll.setWidget(inner)
        vb = QVBoxLayout(inner)

        # Per-architecture axes
        self._cb_arch    = {}   # arch -> {axis_name -> {level -> QCheckBox}}
        self._arch_enabled = {} # arch -> QCheckBox (master toggle)

        for arch in _ARCHITECTURES:
            box = QGroupBox(_ARCH_LABELS[arch])
            bl = QVBoxLayout(box)

            arch_toggle = QCheckBox("Run this architecture")
            arch_toggle.setChecked(True)
            arch_toggle.setStyleSheet("font-weight: bold;")
            bl.addWidget(arch_toggle)
            self._arch_enabled[arch] = arch_toggle

            self._cb_arch[arch] = {}

            # RE basis row
            re_row = QHBoxLayout()
            re_row.addWidget(QLabel("RE basis:"))
            self._cb_arch[arch]["re"] = {}
            for b in _RE_BASES:
                cb = QCheckBox(b)
                cb.setChecked(True)
                cb.setToolTip(_RE_LABELS[b])
                re_row.addWidget(cb)
                self._cb_arch[arch]["re"][b] = cb
            re_row.addStretch()
            bl.addLayout(re_row)

            # By-variable row (only B and C)
            if arch in ("B", "C"):
                by_row = QHBoxLayout()
                by_row.addWidget(QLabel("By-var:"))
                self._cb_arch[arch]["by"] = {}
                for v in _BY_VARIABLES:
                    cb = QCheckBox(v)
                    cb.setChecked(True)
                    cb.setToolTip(_BY_LABELS[v])
                    by_row.addWidget(cb)
                    self._cb_arch[arch]["by"][v] = cb
                by_row.addStretch()
                bl.addLayout(by_row)

                # AR.start row (B and C only)
                ar_row = QHBoxLayout()
                ar_row.addWidget(QLabel("AR.start:"))
                self._cb_arch[arch]["ar"] = {}
                for a in _AR_START_MODES:
                    cb = QCheckBox(a)
                    cb.setChecked(True)
                    cb.setToolTip(_AR_LABELS[a])
                    ar_row.addWidget(cb)
                    self._cb_arch[arch]["ar"][a] = cb
                ar_row.addStretch()
                bl.addLayout(ar_row)
            else:
                # Arch A: by_var fixed to cond_combo, AR.start fixed to animal
                self._cb_arch[arch]["by"] = {"cond_combo": None}
                self._cb_arch[arch]["ar"] = {"animal":     None}
                note = QLabel(
                    "<i>Arch A: by-var fixed to <b>cond_combo</b>, "
                    "AR.start fixed to <b>animal</b> (single-phase models have "
                    "no internal phase boundary)</i>"
                )
                note.setStyleSheet("color: #888888; padding-left: 6px;")
                bl.addWidget(note)

            vb.addWidget(box)

        # Run controls
        ctrl_box = QGroupBox("Run controls")
        cbl = QHBoxLayout(ctrl_box)

        self._import_btn = QPushButton("Import existing results from output dir")
        self._import_btn.setMinimumHeight(36)
        self._import_btn.setToolTip(
            "Scan <output_dir>/ablation/ for prior run_metrics.json and "
            "cv_metrics.json files. Populates the Stage 1 table + Pareto "
            "fronts + (optionally) Stage 3 CV ranking without re-fitting. "
            "Useful for resuming a completed sweep or viewing past results."
        )
        self._import_btn.clicked.connect(self._on_import_existing)

        self._run_stage1_btn = QPushButton("Run Stage 1 (cheap metrics + adequacy)")
        self._run_stage1_btn.setMinimumHeight(36)
        self._run_stage1_btn.setStyleSheet("QPushButton { font-weight: bold; }")
        self._run_stage1_btn.clicked.connect(self._on_run_stage1)

        self._run_stage3_btn = QPushButton("Run Stage 3 CV on Pareto front")
        self._run_stage3_btn.setMinimumHeight(36)
        self._run_stage3_btn.setEnabled(False)
        self._run_stage3_btn.clicked.connect(self._on_run_stage3)

        self._stop_btn = QPushButton("Stop")
        self._stop_btn.setMinimumHeight(36)
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._on_stop)

        self._status_label = QLabel("idle")
        self._status_label.setStyleSheet("color: #aaaaaa;")

        cbl.addWidget(self._import_btn)
        cbl.addWidget(self._run_stage1_btn)
        cbl.addWidget(self._run_stage3_btn)
        cbl.addWidget(self._stop_btn)
        cbl.addWidget(self._status_label, stretch=1)

        vb.addWidget(ctrl_box)
        return scroll

    def _build_results_tabs(self) -> QWidget:
        self._tabs = QTabWidget()

        # 1. Stage 1 metrics matrix
        self._metrics_table = QTableWidget()
        self._metrics_table.setAlternatingRowColors(True)
        self._metrics_table.setSortingEnabled(True)
        self._metrics_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self._metrics_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self._tabs.addTab(self._metrics_table, "Stage 1 metrics")

        # 2. Random-eats-fixed bar chart
        self._eats_canvas = FigureCanvas(Figure(figsize=(8, 5), tight_layout=True))
        self._eats_canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._tabs.addTab(self._eats_canvas, "Random-eats-fixed")

        # 3. Cross-arch comparison
        self._cross_canvas = FigureCanvas(Figure(figsize=(8, 5), tight_layout=True))
        self._cross_canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._tabs.addTab(self._cross_canvas, "Cross-arch comparison")

        # 4. Pareto front
        self._pareto_canvas = FigureCanvas(Figure(figsize=(8, 5), tight_layout=True))
        self._pareto_canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._tabs.addTab(self._pareto_canvas, "Pareto front")

        # 5. Stage 3 CV ranking
        self._cv_table = QTableWidget()
        self._cv_table.setAlternatingRowColors(True)
        self._cv_table.setSortingEnabled(True)
        self._cv_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self._tabs.addTab(self._cv_table, "Stage 3 CV ranking")

        # 6. Per-contrast cross-variant table with RE-sensitivity
        self._contrast_table = QTableWidget()
        self._contrast_table.setAlternatingRowColors(True)
        self._contrast_table.setSortingEnabled(True)
        self._contrast_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self._tabs.addTab(self._contrast_table, "Per-contrast comparison")

        return self._tabs

    def _build_log_panel(self) -> QWidget:
        box = QGroupBox("Combined R log")
        bl = QVBoxLayout(box)
        self._log = QTextEdit()
        self._log.setReadOnly(True)
        self._log.setFont(QFont("Courier New", 9))
        bl.addWidget(self._log)
        return box

    # ── Run-plan generation ──────────────────────────────────────────────

    def _generate_stage1_plans(self) -> list:
        """Returns list of run-spec dicts, one per Stage 1 fit."""
        plans = []
        bam = self._bam_widget
        if bam is None or bam._df is None or bam._df.empty:
            return plans

        groups = sorted(int(g) for g in bam._df["Group"].dropna().unique())

        for arch in _ARCHITECTURES:
            if not self._arch_enabled[arch].isChecked():
                continue
            re_bases = [b for b in _RE_BASES if self._cb_arch[arch]["re"][b].isChecked()]
            if not re_bases:
                continue

            if arch == "A":
                # 5 RE × n_groups separate fits (one per phase group)
                for re in re_bases:
                    for g in groups:
                        plans.append({
                            "arch": arch, "re": re, "by": "cond_combo",
                            "ar":   "animal", "group_index": g,
                        })
            else:
                by_vars = [v for v in _BY_VARIABLES
                           if self._cb_arch[arch]["by"][v] is not None
                           and self._cb_arch[arch]["by"][v].isChecked()]
                ar_modes = [a for a in _AR_START_MODES
                            if self._cb_arch[arch]["ar"][a] is not None
                            and self._cb_arch[arch]["ar"][a].isChecked()]
                for re in re_bases:
                    for by in by_vars:
                        for ar in ar_modes:
                            plans.append({
                                "arch": arch, "re": re, "by": by,
                                "ar":   ar,   "group_index": None,
                            })
        return plans

    @staticmethod
    def _plan_run_id(p) -> str:
        if p["arch"] == "A":
            return f"A_{p['re']}_g{p['group_index']}"
        return f"{p['arch']}_{p['re']}_{p['by']}_{p['ar']}"

    @staticmethod
    def _plan_variant_id(p) -> str:
        """Variant ID groups Arch A's per-phase fits under a single 'variant'."""
        if p["arch"] == "A":
            return f"A_{p['re']}"
        return f"{p['arch']}_{p['re']}_{p['by']}_{p['ar']}"

    # ── Stage 1 run loop ──────────────────────────────────────────────────

    def _on_run_stage1(self):
        if self._bam_widget is None:
            QMessageBox.warning(self, "Not Ready", "BAM widget reference missing.")
            return
        bam = self._bam_widget
        if bam._df is None or bam._df.empty:
            QMessageBox.warning(self, "Not Ready",
                "Load data in the Data Loading tab before running the ablation.")
            return
        if not bam._output_dir:
            QMessageBox.warning(self, "Not Ready",
                "Set an output directory in the Data Loading tab before running.")
            return

        plans = self._generate_stage1_plans()
        if not plans:
            QMessageBox.warning(self, "Nothing to Run",
                "No variants selected. Enable at least one architecture and "
                "axis combination.")
            return

        self._stage1_queue       = plans
        self._stage1_total       = len(plans)
        self._stage1_metrics     = {}
        self._stage1_master      = {}
        self._stage1_fail_reason = {}
        self._pareto_candidates  = {}
        self._pareto_union       = []
        self._stage3_results     = {}
        self._user_stopped       = False
        self._current_stage      = "stage1"
        self._log.clear()
        self._metrics_table.setRowCount(0)
        self._metrics_table.setColumnCount(0)
        self._cv_table.setRowCount(0)
        self._cv_table.setColumnCount(0)
        self._contrast_table.setRowCount(0)
        self._contrast_table.setColumnCount(0)
        for c in (self._eats_canvas, self._cross_canvas, self._pareto_canvas):
            c.figure.clear()
            c.draw_idle()

        self._run_stage1_btn.setEnabled(False)
        self._run_stage3_btn.setEnabled(False)
        self._stop_btn.setEnabled(True)
        self._append_log(f"--- Stage 1 starting: {self._stage1_total} fit(s) ---")
        self._run_next_stage1()

    def _run_next_stage1(self):
        if self._user_stopped or not self._stage1_queue:
            self._finish_stage1()
            return
        plan = self._stage1_queue.pop(0)
        run_id = self._plan_run_id(plan)
        self._current_run_id = run_id

        done = self._stage1_total - len(self._stage1_queue)
        self._status_label.setText(
            f"Stage 1: {done}/{self._stage1_total} — running {run_id}"
        )
        self._append_log(f"\n========== Stage 1 [{done}/{self._stage1_total}]: {run_id} ==========")

        subdir = self._stage1_subdir(plan)
        os.makedirs(subdir, exist_ok=True)

        csv_path = self._export_csv_to(subdir)
        if not csv_path:
            self._stage1_fail_reason[run_id] = "CSV export failed"
            self._stage1_done.emit(run_id, False)
            return
        self._write_sidecars(subdir, skip_posterior=True)

        rscript = BamWidget._find_rscript()
        if not rscript:
            QMessageBox.critical(self, "Rscript Not Found",
                "Could not locate Rscript. Install R or set up the bundled R runtime.")
            self._stage1_fail_reason[run_id] = "Rscript not found"
            self._stage1_done.emit(run_id, False)
            return

        env = BamWidget._r_env(rscript)
        env["ABL_ARCH"]            = plan["arch"]
        env["ABL_BY_VARIABLE"]     = plan["by"]
        env["ABL_AR_START_MODE"]   = plan["ar"]
        env["ABL_RHO_FIXED"]       = _ABLATION_RHO
        env["ABL_SKIP_POSTERIOR"]  = "TRUE"
        env["RANDOM_BASIS"]        = plan["re"]
        if plan["arch"] == "A":
            env["ABL_GROUP_INDEX"] = str(plan["group_index"])

        cmd = self._build_cmd(rscript, csv_path, subdir)
        self._launch(cmd, env, on_done=self._stage1_done)

    def _on_stage1_done(self, run_id, success):
        # Find the matching plan to know whether it was an Arch A per-group fit
        plan = self._infer_plan_from_run_id(run_id)
        subdir = self._stage1_subdir(plan) if plan else None

        if success and subdir:
            metrics_path = os.path.join(subdir, "run_metrics.json")
            if os.path.isfile(metrics_path):
                try:
                    with open(metrics_path, "r", encoding="utf-8") as f:
                        self._stage1_metrics[run_id] = json.load(f)
                    self._append_log(f"Loaded run_metrics.json for {run_id}")
                except Exception as e:
                    self._stage1_fail_reason[run_id] = f"Metrics read error: {e}"
                    self._append_log(f"Could not read run_metrics.json: {e}")
            else:
                self._stage1_fail_reason[run_id] = "run_metrics.json missing"
                self._append_log(f"[Missing] {metrics_path}")
            master_path = os.path.join(subdir, "master_results.csv")
            if os.path.isfile(master_path):
                try:
                    self._stage1_master[run_id] = pd.read_csv(master_path)
                except Exception as e:
                    self._append_log(f"Could not read master_results.csv: {e}")
        elif not success:
            self._stage1_fail_reason.setdefault(run_id, "R exited non-zero")
            self._append_log(f"--- Variant {run_id} FAILED ---")
        self._run_next_stage1()

    def _finish_stage1(self):
        self._run_stage1_btn.setEnabled(True)
        self._stop_btn.setEnabled(False)
        if self._user_stopped:
            self._status_label.setText("stopped")
            return
        n_loaded = len(self._stage1_metrics)
        self._status_label.setText(
            f"Stage 1 done — loaded {n_loaded}/{self._stage1_total} variants"
        )
        self._process_stage1_results()

    def _process_stage1_results(self):
        """Aggregate, filter, populate displays, compute Pareto. Called from
        both Stage 1 completion and from import-existing-results."""
        variants = self._aggregate_variants()
        adequacy_flags = self._apply_adequacy_filters(variants)
        self._populate_metrics_table(variants, adequacy_flags)
        self._populate_random_eats_chart(variants, adequacy_flags)
        self._populate_cross_arch_chart(variants, adequacy_flags)

        # Per-arch Pareto over surviving variants
        self._pareto_candidates = {}
        self._pareto_union = []
        for arch in _ARCHITECTURES:
            arch_vars = [v for v in variants
                         if v["arch"] == arch
                         and not adequacy_flags.get(v["variant_id"], {}).get("failed")]
            front_idx = _pareto_front(arch_vars, _PARETO_METRICS)
            front_ids = [arch_vars[i]["variant_id"] for i in front_idx]
            self._pareto_candidates[arch] = front_ids
            self._pareto_union.extend(front_ids)
        self._populate_pareto_chart(variants, adequacy_flags)
        self._populate_contrast_table(variants)

        if self._pareto_union:
            self._run_stage3_btn.setEnabled(True)
            self._append_log(
                f"\nPareto fronts: " +
                ", ".join(f"{a}={len(v)}" for a, v in self._pareto_candidates.items())
            )
        else:
            self._append_log("\nNo surviving Pareto candidates — Stage 3 disabled.")

    # ── Stage 1 helper: aggregate Arch A per-group fits into one variant ──

    def _aggregate_variants(self) -> list:
        """Build variant records from loaded run_metrics. For Arch A, combine
        per-group fits into one record with summed AIC/BIC/etc."""
        records = []
        # Bucket runs by variant_id
        buckets = {}
        for run_id, metrics in self._stage1_metrics.items():
            plan = self._infer_plan_from_run_id(run_id)
            if plan is None:
                continue
            vid = self._plan_variant_id(plan)
            buckets.setdefault(vid, []).append((run_id, plan, metrics))

        for vid, items in buckets.items():
            arch = items[0][1]["arch"]
            re = items[0][1]["re"]
            by = items[0][1]["by"]
            ar = items[0][1]["ar"]
            # Sum-style for Arch A; single-fit for B/C
            metrics_list = [m for _, _, m in items]
            agg = {
                "variant_id":  vid,
                "arch":        arch,
                "re":          re,
                "by":          by,
                "ar":          ar,
                "n_fits":      len(items),
                # Sums
                "AIC":         self._sum(metrics_list, "fit_quality", "AIC"),
                "BIC":         self._sum(metrics_list, "fit_quality", "BIC"),
                "logLik":      self._sum(metrics_list, "fit_quality", "logLik"),
                "fREML":       self._sum(metrics_list, "fit_quality", "fREML"),
                "total_edf":   self._sum(metrics_list, "fit_quality", "total_edf"),
                # Means (Arch A averages across phases; B/C single value)
                "dev_explained":    self._mean(metrics_list, "fit_quality", "dev_explained"),
                "r_sq_adj":         self._mean(metrics_list, "fit_quality", "r_sq_adj"),
                "tweedie_phi":      self._mean(metrics_list, "fit_quality", "tweedie_phi"),
                "ar1_residual_lag1": self._mean(metrics_list, "adequacy", "ar1_residual_lag1"),
                # Max (worst across phases for Arch A)
                "max_overall_concurvity":  self._max(metrics_list, "concurvity", "max_overall"),
                "max_random_vs_fixed":     self._max(metrics_list, "concurvity", "max_random_vs_fixed"),
                "edf_random_to_fixed":     self._max(metrics_list, "edf",        "random_to_fixed"),
                "random_sd_to_phi":        self._max(metrics_list, "variance_components", "random_sd_to_phi"),
                # k-index pass rate (mean across phase fits)
                "k_index_pass_rate":       self._mean(metrics_list, "adequacy", "k_index", "pass_rate"),
                # Convergence
                "convergence_warnings":    sum(
                    len(_safe_get(m, "adequacy", "convergence_warnings", default=[]) or [])
                    for m in metrics_list
                ),
            }
            # Contrast SE: mean across Test_Families × items
            agg["mean_contrast_SE"] = self._contrast_mean_se(metrics_list)
            records.append(agg)
        # Sort by arch, re, by, ar
        records.sort(key=lambda r: (r["arch"], r["re"], r["by"], r["ar"]))
        return records

    @staticmethod
    def _sum(ms, *keys):
        vals = [_safe_get(m, *keys) for m in ms]
        vals = [v for v in vals if isinstance(v, (int, float)) and not math.isnan(v)]
        return sum(vals) if vals else float("nan")

    @staticmethod
    def _mean(ms, *keys):
        vals = [_safe_get(m, *keys) for m in ms]
        vals = [v for v in vals if isinstance(v, (int, float)) and not math.isnan(v)]
        return (sum(vals) / len(vals)) if vals else float("nan")

    @staticmethod
    def _max(ms, *keys):
        vals = [_safe_get(m, *keys) for m in ms]
        vals = [v for v in vals if isinstance(v, (int, float)) and not math.isnan(v)]
        return max(vals) if vals else float("nan")

    @staticmethod
    def _contrast_mean_se(metrics_list):
        ses = []
        for m in metrics_list:
            cv = _safe_get(m, "contrast_variance", default={}) or {}
            if not isinstance(cv, dict):
                continue
            for fam, info in cv.items():
                if isinstance(info, dict):
                    se = info.get("mean_SE")
                    if isinstance(se, (int, float)) and not math.isnan(se):
                        ses.append(se)
        return (sum(ses) / len(ses)) if ses else float("nan")

    # ── Adequacy filters ─────────────────────────────────────────────────

    def _apply_adequacy_filters(self, variants):
        """Returns dict {variant_id: {"failed": bool, "reasons": [str]}}."""
        out = {}
        # Need min phi across surviving (non-NaN) variants for the relative phi check
        phis = [v["tweedie_phi"] for v in variants
                if isinstance(v["tweedie_phi"], (int, float)) and not math.isnan(v["tweedie_phi"])]
        min_phi = min(phis) if phis else float("nan")

        for v in variants:
            reasons = []
            if v["convergence_warnings"] and v["convergence_warnings"] > 0:
                reasons.append(f"{v['convergence_warnings']} convergence warning(s)")
            k_rate = v.get("k_index_pass_rate")
            if isinstance(k_rate, (int, float)) and not math.isnan(k_rate):
                if k_rate < _ADEQUACY["k_index_pass_rate_min"]:
                    reasons.append(f"k-index pass rate {k_rate:.0%} < 30%")
            oc = v.get("max_overall_concurvity")
            if isinstance(oc, (int, float)) and not math.isnan(oc):
                if oc > _ADEQUACY["max_overall_concurvity"]:
                    reasons.append(f"max concurvity {oc:.2f} > 0.97")
            rf = v.get("max_random_vs_fixed")
            if isinstance(rf, (int, float)) and not math.isnan(rf):
                if rf > _ADEQUACY["max_random_vs_fixed_conc"]:
                    reasons.append(f"random-vs-fixed concurvity {rf:.2f} > 0.95")
            ar = v.get("ar1_residual_lag1")
            if isinstance(ar, (int, float)) and not math.isnan(ar):
                if abs(ar) > _ADEQUACY["max_ar1_residual_lag1"]:
                    reasons.append(f"|AR(1) residual lag-1| {abs(ar):.2f} > 0.7")
            phi = v["tweedie_phi"]
            if (isinstance(phi, (int, float)) and not math.isnan(phi)
                and not math.isnan(min_phi) and min_phi > 0):
                if phi > _ADEQUACY["phi_relative_multiplier"] * min_phi:
                    reasons.append(
                        f"phi {phi:.3g} > {_ADEQUACY['phi_relative_multiplier']:.0f}× min"
                    )
            out[v["variant_id"]] = {"failed": bool(reasons), "reasons": reasons}
        return out

    # ── Table population ─────────────────────────────────────────────────

    def _populate_metrics_table(self, variants, adequacy_flags):
        cols = [
            "variant_id", "arch", "re", "by", "ar",
            "AIC", "BIC", "fREML", "dev_explained", "r_sq_adj",
            "tweedie_phi", "total_edf",
            "k_index_pass_rate", "ar1_residual_lag1",
            "max_overall_concurvity", "max_random_vs_fixed",
            "edf_random_to_fixed", "random_sd_to_phi",
            "mean_contrast_SE", "convergence_warnings",
            "ADEQUATE", "ADEQUACY_NOTE",
        ]
        self._metrics_table.setSortingEnabled(False)
        self._metrics_table.setColumnCount(len(cols))
        self._metrics_table.setHorizontalHeaderLabels(cols)
        self._metrics_table.setRowCount(len(variants))

        for r, v in enumerate(variants):
            flag = adequacy_flags.get(v["variant_id"], {})
            failed = flag.get("failed", False)
            note = "; ".join(flag.get("reasons", []))
            row_color = QColor("#5a2222") if failed else None
            for c, col in enumerate(cols):
                if col == "ADEQUATE":
                    text = "FAIL" if failed else "OK"
                    item = QTableWidgetItem(text)
                    item.setForeground(QColor("#ff8888") if failed else QColor("#88dd88"))
                elif col == "ADEQUACY_NOTE":
                    item = QTableWidgetItem(note or "—")
                else:
                    val = v.get(col)
                    if isinstance(val, (int, float)):
                        item = QTableWidgetItem(_fmt(val))
                    else:
                        item = QTableWidgetItem(str(val) if val is not None else "—")
                if row_color is not None:
                    item.setBackground(row_color)
                self._metrics_table.setItem(r, c, item)
        self._metrics_table.setSortingEnabled(True)
        self._metrics_table.resizeColumnsToContents()

        # Also dump failed variants explicitly to the log so the user has a
        # text record of what was crossed off and why.
        failed_lines = []
        for v in variants:
            f = adequacy_flags.get(v["variant_id"], {})
            if f.get("failed"):
                failed_lines.append(f"  {v['variant_id']}: {'; '.join(f['reasons'])}")
        if failed_lines:
            self._append_log("\n--- Variants failed adequacy filters (excluded from Pareto) ---")
            for l in failed_lines:
                self._append_log(l)

    def _populate_random_eats_chart(self, variants, adequacy_flags):
        def _num(v):
            # Preserve genuine zero values; only coerce None / non-finite to NaN.
            if v is None:
                return float("nan")
            try:
                vf = float(v)
                return vf if math.isfinite(vf) else float("nan")
            except (TypeError, ValueError):
                return float("nan")

        fig = self._eats_canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        labels = [v["variant_id"] for v in variants]
        c_rf  = [_num(v.get("max_random_vs_fixed")) for v in variants]
        ratio = [_num(v.get("edf_random_to_fixed")) for v in variants]
        vcomp = [_num(v.get("random_sd_to_phi"))    for v in variants]
        n = len(variants)
        x = list(range(n))
        w = 0.27
        ax.bar([xi - w for xi in x], c_rf, width=w, label="concurvity(R,F)")
        ax.bar(x,                     ratio, width=w, label="EDF random/fixed")
        ax.bar([xi + w for xi in x], vcomp, width=w, label="random_sd / √φ")
        ax.axhline(_ADEQUACY["max_random_vs_fixed_conc"], color="red", ls="--",
                   lw=0.8, label="concurvity threshold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=80, fontsize=7)
        ax.set_title("Random-eats-fixed diagnostics per variant")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(True, axis="y", alpha=0.3)
        self._eats_canvas.draw_idle()

    def _populate_cross_arch_chart(self, variants, adequacy_flags):
        fig = self._cross_canvas.figure
        fig.clear()
        # For each arch, take the BEST surviving AIC as the representative
        per_arch = {}
        for v in variants:
            if adequacy_flags.get(v["variant_id"], {}).get("failed"):
                continue
            arch = v["arch"]
            if (arch not in per_arch
                or (isinstance(v.get("AIC"), (int, float))
                    and v["AIC"] < per_arch[arch]["AIC"])):
                per_arch[arch] = v
        archs = sorted(per_arch.keys())
        if not archs:
            ax = fig.add_subplot(111)
            ax.set_title("No surviving variants for cross-arch comparison")
            self._cross_canvas.draw_idle()
            return

        # Metrics span very different scales (AIC ~ millions, phi ~ tens,
        # dev_explained ~ [0,1]). Use one subplot per metric so each gets
        # its own y-axis and the small ones don't get crushed.
        metrics = ("AIC", "BIC", "dev_explained", "r_sq_adj",
                   "tweedie_phi", "mean_contrast_SE")
        ncols = 3
        nrows = (len(metrics) + ncols - 1) // ncols
        axes = fig.subplots(nrows, ncols, sharex=False)
        axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

        x = list(range(len(archs)))
        for j, m in enumerate(metrics):
            ax = axes[j]
            vals = []
            for arch in archs:
                raw = per_arch[arch].get(m)
                try:
                    vals.append(float(raw) if raw is not None and math.isfinite(float(raw))
                                else float("nan"))
                except (TypeError, ValueError):
                    vals.append(float("nan"))
            ax.bar(x, vals, width=0.6,
                   color=[("#4488dd","#dd4488","#44dd88")[
                       ("A","B","C").index(a) if a in ("A","B","C") else 0
                   ] for a in archs])
            ax.set_xticks(x)
            ax.set_xticklabels([f"Arch {a}" for a in archs], fontsize=8)
            ax.set_title(m, fontsize=10)
            ax.grid(True, axis="y", alpha=0.3)
            # Annotate each bar with its value
            for xi, vi in zip(x, vals):
                if not math.isnan(vi):
                    ax.text(xi, vi, _fmt(vi, 3),
                            ha="center", va="bottom", fontsize=7)
        # Hide any unused subplot cells
        for k in range(len(metrics), len(axes)):
            axes[k].axis("off")
        fig.suptitle("Cross-architecture comparison (best-AIC variant per arch)",
                     fontsize=11)
        self._cross_canvas.draw_idle()

    def _populate_pareto_chart(self, variants, adequacy_flags):
        fig = self._pareto_canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        # X = AIC (minimize), Y = mean_contrast_SE (minimize) — useful 2D view
        arch_colors = {"A": "#4488dd", "B": "#dd4488", "C": "#44dd88"}
        front_ids = set(self._pareto_union)
        for v in variants:
            fid = v["variant_id"]
            adequate = not adequacy_flags.get(fid, {}).get("failed")
            color = arch_colors.get(v["arch"], "#888888")
            x = v.get("AIC")
            y = v.get("mean_contrast_SE")
            if not (isinstance(x, (int, float)) and isinstance(y, (int, float))):
                continue
            on_front = fid in front_ids
            ax.scatter(
                [x], [y],
                c=[color],
                s=120 if on_front else 40,
                edgecolors="white" if on_front else "none",
                linewidths=1.5 if on_front else 0,
                alpha=1.0 if adequate else 0.25,
                marker="o" if adequate else "x",
            )
            if on_front:
                ax.annotate(fid, (x, y), fontsize=7,
                            textcoords="offset points", xytext=(4, 4))
        ax.set_xlabel("AIC (lower = better)")
        ax.set_ylabel("mean contrast SE (lower = better)")
        ax.set_title("Pareto front overlay — large dots = on Pareto front per arch")
        ax.grid(True, alpha=0.3)
        # Legend hack
        for arch, col in arch_colors.items():
            ax.scatter([], [], c=col, label=f"Arch {arch}", s=60)
        ax.legend(fontsize=9, loc="best")
        self._pareto_canvas.draw_idle()

    def _populate_contrast_table(self, variants):
        # Build per-contrast cross-variant table from stored master_results.csv
        if not self._stage1_master:
            return
        # Merge by (Test_Family, Group, Split_By, Tested_Level)
        merge_keys = ["Test_Family", "Group", "Split_By", "Tested_Level"]
        merged = None
        # Pull master_results per VARIANT (already aggregated per Arch A by run_id reuse:
        # for A we have separate masters per group — concat them per variant_id)
        per_variant_master = {}
        for run_id, df in self._stage1_master.items():
            plan = self._infer_plan_from_run_id(run_id)
            if plan is None:
                continue
            vid = self._plan_variant_id(plan)
            if vid in per_variant_master:
                per_variant_master[vid] = pd.concat(
                    [per_variant_master[vid], df], ignore_index=True
                )
            else:
                per_variant_master[vid] = df.copy()

        # Use the variant ordering from `variants`
        ordered_vids = [v["variant_id"] for v in variants]
        for vid in ordered_vids:
            df = per_variant_master.get(vid)
            if df is None:
                continue
            sub = df[merge_keys + ["estimate", "SE", "raw_pvalue"]].copy()
            sub = sub.rename(columns={
                "estimate":   f"est_{vid}",
                "SE":         f"SE_{vid}",
                "raw_pvalue": f"p_{vid}",
            })
            merged = sub if merged is None else merged.merge(
                sub, on=merge_keys, how="outer"
            )
        if merged is None or merged.empty:
            return

        # RE-sensitivity column: range of estimate across variants for each contrast
        est_cols = [f"est_{vid}" for vid in ordered_vids if f"est_{vid}" in merged.columns]
        if est_cols:
            merged["re_sensitivity_range"] = (
                merged[est_cols].max(axis=1) - merged[est_cols].min(axis=1)
            )
        merged = merged.sort_values(merge_keys).reset_index(drop=True)

        cols = list(merge_keys)
        for vid in ordered_vids:
            cols += [f"p_{vid}", f"SE_{vid}", f"est_{vid}"]
        if "re_sensitivity_range" in merged.columns:
            cols.append("re_sensitivity_range")
        cols = [c for c in cols if c in merged.columns]

        self._contrast_table.setSortingEnabled(False)
        self._contrast_table.setColumnCount(len(cols))
        self._contrast_table.setHorizontalHeaderLabels(cols)
        self._contrast_table.setRowCount(len(merged))
        for r in range(len(merged)):
            row = merged.iloc[r]
            for c, col in enumerate(cols):
                val = row[col]
                if pd.isna(val):
                    item = QTableWidgetItem("—")
                elif col.startswith("p_"):
                    item = QTableWidgetItem(_fmt(val))
                    item.setForeground(QColor("#5fb35f") if val < 0.05 else QColor("#cc6666"))
                elif col.startswith("SE_") or col.startswith("est_"):
                    item = QTableWidgetItem(_fmt(val, 3))
                elif col == "Group":
                    item = QTableWidgetItem(str(int(val)))
                elif col == "re_sensitivity_range":
                    item = QTableWidgetItem(_fmt(val, 3))
                else:
                    item = QTableWidgetItem(str(val))
                self._contrast_table.setItem(r, c, item)
        self._contrast_table.setSortingEnabled(True)
        self._contrast_table.resizeColumnsToContents()

    # ── Import existing results from output dir (skip Stage 1 re-fit) ─────

    def _on_import_existing(self):
        if self._bam_widget is None or not self._bam_widget._output_dir:
            QMessageBox.warning(self, "Not Ready",
                "Set an output directory in the Data Loading tab first.")
            return
        if self._process is not None and self._process.poll() is None:
            QMessageBox.warning(self, "Busy",
                "A run is currently in progress — stop it before importing.")
            return

        output_dir = self._bam_widget._output_dir
        stage1_dir = os.path.join(output_dir, "ablation", "stage1")
        if not os.path.isdir(stage1_dir):
            QMessageBox.information(self, "No Existing Results",
                f"No prior ablation results found at:\n{stage1_dir}\n\n"
                "Run Stage 1 to start a new ablation.")
            return

        # Reset state (mirrors _on_run_stage1's reset block)
        self._stage1_metrics     = {}
        self._stage1_master      = {}
        self._stage1_fail_reason = {}
        self._pareto_candidates  = {}
        self._pareto_union       = []
        self._stage3_results     = {}
        self._user_stopped       = False
        self._metrics_table.setRowCount(0); self._metrics_table.setColumnCount(0)
        self._cv_table.setRowCount(0);      self._cv_table.setColumnCount(0)
        self._contrast_table.setRowCount(0); self._contrast_table.setColumnCount(0)
        for c in (self._eats_canvas, self._cross_canvas, self._pareto_canvas):
            c.figure.clear(); c.draw_idle()
        self._log.clear()

        self._append_log(f"--- Importing Stage 1 results from: {stage1_dir} ---")
        n_metrics, n_master = self._scan_stage1_dir(stage1_dir)

        if n_metrics == 0:
            QMessageBox.information(self, "No Existing Results",
                f"No run_metrics.json files found under:\n{stage1_dir}")
            self._status_label.setText("idle")
            return

        self._stage1_total = n_metrics
        self._append_log(
            f"Loaded {n_metrics} run_metrics.json + {n_master} master_results.csv file(s)."
        )

        # Completeness check against the current UI axis selection.
        completeness = self._check_stage1_completeness(stage1_dir)
        if completeness["missing_runs"] or completeness["missing_metrics"] or completeness["missing_masters"]:
            self._show_completeness_report(completeness)

        self._status_label.setText(
            f"Imported {n_metrics} Stage 1 variant(s) — running Stage 2 analysis…"
        )

        # Stage 2 (Pareto + displays) + Stage 3 button enablement
        self._process_stage1_results()
        self._status_label.setText(
            f"Imported {n_metrics} Stage 1 variant(s); "
            f"{len(self._pareto_union)} Pareto candidate(s)."
        )

        # Optionally pick up any existing Stage 3 CV results in the output dir
        stage3_dir = os.path.join(output_dir, "ablation", "stage3_cv")
        n_cv = self._scan_stage3_dir(stage3_dir) if os.path.isdir(stage3_dir) else 0
        if n_cv > 0:
            self._append_log(f"Loaded {n_cv} Stage 3 CV fold result(s) from {stage3_dir}.")
            self._populate_cv_table()
            self._status_label.setText(
                self._status_label.text() + f" Stage 3 imported ({n_cv} folds)."
            )

    def _check_stage1_completeness(self, stage1_dir: str) -> dict:
        """Compare expected Stage 1 runs (from current checkbox selection) to
        what was actually found on disk under stage1_dir. Returns the
        partition: present-and-complete / subdir-missing-entirely /
        subdir-present-but-metrics-missing / subdir-present-but-master-missing.
        """
        expected_plans = self._generate_stage1_plans()
        expected_run_ids = [self._plan_run_id(p) for p in expected_plans]

        present_subdirs = set()
        if os.path.isdir(stage1_dir):
            for entry in os.listdir(stage1_dir):
                if os.path.isdir(os.path.join(stage1_dir, entry)):
                    present_subdirs.add(entry)

        complete         = []
        missing_runs     = []   # no subdir at all
        missing_metrics  = []   # subdir present, run_metrics.json absent
        missing_masters  = []   # subdir present, master_results.csv absent

        for rid in expected_run_ids:
            if rid not in present_subdirs:
                missing_runs.append(rid)
                continue
            subdir = os.path.join(stage1_dir, rid)
            has_metrics = os.path.isfile(os.path.join(subdir, "run_metrics.json"))
            has_master  = os.path.isfile(os.path.join(subdir, "master_results.csv"))
            if not has_metrics:
                missing_metrics.append(rid)
            if not has_master:
                missing_masters.append(rid)
            if has_metrics and has_master:
                complete.append(rid)

        # Also flag unexpected subdirs (likely from a prior axis selection)
        unexpected = sorted(present_subdirs - set(expected_run_ids))

        return {
            "expected":         expected_run_ids,
            "complete":         complete,
            "missing_runs":     missing_runs,
            "missing_metrics":  missing_metrics,
            "missing_masters":  missing_masters,
            "unexpected":       unexpected,
        }

    def _show_completeness_report(self, completeness: dict):
        """Log the completeness report and surface a dialog if anything is
        missing relative to the current axis selection."""
        n_expected     = len(completeness["expected"])
        n_complete     = len(completeness["complete"])
        n_missing_runs = len(completeness["missing_runs"])
        n_no_metrics   = len(completeness["missing_metrics"])
        n_no_master    = len(completeness["missing_masters"])
        n_unexpected   = len(completeness["unexpected"])

        self._append_log(
            f"\nCompleteness vs current checkbox selection: "
            f"{n_complete}/{n_expected} complete, "
            f"{n_missing_runs} never run, "
            f"{n_no_metrics} missing run_metrics.json, "
            f"{n_no_master} missing master_results.csv."
        )
        if completeness["missing_runs"]:
            self._append_log("  Never-run variants (no subdir):")
            for rid in completeness["missing_runs"]:
                self._append_log(f"    {rid}")
        if completeness["missing_metrics"]:
            self._append_log("  Subdir present but run_metrics.json missing "
                             "(likely crashed mid-fit):")
            for rid in completeness["missing_metrics"]:
                self._append_log(f"    {rid}")
        if completeness["missing_masters"]:
            self._append_log("  Subdir present but master_results.csv missing:")
            for rid in completeness["missing_masters"]:
                self._append_log(f"    {rid}")
        if completeness["unexpected"]:
            self._append_log(
                "  Unexpected subdirs (present on disk but not in current "
                "checkbox selection — leftover from a prior run with "
                "different axes):"
            )
            for rid in completeness["unexpected"]:
                self._append_log(f"    {rid}")

        # Show a dialog with a brief summary so the user can decide whether
        # to proceed or re-run the missing variants
        any_missing = n_missing_runs + n_no_metrics + n_no_master
        if any_missing:
            QMessageBox.warning(
                self, "Incomplete Stage 1 import",
                f"Imported {n_complete} of {n_expected} expected variants "
                f"(based on the current checkbox selection).\n\n"
                f"  {n_missing_runs} variant(s) never run\n"
                f"  {n_no_metrics} variant(s) missing run_metrics.json\n"
                f"  {n_no_master} variant(s) missing master_results.csv\n"
                + (f"  {n_unexpected} extra subdir(s) from a prior axis selection\n"
                   if n_unexpected else "")
                + "\nStage 2 / 3 will use only the loaded variants. See the "
                  "log for the full list. Re-run Stage 1 if you want the "
                  "missing variants populated."
            )

    def _scan_stage1_dir(self, stage1_dir: str):
        """Scan stage1_dir for per-variant run_metrics.json + master_results.csv.
        Returns (n_metrics_loaded, n_master_loaded)."""
        n_metrics = 0
        n_master  = 0
        for entry in sorted(os.listdir(stage1_dir)):
            subdir = os.path.join(stage1_dir, entry)
            if not os.path.isdir(subdir):
                continue
            run_id = entry
            # Validate that the directory name parses as a known run plan
            if self._infer_plan_from_run_id(run_id) is None:
                self._append_log(f"  Skipping unrecognized subdir: {run_id}")
                continue
            metrics_path = os.path.join(subdir, "run_metrics.json")
            master_path  = os.path.join(subdir, "master_results.csv")
            if os.path.isfile(metrics_path):
                try:
                    with open(metrics_path, "r", encoding="utf-8") as f:
                        self._stage1_metrics[run_id] = json.load(f)
                    n_metrics += 1
                except Exception as e:
                    self._append_log(f"  Could not parse {metrics_path}: {e}")
            if os.path.isfile(master_path):
                try:
                    self._stage1_master[run_id] = pd.read_csv(master_path)
                    n_master += 1
                except Exception as e:
                    self._append_log(f"  Could not read {master_path}: {e}")
        return n_metrics, n_master

    def _scan_stage3_dir(self, stage3_dir: str) -> int:
        """Scan stage3_cv subdirs for cv_metrics.json files. Returns count loaded."""
        n_loaded = 0
        for entry in sorted(os.listdir(stage3_dir)):
            subdir = os.path.join(stage3_dir, entry)
            if not os.path.isdir(subdir):
                continue
            cv_path = os.path.join(subdir, "cv_metrics.json")
            if not os.path.isfile(cv_path):
                continue
            plan = self._infer_stage3_plan_from_run_id(entry)
            if plan is None:
                continue
            try:
                with open(cv_path, "r", encoding="utf-8") as f:
                    cv_metrics = json.load(f)
                dev_raw = cv_metrics.get("holdout_deviance")
                dev = float(dev_raw) if isinstance(dev_raw, (int, float)) else float("nan")
            except Exception as e:
                self._append_log(f"  Could not read {cv_path}: {e}")
                continue
            vid = plan["variant_id"]
            self._stage3_results.setdefault(vid, []).append({
                "fold":  plan["fold"],
                "group": plan.get("group_index"),
                "dev":   dev,
            })
            n_loaded += 1
        return n_loaded

    # ── Stage 3 CV ───────────────────────────────────────────────────────

    def _on_run_stage3(self):
        if not self._pareto_union:
            return
        # Build queue: per candidate × _N_CV_FOLDS folds × (n_groups for Arch A)
        bam = self._bam_widget
        groups = sorted(int(g) for g in bam._df["Group"].dropna().unique())
        plans = []
        for vid in self._pareto_union:
            # Recover the canonical plan(s) for this variant
            arch = vid.split("_")[0]
            if arch == "A":
                # vid format: "A_<re>"
                re = vid.split("_", 1)[1]
                for g in groups:
                    for fold in range(_N_CV_FOLDS):
                        plans.append({
                            "variant_id":  vid,
                            "arch":        "A", "re": re, "by": "cond_combo",
                            "ar":          "animal",
                            "group_index": g, "fold": fold,
                        })
            else:
                # vid format: "<arch>_<re>_<by>_<ar>"
                parts = vid.split("_")
                # ar is the last token; by is the second-to-last; re is parts[1]
                # by might itself contain underscores (cond_combo, group_factor)
                # — handle robustly by checking suffixes:
                arch = parts[0]
                ar = parts[-1]
                # by candidates
                for cand_by in _BY_VARIABLES:
                    suffix = f"_{cand_by}_{ar}"
                    if vid.endswith(suffix):
                        by = cand_by
                        re = vid[len(arch) + 1 : -len(suffix)]
                        break
                else:
                    self._append_log(f"Could not parse variant id {vid} — skipping")
                    continue
                for fold in range(_N_CV_FOLDS):
                    plans.append({
                        "variant_id":  vid,
                        "arch":        arch, "re": re, "by": by, "ar": ar,
                        "group_index": None, "fold": fold,
                    })

        self._stage3_queue   = plans
        self._stage3_total   = len(plans)
        self._stage3_results = {}
        self._user_stopped   = False
        self._current_stage  = "stage3"
        self._run_stage1_btn.setEnabled(False)
        self._run_stage3_btn.setEnabled(False)
        self._stop_btn.setEnabled(True)
        self._append_log(f"\n--- Stage 3 CV starting: {self._stage3_total} fit(s) ---")
        self._run_next_stage3()

    def _run_next_stage3(self):
        if self._user_stopped or not self._stage3_queue:
            self._finish_stage3()
            return
        plan = self._stage3_queue.pop(0)
        run_id = self._stage3_run_id(plan)
        self._current_run_id = run_id

        done = self._stage3_total - len(self._stage3_queue)
        self._status_label.setText(
            f"Stage 3 CV: {done}/{self._stage3_total} — {run_id}"
        )
        self._append_log(f"\n========== Stage 3 [{done}/{self._stage3_total}]: {run_id} ==========")

        subdir = self._stage3_subdir(plan)
        os.makedirs(subdir, exist_ok=True)
        csv_path = self._export_csv_to(subdir)
        if not csv_path:
            self._append_log("CSV export failed; skipping fold")
            self._stage3_done.emit(run_id, False)
            return
        self._write_sidecars(subdir, skip_posterior=True)

        rscript = BamWidget._find_rscript()
        if not rscript:
            QMessageBox.critical(self, "Rscript Not Found",
                "Could not locate Rscript.")
            self._stage3_done.emit(run_id, False)
            return

        env = BamWidget._r_env(rscript)
        env["ABL_ARCH"]            = plan["arch"]
        env["ABL_BY_VARIABLE"]     = plan["by"]
        env["ABL_AR_START_MODE"]   = plan["ar"]
        env["ABL_RHO_FIXED"]       = _ABLATION_RHO
        env["ABL_SKIP_POSTERIOR"]  = "TRUE"
        env["RANDOM_BASIS"]        = plan["re"]
        env["ABL_CV_FOLD"]         = str(plan["fold"])
        if plan["arch"] == "A":
            env["ABL_GROUP_INDEX"] = str(plan["group_index"])

        cmd = self._build_cmd(rscript, csv_path, subdir)
        self._launch(cmd, env, on_done=self._stage3_done)

    def _on_stage3_done(self, run_id, success):
        # Parse the fold's cv_metrics.json from the subdir
        plan = self._infer_stage3_plan_from_run_id(run_id)
        if plan is None:
            self._run_next_stage3()
            return
        subdir = self._stage3_subdir(plan)
        cv_path = os.path.join(subdir, "cv_metrics.json")
        dev = float("nan")
        if success and os.path.isfile(cv_path):
            try:
                with open(cv_path, "r", encoding="utf-8") as f:
                    cv_metrics = json.load(f)
                dev_raw = cv_metrics.get("holdout_deviance")
                if isinstance(dev_raw, (int, float)):
                    dev = float(dev_raw)
            except Exception as e:
                self._append_log(f"Could not read cv_metrics.json: {e}")
        vid = plan["variant_id"]
        self._stage3_results.setdefault(vid, []).append({
            "fold":  plan["fold"],
            "group": plan.get("group_index"),
            "dev":   dev,
        })
        self._run_next_stage3()

    def _finish_stage3(self):
        self._run_stage1_btn.setEnabled(True)
        self._run_stage3_btn.setEnabled(bool(self._pareto_union))
        self._stop_btn.setEnabled(False)
        if self._user_stopped:
            self._status_label.setText("stopped")
            return
        self._status_label.setText("Stage 3 CV done")
        self._populate_cv_table()

    def _populate_cv_table(self):
        # Aggregate per variant: total deviance per fold (Arch A sums across groups)
        # then mean ± SE across folds.
        rows = []
        for vid, fold_entries in self._stage3_results.items():
            # Group entries by fold and sum the holdout deviance across phase groups
            per_fold = {}
            for e in fold_entries:
                per_fold.setdefault(e["fold"], 0.0)
                if isinstance(e["dev"], float) and not math.isnan(e["dev"]):
                    per_fold[e["fold"]] += e["dev"]
                else:
                    per_fold[e["fold"]] = float("nan")
            dev_vals = [v for v in per_fold.values()
                        if isinstance(v, float) and not math.isnan(v)]
            if dev_vals:
                mean_dev = sum(dev_vals) / len(dev_vals)
                if len(dev_vals) > 1:
                    var = sum((v - mean_dev) ** 2 for v in dev_vals) / (len(dev_vals) - 1)
                    se   = math.sqrt(var / len(dev_vals))
                else:
                    se = float("nan")
            else:
                mean_dev = float("nan"); se = float("nan")
            arch = vid.split("_")[0]
            rows.append({
                "variant_id":     vid,
                "arch":           arch,
                "mean_cv_dev":    mean_dev,
                "se_cv_dev":      se,
                "n_folds":        len(dev_vals),
            })
        # Sort by mean_cv_dev ascending (lower = better)
        rows.sort(key=lambda r: (
            float("inf") if isinstance(r["mean_cv_dev"], float) and math.isnan(r["mean_cv_dev"]) else r["mean_cv_dev"]
        ))

        cols = ["variant_id", "arch", "mean_cv_dev", "se_cv_dev", "n_folds"]
        self._cv_table.setSortingEnabled(False)
        self._cv_table.setColumnCount(len(cols))
        self._cv_table.setHorizontalHeaderLabels(cols)
        self._cv_table.setRowCount(len(rows))
        for r, row in enumerate(rows):
            for c, col in enumerate(cols):
                val = row.get(col)
                if isinstance(val, (int, float)):
                    item = QTableWidgetItem(_fmt(val))
                else:
                    item = QTableWidgetItem(str(val))
                if r == 0 and col == "variant_id":
                    item.setForeground(QColor("#5fb35f"))
                    f = item.font(); f.setBold(True); item.setFont(f)
                self._cv_table.setItem(r, c, item)
        self._cv_table.setSortingEnabled(True)
        self._cv_table.resizeColumnsToContents()

        if rows and not math.isnan(rows[0].get("mean_cv_dev", float("nan"))):
            self._append_log(
                f"\nStage 3 winner by CV deviance: {rows[0]['variant_id']} "
                f"(mean dev {rows[0]['mean_cv_dev']:.4g} ± {rows[0]['se_cv_dev']:.3g})"
            )

    # ── Subprocess plumbing ──────────────────────────────────────────────

    def _launch(self, cmd, env, on_done):
        kwargs = {}
        if sys.platform == "win32":
            kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW
        else:
            kwargs["start_new_session"] = True
        try:
            self._process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
                env=env,
                **kwargs,
            )
        except Exception as e:
            self._append_log(f"Failed to launch R: {e}")
            on_done.emit(self._current_run_id, False)
            return

        rid = self._current_run_id
        t = threading.Thread(
            target=self._stream_output,
            args=(self._process, rid, on_done),
            daemon=True,
        )
        t.start()

    def _stream_output(self, proc, run_id, on_done):
        for line in iter(proc.stdout.readline, ""):
            self._log_line_ready.emit(line.rstrip())
        proc.wait()
        on_done.emit(run_id, proc.returncode == 0)

    def _on_stop(self):
        self._user_stopped = True
        self._append_log("--- Stop requested ---")
        if not self.request_termination():
            # No process running — finish current stage immediately
            if self._current_stage == "stage1":
                self._finish_stage1()
            elif self._current_stage == "stage3":
                self._finish_stage3()
            else:
                self._stop_btn.setEnabled(False)

    def request_termination(self) -> bool:
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
                try: proc.terminate()
                except Exception: pass
        else:
            try:
                pgid = os.getpgid(proc.pid)
                os.killpg(pgid, signal.SIGTERM)
                try: proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    os.killpg(pgid, signal.SIGKILL)
            except Exception:
                try: proc.terminate()
                except Exception: pass
        return True

    # ── Path helpers ─────────────────────────────────────────────────────

    def _ablation_root(self) -> str:
        return os.path.join(self._bam_widget._output_dir, "ablation")

    def _stage1_subdir(self, plan) -> str:
        run_id = self._plan_run_id(plan)
        return os.path.join(self._ablation_root(), "stage1", run_id)

    def _stage3_subdir(self, plan) -> str:
        run_id = self._stage3_run_id(plan)
        return os.path.join(self._ablation_root(), "stage3_cv", run_id)

    @staticmethod
    def _stage3_run_id(plan) -> str:
        parts = [plan["arch"], plan["re"], plan["by"], plan["ar"],
                 f"fold{plan['fold']}"]
        if plan.get("group_index") is not None:
            parts.append(f"g{plan['group_index']}")
        return "_".join(parts)

    def _infer_plan_from_run_id(self, run_id):
        # Reverse the encoding used in _plan_run_id
        parts = run_id.split("_")
        if not parts:
            return None
        arch = parts[0]
        if arch == "A":
            # A_<re>_g<n>
            re = parts[1]
            gpart = parts[-1]
            if not gpart.startswith("g"):
                return None
            return {
                "arch": "A", "re": re, "by": "cond_combo",
                "ar":   "animal",
                "group_index": int(gpart[1:]),
            }
        # <arch>_<re>_<by>_<ar>
        ar = parts[-1]
        if ar not in _AR_START_MODES:
            return None
        for cand_by in _BY_VARIABLES:
            suffix = f"_{cand_by}_{ar}"
            if run_id.endswith(suffix):
                re = run_id[len(arch) + 1 : -len(suffix)]
                if re not in _RE_BASES:
                    return None
                return {
                    "arch": arch, "re": re, "by": cand_by, "ar": ar,
                    "group_index": None,
                }
        return None

    def _infer_stage3_plan_from_run_id(self, run_id):
        # <stage1_run_id>_foldN(_gG for Arch A)
        # Find "_fold<N>"
        import re as _re
        m = _re.search(r"_fold(\d+)(?:_g(\d+))?$", run_id)
        if not m:
            return None
        fold = int(m.group(1))
        g    = int(m.group(2)) if m.group(2) is not None else None
        prefix_end = m.start()
        s1_run_id = run_id[:prefix_end]
        plan = self._infer_plan_from_run_id(s1_run_id + (f"_g{g}" if g is not None else ""))
        if plan is None:
            return None
        variant_id = self._plan_variant_id(plan)
        plan["variant_id"] = variant_id
        plan["fold"] = fold
        return plan

    # ── CSV + sidecars (duplicated from bam_widget for self-containedness) ──

    def _export_csv_to(self, output_dir: str):
        bam = self._bam_widget
        df = bam._df.copy()
        df["pixel_diff"] = df["pxl_diff"] if "pxl_diff" in df.columns else 0
        df["animal_id"] = df["plate"].astype(str) + "_" + df["location"].astype(str)
        export_cols = [
            "time_sec", "location", "loc_coord", "pixel_diff",
            "Condition", "Phase", "Group", "animal_id", "plate",
        ]
        for var in bam._variable_names:
            if var in df.columns and var not in export_cols:
                export_cols.append(var)
        missing = [c for c in export_cols if c not in df.columns]
        if missing:
            QMessageBox.critical(
                self, "Missing Columns",
                f"Cannot export — missing columns:\n{missing}"
            )
            return None
        path = os.path.join(output_dir, "finomena_pre-processed_data.csv")
        try:
            df[export_cols].to_csv(path, index=False)
        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))
            return None
        return path

    def _write_sidecars(self, output_dir: str, skip_posterior: bool):
        if self._contrast_widget is not None:
            try:
                pe_settings = self._contrast_widget.get_posterior_settings() or {}
                if skip_posterior:
                    # Force-disable posterior equivalence regardless of UI state
                    pe_settings = dict(pe_settings)
                    pe_settings["enabled"] = False
                sidecar = {
                    "kept_pairs":           self._contrast_widget.get_kept_pairs(),
                    "posterior_equivalence": pe_settings,
                }
                with open(os.path.join(output_dir, "contrast_selection.json"),
                          "w", encoding="utf-8") as f:
                    json.dump(sidecar, f, indent=2)
            except Exception as e:
                self._append_log(f"Could not write contrast_selection.json: {e}")
        if self._correction_widget is not None:
            try:
                self._correction_widget.set_output_dir(output_dir)
                self._correction_widget.write_sidecar()
            except Exception as e:
                self._append_log(f"Could not write correction.json: {e}")

    def _build_cmd(self, rscript: str, csv_path: str, output_dir: str) -> list:
        bam = self._bam_widget
        var_names_str  = ",".join(bam._variable_names)
        ref_values_str = ",".join(
            bam._variable_refs.get(v, "") for v in bam._variable_names
        )
        roles_str = ",".join(
            f"{cond}={_ROLE_ABBREV.get(role, role)}"
            for cond, role in bam._roles.items()
            if role
        )
        return [
            rscript, _DEFAULT_R_SCRIPT,
            csv_path, output_dir,
            var_names_str, ref_values_str, bam._ref_condition,
            "BH", "none",
            roles_str,
            bam._family_name, str(bam._shift_val),
        ]

    # ── Log ──────────────────────────────────────────────────────────────

    def _append_log(self, line: str):
        self._log.append(line)
        self._log.verticalScrollBar().setValue(
            self._log.verticalScrollBar().maximum()
        )

"""
Catch22 Clustering Widget
=========================
Ports Cells 5 + 6 of catch22 clustering.ipynb into a PySide6 GUI.

Receives the preprocessed full_df from the Python pipeline,
extracts CATCH24 time-series features per animal per phase,
generates phase-by-phase Chebyshev-distance clustermaps,
and reports the top biological drivers between each role and the reference.
"""

import io
import os
import threading
from typing import Optional

import numpy as np
import pandas as pd

from PySide6.QtCore import Qt, Signal
from PySide6.QtCore import QByteArray
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QGroupBox, QHBoxLayout, QLabel, QMessageBox,
    QProgressBar, QPushButton, QSizePolicy, QSplitter, QTextEdit,
    QVBoxLayout, QWidget
)

from figure_viewer import FigureViewerWidget

# ── CATCH24 feature descriptions (verbatim from Cell 6) ───────────────────────
CATCH24_DESCRIPTIONS = {
    "DN_HistogramMode_5":                           "Distribution mode (5 bins)",
    "DN_HistogramMode_10":                          "Distribution mode (10 bins)",
    "CO_f1ecac":                                    "First 1/e autocorrelation crossing",
    "CO_FirstMin_ac":                               "First autocorrelation minimum",
    "CO_HistogramAMI_even_2_5":                     "Auto-mutual information (lag 2, 5 bins)",
    "CO_trev_1_num":                                "Time-reversibility statistic (lag 1)",
    "CO_Embed2_Dist_tau_d_expfit_meandiff":         "Embedding distance exponential fit",
    "IN_AutoMutualInfoStats_40_gaussian_fmmi":      "First min. of auto-mutual information",
    "MD_hrv_classic_pnn40":                         "HRV: consecutive differences > 40ms",
    "SB_BinaryStats_mean_longstretch1":             "Mean longest above-mean stretch",
    "SB_BinaryStats_diff_longstretch0":             "Longest flat (unchanged) stretch",
    "SB_MotifThree_quantile_hh":                    "High-high motif entropy (3-quantile)",
    "SB_TransitionMatrix_3ac_sumdiagcov":           "Transition matrix diagonal covariance",
    "SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1":       "DFA scaling exponent",
    "SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1":  "Long-range dependence (DFA range)",
    "SP_Summaries_welch_rect_area_5_1":             "Power in lowest 20% of frequencies",
    "SP_Summaries_welch_rect_centroid":             "Welch power spectrum centroid",
    "FC_LocalSimple_mean1_tauresrat":               "Local mean forecast error (lag 1)",
    "FC_LocalSimple_mean3_stderr":                  "Local mean forecast std error (lag 3)",
    "DN_OutlierInclude_p_001_mdrmd":                "Positive outlier median (0.1%)",
    "DN_OutlierInclude_n_001_mdrmd":                "Negative outlier median (0.1%)",
    "PD_PeriodicityWang_th0_01":                    "Periodicity (Wang method)",
    "mean":                                         "Mean (catch24)",
    "variance":                                     "Variance (catch24)",
}


class Catch22Widget(QWidget):
    """
    Tab widget for CATCH24 feature extraction and phase-by-phase clustermap
    visualization.  Compares all non-reference conditions against the reference
    control to identify top biological drivers per role.
    """

    # Public signal
    analysis_complete = Signal()

    # Private thread-safe signals
    _progress_signal = Signal(int, str)           # (percent, message)
    _analysis_done   = Signal(bool, str, list)    # (success, features_text, figures_list)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._df            = None
        self._scaled_df     = None  # kept for re-running top-features
        self._pending_figs  = []
        self._ref_condition = ""
        self._roles         = {}    # {condition_name: role_str}

        self._build_ui()
        self._progress_signal.connect(self._on_progress)
        self._analysis_done.connect(self._on_analysis_done)

    # ── Public API ────────────────────────────────────────────────────────────

    def set_references(self, variable_refs: dict, ref_condition: str):
        """Sets the reference (baseline) condition."""
        self._ref_condition = ref_condition
        self._update_ref_label()

    def set_roles(self, roles: dict):
        """Stores the {condition_name: role_str} mapping from MetadataAssignment."""
        self._roles = dict(roles)
        self._update_ref_label()

    def _update_ref_label(self):
        """Updates the reference condition display label."""
        if self._ref_condition:
            self._ref_label.setText(
                f"Reference control: <b>{self._ref_condition}</b> — "
                f"all other conditions will be compared against it."
            )
        else:
            self._ref_label.setText(
                "No reference condition set. "
                "Assign roles in the Metadata Assignment tab."
            )

    def load_data(self, df: pd.DataFrame):
        """Receives the phase-assigned full_df from the pipeline."""
        self._df = df
        if df is not None and not df.empty:
            conditions = sorted(df['Condition'].dropna().unique()) if 'Condition' in df.columns else []
            n_rows  = len(df)
            n_conds = len(conditions)
            self._data_status_label.setText(
                f"Data loaded: {n_rows:,} rows, {n_conds} condition(s)"
            )
            self._run_button.setEnabled(True)
        else:
            self._data_status_label.setText("No data loaded.")
            self._run_button.setEnabled(False)

    # ── Build UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # ── 1. Data Status ────────────────────────────────────────────────────
        status_group = QGroupBox("1. Data Status")
        sg_layout = QVBoxLayout(status_group)
        self._data_status_label = QLabel(
            "No data loaded. Process data in the Data Loading tab first."
        )
        self._data_status_label.setWordWrap(True)
        sg_layout.addWidget(self._data_status_label)
        layout.addWidget(status_group)

        # ── 2. Reference Info ─────────────────────────────────────────────────
        ref_group = QGroupBox("2. Reference Control")
        rg_layout = QVBoxLayout(ref_group)
        self._ref_label = QLabel(
            "No reference condition set. "
            "Assign roles in the Metadata Assignment tab."
        )
        self._ref_label.setWordWrap(True)
        rg_layout.addWidget(self._ref_label)
        layout.addWidget(ref_group)

        # ── 3. Run Analysis ───────────────────────────────────────────────────
        run_group = QGroupBox("3. Run Analysis")
        rg_outer = QVBoxLayout(run_group)

        self._catch24_checkbox = QCheckBox("Include mean && variance features (CATCH24 — 24 features total)")
        self._catch24_checkbox.setChecked(True)
        self._catch24_checkbox.setToolTip(
            "CATCH22 extracts 22 time-series features.\n"
            "CATCH24 adds the signal mean and variance as features 23 and 24.\n"
            "Uncheck to use only the 22 interpretable dynamics features,\n"
            "which can be useful when mean activity level differences would\n"
            "dominate the clustering and obscure subtler temporal patterns."
        )
        rg_outer.addWidget(self._catch24_checkbox)

        rg_row = QHBoxLayout()
        self._run_button = QPushButton("Run catch22 Analysis")
        self._run_button.setEnabled(False)
        self._run_button.setMinimumHeight(36)
        self._run_button.setStyleSheet("QPushButton { font-weight: bold; }")
        self._run_button.clicked.connect(self._on_run)
        self._progress_bar = QProgressBar()
        self._progress_bar.setRange(0, 100)
        self._progress_bar.setValue(0)
        self._progress_bar.setTextVisible(True)
        self._status_label = QLabel("")
        rg_row.addWidget(self._run_button)
        rg_row.addWidget(self._progress_bar, stretch=1)
        rg_row.addWidget(self._status_label)
        rg_outer.addLayout(rg_row)
        layout.addWidget(run_group)

        # Feature reference note
        from paths import resource_path
        readme_path = resource_path("README.md")
        feature_note = QLabel(
            f'For a description of each CATCH22/24 feature, see '
            f'<a href="file:///{readme_path.replace(os.sep, "/")}">README.md — '
            f'CATCH22/24 Feature Reference</a>.'
        )
        feature_note.setOpenExternalLinks(True)
        feature_note.setStyleSheet("font-style: italic; padding: 2px 4px;")
        layout.addWidget(feature_note)

        # ── 4 + 5: Resizable splitter between clustermaps and top features ────
        splitter = QSplitter(Qt.Vertical)

        clust_group = QGroupBox("4. Phase Clustermaps")
        cl_layout = QVBoxLayout(clust_group)
        self._figure_viewer = FigureViewerWidget(title="Condition Clustermaps by Phase")
        cl_layout.addWidget(self._figure_viewer)
        splitter.addWidget(clust_group)

        feat_group = QGroupBox("5. Top Driving Features (vs Reference)")
        fg_layout = QVBoxLayout(feat_group)
        self._features_text = QTextEdit()
        self._features_text.setReadOnly(True)
        mono_font = QFont("Courier New", 9)
        self._features_text.setFont(mono_font)
        fg_layout.addWidget(self._features_text)
        splitter.addWidget(feat_group)

        # Start with clustermaps ~75% and features ~25%
        splitter.setSizes([450, 150])
        layout.addWidget(splitter, stretch=1)

    # ── Run ───────────────────────────────────────────────────────────────────

    def _on_run(self):
        if self._df is None or self._df.empty:
            QMessageBox.warning(self, "No Data", "No data loaded yet.")
            return

        if not self._ref_condition:
            QMessageBox.warning(
                self, "No Reference",
                "No reference condition is set.\n\n"
                "Go to Experimental Design → Metadata Assignment and assign "
                "a 'Reference Control' role first."
            )
            return

        try:
            import pycatch22  # noqa: F401
        except ImportError:
            QMessageBox.critical(
                self, "Missing Dependency",
                "pycatch22 is not installed.\n\nInstall it with:\n    pip install pycatch22"
            )
            return

        try:
            from sklearn.preprocessing import RobustScaler  # noqa: F401
        except Exception as exc:
            import traceback
            QMessageBox.critical(
                self, "Missing Dependency",
                f"scikit-learn import failed:\n\n{traceback.format_exc()}"
            )
            return

        catch24  = self._catch24_checkbox.isChecked()

        self._run_button.setEnabled(False)
        self._progress_bar.setValue(0)
        self._status_label.setText("Starting…")
        self._features_text.clear()
        self._figure_viewer.clear()

        work_df = self._df.copy()
        # Compute pixel_diff if not already present
        if 'pxl_diff' not in work_df.columns:
            data_cols = [c for c in work_df.columns if c.startswith('data1_')]
            if data_cols:
                work_df['pxl_diff'] = work_df[data_cols].sum(axis=1)
            else:
                work_df['pxl_diff'] = 0

        t = threading.Thread(
            target=self._analysis_worker,
            args=(work_df, self._ref_condition, catch24),
            daemon=True
        )
        t.start()

    # ── Background worker ─────────────────────────────────────────────────────

    def _analysis_worker(self, work_df: pd.DataFrame, ref_condition: str, catch24: bool = True):
        """Runs Cell 6 logic entirely in a background thread."""
        try:
            import pycatch22
            import seaborn as sns
            import matplotlib.pyplot as plt
            from sklearn.preprocessing import RobustScaler

            # ── STEP 1: Feature Extraction ────────────────────────────────────
            self._progress_signal.emit(5, "Extracting CATCH24 features…")

            sorted_df = work_df.sort_values(['loc_coord', 'Group', 'time_sec'])

            # Determine feature names dynamically from the chosen feature set
            dummy = pycatch22.catch22_all([0.0, 0.0, 0.0], catch24=catch24)
            base_feature_names = dummy['names']

            # Determine unique groups (phases) — don't assume exactly 8
            groups = sorted(sorted_df['Group'].dropna().unique())
            n_groups = len(groups)

            col_names = [
                f"{fname}_Phase{int(g)}"
                for g in groups
                for fname in base_feature_names
            ]

            all_features, subject_ids, true_labels = [], [], []

            subjects = sorted_df['loc_coord'].unique()
            n_subjects = len(subjects)

            for s_idx, subject in enumerate(subjects):
                self._progress_signal.emit(
                    5 + int(45 * s_idx / max(n_subjects, 1)),
                    f"Extracting features: animal {s_idx + 1}/{n_subjects}"
                )
                subj_df = sorted_df[sorted_df['loc_coord'] == subject]
                subject_ids.append(subject)
                true_labels.append(subj_df['Condition'].iloc[0])

                subject_features = []
                for g in groups:
                    phase_df = subj_df[subj_df['Group'] == g]
                    if phase_df.empty:
                        # Fill with zeros if this animal has no data for this phase
                        subject_features.extend([0.0] * len(base_feature_names))
                    else:
                        ts = np.log1p(
                            np.clip(phase_df['pxl_diff'].values.astype(float), 0, None)
                        )
                        result = pycatch22.catch22_all(list(ts), catch24=catch24)
                        subject_features.extend(result['values'])

                all_features.append(subject_features)

            feature_matrix = pd.DataFrame(
                all_features, index=subject_ids, columns=col_names
            )

            # Clean matrix
            feature_matrix.replace([np.inf, -np.inf], np.nan, inplace=True)
            feature_matrix.dropna(axis=1, how='all', inplace=True)
            feature_matrix.fillna(0, inplace=True)
            feature_matrix = feature_matrix.loc[
                :, feature_matrix.var(axis=0) >= 1e-6
            ]

            # ── STEP 2: Scaling ───────────────────────────────────────────────
            self._progress_signal.emit(52, "Scaling features (RobustScaler)…")
            scaled_features = RobustScaler().fit_transform(feature_matrix)
            scaled_df = pd.DataFrame(
                scaled_features,
                index=feature_matrix.index,
                columns=feature_matrix.columns,
            )
            scaled_df['Condition'] = true_labels
            self._scaled_df = scaled_df  # store for re-use

            # ── STEP 3: Phase-by-phase clustermaps ────────────────────────────
            figures = []
            for phase_idx, g in enumerate(groups):
                phase_num = int(g)
                self._progress_signal.emit(
                    52 + int(40 * (phase_idx + 1) / max(n_groups, 1)),
                    f"Generating clustermap for Phase {phase_num} ({phase_idx + 1}/{n_groups})…"
                )

                phase_cols = [c for c in feature_matrix.columns
                              if c.endswith(f'_Phase{phase_num}')]
                if not phase_cols:
                    continue

                phase_centroids = (
                    scaled_df[phase_cols + ['Condition']]
                    .groupby('Condition')
                    .median()
                )

                # Translate column names to English descriptions
                phase_centroids.columns = [
                    CATCH24_DESCRIPTIONS.get(
                        c.rsplit('_Phase', 1)[0],
                        c.rsplit('_Phase', 1)[0]
                    )
                    for c in phase_centroids.columns
                ]

                clustermap_kws = dict(
                    method='average',
                    metric='chebyshev',
                    cmap='vlag',
                    center=0,
                    figsize=(9, 13),
                    col_cluster=True,
                    row_cluster=True,
                    cbar_kws={'label': 'Robust Z-Score (Median)'},
                    dendrogram_ratio=(0.15, 0.2),
                    linewidths=0.3,
                    linecolor='white',
                )

                # Handle case with fewer than 2 samples (clustermap requires >= 2)
                if phase_centroids.shape[0] < 2:
                    clustermap_kws['col_cluster'] = False
                    clustermap_kws['row_cluster'] = False

                try:
                    g_plot = sns.clustermap(phase_centroids.T, **clustermap_kws)
                    plt.setp(
                        g_plot.ax_heatmap.get_xticklabels(),
                        rotation=45, ha='right', fontsize=11, fontweight='bold'
                    )
                    plt.setp(
                        g_plot.ax_heatmap.get_yticklabels(),
                        rotation=0, fontsize=9
                    )
                    g_plot.fig.suptitle(
                        f"Phase Group {phase_num} — Condition Clustering",
                        fontsize=14, fontweight='bold', y=1.02
                    )

                    buf = io.BytesIO()
                    g_plot.savefig(buf, format='png', bbox_inches='tight', dpi=150)
                    plt.close(g_plot.fig)
                    buf.seek(0)

                    ba = QByteArray(buf.read())
                    pixmap = QPixmap()
                    pixmap.loadFromData(ba, "PNG")

                    figures.append({
                        'pixmap':   pixmap,
                        'title':    f"Phase {phase_num}",
                        'filepath': None,
                    })
                except Exception as e:
                    self._progress_signal.emit(
                        52 + int(40 * (phase_idx + 1) / max(n_groups, 1)),
                        f"Warning: Phase {phase_num} clustermap failed: {e}"
                    )

            # ── STEP 4: Top driving features (all conditions vs reference) ────
            self._progress_signal.emit(95, "Extracting top driving features…")
            features_report = self._extract_all_vs_reference(
                scaled_df, ref_condition, top_k=5
            )

            self._progress_signal.emit(100, f"Done. {len(figures)} clustermaps generated.")
            self._analysis_done.emit(True, features_report, figures)

        except Exception as e:
            import traceback
            self._analysis_done.emit(False, f"{e}\n\n{traceback.format_exc()}", [])

    # ── Thread-safe slots ─────────────────────────────────────────────────────

    def _on_progress(self, percent: int, message: str):
        self._progress_bar.setValue(percent)
        self._status_label.setText(message)

    def _on_analysis_done(self, success: bool, message: str, figures: list):
        self._run_button.setEnabled(True)
        if success:
            self._figure_viewer.load_figures(figures)
            self._features_text.setPlainText(message)
            self.analysis_complete.emit()
        else:
            self._progress_bar.setValue(0)
            self._status_label.setText("Error — see details.")
            QMessageBox.critical(self, "catch22 Analysis Error", message)

    # ── Top features helper ───────────────────────────────────────────────────

    @staticmethod
    def _extract_all_vs_reference(
        scaled_df: pd.DataFrame,
        ref_condition: str,
        top_k: int = 5,
    ) -> str:
        """
        Returns a formatted string of the top driving features between
        each non-reference condition and the reference.
        """
        lines = []
        conditions = sorted(scaled_df['Condition'].unique())

        if ref_condition not in conditions:
            return f"Reference condition '{ref_condition}' not found in data."

        ref_centroid = (
            scaled_df[scaled_df['Condition'] == ref_condition]
            .drop(columns='Condition')
            .median(axis=0)
        )

        other_conditions = [c for c in conditions if c != ref_condition]

        for cond in other_conditions:
            cond_centroid = (
                scaled_df[scaled_df['Condition'] == cond]
                .drop(columns='Condition')
                .median(axis=0)
            )

            abs_diff = np.abs(cond_centroid - ref_centroid)
            top_features = abs_diff.sort_values(ascending=False).head(top_k)

            lines.append("\n" + "=" * 70)
            lines.append(f"TOP BIOLOGICAL DRIVERS: {cond} vs {ref_condition} (Reference)")
            lines.append("=" * 70)
            for feat, diff in top_features.items():
                base_feat = feat.rsplit('_Phase', 1)[0]
                desc = CATCH24_DESCRIPTIONS.get(base_feat, base_feat)
                phase_part = feat.rsplit('_Phase', 1)
                phase_str  = f"Phase {phase_part[1]}" if len(phase_part) == 2 else ""
                lines.append(f"Difference: {diff:.3f} robust units")
                lines.append(f"Feature:    {feat}  ({phase_str})")
                lines.append(f"Meaning:    {desc}\n")

        return "\n".join(lines)

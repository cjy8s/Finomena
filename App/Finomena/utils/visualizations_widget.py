"""
Visualizations Widget
=====================
Interactive plots over the BAM master-results CSV. Replaces R's static PNG
output for the main contrast results. Renders entirely in Python so users can
change plot type, filters, and layout without re-running R.

Plot types
----------
- Volcano:   log_2 FC (x) vs -log_10 FDR p-value (y). One dot per contrast;
             scales gracefully to hundreds of conditions.
- Top-N forest: classic forest plot truncated to the N most significant (or
             largest-|effect|) contrasts. Replaces the unscalable full forest.
- Beeswarm:  distribution of effect sizes within each Test_Family / Group.
             Shows the population of effects, not just individuals.
- PE density: posterior densities of M (max absolute trajectory difference)
             with the equivalence margin delta marked. Optional — only
             available when posterior equivalence was enabled.
- PE Pr(M<delta): heatmap of the posterior equivalence probability per pair
             x group. Same condition: only when PE was enabled.

Data contract
-------------
Reads ``master_results.csv`` from the chosen output directory. Required cols:
``Test_Family, Group, Split_By, Tested_Level, estimate, SE, raw_pvalue,
adjusted_pvalue``. Optional cols: ``correction_method``, ``tree_level``,
``tree_status`` (NA-filled for flat correction).

Optionally reads ``posterior_equivalence_draws.csv`` (long: ``Group``, ``pair``,
``M_log2``) and ``posterior_equivalence_summary.csv`` (one row per pair x group
with ``Pr_equiv`` and ``delta_log2``) when present. Both are written by
[TweedieAR1 BAM.R] only if posterior equivalence is enabled in the contrast
selection.
"""

import os
import re

import numpy as np
import pandas as pd

from PySide6.QtCore import Qt, QEvent, QTimer
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox, QColorDialog, QComboBox, QDoubleSpinBox, QFileDialog,
    QGroupBox, QHBoxLayout, QLabel, QMessageBox, QPushButton,
    QScrollArea, QSizePolicy, QSpinBox, QSplitter, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as _FigureCanvasBase
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavToolbar
from matplotlib.figure import Figure


class FigureCanvas(_FigureCanvasBase):
    """FigureCanvas that lets wheel events propagate to the parent QScrollArea.

    The default Qt-backend canvas accepts wheel events so it can emit
    matplotlib's ``scroll_event`` for plugin code that wants to zoom/pan. We
    don't subscribe to that, and accepting blocks the QScrollArea from ever
    seeing the wheel event — so the user can't mouse-wheel-scroll tall plots
    like the top-N forest or the PE heatmap.
    """

    def wheelEvent(self, event):
        event.ignore()  # bubble to the QScrollArea wrapping the canvas


# ── Figure styles ─────────────────────────────────────────────────────────────
# Two themes for the matplotlib canvas. The choice is independent of the Qt
# (qdarktheme) chrome — users may want light-mode plots even with a dark app
# UI for publication or printing.
_DARK_STYLE = {
    "figure.facecolor":  "#202124",
    "axes.facecolor":    "#202124",
    "axes.edgecolor":    "#cccccc",
    "axes.labelcolor":   "#e8e8e8",
    "axes.titlecolor":   "#e8e8e8",
    "xtick.color":       "#cccccc",
    "ytick.color":       "#cccccc",
    "text.color":        "#e8e8e8",
    "grid.color":        "#444444",
    "grid.linestyle":    ":",
    "axes.grid":         True,
    "axes.grid.which":   "major",
    "savefig.facecolor": "#202124",
}

_LIGHT_STYLE = {
    "figure.facecolor":  "#ffffff",
    "axes.facecolor":    "#ffffff",
    "axes.edgecolor":    "#000000",
    "axes.labelcolor":   "#000000",
    "axes.titlecolor":   "#000000",
    "xtick.color":       "#000000",
    "ytick.color":       "#000000",
    "text.color":        "#000000",
    "grid.color":        "#dddddd",
    "grid.linestyle":    ":",
    "axes.grid":         True,
    "axes.grid.which":   "major",
    "savefig.facecolor": "#ffffff",
}

# Default colors. Significance color is user-configurable per session via the
# Theme controls; the others follow the theme.
_DEFAULT_SIG_COLOR    = "#ff7043"   # orange — works on both backgrounds
_NOTSIG_COLOR_DARK    = "#9e9e9e"   # light grey, visible on dark
_NOTSIG_COLOR_LIGHT   = "#757575"   # darker grey, visible on white
_COLOR_REF            = "#42a5f5"   # blue — threshold lines stay blue in both themes


class VisualizationsWidget(QWidget):
    """Interactive matplotlib plots of BAM master results."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._output_dir: str | None = None
        self._df: pd.DataFrame | None = None
        self._pe_draws_df: pd.DataFrame | None = None
        self._pe_summary_df: pd.DataFrame | None = None
        self._pe_delta: float | None = None
        self._roles: dict = {}           # {condition_name: role_str}
        self._role_abbrev_map: dict = {
            "Reference Control": "RC",
            "Experimental":      "EXP",
            "Positive Control":  "PC",
            "Potential Rescue":  "RES",
        }
        # Theme state — figure background + line/text colors. Independent of
        # the Qt app theme. Default to dark mode to match the rest of the UI;
        # the user can flip to light mode (white bg, black ink) for export.
        self._theme: str = "dark"
        self._cur_style: dict = _DARK_STYLE
        self._sig_color: str = _DEFAULT_SIG_COLOR
        # Threshold reference lines (volcano ±1 lines, FDR cutoff, y=0, the
        # beeswarm median tick, etc.) share one user-pickable color.
        self._ref_color: str = _COLOR_REF
        # Confidence-interval ribbon color (cross-group lines plot, etc.).
        # Default to a neutral mid-grey readable in either theme.
        self._ci_color: str = "#9e9e9e"
        # M-density histogram fill color (PE density plot). Independent of
        # the significant-point color so users can tune the two separately.
        self._density_color: str = _DEFAULT_SIG_COLOR
        # Heatmap colormap (PE Pr(M < δ) plot). Sequential cmaps make most
        # sense for a probability in [0, 1]; diverging options are useful when
        # users want to anchor at 0.5.
        self._heatmap_cmap: str = "viridis"
        # Global multiple-testing correction method. Used purely for axis
        # labels: drives whether the corrected-p axis reads "FDR" (BH) or
        # "FWER" (Holm). Set via set_correction_method() from app.py.
        self._correction_method: str = "BH"
        # Seed QColorDialog's custom-color slots with the three defaults so
        # users can one-click revert to the originals from inside the picker.
        # Static method — setting once affects every dialog in this app session.
        for slot, hex_default in (
            (0, _DEFAULT_SIG_COLOR),   # significant-point default
            (1, _COLOR_REF),           # threshold-line default
            (2, self._ci_color),       # CI ribbon default
            (3, self._density_color),  # M-density histogram default
        ):
            try:
                QColorDialog.setCustomColor(slot, QColor(hex_default))
            except Exception:
                # On some Qt builds the static API expects QColor(int rgb);
                # fall through silently if the call signature differs.
                pass
        self._suspend_render = False     # set during bulk control updates
        self._build_ui()

    # ── Public API ────────────────────────────────────────────────────────────

    def set_output_dir(self, path: str):
        self._output_dir = path
        self._try_load()

    def set_roles(self, roles: dict):
        """Receives {condition_name: role_str} from metadata assignment."""
        self._roles = dict(roles)
        self._render()

    def set_correction_method(self, method: str):
        """Tell the widget which global correction the user chose ("BH" / "holm").

        Drives axis-label terminology only; the underlying p-values come from
        R's master_results.csv (whose column name is fixed as
        ``adjusted_pvalue`` regardless of method).
        """
        self._correction_method = method or "BH"
        self._render()

    def _correction_label(self) -> str:
        """Plain-language correction tag for axis/legend labels."""
        m = str(self._correction_method).lower()
        if m in ("holm", "holm-gatekeeping"):
            return "Holm-FWER"
        if m in ("gmcp", "graphicalmcp"):
            return "graphicalMCP-FWER"
        if m == "treebh":
            return "TreeBH-FDR"
        # Default / unknown → assume BH-FDR
        return "BH-FDR"

    @staticmethod
    def _pretty(s) -> str:
        """Underscore-to-space conversion for publication-style display.

        Test_Family values arrive as R-friendly identifiers
        ("Drug_Effect", "Full_Interaction"); on plot titles and axis labels
        we want them as natural English ("Drug Effect", "Full Interaction").
        """
        return str(s).replace("_", " ")

    def on_bam_completed(self, output_dir: str):
        """Called when a BAM run finishes — reloads the master CSV."""
        self._output_dir = output_dir
        self._try_load()

    # ── Build UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        outer = QVBoxLayout(self)

        splitter = QSplitter(Qt.Horizontal)
        outer.addWidget(splitter, stretch=1)

        # Left: controls panel
        controls_box = QWidget()
        controls_layout = QVBoxLayout(controls_box)
        controls_layout.setAlignment(Qt.AlignTop)
        controls_box.setMinimumWidth(280)
        controls_box.setMaximumWidth(360)

        # Plot type
        type_group = QGroupBox("Plot type")
        type_layout = QVBoxLayout(type_group)
        self._plot_type = QComboBox()
        # The dropdown is grouped into three conceptual buckets so users can
        # pick the plot that matches their data shape at a glance:
        #   - Single-group views are designed for one phase × many conditions
        #     (volcano / forest / beeswarm scale well across hundreds of dots).
        #   - Cross-group views shine when the same contrast is repeated across
        #     several phase groups (effect-size trajectories with CI ribbons).
        #   - Posterior equivalence plots work in either regime.
        self._add_plot_type_header("Single-group views (many conditions)")
        self._plot_type.addItem("Volcano")
        self._plot_type.addItem("Top-N forest")
        self._plot_type.addItem("Beeswarm")
        self._plot_type.insertSeparator(self._plot_type.count())
        self._add_plot_type_header("Cross-group views (many phases)")
        self._plot_type.addItem("Effect lines across groups")
        self._plot_type.insertSeparator(self._plot_type.count())
        self._add_plot_type_header("Posterior equivalence (either)")
        self._plot_type.addItem("Posterior eq.: M density")
        self._plot_type.addItem("Posterior eq.: Probability of equivalence")
        # Headers are non-selectable, but the very first item must still be
        # the default. Snap to "Volcano" (the first real plot).
        idx_volcano = self._plot_type.findText("Volcano")
        if idx_volcano >= 0:
            self._plot_type.setCurrentIndex(idx_volcano)
        self._plot_type.currentIndexChanged.connect(self._on_plot_type_changed)
        type_layout.addWidget(self._plot_type)
        # Hint shown when a PE plot type is selected but PE files are absent
        self._pe_availability_label = QLabel("")
        self._pe_availability_label.setWordWrap(True)
        self._pe_availability_label.setStyleSheet(
            "color: #ffb74d; font-style: italic;"
        )
        self._pe_availability_label.hide()
        type_layout.addWidget(self._pe_availability_label)
        controls_layout.addWidget(type_group)

        # Filters
        filt_group = QGroupBox("Filters")
        filt_layout = QVBoxLayout(filt_group)

        filt_layout.addWidget(QLabel("Test family"))
        self._family_combo = QComboBox()
        self._family_combo.currentIndexChanged.connect(self._render)
        filt_layout.addWidget(self._family_combo)

        filt_layout.addWidget(QLabel("Group (phase)"))
        self._group_combo = QComboBox()
        self._group_combo.currentIndexChanged.connect(self._render)
        filt_layout.addWidget(self._group_combo)

        filt_layout.addWidget(QLabel("Split-By"))
        self._splitby_combo = QComboBox()
        self._splitby_combo.currentIndexChanged.connect(self._render)
        filt_layout.addWidget(self._splitby_combo)

        sig_row = QHBoxLayout()
        sig_row.addWidget(QLabel("Significant p-value threshold:"))
        self._sig_thresh = QDoubleSpinBox()
        self._sig_thresh.setRange(0.0001, 0.5)
        self._sig_thresh.setSingleStep(0.01)
        self._sig_thresh.setDecimals(4)
        self._sig_thresh.setValue(0.05)
        self._sig_thresh.valueChanged.connect(self._render)
        sig_row.addWidget(self._sig_thresh)
        filt_layout.addLayout(sig_row)

        eff_row = QHBoxLayout()
        eff_row.addWidget(QLabel("Min |log₂ FC|:"))
        self._effect_min = QDoubleSpinBox()
        self._effect_min.setRange(0.0, 20.0)
        self._effect_min.setSingleStep(0.1)
        self._effect_min.setDecimals(2)
        self._effect_min.setValue(0.0)
        self._effect_min.valueChanged.connect(self._render)
        eff_row.addWidget(self._effect_min)
        filt_layout.addLayout(eff_row)

        topn_row = QHBoxLayout()
        topn_row.addWidget(QLabel("Top-N rows:"))
        self._top_n = QSpinBox()
        self._top_n.setRange(1, 500)
        self._top_n.setValue(30)
        self._top_n.valueChanged.connect(self._render)
        topn_row.addWidget(self._top_n)
        filt_layout.addLayout(topn_row)

        self._rank_by_combo = QComboBox()
        self._rank_by_combo.addItem("Rank by significance (corrected p)", "fdr")
        self._rank_by_combo.addItem("Rank by |effect size|",        "effect")
        self._rank_by_combo.currentIndexChanged.connect(self._render)
        filt_layout.addWidget(self._rank_by_combo)

        # tree_status filter — only meaningful when the master CSV came from a
        # Tree-strategy correction. Hidden otherwise via _update_tree_status_visibility.
        filt_layout.addWidget(QLabel("Tree status"))
        self._tree_status_combo = QComboBox()
        self._tree_status_combo.addItem("All",                     "all")
        self._tree_status_combo.addItem("Tested only",             "tested")
        self._tree_status_combo.addItem("Tested + gate-failed",    "tested_plus_gate_failed")
        self._tree_status_combo.addItem("Not reached (gated out)", "not_reached")
        self._tree_status_combo.addItem("Off-tree (outside spec)", "off_tree")
        self._tree_status_combo.currentIndexChanged.connect(self._render)
        filt_layout.addWidget(self._tree_status_combo)

        controls_layout.addWidget(filt_group)

        # Layout / appearance
        layout_group = QGroupBox("Layout")
        layout_layout = QVBoxLayout(layout_group)

        self._show_labels = QCheckBox("Show top contrast labels (volcano + beeswarm)")
        self._show_labels.setChecked(True)
        self._show_labels.toggled.connect(self._render)
        layout_layout.addWidget(self._show_labels)

        # Choose whether labels go only on significant points or on every
        # point. Affects volcano and beeswarm. "Significant only" matches the
        # historical behaviour; "All points" labels every dot up to a
        # defensive cap so the plot doesn't turn into a wall of text.
        vlabel_row = QHBoxLayout()
        vlabel_row.addWidget(QLabel("Data-point labels:"))
        self._label_scope = QComboBox()
        self._label_scope.addItem("Significant only", "sig")
        self._label_scope.addItem("All points",       "all")
        self._label_scope.currentIndexChanged.connect(self._render)
        vlabel_row.addWidget(self._label_scope, stretch=1)
        layout_layout.addLayout(vlabel_row)

        self._show_role_tags = QCheckBox("Show role tags ([RC], [EXP], …)")
        self._show_role_tags.setChecked(True)
        self._show_role_tags.toggled.connect(self._render)
        layout_layout.addWidget(self._show_role_tags)

        self._show_group_in_label = QCheckBox("Include phase group in label (e.g. [G1])")
        self._show_group_in_label.setChecked(False)
        self._show_group_in_label.toggled.connect(self._render)
        layout_layout.addWidget(self._show_group_in_label)

        fs_row = QHBoxLayout()
        fs_row.addWidget(QLabel("Font size:"))
        self._font_size = QSpinBox()
        self._font_size.setRange(6, 24)
        self._font_size.setValue(10)
        self._font_size.valueChanged.connect(self._render)
        fs_row.addWidget(self._font_size)
        layout_layout.addLayout(fs_row)

        controls_layout.addWidget(layout_group)

        # ── Theme ────────────────────────────────────────────────────────────
        # Figure-only theme. Independent of the Qt app theme so users can
        # publish/print light-mode figures while keeping the dark UI.
        theme_group = QGroupBox("Theme")
        theme_layout = QVBoxLayout(theme_group)

        # Use a uniform label width so the three color-picker buttons line
        # up vertically with each other (and with the Mode dropdown).
        _theme_label_w = 170

        theme_row = QHBoxLayout()
        mode_lbl = QLabel("Mode:")
        mode_lbl.setMinimumWidth(_theme_label_w)
        theme_row.addWidget(mode_lbl)
        self._theme_combo = QComboBox()
        self._theme_combo.addItem("Dark mode",  "dark")
        self._theme_combo.addItem("Light mode", "light")
        self._theme_combo.currentIndexChanged.connect(self._on_theme_changed)
        theme_row.addWidget(self._theme_combo, stretch=1)
        theme_layout.addLayout(theme_row)

        sig_color_row = QHBoxLayout()
        sig_lbl = QLabel("Significant-point color:")
        sig_lbl.setMinimumWidth(_theme_label_w)
        sig_color_row.addWidget(sig_lbl)
        self._sig_color_btn = QPushButton()
        self._sig_color_btn.setFixedWidth(80)
        self._sig_color_btn.setToolTip("Click to pick a color for FDR-significant points")
        self._sig_color_btn.clicked.connect(self._pick_sig_color)
        self._update_sig_color_btn_style()
        sig_color_row.addWidget(self._sig_color_btn)
        sig_color_row.addStretch()
        theme_layout.addLayout(sig_color_row)

        ref_color_row = QHBoxLayout()
        ref_lbl = QLabel("Threshold-line color:")
        ref_lbl.setMinimumWidth(_theme_label_w)
        ref_color_row.addWidget(ref_lbl)
        self._ref_color_btn = QPushButton()
        self._ref_color_btn.setFixedWidth(80)
        self._ref_color_btn.setToolTip(
            "Click to pick a color for reference / threshold lines "
            "(volcano ±1, FDR cutoff, y=0, beeswarm median tick)"
        )
        self._ref_color_btn.clicked.connect(self._pick_ref_color)
        self._update_ref_color_btn_style()
        ref_color_row.addWidget(self._ref_color_btn)
        ref_color_row.addStretch()
        theme_layout.addLayout(ref_color_row)

        ci_color_row = QHBoxLayout()
        ci_lbl = QLabel("CI ribbon color:")
        ci_lbl.setMinimumWidth(_theme_label_w)
        ci_color_row.addWidget(ci_lbl)
        self._ci_color_btn = QPushButton()
        self._ci_color_btn.setFixedWidth(80)
        self._ci_color_btn.setToolTip(
            "Click to pick a color for confidence-interval ribbons "
            "(cross-group lines plot)"
        )
        self._ci_color_btn.clicked.connect(self._pick_ci_color)
        self._update_ci_color_btn_style()
        ci_color_row.addWidget(self._ci_color_btn)
        ci_color_row.addStretch()
        theme_layout.addLayout(ci_color_row)

        density_color_row = QHBoxLayout()
        density_lbl = QLabel("M-density color:")
        density_lbl.setMinimumWidth(_theme_label_w)
        density_color_row.addWidget(density_lbl)
        self._density_color_btn = QPushButton()
        self._density_color_btn.setFixedWidth(80)
        self._density_color_btn.setToolTip(
            "Click to pick a fill color for the M-density histograms "
            "(Posterior eq.: M density plot)"
        )
        self._density_color_btn.clicked.connect(self._pick_density_color)
        self._update_density_color_btn_style()
        density_color_row.addWidget(self._density_color_btn)
        density_color_row.addStretch()
        theme_layout.addLayout(density_color_row)

        cmap_row = QHBoxLayout()
        cmap_lbl = QLabel("Heatmap palette:")
        cmap_lbl.setMinimumWidth(_theme_label_w)
        cmap_row.addWidget(cmap_lbl)
        self._heatmap_cmap_combo = QComboBox()
        self._invert_cmap = QCheckBox("invert")
        self._invert_cmap.setToolTip(
            "Reverse the colormap direction — e.g. viridis flips so low values"
            " become bright and high values become dark."
        )
        self._invert_cmap.toggled.connect(self._render)
        for label, key in (
            ("Viridis (default)",   "viridis"),
            ("Plasma",              "plasma"),
            ("Magma",               "magma"),
            ("Inferno",             "inferno"),
            ("Cividis",             "cividis"),
            ("Turbo",               "turbo"),
            ("Blues",               "Blues"),
            ("Greens",              "Greens"),
            ("Reds",                "Reds"),
            ("Purples",             "Purples"),
            ("Oranges",             "Oranges"),
            ("RdBu (diverging)",    "RdBu"),
            ("Coolwarm (diverging)", "coolwarm"),
        ):
            self._heatmap_cmap_combo.addItem(label, key)
        self._heatmap_cmap_combo.setToolTip(
            "Colormap used for the Posterior eq.: Probability of equivalence heatmap"
        )
        self._heatmap_cmap_combo.currentIndexChanged.connect(self._on_cmap_changed)
        cmap_row.addWidget(self._heatmap_cmap_combo, stretch=1)
        cmap_row.addWidget(self._invert_cmap)
        theme_layout.addLayout(cmap_row)

        controls_layout.addWidget(theme_group)

        # Export
        export_row = QHBoxLayout()
        self._export_btn = QPushButton("Export figure…")
        self._export_btn.clicked.connect(self._on_export)
        export_row.addWidget(self._export_btn)
        controls_layout.addLayout(export_row)

        controls_layout.addStretch()

        # Status label at bottom of controls column
        self._status_label = QLabel("No results loaded.")
        self._status_label.setWordWrap(True)
        self._status_label.setStyleSheet("color: #aaaaaa; font-style: italic;")
        controls_layout.addWidget(self._status_label)

        splitter.addWidget(controls_box)

        # Right: figure canvas
        # The canvas lives inside a QScrollArea with widgetResizable=False so
        # tall figures (many-panel forest / lines / heatmap) scroll vertically
        # while always matching the viewport width. The previous design called
        # figure.set_size_inches(..., forward=True) in each render which both
        # fought the Qt layout and persisted oversize state across renders.
        fig_box = QWidget()
        fig_layout = QVBoxLayout(fig_box)
        fig_layout.setContentsMargins(0, 0, 0, 0)
        self._figure = Figure(figsize=(8, 6))
        self._figure.patch.set_facecolor(self._cur_style["figure.facecolor"])
        self._canvas = FigureCanvas(self._figure)
        self._canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self._fig_scroll = QScrollArea()
        self._fig_scroll.setWidgetResizable(False)
        self._fig_scroll.setWidget(self._canvas)
        self._fig_scroll.installEventFilter(self)
        self._toolbar = NavToolbar(self._canvas, fig_box)
        fig_layout.addWidget(self._toolbar)
        fig_layout.addWidget(self._fig_scroll, stretch=1)
        splitter.addWidget(fig_box)

        # Re-render with a small debounce after window resize so dragging the
        # window doesn't fire dozens of redraws.
        self._resize_timer = QTimer(self)
        self._resize_timer.setSingleShot(True)
        self._resize_timer.timeout.connect(self._render)

        splitter.setSizes([320, 900])

    # ── Plot-type dropdown structure ──────────────────────────────────────────

    def _add_plot_type_header(self, label: str):
        """Append a non-selectable, italicised section header to the plot-type combo."""
        self._plot_type.addItem(f"── {label} ──")
        idx = self._plot_type.count() - 1
        item = self._plot_type.model().item(idx)
        # Strip selectable/enabled flags so the header behaves like a label.
        item.setFlags(item.flags() & ~Qt.ItemIsSelectable & ~Qt.ItemIsEnabled)
        f = item.font()
        f.setItalic(True)
        item.setFont(f)

    def _set_plot_type_enabled(self, text: str, enabled: bool, tooltip: str = ""):
        """Toggle whether a particular plot-type entry can be selected."""
        idx = self._plot_type.findText(text)
        if idx < 0:
            return
        item = self._plot_type.model().item(idx)
        if item is None:
            return
        flags = item.flags()
        if enabled:
            flags |= Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            flags &= ~Qt.ItemIsEnabled & ~Qt.ItemIsSelectable
        item.setFlags(flags)
        item.setToolTip(tooltip)

    def _is_plot_type_enabled(self, text: str) -> bool:
        idx = self._plot_type.findText(text)
        if idx < 0:
            return False
        item = self._plot_type.model().item(idx)
        return bool(item and (item.flags() & Qt.ItemIsEnabled))

    def _update_plot_type_availability(self):
        """Disable plot-type entries that don't apply to the loaded data."""
        if self._df is None or self._df.empty:
            return
        n_groups = (int(self._df["Group"].nunique())
                    if "Group" in self._df.columns else 0)

        # Cross-group lines need at least two phase groups
        self._set_plot_type_enabled(
            "Effect lines across groups",
            enabled=(n_groups >= 2),
            tooltip="" if n_groups >= 2
                    else "Requires 2+ phase groups (only 1 present)."
        )

        # If the currently-selected plot just got disabled, fall back to a
        # sensible default so the user isn't stuck on a greyed-out entry.
        if not self._is_plot_type_enabled(self._plot_type.currentText()):
            for fallback in ("Volcano", "Top-N forest", "Beeswarm"):
                if self._is_plot_type_enabled(fallback):
                    idx = self._plot_type.findText(fallback)
                    if idx >= 0:
                        self._plot_type.setCurrentIndex(idx)
                        break

    def _on_plot_type_changed(self, _idx: int):
        # Some controls only matter for specific plot types — toggle enabled state.
        ptype = self._plot_type.currentText()
        is_pe    = ptype.startswith("Posterior eq.")
        is_lines = ptype == "Effect lines across groups"
        self._top_n.setEnabled(
            ptype == "Top-N forest" or is_pe or is_lines
        )
        self._rank_by_combo.setEnabled(ptype == "Top-N forest" or is_lines)
        self._show_labels.setEnabled(ptype == "Volcano" or ptype == "Beeswarm")
        # The scope combo applies to both volcano and beeswarm; it's only
        # meaningful when point labels are turned on.
        self._label_scope.setEnabled(
            (ptype == "Volcano" or ptype == "Beeswarm")
            and self._show_labels.isChecked()
        )
        # PE plots don't use Test_Family or Split_By filters.
        self._family_combo.setEnabled(not is_pe)
        self._splitby_combo.setEnabled(not is_pe)
        self._effect_min.setEnabled(not is_pe)
        # The cross-group lines plot integrates across groups; the per-group
        # filter would only collapse it to a single point per panel.
        if is_lines:
            self._group_combo.setEnabled(False)
        else:
            self._group_combo.setEnabled(True)

        # Warn if PE was selected but the draws/summary aren't on disk
        if is_pe and not self._pe_files_loaded():
            self._pe_availability_label.setText(
                "Posterior equivalence files not found. Enable PE in Contrast "
                "Selection and re-run BAM."
            )
            self._pe_availability_label.show()
        else:
            self._pe_availability_label.hide()
        self._render()

    def _pe_files_loaded(self) -> bool:
        return self._pe_draws_df is not None and self._pe_summary_df is not None

    # ── Theme helpers ────────────────────────────────────────────────────────

    def _notsig_color(self) -> str:
        """Color for non-significant points; theme-dependent for contrast."""
        return _NOTSIG_COLOR_DARK if self._theme == "dark" else _NOTSIG_COLOR_LIGHT

    def _on_theme_changed(self, _idx: int):
        mode = self._theme_combo.currentData()
        self._theme = "light" if mode == "light" else "dark"
        self._cur_style = _LIGHT_STYLE if self._theme == "light" else _DARK_STYLE
        # Match the canvas background to the new theme; otherwise the area
        # around the axes briefly shows the old theme until the next render.
        self._figure.patch.set_facecolor(self._cur_style["figure.facecolor"])
        self._render()

    def _pick_sig_color(self):
        """Open a color picker for the significant-point color."""
        initial = QColor(self._sig_color)
        chosen  = QColorDialog.getColor(
            initial, self, "Pick significant-point color"
        )
        if chosen.isValid():
            self._sig_color = chosen.name()  # "#rrggbb"
            self._update_sig_color_btn_style()
            self._render()

    def _update_sig_color_btn_style(self):
        self._restyle_color_btn(self._sig_color_btn, self._sig_color)

    def _pick_ref_color(self):
        """Open a color picker for the threshold/reference-line color."""
        initial = QColor(self._ref_color)
        chosen  = QColorDialog.getColor(
            initial, self, "Pick threshold-line color"
        )
        if chosen.isValid():
            self._ref_color = chosen.name()
            self._update_ref_color_btn_style()
            self._render()

    def _update_ref_color_btn_style(self):
        self._restyle_color_btn(self._ref_color_btn, self._ref_color)

    def _pick_ci_color(self):
        """Open a color picker for the confidence-interval ribbon color."""
        initial = QColor(self._ci_color)
        chosen  = QColorDialog.getColor(
            initial, self, "Pick confidence-interval color"
        )
        if chosen.isValid():
            self._ci_color = chosen.name()
            self._update_ci_color_btn_style()
            self._render()

    def _update_ci_color_btn_style(self):
        self._restyle_color_btn(self._ci_color_btn, self._ci_color)

    def _pick_density_color(self):
        """Open a color picker for the M-density histogram fill color."""
        initial = QColor(self._density_color)
        chosen  = QColorDialog.getColor(
            initial, self, "Pick M-density fill color"
        )
        if chosen.isValid():
            self._density_color = chosen.name()
            self._update_density_color_btn_style()
            self._render()

    def _update_density_color_btn_style(self):
        self._restyle_color_btn(self._density_color_btn, self._density_color)

    def _on_cmap_changed(self, _idx: int):
        key = self._heatmap_cmap_combo.currentData()
        if key:
            self._heatmap_cmap = key
            self._render()

    @staticmethod
    def _restyle_color_btn(btn, hex_color: str):
        """Repaint a color-picker button so it previews its current color.

        Auto-flips the text color between black/white based on perceived
        luminance so the hex label stays readable.
        """
        if btn is None:
            return
        btn.setText(hex_color)
        c = QColor(hex_color)
        luminance = 0.299 * c.red() + 0.587 * c.green() + 0.114 * c.blue()
        text_color = "#000000" if luminance > 140 else "#ffffff"
        btn.setStyleSheet(
            f"QPushButton {{"
            f" background-color: {hex_color};"
            f" color: {text_color};"
            f" border: 1px solid #555;"
            f"}}"
        )

    # ── Data loading ──────────────────────────────────────────────────────────

    def _try_load(self):
        if not self._output_dir:
            return
        path = os.path.join(self._output_dir, "master_results.csv")
        if not os.path.isfile(path):
            self._status_label.setText(
                "master_results.csv not found in output directory. Run BAM first."
            )
            self._df = None
            self._clear_axes("Run BAM to produce master_results.csv")
            return

        try:
            df = pd.read_csv(path)
        except Exception as exc:
            self._status_label.setText(f"Could not read CSV: {exc}")
            return

        required = {"Test_Family", "Group", "Split_By", "Tested_Level",
                    "estimate", "SE", "raw_pvalue", "adjusted_pvalue"}
        missing = required - set(df.columns)
        if missing:
            self._status_label.setText(
                f"master_results.csv missing required columns: {sorted(missing)}"
            )
            return

        # Add a "minus-log10 adjusted-p" once for reuse. Label annotation is
        # done lazily inside the render functions so the role-tag toggle takes
        # effect live without reloading the CSV.
        df = df.copy()
        df["neg_log10_fdr"] = -np.log10(df["adjusted_pvalue"].clip(lower=1e-300))
        df["abs_estimate"] = df["estimate"].abs()

        # Pick up the correction method R recorded, if available — keeps axis
        # labels in sync without needing the BAM widget to broadcast it.
        if "correction_method" in df.columns and not df["correction_method"].dropna().empty:
            method = str(df["correction_method"].dropna().iloc[0])
            if method.lower() in ("bh", "holm", "treebh",
                                  "holm-gatekeeping", "gmcp",
                                  "graphicalmcp"):
                self._correction_method = method

        # Tree-status filter visibility: only useful when the CSV actually
        # has non-NA tree_status (i.e. came from a Tree-strategy run).
        has_tree_info = (
            "tree_status" in df.columns and
            df["tree_status"].notna().any()
        )
        if hasattr(self, "_tree_status_combo"):
            self._tree_status_combo.setEnabled(has_tree_info)

        self._df = df
        self._status_label.setText(
            f"Loaded {len(df)} contrasts from {os.path.basename(path)}."
        )
        self._populate_filter_combos()
        self._try_load_pe()
        # Grey out plot types that don't apply to this dataset's shape (e.g.
        # cross-group lines when only one phase group is present).
        self._update_plot_type_availability()
        # Re-evaluate plot-type-dependent control state with new data on hand
        self._on_plot_type_changed(self._plot_type.currentIndex())

    def _try_load_pe(self):
        """Optional posterior-equivalence files. Absent when PE wasn't enabled."""
        if not self._output_dir:
            return
        draws_path   = os.path.join(self._output_dir, "posterior_equivalence_draws.csv")
        summary_path = os.path.join(self._output_dir, "posterior_equivalence_summary.csv")

        self._pe_draws_df   = None
        self._pe_summary_df = None
        self._pe_delta      = None

        if os.path.isfile(draws_path):
            try:
                self._pe_draws_df = pd.read_csv(draws_path)
            except Exception as exc:
                self._status_label.setText(
                    f"PE draws file present but unreadable: {exc}"
                )
        if os.path.isfile(summary_path):
            try:
                pe_sum = pd.read_csv(summary_path)
                self._pe_summary_df = pe_sum
                if "delta_log2" in pe_sum.columns and not pe_sum.empty:
                    self._pe_delta = float(pe_sum["delta_log2"].iloc[0])
            except Exception as exc:
                self._status_label.setText(
                    f"PE summary file present but unreadable: {exc}"
                )

    def _populate_filter_combos(self):
        """Refill the family / group / split-by combos from the loaded data."""
        if self._df is None:
            return
        self._suspend_render = True
        try:
            for combo, col, includes_all in (
                (self._family_combo,  "Test_Family", True),
                (self._group_combo,   "Group",       True),
                (self._splitby_combo, "Split_By",    True),
            ):
                current = combo.currentText()
                combo.blockSignals(True)
                combo.clear()
                if includes_all:
                    combo.addItem("(all)")
                values = sorted(
                    self._df[col].dropna().unique().tolist(),
                    key=lambda x: str(x),
                )
                combo.addItems([str(v) for v in values])
                # Restore prior selection if still present, otherwise default to (all)
                idx = combo.findText(current) if current else -1
                combo.setCurrentIndex(idx if idx >= 0 else 0)
                combo.blockSignals(False)
        finally:
            self._suspend_render = False

    # ── Label annotation ──────────────────────────────────────────────────────

    @staticmethod
    def _split_label(label: str) -> tuple[str, str | None]:
        """emmeans labels for pairs use ' - ' as a separator."""
        parts = re.split(r"\s+-\s+", str(label), maxsplit=1)
        if len(parts) == 2:
            return parts[0].strip(), parts[1].strip()
        return str(label).strip(), None

    def _annotate_one(self, cond_label: str, split_by: str | None = None) -> str:
        """Adds [RC] / [EXP] / etc. role tag if known and the user wants it.

        The role map is keyed by *full* condition strings like "WT+DMSO". For
        per-variable contrast labels (e.g. just "KO" out of Genotype_Effect),
        the level alone won't match a role key. ``split_by`` carries the
        other-variable context (e.g. "DMSO") so we can rebuild the full
        condition and resolve a role on those rows too.
        """
        if not self._roles:
            return cond_label
        # Layout toggle: user may prefer raw condition names
        if hasattr(self, "_show_role_tags") and not self._show_role_tags.isChecked():
            return cond_label

        # Direct match: emmeans returns levels space-joined for combined
        # factors; the role map keys are "+"-joined.
        cond_plus = cond_label.replace(" ", "+")
        role = self._roles.get(cond_plus)

        # Reconstruction match: per-variable contrasts (e.g. cond_label="KO",
        # split_by="DMSO") need to be re-glued into "KO+DMSO" or "DMSO+KO".
        # Try both orders since we don't know which variable came first in
        # the design.
        if role is None and split_by:
            # Split_By can be " + " joined when multiple split variables exist.
            split_plus = (
                str(split_by).replace(" + ", "+").replace(" ", "+")
            )
            for combined in (f"{cond_plus}+{split_plus}",
                             f"{split_plus}+{cond_plus}"):
                role = self._roles.get(combined)
                if role:
                    break

        if role:
            tag = self._role_abbrev_map.get(role, role)
            return f"{cond_label} [{tag}]"
        return cond_label

    def _annotate_label(self, label: str, group=None,
                        split_by: str | None = None) -> str:
        lhs, rhs = self._split_label(label)
        if rhs is None:
            base = self._annotate_one(lhs, split_by=split_by)
        else:
            base = (f"{self._annotate_one(lhs, split_by=split_by)} - "
                    f"{self._annotate_one(rhs, split_by=split_by)}")

        # Optional phase-group prefix. The checkbox lives in the Layout group
        # and is unrelated to which plot type is active — callers that don't
        # have a group context simply omit the argument.
        if (group is not None
                and hasattr(self, "_show_group_in_label")
                and self._show_group_in_label.isChecked()):
            try:
                gtxt = f"G{int(group)}"
            except (TypeError, ValueError):
                gtxt = f"G{group}"
            base = f"[{gtxt}] {base}"
        return base

    # ── Rendering ─────────────────────────────────────────────────────────────

    # ── Canvas sizing ────────────────────────────────────────────────────────

    def eventFilter(self, obj, event):
        """Debounced re-render when the canvas scroll area is resized.

        This keeps the figure width matched to the viewport width so plots
        don't end up wedged into a too-small or too-large area when the user
        resizes the window.
        """
        if obj is getattr(self, "_fig_scroll", None) and event.type() == QEvent.Resize:
            self._resize_timer.start(150)
        return super().eventFilter(obj, event)

    def _sync_canvas_size(self, content_height_in: float | None = None):
        """Resize the figure + canvas to match the scroll viewport.

        Width always tracks the viewport. Height fits the viewport unless
        ``content_height_in`` is provided AND larger than the viewport — in
        which case the canvas grows tall and the scroll area shows a vertical
        scrollbar. Called at the start of every render so figure size is
        always pinned by the current plot type, never inherited from the
        previous one.
        """
        if not hasattr(self, "_fig_scroll"):
            return
        # Trim a few pixels to account for scrollbar / viewport frame
        vp_w = max(360, self._fig_scroll.viewport().width()  - 4)
        vp_h = max(280, self._fig_scroll.viewport().height() - 4)
        dpi  = self._figure.dpi

        w_in = vp_w / dpi
        if content_height_in is None:
            h_in = vp_h / dpi
        else:
            h_in = max(vp_h / dpi, float(content_height_in))

        w_px = int(round(w_in * dpi))
        h_px = int(round(h_in * dpi))
        self._canvas.setFixedSize(w_px, h_px)
        # forward=False because we just set the canvas size ourselves;
        # forward=True would cause Qt-layout fights and was the source of
        # the inconsistent-sizing symptom.
        self._figure.set_size_inches(w_in, h_in, forward=False)

    def _clear_axes(self, message: str | None = None):
        # Empty-state messages should fit the viewport, not inherit a tall size
        # from a previous many-panel render.
        self._sync_canvas_size(None)
        self._figure.clear()
        ax = self._figure.add_subplot(111)
        ax.set_facecolor(self._cur_style["axes.facecolor"])
        if message:
            ax.text(0.5, 0.5, message,
                    ha="center", va="center",
                    color=self._cur_style["text.color"],
                    fontsize=12, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        self._canvas.draw_idle()

    def _filtered_df(self) -> pd.DataFrame | None:
        if self._df is None:
            return None
        df = self._df

        fam = self._family_combo.currentText()
        if fam and fam != "(all)":
            df = df[df["Test_Family"].astype(str) == fam]

        # The cross-group lines plot integrates across groups, so the per-group
        # selector is ignored even if a value lingers from another plot type.
        # The effect-min filter is also skipped — a missing point would leave
        # gaps in the line trajectory.
        ptype = self._plot_type.currentText() if hasattr(self, "_plot_type") else ""
        is_lines = ptype == "Effect lines across groups"

        grp = self._group_combo.currentText()
        if grp and grp != "(all)" and not is_lines:
            df = df[df["Group"].astype(str) == grp]

        sb = self._splitby_combo.currentText()
        if sb and sb != "(all)":
            df = df[df["Split_By"].astype(str) == sb]

        # Tree-status filter — only acts when the CSV has tree_status info.
        if "tree_status" in df.columns and hasattr(self, "_tree_status_combo"):
            scope = self._tree_status_combo.currentData()
            if scope == "tested":
                df = df[df["tree_status"].astype(str) == "tested"]
            elif scope == "tested_plus_gate_failed":
                df = df[df["tree_status"].astype(str).isin(["tested", "gate_failed"])]
            elif scope == "not_reached":
                df = df[df["tree_status"].astype(str) == "not_reached"]
            elif scope == "off_tree":
                df = df[df["tree_status"].astype(str) == "off_tree"]
            # "all" or any other value → no filter

        if not is_lines:
            df = df[df["abs_estimate"] >= self._effect_min.value()]
        return df

    def _render(self):
        if self._suspend_render:
            return

        # Apply dark style to all matplotlib operations in this redraw
        with matplotlib.rc_context(rc=self._cur_style):
            matplotlib.rcParams["font.size"] = self._font_size.value()
            ptype = self._plot_type.currentText()
            self._figure.clear()

            if ptype.startswith("Posterior eq."):
                if not self._pe_files_loaded():
                    self._clear_axes(
                        "Posterior equivalence files not found.\n"
                        "Enable PE in Contrast Selection and re-run BAM."
                    )
                    return
                if ptype == "Posterior eq.: M density":
                    self._render_pe_density()
                elif ptype == "Posterior eq.: Probability of equivalence":
                    self._render_pe_heatmap()
                else:
                    self._clear_axes(f"Unknown plot type: {ptype}")
                    return
            else:
                df = self._filtered_df()
                if df is None or df.empty:
                    self._clear_axes("No data after current filters.")
                    return

                if ptype == "Volcano":
                    self._render_volcano(df)
                elif ptype == "Top-N forest":
                    self._render_topn_forest(df)
                elif ptype == "Beeswarm":
                    self._render_beeswarm(df)
                elif ptype == "Effect lines across groups":
                    self._render_cross_group_lines(df)
                else:
                    self._clear_axes(f"Unknown plot type: {ptype}")
                    return

            # Reserve a thin strip at the top of the figure for suptitles
            # so they don't collide with the first row of subplot titles.
            try:
                self._figure.tight_layout(rect=[0, 0, 1, 0.99])
            except Exception:
                try:
                    self._figure.tight_layout()
                except Exception:
                    pass
            self._canvas.draw_idle()

    # ---- Volcano -----------------------------------------------------------------

    def _render_volcano(self, df: pd.DataFrame):
        # Volcano fits the viewport — no scroll, no inherited tall size.
        self._sync_canvas_size(None)

        sig_thresh = self._sig_thresh.value()
        grp_sel    = self._group_combo.currentText()
        fam_sel    = self._family_combo.currentText()

        # Facet when no specific family is selected
        if fam_sel == "(all)":
            families = sorted(df["Test_Family"].dropna().unique().tolist(),
                              key=lambda x: str(x))
        else:
            families = [fam_sel]

        n_panels = len(families)
        ncols = min(2, n_panels)
        nrows = int(np.ceil(n_panels / ncols))

        axes = self._figure.subplots(nrows, ncols, squeeze=False)
        for ax in axes.flatten():
            ax.set_visible(False)

        # Shared symmetric x-limit across all panels — centers every facet on 0
        # using the global max |log2 FC| so the +1 / -1 reference lines stay in
        # the same visual position across families.
        x_pad = 1.05
        global_max_abs = float(df["abs_estimate"].max()) if not df.empty else 1.0
        # Always keep at least ±1.2 visible so the ±1 reference lines aren't
        # flush against the axis edges when all effects are tiny.
        x_extent = max(global_max_abs * x_pad, 1.2)

        for i, fam in enumerate(families):
            r, c = divmod(i, ncols)
            ax = axes[r][c]
            ax.set_visible(True)
            ax.set_facecolor(self._cur_style["axes.facecolor"])

            sub = df[df["Test_Family"].astype(str) == str(fam)]
            if sub.empty:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes,
                        color=self._cur_style["text.color"])
                ax.set_xlim(-x_extent, x_extent)
                continue

            is_sig    = sub["adjusted_pvalue"] < sig_thresh
            sig_df    = sub[is_sig]
            notsig_df = sub[~is_sig]

            ax.scatter(notsig_df["estimate"], notsig_df["neg_log10_fdr"],
                       s=22, c=self._notsig_color(), alpha=0.55, edgecolors="none",
                       label=f"ns (n={len(notsig_df)})")
            ax.scatter(sig_df["estimate"], sig_df["neg_log10_fdr"],
                       s=32, c=self._sig_color, alpha=0.85, edgecolors="none",
                       label=f"{self._correction_label()}<{sig_thresh:g} "
                             f"(n={len(sig_df)})")

            # Reference lines: 2× fold change in either direction (log_2 = ±1).
            ax.axvline( 1, color=self._ref_color, linestyle="--",
                       linewidth=0.8, alpha=0.7)
            ax.axvline(-1, color=self._ref_color, linestyle="--",
                       linewidth=0.8, alpha=0.7)
            ax.axhline(-np.log10(sig_thresh), color=self._ref_color,
                       linestyle=":", linewidth=0.8, alpha=0.7)

            # Symmetric x-axis around 0 — same range in every panel.
            ax.set_xlim(-x_extent, x_extent)

            ax.set_xlabel(r"log$_2$ fold change (vs reference)")
            ax.set_ylabel(rf"$-\log_{{10}}$({self._correction_label()} "
                          rf"p-value)")

            panel_title = self._pretty(fam)
            if grp_sel and grp_sel != "(all)":
                panel_title += f"  ·  Group {grp_sel}"
            ax.set_title(panel_title,
                         fontsize=max(8, self._font_size.value() - 1))

            if self._show_labels.isChecked():
                scope = self._label_scope.currentData()
                if scope == "all":
                    # Label every point, capped so a 1000-pair volcano
                    # doesn't turn into a wall of overlapping text. Rank by
                    # -log10(FDR) so the most significant are guaranteed
                    # to appear if the cap kicks in.
                    cap = 200
                    label_df = sub.nlargest(min(cap, len(sub)), "neg_log10_fdr")
                else:  # "sig" → significant-only (historical behaviour)
                    label_df = sig_df.nlargest(min(20, len(sig_df)), "neg_log10_fdr")
                for _, row in label_df.iterrows():
                    ax.annotate(
                        self._annotate_label(row["Tested_Level"],
                                             group=row.get("Group"),
                                             split_by=row.get("Split_By")),
                        xy=(row["estimate"], row["neg_log10_fdr"]),
                        xytext=(4, 2), textcoords="offset points",
                        fontsize=max(6, self._font_size.value() - 2),
                        color=self._cur_style["text.color"], alpha=0.85,
                    )

            ax.legend(loc="best", framealpha=0.5,
                      facecolor=self._cur_style["axes.facecolor"],
                      edgecolor=self._cur_style["axes.edgecolor"],
                      labelcolor=self._cur_style["text.color"],
                      fontsize=max(7, self._font_size.value() - 2))

        # Suptitle reflects the global filter context
        sup = "Volcano"
        if grp_sel and grp_sel != "(all)":
            sup += f"  ·  Group {grp_sel}"
        self._figure.suptitle(sup, color=self._cur_style["text.color"])

    # ---- Top-N forest ------------------------------------------------------------

    def _render_topn_forest(self, df: pd.DataFrame):
        rank_by    = self._rank_by_combo.currentData()
        n_req      = int(self._top_n.value())
        sig_thresh = self._sig_thresh.value()
        rank_col   = "neg_log10_fdr" if rank_by == "fdr" else "abs_estimate"
        fam_sel    = self._family_combo.currentText()

        # Facet when no specific family is selected — each panel gets its own
        # top-N within its family (effect-size scales differ across families).
        if fam_sel == "(all)":
            families = sorted(df["Test_Family"].dropna().unique().tolist(),
                              key=lambda x: str(x))
        else:
            families = [fam_sel]

        panels = []
        for fam in families:
            sub = df[df["Test_Family"].astype(str) == str(fam)]
            if sub.empty:
                continue
            ranked = (sub.nlargest(n_req, rank_col)
                         .sort_values("estimate", ascending=True))
            if ranked.empty:
                continue
            panels.append((str(fam), ranked))

        if not panels:
            self._clear_axes("No rows match the current filters.")
            return

        ncols = min(2, len(panels))
        nrows = int(np.ceil(len(panels) / ncols))

        # Height scales with the tallest panel
        max_rows = max(len(r) for _, r in panels)
        per_panel_h = max(2.5, min(0.3 * max_rows + 1.0, 18.0))
        self._sync_canvas_size(per_panel_h * nrows)

        axes = self._figure.subplots(nrows, ncols, squeeze=False)
        for ax in axes.flatten():
            ax.set_visible(False)

        for i, (fam, ranked) in enumerate(panels):
            r, c = divmod(i, ncols)
            ax = axes[r][c]
            ax.set_visible(True)
            ax.set_facecolor(self._cur_style["axes.facecolor"])

            n_rows_p = len(ranked)
            ys      = np.arange(n_rows_p)
            ci_half = 1.96 * ranked["SE"].values
            x       = ranked["estimate"].values
            sig_mask = (ranked["adjusted_pvalue"].values < sig_thresh)
            colors   = np.where(sig_mask, self._sig_color, self._notsig_color())

            ax.hlines(ys, x - ci_half, x + ci_half,
                      colors=colors, linewidth=1.6, alpha=0.85)
            ax.scatter(x, ys, c=colors, s=36, zorder=3)
            ax.axvline(0, color=self._ref_color, linestyle="--",
                       linewidth=0.8, alpha=0.7)

            labels = [
                self._annotate_label(t, group=g, split_by=s)
                for t, g, s in zip(
                    ranked["Tested_Level"].astype(str),
                    ranked["Group"] if "Group" in ranked.columns
                                    else [None] * len(ranked),
                    ranked["Split_By"] if "Split_By" in ranked.columns
                                       else [None] * len(ranked),
                )
            ]
            ax.set_yticks(ys)
            ax.set_yticklabels(labels,
                               fontsize=max(6, self._font_size.value() - 1))
            ax.set_xlabel(r"log$_2$ fold change (95% CI)")

            panel_title = self._pretty(fam)
            if len(panels) > 1:
                panel_title += f"  ·  top {n_rows_p}"
            ax.set_title(panel_title,
                         fontsize=max(8, self._font_size.value() - 1))

            ax.invert_yaxis()
            # Small vertical padding so the topmost and bottommost rows
            # aren't crammed against the axes border.
            ax.margins(y=0.03)

        rank_label = (f"{self._correction_label()} p-value"
                      if rank_by == "fdr" else "|effect size|")
        sup = f"Top-{n_req} forest  ·  ranked by {rank_label}"
        grp_sel = self._group_combo.currentText()
        if grp_sel and grp_sel != "(all)":
            sup += f"  ·  Group {grp_sel}"
        self._figure.suptitle(sup, color=self._cur_style["text.color"])

    # ---- Beeswarm ----------------------------------------------------------------

    def _render_beeswarm(self, df: pd.DataFrame):
        # Beeswarm fits the viewport; reset canvas so a prior tall plot doesn't
        # leave the figure oversized.
        self._sync_canvas_size(None)

        # Group along x by Test_Family (or by Group if a single family is active)
        if self._family_combo.currentText() == "(all)":
            group_col = "Test_Family"
        else:
            group_col = "Group"

        ax = self._figure.add_subplot(111)
        ax.set_facecolor(self._cur_style["axes.facecolor"])

        sig_thresh   = self._sig_thresh.value()
        show_labels  = self._show_labels.isChecked()
        label_scope  = self._label_scope.currentData() if show_labels else None
        cats = sorted(df[group_col].dropna().unique().tolist(), key=lambda x: str(x))
        if not cats:
            self._clear_axes("No categories to plot.")
            return

        # Per-category label caps. "Significant only" stays terse so the
        # callouts focus attention on the few hits; "All points" goes wider
        # but still bounds the total to avoid total overlap.
        cap_sig = 10
        cap_all = 30

        # Deterministic horizontal jitter for swarm-like effect
        rng = np.random.default_rng(seed=0)
        for i, cat in enumerate(cats):
            sub = df[df[group_col] == cat]
            if sub.empty:
                continue

            est = sub["estimate"].values
            fdr = sub["adjusted_pvalue"].values
            xs  = i + rng.uniform(-0.25, 0.25, size=len(est))
            sig_mask = (fdr < sig_thresh)

            ax.scatter(xs[~sig_mask], est[~sig_mask],
                       s=20, c=self._notsig_color(), alpha=0.55, edgecolors="none")
            ax.scatter(xs[sig_mask], est[sig_mask],
                       s=30, c=self._sig_color, alpha=0.9, edgecolors="none")

            # Median tick — anchor for the distribution at a glance
            med = np.median(est)
            ax.hlines(med, i - 0.32, i + 0.32, colors=self._ref_color,
                      linewidth=2.0, alpha=0.85)
            # Mean tick (smaller, dashed) — gives sense of skew vs the median
            mean_val = float(np.mean(est))
            ax.hlines(mean_val, i - 0.22, i + 0.22, colors="#80cbc4",
                      linewidth=1.0, alpha=0.85, linestyles="dashed")

            # Annotate dots according to the user-chosen scope.
            if show_labels:
                sub_local = sub.assign(_x=xs, _abs=np.abs(est))
                if label_scope == "all":
                    label_rows = sub_local.sort_values(
                        "_abs", ascending=False
                    ).head(cap_all)
                else:  # "sig"
                    label_rows = (
                        sub_local[sub_local["adjusted_pvalue"] < sig_thresh]
                        .sort_values("_abs", ascending=False)
                        .head(cap_sig)
                    )
                for _, row in label_rows.iterrows():
                    ax.annotate(
                        self._annotate_label(row["Tested_Level"],
                                             group=row.get("Group"),
                                             split_by=row.get("Split_By")),
                        xy=(row["_x"], row["estimate"]),
                        xytext=(4, 2), textcoords="offset points",
                        fontsize=max(6, self._font_size.value() - 3),
                        color=self._cur_style["text.color"], alpha=0.85,
                    )

        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels([self._pretty(c) for c in cats],
                           rotation=30, ha="right")
        ax.axhline(0, color=self._ref_color, linestyle="--", linewidth=0.8, alpha=0.7)
        ax.set_ylabel(r"log$_2$ fold change")
        ax.set_xlabel(self._pretty(group_col))

        # Small legend explaining the central tick marks
        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], color=self._ref_color, lw=2, label="median"),
            Line2D([0], [0], color="#80cbc4", lw=1, linestyle="--", label="mean"),
            Line2D([0], [0], marker="o", color="none",
                   markerfacecolor=self._sig_color, markersize=7,
                   label=f"{self._correction_label()}<{sig_thresh:g}"),
            Line2D([0], [0], marker="o", color="none",
                   markerfacecolor=self._notsig_color(), markersize=7, label="ns"),
        ]
        ax.legend(handles=legend_handles, loc="best", framealpha=0.5,
                  facecolor=self._cur_style["axes.facecolor"],
                  edgecolor=self._cur_style["axes.edgecolor"],
                  labelcolor=self._cur_style["text.color"],
                  fontsize=max(7, self._font_size.value() - 2))

        ax.set_title(f"Effect-size distribution by {self._pretty(group_col)}")

    # ---- Effect lines across groups ---------------------------------------------

    def _render_cross_group_lines(self, df: pd.DataFrame):
        """Small-multiples line plot of one contrast's effect-size trajectory
        across phase groups, with a 95% CI ribbon and significance-colored points.

        Replaces the old "rescue assessment" R plot but works for any contrast
        — pick a Test_Family with a single dropdown to narrow the set, and use
        Top-N + the rank-by dropdown to limit how many panels are drawn.
        """
        # Need at least two groups in the data for a line plot to be meaningful.
        if df["Group"].nunique() < 2:
            self._clear_axes(
                "This plot needs multiple phase groups.\n"
                "Only one group is present in master_results.csv."
            )
            return

        sig_thresh = self._sig_thresh.value()
        n_req      = int(self._top_n.value())
        rank_by    = self._rank_by_combo.currentData()

        # Per-variable contrasts (e.g. "KO - WT" out of Genotype_Effect) have
        # multiple rows per Group — one per Split_By context (e.g. one for
        # DMSO, one for compound). Grouping by Tested_Level alone stitches
        # those into a single zig-zagging line. Use the (Tested_Level,
        # Split_By) pair as the trajectory identifier so each context
        # becomes its own smooth panel.
        df = df.copy()
        df["_traj_key"] = (df["Tested_Level"].astype(str) + "  |  "
                           + df["Split_By"].astype(str))

        # Compute a per-trajectory score across groups, then keep the top N.
        if rank_by == "fdr":
            # Most significant in ANY single group = lowest min FDR
            scores = df.groupby("_traj_key")["adjusted_pvalue"].min()
            picked = scores.nsmallest(n_req).index.tolist()
        else:
            # Largest effect magnitude in ANY single group
            scores = df.groupby("_traj_key")["abs_estimate"].max()
            picked = scores.nlargest(n_req).index.tolist()

        if not picked:
            self._clear_axes("No contrasts to plot.")
            return

        sub = df[df["_traj_key"].astype(str).isin(picked)]
        if sub.empty:
            self._clear_axes("No data to plot for the selected contrasts.")
            return

        # ── Layout: per-variable effects in the left column, full-interaction
        # contrasts in the right column. Within the per-variable column,
        # families are ordered by their first appearance in the CSV — which
        # matches the order R wrote them, which in turn matches the design
        # variable order (Genotype before Drug, for the typical setup).
        family_first = {}
        for fam in df["Test_Family"].astype(str):
            if fam not in family_first:
                family_first[fam] = len(family_first)

        def _family_of(key: str) -> str:
            rows_k = sub[sub["_traj_key"] == key]
            return str(rows_k["Test_Family"].iloc[0]) if not rows_k.empty else ""

        left_keys, right_keys = [], []
        for key in picked:
            fam = _family_of(key)
            if fam == "Full_Interaction":
                right_keys.append(key)
            else:
                left_keys.append((family_first.get(fam, 999), key))
        # Stable sort by family-appearance order; rank-within-family preserved
        left_keys.sort(key=lambda x: x[0])
        left_keys = [k for _, k in left_keys]

        # Decide grid shape and where each key lands. If one column is empty
        # collapse to a single-column layout so the remaining column doesn't
        # get squished by an absent sibling.
        keys_in_grid: list[tuple[str, int, int]] = []  # (key, row, col)
        if left_keys and right_keys:
            nrows = max(len(left_keys), len(right_keys))
            ncols = 2
            keys_in_grid += [(k, i, 0) for i, k in enumerate(left_keys)]
            keys_in_grid += [(k, i, 1) for i, k in enumerate(right_keys)]
        elif left_keys:
            ncols = 1
            nrows = len(left_keys)
            keys_in_grid += [(k, i, 0) for i, k in enumerate(left_keys)]
        else:
            ncols = 1
            nrows = len(right_keys)
            keys_in_grid += [(k, i, 0) for i, k in enumerate(right_keys)]

        # Height scales with the tallest column
        per_row_h = 2.5
        self._sync_canvas_size(max(3.5, min(per_row_h * nrows + 1.2, 36.0)))

        axes = self._figure.subplots(nrows, ncols, squeeze=False, sharex=True)
        for ax in axes.flatten():
            ax.set_visible(False)

        all_groups_sorted = sorted(df["Group"].dropna().unique().tolist())

        # Last row used in each column — needed so the x-label appears at the
        # bottom of EACH column rather than only at nrows-1 (which would skip
        # the shorter column's bottom panel when the two columns differ in
        # length).
        last_row_in_col: dict[int, int] = {}
        for _key, _r, _c in keys_in_grid:
            last_row_in_col[_c] = max(last_row_in_col.get(_c, -1), _r)

        for traj_key, r, c in keys_in_grid:
            ax = axes[r][c]
            ax.set_visible(True)
            ax.set_facecolor(self._cur_style["axes.facecolor"])

            rows = sub[sub["_traj_key"] == traj_key].sort_values("Group")
            if rows.empty:
                continue

            x      = rows["Group"].to_numpy()
            y      = rows["estimate"].to_numpy()
            ci_lo  = y - 1.96 * rows["SE"].to_numpy()
            ci_hi  = y + 1.96 * rows["SE"].to_numpy()
            sig    = (rows["adjusted_pvalue"].to_numpy() < sig_thresh)
            colors = np.where(sig, self._sig_color, self._notsig_color())

            ax.fill_between(x, ci_lo, ci_hi, color=self._ci_color, alpha=0.18,
                            linewidth=0)
            ax.plot(x, y, color="#cfd8dc", linewidth=1.3, alpha=0.85)
            ax.scatter(x, y, c=colors, s=42, zorder=3, edgecolors="none")
            ax.axhline(0, color=self._ref_color, linestyle="--",
                       linewidth=0.8, alpha=0.7)

            # Pull contrast + split_by back out of the key for the panel
            # title. Skip the split label when it's the placeholder "None"
            # used by Full_Interaction rows so those titles stay clean.
            contrast = rows["Tested_Level"].iloc[0]
            split_by = rows["Split_By"].iloc[0]
            title = self._annotate_label(contrast, split_by=split_by)
            if str(split_by) not in ("None", "All", "nan"):
                title += f"  ·  in {self._pretty(split_by)}"
            ax.set_title(title, fontsize=max(8, self._font_size.value() - 1))
            # X-label on the LAST USED row of THIS column (the two columns
            # may have different heights); y-label only on column 0.
            if r == last_row_in_col.get(c, -1):
                ax.set_xlabel("Group (Experimental Phase)")
            if c == 0:
                ax.set_ylabel(r"Effect size (log$_2$ FC, ±95% CI)")

            ax.set_xticks(all_groups_sorted)

        # Legend lives on the figure so it doesn't crowd any one panel
        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], marker="o", color="none",
                   markerfacecolor=self._notsig_color(), markersize=7, label="ns"),
            Line2D([0], [0], marker="o", color="none",
                   markerfacecolor=self._sig_color, markersize=7,
                   label=f"{self._correction_label()} < {sig_thresh:g}"),
        ]
        self._figure.legend(
            handles=legend_handles, loc="lower center", ncol=2,
            facecolor=self._cur_style["axes.facecolor"],
            edgecolor=self._cur_style["axes.edgecolor"],
            labelcolor=self._cur_style["text.color"],
            fontsize=max(8, self._font_size.value() - 1),
            bbox_to_anchor=(0.5, 0.0),
        )

        # The top-N count and rank-by criterion already live in the Filters
        # panel, so the suptitle stays minimal — just what the plot shows.
        sup = "Effect-size trajectory across phase groups"
        fam_sel = self._family_combo.currentText()
        if fam_sel and fam_sel != "(all)":
            sup += f"  ·  {self._pretty(fam_sel)}"
        self._figure.suptitle(sup, color=self._cur_style["text.color"])

    # ---- Posterior equivalence: M density ---------------------------------------

    def _pe_filtered_pairs(self) -> list[str]:
        """Top-N pairs (by Pr_equiv descending), respecting the current group filter."""
        if self._pe_summary_df is None:
            return []
        sub = self._pe_summary_df
        grp = self._group_combo.currentText()
        if grp and grp != "(all)":
            sub = sub[sub["Group"].astype(str) == grp]
        if sub.empty:
            return []
        n = max(1, int(self._top_n.value()))
        ranked = sub.sort_values("Pr_equiv", ascending=False).head(n)
        return ranked["pair"].astype(str).tolist()

    def _render_pe_density(self):
        draws = self._pe_draws_df
        if draws is None or draws.empty:
            self._clear_axes("No posterior draws available.")
            return

        pairs = self._pe_filtered_pairs()
        grp_sel = self._group_combo.currentText()

        sub = draws[draws["pair"].astype(str).isin(pairs)]
        if grp_sel and grp_sel != "(all)":
            sub = sub[sub["Group"].astype(str) == grp_sel]
        if sub.empty:
            self._clear_axes("No PE draws match the current filters.")
            return

        # One small subplot per (pair, group) combination — facet-wrap style.
        keys = list(sub.groupby(["pair", "Group"]).groups.keys())
        n_panels = len(keys)
        ncols = min(3, n_panels)
        nrows = int(np.ceil(n_panels / ncols))

        # Scale figure height with panels
        self._sync_canvas_size(max(4.0, min(2.2 * nrows + 1.2, 36.0)))

        axes = self._figure.subplots(nrows, ncols, squeeze=False)
        delta = self._pe_delta

        # Annotated label cache so each panel shows role tags
        for ax in axes.flatten():
            ax.set_visible(False)

        for i, (pair, group) in enumerate(keys):
            r, c = divmod(i, ncols)
            ax = axes[r][c]
            ax.set_visible(True)
            ax.set_facecolor(self._cur_style["axes.facecolor"])

            vals = sub[(sub["pair"] == pair) & (sub["Group"] == group)]["M_log2"].values
            if vals.size < 5:
                ax.text(0.5, 0.5, "n < 5", ha="center", va="center",
                        transform=ax.transAxes, color=self._cur_style["text.color"])
                ax.set_xticks([]); ax.set_yticks([])
                continue

            # Histogram with KDE-like smoothing via matplotlib density=True
            ax.hist(vals, bins=40, color=self._density_color, alpha=0.55,
                    density=True, edgecolor="none")
            if delta is not None:
                # Equivalence-margin line follows the user-pickable
                # threshold-line color so the choice is consistent across
                # every plot type that draws thresholds.
                ax.axvline(delta, color=self._ref_color, linestyle="--",
                           linewidth=1.2, alpha=0.9)
                # Annotate Pr(M < delta) in-panel. The current text-color
                # already contrasts with the panel background, so no bbox.
                pr_equiv = float(np.mean(vals < delta))
                ax.text(0.97, 0.93, f"Pr(M<δ)={pr_equiv:.2f}",
                        ha="right", va="top", transform=ax.transAxes,
                        color=self._cur_style["text.color"],
                        fontsize=max(7, self._font_size.value() - 2))

            ax.set_title(f"{self._annotate_pe_pair(pair)}  ·  G{group}",
                         fontsize=max(7, self._font_size.value() - 2))
            ax.set_xlabel(r"M (log$_2$ FC)",
                          fontsize=max(7, self._font_size.value() - 2))
            ax.set_ylabel("density",
                          fontsize=max(7, self._font_size.value() - 2))

        suptitle = "Posterior densities of M (max absolute trajectory difference)"
        if delta is not None:
            suptitle += f"  ·  δ = {delta:.2f} log₂"
        # No explicit y — matplotlib's default suptitle position works fine
        # now that the dispatcher's tight_layout reserves a small top sliver.
        self._figure.suptitle(suptitle, color=self._cur_style["text.color"])

    # ---- Posterior equivalence: Pr(M < δ) heatmap -------------------------------

    def _render_pe_heatmap(self):
        summary = self._pe_summary_df
        if summary is None or summary.empty:
            self._clear_axes("No PE summary available.")
            return

        sub = summary
        grp_sel = self._group_combo.currentText()
        if grp_sel and grp_sel != "(all)":
            sub = sub[sub["Group"].astype(str) == grp_sel]
        if sub.empty:
            self._clear_axes("No PE rows match the current filters.")
            return

        # Top-N rows by Pr_equiv descending so the most-equivalent pairs sit at top.
        n = max(1, int(self._top_n.value()))
        pairs_ranked = (
            sub.groupby("pair")["Pr_equiv"]
               .mean()
               .sort_values(ascending=False)
               .head(n)
               .index.tolist()
        )
        sub = sub[sub["pair"].isin(pairs_ranked)]

        pivot = (
            sub.pivot_table(index="pair", columns="Group", values="Pr_equiv")
               .reindex(pairs_ranked)
        )
        if pivot.empty:
            self._clear_axes("No PE matrix to plot.")
            return

        # Height scales with rows
        self._sync_canvas_size(max(4.0, min(0.3 * len(pivot) + 1.5, 36.0)))

        ax = self._figure.add_subplot(111)
        ax.set_facecolor(self._cur_style["axes.facecolor"])
        # Suppress the dark-style major grid — it shows up as dotted lines
        # piercing the heatmap cells and is purely visual noise here.
        ax.grid(False)
        ax.set_axisbelow(False)
        for ax2 in [ax]:
            ax2.xaxis.grid(False)
            ax2.yaxis.grid(False)
        # Optional cmap reversal — every matplotlib cmap has an "_r" twin.
        cmap_name = self._heatmap_cmap
        if self._invert_cmap.isChecked():
            cmap_name = cmap_name + "_r"
        cmap_obj = matplotlib.colormaps.get_cmap(cmap_name)

        im = ax.imshow(pivot.values, aspect="auto", cmap=cmap_obj,
                       vmin=0.0, vmax=1.0, origin="upper")

        # Annotate each cell with the probability. Pick the text color from
        # the underlying cell's luminance so it stays readable regardless of
        # which cmap is active (and whether the user inverted it). Rec. 601
        # weighting: Y = 0.299 R + 0.587 G + 0.114 B; the 0.6 threshold
        # leaves a small margin so the swap happens before contrast gets
        # uncomfortable.
        n_rows, n_cols = pivot.values.shape
        for r in range(n_rows):
            for c in range(n_cols):
                val = pivot.values[r, c]
                if np.isnan(val):
                    continue
                rgba = cmap_obj(float(val))   # rgba in [0, 1]
                lum  = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
                text_color = "black" if lum > 0.6 else "white"
                ax.text(c, r, f"{val:.2f}", ha="center", va="center",
                        color=text_color,
                        fontsize=max(7, self._font_size.value() - 2))

        ax.set_xticks(range(n_cols))
        ax.set_xticklabels([f"G{g}" for g in pivot.columns])
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels([self._annotate_pe_pair(p) for p in pivot.index],
                           fontsize=max(7, self._font_size.value() - 1))

        cbar = self._figure.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
        cbar.set_label("Probability of equivalence",
                       color=self._cur_style["text.color"])
        cbar.ax.yaxis.set_tick_params(color=self._cur_style["text.color"])
        for t in cbar.ax.get_yticklabels():
            t.set_color(self._cur_style["text.color"])

        title = "Probability of equivalence by pair × phase group"
        if self._pe_delta is not None:
            title += (f"  ·  equivalence margin δ = "
                      f"{self._pe_delta:.2f} log₂ fold change")
        ax.set_title(title)
        ax.set_xlabel("Phase group")

    def _annotate_pe_pair(self, pair_label: str) -> str:
        """PE pair strings look like 'WT+DMSO vs WT+DrugX'. Annotate both sides."""
        parts = re.split(r"\s+vs\s+", str(pair_label), maxsplit=1)
        if len(parts) != 2:
            return pair_label
        lhs, rhs = parts
        return f"{self._annotate_one(lhs)} vs {self._annotate_one(rhs)}"

    # ── Export ────────────────────────────────────────────────────────────────

    def _on_export(self):
        if self._df is None:
            QMessageBox.information(self, "Nothing to export",
                                    "Load a master_results.csv first.")
            return
        default_name = self._suggest_export_name()
        default_path = (
            os.path.join(self._output_dir, default_name)
            if self._output_dir else default_name
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Export figure", default_path,
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)"
        )
        if not path:
            return
        try:
            self._figure.savefig(path, dpi=300, bbox_inches="tight",
                                 facecolor=self._cur_style["figure.facecolor"])
            # Stay consistent with the current theme on save
            self._status_label.setText(f"Saved: {os.path.basename(path)}")
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))

    def _suggest_export_name(self) -> str:
        ptype = self._plot_type.currentText().lower().replace(" ", "_").replace("-", "_")
        fam = self._family_combo.currentText()
        grp = self._group_combo.currentText()
        parts = [ptype]
        if fam and fam != "(all)":
            parts.append(re.sub(r"[^A-Za-z0-9]+", "_", fam))
        if grp and grp != "(all)":
            parts.append(f"g{grp}")
        return "_".join(parts) + ".png"

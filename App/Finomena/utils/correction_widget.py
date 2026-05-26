"""
Correction Widget
=================
Configures the multiple-testing correction strategy for the BAM pipeline and
writes a ``correction.json`` sidecar that the R script reads at run time.

Three top-level choices:

    Strategy:     Flat | Tree
    Error rate:   FDR  | FWER
    Contrast set: All-pairs | Reference-only

If Tree is chosen, the user designs an ordered list of levels (drag to
reorder). Each level specifies:

    Test_Family   – which contrast family it draws from
    Split_By      – optional filter on the Split_By column
    Pair filter   – optional regex on Tested_Level
    Linkage key   – how identity carries from this level's parent

The correction.json schema mirrors the R-side dispatcher's expectation:

    {
      "strategy":     "flat" | "tree",
      "error_rate":   "FDR" | "FWER",
      "contrast_set": "all_pairs" | "ref_only",
      "threshold":    0.05,
      "scope":        "per_family" | "pooled",   (flat only)
      "tree": [
        {
          "name":        "L1 phase gate",
          "Test_Family": "Genotype_Effect",
          "split_by":    "DMSO",
          "pair_filter": null,
          "linkage_key": null
        },
        ...
      ]
    }
"""

import json
import os

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFrame,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QPushButton, QRadioButton, QScrollArea,
    QSizePolicy, QSplitter, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch


# Default Test_Family options. The widget can extend this list dynamically
# from the experimental-design variables — every <Var>_Effect becomes an
# eligible Test_Family, plus the special "Full_Interaction".
_DEFAULT_TEST_FAMILIES = ["Full_Interaction"]


class CorrectionWidget(QWidget):
    """Top-level Correction-tab widget. Owns the strategy state and writes
    correction.json into the output directory at run time."""

    # Emitted when the user changes anything that affects the sidecar.
    spec_changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._variable_names: list[str] = ["Genotype", "Drug"]
        self._variable_refs:  dict      = {}
        self._output_dir:     str       = ""
        self._build_ui()
        self._update_visibility()

    # ── Public API ────────────────────────────────────────────────────────────

    def set_variable_names(self, names: list):
        self._variable_names = list(names)
        self._refresh_test_family_choices()

    def set_references(self, variable_refs: dict, ref_condition: str):
        self._variable_refs = dict(variable_refs)
        # Reset-to-default behaves more sensibly once we know the refs.

    def set_output_dir(self, path: str):
        self._output_dir = path

    def get_spec(self) -> dict:
        """Build the correction.json spec from current control state."""
        if self._flat_radio.isChecked():
            return {
                "strategy":     "flat",
                "error_rate":   "FDR" if self._fdr_radio.isChecked() else "FWER",
                "contrast_set": "all_pairs" if self._allpairs_radio.isChecked()
                                            else "ref_only",
                "threshold":    float(self._threshold_spin.value()),
                "scope":        ("per_family" if self._per_family_radio.isChecked()
                                              else "pooled"),
            }
        # Tree
        return {
            "strategy":     "tree",
            "error_rate":   "FDR" if self._fdr_radio.isChecked() else "FWER",
            "contrast_set": "all_pairs" if self._allpairs_radio.isChecked()
                                        else "ref_only",
            "threshold":    float(self._threshold_spin.value()),
            "tree":         self._collect_tree_levels(),
        }

    def write_sidecar(self) -> str:
        """Serialize the current spec to <output_dir>/correction.json.
        Returns the absolute path written, or empty string on failure."""
        if not self._output_dir:
            return ""
        try:
            spec = self.get_spec()
            path = os.path.join(self._output_dir, "correction.json")
            with open(path, "w", encoding="utf-8") as f:
                json.dump(spec, f, indent=2)
            return path
        except Exception:
            return ""

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        outer = QVBoxLayout(self)

        # ── Top: strategy switches ──────────────────────────────────────────
        strat_group = QGroupBox("Correction Strategy")
        sg_layout = QHBoxLayout(strat_group)

        # Strategy: Flat / Tree
        strat_col = QVBoxLayout()
        strat_col.addWidget(QLabel("<b>Structure</b>"))
        self._flat_radio = QRadioButton("Flat (no hierarchy)")
        self._tree_radio = QRadioButton("Tree (gatekeeping hierarchy)")
        self._flat_radio.setChecked(True)
        self._flat_radio.toggled.connect(self._on_any_change)
        self._tree_radio.toggled.connect(self._on_any_change)
        strat_col.addWidget(self._flat_radio)
        strat_col.addWidget(self._tree_radio)
        strat_col.addStretch()
        sg_layout.addLayout(strat_col)

        # Error rate: FDR / FWER
        err_col = QVBoxLayout()
        err_col.addWidget(QLabel("<b>Error rate</b>"))
        self._fdr_radio  = QRadioButton("FDR (Benjamini-Hochberg / TreeBH)")
        self._fwer_radio = QRadioButton("FWER (Holm / graphicalMCP / Holm-gatekeeping)")
        self._fdr_radio.setChecked(True)
        self._fdr_radio.toggled.connect(self._on_any_change)
        self._fwer_radio.toggled.connect(self._on_any_change)
        err_col.addWidget(self._fdr_radio)
        err_col.addWidget(self._fwer_radio)
        err_col.addStretch()
        sg_layout.addLayout(err_col)

        # Contrast set: all-pairs / ref-only
        cs_col = QVBoxLayout()
        cs_col.addWidget(QLabel("<b>Contrast set</b>"))
        self._allpairs_radio = QRadioButton("All comparisons (pairwise)")
        self._refonly_radio  = QRadioButton("Reference-only (treatment vs ref)")
        self._refonly_radio.setChecked(True)
        self._allpairs_radio.toggled.connect(self._on_any_change)
        self._refonly_radio.toggled.connect(self._on_any_change)
        cs_col.addWidget(self._allpairs_radio)
        cs_col.addWidget(self._refonly_radio)
        cs_col.addStretch()
        sg_layout.addLayout(cs_col)

        # Threshold + flat-scope
        thr_col = QVBoxLayout()
        thr_col.addWidget(QLabel("<b>Threshold</b>"))
        thr_row = QHBoxLayout()
        thr_row.addWidget(QLabel("q (FDR) / α (FWER):"))
        self._threshold_spin = QDoubleSpinBox()
        self._threshold_spin.setRange(0.001, 0.5)
        self._threshold_spin.setSingleStep(0.01)
        self._threshold_spin.setDecimals(3)
        self._threshold_spin.setValue(0.05)
        self._threshold_spin.setMinimumWidth(120)
        self._threshold_spin.valueChanged.connect(self._on_any_change)
        thr_row.addWidget(self._threshold_spin)
        thr_row.addStretch()
        thr_col.addLayout(thr_row)

        self._scope_label = QLabel("<b>Flat scope</b>")
        thr_col.addWidget(self._scope_label)
        self._per_family_radio = QRadioButton("Per Test_Family")
        self._pooled_radio     = QRadioButton("Pooled across all")
        self._per_family_radio.setChecked(True)
        self._per_family_radio.toggled.connect(self._on_any_change)
        self._pooled_radio.toggled.connect(self._on_any_change)
        thr_col.addWidget(self._per_family_radio)
        thr_col.addWidget(self._pooled_radio)
        thr_col.addStretch()
        sg_layout.addLayout(thr_col)

        # Group them so only one strategy/error/contrast is selected
        self._strat_grp  = QButtonGroup(self)
        self._strat_grp.addButton(self._flat_radio)
        self._strat_grp.addButton(self._tree_radio)
        self._err_grp    = QButtonGroup(self)
        self._err_grp.addButton(self._fdr_radio)
        self._err_grp.addButton(self._fwer_radio)
        self._cs_grp     = QButtonGroup(self)
        self._cs_grp.addButton(self._allpairs_radio)
        self._cs_grp.addButton(self._refonly_radio)
        self._scope_grp  = QButtonGroup(self)
        self._scope_grp.addButton(self._per_family_radio)
        self._scope_grp.addButton(self._pooled_radio)

        outer.addWidget(strat_group)

        # ── Below: split between tree designer (left) and preview (right) ──
        self._splitter = QSplitter(Qt.Horizontal)

        # Designer
        designer = QWidget()
        d_layout = QVBoxLayout(designer)
        d_layout.addWidget(QLabel("<b>Tree levels</b> (drag to reorder)"))

        self._level_list = QListWidget()
        self._level_list.setDragDropMode(QListWidget.InternalMove)
        self._level_list.setDefaultDropAction(Qt.MoveAction)
        self._level_list.setSelectionMode(QListWidget.SingleSelection)
        self._level_list.model().rowsMoved.connect(self._on_any_change)
        d_layout.addWidget(self._level_list, stretch=1)

        btn_row = QHBoxLayout()
        self._add_btn   = QPushButton("Add level")
        self._add_btn.clicked.connect(self._on_add_level)
        self._reset_btn = QPushButton("Reset to default rescue tree")
        self._reset_btn.clicked.connect(self._on_reset)
        btn_row.addWidget(self._add_btn)
        btn_row.addWidget(self._reset_btn)
        btn_row.addStretch()
        d_layout.addLayout(btn_row)

        self._splitter.addWidget(designer)

        # Preview / visualization — schematic tree diagram showing the
        # structure with variable-name placeholders rather than expanded
        # individual contrasts (a literal rendering of 8 phases × 82 drugs
        # would be unreadable).
        preview = QWidget()
        p_layout = QVBoxLayout(preview)
        p_layout.addWidget(QLabel("<b>Tree preview</b>"))
        self._preview_figure = Figure(figsize=(6, 5))
        self._preview_figure.patch.set_facecolor("#202124")
        self._preview_canvas = FigureCanvas(self._preview_figure)
        self._preview_canvas.setSizePolicy(QSizePolicy.Expanding,
                                           QSizePolicy.Expanding)
        p_layout.addWidget(self._preview_canvas, stretch=1)
        self._splitter.addWidget(preview)

        # Default split: ~55% designer / ~45% preview (5% right-shift from
        # 50/50). User can drag the divider freely after that.
        self._splitter.setSizes([550, 450])
        outer.addWidget(self._splitter, stretch=1)

        # Status line at bottom
        self._status_label = QLabel("")
        self._status_label.setStyleSheet("color: #aaaaaa; font-style: italic; padding: 4px;")
        outer.addWidget(self._status_label)

        # Seed with default rescue tree + render preview
        self._populate_default_tree()
        self._update_status()
        self._render_preview()

    # ── State helpers ─────────────────────────────────────────────────────────

    def _on_any_change(self, *_):
        self._update_visibility()
        self.spec_changed.emit()
        self._update_status()
        self._render_preview()

    def _update_visibility(self):
        # Tree designer only visible when Tree strategy is active
        tree_mode = self._tree_radio.isChecked()
        self._splitter.setVisible(tree_mode)
        # Flat scope radios only visible when Flat strategy is active
        flat_mode = self._flat_radio.isChecked()
        self._scope_label.setVisible(flat_mode)
        self._per_family_radio.setVisible(flat_mode)
        self._pooled_radio.setVisible(flat_mode)

    def _update_status(self):
        spec = self.get_spec()
        strat   = spec["strategy"]
        err     = spec["error_rate"]
        cs      = spec["contrast_set"]
        thr     = spec["threshold"]
        if strat == "flat":
            method = "BH" if err == "FDR" else "Holm"
        else:
            if err == "FDR":
                method = "TreeBH"
            else:
                method = "graphicalMCP" if cs == "ref_only" else "Holm-gatekeeping"
        self._status_label.setText(
            f"Selected: {strat} · {err} · {cs} · {method} at {thr:g}"
        )

    # ── Tree level rows ──────────────────────────────────────────────────────

    def _refresh_test_family_choices(self):
        """Per-variable Test_Family options come from the design variables."""
        choices = [f"{v}_Effect" for v in self._variable_names] + _DEFAULT_TEST_FAMILIES
        for i in range(self._level_list.count()):
            item = self._level_list.item(i)
            row  = self._level_list.itemWidget(item)
            if isinstance(row, _LevelRow):
                row.set_test_family_choices(choices)

    def _on_add_level(self):
        self._append_level({
            "name":        f"Level {self._level_list.count() + 1}",
            "Test_Family": "Full_Interaction",
            "split_by":    "",
            "pair_filter": "",
            "linkage_key": "Group",
        })
        self._on_any_change()

    def _on_reset(self):
        self._populate_default_tree()
        self._on_any_change()

    def _populate_default_tree(self):
        self._level_list.clear()
        if len(self._variable_names) >= 2:
            gate_var   = self._variable_names[0]
            screen_var = self._variable_names[1]
            ref_screen = self._variable_refs.get(screen_var, "")
            defaults = [
                {
                    "name":        f"L1 {gate_var} gate",
                    "Test_Family": f"{gate_var}_Effect",
                    "split_by":    ref_screen,
                    "pair_filter": "",
                    "linkage_key": "",
                },
                {
                    "name":        "L2 rescue claim",
                    "Test_Family": "Full_Interaction",
                    "split_by":    "",
                    "pair_filter": "",
                    "linkage_key": "Group",
                },
                {
                    "name":        f"L3 {screen_var} specificity",
                    "Test_Family": f"{screen_var}_Effect",
                    "split_by":    "",
                    "pair_filter": "",
                    "linkage_key": "Tested_Level",
                },
            ]
        else:
            defaults = []
        for spec in defaults:
            self._append_level(spec)

    def _append_level(self, spec: dict):
        row = _LevelRow(self._variable_names, spec, parent=self._level_list)
        row.changed.connect(self._on_any_change)
        row.remove_requested.connect(self._on_remove_level)
        item = QListWidgetItem(self._level_list)
        item.setSizeHint(row.sizeHint())
        self._level_list.addItem(item)
        self._level_list.setItemWidget(item, row)

    def _on_remove_level(self, row_widget):
        for i in range(self._level_list.count()):
            item = self._level_list.item(i)
            if self._level_list.itemWidget(item) is row_widget:
                self._level_list.takeItem(i)
                break
        self._on_any_change()

    # ── Tree preview rendering ────────────────────────────────────────────

    def _render_preview(self):
        """Draw a schematic of the current spec on the preview canvas."""
        fig = self._preview_figure
        fig.clear()
        ax  = fig.add_subplot(111)
        ax.set_facecolor("#202124")
        ax.set_xlim(0, 10)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        spec = self.get_spec()
        strat = spec["strategy"]
        err   = spec["error_rate"]
        cs    = spec["contrast_set"]
        thr   = spec["threshold"]

        if strat == "flat":
            method = "BH" if err == "FDR" else "Holm"
            ax.text(5, 5.5,
                    f"Flat correction\n{method}  ·  threshold {thr:g}",
                    ha="center", va="center",
                    color="#e8e8e8", fontsize=12,
                    bbox=dict(boxstyle="round,pad=0.6",
                              facecolor="#37474f", edgecolor="#cfd8dc"))
            ax.text(5, 3,
                    f"Scope: {spec.get('scope', 'per_family')}\n"
                    f"Contrasts: {cs}",
                    ha="center", va="center",
                    color="#cfd8dc", fontsize=10)
            ax.set_ylim(0, 10)
            self._preview_canvas.draw_idle()
            return

        # Tree: stack boxes top-to-bottom with linkage arrows
        levels = spec.get("tree", [])
        n = len(levels)
        if n == 0:
            ax.text(5, 5, "No levels defined yet.\nClick 'Add level' to begin.",
                    ha="center", va="center", color="#aaaaaa",
                    fontsize=11, fontstyle="italic")
            ax.set_ylim(0, 10)
            self._preview_canvas.draw_idle()
            return

        # Lay out one box per level, top to bottom
        ax.set_ylim(0, n + 1.5)
        # Algorithm label
        if err == "FDR":
            algo = "TreeBH"
        else:
            algo = "graphicalMCP (correlation-aware)" if cs == "ref_only" else "Holm-gatekeeping"
        ax.text(5, n + 1.0,
                f"Tree · {err} · {cs} · {algo} at {thr:g}",
                ha="center", va="bottom",
                color="#80cbc4", fontsize=11, fontweight="bold")

        for i, lvl in enumerate(levels):
            y = n - i  # top → bottom
            label = lvl.get("name", f"Level {i+1}") or f"Level {i+1}"
            fam   = lvl.get("Test_Family", "—") or "—"
            sb    = lvl.get("split_by") or "(all)"
            pf    = lvl.get("pair_filter") or "(none)"
            lk    = lvl.get("linkage_key")
            if i == 0:
                lk_str = "(root)"
            else:
                if isinstance(lk, dict):
                    lk_str = lk.get("parent_var", "Group")
                else:
                    lk_str = lk or "Group"

            text = (
                f"{label}\n"
                f"family: {fam}    Split_By: {sb}    pair: {pf}\n"
                f"linkage from parent: {lk_str}"
            )
            ax.text(5, y, text,
                    ha="center", va="center",
                    color="#e8e8e8", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.4",
                              facecolor="#37474f", edgecolor="#cfd8dc"))
            # Arrow to next level (if any)
            if i < n - 1:
                ax.annotate("",
                            xy=(5, y - 0.55), xytext=(5, y - 0.35),
                            arrowprops=dict(arrowstyle="->", color="#80cbc4",
                                            lw=1.2))

        try:
            fig.tight_layout()
        except Exception:
            pass
        self._preview_canvas.draw_idle()

    def _collect_tree_levels(self) -> list:
        levels = []
        for i in range(self._level_list.count()):
            item = self._level_list.item(i)
            row  = self._level_list.itemWidget(item)
            if not isinstance(row, _LevelRow):
                continue
            spec = row.get_spec()
            # Normalize empty strings to None for cleaner JSON
            for k in ("split_by", "pair_filter"):
                if spec.get(k) == "":
                    spec[k] = None
            # First level has no linkage key
            if i == 0:
                spec["linkage_key"] = None
            else:
                lk = spec.get("linkage_key") or "Group"
                # JSON shape mirrors what tree_metadata.R expects
                spec["linkage_key"] = {"parent_var": lk, "child_var": lk}
            levels.append(spec)
        return levels


class _LevelRow(QFrame):
    """One row in the tree designer — represents one level's spec."""

    changed         = Signal()
    remove_requested = Signal(object)

    def __init__(self, variable_names: list, spec: dict, parent=None):
        super().__init__(parent)
        self.setFrameShape(QFrame.StyledPanel)
        self._variable_names = list(variable_names)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(6, 4, 6, 4)
        layout.setSpacing(8)

        # Drag handle marker (purely visual)
        layout.addWidget(QLabel("≡"))

        # Editable level name
        self._name_edit = QLineEdit(spec.get("name", "Level"))
        self._name_edit.setMaximumWidth(160)
        self._name_edit.textChanged.connect(self.changed)
        layout.addWidget(self._name_edit)

        layout.addWidget(QLabel("Family:"))
        self._family_combo = QComboBox()
        self.set_test_family_choices(
            [f"{v}_Effect" for v in variable_names] + _DEFAULT_TEST_FAMILIES
        )
        idx = self._family_combo.findText(spec.get("Test_Family", ""))
        if idx >= 0:
            self._family_combo.setCurrentIndex(idx)
        self._family_combo.currentIndexChanged.connect(self.changed)
        layout.addWidget(self._family_combo)

        layout.addWidget(QLabel("Split_By:"))
        self._split_edit = QLineEdit(spec.get("split_by", "") or "")
        self._split_edit.setMaximumWidth(120)
        self._split_edit.setPlaceholderText("(all)")
        self._split_edit.textChanged.connect(self.changed)
        layout.addWidget(self._split_edit)

        layout.addWidget(QLabel("Pair filter:"))
        self._pair_edit = QLineEdit(spec.get("pair_filter", "") or "")
        self._pair_edit.setMaximumWidth(140)
        self._pair_edit.setPlaceholderText("(none, regex on Tested_Level)")
        self._pair_edit.textChanged.connect(self.changed)
        layout.addWidget(self._pair_edit)

        layout.addWidget(QLabel("Linkage:"))
        self._linkage_combo = QComboBox()
        link_options = ["Group", "Tested_Level"] + list(variable_names)
        for opt in link_options:
            self._linkage_combo.addItem(opt)
        idx = self._linkage_combo.findText(spec.get("linkage_key", "Group"))
        if idx >= 0:
            self._linkage_combo.setCurrentIndex(idx)
        self._linkage_combo.currentIndexChanged.connect(self.changed)
        layout.addWidget(self._linkage_combo)

        layout.addStretch()

        self._remove_btn = QPushButton("✕")
        self._remove_btn.setMaximumWidth(28)
        self._remove_btn.clicked.connect(lambda: self.remove_requested.emit(self))
        layout.addWidget(self._remove_btn)

    def set_test_family_choices(self, choices: list):
        current = self._family_combo.currentText() if hasattr(self, "_family_combo") else ""
        self._family_combo.blockSignals(True)
        self._family_combo.clear()
        # Deduplicate while preserving order
        seen = set()
        for c in choices:
            if c not in seen:
                self._family_combo.addItem(c)
                seen.add(c)
        idx = self._family_combo.findText(current) if current else -1
        if idx >= 0:
            self._family_combo.setCurrentIndex(idx)
        self._family_combo.blockSignals(False)

    def get_spec(self) -> dict:
        return {
            "name":        self._name_edit.text().strip(),
            "Test_Family": self._family_combo.currentText(),
            "split_by":    self._split_edit.text().strip(),
            "pair_filter": self._pair_edit.text().strip(),
            "linkage_key": self._linkage_combo.currentText(),
        }

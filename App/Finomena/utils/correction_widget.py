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
    Pair filter   – optional regex on Passed_Contrasts
    Linkage key   – how identity carries from this level's parent

The correction.json schema mirrors the R-side dispatcher's expectation:

    {
      "strategy":     "flat" | "tree",
      "error_rate":   "FDR" | "FWER",
      "contrast_set": "all_pairs" | "ref_only",
      "threshold":    0.05,
      "scope":        "per_family",               (flat only; pooled removed)
      "tree": [
        {
          "name":        "L1 phase gate",
          "Test_Family": "Genotype_Effect",
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
        self._conditions:     dict      = {}    # {full_interaction_name: hex_color}
        self._output_dir:     str       = ""
        self._build_ui()
        self._update_visibility()

    # ── Public API ────────────────────────────────────────────────────────────

    def set_variable_names(self, names: list):
        self._variable_names = list(names)
        # Rebuild the default tree from the current variable list:
        #   L1 = <var0>_Effect (gate)
        #   L2 = <var1>_Effect (effect)
        #   L3 = Full_Interaction (specificity)
        # This is the canonical rescue-paradigm tree, simpler→complex order.
        # Rebuilding here guarantees the default reflects the current design
        # variables.
        self._populate_default_tree()
        self._refresh_test_family_choices()
        self._render_preview()

    def set_references(self, variable_refs: dict, ref_condition: str):
        self._variable_refs = dict(variable_refs)
        # Reset-to-default behaves more sensibly once we know the refs.

    def set_conditions(self, conditions: dict):
        """Wired to ExperimentalConditionsWidget.conditions_updated. The dict
        is {full_interaction_name (e.g. "WT+DMSO"): hex_color}. We don't
        care about colors; we use the keys to build the condition_lookup
        that gets shipped to R via correction.json (so the design-variable
        parser can decode multi-variable contrast strings by dict access)."""
        self._conditions = dict(conditions)

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
                # scope is fixed to per_family — pooled across all was removed
                # because each Test_Family represents a distinct scientific
                # question (Genotype / Drug / Full_Interaction).
                "scope":        "per_family",
            }
        # Tree
        return {
            "strategy":          "tree",
            "error_rate":        "FDR" if self._fdr_radio.isChecked() else "FWER",
            "contrast_set":      "all_pairs" if self._allpairs_radio.isChecked()
                                             else "ref_only",
            "threshold":         float(self._threshold_spin.value()),
            "tree":              self._collect_tree_levels(),
            # Per-variable decoding for every condition name. Python already
            # knows these — passing them to R lets tree_metadata.R decode
            # Full_Interaction contrast strings via dict access instead of
            # parsing "+"-joined condition names.
            "condition_lookup":  self._build_condition_lookup(),
        }

    def _build_condition_lookup(self) -> dict:
        """{condition_name: {var_name: value, ...}} for every loaded condition."""
        out = {}
        for cond in self._conditions.keys():
            parts = cond.split("+")
            if len(parts) != len(self._variable_names):
                continue
            out[cond] = {v: p.strip() for v, p in zip(self._variable_names, parts)}
        return out

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

        # "Pooled across all" was removed; flat correction is always
        # per Test_Family because each family represents a distinct
        # scientific question. See ABLATION_METHODOLOGY.md / earlier
        # discussion of scope choice.
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
        # When a row is dragged: refresh whichever row is now at position 0
        # (root) so its linkage combo is greyed/None, and any previously-root
        # row that moved DOWN gets its combo re-enabled.
        self._level_list.model().rowsMoved.connect(self._refresh_root_level_state)
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

    def _refresh_root_level_state(self):
        """Whichever row sits at position 0 is the L1 root; its linkage_key
        is structurally ignored by the R-side correction algorithms (BH/Holm
        at L1 corrects across all rows with one shared alpha). Reflect that
        in the UI by forcing the linkage combo to "None" and disabling it
        for the row at index 0. Call this whenever the row order changes."""
        for i in range(self._level_list.count()):
            item = self._level_list.item(i)
            row  = self._level_list.itemWidget(item)
            if isinstance(row, _LevelRow):
                row.set_root_level(i == 0)

    def _on_add_level(self):
        self._append_level({
            "name":        f"Level {self._level_list.count() + 1}",
            "Test_Family": "Full_Interaction",
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
            # split_by was removed — every row of a Test_Family is covered by
            # exactly one tree level, so every computed p-value gets corrected.
            # Default tree (rescue paradigm), simpler → more complex hypothesis:
            #   L1 = gate     (Genotype_Effect — all contexts)
            #     "Does the disease phenotype exist?"
            #   L2 = effect   (Drug_Effect — both genotype contexts)
            #     "Does the drug do anything?"
            #   L3 = specificity (Full_Interaction)
            #     "Is the drug's effect specifically rescue (differs by genotype)?"
            defaults = [
                {
                    "name":        f"L1 {gate_var} gate",
                    "Test_Family": f"{gate_var} Effect",
                    "linkage_key": None,   # root — auto-forced to None at i==0
                },
                {
                    "name":        f"L2 {screen_var} effect",
                    "Test_Family": f"{screen_var} Effect",
                    "linkage_key": f"{gate_var}",
                },
                {
                    "name":        f"L3 Interaction specificity",
                    "Test_Family": "Full_Interaction",
                    "linkage_key": f"{screen_var}",
                },
            ]
        else:
            defaults = []
        for spec in defaults:
            self._append_level(spec)

    def _append_level(self, spec: dict):
        row = _LevelRow(
            self._variable_names, spec,
            parent=self._level_list,
        )
        row.changed.connect(self._on_any_change)
        row.remove_requested.connect(self._on_remove_level)
        item = QListWidgetItem(self._level_list)
        item.setSizeHint(row.sizeHint())
        self._level_list.addItem(item)
        self._level_list.setItemWidget(item, row)
        self._refresh_root_level_state()

    def _on_remove_level(self, row_widget):
        for i in range(self._level_list.count()):
            item = self._level_list.item(i)
            if self._level_list.itemWidget(item) is row_widget:
                self._level_list.takeItem(i)
                break
        self._refresh_root_level_state()
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
                    f"Scope: per Test_Family\n"
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

        # Lay out: header at top, then a "Group root" indicator, then one
        # box per level descending. Total vertical span = n + 2.5 units.
        ax.set_ylim(0, n + 2.5)

        # Algorithm label, with a one-line note about alpha-redistribution
        # behaviour so the user knows what the method actually does.
        if err == "FDR":
            algo = "TreeBH"
            redist_note = "BH per parent-rejected family · no alpha redistribution"
        else:
            if cs == "ref_only":
                algo = "graphicalMCP (correlation-aware)"
                redist_note = "alpha redistributes via graph edges (children + siblings)"
            else:
                algo = "Holm-gatekeeping"
                redist_note = "Holm per parent-rejected family · no alpha redistribution"
        ax.text(5, n + 2.2,
                f"Tree · {err} · {cs} · {algo} at {thr:g}",
                ha="center", va="bottom",
                color="#80cbc4", fontsize=11, fontweight="bold")
        ax.text(5, n + 1.95,
                redist_note,
                ha="center", va="top",
                color="#9eaeb6", fontsize=8, fontstyle="italic")

        # Group root indicator — visualises that L1 has one node per phase,
        # all siblings of an implicit Group root. Important for understanding
        # alpha flow at L1 (especially under graphicalMCP).
        ax.text(5, n + 1.0,
                "Group root\nL1 phase gate siblings share this level's alpha pool",
                ha="center", va="center",
                color="#cfd8dc", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.4",
                          facecolor="#26323a", edgecolor="#546e7a",
                          linestyle="dashed"))
        # Arrow from root down to L1
        ax.annotate("",
                    xy=(5, n + 0.45), xytext=(5, n + 0.75),
                    arrowprops=dict(arrowstyle="->", color="#80cbc4", lw=1.2))

        for i, lvl in enumerate(levels):
            y = n - i  # top → bottom
            label = lvl.get("name", f"Level {i+1}") or f"Level {i+1}"
            fam   = lvl.get("Test_Family", "—") or "—"
            lk    = lvl.get("linkage_key")
            if i == 0:
                lk_str = "(from Group root, per-phase siblings)"
            else:
                if isinstance(lk, dict):
                    lk_str = lk.get("parent_var", "Group")
                else:
                    lk_str = lk or "Group"

            text = (
                f"{label}\n"
                f"family: {fam}\n"
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
            # First level has no linkage key
            if i == 0:
                spec["linkage_key"] = None
            else:
                lk = spec.get("linkage_key") or "Group"
                # Safety net: the "None" option is meant to live only on the
                # L1 row. If it somehow leaks into a non-root spec (e.g. timing
                # during reorder), fall back to "Group".
                if lk == "None":
                    lk = "Group"
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

        # split_by + pair_filter were both removed — narrowing tree levels
        # would leave some computed p-values uncorrected (off_tree). p-value
        # correction must be all-or-nothing for the family of computed tests.
        # Pair narrowing belongs in Contrast Selection / plate layout.

        layout.addWidget(QLabel("Linkage:"))
        self._linkage_combo = QComboBox()
        link_options = ["Group", "Passed_Contrasts"] + list(variable_names)
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
            "linkage_key": self._linkage_combo.currentText(),
        }

    def set_root_level(self, is_root: bool):
        """L1 (position 0) has no parent — corrections.R applies BH/Holm
        across ALL its rows with one shared alpha, ignoring whatever
        linkage_key the JSON carries. Reflect that in the UI by forcing
        the combo to "None" and disabling it for the root row. If this
        row is later moved out of position 0 (or another row takes its
        place), the combo re-enables with its prior options."""
        self._linkage_combo.blockSignals(True)
        none_idx = self._linkage_combo.findText("None")
        if is_root:
            if none_idx < 0:
                self._linkage_combo.insertItem(0, "None")
            self._linkage_combo.setCurrentIndex(self._linkage_combo.findText("None"))
            self._linkage_combo.setEnabled(False)
        else:
            if none_idx >= 0:
                self._linkage_combo.removeItem(none_idx)
            self._linkage_combo.setEnabled(True)
        self._linkage_combo.blockSignals(False)

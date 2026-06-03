"""
Contrast Selection Widget
=========================
Lets the user choose which pairwise Full_Interaction contrasts to include
in statistical testing.

Excluded pairs are filtered BEFORE both the within-family (Sidak) and
cross-phase (FDR/Holm) corrections, so kept comparisons get full power.
"""

from itertools import combinations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox, QDoubleSpinBox, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QScrollArea, QSpinBox, QVBoxLayout, QWidget,
)


class ContrastSelectionWidget(QWidget):
    """Tab: pick which pairwise condition contrasts to keep for BAM analysis."""

    selection_changed = Signal()  # fired on any checkbox state change

    def __init__(self, parent=None):
        super().__init__(parent)
        self._conditions: dict = {}         # {cond_name: hex_color}
        self._roles: dict = {}              # {cond_name: role_string}
        self._pair_checkboxes: dict = {}    # {(condA, condB): QCheckBox}
        self._pending_excluded: set = set() # applied when pair list is (re)built
        # Tracks whether the user has manually toggled a checkbox. Until True,
        # arriving role info is allowed to re-apply role-based defaults.
        self._user_has_interacted: bool = False
        self._build_ui()
        self._update_summary()

    # ── Public API ─────────────────────────────────────────────────────────

    def set_roles(self, roles: dict):
        """
        Receives {condition_name: role_string} from the metadata assignment widget.
        If the user hasn't manually toggled anything yet, re-applies the role-based
        defaults to all existing pairs (Experimental vs Reference/Positive Control
        are checked; everything else is unchecked).
        """
        self._roles = dict(roles)
        if self._user_has_interacted or not self._pair_checkboxes:
            return
        for pair, cb in self._pair_checkboxes.items():
            cb.blockSignals(True)
            cb.setChecked(self._is_default_pair(pair))
            cb.blockSignals(False)
        self._update_summary()
        self.selection_changed.emit()

    def on_bam_completed(self):
        """
        Called after BAM analysis completes successfully. Shows a persistent
        warning so the user knows that any contrast change here will
        invalidate the BAM results and require a re-run.
        """
        self._bam_run_warning_label.show()

    # ── UI ─────────────────────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)

        info = QLabel(
            "<b>Pairwise comparisons to include in statistical testing</b><br>"
            "Uncheck any comparison you don't care about. Excluded pairs are "
            "dropped <i>before</i> both the within-family (Sidak) and the "
            "cross-phase (FDR/Holm) corrections, so your kept comparisons "
            "keep full statistical power."
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        # Shown after a successful BAM run so the user knows that any contrast
        # change here will invalidate those results and require a re-run.
        self._bam_run_warning_label = QLabel(
            "⚠ BAM has been run. Changing contrasts here will require re-running the BAM analysis."
        )
        self._bam_run_warning_label.setStyleSheet(
            "color: #c0392b; font-weight: bold; padding: 4px;"
        )
        self._bam_run_warning_label.setWordWrap(True)
        self._bam_run_warning_label.hide()
        layout.addWidget(self._bam_run_warning_label)

        btn_row = QHBoxLayout()
        self._select_all_btn = QPushButton("Select All")
        self._select_all_btn.clicked.connect(self._select_all)
        self._deselect_all_btn = QPushButton("Deselect All")
        self._deselect_all_btn.clicked.connect(self._deselect_all)
        btn_row.addWidget(self._select_all_btn)
        btn_row.addWidget(self._deselect_all_btn)
        btn_row.addSpacing(12)
        btn_row.addWidget(QLabel("Filter:"))
        self._filter_edit = QLineEdit()
        self._filter_edit.setPlaceholderText("substring match on either side of 'vs'")
        self._filter_edit.textChanged.connect(self._apply_filter)
        btn_row.addWidget(self._filter_edit, 1)
        layout.addLayout(btn_row)

        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll_content = QWidget()
        self._scroll_layout = QVBoxLayout(self._scroll_content)
        self._scroll_layout.setAlignment(Qt.AlignTop)
        self._scroll_layout.setSpacing(2)
        self._scroll.setWidget(self._scroll_content)
        layout.addWidget(self._scroll, 1)

        self._summary_label = QLabel("")
        self._summary_label.setWordWrap(True)
        layout.addWidget(self._summary_label)

        # ── Posterior equivalence section ─────────────────────────────────
        pe_group = QGroupBox("Posterior Equivalence Testing (Bayesian)")
        pe_layout = QVBoxLayout(pe_group)

        self._pe_enabled_cb = QCheckBox(
            "Run posterior equivalence test for the kept pairs (during BAM run)"
        )
        self._pe_enabled_cb.stateChanged.connect(self._on_pe_toggle)
        pe_layout.addWidget(self._pe_enabled_cb)

        pe_help = QLabel(
            "<i>For every kept pair, draws from the BAM posterior to compute "
            "<b>M</b> = the maximum absolute log<sub>2</sub> fold change between "
            "the two condition trajectories across each phase's time course. "
            "Reports the posterior distribution of <b>M</b> and "
            "<b>Pr(M &lt; δ)</b> — a direct Bayesian credibility statement of "
            "equivalence within tolerance δ.</i>"
        )
        pe_help.setWordWrap(True)
        pe_help.setStyleSheet("color: #888; padding-left: 22px;")
        pe_layout.addWidget(pe_help)

        pe_settings = QHBoxLayout()
        pe_settings.addSpacing(22)
        pe_settings.addWidget(QLabel("δ (log<sub>2</sub> fold change):"))
        self._pe_delta_spin = QDoubleSpinBox()
        self._pe_delta_spin.setRange(0.01, 10.0)
        self._pe_delta_spin.setSingleStep(0.1)
        self._pe_delta_spin.setDecimals(2)
        self._pe_delta_spin.setValue(1.0)
        self._pe_delta_spin.setToolTip(
            "Equivalence margin on log₂ scale. Default 1.0 = trajectories never differ "
            "by more than a factor of 2 at any time point."
        )
        pe_settings.addWidget(self._pe_delta_spin)
        pe_settings.addSpacing(20)
        pe_settings.addWidget(QLabel("Posterior draws:"))
        self._pe_draws_spin = QSpinBox()
        self._pe_draws_spin.setRange(100, 100000)
        self._pe_draws_spin.setSingleStep(1000)
        self._pe_draws_spin.setValue(10000)
        self._pe_draws_spin.setToolTip(
            "Number of posterior samples used to estimate the M distribution and Pr(M < δ). "
            "Default 10,000 gives ~0.5% Monte Carlo precision on Pr. Lower values run faster "
            "but produce a noisier estimate."
        )
        pe_settings.addWidget(self._pe_draws_spin)
        pe_settings.addStretch()
        pe_layout.addLayout(pe_settings)

        pe_draws_note = QLabel(
            "<i>10,000 draws is recommended for final results. Reduce for "
            "quicker runs at lower precision.</i>"
        )
        pe_draws_note.setWordWrap(True)
        pe_draws_note.setStyleSheet("color: #888; padding-left: 22px;")
        pe_layout.addWidget(pe_draws_note)

        self._pe_settings_widgets = [self._pe_delta_spin, self._pe_draws_spin]
        self._update_pe_settings_enabled()

        layout.addWidget(pe_group)

    # ── Public API ─────────────────────────────────────────────────────────

    def update_conditions(self, conditions: dict):
        """Slot: receives {cond_name: hex_color} from ExperimentalConditionsWidget."""
        self._conditions = dict(conditions)
        self._rebuild_pair_list()

    def get_kept_pairs(self) -> list:
        """Returns [[lhs, rhs], ...] for every checked pair."""
        return [list(pair) for pair, cb in self._pair_checkboxes.items() if cb.isChecked()]

    def get_excluded_pairs(self) -> list:
        return [list(pair) for pair, cb in self._pair_checkboxes.items() if not cb.isChecked()]

    def get_posterior_settings(self) -> dict:
        """Returns the posterior-equivalence settings (enabled / delta / n_draws)."""
        return {
            "enabled": self._pe_enabled_cb.isChecked(),
            "delta":   float(self._pe_delta_spin.value()),
            "n_draws": int(self._pe_draws_spin.value()),
        }

    def save_config(self) -> dict:
        """Config persists EXCLUDED pairs (defaults to 'all kept') + posterior settings."""
        return {
            "excluded_pairs":       self.get_excluded_pairs(),
            "posterior_equivalence": self.get_posterior_settings(),
        }

    def load_config(self, data: dict):
        """
        Restores excluded-pair set + posterior settings.

        If the saved config has no exclusions (empty list or missing key), the
        role-based defaults stay in place — empty exclusions usually means
        "user never customized contrasts", not "user explicitly wanted every
        pair included". An explicit non-empty exclusion list is treated as a
        deliberate user choice and locks future role updates out.
        """
        excluded = {tuple(self._canon(p)) for p in data.get("excluded_pairs", [])}
        if excluded:
            # Explicit user customization — apply exactly and lock defaults.
            self._user_has_interacted = True
            self._pending_excluded = excluded
            for pair, cb in self._pair_checkboxes.items():
                cb.blockSignals(True)
                cb.setChecked(pair not in excluded)
                cb.blockSignals(False)
            self._update_summary()
        else:
            # No saved exclusions — leave role-based defaults from set_roles in place.
            self._pending_excluded = set()

        pe = data.get("posterior_equivalence") or {}
        self._pe_enabled_cb.setChecked(bool(pe.get("enabled", False)))
        if "delta" in pe:
            try:
                self._pe_delta_spin.setValue(float(pe["delta"]))
            except (TypeError, ValueError):
                pass
        if "n_draws" in pe:
            try:
                self._pe_draws_spin.setValue(int(pe["n_draws"]))
            except (TypeError, ValueError):
                pass
        self._update_pe_settings_enabled()

        self._update_summary()
        self.selection_changed.emit()

    # ── Internal ───────────────────────────────────────────────────────────

    @staticmethod
    def _canon(pair):
        """Canonicalize a pair so (A, B) and (B, A) are treated as equal."""
        a, b = pair
        return (a, b) if a <= b else (b, a)

    # Roles considered "controls" for the default-checked logic.
    _CONTROL_ROLES = ("Reference Control", "Positive Control")

    def _is_default_pair(self, pair) -> bool:
        """
        Default-checked iff one side is Experimental and the other is a control
        (Reference or Positive). When no role info is available, falls back to
        True so behavior matches the old "all checked by default".
        """
        if not self._roles:
            return True
        r1 = self._roles.get(pair[0])
        r2 = self._roles.get(pair[1])
        return (
            (r1 == "Experimental" and r2 in self._CONTROL_ROLES) or
            (r2 == "Experimental" and r1 in self._CONTROL_ROLES)
        )

    def _rebuild_pair_list(self):
        # Preserve current check state across rebuild (only for pairs that still exist)
        prior_state = {pair: cb.isChecked() for pair, cb in self._pair_checkboxes.items()}

        # Clear widgets
        for i in reversed(range(self._scroll_layout.count())):
            w = self._scroll_layout.itemAt(i).widget()
            if w is not None:
                w.deleteLater()
        self._pair_checkboxes.clear()

        # Enumerate new pairs (canonicalized, sorted for stable UI order)
        conds = sorted(self._conditions.keys())
        pairs = [self._canon(p) for p in combinations(conds, 2)]

        for pair in pairs:
            cb = QCheckBox(f"{pair[0]}   vs   {pair[1]}")
            # Precedence for initial state:
            #   1. prior state (if pair existed before rebuild)
            #   2. pending excluded (from load_config before conditions arrived)
            #   3. role-based default: experimental vs control → checked
            if pair in prior_state:
                cb.setChecked(prior_state[pair])
            elif pair in self._pending_excluded:
                cb.setChecked(False)
            else:
                cb.setChecked(self._is_default_pair(pair))
            cb.stateChanged.connect(self._on_checkbox_changed)
            self._scroll_layout.addWidget(cb)
            self._pair_checkboxes[pair] = cb

        # Pending excluded was consumed (or stays for pairs not yet present)
        # Clearing on every rebuild is fine: if load_config is re-called, it repopulates.
        self._pending_excluded = {e for e in self._pending_excluded if e not in set(pairs)}

        self._apply_filter()
        self._update_summary()

    def _on_checkbox_changed(self):
        self._user_has_interacted = True
        self._update_summary()
        self.selection_changed.emit()

    def _update_summary(self):
        total   = len(self._pair_checkboxes)
        checked = sum(1 for cb in self._pair_checkboxes.values() if cb.isChecked())
        if total == 0:
            self._summary_label.setText(
                "<i>Define conditions in the Experimental Conditions tab to populate "
                "pairwise comparisons.</i>"
            )
            return
        self._summary_label.setText(
            f"<b>{checked}</b> of <b>{total}</b> comparisons selected "
            f"&nbsp;·&nbsp; Sidak family size per phase group: <b>{checked}</b>"
        )

    def _select_all(self):
        self._user_has_interacted = True
        for cb in self._pair_checkboxes.values():
            cb.setChecked(True)

    def _deselect_all(self):
        self._user_has_interacted = True
        for cb in self._pair_checkboxes.values():
            cb.setChecked(False)

    def _apply_filter(self):
        query = self._filter_edit.text().strip().lower()
        for cb in self._pair_checkboxes.values():
            cb.setVisible(not query or query in cb.text().lower())

    def _on_pe_toggle(self, _state):
        self._update_pe_settings_enabled()
        self.selection_changed.emit()

    def _update_pe_settings_enabled(self):
        on = self._pe_enabled_cb.isChecked()
        for w in self._pe_settings_widgets:
            w.setEnabled(on)

import sys
import os
import pandas as pd
from functools import partial
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QFileDialog, QLabel, QMessageBox, QGroupBox,
    QTreeWidget, QTreeWidgetItem, QHeaderView, QSplitter, QTableWidget, QTableWidgetItem,
    QAbstractItemView, QTabWidget, QLineEdit, QCheckBox, QSpinBox, QColorDialog, QComboBox
)
from PySide6.QtCore import QDir, QDirIterator, Qt, Signal, QSize, QObject, QMetaMethod
from PySide6.QtGui import QIntValidator, QColor, QPalette, QPixmap, QPainter, QIcon

# =============================================================================
#  EXPERIMENTAL DESIGN WIDGETS
# =============================================================================

class ExperimentalSetupWidget(QWidget):
    """Sub-tab for defining experimental phases."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Experimental Setup")
        self._is_recalculating = False
        self.data_for_setup = None

        layout = QVBoxLayout(self)

        # Add a label to display the maximum time
        self.max_time_label = QLabel("Max time from data: N/A")
        layout.addWidget(self.max_time_label)

        control_layout = QHBoxLayout()
        control_layout.addWidget(QLabel("Number of Phases:"))
        self.num_phases_spinbox = QSpinBox()
        self.num_phases_spinbox.setRange(1, 100)
        self.num_phases_spinbox.setValue(7)
        control_layout.addWidget(self.num_phases_spinbox)
        control_layout.addStretch()
        layout.addLayout(control_layout)

        self.table = QTableWidget()
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels(["Group", "Phase Name", "Length (s)", "Start (s)", "End (s)"])
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        layout.addWidget(self.table)

        # --- Connections ---
        self.num_phases_spinbox.valueChanged.connect(self.generate_phase_rows)

        # --- Initial State ---
        self.generate_phase_rows()

    def load_data(self, data_tuple):
        # Unpack the tuple, now expecting max_time_sec as the fourth element
        final_df_with_notation, file_dataframes, suffix_map, max_time_sec = data_tuple
        self.max_time_label.setText(f"Max time from data: {max_time_sec:.2f} seconds")

    def generate_phase_rows(self):
        """
        Adjusts the number of rows while preserving all existing data.
        Existing rows are left untouched; only genuinely new rows get
        fresh default widgets.
        """
        num_rows     = self.num_phases_spinbox.value()
        current_rows = self.table.rowCount()

        self.table.cellChanged.disconnect() if self.isSignalConnected(self.table, "cellChanged") else None

        if num_rows > current_rows:
            # ── Adding rows: only touch the new ones ──────────────────────────
            self.table.setRowCount(num_rows)
            for row in range(current_rows, num_rows):
                group_spin = QSpinBox()
                group_spin.setRange(1, 99)
                self.table.setCellWidget(row, 0, group_spin)

                self.table.setItem(row, 1, QTableWidgetItem(f"Phase {row + 1}"))

                length_spin = QSpinBox()
                length_spin.setRange(0, 100000)
                length_spin.setValue(300)
                length_spin.valueChanged.connect(self.recalculate_phases)
                self.table.setCellWidget(row, 2, length_spin)

                self.table.setItem(row, 3, QTableWidgetItem("0"))
                self.table.setItem(row, 4, QTableWidgetItem("300"))

        elif num_rows < current_rows:
            # ── Removing rows: Qt drops the last N rows automatically ─────────
            self.table.setRowCount(num_rows)

        # Reconnect signals for any length spinboxes that lost their connection
        # (only needed for pre-existing rows if Qt recreated any widgets)
        for row in range(self.table.rowCount()):
            spin = self.table.cellWidget(row, 2)
            if spin is not None:
                try:
                    spin.valueChanged.disconnect(self.recalculate_phases)
                except RuntimeError:
                    pass
                spin.valueChanged.connect(self.recalculate_phases)

        self.recalculate_phases()
        self.table.cellChanged.connect(self.manual_time_edit)

    def recalculate_phases(self, _=None):
        if self._is_recalculating:
            return
        self._is_recalculating = True

        current_start_time = 0
        for i in range(self.table.rowCount()):
            length = self.table.cellWidget(i, 2).value()

            start_item = self.table.item(i, 3)
            if start_item is None:
                start_item = QTableWidgetItem()
                self.table.setItem(i, 3, start_item)
            start_item.setText(str(current_start_time))

            end_time = current_start_time + length
            end_item = self.table.item(i, 4)
            if end_item is None:
                end_item = QTableWidgetItem()
                self.table.setItem(i, 4, end_item)
            end_item.setText(str(end_time))

            current_start_time = end_time

        self._is_recalculating = False

    def manual_time_edit(self, row, column):
        if self._is_recalculating or (column != 3 and column != 4): # Adjusted for new column indices
            return

        self._is_recalculating = True

        try:
            if column == 4: # End time was changed
                end_time = int(self.table.item(row, 4).text())
                start_time = int(self.table.item(row, 3).text())
                if end_time < start_time:
                    end_time = start_time
                    self.table.item(row, 4).setText(str(end_time))

                new_length = end_time - start_time
                self.table.cellWidget(row, 2).setValue(new_length)

            elif column == 3: # Start time was changed
                start_time = int(self.table.item(row, 3).text())
                end_time_item = self.table.item(row, 4)
                if end_time_item:
                    end_time = int(end_time_item.text())
                    if start_time > end_time:
                        start_time = end_time
                        self.table.item(row, 3).setText(str(start_time))

                    new_length = end_time - start_time
                    self.table.cellWidget(row, 2).setValue(new_length)

        except (ValueError, AttributeError):
            pass
        finally:
            self._is_recalculating = False
            self.recalculate_phases()

    def get_phases_data(self):
        """Returns a list of dictionaries with phase information."""
        phases_data = []
        for row in range(self.table.rowCount()):
            try:
                group = self.table.cellWidget(row, 0).value()
                phase_name = self.table.item(row, 1).text()
                length = self.table.cellWidget(row, 2).value()
                start_time = int(self.table.item(row, 3).text())
                end_time = int(self.table.item(row, 4).text())

                phases_data.append({
                    "Group": group,
                    "Phase": phase_name,
                    "Length_s": length,
                    "Start_s": start_time,
                    "End_s": end_time,
                })
            except (AttributeError, ValueError) as e:
                print(f"Could not get data for row {row}: {e}")
        return phases_data

    def load_phases_config(self, phases_data: list):
        """Restores the phase table from a list of phase dicts (as returned by get_phases_data)."""
        if not phases_data:
            return
        self.num_phases_spinbox.setValue(len(phases_data))
        # generate_phase_rows is triggered by the spinbox; wait for rows to exist
        for row, phase in enumerate(phases_data):
            try:
                self.table.cellWidget(row, 0).setValue(phase.get("Group", 1))
                self.table.item(row, 1).setText(phase.get("Phase", f"Phase {row + 1}"))
                self.table.cellWidget(row, 2).setValue(phase.get("Length_s", 300))
            except (AttributeError, TypeError):
                pass
        self.recalculate_phases()

    def isSignalConnected(self, obj, name):
        """Check if a signal is connected."""
        if isinstance(obj, QObject):
            meta_obj = obj.metaObject()
            for i in range(meta_obj.methodCount()):
                meta_method = meta_obj.method(i)
                if meta_method.methodType() == QMetaMethod.Signal and meta_method.name() == name:
                    return obj.isSignalConnected(meta_method)
        return False


# =============================================================================
#  EXPERIMENTAL CONDITIONS WIDGET (multi-variable descriptors)
# =============================================================================

class ExperimentalConditionsWidget(QWidget):
    """Sub-tab for defining experimental conditions with multiple descriptor variables.

    Each condition is defined by values for N input variables (e.g. Genotype, Drug,
    Antidote).  All combinatorial interaction columns are auto-computed:
      - For 2 vars (V1, V2): one combo column V1+V2
      - For 3 vars (V1, V2, V3): combos V1+V2, V1+V3, V2+V3, V1+V2+V3
      - For N vars: all C(N,2) + C(N,3) + ... + C(N,N) combination columns

    The full interaction column (all variables combined) is the "Condition" used
    downstream by the plate format, GAMM, and catch22 widgets.
    """
    conditions_updated = Signal(dict)   # {full_interaction_name: hex_color}
    variables_updated  = Signal(list)   # [var_name_str, ...]

    DEFAULT_VAR_NAMES = ["Variable 1", "Variable 2"]

    def __init__(self, parent=None):
        super().__init__(parent)
        self._num_vars = 2
        self._var_name_edits = []
        self._updating = False          # guard against recursive signals
        # Column layout tracking (set by _setup_table_columns)
        self._combo_specs = []          # list of tuples of var indices for each combo column
        self._combo_col_start = 0       # first table column index for combo columns
        self._color_col = 0             # table column index for the Color button

        layout = QVBoxLayout(self)

        self._load_palettes()

        # ── Row 1: Variable controls ─────────────────────────────────────────
        var_row = QHBoxLayout()
        var_row.addWidget(QLabel("Number of Variables:"))
        self.num_vars_spinbox = QSpinBox()
        self.num_vars_spinbox.setRange(1, 10)
        self.num_vars_spinbox.setValue(2)
        var_row.addWidget(self.num_vars_spinbox)
        var_row.addSpacing(20)
        var_row.addWidget(QLabel("Variable Names:"))
        self._var_names_layout = QHBoxLayout()
        var_row.addLayout(self._var_names_layout)
        var_row.addStretch()
        layout.addLayout(var_row)

        # ── Row 2: Condition count + palette ─────────────────────────────────
        control_layout = QHBoxLayout()
        control_layout.addWidget(QLabel("Number of Conditions:"))
        self.num_conditions_spinbox = QSpinBox()
        self.num_conditions_spinbox.setRange(1, 100)
        self.num_conditions_spinbox.setValue(4)
        control_layout.addWidget(self.num_conditions_spinbox)

        control_layout.addSpacing(20)

        control_layout.addWidget(QLabel("Select Palette:"))
        self.palette_combo = QComboBox()

        swatch_width, swatch_height = 120, 18
        self.palette_combo.setIconSize(QSize(swatch_width, swatch_height))

        for name in sorted(self.palettes.keys()):
            colors = self.palettes[name]
            icon = self._create_palette_swatch(colors, width=swatch_width, height=swatch_height)
            self.palette_combo.addItem(icon, name)

        control_layout.addWidget(self.palette_combo)
        control_layout.addStretch()
        layout.addLayout(control_layout)

        # ── Conditions table ─────────────────────────────────────────────────
        self.conditions_table = QTableWidget()
        layout.addWidget(self.conditions_table)

        # ── Build initial state ──────────────────────────────────────────────
        self._build_var_name_editors()
        self._setup_table_columns()

        # ── Connections ──────────────────────────────────────────────────────
        self.num_vars_spinbox.valueChanged.connect(self._on_var_count_changed)
        self.num_conditions_spinbox.valueChanged.connect(self.generate_condition_rows)
        self.palette_combo.currentTextChanged.connect(self.generate_condition_rows)
        self.conditions_table.itemChanged.connect(self._on_cell_changed)

        # ── Initial state ────────────────────────────────────────────────────
        self.generate_condition_rows()
        self.emit_conditions_data()

    # ── Combinatorial helpers ────────────────────────────────────────────────

    def _get_combo_specs(self):
        """Returns a list of index-tuples for all combination columns (size 2..N).

        For N=1: no combos (only one variable).
        For N=2: [(0,1)]
        For N=3: [(0,1), (0,2), (1,2), (0,1,2)]
        """
        specs = []
        for k in range(2, self._num_vars + 1):
            for combo in combinations(range(self._num_vars), k):
                specs.append(combo)
        return specs

    # ── Variable name editors ────────────────────────────────────────────────

    def _build_var_name_editors(self):
        """Creates QLineEdit widgets for each variable name."""
        while self._var_names_layout.count():
            item = self._var_names_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        self._var_name_edits.clear()

        for i in range(self._num_vars):
            edit = QLineEdit()
            edit.setMaximumWidth(120)
            if i < len(self.DEFAULT_VAR_NAMES):
                edit.setText(self.DEFAULT_VAR_NAMES[i])
            else:
                edit.setText(f"Variable {i + 1}")
            edit.textChanged.connect(self._on_var_name_changed)
            self._var_names_layout.addWidget(edit)
            self._var_name_edits.append(edit)

    def _setup_table_columns(self):
        """Sets up table columns: [V1] [V2] ... [VN] [combo cols...] [Color]."""
        self._combo_specs = self._get_combo_specs()
        n_combos = len(self._combo_specs)
        self._combo_col_start = self._num_vars
        self._color_col = self._num_vars + n_combos

        # Fully clear old table contents and widgets before resizing
        self._clear_all_cell_widgets()
        self.conditions_table.setRowCount(0)

        col_count = self._num_vars + n_combos + 1  # vars + combos + Color
        self.conditions_table.setColumnCount(col_count)

        var_names = self.get_variable_names()
        headers = list(var_names)
        for spec in self._combo_specs:
            combo_name = "+".join(var_names[i] for i in spec)
            headers.append(combo_name)
        headers.append("Color")
        self.conditions_table.setHorizontalHeaderLabels(headers)

        # Stretch the full interaction column (last combo) if there are combos
        if n_combos > 0:
            self.conditions_table.horizontalHeader().setSectionResizeMode(
                self._color_col - 1, QHeaderView.Stretch
            )

    # ── Variable count / name changes ────────────────────────────────────────

    def _on_var_count_changed(self, new_count):
        """Called when the Number of Variables spinbox changes."""
        if new_count == self._num_vars:
            return

        saved = self._snapshot_table()
        old_count = self._num_vars

        self._num_vars = new_count
        self._build_var_name_editors()
        self._setup_table_columns()

        self._restore_and_generate(saved, old_count)

        self.variables_updated.emit(self.get_variable_names())
        self.emit_conditions_data()

    def _on_var_name_changed(self):
        """Called when any variable name QLineEdit changes."""
        var_names = self.get_variable_names()
        headers = list(var_names)
        for spec in self._combo_specs:
            combo_name = "+".join(var_names[i] for i in spec)
            headers.append(combo_name)
        headers.append("Color")
        for i, h in enumerate(headers):
            item = QTableWidgetItem(h)
            self.conditions_table.setHorizontalHeaderItem(i, item)
        self.variables_updated.emit(var_names)

    def _snapshot_table(self):
        """Saves variable values and colors for all rows."""
        rows = []
        for row in range(self.conditions_table.rowCount()):
            vals = []
            for col in range(self._num_vars):
                item = self.conditions_table.item(row, col)
                vals.append(item.text() if item else "")
            btn = self.conditions_table.cellWidget(row, self._color_col)
            color = "#ffffff"
            if btn:
                ss = btn.styleSheet()
                if "background-color:" in ss:
                    color = ss.split("background-color:")[1].split(";")[0].strip()
            rows.append({"vals": vals, "color": color})
        return rows

    def _restore_and_generate(self, saved, old_var_count):
        """Rebuilds rows, restoring saved variable values where possible."""
        self.generate_condition_rows()

        self._updating = True
        for row in range(min(len(saved), self.conditions_table.rowCount())):
            old_vals = saved[row]["vals"]
            for col in range(min(old_var_count, self._num_vars)):
                item = self.conditions_table.item(row, col)
                if item and col < len(old_vals) and old_vals[col]:
                    item.setText(old_vals[col])
            btn = self.conditions_table.cellWidget(row, self._color_col)
            if btn:
                btn.setStyleSheet(
                    f"background-color: {saved[row]['color']}; border: 1px solid #888;"
                )
            self._update_combo_cells(row)
        self._updating = False

    # ── Table row generation ─────────────────────────────────────────────────

    def _clear_all_cell_widgets(self):
        """Removes every cell widget in the table to prevent stale buttons."""
        for r in range(self.conditions_table.rowCount()):
            for c in range(self.conditions_table.columnCount()):
                self.conditions_table.removeCellWidget(r, c)

    def generate_condition_rows(self):
        """Generates/refreshes condition rows with palette colors."""
        num_rows = self.num_conditions_spinbox.value()
        if num_rows <= 0:
            self.conditions_table.setRowCount(0)
            return

        current_palette_name = self.palette_combo.currentText()
        anchor_colors = self.palettes.get(current_palette_name, [])

        if not anchor_colors:
            hex_colors = [mcolors.to_hex(c) for c in np.linspace(0.2, 0.8, num_rows)]
        elif num_rows <= len(anchor_colors):
            indices = np.linspace(0, len(anchor_colors) - 1, num_rows, dtype=int)
            hex_colors = [anchor_colors[i] for i in indices]
        else:
            cmap = mcolors.LinearSegmentedColormap.from_list("custom_gradient", anchor_colors)
            generated_colors = cmap(np.linspace(0, 1, num_rows))
            hex_colors = [mcolors.to_hex(c) for c in generated_colors]

        self._updating = True
        # Clear stale cell widgets (e.g. old color buttons at previous column positions)
        self._clear_all_cell_widgets()
        self.conditions_table.setRowCount(num_rows)

        var_names = self.get_variable_names()
        # Column-index letters for defaults: Variable 1→A, Variable 2→B, ...
        _LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        for row in range(num_rows):
            # Variable value columns — set defaults only for new (empty) cells
            for col in range(self._num_vars):
                existing = self.conditions_table.item(row, col)
                if not existing:
                    prefix = _LETTERS[col] if col < len(_LETTERS) else f"V{col+1}_"
                    default = f"{prefix}{row + 1}"
                    self.conditions_table.setItem(row, col, QTableWidgetItem(default))

            # Combo columns (auto-computed, read-only)
            self._update_combo_cells(row)

            # Color button
            color = hex_colors[row]
            color_button = QPushButton()
            color_button.setStyleSheet(f"background-color: {color}; border: 1px solid #888;")
            color_button.clicked.connect(partial(self.open_color_picker, row))
            self.conditions_table.setCellWidget(row, self._color_col, color_button)

        self._updating = False
        self.emit_conditions_data()

    def _update_combo_cells(self, row):
        """Recomputes ALL auto-generated combination columns for a row."""
        # Read variable values
        var_vals = []
        for col in range(self._num_vars):
            item = self.conditions_table.item(row, col)
            var_vals.append(item.text() if item else "")

        # Fill each combination column
        for ci, spec in enumerate(self._combo_specs):
            combo_val = "+".join(var_vals[i] for i in spec)
            table_col = self._combo_col_start + ci
            item = self.conditions_table.item(row, table_col)
            if not item:
                item = QTableWidgetItem()
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.conditions_table.setItem(row, table_col, item)
            else:
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            item.setText(combo_val)

    def _on_cell_changed(self, item):
        """Called when any table cell is edited (itemChanged emits QTableWidgetItem)."""
        if self._updating or item is None:
            return
        row = item.row()
        col = item.column()
        if col < self._num_vars:
            self._updating = True
            self._update_combo_cells(row)
            self._updating = False
            self.emit_conditions_data()

    # ── Data emission ────────────────────────────────────────────────────────

    def _full_interaction_col(self):
        """Returns the table column index of the full interaction (last combo)."""
        if self._combo_specs:
            return self._combo_col_start + len(self._combo_specs) - 1
        # If only 1 variable (no combos), treat column 0 as the condition
        return 0

    def emit_conditions_data(self):
        """Gathers {full_interaction_name: hex_color} and emits conditions_updated."""
        if self._updating:
            return
        conditions = {}
        fi_col = self._full_interaction_col()
        for row in range(self.conditions_table.rowCount()):
            cond_item = self.conditions_table.item(row, fi_col)
            color_btn = self.conditions_table.cellWidget(row, self._color_col)
            if cond_item and color_btn:
                name = cond_item.text()
                ss = color_btn.styleSheet()
                color = ss.split("background-color:")[1].split(";")[0].strip()
                conditions[name] = color
        self.conditions_updated.emit(conditions)

    # ── Public getters ───────────────────────────────────────────────────────

    def get_variable_names(self):
        """Returns list of variable names, e.g. ['Variable 1', 'Variable 2']."""
        return [edit.text() for edit in self._var_name_edits]

    def get_condition_variables(self):
        """Returns {condition_name: {var_name: value}} for all conditions."""
        var_names = self.get_variable_names()
        fi_col = self._full_interaction_col()
        result = {}
        for row in range(self.conditions_table.rowCount()):
            cond_item = self.conditions_table.item(row, fi_col)
            if not cond_item:
                continue
            cond_name = cond_item.text()
            var_dict = {}
            for col, var_name in enumerate(var_names):
                item = self.conditions_table.item(row, col)
                var_dict[var_name] = item.text() if item else ""
            result[cond_name] = var_dict
        return result

    def get_conditions(self):
        """Returns {condition_name: hex_color}."""
        conditions = {}
        fi_col = self._full_interaction_col()
        for row in range(self.conditions_table.rowCount()):
            cond_item = self.conditions_table.item(row, fi_col)
            color_btn = self.conditions_table.cellWidget(row, self._color_col)
            if cond_item and color_btn:
                name = cond_item.text()
                ss = color_btn.styleSheet()
                color = ss.split("background-color:")[1].split(";")[0].strip()
                conditions[name] = color
        return conditions

    # ── Save / Load config ───────────────────────────────────────────────────

    def save_conditions_config(self):
        """Returns a dict for JSON serialisation (all user inputs)."""
        var_names = self.get_variable_names()
        rows = []
        for row in range(self.conditions_table.rowCount()):
            vals = []
            for col in range(self._num_vars):
                item = self.conditions_table.item(row, col)
                vals.append(item.text() if item else "")
            btn = self.conditions_table.cellWidget(row, self._color_col)
            color = "#ffffff"
            if btn:
                ss = btn.styleSheet()
                if "background-color:" in ss:
                    color = ss.split("background-color:")[1].split(";")[0].strip()
            rows.append({"values": vals, "color": color})
        return {
            "variable_names": var_names,
            "num_conditions": self.num_conditions_spinbox.value(),
            "palette": self.palette_combo.currentText(),
            "rows": rows,
        }

    def load_conditions_config(self, config):
        """Restores from config dict.  Supports new format and legacy {name: color} format."""
        if isinstance(config, dict) and "variable_names" in config:
            var_names = config["variable_names"]
            rows = config.get("rows", [])

            self.num_vars_spinbox.setValue(len(var_names))
            for i, name in enumerate(var_names):
                if i < len(self._var_name_edits):
                    self._var_name_edits[i].setText(name)

            palette = config.get("palette", "")
            if palette:
                idx = self.palette_combo.findText(palette)
                if idx >= 0:
                    self.palette_combo.setCurrentIndex(idx)

            self.num_conditions_spinbox.setValue(config.get("num_conditions", len(rows)))

            self._updating = True
            for row_idx, row_data in enumerate(rows):
                vals = row_data.get("values", [])
                color = row_data.get("color", "#ffffff")
                for col, val in enumerate(vals):
                    if col < self._num_vars:
                        item = self.conditions_table.item(row_idx, col)
                        if item:
                            item.setText(val)
                self._update_combo_cells(row_idx)
                btn = self.conditions_table.cellWidget(row_idx, self._color_col)
                if btn:
                    btn.setStyleSheet(f"background-color: {color}; border: 1px solid #888;")
            self._updating = False
        elif isinstance(config, dict):
            # Legacy format: plain {condition_name: hex_color}
            items = list(config.items())
            self.num_conditions_spinbox.setValue(len(items))

            self._updating = True
            for row, (name, color) in enumerate(items):
                parts = name.split("+")
                if len(parts) > self._num_vars:
                    self.num_vars_spinbox.setValue(len(parts))
                for col, val in enumerate(parts):
                    if col < self._num_vars:
                        item = self.conditions_table.item(row, col)
                        if item:
                            item.setText(val)
                self._update_combo_cells(row)
                btn = self.conditions_table.cellWidget(row, self._color_col)
                if btn:
                    btn.setStyleSheet(f"background-color: {color}; border: 1px solid #888;")
            self._updating = False

        self.emit_conditions_data()

    # ── Color picker ─────────────────────────────────────────────────────────

    def open_color_picker(self, row):
        """Opens a color dialog and updates the button for the specified row."""
        button = self.conditions_table.cellWidget(row, self._color_col)
        if not button:
            return
        stylesheet = button.styleSheet()
        current_color_hex = "#ffffff"
        if "background-color:" in stylesheet:
            current_color_hex = stylesheet.split("background-color:")[1].split(";")[0].strip()
        current_color = QColor(current_color_hex)
        color = QColorDialog.getColor(current_color, self, "Select Condition Color")
        if color.isValid():
            button.setStyleSheet(f"background-color: {color.name()}; border: 1px solid #888;")
            self.emit_conditions_data()

    # ── Palette swatch helper ────────────────────────────────────────────────

    def _create_palette_swatch(self, colors, width=120, height=18):
        """Creates a QIcon swatch from a list of hex color strings."""
        pixmap = QPixmap(width, height)
        pixmap.fill(Qt.transparent)
        painter = QPainter(pixmap)
        num_colors = len(colors)
        if num_colors == 0:
            painter.end()
            return QIcon(pixmap)
        rect_width = width / num_colors
        for i, hex_color in enumerate(colors):
            color = QColor(hex_color)
            painter.setBrush(color)
            painter.setPen(Qt.NoPen)
            painter.drawRect(int(i * rect_width), 0, int(rect_width + 1), height)
        painter.end()
        return QIcon(pixmap)

    # ── Palette loading (unchanged) ──────────────────────────────────────────

    def _load_palettes(self):
        """Loads a comprehensive set of color palettes for up to 30+ conditions.

        Sources: seaborn, matplotlib, distinctipy, glasbey, colorcet, bokeh,
        palettable, cmocean, cmasher, husl, vapeplot, plotly (hardcoded),
        and cycler.
        """
        self.palettes = {}
        N = 30  # pre-generate this many colors for continuous palettes

        # ══════════════════════════════════════════════════════════════════════
        #  CATEGORICAL — PURPOSE-BUILT FOR 30+ DISTINCT GROUPS
        # ══════════════════════════════════════════════════════════════════════

        # distinctipy — maximally distinct colors in perceptual space
        try:
            import distinctipy
            cols = distinctipy.get_colors(N)
            self.palettes["Distinct: distinctipy (30)"] = [mcolors.to_hex(c) for c in cols]
        except Exception:
            pass

        # glasbey — optimised categorical palettes
        try:
            import glasbey
            cols = glasbey.create_palette(palette_size=N)
            self.palettes["Distinct: glasbey (30)"] = list(cols)
        except Exception:
            pass

        # colorcet — perceptually uniform categorical
        try:
            import colorcet as cc
            for name, attr in [("glassbey", "glasbey"), ("glassbey_light", "glasbey_light"),
                                ("glassbey_dark", "glasbey_dark"), ("glassbey_warm", "glasbey_warm"),
                                ("glassbey_cool", "glasbey_cool")]:
                pal = getattr(cc, attr, None)
                if pal:
                    self.palettes[f"Colorcet: {name}"] = [c for c in pal[:N]]
        except Exception:
            pass

        # ══════════════════════════════════════════════════════════════════════
        #  GENERAL PURPOSE — LARGE BUILT-IN QUALITATIVE SETS
        # ══════════════════════════════════════════════════════════════════════

        # Matplotlib qualitative
        for name in ['tab10', 'tab20', 'tab20b', 'tab20c',
                      'Set1', 'Set2', 'Set3', 'Paired', 'Accent',
                      'Dark2', 'Pastel1', 'Pastel2']:
            try:
                cmap = plt.get_cmap(name)
                colors = cmap.colors if hasattr(cmap, 'colors') else cmap(np.linspace(0, 1, 20))
                self.palettes[f"Qualitative: {name}"] = [mcolors.to_hex(c) for c in colors]
            except Exception:
                pass

        # Matplotlib HSV (continuous, sample N)
        try:
            cmap = plt.get_cmap('hsv')
            colors = cmap(np.linspace(0, 0.92, N))  # avoid wrapping back to red
            self.palettes["Qualitative: hsv"] = [mcolors.to_hex(c) for c in colors]
        except Exception:
            pass

        # Plotly qualitative (hardcoded — avoids 50 MB plotly dependency)
        self.palettes["Plotly: Alphabet (26)"] = [
            '#AA0DFE', '#3283FE', '#85660D', '#782AB6', '#565656',
            '#1C8356', '#16FF32', '#F7E1A0', '#E2E2E2', '#1CBE4F',
            '#C4451C', '#DEA0FD', '#FE00FA', '#325A9B', '#FEAF16',
            '#F8A19F', '#90AD1C', '#F6222E', '#1CFFCE', '#2ED9FF',
            '#B10DA1', '#C075A6', '#FC1CBF', '#B00068', '#FBE426',
            '#FA0087',
        ]
        self.palettes["Plotly: Dark24"] = [
            '#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF',
            '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB',
            '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32',
            '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7',
            '#DA60CA', '#6C4516', '#0D2A63', '#AF0038',
        ]
        self.palettes["Plotly: Light24"] = [
            '#FD3216', '#00FE35', '#6A76FC', '#FED4C4', '#FE00CE',
            '#0DF9FF', '#F6F926', '#FF9616', '#DEA0FD', '#922B21',
            '#9467BD', '#563098', '#FC5A50', '#B00068', '#FF00FF',
            '#5A5156', '#1CBE4F', '#C4451C', '#FBE426', '#FA0087',
            '#325A9B', '#B10DA1', '#16FF32', '#3283FE',
        ]

        # Bokeh Category20
        try:
            from bokeh.palettes import Category20_20
            self.palettes["Bokeh: Category20"] = list(Category20_20)
        except Exception:
            pass

        # ══════════════════════════════════════════════════════════════════════
        #  SEABORN QUALITATIVE & GENERATED
        # ══════════════════════════════════════════════════════════════════════

        for name in ['deep', 'muted', 'pastel', 'bright', 'dark', 'colorblind']:
            palette = sns.color_palette(name)
            self.palettes[f"Seaborn: {name.capitalize()}"] = [mcolors.to_hex(c) for c in palette]

        # HUSL and HLS — scale to any N with evenly spaced hues
        for name in ['husl', 'hls']:
            palette = sns.color_palette(name, n_colors=N)
            self.palettes[f"Seaborn: {name.upper()} ({N})"] = [mcolors.to_hex(c) for c in palette]

        # Seaborn cubehelix
        try:
            palette = sns.cubehelix_palette(n_colors=N)
            self.palettes[f"Seaborn: Cubehelix ({N})"] = [mcolors.to_hex(c) for c in palette]
        except Exception:
            pass

        # ══════════════════════════════════════════════════════════════════════
        #  PALETTABLE — ColorBrewer, Cartocolors, Scientific, etc.
        # ══════════════════════════════════════════════════════════════════════

        try:
            from palettable.cartocolors.qualitative import Bold_10, Prism_10, Vivid_10, Safe_10, Antique_10
            for pal, label in [(Bold_10, "Bold"), (Prism_10, "Prism"),
                                (Vivid_10, "Vivid"), (Safe_10, "Safe"), (Antique_10, "Antique")]:
                self.palettes[f"Palettable: {label}"] = [mcolors.to_hex([c/255 for c in rgb]) for rgb in pal.colors]
        except Exception:
            pass

        try:
            from palettable.scientific.diverging import Roma_20, Vik_20, Berlin_20
            for pal, label in [(Roma_20, "Roma"), (Vik_20, "Vik"), (Berlin_20, "Berlin")]:
                self.palettes[f"Palettable: {label}"] = [mcolors.to_hex([c/255 for c in rgb]) for rgb in pal.colors]
        except Exception:
            pass

        # ══════════════════════════════════════════════════════════════════════
        #  SCIENTIFIC COLORMAPS — cmocean, cmasher
        # ══════════════════════════════════════════════════════════════════════

        try:
            import cmocean
            for name in ['thermal', 'haline', 'solar', 'ice', 'deep',
                          'dense', 'algae', 'matter', 'turbid', 'speed',
                          'amp', 'tempo', 'rain', 'phase', 'balance', 'delta', 'curl']:
                cmap = getattr(cmocean.cm, name, None)
                if cmap:
                    colors = cmap(np.linspace(0.05, 0.95, N))
                    self.palettes[f"cmocean: {name}"] = [mcolors.to_hex(c) for c in colors]
        except Exception:
            pass

        try:
            import cmasher as cmr
            for name in ['ember', 'flamingo', 'gothic', 'jungle',
                          'lavender', 'lilac', 'neon', 'ocean', 'pepper',
                          'rainforest', 'savanna', 'sepia', 'sunburst',
                          'torch', 'tropical', 'voltage', 'wildfire',
                          'chroma', 'pride', 'cosmic', 'freeze', 'fusion',
                          'guppy', 'holly', 'iceburn', 'redshift', 'watermelon']:
                try:
                    cmap = cmr.get_sub_cmap(f"cmr.{name}", 0.05, 0.95)
                    colors = cmap(np.linspace(0, 1, N))
                    self.palettes[f"cmasher: {name}"] = [mcolors.to_hex(c) for c in colors]
                except Exception:
                    pass
        except Exception:
            pass

        # ══════════════════════════════════════════════════════════════════════
        #  AESTHETIC / THEMED
        # ══════════════════════════════════════════════════════════════════════

        try:
            import vapeplot
            for name, pal in vapeplot.palettes.items():
                if pal:
                    self.palettes[f"Vapeplot: {name}"] = [mcolors.to_hex(c) for c in pal]
        except Exception:
            pass

        # ══════════════════════════════════════════════════════════════════════
        #  PERCEPTUALLY UNIFORM SEQUENTIAL
        # ══════════════════════════════════════════════════════════════════════

        for name in ['viridis', 'plasma', 'inferno', 'magma', 'cividis',
                      'turbo', 'rainbow', 'jet']:
            try:
                cmap = plt.get_cmap(name)
                colors = cmap(np.linspace(0.05, 0.95, N))
                self.palettes[f"Sequential: {name.capitalize()}"] = [mcolors.to_hex(c) for c in colors]
            except Exception:
                pass

        # Bokeh Turbo256 and Viridis256
        try:
            from bokeh.palettes import Turbo256, Viridis256
            step = max(1, len(Turbo256) // N)
            self.palettes["Bokeh: Turbo256"] = list(Turbo256[::step])[:N]
            self.palettes["Bokeh: Viridis256"] = list(Viridis256[::step])[:N]
        except Exception:
            pass

        # ══════════════════════════════════════════════════════════════════════
        #  DIVERGING
        # ══════════════════════════════════════════════════════════════════════

        for name in ['coolwarm', 'RdYlBu', 'RdYlGn', 'Spectral',
                      'PiYG', 'PRGn', 'BrBG', 'RdBu']:
            try:
                cmap = plt.get_cmap(name)
                colors = cmap(np.linspace(0.05, 0.95, N))
                self.palettes[f"Diverging: {name}"] = [mcolors.to_hex(c) for c in colors]
            except Exception:
                pass


# =============================================================================
#  METADATA ASSIGNMENT WIDGET (replaces ReferenceControlsWidget)
# =============================================================================

class MetadataAssignmentWidget(QWidget):
    """Sub-tab for assigning experimental roles to conditions and auto-deriving reference levels.

    Roles
    -----
    Reference Control — baseline for all statistical comparisons
    Experimental      — treatments being tested
    Positive Control  — known responders (validation)
    Potential Rescue  — rescue / combination conditions
    """
    roles_updated      = Signal(dict)        # {condition_name: role_str}
    references_updated = Signal(dict, str)   # ({var_name: ref_value}, ref_condition)

    ROLES = ["Reference Control", "Experimental", "Positive Control", "Potential Rescue"]
    ROLE_ABBREVIATIONS = {
        "Reference Control": "RC",
        "Experimental":      "EXP",
        "Positive Control":  "PC",
        "Potential Rescue":  "RES",
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self._variable_names       = ["Genotype", "Drug"]
        self._conditions           = {}    # {name: color}
        self._condition_variables  = {}    # {name: {var: val}}

        layout = QVBoxLayout(self)
        layout.setAlignment(Qt.AlignTop)

        layout.addWidget(QLabel(
            "<b>Assign a role to each condition for statistical analysis.</b><br>"
            "The <i>Reference Control</i> condition determines the baseline for all comparisons."
        ))
        layout.addSpacing(8)

        # ── Role assignment table ────────────────────────────────────────────
        self._role_table = QTableWidget()
        self._role_table.setColumnCount(2)
        self._role_table.setHorizontalHeaderLabels(["Condition", "Role"])
        self._role_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self._role_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        layout.addWidget(self._role_table)

        layout.addSpacing(12)

        # ── Auto-derived reference levels ────────────────────────────────────
        self._ref_group = QGroupBox("Auto-Derived Reference Levels (from Reference Control)")
        ref_layout = QVBoxLayout(self._ref_group)
        self._ref_info_label = QLabel("No Reference Control assigned yet.")
        self._ref_info_label.setWordWrap(True)
        ref_layout.addWidget(self._ref_info_label)
        layout.addWidget(self._ref_group)

        layout.addStretch()

    # ── Public API ───────────────────────────────────────────────────────────

    def set_variable_names(self, var_names: list):
        """Updates the known descriptor variable names."""
        self._variable_names = list(var_names)
        self._update_ref_display()

    def update_from_conditions(self, conditions: dict):
        """Receives {condition_name: hex_color} and rebuilds the role table."""
        old_roles = self.get_roles()
        self._conditions = dict(conditions)

        cond_names = list(conditions.keys())
        self._role_table.setRowCount(len(cond_names))

        for row, name in enumerate(cond_names):
            # Condition name (read-only)
            item = QTableWidgetItem(name)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            self._role_table.setItem(row, 0, item)

            # Role combo
            combo = QComboBox()
            combo.addItems(self.ROLES)
            prev_role = old_roles.get(name, "Experimental")
            idx = combo.findText(prev_role)
            combo.setCurrentIndex(idx if idx >= 0 else 1)
            combo.currentTextChanged.connect(self._on_role_changed)
            self._role_table.setCellWidget(row, 1, combo)

        # Ensure at least one Reference Control exists
        roles = self.get_roles()
        if "Reference Control" not in roles.values() and cond_names:
            combo = self._role_table.cellWidget(0, 1)
            if combo:
                combo.setCurrentText("Reference Control")

        self._emit()

    def update_condition_variables(self, cond_vars: dict):
        """Receives {condition_name: {var_name: value}} from ExperimentalConditionsWidget."""
        self._condition_variables = dict(cond_vars)
        self._update_ref_display()

    # ── Save / Load ──────────────────────────────────────────────────────────

    def save_roles_config(self):
        """Returns a serialisable dict of current role assignments."""
        return self.get_roles()

    def load_roles_config(self, roles: dict):
        """Restores role assignments from a {condition_name: role} dict."""
        for row in range(self._role_table.rowCount()):
            item = self._role_table.item(row, 0)
            combo = self._role_table.cellWidget(row, 1)
            if item and combo:
                role = roles.get(item.text(), "Experimental")
                idx = combo.findText(role)
                if idx >= 0:
                    combo.setCurrentIndex(idx)
        self._emit()

    # ── Internal ─────────────────────────────────────────────────────────────

    def _on_role_changed(self):
        self._emit()

    def _emit(self):
        roles = self.get_roles()
        self.roles_updated.emit(roles)

        # Derive per-variable reference values from the Reference Control condition
        ref_condition = ""
        ref_vars = {}
        for cond, role in roles.items():
            if role == "Reference Control":
                ref_condition = cond
                if cond in self._condition_variables:
                    ref_vars = dict(self._condition_variables[cond])
                else:
                    # Fallback: split on "+"
                    parts = cond.split("+")
                    for i, var in enumerate(self._variable_names):
                        ref_vars[var] = parts[i] if i < len(parts) else ""
                break

        self.references_updated.emit(ref_vars, ref_condition)
        self._update_ref_display()

    def _update_ref_display(self):
        """Updates the reference level info label."""
        roles = self.get_roles()
        ref_cond = None
        for cond, role in roles.items():
            if role == "Reference Control":
                ref_cond = cond
                break

        if not ref_cond:
            self._ref_info_label.setText("No Reference Control assigned yet.")
            return

        if ref_cond in self._condition_variables:
            ref_vars = self._condition_variables[ref_cond]
        else:
            parts = ref_cond.split("+")
            ref_vars = {}
            for i, var in enumerate(self._variable_names):
                ref_vars[var] = parts[i] if i < len(parts) else ""

        lines = [f"<b>Reference Condition:</b> {ref_cond}"]
        for var, val in ref_vars.items():
            lines.append(f"&nbsp;&nbsp;<b>Reference {var}:</b> {val}")
        self._ref_info_label.setText("<br>".join(lines))

    # ── Public getters ───────────────────────────────────────────────────────

    def get_roles(self) -> dict:
        """Returns {condition_name: role_string}."""
        roles = {}
        for row in range(self._role_table.rowCount()):
            item = self._role_table.item(row, 0)
            combo = self._role_table.cellWidget(row, 1)
            if item and combo:
                roles[item.text()] = combo.currentText()
        return roles

    def get_references(self) -> tuple:
        """Returns ({var_name: ref_value}, ref_condition_name)."""
        roles = self.get_roles()
        for cond, role in roles.items():
            if role == "Reference Control":
                if cond in self._condition_variables:
                    return (dict(self._condition_variables[cond]), cond)
                else:
                    parts = cond.split("+")
                    ref_vars = {}
                    for i, var in enumerate(self._variable_names):
                        ref_vars[var] = parts[i] if i < len(parts) else ""
                    return (ref_vars, cond)
        return ({}, "")

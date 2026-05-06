from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QMessageBox, QTableWidget, QTableWidgetItem, QHeaderView, QSplitter,
    QComboBox, QAbstractItemView, QGridLayout, QInputDialog,
)
print("[DBG] plate_format: imports done", flush=True)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor

# =============================================================================
#  EXPERIMENTAL DESIGN WIDGET (PLATE FORMAT)
# =============================================================================

class ExperimentalPlateWidget(QWidget):
    """Sub-tab for assigning conditions to wells on a plate."""
    plate_layout_applied = Signal(object) # Emits the updated DataFrame
    layouts_changed      = Signal(list)   # Emits list of layout names when the library changes

    PLATE_FORMATS = {
        "6-well": (2, 3),
        "12-well": (3, 4),
        "24-well": (4, 6),
        "48-well": (6, 8),
        "96-well": (8, 12),
        "384-well": (16, 24)
    }

    def __init__(self, parent=None):
        print("[DBG] ExperimentalPlateWidget.__init__ START", flush=True)
        super().__init__(parent)
        self.df = None
        self.conditions = {}
        self.active_condition = None
        self.active_button = None
        self._roles = {}  # {condition_name: role_string}

        # Multi-layout library: {name: {"format": str, "map": {well_coord: condition}}}
        self._layouts     = {"Default": {"format": "96-well", "map": {}}}
        self._active_name = "Default"
        print("[DBG] plate_format: state initialized", flush=True)

        # --- Main Layout ---
        self.main_layout = QHBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        self.main_layout.addWidget(splitter)

        # --- Left Panel (Controls) ---
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setAlignment(Qt.AlignTop)

        self.info_label = QLabel("1. Define conditions in the previous tab.\n2. Click a condition button below.\n3. Click on wells, rows, or columns to paint the condition.")
        self.info_label.setWordWrap(True)

        # Layout library selector
        print("[DBG] plate_format: building layout library row", flush=True)
        layout_lib_row = QHBoxLayout()
        layout_lib_row.addWidget(QLabel("<b>Layout:</b>"))
        self.layout_combo = QComboBox()
        self.layout_combo.setMinimumWidth(140)
        self.layout_combo.addItem("Default")
        self.layout_combo.currentTextChanged.connect(self._on_layout_switched)
        layout_lib_row.addWidget(self.layout_combo, 1)

        self.add_layout_btn = QPushButton("+")
        self.add_layout_btn.setFixedWidth(28)
        self.add_layout_btn.setToolTip("Add a new plate layout")
        self.add_layout_btn.clicked.connect(self._on_add_layout)

        self.rename_layout_btn = QPushButton("Rename")
        self.rename_layout_btn.setToolTip("Rename the current layout")
        self.rename_layout_btn.clicked.connect(self._on_rename_layout)

        self.duplicate_layout_btn = QPushButton("Dup")
        self.duplicate_layout_btn.setToolTip("Duplicate the current layout")
        self.duplicate_layout_btn.clicked.connect(self._on_duplicate_layout)

        self.delete_layout_btn = QPushButton("-")
        self.delete_layout_btn.setFixedWidth(28)
        self.delete_layout_btn.setToolTip("Delete the current layout")
        self.delete_layout_btn.clicked.connect(self._on_delete_layout)

        for b in (self.add_layout_btn, self.rename_layout_btn,
                  self.duplicate_layout_btn, self.delete_layout_btn):
            layout_lib_row.addWidget(b)
        print("[DBG] plate_format: layout library row done", flush=True)

        # Plate format selection
        plate_format_layout = QHBoxLayout()
        plate_format_layout.addWidget(QLabel("<b>Plate Format:</b>"))
        self.plate_format_combo = QComboBox()
        self.plate_format_combo.addItems(self.PLATE_FORMATS.keys())
        self.plate_format_combo.currentTextChanged.connect(self._on_format_changed)
        plate_format_layout.addWidget(self.plate_format_combo)

        # Conditions container
        self.conditions_label = QLabel("<b>Click a Condition to Activate:</b>")
        self.conditions_widget = QWidget() # A container for the buttons
        self.conditions_layout = QVBoxLayout(self.conditions_widget)
        self.conditions_layout.setAlignment(Qt.AlignTop)

        # Action buttons
        action_layout = QGridLayout()
        self.clear_selection_button = QPushButton("Clear Plate Assignments")
        self.clear_selection_button.clicked.connect(self.clear_plate)
        self.apply_button = QPushButton("✔ Apply Plate Layout to Data")
        self.apply_button.clicked.connect(self.apply_layout_to_dataframe)
        self.apply_button.setEnabled(False)
        self.apply_button.setStyleSheet("background-color: #2a82da; color: white; font-weight: bold;")
        action_layout.addWidget(self.clear_selection_button, 0, 0)
        action_layout.addWidget(self.apply_button, 0, 1)


        left_layout.addWidget(self.info_label)
        left_layout.addSpacing(15)
        left_layout.addLayout(layout_lib_row)
        left_layout.addSpacing(8)
        left_layout.addLayout(plate_format_layout)
        left_layout.addSpacing(15)
        left_layout.addWidget(self.conditions_label)
        left_layout.addWidget(self.conditions_widget)
        left_layout.addStretch()
        left_layout.addLayout(action_layout)

        # --- Right Panel (Plate Grid) ---
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        self.plate_table = QTableWidget()
        self.plate_table.setSelectionMode(QAbstractItemView.NoSelection)
        right_layout.addWidget(self.plate_table)

        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setSizes([300, 700])

        # --- Connections for painting ---
        self.plate_table.itemClicked.connect(self.paint_well)
        self.plate_table.horizontalHeader().sectionClicked.connect(self.paint_column)
        self.plate_table.verticalHeader().sectionClicked.connect(self.paint_row)

        # Initialize
        print("[DBG] plate_format: calling setup_plate_grid", flush=True)
        self.setup_plate_grid()
        print("[DBG] ExperimentalPlateWidget.__init__ END", flush=True)

    def load_data(self, df):
        """Receives the DataFrame from the file/plate assignment step."""
        self.df = df.copy()
        self.info_label.setText(f"{len(df)} data rows loaded. Assign conditions to the plate wells.")
        self.apply_button.setEnabled(True)
        print("ExperimentalPlateWidget received DataFrame.")

    def update_roles(self, roles: dict):
        """Receives {condition_name: role_string} and refreshes button labels."""
        self._roles = dict(roles)
        self._refresh_button_labels()

    def _refresh_button_labels(self):
        """Updates button text to include role abbreviations."""
        _ABBREVS = {
            "Reference Control": "RC",
            "Experimental":      "EXP",
            "Positive Control":  "PC",
            "Potential Rescue":  "RES",
        }
        for i in range(self.conditions_layout.count()):
            widget = self.conditions_layout.itemAt(i).widget()
            if widget and isinstance(widget, QPushButton):
                cond_name = widget.property("condition_name")
                if cond_name:
                    role = self._roles.get(cond_name, "")
                    abbrev = _ABBREVS.get(role, "")
                    label = f"{cond_name}  [{abbrev}]" if abbrev else cond_name
                    widget.setText(label)

    def update_conditions(self, conditions):
        """Receives the dictionary of conditions and creates sorted buttons.
        Also repaints any already-painted wells on the grid to reflect new colors."""
        self.conditions = conditions  # Store for later use
        # Reset active state before deleting buttons to prevent runtime errors
        self.active_condition = None
        self.active_button = None

        # Clear existing buttons
        for i in reversed(range(self.conditions_layout.count())):
            widget = self.conditions_layout.itemAt(i).widget()
            if widget:
                widget.deleteLater()

        # --- Custom Sorting ---
        def sort_key(item):
            # Extracts the number from "Condition X" and returns it as an integer
            name = item[0]
            try:
                return int(''.join(filter(str.isdigit, name)))
            except (ValueError, TypeError):
                return float('inf') # Put non-numeric names at the end

        sorted_conditions = sorted(conditions.items(), key=sort_key)
        # --- End Custom Sorting ---

        _ABBREVS = {
            "Reference Control": "RC",
            "Experimental":      "EXP",
            "Positive Control":  "PC",
            "Potential Rescue":  "RES",
        }

        for name, color in sorted_conditions:
            role = self._roles.get(name, "")
            abbrev = _ABBREVS.get(role, "")
            label = f"{name}  [{abbrev}]" if abbrev else name

            btn = QPushButton(label)
            btn.setProperty("condition_name", name)  # store original name
            btn.setCheckable(True)
            btn.setStyleSheet(f"""
                QPushButton {{
                    background-color: {color};
                    color: {self._get_contrasting_text_color(color)};
                    border: 1px solid #555;
                    padding: 5px;
                }}
                QPushButton:checked {{
                    border: 3px solid #0078d7;
                }}
            """)
            btn.clicked.connect(lambda checked, b=btn, n=name, c=color: self.set_active_condition(n, c, b) if checked else self.deactivate_condition(b))
            self.conditions_layout.addWidget(btn)

        # Repaint existing wells on the grid to reflect updated colors
        self._repaint_grid_colors()

    def _repaint_grid_colors(self):
        """Updates the colors of all painted wells to match the current conditions dict."""
        for r in range(self.plate_table.rowCount()):
            for c in range(self.plate_table.columnCount()):
                item = self.plate_table.item(r, c)
                if item is None:
                    continue
                cond_name = item.data(Qt.UserRole)
                if not cond_name:
                    continue
                new_color = self.conditions.get(cond_name)
                if new_color:
                    item.setBackground(QColor(new_color))
                    item.setForeground(QColor(self._get_contrasting_text_color(new_color)))

    def _get_contrasting_text_color(self, hex_color):
        """Determines a contrasting text color (black or white) based on the background color."""
        try:
            # Convert hex to RGB
            hex_color = hex_color.lstrip('#')
            r = int(hex_color[0:2], 16)
            g = int(hex_color[2:4], 16)
            b = int(hex_color[4:6], 16)

            # Calculate brightness (YIQ color space)
            brightness = (r * 299 + g * 587 + b * 114) / 1000

            # Return black for light colors, white for dark colors
            return '#000000' if brightness > 128 else '#FFFFFF'
        except (ValueError, IndexError):
            return '#000000' # Default to black on error

    def set_active_condition(self, name, color, button_pressed):
        """Sets the currently active condition for 'painting' onto the grid."""
        # Deactivate any other active button
        if self.active_button and self.active_button is not button_pressed:
            self.active_button.setChecked(False)
        
        self.active_button = button_pressed
        self.active_condition = {'name': name, 'color': color}

    def deactivate_condition(self, button_pressed):
        """Deactivates the current condition if the active button is clicked again."""
        if self.active_button is button_pressed:
            self.active_condition = None
            self.active_button = None
            button_pressed.setChecked(False)

    def setup_plate_grid(self):
        """Configures the QTableWidget based on the selected plate format."""
        format_name = self.plate_format_combo.currentText()
        rows, cols = self.PLATE_FORMATS[format_name]

        self.plate_table.setRowCount(rows)
        self.plate_table.setColumnCount(cols)

        # Set headers
        row_headers = [chr(ord('A') + i) for i in range(rows)]
        col_headers = [str(i + 1) for i in range(cols)]
        self.plate_table.setVerticalHeaderLabels(row_headers)
        self.plate_table.setHorizontalHeaderLabels(col_headers)

        # Populate with empty items and resize
        for r in range(rows):
            self.plate_table.setRowHeight(r, 40)
            for c in range(cols):
                self.plate_table.setColumnWidth(c, 40)
                item = QTableWidgetItem()
                item.setTextAlignment(Qt.AlignCenter)
                self.plate_table.setItem(r, c, item)

        self.plate_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.plate_table.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)

    def _apply_active_to_cell(self, item):
        """Paints the active condition onto a cell unconditionally."""
        if item is None or not self.active_condition:
            return
        name  = self.active_condition['name']
        color = self.active_condition['color']
        item.setBackground(QColor(color))
        item.setForeground(QColor(self._get_contrasting_text_color(color)))
        item.setText(name)
        item.setData(Qt.UserRole, name)

    def _clear_cell(self, item):
        """Removes any condition assignment from a cell, restoring the default
        theme background by clearing the explicit background/foreground roles."""
        if item is None:
            return
        item.setData(Qt.BackgroundRole, None)
        item.setData(Qt.ForegroundRole, None)
        item.setText("")
        item.setData(Qt.UserRole, None)

    def paint_well(self, item):
        """Single-cell click handler.

        - If no active condition is selected: clicking a painted cell clears it
          (empty cell does nothing).
        - If an active condition is selected: clicking a cell that already has
          that condition clears it (toggle); otherwise paints it with the
          active condition.
        """
        if item is None:
            return
        current = item.data(Qt.UserRole)
        if not self.active_condition:
            if current:
                self._clear_cell(item)
            return
        if current == self.active_condition['name']:
            self._clear_cell(item)
        else:
            self._apply_active_to_cell(item)

    def paint_column(self, col_index):
        """Applies the active condition to all wells in a clicked column."""
        if not self.active_condition:
            return
        for row in range(self.plate_table.rowCount()):
            self._apply_active_to_cell(self.plate_table.item(row, col_index))

    def paint_row(self, row_index):
        """Applies the active condition to all wells in a clicked row."""
        if not self.active_condition:
            return
        for col in range(self.plate_table.columnCount()):
            self._apply_active_to_cell(self.plate_table.item(row_index, col))

    def apply_condition_to_selection(self, name, color):
        """Applies the selected condition's name and color to the selected wells."""
        selected_items = self.plate_table.selectedItems()
        if not selected_items:
            QMessageBox.information(self, "No Selection", "Please select one or more wells in the grid first.")
            return

        q_color = QColor(color)
        for item in selected_items:
            item.setBackground(q_color)
            item.setData(Qt.UserRole, name) # Store the condition name in the item
        
        self.plate_table.clearSelection() # Deselect after applying

    def load_plate_config(self, format_name: str, well_condition_map: dict):
        """
        Restores plate format and well assignments from saved config.

        Parameters
        ----------
        format_name : str
            One of the keys in PLATE_FORMATS (e.g., "96-well").
        well_condition_map : dict
            {well_coord: condition_name}, e.g., {"A1": "WT+DMSO", ...}
        """
        # Set plate format (triggers setup_plate_grid)
        idx = self.plate_format_combo.findText(format_name)
        if idx >= 0:
            self.plate_format_combo.setCurrentIndex(idx)

        row_headers = [self.plate_table.verticalHeaderItem(r).text()
                       for r in range(self.plate_table.rowCount())]

        for well, condition_name in well_condition_map.items():
            color = self.conditions.get(condition_name, "#cccccc")
            # Parse well coord: letter(s) + number
            # e.g., "A1" → row letter A, col 1
            letter_part = ''.join(c for c in well if c.isalpha())
            num_part    = ''.join(c for c in well if c.isdigit())
            if not letter_part or not num_part:
                continue
            try:
                row_idx = row_headers.index(letter_part)
                col_idx = int(num_part) - 1
            except (ValueError, IndexError):
                continue
            item = self.plate_table.item(row_idx, col_idx)
            if item:
                bg = QColor(color)
                item.setBackground(bg)
                item.setForeground(QColor(self._get_contrasting_text_color(color)))
                item.setText(condition_name)
                item.setData(Qt.UserRole, condition_name)

    def get_well_condition_map(self, name: str = None):
        """Returns {well_coord: condition_name}. Active layout reads from the grid;
        named inactive layout returns its stored map."""
        if name is not None and name != self._active_name:
            data = self._layouts.get(name)
            return dict(data["map"]) if data else {}

        well_to_condition = {}
        row_headers = [self.plate_table.verticalHeaderItem(r).text()
                       for r in range(self.plate_table.rowCount())]
        for r in range(self.plate_table.rowCount()):
            for c in range(self.plate_table.columnCount()):
                item = self.plate_table.item(r, c)
                if item:
                    condition_name = item.data(Qt.UserRole)
                    if condition_name:
                        well_name = f"{row_headers[r]}{c + 1}"
                        well_to_condition[well_name] = condition_name
        return well_to_condition

    def get_plate_format(self, name: str = None):
        """Returns (rows, cols) for the named layout, or active if ``name`` is None."""
        if name is not None and name != self._active_name:
            data = self._layouts.get(name)
            format_name = data["format"] if data else "96-well"
        else:
            format_name = self.plate_format_combo.currentText()
        return self.PLATE_FORMATS.get(format_name, (8, 12))

    def clear_plate(self):
        """Resets all well assignments in the table."""
        for r in range(self.plate_table.rowCount()):
            for c in range(self.plate_table.columnCount()):
                self._clear_cell(self.plate_table.item(r, c))

    def apply_layout_to_dataframe(self):
        """Maps the plate layout to the DataFrame and emits the result."""
        if self.df is None:
            QMessageBox.warning(self, "No Data", "No data has been loaded to apply the layout to.")
            return

        if 'plate_notation' not in self.df.columns:
            QMessageBox.critical(self, "Missing Column", "The required 'plate_notation' column was not found in the data. Please check the data binning step.")
            return

        well_to_condition_map = {}
        well_to_color_map = {}
        row_headers = [self.plate_table.verticalHeaderItem(r).text() for r in range(self.plate_table.rowCount())]
        
        for r in range(self.plate_table.rowCount()):
            for c in range(self.plate_table.columnCount()):
                item = self.plate_table.item(r, c)
                condition_name = item.data(Qt.UserRole)
                if condition_name:
                    well_name = f"{row_headers[r]}{c + 1}"
                    well_to_condition_map[well_name] = condition_name
                    # Also get the color from our stored conditions dictionary
                    well_to_color_map[well_name] = self.conditions.get(condition_name, "#ffffff") # Default to white if not found

        if not well_to_condition_map:
            QMessageBox.warning(self, "No Assignments", "No conditions have been assigned to any wells.")
            return
            
        # Create the new columns by mapping the 'plate_notation' column
        self.df['Condition'] = self.df['plate_notation'].map(well_to_condition_map)
        self.df['Color'] = self.df['plate_notation'].map(well_to_color_map)
        
        # Fill any unmapped wells with default values
        self.df['Condition'] = self.df['Condition'].fillna('Unassigned')
        self.df['Color'] = self.df['Color'].fillna('#ffffff')

        # Check how many rows were unassigned
        unassigned_count = (self.df['Condition'] == 'Unassigned').sum()
        
        msg = (f"Successfully applied conditions to the data.\n\n"
            f"{len(well_to_condition_map)} well assignments were mapped.\n"
            f"{len(self.df) - unassigned_count} data rows were assigned a condition.\n"
            f"{unassigned_count} rows had no matching well assignment and were marked as 'Unassigned'.")

        QMessageBox.information(self, "Success", msg)

        # Emit the updated dataframe for the next step in the pipeline
        self.plate_layout_applied.emit(self.df)
        print("Emitting DataFrame with 'Condition' and 'Color' columns.")

    # ── Multi-layout library ────────────────────────────────────────────────

    def _flush_active_to_layouts(self):
        if self._active_name is None:
            return
        self._layouts[self._active_name] = {
            "format": self.plate_format_combo.currentText(),
            "map":    self.get_well_condition_map(None),
        }

    def _paint_map_on_grid(self, well_condition_map: dict):
        row_headers = [self.plate_table.verticalHeaderItem(r).text()
                       for r in range(self.plate_table.rowCount())]
        for well, condition_name in well_condition_map.items():
            color = self.conditions.get(condition_name, "#cccccc")
            letter_part = ''.join(c for c in well if c.isalpha())
            num_part    = ''.join(c for c in well if c.isdigit())
            if not letter_part or not num_part:
                continue
            try:
                row_idx = row_headers.index(letter_part)
                col_idx = int(num_part) - 1
            except (ValueError, IndexError):
                continue
            item = self.plate_table.item(row_idx, col_idx)
            if item:
                item.setBackground(QColor(color))
                item.setForeground(QColor(self._get_contrasting_text_color(color)))
                item.setText(condition_name)
                item.setData(Qt.UserRole, condition_name)

    def _load_layout_into_grid(self, name: str):
        data = self._layouts.get(name)
        if not data:
            return
        self.plate_format_combo.blockSignals(True)
        idx = self.plate_format_combo.findText(data["format"])
        if idx >= 0:
            self.plate_format_combo.setCurrentIndex(idx)
        self.plate_format_combo.blockSignals(False)
        self.setup_plate_grid()
        self._paint_map_on_grid(data["map"])

    def _on_layout_switched(self, new_name: str):
        if not new_name or new_name == self._active_name:
            return
        if new_name not in self._layouts:
            return
        self._flush_active_to_layouts()
        self._active_name = new_name
        self._load_layout_into_grid(new_name)

    def _on_format_changed(self, new_format: str):
        if self._active_name:
            self._layouts[self._active_name] = {
                "format": new_format,
                "map":    {},
            }
        self.setup_plate_grid()

    def _refresh_layout_combo(self):
        self.layout_combo.blockSignals(True)
        self.layout_combo.clear()
        for name in self._layouts:
            self.layout_combo.addItem(name)
        if self._active_name:
            idx = self.layout_combo.findText(self._active_name)
            if idx >= 0:
                self.layout_combo.setCurrentIndex(idx)
        self.layout_combo.blockSignals(False)

    def _on_add_layout(self):
        name, ok = QInputDialog.getText(self, "New Plate Layout",
                                        "Name for the new layout:")
        name = name.strip() if name else ""
        if not ok or not name:
            return
        if name in self._layouts:
            QMessageBox.warning(self, "Name Exists",
                                f"A layout named '{name}' already exists.")
            return
        self._flush_active_to_layouts()
        self._layouts[name] = {"format": "96-well", "map": {}}
        self._active_name = name
        self._refresh_layout_combo()
        self._load_layout_into_grid(name)
        self.layouts_changed.emit(list(self._layouts.keys()))

    def _on_rename_layout(self):
        if not self._active_name:
            return
        current = self._active_name
        new_name, ok = QInputDialog.getText(self, "Rename Layout",
                                            "New name:", text=current)
        new_name = new_name.strip() if new_name else ""
        if not ok or not new_name or new_name == current:
            return
        if new_name in self._layouts:
            QMessageBox.warning(self, "Name Exists",
                                f"A layout named '{new_name}' already exists.")
            return
        self._flush_active_to_layouts()
        self._layouts = {
            (new_name if k == current else k): v
            for k, v in self._layouts.items()
        }
        self._active_name = new_name
        self._refresh_layout_combo()
        self.layouts_changed.emit(list(self._layouts.keys()))

    def _on_duplicate_layout(self):
        if not self._active_name:
            return
        base = self._active_name
        new_name, ok = QInputDialog.getText(self, "Duplicate Layout",
                                            "Name for the duplicate:",
                                            text=f"{base} (copy)")
        new_name = new_name.strip() if new_name else ""
        if not ok or not new_name:
            return
        if new_name in self._layouts:
            QMessageBox.warning(self, "Name Exists",
                                f"A layout named '{new_name}' already exists.")
            return
        self._flush_active_to_layouts()
        src = self._layouts[base]
        self._layouts[new_name] = {
            "format": src["format"],
            "map":    dict(src["map"]),
        }
        self._active_name = new_name
        self._refresh_layout_combo()
        self._load_layout_into_grid(new_name)
        self.layouts_changed.emit(list(self._layouts.keys()))

    def _on_delete_layout(self):
        if not self._active_name:
            return
        if len(self._layouts) == 1:
            QMessageBox.warning(self, "Cannot Delete",
                                "At least one plate layout must exist.")
            return
        current = self._active_name
        reply = QMessageBox.question(
            self, "Delete Layout",
            f"Delete layout '{current}'? This cannot be undone.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No
        )
        if reply != QMessageBox.Yes:
            return
        del self._layouts[current]
        self._active_name = next(iter(self._layouts))
        self._refresh_layout_combo()
        self._load_layout_into_grid(self._active_name)
        self.layouts_changed.emit(list(self._layouts.keys()))

    def get_layout_names(self) -> list:
        return list(self._layouts.keys())

    def get_all_layouts(self) -> dict:
        self._flush_active_to_layouts()
        return {
            name: {"format": data["format"], "map": dict(data["map"])}
            for name, data in self._layouts.items()
        }

    def load_layouts_config(self, layouts: dict, active: str = None):
        if not layouts:
            return
        self._layouts = {
            name: {
                "format": data.get("format", "96-well"),
                "map":    dict(data.get("map", {})),
            }
            for name, data in layouts.items()
        }
        chosen = active if (active and active in self._layouts) else next(iter(self._layouts))
        self._active_name = chosen
        self._refresh_layout_combo()
        self._load_layout_into_grid(chosen)
        self.layouts_changed.emit(list(self._layouts.keys()))

    def emit_layouts_state(self):
        print("[DBG] plate_format.emit_layouts_state() called", flush=True)
        self.layouts_changed.emit(list(self._layouts.keys()))
        print("[DBG] plate_format.emit_layouts_state() done", flush=True)
"""
Data Loader Widget
==================
Handles multi-directory data ingestion and preprocessing pipeline.

Workflow:
  1. User configures plate layout in ExperimentalPlateWidget (painted wells).
  2. User configures phases in ExperimentalSetupWidget (phase table).
  3. User adds one or more directories here — each directory = one plate replicate.
  4. Click "Load & Process Data" → runs the preprocessing pipeline → emits data_ready(full_df).
"""

import os
import threading
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
    QComboBox, QGroupBox, QProgressBar,
    QMessageBox, QFileDialog, QScrollArea, QSizePolicy, QDoubleSpinBox
)
print("[DBG] data_loader: imports done", flush=True)
from PySide6.QtCore import Qt, QEvent, Signal
from PySide6.QtGui import QPixmap, QImage


class DataLoaderWidget(QWidget):
    """
    Multi-directory data loading and preprocessing widget.

    Receives references to ExperimentalPlateWidget and ExperimentalSetupWidget
    at construction time so it can call get_well_condition_map() and
    get_phases_data() at the moment the user clicks "Load & Process Data".

    Emits
    -----
    data_ready : Signal(pd.DataFrame)
        Emits the fully preprocessed full_df once processing is complete.
        Columns: plate, Condition, loc_coord, location, time_sec, pxl_diff,
                 Phase, Group, Genotype, Drug
    """

    data_ready = Signal(object)   # pd.DataFrame

    # Private thread-safe signals
    _progress_signal = Signal(int, str)
    _done_signal     = Signal(bool, str)   # (success, error_message)

    def __init__(self, plate_widget=None, setup_widget=None, parent=None):
        """
        Parameters
        ----------
        plate_widget : ExperimentalPlateWidget
            Used to call get_well_condition_map() and get_plate_format().
        setup_widget : ExperimentalSetupWidget
            Used to call get_phases_data().
        """
        print("[DBG] DataLoaderWidget.__init__ START", flush=True)
        super().__init__(parent)
        self._plate_widget = plate_widget
        self._setup_widget = setup_widget
        # Each entry: {"path": str, "layout": str | None}
        self._directories      = []
        self._conditions       = {}   # {condition_name: hex_color} from Experimental Design
        self._variable_names   = ["Genotype", "Drug"]  # descriptor variable names
        self._available_layouts = []  # list of plate layout names (from plate widget)
        self._result_df        = None
        print("[DBG] data_loader: calling _build_ui()", flush=True)
        self._build_ui()
        print("[DBG] data_loader: _build_ui() done", flush=True)
        self._progress_signal.connect(self._on_progress)
        self._done_signal.connect(self._on_done)
        print("[DBG] DataLoaderWidget.__init__ END", flush=True)

    # ── Public API ────────────────────────────────────────────────────────────

    def set_plate_widget(self, widget):
        """Set the ExperimentalPlateWidget reference (can be called after init)."""
        self._plate_widget = widget

    def set_setup_widget(self, widget):
        """Set the ExperimentalSetupWidget reference (can be called after init)."""
        self._setup_widget = widget

    def set_variable_names(self, names: list):
        """Stores the descriptor variable names for condition splitting."""
        self._variable_names = list(names)

    # ── Build UI ──────────────────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # ── Directory list group ──────────────────────────────────────────────
        print("[DBG] data_loader: building dir_group (QTableWidget)", flush=True)
        dir_group = QGroupBox("Input Directories  (each directory = one plate replicate, in order)")
        dir_group_layout = QVBoxLayout(dir_group)

        self._dir_table = QTableWidget(0, 2)
        self._dir_table.setHorizontalHeaderLabels(["Directory", "Plate Layout"])
        self._dir_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        # Plate Layout column is sized to ~1/4 of the viewport width and kept
        # proportional via the eventFilter installed below.
        self._dir_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Fixed)
        self._dir_table.setColumnWidth(1, 240)
        self._dir_table.installEventFilter(self)
        self._dir_table.verticalHeader().setVisible(False)
        self._dir_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._dir_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self._dir_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._dir_table.setAlternatingRowColors(True)
        self._dir_table.setMaximumHeight(140)
        self._dir_table.setToolTip(
            "Each directory is treated as a separate plate replicate.\n"
            "The topmost row = Plate 1, next = Plate 2, etc.\n"
            "Pick a Plate Layout for each plate from the dropdown in the right column."
        )
        dir_group_layout.addWidget(self._dir_table)
        print("[DBG] data_loader: dir_group built", flush=True)

        btn_row = QHBoxLayout()
        self._add_btn    = QPushButton("Add Directory…")
        self._remove_btn = QPushButton("Remove Selected")
        self._up_btn     = QPushButton("▲ Move Up")
        self._down_btn   = QPushButton("▼ Move Down")

        self._add_btn.clicked.connect(self._on_add)
        self._remove_btn.clicked.connect(self._on_remove)
        self._up_btn.clicked.connect(self._on_move_up)
        self._down_btn.clicked.connect(self._on_move_down)

        for b in (self._add_btn, self._remove_btn, self._up_btn, self._down_btn):
            btn_row.addWidget(b)
        btn_row.addStretch()
        dir_group_layout.addLayout(btn_row)
        layout.addWidget(dir_group)

        # ── Processing group ──────────────────────────────────────────────────
        proc_group = QGroupBox("Processing")
        proc_layout = QVBoxLayout(proc_group)

        run_row = QHBoxLayout()
        self._run_btn = QPushButton("Load & Process Data")
        self._run_btn.setEnabled(False)
        self._run_btn.setMinimumHeight(36)
        self._run_btn.setStyleSheet("QPushButton { font-weight: bold; }")
        self._run_btn.clicked.connect(self._on_run)
        run_row.addWidget(self._run_btn)
        run_row.addStretch()
        proc_layout.addLayout(run_row)

        self._progress_bar = QProgressBar()
        self._progress_bar.setRange(0, 100)
        self._progress_bar.setValue(0)
        self._progress_bar.setTextVisible(True)
        proc_layout.addWidget(self._progress_bar)

        self._status_label = QLabel("Ready. Add directories and configure plate layout / phases first.")
        self._status_label.setWordWrap(True)
        proc_layout.addWidget(self._status_label)
        layout.addWidget(proc_group)

        # ── Data preview figures ──────────────────────────────────────────────
        fig_group = QGroupBox("Data Preview")
        fig_layout = QVBoxLayout(fig_group)

        fig_btn_row = QHBoxLayout()
        self._update_fig_btn = QPushButton("Update Figures")
        self._update_fig_btn.setEnabled(False)
        self._update_fig_btn.setMinimumHeight(32)
        self._update_fig_btn.setStyleSheet("QPushButton { font-weight: bold; }")
        self._update_fig_btn.clicked.connect(self._update_figures)
        self._save_fig_btn = QPushButton("Save Figures…")
        self._save_fig_btn.setEnabled(False)
        self._save_fig_btn.clicked.connect(self._save_figures)
        fig_btn_row.addWidget(self._update_fig_btn)
        fig_btn_row.addWidget(self._save_fig_btn)
        fig_btn_row.addSpacing(20)
        fig_btn_row.addWidget(QLabel("Line width:"))
        self._line_width_spin = QDoubleSpinBox()
        self._line_width_spin.setRange(0.1, 5.0)
        self._line_width_spin.setSingleStep(0.1)
        self._line_width_spin.setValue(0.4)
        self._line_width_spin.setFixedWidth(70)
        fig_btn_row.addWidget(self._line_width_spin)
        fig_btn_row.addStretch()
        fig_layout.addLayout(fig_btn_row)

        # Navigation bar (prev / page label / next)
        nav_bar = QHBoxLayout()
        self._fig_prev_btn = QPushButton("◀")
        self._fig_prev_btn.setFixedWidth(40)
        self._fig_prev_btn.setEnabled(False)
        self._fig_prev_btn.clicked.connect(self._fig_prev)
        self._fig_page_label = QLabel("")
        self._fig_page_label.setAlignment(Qt.AlignCenter)
        self._fig_page_label.setMinimumWidth(120)
        self._fig_next_btn = QPushButton("▶")
        self._fig_next_btn.setFixedWidth(40)
        self._fig_next_btn.setEnabled(False)
        self._fig_next_btn.clicked.connect(self._fig_next)
        nav_bar.addWidget(self._fig_prev_btn)
        nav_bar.addWidget(self._fig_page_label)
        nav_bar.addWidget(self._fig_next_btn)
        nav_bar.addStretch()
        fig_layout.addLayout(nav_bar)

        # Single figure display area
        self._fig_pixmaps = []     # list of QPixmap
        self._fig_titles  = []     # list of str
        self._fig_index   = 0
        self._fig_scroll = QScrollArea()
        self._fig_scroll.setWidgetResizable(True)
        self._fig_scroll.setAlignment(Qt.AlignCenter)
        self._fig_image_label = QLabel("No figures loaded.")
        self._fig_image_label.setAlignment(Qt.AlignCenter)
        self._fig_image_label.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored)
        self._fig_scroll.setWidget(self._fig_image_label)
        fig_layout.addWidget(self._fig_scroll, stretch=1)

        layout.addWidget(fig_group, stretch=1)

    # ── Condition colors ────────────────────────────────────────────────────────

    def update_conditions(self, conditions: dict):
        """Receives {condition_name: hex_color} from ExperimentalConditionsWidget."""
        self._conditions = dict(conditions)

    # ── Figure generation ──────────────────────────────────────────────────────

    def _update_figures(self):
        """Generates both preview figures using the current data and condition colors."""
        if self._result_df is None or self._result_df.empty:
            QMessageBox.information(self, "No Data", "Load data first before generating figures.")
            return

        df = self._result_df
        colors = self._conditions

        self._fig_pixmaps = []
        self._fig_titles = []

        # Build figure 1: all individual data points
        fig1 = self._make_raw_figure(df, colors)
        self._fig_pixmaps.append(self._fig_to_pixmap(fig1))
        self._fig_titles.append("All Individual Traces")

        # Build figure 2: mean per animal with 95% CI
        fig2 = self._make_mean_ci_figure(df, colors)
        self._fig_pixmaps.append(self._fig_to_pixmap(fig2))
        self._fig_titles.append("Mean Movement per Condition (95% CI)")

        self._fig_index = 0
        self._show_current_figure()
        self._save_fig_btn.setEnabled(True)

    def _show_current_figure(self):
        """Displays the figure at _fig_index and updates nav buttons."""
        n = len(self._fig_pixmaps)
        if n == 0:
            return
        idx = self._fig_index
        pixmap = self._fig_pixmaps[idx]
        available = self._fig_scroll.viewport().size()
        scaled = pixmap.scaled(available, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self._fig_image_label.setPixmap(scaled)
        self._fig_page_label.setText(f"{self._fig_titles[idx]}  ({idx + 1}/{n})")
        self._fig_prev_btn.setEnabled(idx > 0)
        self._fig_next_btn.setEnabled(idx < n - 1)

    def _fig_prev(self):
        if self._fig_index > 0:
            self._fig_index -= 1
            self._show_current_figure()

    def _fig_next(self):
        if self._fig_index < len(self._fig_pixmaps) - 1:
            self._fig_index += 1
            self._show_current_figure()

    def _make_raw_figure(self, df, colors):
        """All individual data points, one line per animal, colored by condition."""
        lw = self._line_width_spin.value()
        fig, ax = plt.subplots(figsize=(7, 4.5), dpi=120)

        conditions = df['Condition'].unique()
        for cond in sorted(conditions):
            color = colors.get(cond, '#888888')
            sub = df[df['Condition'] == cond]
            for (plate, loc), animal_df in sub.groupby(['plate', 'location']):
                animal_df = animal_df.sort_values('time_sec')
                ax.plot(animal_df['time_sec'], animal_df['pxl_diff'],
                        color=color, alpha=0.15, linewidth=lw, rasterized=True)
            # Invisible line for legend entry
            ax.plot([], [], color=color, linewidth=4, label=cond)

        ax.set_xlabel('Time (s)', fontsize=10)
        ax.set_ylabel('Movement (px diff)', fontsize=10)
        ax.set_title('All Individual Traces', fontsize=11, fontweight='bold')
        leg = ax.legend(fontsize=8, loc='upper right', framealpha=0.8)
        for line in leg.get_lines():
            line.set_linewidth(4)
        ax.tick_params(labelsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.tight_layout()
        return fig

    def _make_mean_ci_figure(self, df, colors):
        """Mean movement per condition."""
        lw = self._line_width_spin.value()
        fig, ax = plt.subplots(figsize=(7, 4.5), dpi=120)

        conditions = df['Condition'].unique()
        for cond in sorted(conditions):
            color = colors.get(cond, '#888888')
            sub = df[df['Condition'] == cond]

            animal_means = (
                sub.groupby(['plate', 'location', 'time_sec'])['pxl_diff']
                .sum()
                .reset_index()
            )
            stats = animal_means.groupby('time_sec')['pxl_diff'].agg(['mean', 'sem', 'count']).reset_index()
            stats.sort_values('time_sec', inplace=True)

            ax.plot(stats['time_sec'], stats['mean'],
                    color=color, linewidth=lw, label=cond)

        ax.set_xlabel('Time (s)', fontsize=10)
        ax.set_ylabel('Mean Movement (px diff)', fontsize=10)
        ax.set_title('Mean Movement per Condition', fontsize=11, fontweight='bold')
        leg = ax.legend(fontsize=8, loc='upper right', framealpha=0.8)
        for line in leg.get_lines():
            line.set_linewidth(4)
        ax.tick_params(labelsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fig.tight_layout()
        return fig

    @staticmethod
    def _fig_to_pixmap(fig):
        """Renders a matplotlib Figure to a QPixmap."""
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        buf = canvas.buffer_rgba()
        w, h = canvas.get_width_height()
        img = QImage(buf, w, h, QImage.Format_RGBA8888)
        pixmap = QPixmap.fromImage(img)
        plt.close(fig)
        return pixmap

    def _save_figures(self):
        """Saves both preview figures to a user-selected directory."""
        dest_dir = QFileDialog.getExistingDirectory(self, "Select Folder to Save Figures")
        if not dest_dir:
            return

        # Regenerate fresh figures for saving at higher quality
        df = self._result_df
        colors = self._conditions
        saved = 0
        for name, make_fn in [("all_traces", self._make_raw_figure),
                               ("mean_95ci", self._make_mean_ci_figure)]:
            try:
                fig = make_fn(df, colors)
                path = os.path.join(dest_dir, f"{name}.png")
                fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
                plt.close(fig)
                saved += 1
            except Exception as e:
                QMessageBox.warning(self, "Save Error", f"Could not save {name}: {e}")

        if saved:
            QMessageBox.information(self, "Saved", f"Saved {saved} figure(s) to:\n{dest_dir}")

    # ── Directory list management ─────────────────────────────────────────────

    def update_layouts(self, names: list):
        """Slot: called when the plate widget's layout library changes."""
        print(f"[DBG] data_loader.update_layouts({names}) called", flush=True)
        self._available_layouts = list(names)
        for entry in self._directories:
            if entry["layout"] and entry["layout"] not in self._available_layouts:
                entry["layout"] = None
        self._rebuild_table()
        print("[DBG] data_loader.update_layouts done", flush=True)

    def _on_add(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Input Directory")
        if not directory:
            return
        self._directories.append({"path": directory, "layout": None})
        self._rebuild_table()
        self._run_btn.setEnabled(True)

    def _on_remove(self):
        row = self._dir_table.currentRow()
        if row < 0:
            return
        del self._directories[row]
        self._rebuild_table()
        self._run_btn.setEnabled(bool(self._directories))

    def _on_move_up(self):
        row = self._dir_table.currentRow()
        if row <= 0:
            return
        self._directories[row], self._directories[row - 1] = (
            self._directories[row - 1], self._directories[row]
        )
        self._rebuild_table()
        self._dir_table.setCurrentCell(row - 1, 0)

    def _on_move_down(self):
        row = self._dir_table.currentRow()
        if row < 0 or row >= len(self._directories) - 1:
            return
        self._directories[row], self._directories[row + 1] = (
            self._directories[row + 1], self._directories[row]
        )
        self._rebuild_table()
        self._dir_table.setCurrentCell(row + 1, 0)

    def _rebuild_table(self):
        """Rebuild the directory table rows from self._directories."""
        self._dir_table.setRowCount(len(self._directories))
        for i, entry in enumerate(self._directories):
            label_item = QTableWidgetItem(f"Plate {i + 1}:  {entry['path']}")
            label_item.setToolTip(entry["path"])
            self._dir_table.setItem(i, 0, label_item)

            combo = QComboBox()
            combo.addItem("- pick a layout -", None)
            for name in self._available_layouts:
                combo.addItem(name, name)
            if entry["layout"]:
                idx = combo.findData(entry["layout"])
                if idx >= 0:
                    combo.setCurrentIndex(idx)
            # Let the dropdown breathe so it doesn't collapse to the longest
            # layout name; the column itself is sized to ~1/4 of the table.
            combo.setMinimumWidth(200)
            combo.currentIndexChanged.connect(
                lambda _idx, row=i: self._on_layout_combo_changed(row)
            )
            self._dir_table.setCellWidget(i, 1, combo)

    def _on_layout_combo_changed(self, row: int):
        if row < 0 or row >= len(self._directories):
            return
        combo = self._dir_table.cellWidget(row, 1)
        if combo is None:
            return
        self._directories[row]["layout"] = combo.currentData()
        # Clear the missing-layout error styling now that a layout is picked.
        if combo.currentData() is not None:
            self._set_combo_error(combo, False)

    # ── Plate-layout column sizing + error highlighting ─────────────────────────

    def eventFilter(self, obj, event):
        """Keep the Plate Layout column at ~1/4 of the table viewport width."""
        if obj is self._dir_table and event.type() == QEvent.Resize:
            viewport_w = self._dir_table.viewport().width()
            if viewport_w > 0:
                self._dir_table.setColumnWidth(1, max(200, viewport_w // 4))
        return super().eventFilter(obj, event)

    @staticmethod
    def _set_combo_error(combo: QComboBox, on: bool):
        """Toggle a red-border highlight on a layout combo to flag a missing pick."""
        if on:
            combo.setStyleSheet(
                "QComboBox {"
                " border: 2px solid #e53935;"
                " background-color: rgba(229, 57, 53, 0.18);"
                "}"
            )
        else:
            combo.setStyleSheet("")

    # ── Run pipeline ──────────────────────────────────────────────────────────

    def _on_run(self):
        if self._plate_widget is None or self._setup_widget is None:
            QMessageBox.critical(self, "Configuration Error",
                                 "Plate widget or setup widget not configured.")
            return

        if not self._directories:
            QMessageBox.warning(self, "No Directories",
                                "Add at least one input directory first.")
            return

        missing_rows = [i for i, d in enumerate(self._directories) if not d["layout"]]
        if missing_rows:
            # Highlight every offending combo in red so the user can spot them
            # at a glance instead of cross-referencing plate numbers from the dialog.
            for row in missing_rows:
                combo = self._dir_table.cellWidget(row, 1)
                if combo is not None:
                    self._set_combo_error(combo, True)
            QMessageBox.warning(
                self, "Plate Layout Required",
                "Please choose a plate layout for plate(s): "
                f"{', '.join(str(r + 1) for r in missing_rows)}.\n"
                "Every directory must have a plate layout assigned before processing."
            )
            return

        per_directory_layouts = []
        for i, entry in enumerate(self._directories, start=1):
            name = entry["layout"]
            wc_map     = self._plate_widget.get_well_condition_map(name)
            rows, cols = self._plate_widget.get_plate_format(name)
            if not wc_map:
                QMessageBox.warning(
                    self, "Empty Plate Layout",
                    f"Plate layout '{name}' (assigned to Plate {i}) has no wells painted.\n"
                    "Please go to the Plate Format tab and paint wells first."
                )
                return
            per_directory_layouts.append({
                "name": name,
                "rows": rows,
                "cols": cols,
                "map":  wc_map,
            })

        phases_data = self._setup_widget.get_phases_data()
        if not phases_data:
            QMessageBox.warning(
                self, "No Phase Data",
                "No phases have been defined.\n"
                "Please go to the Experimental Design tab and configure the phases first."
            )
            return

        time_frame = pd.DataFrame([
            {'End': p['End_s'], 'Phase': p['Phase'], 'Group': p['Group']}
            for p in phases_data
        ])

        self._run_btn.setEnabled(False)
        self._progress_bar.setValue(0)
        self._status_label.setText("Processing…")

        paths = [e["path"] for e in self._directories]
        t = threading.Thread(
            target=self._worker,
            args=(paths, per_directory_layouts, time_frame),
            daemon=True
        )
        t.start()

    def _worker(self, directories, per_directory_layouts, time_frame):
        """Background thread: runs the preprocessing pipeline."""
        try:
            self._progress_signal.emit(0, "Scanning directories…")
            plate_file_lists = []
            total_files = 0
            for directory in directories:
                base_path = Path(directory)
                folders = [base_path] + [f for f in base_path.rglob('*') if f.is_dir()]
                folder_files = []
                for folder in folders:
                    files = sorted(
                        f for f in folder.glob('*')
                        if f.suffix.lower() in ('.csv', '.xls', '.xlsx')
                    )
                    if files:
                        folder_files.append((folder, files))
                        total_files += len(files)
                plate_file_lists.append(folder_files)

            if total_files == 0:
                self._done_signal.emit(False, "No .csv / .xls / .xlsx files were found in the selected directories.")
                return

            all_dfs    = []
            n_dirs     = len(directories)
            files_done = 0

            for plate_idx, (directory, folder_files, layout_info) in enumerate(
                zip(directories, plate_file_lists, per_directory_layouts), start=1
            ):
                rows_p, cols_p = layout_info["rows"], layout_info["cols"]
                row_letters = [chr(ord('A') + i) for i in range(rows_p)]
                n_wells = rows_p * cols_p
                coords = [f"{row_letters[i // cols_p]}{(i % cols_p) + 1}"
                          for i in range(n_wells)]
                loc_to_coord     = {i + 1: coords[i] for i in range(n_wells)}
                loc_to_condition = {i + 1: layout_info["map"].get(coords[i])
                                    for i in range(n_wells)}
                layout_name = layout_info["name"]

                plate_dfs = []

                for folder, files in folder_files:
                    file_dfs = []
                    for f in files:
                        files_done += 1
                        pct = int(files_done / total_files * 78)
                        self._progress_signal.emit(
                            pct,
                            f"Plate {plate_idx}/{n_dirs} [{layout_name}] - "
                            f"file {files_done}/{total_files}: {f.name}"
                        )

                        try:
                            if f.suffix.lower() == '.xlsx':
                                df_f = pd.read_excel(
                                    f, dtype={'data1': str}, keep_default_na=False
                                )
                            else:
                                df_f = pd.read_csv(
                                    f, sep='\t', low_memory=False,
                                    dtype={'data1': str}, keep_default_na=False
                                )
                            file_dfs.append(df_f)
                        except Exception as e:
                            print(f"Warning: could not read {f}: {e}")

                    if not file_dfs:
                        continue

                    temp_df = pd.concat(file_dfs, ignore_index=True)

                    needed = [c for c in ('time', 'location', 'data1') if c in temp_df.columns]
                    if not needed or 'time' not in needed:
                        print(f"Warning: missing required columns in {directory}. Skipping.")
                        continue

                    temp_df = temp_df[temp_df['time'] > 0][['time', 'location', 'data1']].copy()
                    temp_df['loc_id']    = temp_df['location'].str.extract(r'(\d+)').astype(int)
                    temp_df['Condition'] = temp_df['loc_id'].map(loc_to_condition)
                    temp_df['loc_coord'] = temp_df['loc_id'].map(loc_to_coord)
                    temp_df['time_sec']  = (temp_df['time'] / 1_000_000).astype(int)
                    plate_dfs.append(temp_df)

                if plate_dfs:
                    plate_combined = pd.concat(plate_dfs, ignore_index=True)
                    plate_combined['plate']        = plate_idx
                    plate_combined['plate_layout'] = layout_name
                    all_dfs.append(plate_combined)

            if not all_dfs:
                self._done_signal.emit(False, "Files were found but none contained the required columns (time, location, data1).")
                return

            self._progress_signal.emit(80, f"Aggregating {len(all_dfs)} plate(s)…")

            merged_df = pd.concat(all_dfs, ignore_index=True)
            merged_df['pxl_diff'] = pd.to_numeric(merged_df['data1'], errors='coerce')

            binned_df = (
                merged_df
                .groupby(['plate', 'plate_layout', 'Condition', 'loc_coord', 'location', 'time_sec'])['pxl_diff']
                .sum()
                .reset_index()
                .sort_values(['plate', 'time_sec', 'location'])
            )

            self._progress_signal.emit(90, "Assigning phases (temporal join)…")

            # Causal pivot: merge_asof to assign Phase and Group
            binned_sorted = binned_df.sort_values('time_sec').copy()
            binned_sorted['time_sec'] = binned_sorted['time_sec'].astype(np.int64)

            tf_sorted = time_frame[['End', 'Phase', 'Group']].sort_values('End').copy()
            tf_sorted['End'] = tf_sorted['End'].astype(np.int64)

            binned_merged = pd.merge_asof(
                binned_sorted,
                tf_sorted,
                left_on='time_sec',
                right_on='End',
                direction='forward'
            ).drop(columns='End')

            # Remove empty wells
            full_df = binned_merged[
                binned_merged['Condition'].notna() &
                (binned_merged['Condition'] != '')
            ].sort_values(['plate', 'time_sec', 'location']).reset_index(drop=True)

            # Split Condition → individual descriptor variables (dynamic)
            split = full_df['Condition'].str.split('+', expand=True)
            for i, var_name in enumerate(self._variable_names):
                full_df[var_name] = split[i] if i in split.columns else ''

            self._progress_signal.emit(100, f"Done. {len(full_df):,} rows loaded across {len(directories)} plate(s).")
            self._done_signal.emit(True, "")
            # Store result for _on_done to emit
            self._result_df = full_df

        except Exception as e:
            import traceback
            self._done_signal.emit(False, f"{e}\n\n{traceback.format_exc()}")

    # ── Thread-safe slots ─────────────────────────────────────────────────────

    def _on_progress(self, percent: int, message: str):
        self._progress_bar.setValue(percent)
        self._status_label.setText(message)

    def _on_done(self, success: bool, error_msg: str):
        self._run_btn.setEnabled(True)
        if success:
            df = self._result_df
            n_plates     = df['plate'].nunique()
            n_conditions = df['Condition'].nunique()
            n_phases     = df['Phase'].nunique() if 'Phase' in df.columns else '?'
            self._status_label.setText(
                f"✓ Loaded {len(df):,} rows  |  "
                f"{n_plates} plate(s)  |  "
                f"{n_conditions} condition(s)  |  "
                f"{n_phases} phase(s)"
            )
            self._update_fig_btn.setEnabled(True)
            self.data_ready.emit(df)
            # Auto-generate figures on first load
            self._update_figures()
        else:
            self._progress_bar.setValue(0)
            self._status_label.setText("Error during processing. See details below.")
            QMessageBox.critical(self, "Processing Error", error_msg)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self._fig_pixmaps:
            self._show_current_figure()

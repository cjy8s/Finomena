"""
Finomena: Zebrafish Time-Series Behavioral Analysis Suite
===============================================
Main entry point for the rebuilt app.

Tab layout:
  1. Experimental Design  — plate layout + phases + conditions + reference controls
  2. Data Loading         — multi-directory input, Cell-4 preprocessing pipeline
  3. Analysis in R        — optional family selection + BAM analysis
"""

import os
import signal

# Workarounds for VS Code integrated-terminal shell-integration escape
# sequences being misinterpreted as Ctrl+C on Windows:
#   1. Disable Intel Fortran's console handler (must be set before
#      libifcoremd.dll loads via numpy/scipy), or the process aborts
#      with forrtl error 200.
#   2. Ignore SIGINT in the Python main thread, or imports get
#      interrupted with KeyboardInterrupt mid-decode.
# Closing the GUI window remains the normal way to exit.
os.environ["FOR_DISABLE_CONSOLE_CTRL_HANDLER"] = "1"
signal.signal(signal.SIGINT, signal.SIG_IGN)

import json
import sys

print("[DBG] app.py: script start", flush=True)

# ── Path setup ─────────────────────────────────────────────────────────────────
# Must be self-contained (cannot import paths.py yet — it lives inside utils/).
if getattr(sys, 'frozen', False):
    _APP_DIR = getattr(sys, '_MEIPASS', os.path.dirname(sys.executable))
else:
    _APP_DIR = os.path.dirname(os.path.abspath(__file__))

_UTILS_DIR = os.path.join(_APP_DIR, "Finomena", "utils")
if _UTILS_DIR not in sys.path:
    sys.path.insert(0, _UTILS_DIR)

from PySide6.QtWidgets import (
    QApplication, QFileDialog, QHBoxLayout, QLabel, QMainWindow,
    QMessageBox, QPushButton, QTabWidget, QVBoxLayout, QWidget
)
from PySide6.QtCore import Qt

# ── Widget imports ─────────────────────────────────────────────────────────────
from experiment_design import (
    ExperimentalSetupWidget,
    ExperimentalConditionsWidget,
    MetadataAssignmentWidget,
)
print("[DBG] app.py: about to import plate_format", flush=True)
from plate_format     import ExperimentalPlateWidget
print("[DBG] app.py: plate_format imported", flush=True)
from dataframe_viewer import DataFrameViewerWidget
print("[DBG] app.py: about to import data_loader", flush=True)
from data_loader      import DataLoaderWidget
print("[DBG] app.py: data_loader imported", flush=True)
from bam_widget                import BamWidget
from family_selection_widget   import FamilySelectionWidget
from contrast_selection_widget import ContrastSelectionWidget
print("[DBG] app.py: all widget imports done", flush=True)


# =============================================================================
#  MAIN WINDOW
# =============================================================================

class MainWindow(QMainWindow):
    def __init__(self):
        print("[DBG] MainWindow.__init__ START", flush=True)
        super().__init__()
        self.setWindowTitle("Finomena: Zebrafish Time-Series Behavioral Analysis Suite")
        self.setGeometry(100, 100, 1280, 860)
        self.setStyleSheet(_SPINBOX_STYLESHEET)

        self.main_tabs = QTabWidget()
        self.setCentralWidget(self.main_tabs)

        # ── Tab 1: Experimental Design ─────────────────────────────────────────
        exp_design_tab = QWidget()
        exp_design_layout = QVBoxLayout(exp_design_tab)

        # Save / Load config bar
        config_bar = QHBoxLayout()
        config_bar.addWidget(QLabel("<b>Experimental Design Config:</b>"))
        self._save_config_btn = QPushButton("💾  Save Config…")
        self._save_config_btn.setToolTip("Save phases, conditions, and plate layout to a JSON file")
        self._save_config_btn.clicked.connect(self._save_config)
        self._load_config_btn = QPushButton("📂  Load Config…")
        self._load_config_btn.setToolTip("Load a previously saved config JSON file")
        self._load_config_btn.clicked.connect(self._load_config)
        config_bar.addWidget(self._save_config_btn)
        config_bar.addWidget(self._load_config_btn)
        config_bar.addStretch()
        exp_design_layout.addLayout(config_bar)

        exp_sub_tabs = QTabWidget()
        exp_design_layout.addWidget(exp_sub_tabs)

        print("[DBG] building ExperimentalSetupWidget", flush=True)
        self.experiment_format    = ExperimentalSetupWidget()
        print("[DBG] building ExperimentalConditionsWidget", flush=True)
        self.conditions_format    = ExperimentalConditionsWidget()
        print("[DBG] building MetadataAssignmentWidget", flush=True)
        self.metadata_assignment  = MetadataAssignmentWidget()
        print("[DBG] building ExperimentalPlateWidget", flush=True)
        self.plate_format_widget  = ExperimentalPlateWidget()
        print("[DBG] ExperimentalPlateWidget built", flush=True)

        exp_sub_tabs.addTab(self.conditions_format,   "1. Experimental Conditions")
        exp_sub_tabs.addTab(self.metadata_assignment,  "2. Metadata Assignment")
        exp_sub_tabs.addTab(self.plate_format_widget, "3. Plate Format")
        exp_sub_tabs.addTab(self.experiment_format,   "4. Experimental Setup")

        self.main_tabs.addTab(exp_design_tab, "Experimental Design")

        # ── Tab 2: Data Loading ────────────────────────────────────────────────
        data_loading_tab = QWidget()
        dl_layout = QVBoxLayout(data_loading_tab)
        dl_sub_tabs = QTabWidget()
        dl_layout.addWidget(dl_sub_tabs)

        print("[DBG] building DataLoaderWidget", flush=True)
        self.data_loader = DataLoaderWidget(
            plate_widget=self.plate_format_widget,
            setup_widget=self.experiment_format,
        )
        print("[DBG] DataLoaderWidget built", flush=True)

        # Output directory chooser
        self._output_dir = None
        output_dir_tab = QWidget()
        od_layout = QVBoxLayout(output_dir_tab)
        od_layout.setAlignment(Qt.AlignTop)
        od_layout.addWidget(QLabel(
            "<b>Set the output directory for all analysis results.</b><br>"
            "The pre-processed CSV, BAM results, and figures will be saved here."
        ))
        od_row = QHBoxLayout()
        self._output_dir_btn = QPushButton("Choose Output Directory…")
        self._output_dir_btn.setMinimumHeight(36)
        self._output_dir_btn.setStyleSheet("QPushButton { font-weight: bold; }")
        self._output_dir_btn.clicked.connect(self._choose_output_dir)
        self._output_dir_label = QLabel("No output directory selected.")
        self._output_dir_label.setWordWrap(True)
        od_row.addWidget(self._output_dir_btn)
        od_row.addWidget(self._output_dir_label, stretch=1)
        od_layout.addLayout(od_row)

        self.sanity_check_viewer = DataFrameViewerWidget()

        dl_sub_tabs.addTab(self.data_loader,         "1. Load Data")
        dl_sub_tabs.addTab(output_dir_tab,           "2. Set Output Directory")
        dl_sub_tabs.addTab(self.sanity_check_viewer, "3. Sanity Check")

        self.main_tabs.addTab(data_loading_tab, "Data Loading")

        # ── Tab 3: Analysis in R ───────────────────────────────────────────────
        r_analysis_tab = QWidget()
        r_layout       = QVBoxLayout(r_analysis_tab)
        r_sub_tabs     = QTabWidget()
        r_layout.addWidget(r_sub_tabs)

        self.family_selection_widget  = FamilySelectionWidget()
        self.contrast_selection_widget = ContrastSelectionWidget()
        self.bam_widget               = BamWidget()
        self.bam_widget.set_contrast_selection_widget(self.contrast_selection_widget)

        r_sub_tabs.addTab(self.family_selection_widget,  "1. Family Selection (Optional)")
        r_sub_tabs.addTab(self.contrast_selection_widget, "2. Contrast Selection")
        r_sub_tabs.addTab(self.bam_widget,               "3. BAM Analysis")

        self.main_tabs.addTab(r_analysis_tab, "Analysis in R")

        # ── Signal connections ─────────────────────────────────────────────────

        # Condition colors → plate format + data loader + metadata assignment
        self.conditions_format.conditions_updated.connect(
            self.plate_format_widget.update_conditions
        )
        self.conditions_format.conditions_updated.connect(
            self.data_loader.update_conditions
        )
        self.conditions_format.conditions_updated.connect(
            self.metadata_assignment.update_from_conditions
        )
        self.conditions_format.conditions_updated.connect(
            self.contrast_selection_widget.update_conditions
        )
        # Also push the per-variable breakdown to metadata for reference derivation
        self.conditions_format.conditions_updated.connect(
            self._push_condition_variables
        )

        # Variable names → metadata, bam, family selection, data loader
        self.conditions_format.variables_updated.connect(
            self.metadata_assignment.set_variable_names
        )
        self.conditions_format.variables_updated.connect(
            self.bam_widget.set_variable_names
        )
        self.conditions_format.variables_updated.connect(
            self.family_selection_widget.set_variable_names
        )
        self.conditions_format.variables_updated.connect(
            self.data_loader.set_variable_names
        )

        # Roles → plate format (for button labels)
        self.metadata_assignment.roles_updated.connect(
            self.plate_format_widget.update_roles
        )
        self.metadata_assignment.roles_updated.connect(
            self.bam_widget.set_roles
        )

        # Plate layout library → data loader
        self.plate_format_widget.layouts_changed.connect(
            self.data_loader.update_layouts
        )

        # Reference controls → analysis widgets
        self.metadata_assignment.references_updated.connect(self._on_references_updated)

        # Data loading → sanity check → fan-out to analysis widgets
        self.data_loader.data_ready.connect(self.sanity_check_viewer.load_data)
        self.sanity_check_viewer.data_passed_through.connect(self._on_data_ready)

        # Family selection → BAM analysis
        self.family_selection_widget.family_changed.connect(
            self.bam_widget.set_family
        )

        # Contrast changes → mark BAM results stale
        self.contrast_selection_widget.selection_changed.connect(
            self.bam_widget.on_contrasts_changed
        )

        # BAM completion → warn on contrast tab that re-running is required
        self.bam_widget.analysis_complete.connect(
            lambda _output_dir: self.contrast_selection_widget.on_bam_completed()
        )

        # Seed with initial conditions
        print("[DBG] seeding conditions", flush=True)
        self.conditions_format.emit_conditions_data()
        print("[DBG] seeding plate layouts", flush=True)
        self.plate_format_widget.emit_layouts_state()
        print("[DBG] MainWindow.__init__ END", flush=True)

    # ── Output directory ────────────────────────────────────────────────────────

    def _choose_output_dir(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if not directory:
            return
        self._output_dir = directory
        self._output_dir_label.setText(directory)
        self.bam_widget.set_output_dir(directory)
        self.family_selection_widget.set_output_dir(directory)

    # ── Condition variable push ──────────────────────────────────────────────

    def _push_condition_variables(self, _conditions_dict):
        """Pushes per-variable breakdown to metadata widget whenever conditions change."""
        cond_vars = self.conditions_format.get_condition_variables()
        self.metadata_assignment.update_condition_variables(cond_vars)

    # ── Reference controls ────────────────────────────────────────────────────

    def _on_references_updated(self, variable_refs: dict, ref_condition: str):
        self.bam_widget.set_references(variable_refs, ref_condition)
        self.family_selection_widget.set_references(variable_refs, ref_condition)

    # ── Data fan-out ──────────────────────────────────────────────────────────

    def _on_data_ready(self, df):
        """Fans out the processed full_df to all analysis tabs."""
        if df is None or df.empty:
            return
        self.bam_widget.load_data(df)
        self.family_selection_widget.load_data(df)

    # ── Save / Load config ────────────────────────────────────────────────────

    def _save_config(self):
        """Saves ALL experimental design inputs to a JSON file."""
        default_path = os.path.join(os.path.expanduser("~"), "Documents",
                                    "experiment_config.json")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Experimental Design Config",
            default_path, "JSON Files (*.json)"
        )
        if not path:
            return

        config = {
            "phases":              self.experiment_format.get_phases_data(),
            "conditions":          self.conditions_format.save_conditions_config(),
            "roles":               self.metadata_assignment.save_roles_config(),
            "plate_layouts":       self.plate_format_widget.get_all_layouts(),
            "contrast_selection":  self.contrast_selection_widget.save_config(),
        }

        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(config, f, indent=2)
            QMessageBox.information(self, "Saved",
                                    f"Experimental design config saved to:\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Save Error", str(e))

    def _load_config(self):
        """Loads a previously saved JSON config and restores all design widgets."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Experimental Design Config",
            os.path.expanduser("~"), "JSON Files (*.json)"
        )
        if not path:
            return

        try:
            with open(path, "r", encoding="utf-8") as f:
                config = json.load(f)
        except Exception as e:
            QMessageBox.critical(self, "Load Error", f"Could not read config file:\n{e}")
            return

        errors = []

        # 1. Restore phases
        phases = config.get("phases")
        if phases:
            try:
                self.experiment_format.load_phases_config(phases)
            except Exception as e:
                errors.append(f"Phases: {e}")

        # 2. Restore conditions (must happen before plate/roles so colors propagate)
        conditions = config.get("conditions")
        if conditions:
            try:
                self.conditions_format.load_conditions_config(conditions)
            except Exception as e:
                errors.append(f"Conditions: {e}")

        # 3. Restore roles (must happen after conditions populate the table)
        roles = config.get("roles")
        if roles:
            try:
                self.metadata_assignment.load_roles_config(roles)
            except Exception as e:
                errors.append(f"Roles: {e}")

        # 3b. Restore contrast selection (must happen after conditions populate
        #     so the pair list exists to apply excluded_pairs to)
        cs = config.get("contrast_selection")
        if cs:
            try:
                self.contrast_selection_widget.load_config(cs)
            except Exception as e:
                errors.append(f"Contrast selection: {e}")

        # 4. Restore plate layout library
        layouts = config.get("plate_layouts")
        if layouts:
            try:
                self.plate_format_widget.load_layouts_config(layouts)
            except Exception as e:
                errors.append(f"Plate layouts: {e}")
        else:
            fmt    = config.get("plate_format")
            layout = config.get("plate_layout")
            if fmt and layout:
                try:
                    self.plate_format_widget.load_layouts_config(
                        {"Default": {"format": fmt, "map": layout}}
                    )
                except Exception as e:
                    errors.append(f"Plate layout (legacy): {e}")

        # 5. Legacy support: old configs with ref_genotype / ref_drug / ref_condition
        if not roles:
            ref_condition = config.get("ref_condition", "")
            if ref_condition:
                try:
                    # Mark that condition as Reference Control
                    legacy_roles = {ref_condition: "Reference Control"}
                    self.metadata_assignment.load_roles_config(legacy_roles)
                except Exception as e:
                    errors.append(f"Legacy reference controls: {e}")

        if errors:
            QMessageBox.warning(
                self, "Loaded with Warnings",
                "Config loaded but some sections had errors:\n\n" + "\n".join(errors)
            )
        else:
            QMessageBox.information(self, "Loaded",
                                    "Experimental design config loaded successfully.")


# =============================================================================
#  SPINBOX STYLESHEET
# =============================================================================

_SPINBOX_STYLESHEET = """
    QSpinBox {
        width: 60px;
        padding-right: 15px;
    }
    QSpinBox::up-button {
        subcontrol-origin: border;
        subcontrol-position: center right;
        width: 16px;
        border: none;
        background: transparent;
    }
    QSpinBox::down-button {
        subcontrol-origin: border;
        subcontrol-position: center right;
        width: 16px;
        right: 17px;
        border: none;
        background: transparent;
    }
    QSpinBox::up-button:hover, QSpinBox::down-button:hover {
        background-color: #3C3C3C;
        border-radius: 3px;
        border: 1px solid #888;
    }
"""


# =============================================================================
#  ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    print("[DBG] __main__: creating QApplication", flush=True)
    app = QApplication(sys.argv)
    print("[DBG] __main__: QApplication created", flush=True)

    try:
        import qdarktheme
        app.setStyleSheet(qdarktheme.load_stylesheet())
    except ImportError:
        print("qdarktheme not found — using system theme.")
        print("Install with: pip install pyqtdarktheme")

    print("[DBG] about to construct MainWindow", flush=True)
    window = MainWindow()
    print("[DBG] about to call window.show()", flush=True)
    window.show()
    print("[DBG] about to enter app.exec()", flush=True)
    sys.exit(app.exec())

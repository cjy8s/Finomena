import os
import shutil

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QScrollArea, QFileDialog, QMessageBox, QSizePolicy
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap, QFont


class FigureViewerWidget(QWidget):
    """
    A reusable widget that displays a sequence of figures (QPixmaps) with
    navigation arrows and a Save All button.

    Usage:
        viewer = FigureViewerWidget(title="BAM Results")
        viewer.load_figures([
            {'pixmap': qpixmap, 'title': 'Figure 1', 'filepath': '/path/to/file.pdf'},
            ...
        ])
    """

    def __init__(self, title: str = "Figures", parent=None):
        super().__init__(parent)
        self._title_text = title
        self._figures = []       # list of {'pixmap', 'title', 'filepath'}
        self._current_index = 0
        self._build_ui()

    # ── Public API ────────────────────────────────────────────────────────────

    def load_figures(self, figures: list):
        """
        Replaces the current figure set.

        Parameters
        ----------
        figures : list of dict
            Each dict must contain:
              'pixmap'  : QPixmap  — the rendered image
              'title'   : str      — short label shown below the image
              'filepath': str|None — original file path (used by Save All)
        """
        self._figures = figures
        self._current_index = 0
        self._save_button.setEnabled(bool(figures))
        if figures:
            self._show_figure(0)
        else:
            self._image_label.clear()
            self._image_label.setText("No figures loaded.")
            self._caption_label.clear()
            self._page_label.setText("No figures")
        self._update_nav_buttons()

    def clear(self):
        """Removes all figures and resets the widget."""
        self._figures = []
        self._current_index = 0
        self._image_label.clear()
        self._image_label.setText("No figures loaded.")
        self._caption_label.clear()
        self._page_label.setText("")
        self._save_button.setEnabled(False)
        self._update_nav_buttons()

    # ── Private: build UI ─────────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        # Title label
        title_label = QLabel(self._title_text)
        title_font = QFont()
        title_font.setBold(True)
        title_label.setFont(title_font)
        layout.addWidget(title_label)

        # Navigation bar
        nav_bar = QHBoxLayout()
        self._prev_button = QPushButton("◀")
        self._prev_button.setFixedWidth(40)
        self._prev_button.clicked.connect(self._on_prev)
        self._prev_button.setEnabled(False)

        self._page_label = QLabel("")
        self._page_label.setAlignment(Qt.AlignCenter)
        self._page_label.setMinimumWidth(120)

        self._next_button = QPushButton("▶")
        self._next_button.setFixedWidth(40)
        self._next_button.clicked.connect(self._on_next)
        self._next_button.setEnabled(False)

        self._save_button = QPushButton("Save All")
        self._save_button.setEnabled(False)
        self._save_button.clicked.connect(self._on_save_all)

        nav_bar.addWidget(self._prev_button)
        nav_bar.addWidget(self._page_label)
        nav_bar.addWidget(self._next_button)
        nav_bar.addStretch()
        nav_bar.addWidget(self._save_button)
        layout.addLayout(nav_bar)

        # Scroll area for the figure
        self._scroll_area = QScrollArea()
        self._scroll_area.setWidgetResizable(True)
        self._scroll_area.setAlignment(Qt.AlignCenter)

        self._image_label = QLabel("No figures loaded.")
        self._image_label.setAlignment(Qt.AlignCenter)
        self._image_label.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored)
        self._scroll_area.setWidget(self._image_label)
        layout.addWidget(self._scroll_area, stretch=1)

        # Caption below the figure
        self._caption_label = QLabel("")
        caption_font = QFont()
        caption_font.setItalic(True)
        caption_font.setPointSize(9)
        self._caption_label.setFont(caption_font)
        self._caption_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self._caption_label)

    # ── Private: navigation ───────────────────────────────────────────────────

    def _show_figure(self, index: int):
        if not self._figures or index < 0 or index >= len(self._figures):
            return
        self._current_index = index
        fig = self._figures[index]
        pixmap: QPixmap = fig['pixmap']

        # Scale to fit the scroll area while keeping aspect ratio
        available = self._scroll_area.viewport().size()
        scaled = pixmap.scaled(
            available,
            Qt.KeepAspectRatio,
            Qt.SmoothTransformation
        )
        self._image_label.setPixmap(scaled)
        self._caption_label.setText(fig.get('title', ''))
        self._page_label.setText(f"Figure {index + 1} / {len(self._figures)}")
        self._update_nav_buttons()

    def _on_prev(self):
        if self._current_index > 0:
            self._show_figure(self._current_index - 1)

    def _on_next(self):
        if self._current_index < len(self._figures) - 1:
            self._show_figure(self._current_index + 1)

    def _update_nav_buttons(self):
        n = len(self._figures)
        self._prev_button.setEnabled(self._current_index > 0)
        self._next_button.setEnabled(self._current_index < n - 1)

    # ── Private: save ─────────────────────────────────────────────────────────

    def _on_save_all(self):
        dest_dir = QFileDialog.getExistingDirectory(
            self, "Select Folder to Save Figures"
        )
        if not dest_dir:
            return

        saved = 0
        errors = []
        for fig in self._figures:
            filepath = fig.get('filepath')
            title = fig.get('title', f'figure_{saved + 1}')
            try:
                if filepath and os.path.isfile(filepath):
                    dest = os.path.join(dest_dir, os.path.basename(filepath))
                    shutil.copy2(filepath, dest)
                else:
                    # Save pixmap as PNG
                    safe_title = "".join(c if c.isalnum() or c in (' ', '-', '_') else '_' for c in title)
                    dest = os.path.join(dest_dir, f"{safe_title}.png")
                    fig['pixmap'].save(dest, "PNG")
                saved += 1
            except Exception as e:
                errors.append(f"{title}: {e}")

        if errors:
            QMessageBox.warning(
                self, "Save Complete with Errors",
                f"Saved {saved} figure(s).\n\nErrors:\n" + "\n".join(errors)
            )
        else:
            QMessageBox.information(
                self, "Saved",
                f"All {saved} figure(s) saved to:\n{dest_dir}"
            )

    # ── Resize event ──────────────────────────────────────────────────────────

    def resizeEvent(self, event):
        super().resizeEvent(event)
        # Re-scale the current figure when the window is resized
        if self._figures and 0 <= self._current_index < len(self._figures):
            self._show_figure(self._current_index)

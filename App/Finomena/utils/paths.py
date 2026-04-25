"""
paths.py — Centralized path resolution for development and PyInstaller frozen mode.

In development:
    __file__ resolves normally; APP_DIR = App/ (two parents up from utils/).

When frozen by PyInstaller (--onedir):
    sys.executable is inside dist/Finomena/; APP_DIR = that directory.
    Bundled data files (R script, README) are extracted to sys._MEIPASS.
    External data (R library, bundled R binary) lives alongside the executable.
"""

import os
import sys


def _get_app_dir() -> str:
    """Returns the App/ root for resolving source-tree resources."""
    if getattr(sys, 'frozen', False):
        return getattr(sys, '_MEIPASS', os.path.dirname(sys.executable))
    # Development: this file is at App/Finomena/utils/paths.py
    this_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.dirname(os.path.dirname(this_dir))


APP_DIR = _get_app_dir()


def resource_path(*parts: str) -> str:
    """
    Resolves a path relative to the App root.
    Works in both development and frozen mode.

    Usage:
        resource_path("R", "scripts", "TweedieAR1 BAM.R")
        resource_path("README.md")
    """
    return os.path.join(APP_DIR, *parts)


def external_data_dir() -> str:
    """
    Returns the directory where large external data lives (R library, R binary).
    In frozen mode this is the executable's directory (not _MEIPASS, which is
    a temp extraction folder). In development mode it equals APP_DIR.
    """
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return APP_DIR

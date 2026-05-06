# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller spec file for Finomena: Zebrafish Behavioral Analysis Suite.

Build with the platform build script (build_windows.bat, build_mac.sh, or
build_linux.sh), which runs pyinstaller and bundles R automatically.

Manual build:
    pyinstaller finomena.spec --noconfirm
"""

import os
import sys
import glob as _glob
from PyInstaller.utils.hooks import collect_submodules, collect_data_files

block_cipher = None

# ── Collect ALL submodules for packages with complex internal structure ────────
# PyInstaller's auto-detection misses many submodules for these packages
# (namespace packages, Cython extensions, lazy imports, internal compat shims).
# collect_submodules() walks the full package tree so nothing is left behind.
_all_hidden = []
_all_datas  = []

for _pkg in [
    'jaraco',          # namespace pkg (no __init__.py), pulled in by pkg_resources
    'sklearn',         # 680 submodules, Cython extensions
    'scipy',           # 1058 submodules, array_api_compat.numpy.fft etc.
    'numpy',           # 470 submodules, MKL backends
    'pandas',          # 1440 submodules, Cython internals
    'matplotlib',      # 162 submodules, backends
    'seaborn',         # 53 submodules
    'joblib',          # 72 submodules (used by sklearn)
    'threadpoolctl',   # used by sklearn
    'PIL',             # 103 submodules (Pillow image codecs)
    'qdarktheme',      # 37 submodules
    'fitz',            # PyMuPDF C extension
    'PySide6',         # Qt GUI framework
    # ── Color palette packages ──
    'distinctipy',     # maximally distinct categorical colors
    'glasbey',         # optimised categorical palettes
    'colorcet',        # perceptually uniform colormaps (256-color glassbey etc.)
    'palettable',      # ColorBrewer, CartoCColors, scientific palettes
    'cmocean',         # oceanographic scientific colormaps
    'cmasher',         # scientific colormaps (ember, neon, pride, etc.)
    'husl',            # HSLuv perceptual color space
    'vapeplot',        # aesthetic/themed palettes
    'bokeh',           # Category20, Turbo256, Viridis256
    'colorspacious',   # dependency of cmasher
    'cycler',          # color cycling utility
]:
    try:
        _all_hidden += collect_submodules(_pkg)
    except Exception:
        pass  # package not installed — skip

for _pkg in ['jaraco', 'colorcet', 'palettable', 'cmocean', 'cmasher', 'vapeplot']:
    try:
        _all_datas += collect_data_files(_pkg)
    except Exception:
        pass

# ── Conda fix: many DLLs live in Library/bin/, not alongside .pyd files ───────
# This covers PySide6/Qt6, Pillow image codecs, MKL/BLAS, yaml, zmq, etc.
_conda_bin = os.path.join(sys.prefix, 'Library', 'bin')
_extra_binaries = []
if os.path.isdir(_conda_bin):
    for pattern in [
        # PySide6 / Qt
        'Qt6Core.dll', 'Qt6Gui.dll', 'Qt6Widgets.dll', 'Qt6Network.dll',
        'Qt6Svg.dll', 'Qt6OpenGL.dll', 'Qt6OpenGLWidgets.dll', 'Qt6Pdf.dll',
        'shiboken6*.dll', 'pyside6*.dll',
        # Pillow image codec DLLs
        'freetype.dll', 'libjpeg*.dll', 'tiff.dll', 'openjp2.dll',
        'libwebp*.dll', 'avif.dll', 'lcms2.dll', 'zlib*.dll',
        'libpng*.dll', 'liblzma*.dll',
        # MKL / BLAS (numpy, scipy)
        'mkl_rt*.dll', 'libgfortran*.dll', 'libifcoremd.dll', 'libmmd.dll',
        'libblas.dll', 'liblapack.dll', 'tbb12.dll',
        # Other conda libs
        'yaml.dll', 'libzmq*.dll',
    ]:
        for dll in _glob.glob(os.path.join(_conda_bin, pattern)):
            _extra_binaries.append((dll, '.'))

a = Analysis(
    ['app.py'],
    pathex=['Finomena/utils'],
    binaries=_extra_binaries,
    datas=[
        # R scripts (small, bundled inside the app)
        ('R/scripts/TweedieAR1 BAM.R',  'R/scripts'),
        ('R/scripts/FamilySelection.R', 'R/scripts'),
        # README
        ('README.md', '.'),
    ] + _all_datas,
    hiddenimports=_all_hidden,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        # Qt bindings we don't use
        'PyQt6', 'PyQt5', 'PyQt6.sip', 'PyQt5.sip',
        # ── Packages installed in the conda env but NOT used by Finomena ──
        # Excluding these prevents DLL-load errors (pyarrow/arrow.dll,
        # numba/llvmlite LLVM) and dramatically shrinks the bundle.
        'pyarrow',
        'numba', 'llvmlite',
        'dask', 'distributed',
        'bokeh',
        'IPython', 'ipykernel', 'jupyter', 'notebook', 'tornado',
        'pyts', 'sktime', 'tslearn', 'stumpy', 'tsfresh',
        'hdbscan', 'umap',
        'statsmodels', 'pingouin', 'patsy',
        'rpy2',
        'xarray', 'polars',
    ],
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='Finomena',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,  # TODO: add app icon (.ico on Windows, .icns on Mac)
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='Finomena',
)

# macOS .app bundle (ignored on Windows)
if sys.platform == 'darwin':
    app = BUNDLE(
        coll,
        name='Finomena.app',
        icon=None,  # TODO: add .icns file
        bundle_identifier='com.finomena.app',
        info_plist={
            'NSHighResolutionCapable': True,
            'CFBundleShortVersionString': '1.0.0',
        },
    )

#!/bin/bash
# Build Finomena for macOS
# Prerequisites: pip install pyinstaller
#
# This creates dist/Finomena.app with a bundled portable R — fully
# self-contained, no system R required on the target Mac.

set -e

echo "Building Finomena..."
pyinstaller finomena.spec --noconfirm

# ── Detect system R installation ──────────────────────────────────────────────
R_SOURCE=""

# CRAN .pkg installer (most common)
if [ -d "/Library/Frameworks/R.framework/Resources" ]; then
    R_SOURCE="/Library/Frameworks/R.framework/Resources"
# Homebrew Apple Silicon
elif [ -d "/opt/homebrew/Cellar/r" ]; then
    R_SOURCE="$(ls -d /opt/homebrew/Cellar/r/*/lib/R 2>/dev/null | sort -V | tail -1)"
# Homebrew Intel
elif [ -d "/usr/local/Cellar/r" ]; then
    R_SOURCE="$(ls -d /usr/local/Cellar/r/*/lib/R 2>/dev/null | sort -V | tail -1)"
fi

if [ -z "$R_SOURCE" ] || [ ! -d "$R_SOURCE" ]; then
    echo "ERROR: Could not find an R installation."
    echo "       Install R from https://cran.r-project.org or: brew install r"
    exit 1
fi

echo "Found R at: $R_SOURCE"
echo "Copying full R installation into dist/Finomena/R/ ..."

# Copy R binaries, base packages, and configuration
mkdir -p "dist/Finomena/R"
for subdir in bin etc lib library modules share; do
    if [ -d "$R_SOURCE/$subdir" ]; then
        cp -R "$R_SOURCE/$subdir" "dist/Finomena/R/$subdir"
    fi
done

echo "R installation bundled."

# ── Install app-specific R packages into the bundled library ──────────────────
echo "Installing required R packages into bundled library..."
"dist/Finomena/R/bin/Rscript" R/install_packages.R

echo ""
echo "========================================"
echo "Build complete: dist/Finomena.app"
echo "========================================"
echo ""
echo "To run:  open dist/Finomena.app"
echo ""
echo "The app includes a bundled R installation — no system R required."
echo ""
echo "If macOS blocks the app (Gatekeeper), right-click > Open,"
echo "or run: xattr -d com.apple.quarantine dist/Finomena.app"

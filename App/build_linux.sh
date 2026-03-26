#!/bin/bash
# Build Finomena for Linux
# Prerequisites: pip install pyinstaller
#
# This creates dist/Finomena/ with the executable and a bundled portable R —
# fully self-contained, no system R required on the target machine.

set -e

echo "Building Finomena..."
pyinstaller finomena.spec --noconfirm

# ── Detect system R installation ──────────────────────────────────────────────
R_SOURCE=""

# R installed via package manager (apt, dnf, etc.)
if [ -d "/usr/lib/R" ]; then
    R_SOURCE="/usr/lib/R"
elif [ -d "/usr/lib64/R" ]; then
    R_SOURCE="/usr/lib64/R"
elif [ -d "/usr/local/lib/R" ]; then
    R_SOURCE="/usr/local/lib/R"
fi

if [ -z "$R_SOURCE" ] || [ ! -d "$R_SOURCE" ]; then
    echo "ERROR: Could not find an R installation."
    echo "       Install with: sudo apt install r-base  (Debian/Ubuntu)"
    echo "                 or: sudo dnf install R       (Fedora/RHEL)"
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

# On some Linux distros, Rscript is in /usr/bin, not inside the R directory.
# Ensure it's present in the bundled R/bin/.
if [ ! -f "dist/Finomena/R/bin/Rscript" ]; then
    RSCRIPT_PATH="$(which Rscript 2>/dev/null || true)"
    if [ -n "$RSCRIPT_PATH" ]; then
        cp "$RSCRIPT_PATH" "dist/Finomena/R/bin/Rscript"
    fi
fi

echo "R installation bundled."

# ── Install app-specific R packages into the bundled library ──────────────────
echo "Installing required R packages into bundled library..."
"dist/Finomena/R/bin/Rscript" R/install_packages.R

echo ""
echo "========================================"
echo "Build complete: dist/Finomena/"
echo "========================================"
echo ""
echo "To run:  ./dist/Finomena/Finomena"
echo ""
echo "The app includes a bundled R installation — no system R required."

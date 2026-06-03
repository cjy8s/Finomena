#!/bin/bash
# ============================================================
#  Finomena — macOS build script
#  Output: Finomena_macOS.dmg
#
#  Prerequisites (install once on your build Mac):
#    pip install pyinstaller
#    brew install create-dmg
#
#  R 4.5.3 is downloaded automatically from CRAN on first run
#  and cached in /tmp for subsequent builds.
#
#  Run from the App/ directory:
#    bash build_mac.sh
# ============================================================

set -euo pipefail

R_VERSION="4.5.3"

# Detect architecture and choose the correct CRAN package
ARCH="$(uname -m)"
if [ "$ARCH" = "arm64" ]; then
    R_PKG="R-${R_VERSION}-arm64.pkg"
    R_URL="https://cran.r-project.org/bin/macosx/big-sur-arm64/base/${R_PKG}"
else
    R_PKG="R-${R_VERSION}-x86_64.pkg"
    R_URL="https://cran.r-project.org/bin/macosx/big-sur-x86_64/base/${R_PKG}"
fi

R_INSTALLER="/tmp/${R_PKG}"
R_FRAMEWORK_VERSION="4.5"
R_FRAMEWORK_RESOURCES="/Library/Frameworks/R.framework/Versions/${R_FRAMEWORK_VERSION}/Resources"

APP_BUNDLE="dist/Finomena.app"
# external_data_dir() in frozen .app = Contents/MacOS/
R_DEST="${APP_BUNDLE}/Contents/MacOS/R"
R_RUNTIME="${R_DEST}/runtime"
R_LIBRARY="${R_DEST}/library"

echo "============================================================"
echo " Building Finomena for macOS (${ARCH})   R ${R_VERSION}"
echo "============================================================"

# ── Step 1: PyInstaller ───────────────────────────────────────
echo ""
echo "[1/4] Packaging Python app with PyInstaller..."
pyinstaller finomena.spec --noconfirm

# ── Step 2: Download & install R (cached after first run) ─────
echo ""
echo "[2/4] Setting up R ${R_VERSION}..."

if [ -f "$R_INSTALLER" ]; then
    echo "  Using cached installer: ${R_INSTALLER}"
else
    echo "  Downloading R ${R_VERSION} from CRAN (${ARCH})..."
    curl -fL --progress-bar -o "$R_INSTALLER" "$R_URL"
fi

# Install R to the system framework location so we can copy it
# (requires sudo — R's pkg installer needs /Library/Frameworks write access)
echo "  Installing R ${R_VERSION} to /Library/Frameworks/R.framework..."
echo "  (sudo required for pkg installation)"
sudo installer -pkg "$R_INSTALLER" -target /

# Copy R runtime into the app bundle (bin, etc, lib, library, modules, share)
echo "  Copying R runtime into app bundle..."
rm -rf "$R_RUNTIME"
mkdir -p "$R_RUNTIME"
for subdir in bin etc lib library modules share; do
    if [ -d "${R_FRAMEWORK_RESOURCES}/${subdir}" ]; then
        cp -R "${R_FRAMEWORK_RESOURCES}/${subdir}" "${R_RUNTIME}/${subdir}"
    fi
done

# ── Step 3: Install app R packages into separate library ──────
echo ""
echo "[3/4] Installing R packages into app bundle library..."
mkdir -p "$R_LIBRARY"
R_HOME="$R_RUNTIME" "${R_RUNTIME}/bin/Rscript" \
    R/install_packages.R "$(pwd)/${R_LIBRARY}"

# ── Step 4: Create .dmg ───────────────────────────────────────
echo ""
echo "[4/4] Creating macOS disk image..."

if ! command -v create-dmg &>/dev/null; then
    echo ""
    echo "  WARNING: create-dmg not found — skipping .dmg creation."
    echo "  Install with: brew install create-dmg"
    echo "  The .app bundle is at ${APP_BUNDLE}"
    exit 0
fi

rm -f Finomena_macOS.dmg

create-dmg \
    --volname "Finomena ${R_VERSION}" \
    --window-pos 200 120 \
    --window-size 600 320 \
    --icon-size 100 \
    --icon "Finomena.app" 175 140 \
    --hide-extension "Finomena.app" \
    --app-drop-link 425 140 \
    --no-internet-enable \
    "Finomena_macOS.dmg" \
    "$APP_BUNDLE"

echo ""
echo "============================================================"
echo " Build complete!"
echo "   Installer : Finomena_macOS.dmg"
echo "   App bundle: ${APP_BUNDLE}"
echo "============================================================"
echo ""
echo "Upload Finomena_macOS.dmg to a GitHub Release for direct"
echo "download by end users."
echo ""

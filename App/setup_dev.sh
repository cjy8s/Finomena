#!/bin/bash
# ============================================================
#  Finomena — macOS dev environment setup
#
#  Run this ONCE after cloning the repo. It downloads R 4.5.3
#  from CRAN and installs it into App/R/runtime/ so the app
#  works without a system R installation.
#
#  Run from the repo root:
#    bash App/setup_dev.sh
# ============================================================

set -euo pipefail

R_VERSION="4.5.3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_RUNTIME="${SCRIPT_DIR}/R/runtime"
R_LIBRARY="${SCRIPT_DIR}/R/library"

# Detect architecture
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

echo "============================================================"
echo " Finomena dev setup   (macOS ${ARCH}, R ${R_VERSION})"
echo "============================================================"

if [ -f "${R_RUNTIME}/bin/Rscript" ]; then
    echo ""
    echo "R runtime already present at: ${R_RUNTIME}"
    read -rp "Reinstall? [y/N] " REINSTALL
    if [ "${REINSTALL,,}" != "y" ]; then
        echo "Skipping — dev environment is already set up."
        goto_packages=1
    fi
fi

if [ -z "${goto_packages:-}" ]; then
    # ── Download R installer ──────────────────────────────────
    echo ""
    echo "[1/2] Downloading R ${R_VERSION} from CRAN (${ARCH})..."
    if [ -f "$R_INSTALLER" ]; then
        echo "  Using cached installer: ${R_INSTALLER}"
    else
        curl -fL --progress-bar -o "$R_INSTALLER" "$R_URL"
    fi

    # ── Install R to system framework, then copy into runtime/ ─
    echo "  Installing R (sudo required for pkg installation)..."
    sudo installer -pkg "$R_INSTALLER" -target /

    echo "  Copying R runtime into ${R_RUNTIME}/ ..."
    rm -rf "$R_RUNTIME"
    mkdir -p "$R_RUNTIME"
    for subdir in bin etc lib library modules share; do
        if [ -d "${R_FRAMEWORK_RESOURCES}/${subdir}" ]; then
            cp -R "${R_FRAMEWORK_RESOURCES}/${subdir}" "${R_RUNTIME}/${subdir}"
        fi
    done
fi

# ── Install app R packages into App/R/library/ ───────────────
echo ""
echo "[2/2] Installing R packages into ${R_LIBRARY}/ ..."
mkdir -p "$R_LIBRARY"
R_HOME="$R_RUNTIME" "${R_RUNTIME}/bin/Rscript" \
    "${SCRIPT_DIR}/R/install_packages.R"

echo ""
echo "============================================================"
echo " Dev setup complete!"
echo "  Bundled R : ${R_RUNTIME}"
echo "  R packages: ${R_LIBRARY}"
echo "  Run the app with:  python App/app.py"
echo "============================================================"
echo ""

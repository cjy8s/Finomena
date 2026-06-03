@echo off
setlocal EnableDelayedExpansion

REM ============================================================
REM  Finomena — Windows build script
REM  Output: Output\Finomena_Setup_Windows.exe
REM
REM  Prerequisites (install once on your build machine):
REM    pip install pyinstaller
REM    Inno Setup 6  https://jrsoftware.org/isinfo.php
REM
REM  R 4.5.3 is downloaded automatically from CRAN on first run
REM  and cached in %TEMP% for subsequent builds.
REM ============================================================

set R_VERSION=4.5.3
set R_URL=https://cran.r-project.org/bin/windows/base/R-%R_VERSION%-win.exe
set R_INSTALLER=%TEMP%\R-%R_VERSION%-win.exe
set R_RUNTIME=dist\Finomena\R\runtime
set R_LIBRARY=dist\Finomena\R\library

echo ============================================================
echo  Building Finomena for Windows   (R %R_VERSION%)
echo ============================================================

REM ── Step 1: PyInstaller ─────────────────────────────────────
echo.
echo [1/4] Packaging Python app with PyInstaller...
pyinstaller finomena.spec --noconfirm
if errorlevel 1 ( echo. & echo ERROR: PyInstaller failed. & exit /b 1 )

REM ── Step 2: Download R installer (cached after first run) ───
echo.
echo [2/4] Setting up R %R_VERSION%...
if exist "%R_INSTALLER%" (
    echo   Using cached installer: %R_INSTALLER%
) else (
    echo   Downloading R %R_VERSION% from CRAN...
    powershell -NoProfile -Command ^
        "Invoke-WebRequest -Uri '%R_URL%' -OutFile '%R_INSTALLER%' -UseBasicParsing"
    if errorlevel 1 ( echo. & echo ERROR: Download failed. & exit /b 1 )
    echo   Saved to: %R_INSTALLER%
)

REM Extract R into dist\Finomena\R\runtime\
echo   Extracting R into %R_RUNTIME%\...
if exist "%R_RUNTIME%" rmdir /s /q "%R_RUNTIME%"
"%R_INSTALLER%" /VERYSILENT /SUPPRESSMSGBOXES /NORESTART /DIR="%CD%\%R_RUNTIME%"
if errorlevel 1 ( echo. & echo ERROR: R extraction failed. & exit /b 1 )

REM ── Step 3: Install app R packages into separate library ────
echo.
echo [3/4] Installing R packages into %R_LIBRARY%\...
if not exist "%R_LIBRARY%" mkdir "%R_LIBRARY%"
"%R_RUNTIME%\bin\Rscript.exe" R\install_packages.R "%CD%\%R_LIBRARY%"
if errorlevel 1 ( echo. & echo ERROR: R package installation failed. & exit /b 1 )

REM ── Step 4: Create installer with Inno Setup ────────────────
echo.
echo [4/4] Creating Windows installer with Inno Setup...
set "ISCC="
for %%P in (
    "%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe"
    "%ProgramFiles%\Inno Setup 6\ISCC.exe"
) do if exist %%P set "ISCC=%%P"

if not defined ISCC (
    echo.
    echo   WARNING: Inno Setup 6 not found — skipping installer creation.
    echo   Install from https://jrsoftware.org/isinfo.php then re-run.
    echo   The unpackaged app is at dist\Finomena\Finomena.exe
    goto :done
)

if not exist Output mkdir Output
"%ISCC%" finomena.iss
if errorlevel 1 ( echo. & echo ERROR: Inno Setup failed. & exit /b 1 )

:done
echo.
echo ============================================================
echo  Build complete!
if exist "Output\Finomena_Setup_Windows.exe" (
    echo   Installer : Output\Finomena_Setup_Windows.exe
)
echo   App folder: dist\Finomena\Finomena.exe
echo ============================================================
echo.
echo Upload Output\Finomena_Setup_Windows.exe to a GitHub Release
echo for direct download by end users.
echo.

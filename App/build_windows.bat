@echo off
REM Build Finomena for Windows
REM Prerequisites: pip install pyinstaller
REM
REM This creates dist\Finomena\ with the .exe, a bundled portable R,
REM and all dependencies — fully self-contained, no system R required.

echo Building Finomena...
pyinstaller finomena.spec --noconfirm

if errorlevel 1 (
    echo Build failed.
    exit /b 1
)

REM ── Bundle the full R installation ──────────────────────────────────────────
REM Detect the system R installation directory.
REM Prefer the newest version found in Program Files.

set "R_SOURCE="
for /d %%D in ("C:\Program Files\R\R-*") do set "R_SOURCE=%%D"
if not defined R_SOURCE (
    for /d %%D in ("C:\Program Files (x86)\R\R-*") do set "R_SOURCE=%%D"
)

if not defined R_SOURCE (
    echo ERROR: Could not find an R installation in Program Files.
    echo        Install R from https://cran.r-project.org and rebuild.
    exit /b 1
)

echo Found R at: %R_SOURCE%
echo Copying full R installation into dist\Finomena\R\ ...

REM Copy the R binaries, base packages, and configuration
xcopy /E /I /Y "%R_SOURCE%\bin"     "dist\Finomena\R\bin"     >nul
xcopy /E /I /Y "%R_SOURCE%\etc"     "dist\Finomena\R\etc"     >nul
xcopy /E /I /Y "%R_SOURCE%\modules" "dist\Finomena\R\modules" >nul
xcopy /E /I /Y "%R_SOURCE%\share"   "dist\Finomena\R\share"   >nul
xcopy /E /I /Y "%R_SOURCE%\library" "dist\Finomena\R\library" >nul

echo R installation bundled.

REM ── Install app-specific R packages into the bundled library ────────────────
echo Installing required R packages into bundled library...
"dist\Finomena\R\bin\Rscript.exe" R\install_packages.R

echo.
echo ========================================
echo Build complete: dist\Finomena\
echo ========================================
echo.
echo To run:  dist\Finomena\Finomena.exe
echo.
echo The app includes a bundled R installation — no system R required.

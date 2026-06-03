@echo off
setlocal EnableDelayedExpansion

REM ============================================================
REM  Finomena — Windows dev environment setup
REM
REM  Run this ONCE after cloning the repo. It downloads R 4.5.3
REM  from CRAN and installs it into App\R\runtime\ so the app
REM  works without a system R installation.
REM
REM  Run from the repo root:
REM    App\setup_dev.bat
REM ============================================================

set R_VERSION=4.5.3
set R_URL=https://cran.r-project.org/bin/windows/base/R-%R_VERSION%-win.exe
set R_INSTALLER=%TEMP%\R-%R_VERSION%-win.exe
set R_RUNTIME=%~dp0R\runtime

echo ============================================================
echo  Finomena dev setup   (R %R_VERSION%)
echo ============================================================

if exist "%R_RUNTIME%\bin\Rscript.exe" (
    echo.
    echo R runtime already present at:
    echo   %R_RUNTIME%
    echo.
    set /p REINSTALL="Reinstall? [y/N] "
    if /i not "!REINSTALL!"=="y" (
        echo Skipping — dev environment is already set up.
        goto :packages
    )
)

REM ── Download R installer ─────────────────────────────────────
echo.
echo [1/2] Downloading R %R_VERSION% from CRAN...
if exist "%R_INSTALLER%" (
    echo   Using cached installer: %R_INSTALLER%
) else (
    powershell -NoProfile -Command ^
        "Invoke-WebRequest -Uri '%R_URL%' -OutFile '%R_INSTALLER%' -UseBasicParsing"
    if errorlevel 1 ( echo. & echo ERROR: Download failed. & exit /b 1 )
)

REM ── Extract R into App\R\runtime\ ────────────────────────────
echo   Extracting R into %R_RUNTIME%\...
if exist "%R_RUNTIME%" rmdir /s /q "%R_RUNTIME%"
"%R_INSTALLER%" /VERYSILENT /SUPPRESSMSGBOXES /NORESTART /DIR="%R_RUNTIME%"
if errorlevel 1 ( echo. & echo ERROR: R extraction failed. & exit /b 1 )

:packages
REM ── Install app R packages into App\R\library\ ───────────────
echo.
echo [2/2] Installing R packages into App\R\library\...
"%R_RUNTIME%\bin\Rscript.exe" "%~dp0R\install_packages.R"
if errorlevel 1 ( echo. & echo ERROR: R package installation failed. & exit /b 1 )

echo.
echo ============================================================
echo  Dev setup complete!
echo  The app will now use the bundled R at:
echo    %R_RUNTIME%
echo  Run the app with:
echo    python App\app.py
echo ============================================================
echo.

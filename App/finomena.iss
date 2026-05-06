; ============================================================
;  Finomena — Inno Setup installer script
;  Builds: Output\Finomena_Setup_Windows.exe
;
;  Run via build_windows.bat, or manually:
;    "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" finomena.iss
; ============================================================

#define AppName      "Finomena"
#define AppVersion   "1.0.0"
#define AppPublisher "Baylor College of Medicine"
#define AppURL       "https://github.com/cjy8s/Finomena"
#define AppExeName   "Finomena.exe"

[Setup]
AppId={{E4A2F3B1-7C6D-4E8F-9A0B-1C2D3E4F5A6B}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher={#AppPublisher}
AppPublisherURL={#AppURL}
AppSupportURL={#AppURL}
AppUpdatesURL={#AppURL}/releases
DefaultDirName={localappdata}\{#AppName}
DisableProgramGroupPage=yes
; No admin rights needed — installs per-user to LocalAppData
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
OutputDir=Output
OutputBaseFilename=Finomena_Setup_Windows
SetupIconFile=
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
; 64-bit Windows only (R 4.5.3 is 64-bit)
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
; Minimum Windows 10
MinVersion=10.0
UninstallDisplayName={#AppName}
UninstallDisplayIcon={app}\{#AppExeName}

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Additional icons:"; Flags: unchecked

[Files]
; Bundle everything from the PyInstaller output directory
Source: "dist\Finomena\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{app}\{#AppExeName}"
Name: "{autodesktop}\{#AppName}";  Filename: "{app}\{#AppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExeName}"; \
    Description: "Launch {#AppName}"; \
    Flags: nowait postinstall skipifsilent

[UninstallDelete]
; Clean up any files written by the app at runtime
Type: filesandordirs; Name: "{app}\Output"
Type: filesandordirs; Name: "{app}\R\library"

@echo off
setlocal
cd /d "%~dp0"
where py >nul 2>nul
if not errorlevel 1 (
  py -3w "%~dp0gui.py"
) else (
  where pythonw >nul 2>nul
  if not errorlevel 1 (
    pythonw "%~dp0gui.py"
  ) else (
    python "%~dp0gui.py"
  )
)

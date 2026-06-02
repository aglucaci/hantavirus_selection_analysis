@echo off
setlocal

set "ROOT_DIR=%~dp0.."
cd /d "%ROOT_DIR%"

if not "%PYTHON%"=="" (
    set "PYTHON_CMD=%PYTHON%"
) else if exist "C:\Users\agluc\anaconda3\python.exe" (
    set "PYTHON_CMD=C:\Users\agluc\anaconda3\python.exe"
) else (
    set "PYTHON_CMD=python"
)

echo [dashboard] Building dashboard-ready TSV/JSON files...
"%PYTHON_CMD%" scripts\build_dashboard_data.py
if errorlevel 1 exit /b %errorlevel%

echo [dashboard] Starting Streamlit dashboard...
echo [dashboard] URL: http://localhost:8501
"%PYTHON_CMD%" -m streamlit run dashboard\andv_gpc_dashboard.py --server.address 0.0.0.0 --server.port 8501

endlocal

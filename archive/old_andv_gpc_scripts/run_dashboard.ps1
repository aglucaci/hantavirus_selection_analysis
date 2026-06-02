$ErrorActionPreference = "Stop"

$Root = Resolve-Path (Join-Path $PSScriptRoot "..")
Set-Location $Root

if ($env:PYTHON) {
    $Python = $env:PYTHON
} elseif (Test-Path "C:\Users\agluc\anaconda3\python.exe") {
    $Python = "C:\Users\agluc\anaconda3\python.exe"
} else {
    $Python = "python"
}

Write-Host "[dashboard] Building dashboard-ready TSV/JSON files..."
& $Python scripts/build_dashboard_data.py

Write-Host "[dashboard] Starting Streamlit dashboard..."
Write-Host "[dashboard] URL: http://localhost:8501"
& $Python -m streamlit run dashboard/andv_gpc_dashboard.py --server.address 0.0.0.0 --server.port 8501

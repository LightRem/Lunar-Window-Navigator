# Скачивает карты высот LOLA (LDEM_*.IMG) в data\terrain\.
# Источник: NASA PDS Geosciences Node, набор LRO-L-LOLA-4-GDR-V1.0; запасной сервер - imbrium.mit.edu.
#
# Запуск (PowerShell, из папки проекта):
#   powershell -ExecutionPolicy Bypass -File scripts\download_terrain.ps1               # все карты (~1.2 ГБ)
#   powershell -ExecutionPolicy Bypass -File scripts\download_terrain.ps1 global south  # только выбранные
# Карты: global (вся Луна, 474 м), south (от 80° ю. ш., 80 м), north (от 80° с. ш., 80 м),
#        south10m (от 87.5° ю. ш., 10 м), all.
# Другая папка: -Dir D:\terrain
# Уже скачанные карты с верной контрольной суммой пропускаются; прерванная загрузка продолжается.
# Нужен curl.exe (есть в Windows 10 1803 и новее).
param(
    [Parameter(Position = 0, ValueFromRemainingArguments = $true)]
    [string[]]$Maps = @('all'),
    [string]$Dir = (Join-Path (Split-Path -Parent $PSScriptRoot) 'data\terrain')
)
$ErrorActionPreference = 'Stop'

$PDS = 'https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr'
$MIT = 'https://imbrium.mit.edu/DATA/LOLA_GDR'

$Catalog = @(
    @{ Key = 'global';   Name = 'LDEM_64';       Bytes = 530841600; Pds = 'cylindrical/img/ldem_64';  Mit = 'CYLINDRICAL/IMG/LDEM_64'
       Sha = '1c4958699d4cffd7e777d51421b6044bba300aa0bef286a544309efb082cdbd6' }
    @{ Key = 'south';    Name = 'LDEM_80S_80M';  Bytes = 115520000; Pds = 'polar/img/ldem_80s_80m';   Mit = 'POLAR/IMG/LDEM_80S_80M'
       Sha = 'e6591eeea35402ca5e124ca00f5bf0774381a2d4fa919144b455d82b8ff65fdc' }
    @{ Key = 'north';    Name = 'LDEM_80N_80M';  Bytes = 115520000; Pds = 'polar/img/ldem_80n_80m';   Mit = 'POLAR/IMG/LDEM_80N_80M'
       Sha = 'dfaa75f2d1337258eaf5e13c0195fcad0d70794911b2daf689026ed12dce2233' }
    @{ Key = 'south10m'; Name = 'LDEM_875S_10M'; Bytes = 460136448; Pds = 'polar/img/ldem_875s_10m';  Mit = 'POLAR/IMG/LDEM_875S_10M'
       Sha = '3c774fd3d8fcab8a9f2e40160821ec9ca83d44ea616f5b7cd13fa78a89f39c53' }
)

$known = @('global', 'south', 'north', 'south10m', 'all')
foreach ($m in $Maps) {
    if ($known -notcontains $m) { [Console]::Error.WriteLine("неизвестная карта '$m' (global, south, north, south10m, all)"); exit 2 }
}
if (-not (Get-Command curl.exe -ErrorAction SilentlyContinue)) { [Console]::Error.WriteLine('нужен curl.exe (Windows 10 1803+)'); exit 1 }
New-Item -ItemType Directory -Force -Path $Dir | Out-Null

# Скачать url в файл, продолжая прерванную загрузку; если докачка не вышла - заново с нуля.
function Get-Url([string]$Url, [string]$Out) {
    & curl.exe -fL --retry 3 --retry-delay 5 -C - -o $Out $Url
    if ($LASTEXITCODE -eq 0) { return $true }
    if (Test-Path $Out) { Remove-Item -Force $Out }
    & curl.exe -fL --retry 3 --retry-delay 5 -o $Out $Url
    return ($LASTEXITCODE -eq 0)
}
function Get-Sha([string]$Path) { (Get-FileHash -Algorithm SHA256 -Path $Path).Hash.ToLower() }
function Get-Size([string]$Path) { (Get-Item $Path).Length }

$failed = $false
foreach ($m in $Catalog) {
    if ($Maps -notcontains 'all' -and $Maps -notcontains $m.Key) { continue }
    $img = Join-Path $Dir "$($m.Name).IMG"
    $lbl = Join-Path $Dir "$($m.Name).LBL"
    $part = "$img.part"

    if (-not (Test-Path $lbl)) {
        Write-Host "== $($m.Name).LBL: скачиваю описание"
        if (-not ((Get-Url "$PDS/$($m.Pds).lbl" $lbl) -or (Get-Url "$MIT/$($m.Mit).LBL" $lbl))) {
            Write-Warning "не удалось скачать $($m.Name).LBL"; $failed = $true
        }
    }

    if ((Test-Path $img) -and (Get-Size $img) -eq $m.Bytes) {
        Write-Host "== $($m.Name).IMG: уже есть, проверяю контрольную сумму"
        if ((Get-Sha $img) -eq $m.Sha) { Write-Host '   в порядке'; continue }
        Write-Host '   сумма не совпала, скачиваю заново'
        Remove-Item -Force $img
    }

    Write-Host "== $($m.Name).IMG: скачиваю $([math]::Floor($m.Bytes / 1MB)) МБ"
    if ((Test-Path $part) -and (Get-Size $part) -gt $m.Bytes) { Remove-Item -Force $part }
    if (-not (Get-Url "$PDS/$($m.Pds).img" $part)) {
        Write-Host '   PDS недоступен, пробую imbrium.mit.edu'
        if (-not (Get-Url "$MIT/$($m.Mit).IMG" $part)) { Write-Warning "не удалось скачать $($m.Name).IMG"; $failed = $true; continue }
    }
    if ((Get-Size $part) -ne $m.Bytes -or (Get-Sha $part) -ne $m.Sha) {
        Write-Warning "$($m.Name).IMG: размер или контрольная сумма не совпали; удалите $part и запустите снова"
        $failed = $true; continue
    }
    Move-Item -Force $part $img
    Write-Host "   готово: $img"
}

if ($failed) { Write-Warning 'Часть карт не скачалась, см. сообщения выше.'; exit 1 }
Write-Host "Карты в $Dir готовы."

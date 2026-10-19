# Compiles the shader sources with the toolchains found on this machine (a manual check, not part of the build):
#   HLSL: dxc -T cs_6_0, every [numthreads] entry of examples/*.hlsl and tests/shaders/cores_check.hlsl
#   GLSL: glslangValidator -V, examples/*.comp and tests/shaders/cores_check.comp
# dxc and glslangValidator are taken from PATH, else from $env:VULKAN_SDK\Bin. A missing toolchain is reported and skipped.
# Exit code: the number of sources or entries that failed.

$root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$out = Join-Path ([System.IO.Path]::GetTempPath()) "qaws_shader_check"
New-Item -ItemType Directory -Force $out | Out-Null

function Find-Tool($name) {
    $cmd = Get-Command $name -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    if ($env:VULKAN_SDK) {
        $p = Join-Path $env:VULKAN_SDK "Bin\$name.exe"
        if (Test-Path $p) { return $p }
    }
    return $null
}

$failed = 0
$dxc = Find-Tool "dxc"
$glslang = Find-Tool "glslangValidator"

if ($dxc) {
    $hlsl = @(Get-ChildItem (Join-Path $root "examples\*.hlsl")) + @(Get-Item (Join-Path $PSScriptRoot "cores_check.hlsl"))
    foreach ($f in $hlsl) {
        $text = Get-Content $f.FullName -Raw
        $entries = [regex]::Matches($text, '\[numthreads[^\]]*\]\s*void\s+(\w+)') | ForEach-Object { $_.Groups[1].Value }
        foreach ($e in $entries) {
            $log = & $dxc -T cs_6_0 -E $e -Fo (Join-Path $out "$($f.BaseName)_$e.dxil") $f.FullName 2>&1
            if ($LASTEXITCODE -ne 0) { $failed++; Write-Output "FAIL hlsl $($f.Name):$e"; $log | Select-Object -First 5 | ForEach-Object { Write-Output "    $_" } }
            else { Write-Output "ok   hlsl $($f.Name):$e" }
        }
    }
} else { Write-Output "skip hlsl: dxc not found" }

if ($glslang) {
    $glsl = @(Get-ChildItem (Join-Path $root "examples\*.comp")) + @(Get-Item (Join-Path $PSScriptRoot "cores_check.comp"))
    foreach ($f in $glsl) {
        $log = & $glslang -V --target-env vulkan1.2 -o (Join-Path $out "$($f.BaseName).spv") $f.FullName 2>&1
        if ($LASTEXITCODE -ne 0) { $failed++; Write-Output "FAIL glsl $($f.Name)"; $log | Select-String "ERROR" | Select-Object -First 5 | ForEach-Object { Write-Output "    $_" } }
        else { Write-Output "ok   glsl $($f.Name)" }
    }
} else { Write-Output "skip glsl: glslangValidator not found" }

Write-Output "$failed failed"
exit $failed

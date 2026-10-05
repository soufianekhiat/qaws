# Converts a photo (jpg, png, bmp, ...) to a binary PPM (P6) the examples
# can read, plus a PNG of the same resized image for figure backgrounds.
#
#   powershell -File photo_to_ppm.ps1 -In photo.jpg -Out photo.ppm [-Width 480]
#
param(
    [Parameter(Mandatory = $true)][string]$In,
    [Parameter(Mandatory = $true)][string]$Out,
    [int]$Width = 480
)

Add-Type -AssemblyName System.Drawing

$src = [System.Drawing.Image]::FromFile((Resolve-Path $In))
$height = [int][math]::Round($src.Height * $Width / $src.Width)
$bmp = New-Object System.Drawing.Bitmap $Width, $height, ([System.Drawing.Imaging.PixelFormat]::Format24bppRgb)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::HighQualityBicubic
$g.DrawImage($src, 0, 0, $Width, $height)
$g.Dispose()
$src.Dispose()

$rect = New-Object System.Drawing.Rectangle 0, 0, $Width, $height
$data = $bmp.LockBits($rect, [System.Drawing.Imaging.ImageLockMode]::ReadOnly, $bmp.PixelFormat)
$stride = $data.Stride
$raw = New-Object byte[] ($stride * $height)
[System.Runtime.InteropServices.Marshal]::Copy($data.Scan0, $raw, 0, $raw.Length)
$bmp.UnlockBits($data)

$header = [System.Text.Encoding]::ASCII.GetBytes("P6`n$Width $height`n255`n")
$pixels = New-Object byte[] ($Width * $height * 3)
for ($y = 0; $y -lt $height; $y++) {
    for ($x = 0; $x -lt $Width; $x++) {
        $s = $y * $stride + $x * 3
        $d = ($y * $Width + $x) * 3
        $pixels[$d] = $raw[$s + 2]      # BGR -> RGB
        $pixels[$d + 1] = $raw[$s + 1]
        $pixels[$d + 2] = $raw[$s]
    }
}
$fs = [System.IO.File]::Create($Out)
$fs.Write($header, 0, $header.Length)
$fs.Write($pixels, 0, $pixels.Length)
$fs.Close()

$bmp.Save([System.IO.Path]::ChangeExtension($Out, ".png"), [System.Drawing.Imaging.ImageFormat]::Png)
$bmp.Dispose()
Write-Output "$Out ($Width x $height)"

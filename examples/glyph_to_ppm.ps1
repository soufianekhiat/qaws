# Renders text with an installed font to a binary PPM (P6), black on white,
# antialiased, for the font outline fitting application.
#
#   powershell -File glyph_to_ppm.ps1 -Text "&" -Out glyph.ppm [-Font Georgia] [-Size 400] [-Style Regular]
#
# -Size is the image height in pixels; the glyph em size is 70 % of it.
#
param(
    [Parameter(Mandatory = $true)][string]$Text,
    [Parameter(Mandatory = $true)][string]$Out,
    [string]$Font = "Georgia",
    [int]$Size = 400,
    [string]$Style = "Regular"
)

Add-Type -AssemblyName System.Drawing

$fontStyle = [System.Drawing.FontStyle]::$Style
$glyphFont = New-Object System.Drawing.Font $Font, ([single]($Size * 0.7)), $fontStyle, ([System.Drawing.GraphicsUnit]::Pixel)
$probe = New-Object System.Drawing.Bitmap 1, 1
$pg = [System.Drawing.Graphics]::FromImage($probe)
$format = [System.Drawing.StringFormat]::GenericTypographic
$extent = $pg.MeasureString($Text, $glyphFont, 100000, $format)
$pg.Dispose()
$probe.Dispose()

$width = [int][math]::Ceiling($extent.Width + 0.3 * $Size)
$height = $Size
$bmp = New-Object System.Drawing.Bitmap $width, $height, ([System.Drawing.Imaging.PixelFormat]::Format24bppRgb)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.Clear([System.Drawing.Color]::White)
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::AntiAliasGridFit
$g.DrawString($Text, $glyphFont, [System.Drawing.Brushes]::Black, [single](0.15 * $Size), [single](($height - $extent.Height) / 2), $format)
$g.Dispose()
$glyphFont.Dispose()

$rect = New-Object System.Drawing.Rectangle 0, 0, $width, $height
$data = $bmp.LockBits($rect, [System.Drawing.Imaging.ImageLockMode]::ReadOnly, $bmp.PixelFormat)
$stride = $data.Stride
$raw = New-Object byte[] ($stride * $height)
[System.Runtime.InteropServices.Marshal]::Copy($data.Scan0, $raw, 0, $raw.Length)
$bmp.UnlockBits($data)
$bmp.Dispose()

$header = [System.Text.Encoding]::ASCII.GetBytes("P6`n$width $height`n255`n")
$pixels = New-Object byte[] ($width * $height * 3)
for ($y = 0; $y -lt $height; $y++) {
    for ($x = 0; $x -lt $width; $x++) {
        $s = $y * $stride + $x * 3
        $d = ($y * $width + $x) * 3
        $pixels[$d] = $raw[$s + 2]
        $pixels[$d + 1] = $raw[$s + 1]
        $pixels[$d + 2] = $raw[$s]
    }
}
$stream = [System.IO.File]::Create($Out)
$stream.Write($header, 0, $header.Length)
$stream.Write($pixels, 0, $pixels.Length)
$stream.Close()
Write-Output "$Out : $width x $height"

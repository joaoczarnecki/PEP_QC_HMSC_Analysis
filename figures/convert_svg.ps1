# Caminho do executável do Inkscape
$ink = "C:\Program Files\Inkscape\bin\inkscape.exe"

# Diretório de saída para os JPG
$outDir = "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\jpg"

# Cria o diretório se não existir
if (!(Test-Path $outDir)) {
    New-Item -ItemType Directory -Path $outDir | Out-Null
}

# Lista de arquivos SVG
$files = @(
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\variance_partitioning_ordered_colorblind.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\gradient_plot_slope.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\species_network_plot.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\species_associations_ggcorrplot.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\gradient_plot_twi.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\gradient_plot_rtp.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\gradient_plot_pH.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\gradient_plot_OM.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\gradient_plot_Clay.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\gradient_plot_MAT.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\gradient_plot_MAP.svg",
  "F:\Thesis\3rdChapter\PEP_QC\results\undisturbed\hmsc_figures\beta_coefficients_heatmap.svg"
)

Write-Host "==============================="
Write-Host " CONVERTENDO SVG → JPEG..."
Write-Host "==============================="

foreach ($f in $files) {
    Write-Host "Convertendo: $f"
    $name = [System.IO.Path]::GetFileNameWithoutExtension($f)
    $outFile = Join-Path $outDir ($name + ".jpg")

    & $ink $f --export-type=jpg --export-filename="$outFile"
}

Write-Host ""
Write-Host "Conversão concluída."
Write-Host "Arquivos salvos em: $outDir"

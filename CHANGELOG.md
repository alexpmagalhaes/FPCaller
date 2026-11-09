# Changelog

## 1.0.0 (2026-11-09)

First release as a Python package. It is based on the `FPCaller.py` script used for the NuFANCI study.

### Added
- Installable package (`pip install .`) with an `fpcaller` command and a Python API (`run_fpcaller`).
- `--legacy` flag that reproduces the original script byte for byte. This was checked on the full-genome HCT116 NPM1 Ab-sort replicate 1 and on the chr19 test data.
- `--genome chm13|hg38|none` to choose the centromere filter.
- Options that were fixed in the script: `--chrom-regex`, `--poly-degree`, `--smooth-bins`, `--min-length`, `--plot-chrom`, `--outdir`.
- `--write-signal`, which writes per-bin z-score and corrected-ratio bedGraphs.
- Summary statistics are written as a TSV (`<sample>.stats.tsv`).
- The QC plot shows called domains and the excluded centromere, uses editable text (PDF fonttype 42) and Arial, and is saved as SVG, PNG and PDF.
- Test suite with HCT116 chr19 test data.

### Fixed (default mode; `--legacy` keeps the original behaviour)
- **First bin dropped:** the first line of each bedGraph was always skipped as if it were a header. Only real header lines are skipped now.
- **Rolling mean across chromosomes:** the 5-bin rolling mean ran over all chromosomes concatenated, so the first bins of each chromosome were averaged with the end of the previous one. It now runs per chromosome.
- **Domains spanning two chromosomes:** a run of positive bins was not split at a chromosome change. When the last bins of one chromosome and the first bins of the next were both positive, the resulting "domain" had an end coordinate on the next chromosome and was then lost by the size filter. Domains reaching chromosome ends are now called correctly.
- **hg38 centromeres:** the table was a Python dict with repeated keys, so only the last interval per chromosome was used. All intervals are now checked.
- **Swapped labels:** the start-up message printed the NuFANCI and input file names swapped.
- **Plot legend labels:** the input track and the log2 ratio were mislabelled.

### Effect of the fixes on real data
HCT116 NPM1 Ab sort replicate 1, full genome: 3,747 of 3,764 legacy domains are unchanged. Default mode adds 10 chromosome-end domains (about 18 Mb in total) and shifts a few chromosome-start and borderline domains by one bin. Default mode calls 3,778 domains in total.

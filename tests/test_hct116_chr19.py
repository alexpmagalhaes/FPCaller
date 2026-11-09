"""End-to-end tests on HCT116 chr19 test data.

* default mode must match the stored expected output;
* legacy mode must reproduce the original FPCaller.py output byte-for-byte
  (``expected/HCT116_chr19.original_FPCaller.bed`` was produced by the
  original script on the same files).
"""

import pandas as pd

from fpcaller import run_fpcaller
from fpcaller.cli import main


def read_bed(path):
    return pd.read_csv(path, sep="\t", header=None, names=["Chromosome", "Start", "End"])


def test_default_mode_matches_expected(chr19_fancy, chr19_input, expected):
    result = run_fpcaller(chr19_fancy, chr19_input)
    want = read_bed(expected / "HCT116_chr19.bed")
    pd.testing.assert_frame_equal(result.domains.reset_index(drop=True), want,
                                  check_dtype=False)


def test_legacy_mode_reproduces_original_script(chr19_fancy, chr19_input, expected, tmp_path):
    main(["-f", chr19_fancy, "-i", chr19_input, "-s", "legacy",
          "-o", str(tmp_path), "--legacy", "--plot-chrom", "none"])
    got = (tmp_path / "legacy.bed").read_bytes()
    want = (expected / "HCT116_chr19.original_FPCaller.bed").read_bytes()
    assert got == want


def test_cli_writes_all_outputs(chr19_fancy, chr19_input, expected, tmp_path):
    rc = main(["-f", chr19_fancy, "-i", chr19_input, "-s", "HCT116_chr19",
               "-o", str(tmp_path), "--plot-chrom", "chr19", "--write-signal"])
    assert rc == 0
    for suffix in [".bed", ".stats.tsv", ".zscore.bedgraph",
                   ".corrected_log2ratio.bedgraph",
                   ".chr19_qc.png", ".chr19_qc.svg", ".chr19_qc.pdf"]:
        assert (tmp_path / f"HCT116_chr19{suffix}").exists(), suffix
    assert (tmp_path / "HCT116_chr19.stats.tsv").read_text() == \
        (expected / "HCT116_chr19.stats.tsv").read_text()


def test_no_domain_overlaps_chm13_centromere(chr19_fancy, chr19_input):
    d = run_fpcaller(chr19_fancy, chr19_input).domains
    cen_start, cen_end = 19_655_575, 34_768_171       # CHM13 chr19
    assert not ((d["End"] >= cen_start) & (d["Start"] <= cen_end)).any()

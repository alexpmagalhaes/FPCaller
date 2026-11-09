"""Unit tests on small synthetic data (fast, no test files needed)."""

import gzip

import numpy as np
import pandas as pd
import pytest

from fpcaller import (call_domains, domain_size_stats, filter_centromeres,
                      mean_inter_domain_distance, read_bedgraph, zscore)


def bins(chrom, z_values, size=1000):
    """Build a per-bin table with the given z-scores on one chromosome."""
    n = len(z_values)
    return pd.DataFrame({"Chromosome": chrom,
                         "Start": np.arange(n) * size,
                         "End": np.arange(1, n + 1) * size,
                         "zscore": z_values})


# --- io -------------------------------------------------------------------
def test_read_bedgraph_keeps_first_bin_and_skips_headers(tmp_path):
    path = tmp_path / "x.bedgraph.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("track type=bedGraph\nchr1\t0\t1000\t1.5\nchr1\t1000\t2000\t2\n"
                 "chrX\t0\t1000\t3\n")
    df = read_bedgraph(str(path), "FANCY")
    assert list(df["Start"]) == [0, 1000]          # chrX dropped, first bin kept
    legacy = read_bedgraph(str(path), "FANCY", legacy=True)
    assert len(legacy) == 2                        # header line skipped either way


# --- domain calling -------------------------------------------------------
def test_call_domains_merges_positive_runs_and_applies_min_length():
    z = [-1, 1, 1, 1, -1, 1, 1, -1, 2, 2, 2, 2]
    d = call_domains(bins("chr1", z))
    # run of 2 bins (2 kb) is dropped; runs of 3 and 4 bins are kept
    assert d.values.tolist() == [["chr1", 1000, 4000], ["chr1", 8000, 12000]]


def test_call_domains_nan_is_negative():
    z = [1, 1, 1, np.nan, 1, 1, 1]
    assert len(call_domains(bins("chr1", z))) == 2


def test_call_domains_splits_at_chromosome_change():
    df = pd.concat([bins("chr1", [-1, 1, 1, 1]), bins("chr2", [1, 1, 1, -1])],
                   ignore_index=True)
    d = call_domains(df)
    assert d.values.tolist() == [["chr1", 1000, 4000], ["chr2", 0, 3000]]
    # legacy mode joins them into one run labelled with the first chromosome;
    # End < Start, so the size filter removes it (original behaviour)
    assert call_domains(df, legacy=True).empty


# --- z-score --------------------------------------------------------------
def test_zscore_uses_population_sd():
    df = pd.DataFrame({"RM_CorrectedLog2CPMRatio": [1.0, 2.0, 3.0, 4.0]})
    z = zscore(df)["zscore"].to_numpy()
    v = df["RM_CorrectedLog2CPMRatio"].to_numpy()
    assert np.allclose(z, (v - v.mean()) / v.std(ddof=0))


# --- centromeres ----------------------------------------------------------
def test_filter_centromeres_chm13():
    d = pd.DataFrame({"Chromosome": ["chr1", "chr1", "chr13", "chrUn"],
                      "Start": [1000, 120_000_000, 5_000_000, 0],
                      "End": [50_000, 121_000_000, 6_000_000, 10]})
    kept = filter_centromeres(d, "chm13")
    assert kept["Start"].tolist() == [1000, 0]     # chr1 cen and chr13 p-arm removed


def test_filter_centromeres_hg38_checks_all_intervals():
    # chr5 has several hg38 centromere intervals; this hits one that is not the last
    d = pd.DataFrame({"Chromosome": ["chr5"], "Start": [47_110_000], "End": [47_120_000]})
    assert filter_centromeres(d, "hg38").empty


def test_unknown_genome_raises():
    with pytest.raises(ValueError):
        filter_centromeres(pd.DataFrame({"Chromosome": ["chr1"], "Start": [0], "End": [1]}), "mm10")


# --- statistics -----------------------------------------------------------
def test_size_stats_and_distance():
    d = pd.DataFrame({"Chromosome": ["chr1"] * 3,
                      "Start": [0, 10_000, 100_000], "End": [6_000, 13_000, 101_000]})
    s = domain_size_stats(d)
    assert s["n_domains"] == 3 and s["total_bp"] == 10_000
    assert s["N50"] == 6_000 and s["N90"] == 3_000
    assert mean_inter_domain_distance(d) == (4_000 + 87_000) / 2

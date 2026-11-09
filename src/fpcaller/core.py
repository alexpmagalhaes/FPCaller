"""Core FPCaller algorithm.

Pipeline (one NuFANCI sample against its matched input):

1. ``log2_cpm_ratio``    per-chromosome log2 CPM ratio NuFANCI / input
2. ``baseline_correct``  remove broad trends per chromosome with a polynomial
                         fit, then smooth with a rolling mean
3. ``zscore``            genome-wide z-score of the smoothed, corrected ratio
4. ``call_domains``      contiguous runs of bins with z >= 0, longer than a
                         minimum length
5. ``filter_centromeres`` (see centromeres.py) drop domains overlapping
                         centromeres

``run_fpcaller`` chains all steps.

Legacy mode
-----------
``legacy=True`` reproduces the original FPCaller.py used for the NuFANCI paper
exactly. It differs from the default in three ways that only affect bins at
chromosome edges:

* the first line of each bedGraph is always skipped (see io.read_bedgraph);
* the rolling mean runs across chromosome boundaries (the last bins of one
  chromosome are averaged into the first bins of the next);
* a run of positive bins is not split when the chromosome changes.

The default mode fixes these three behaviours.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .io import COORDS, DEFAULT_CHROM_REGEX, natural_chrom_key, read_bedgraph
from .centromeres import filter_centromeres
from .stats import domain_size_stats, mean_inter_domain_distance


# ---------------------------------------------------------------------------
# Step 1: log2 CPM ratio
# ---------------------------------------------------------------------------
def log2_cpm_ratio(fancy, inp, legacy=False):
    """Per-chromosome log2 CPM ratio of NuFANCI over input.

    For every chromosome separately::

        log2( (FANCY + 1) / sum(FANCY + 1) ) - log2( (INPUT + 1) / sum(INPUT + 1) )

    A pseudocount of 1 avoids log(0). Normalising within each chromosome
    removes differences in chromosome copy number between samples.

    Parameters
    ----------
    fancy, inp : DataFrame
        Output of ``read_bedgraph`` with value columns ``FANCY`` (NuFANCI sample) and ``INPUT``.
        Only bins present in both files are kept (inner join on coordinates).
    legacy : bool
        Use Python's built-in ``sum`` (as the original script did) so floating
        point results are bit-identical to the original.

    Returns
    -------
    DataFrame with Chromosome, Start, End, FANCY, INPUT, log2CPMRatio.
    Chromosomes are in lexicographic order (as in the original script);
    bins keep their file order within each chromosome.
    """
    merged = fancy.merge(inp, on=COORDS)
    total = sum if legacy else np.sum

    parts = []
    # groupby(sort=True) gives lexicographic chromosome order, like the original.
    for _, group in merged.groupby("Chromosome", sort=True):
        group = group.copy()
        f = group["FANCY"] + 1
        i = group["INPUT"] + 1
        ratio = np.log2(f / total(f)) - np.log2(i / total(i))
        group["log2CPMRatio"] = ratio.fillna(0)
        parts.append(group)

    return pd.concat(parts, ignore_index=True)


# ---------------------------------------------------------------------------
# Step 2: baseline correction and smoothing
# ---------------------------------------------------------------------------
def baseline_correct(df, degree=3, smooth_bins=5, legacy=False):
    """Remove chromosome-wide trends from the log2 ratio, then smooth.

    For each chromosome a polynomial of ``degree`` is fitted to
    log2CPMRatio as a function of bin start position and subtracted
    (column ``CorrectedLog2CPMRatio``). This removes GC/mappability and
    copy-number trends that span whole chromosome arms.

    The corrected ratio is then smoothed with a trailing rolling mean over
    ``smooth_bins`` bins (column ``RM_CorrectedLog2CPMRatio``).

    legacy=True applies the rolling mean to all chromosomes concatenated
    (original behaviour); otherwise it is applied per chromosome.
    """
    required = COORDS + ["log2CPMRatio"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")

    parts = []
    for chrom in df["Chromosome"].unique():
        chrom_df = df[df["Chromosome"] == chrom].copy()
        x = chrom_df["Start"].values
        y = chrom_df["log2CPMRatio"].values

        # Fit and subtract the polynomial baseline.
        coefs = np.polyfit(x, y, degree)
        chrom_df["CorrectedLog2CPMRatio"] = y - np.polyval(coefs, x)

        if not legacy:
            chrom_df["RM_CorrectedLog2CPMRatio"] = (
                chrom_df["CorrectedLog2CPMRatio"]
                .rolling(window=smooth_bins, min_periods=1).mean()
            )
        parts.append(chrom_df)

    out = pd.concat(parts, ignore_index=True)

    if legacy:
        out["RM_CorrectedLog2CPMRatio"] = (
            out["CorrectedLog2CPMRatio"]
            .rolling(window=smooth_bins, min_periods=1).mean()
        )
    return out


# ---------------------------------------------------------------------------
# Step 3: z-score
# ---------------------------------------------------------------------------
def zscore(df, column="RM_CorrectedLog2CPMRatio"):
    """Genome-wide z-score of ``column`` (added as column ``zscore``).

    Uses the population standard deviation (ddof=0), computed as in the
    original script so values are identical.
    """
    if column not in df.columns:
        raise ValueError(f"The DataFrame must contain a '{column}' column.")

    values = df[column]
    n = len(values)
    if n > 1:
        # Equal to np.std(values, ddof=0); kept in this form for exact
        # numerical agreement with the original implementation.
        pop_sd = np.std(values, ddof=1) * np.sqrt(n - 1) / np.sqrt(n)
    else:
        pop_sd = 0
    pop_mean = np.mean(values)

    out = df.copy()
    out["zscore"] = (values - pop_mean) / pop_sd
    return out


# ---------------------------------------------------------------------------
# Step 4: domain calling
# ---------------------------------------------------------------------------
def call_domains(df, min_length=2001, legacy=False):
    """Merge contiguous bins with z-score >= 0 into domains.

    A domain starts at the first positive bin and ends at the end of the
    last consecutive positive bin. Bins with a NaN z-score count as negative.
    Domains of ``min_length`` bp or shorter are discarded (default 2001,
    i.e. at 1 kb resolution a domain needs at least 3 bins).

    legacy=True does not split a run of positive bins at a chromosome
    change (original behaviour). Such a run gets the chromosome of its
    first bin.

    Returns
    -------
    DataFrame with Chromosome, Start, End (one row per domain).
    """
    if df.empty:
        return pd.DataFrame(columns=COORDS)

    positive = (df["zscore"] >= 0).to_numpy()   # NaN >= 0 is False
    chrom = df["Chromosome"].to_numpy()

    # A new run starts at a positive bin whose previous bin is negative
    # (or, outside legacy mode, lies on another chromosome).
    prev_positive = np.r_[False, positive[:-1]]
    starts = positive & ~prev_positive
    if not legacy:
        chrom_changed = np.r_[True, chrom[1:] != chrom[:-1]]
        starts |= positive & chrom_changed

    # Label every positive bin with the index of its run.
    run_id = np.cumsum(starts)
    run_id = np.where(positive, run_id, 0)

    pos_df = df.loc[positive, COORDS].copy()
    pos_df["run"] = run_id[positive]

    domains = pos_df.groupby("run", sort=True).agg(
        Chromosome=("Chromosome", "first"),
        Start=("Start", "first"),
        End=("End", "last"),
    ).reset_index(drop=True)

    size = domains["End"] - domains["Start"]
    return domains[size > min_length].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------------------
@dataclass
class FPCallerResult:
    """Everything produced by one FPCaller run."""
    domains: pd.DataFrame   # final domains (Chromosome, Start, End)
    signal: pd.DataFrame    # per-bin table with all intermediate columns
    stats: dict             # summary statistics of the final domains


def run_fpcaller(fancy_path, input_path, genome="chm13",
                 chrom_regex=DEFAULT_CHROM_REGEX, degree=3, smooth_bins=5,
                 min_length=2001, legacy=False):
    """Run the complete FPCaller pipeline.

    Parameters
    ----------
    fancy_path, input_path : str
        Binned coverage bedGraphs of the NuFANCI sample and its input
        (same binning, e.g. deepTools ``bamCoverage`` output).
    genome : {"chm13", "hg38", None}
        Centromere set used to remove domains; ``None`` disables the filter.
    chrom_regex : str or None
        Chromosomes to analyse (default: autosomes ``chr1``-``chrN``).
    degree : int
        Degree of the per-chromosome baseline polynomial.
    smooth_bins : int
        Window (in bins) of the rolling mean.
    min_length : int
        Domains must be longer than this (bp).
    legacy : bool
        Reproduce the original FPCaller.py exactly (see module docstring).

    Returns
    -------
    FPCallerResult
    """
    fancy = read_bedgraph(fancy_path, "FANCY", chrom_regex, legacy=legacy)
    inp = read_bedgraph(input_path, "INPUT", chrom_regex, legacy=legacy)
    if fancy.empty or inp.empty:
        raise ValueError("No bins left after reading/filtering the bedGraphs. "
                         "Check the files and --chrom-regex.")

    signal = log2_cpm_ratio(fancy, inp, legacy=legacy)
    signal = baseline_correct(signal, degree=degree, smooth_bins=smooth_bins,
                              legacy=legacy)
    signal = zscore(signal)

    domains = call_domains(signal, min_length=min_length, legacy=legacy)
    if genome is not None:
        domains = filter_centromeres(domains, genome)

    if not legacy:
        # Natural chromosome order (chr1, chr2, ..., chr10) for output.
        key = domains["Chromosome"].map(natural_chrom_key)
        domains = (domains.assign(_k=key.map(tuple))
                   .sort_values(["_k", "Start"]).drop(columns="_k")
                   .reset_index(drop=True))

    stats = domain_size_stats(domains)
    stats["mean_inter_domain_distance"] = mean_inter_domain_distance(domains)

    return FPCallerResult(domains=domains, signal=signal, stats=stats)

"""Reading and writing of bedGraph / BED files."""

import re

import pandas as pd

# Column names used throughout the package.
COORDS = ["Chromosome", "Start", "End"]

# Default chromosome filter: autosomes only (chr1, chr2, ...).
# Sex chromosomes, chrM, unplaced contigs and the rDNA contig are dropped.
DEFAULT_CHROM_REGEX = r"^chr\d+$"


def read_bedgraph(path, value_name, chrom_regex=DEFAULT_CHROM_REGEX, legacy=False):
    """Read a 4-column bedGraph (chrom, start, end, value) into a DataFrame.

    Parameters
    ----------
    path : str
        bedGraph file, plain text or gzip-compressed (``.gz``).
    value_name : str
        Name for the signal column, e.g. ``"FANCY"`` or ``"INPUT"``.
    chrom_regex : str or None
        Keep only chromosomes whose name matches this regular expression.
        ``None`` keeps all chromosomes.
    legacy : bool
        If True, reproduce the original FPCaller behaviour of skipping the
        first line of the file unconditionally (it assumed a header line).
        If False, only real header lines (``track``, ``browser``, ``#``) are
        skipped, so the first data bin is kept.

    Returns
    -------
    pandas.DataFrame with columns Chromosome, Start, End, <value_name>.
    """
    names = COORDS + [value_name]

    if legacy:
        # Original behaviour: always drop line 1.
        df = pd.read_csv(path, sep="\t", skiprows=1, names=names, low_memory=False)
    else:
        df = pd.read_csv(path, sep="\t", header=None, names=names,
                         comment="#", low_memory=False)
        # Drop UCSC header lines ("track ..." / "browser ...") if present.
        is_header = df["Chromosome"].astype(str).str.match(r"^(track|browser)\b")
        df = df[~is_header]

    df = df.astype({"Chromosome": str, "Start": "int64", "End": "int64",
                    value_name: "float64"})

    if chrom_regex is not None:
        df = df[df["Chromosome"].str.match(chrom_regex)]

    return df.reset_index(drop=True)


def write_bed(domains, path):
    """Write domains (Chromosome, Start, End) as a headerless 3-column BED file."""
    domains[COORDS].to_csv(path, sep="\t", header=False, index=False)


def write_signal_bedgraph(signal, column, path):
    """Write one per-bin signal column (e.g. ``zscore``) as a bedGraph track."""
    signal[COORDS + [column]].to_csv(path, sep="\t", header=False, index=False)


def natural_chrom_key(chrom):
    """Sort key giving chr1, chr2, ..., chr10 instead of chr1, chr10, chr2."""
    parts = re.split(r"(\d+)", chrom)
    return [int(p) if p.isdigit() else p for p in parts]

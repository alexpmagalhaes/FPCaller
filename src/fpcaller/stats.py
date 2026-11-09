"""Summary statistics of called domains."""

import numpy as np


def domain_size_stats(domains):
    """Number, total length, mean size and N10/N50/N90 of domain sizes.

    Nx is the size of the domain at which the cumulative length of domains,
    sorted from largest to smallest, first reaches x % of the total length.
    For example, half of all domain base pairs lie in domains >= N50.
    """
    sizes = (domains["End"] - domains["Start"]).to_numpy()
    if sizes.size == 0:
        return {"n_domains": 0, "total_bp": 0, "mean_size": 0,
                "N10": 0, "N50": 0, "N90": 0}

    sorted_sizes = np.sort(sizes)[::-1]          # largest first
    cumsum = np.cumsum(sorted_sizes)
    total = cumsum[-1]

    def nx(fraction):
        # first domain whose cumulative sum reaches the threshold
        return int(sorted_sizes[np.argmax(cumsum >= fraction * total)])

    return {
        "n_domains": int(sizes.size),
        "total_bp": int(total),
        "mean_size": float(sizes.mean()),
        "N10": nx(0.10),
        "N50": nx(0.50),
        "N90": nx(0.90),
    }


def mean_inter_domain_distance(domains):
    """Mean gap (bp) between consecutive domains on the same chromosome.

    Gap = start of a domain minus end of the previous domain. Returns 0 when
    no chromosome has two or more domains.
    """
    if len(domains) < 2:
        return 0.0
    d = domains.sort_values(["Chromosome", "Start"])
    same_chrom = d["Chromosome"].to_numpy()[1:] == d["Chromosome"].to_numpy()[:-1]
    gaps = d["Start"].to_numpy()[1:] - d["End"].to_numpy()[:-1]
    gaps = gaps[same_chrom]
    return float(gaps.mean()) if gaps.size else 0.0

"""FPCaller: call FANCY-enriched chromatin domains from binned DNA-seq coverage.

FPCaller compares a FANCY-sorted DNA-seq sample with its matched input and
reports broad enriched domains, for example nucleolus-associated domains (NADs)
or speckle-associated domains (SADs).

Typical use is through the command line (``fpcaller --help``), but every step
is also available as a Python function::

    from fpcaller import run_fpcaller
    result = run_fpcaller("fancy.bedgraph", "input.bedgraph")
    result.domains.head()
"""

__version__ = "1.0.0"

from .io import read_bedgraph, write_bed, write_signal_bedgraph
from .core import (
    log2_cpm_ratio,
    baseline_correct,
    zscore,
    call_domains,
    run_fpcaller,
    FPCallerResult,
)
from .centromeres import filter_centromeres, available_genomes
from .stats import domain_size_stats, mean_inter_domain_distance

__all__ = [
    "__version__",
    "read_bedgraph",
    "write_bed",
    "write_signal_bedgraph",
    "log2_cpm_ratio",
    "baseline_correct",
    "zscore",
    "call_domains",
    "run_fpcaller",
    "FPCallerResult",
    "filter_centromeres",
    "available_genomes",
    "domain_size_stats",
    "mean_inter_domain_distance",
]

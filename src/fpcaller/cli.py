"""Command-line interface: ``fpcaller --fancy F.bedgraph --input I.bedgraph --sample NAME``."""

import argparse
import os
import sys

from . import __version__
from .centromeres import available_genomes
from .core import run_fpcaller
from .io import DEFAULT_CHROM_REGEX, write_bed, write_signal_bedgraph


def build_parser():
    p = argparse.ArgumentParser(
        prog="fpcaller",
        description="Call FANCY-enriched domains (e.g. NADs, SADs) from binned "
                    "DNA-seq coverage of a FANCY sample and its matched input.",
    )
    req = p.add_argument_group("required")
    req.add_argument("-f", "--fancy", required=True,
                     help="bedGraph of the FANCY sample (plain or .gz)")
    req.add_argument("-i", "--input", required=True,
                     help="bedGraph of the matched input (same binning)")
    req.add_argument("-s", "--sample", required=True,
                     help="sample name, used as output file prefix")

    p.add_argument("-o", "--outdir", default=".",
                   help="output directory (default: current directory)")
    p.add_argument("-g", "--genome", default="chm13",
                   choices=available_genomes() + ["none"],
                   help="centromere set for filtering domains; 'none' disables "
                        "the filter (default: chm13)")
    p.add_argument("--chrom-regex", default=DEFAULT_CHROM_REGEX,
                   help=f"chromosomes to analyse (default: '{DEFAULT_CHROM_REGEX}', autosomes)")
    p.add_argument("--poly-degree", type=int, default=3,
                   help="degree of the per-chromosome baseline polynomial (default: 3)")
    p.add_argument("--smooth-bins", type=int, default=5,
                   help="rolling-mean window in bins (default: 5)")
    p.add_argument("--min-length", type=int, default=2001,
                   help="keep domains longer than this many bp (default: 2001)")
    p.add_argument("--plot-chrom", default="chr14",
                   help="chromosome for the QC plot, or 'none' (default: chr14)")
    p.add_argument("--write-signal", action="store_true",
                   help="also write per-bin z-score and corrected ratio bedGraphs")
    p.add_argument("--legacy", action="store_true",
                   help="reproduce the original FPCaller.py exactly "
                        "(see README, 'Legacy mode')")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, args.sample)

    print(f"FPCaller {__version__}{' (legacy mode)' if args.legacy else ''}")
    print(f"  FANCY : {args.fancy}")
    print(f"  Input : {args.input}")
    print(f"  Sample: {args.sample}")

    genome = None if args.genome == "none" else args.genome
    result = run_fpcaller(
        args.fancy, args.input, genome=genome, chrom_regex=args.chrom_regex,
        degree=args.poly_degree, smooth_bins=args.smooth_bins,
        min_length=args.min_length, legacy=args.legacy,
    )

    # Domains
    write_bed(result.domains, f"{prefix}.bed")

    # Summary statistics (key <tab> value)
    with open(f"{prefix}.stats.tsv", "w") as fh:
        for key, value in result.stats.items():
            fh.write(f"{key}\t{value}\n")

    # Optional per-bin tracks
    if args.write_signal:
        write_signal_bedgraph(result.signal, "zscore", f"{prefix}.zscore.bedgraph")
        write_signal_bedgraph(result.signal, "RM_CorrectedLog2CPMRatio",
                              f"{prefix}.corrected_log2ratio.bedgraph")

    # QC plot
    if args.plot_chrom.lower() != "none":
        from .plotting import plot_chromosome  # imported here: matplotlib is slow to load
        written = plot_chromosome(result.signal, result.domains, args.plot_chrom,
                                  f"{prefix}.{args.plot_chrom}_qc", genome=genome)
        if not written:
            print(f"  Note: {args.plot_chrom} not in data, QC plot skipped.",
                  file=sys.stderr)

    s = result.stats
    print(f"  Domains: {s['n_domains']}  total {s['total_bp']:,} bp  N50 {s['N50']:,} bp")
    print(f"  Written: {prefix}.bed, {prefix}.stats.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

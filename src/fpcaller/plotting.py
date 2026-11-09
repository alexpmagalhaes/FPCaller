"""QC plot of one chromosome: raw signal, ratio, z-score and called domains."""

import matplotlib

matplotlib.use("Agg")  # no display needed; works on clusters
import matplotlib as mpl
import matplotlib.pyplot as plt

# Editable text in vector output (PDF/SVG) and Arial for all text.
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]


def plot_chromosome(signal, domains, chrom, out_prefix, genome=None,
                    formats=("svg", "png", "pdf")):
    """Plot NuFANCI/input coverage, log2 ratio, z-score and domains for one chromosome.

    Left axis: NuFANCI and input coverage.
    Right axis: log2 CPM ratio and z-score (zero line in red).
    Grey bands: called domains.
    Hatched band: centromere interval(s) excluded by the filter (if ``genome``).

    Writes ``<out_prefix>.<fmt>`` for each format and returns the list of
    written paths (empty if ``chrom`` is not in ``signal``).
    """
    sub = signal[signal["Chromosome"] == chrom]
    if sub.empty:
        return []

    x = sub["Start"] / 1e6  # Mb
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.plot(x, sub["FANCY"], color="#1f4e9c", lw=0.8, label="NuFANCI")
    ax.plot(x, sub["INPUT"], color="#e377c2", lw=0.8, label="Input")
    ax.set_xlabel(f"{chrom} position (Mb)")
    ax.set_ylabel("Coverage")

    ax2 = ax.twinx()
    ax2.plot(x, sub["log2CPMRatio"], color="#17becf", lw=0.6, label="log2 CPM ratio")
    ax2.plot(x, sub["zscore"], color="#ff7f0e", lw=0.8, label="z-score")
    ax2.axhline(0, color="red", lw=0.8)
    ax2.set_ylabel("log2 ratio / z-score")

    # Shade called domains.
    for _, d in domains[domains["Chromosome"] == chrom].iterrows():
        ax.axvspan(d["Start"] / 1e6, d["End"] / 1e6, color="grey", alpha=0.2, lw=0)

    # Mark centromere intervals removed by the filter.
    if genome is not None:
        from .centromeres import CENTROMERES
        first = True
        for c, s, e in CENTROMERES[genome]:
            if c == chrom:
                ax.axvspan(s / 1e6, e / 1e6, facecolor="none", edgecolor="#999999",
                           hatch="//", lw=0, label="centromere (excluded)" if first else None)
                first = False

    # One legend for both axes.
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", bbox_to_anchor=(1.07, 1.0),
              frameon=False, fontsize=8)
    ax.set_title(f"FPCaller QC: {chrom} (grey = called domains)")
    ax.set_xlim(x.min(), x.max())
    fig.tight_layout()

    paths = []
    for fmt in formats:
        path = f"{out_prefix}.{fmt}"
        fig.savefig(path, format=fmt, dpi=200 if fmt == "png" else None)
        paths.append(path)
    plt.close(fig)
    return paths

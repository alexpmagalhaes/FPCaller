"""Centromere coordinates and the centromere filter.

Domains overlapping a centromere are removed because centromeric satellite
arrays are poorly mappable and give unreliable NuFANCI/input ratios.

Coordinates are 0-based, half-open intervals (BED convention), taken from the
original FPCaller.py:

* ``chm13``: one interval per chromosome for T2T-CHM13 v2.0. For the
  acrocentric chromosomes (chr13, 14, 15, 21, 22) the interval runs from 0 to
  the end of the centromere, so the short arm (including the rDNA arrays) is
  excluded as well.
* ``hg38``: GRCh38 centromere model intervals (several per chromosome).

Add a new genome by adding an entry to ``CENTROMERES``.
"""

import numpy as np
import pandas as pd

# (chromosome, start, end)
CHM13_CENTROMERES = [
    ("chr1", 116796047, 147241659),
    ("chr2", 85991677, 99673021),
    ("chr3", 85804701, 101415026),
    ("chr4", 44705153, 59870510),
    ("chr5", 42077201, 54596624),
    ("chr6", 53286688, 66058390),
    ("chr7", 55414371, 68714499),
    ("chr8", 39243545, 51325080),
    ("chr9", 39952797, 81694041),
    ("chr10", 34633793, 46664589),
    ("chr11", 46061947, 59413484),
    ("chr12", 29620498, 42202490),
    ("chr13", 0, 22498291),
    ("chr14", 0, 17708240),
    ("chr15", 0, 22694129),
    ("chr16", 30848286, 57219471),
    ("chr17", 18892419, 32486939),
    ("chr18", 10965698, 25933550),
    ("chr19", 19655575, 34768171),
    ("chr20", 21383660, 37969538),
    ("chr21", 0, 16341849),
    ("chr22", 0, 20711065),
    ("chrX", 52820106, 65927025),
    ("chrY", 21838, 62453637),
]

HG38_CENTROMERES = [
    ("chr1", 122026459, 122224535),
    ("chr1", 122224635, 122503147),
    ("chr1", 122503247, 124785432),
    ("chr1", 124785532, 124849129),
    ("chr1", 124849229, 124932724),
    ("chr10", 39686682, 39935900),
    ("chr10", 39936000, 41497440),
    ("chr10", 41497540, 41545720),
    ("chr10", 41545820, 41593521),
    ("chr11", 51078348, 51090317),
    ("chr11", 51090417, 54342399),
    ("chr11", 54342499, 54425074),
    ("chr12", 34769407, 34816611),
    ("chr12", 34835295, 37185252),
    ("chr13", 16000000, 16022537),
    ("chr13", 16022637, 16110659),
    ("chr13", 16110759, 16164892),
    ("chr13", 16164992, 16228527),
    ("chr13", 16228627, 16249297),
    ("chr13", 16249397, 16256067),
    ("chr13", 16256167, 16259412),
    ("chr13", 16259512, 16282073),
    ("chr13", 16282173, 17416384),
    ("chr13", 17416484, 17416824),
    ("chr13", 17416924, 17417264),
    ("chr13", 17417364, 17418562),
    ("chr13", 17418662, 18051248),
    ("chr14", 15835778, 16307609),
    ("chr14", 16000000, 16022537),
    ("chr14", 16140627, 16228649),
    ("chr14", 16228749, 16282882),
    ("chr14", 16282982, 16346517),
    ("chr14", 16346617, 16367287),
    ("chr14", 16367387, 16374057),
    ("chr14", 16374157, 16377402),
    ("chr14", 16377502, 16400063),
    ("chr14", 16404448, 17538659),
    ("chr14", 17538759, 17539099),
    ("chr14", 17539199, 17539539),
    ("chr14", 17539639, 17540837),
    ("chr14", 17540937, 18173523),
    ("chr15", 17083673, 17498951),
    ("chr15", 17499051, 18355008),
    ("chr15", 18355108, 19725254),
    ("chr16", 36311158, 36334460),
    ("chr16", 36337666, 38265669),
    ("chr17", 22813679, 23194918),
    ("chr17", 23195018, 26566633),
    ("chr17", 26566733, 26616164),
    ("chr18", 15460899, 15780377),
    ("chr18", 15797855, 20561439),
    ("chr18", 20603247, 20696289),
    ("chr18", 20696389, 20736025),
    ("chr18", 20736125, 20813083),
    ("chr18", 20839797, 20861206),
    ("chr19", 24498980, 24552652),
    ("chr19", 24552752, 24891256),
    ("chr19", 24908689, 27190874),
    ("chr2", 92188145, 94090557),
    ("chr20", 26436232, 26586955),
    ("chr20", 26608145, 28494539),
    ("chr20", 28508997, 28556953),
    ("chr20", 28648108, 28728874),
    ("chr20", 29125793, 29204668),
    ("chr20", 29917404, 30038348),
    ("chr21", 10864560, 10887097),
    ("chr21", 10887197, 10975219),
    ("chr21", 10975319, 11029452),
    ("chr21", 11029552, 11093087),
    ("chr21", 11093187, 11113857),
    ("chr21", 11113957, 11120627),
    ("chr21", 11120727, 11123972),
    ("chr21", 11124072, 11146633),
    ("chr21", 11146733, 12280944),
    ("chr21", 12281044, 12281384),
    ("chr21", 12281484, 12281824),
    ("chr21", 12281924, 12283122),
    ("chr21", 12283222, 12915808),
    ("chr22", 12954788, 12977325),
    ("chr22", 13021422, 13109444),
    ("chr22", 13109544, 13163677),
    ("chr22", 13163777, 13227312),
    ("chr22", 13227412, 13248082),
    ("chr22", 13248182, 13254852),
    ("chr22", 13254952, 13258197),
    ("chr22", 13258297, 13280858),
    ("chr22", 13285243, 14419454),
    ("chr22", 14419554, 14419894),
    ("chr22", 14419994, 14420334),
    ("chr22", 14420434, 14421632),
    ("chr22", 14421732, 15054318),
    ("chr3", 90772458, 91233586),
    ("chr3", 91233686, 91247622),
    ("chr3", 91553419, 93655574),
    ("chr4", 49712061, 51743951),
    ("chr5", 46485900, 46569062),
    ("chr5", 46569162, 46796725),
    ("chr5", 46796825, 47061288),
    ("chr5", 47106994, 47153339),
    ("chr5", 47153439, 47296069),
    ("chr5", 47309184, 49591369),
    ("chr5", 49667531, 49721203),
    ("chr5", 49721303, 50059807),
    ("chr6", 58553888, 59829934),
    ("chr7", 58169653, 60828234),
    ("chr7", 61377788, 61528020),
    ("chr8", 44033744, 45877265),
    ("chr9", 43389635, 45518558),
    ("chrX", 58605579, 62412542),
    ("chrY", 10316944, 10544039),
]

CENTROMERES = {
    "chm13": CHM13_CENTROMERES,
    "hg38": HG38_CENTROMERES,
}


def available_genomes():
    """Names of genomes with built-in centromere coordinates."""
    return sorted(CENTROMERES)


def filter_centromeres(domains, genome="chm13"):
    """Remove domains that overlap any centromere interval of ``genome``.

    Overlap is tested as in the original script: a domain is kept only if it
    ends before the centromere starts or starts after it ends
    (``End < cen_start or Start > cen_end``). Chromosomes without a
    centromere entry are kept unchanged.

    Parameters
    ----------
    domains : DataFrame with Chromosome, Start, End.
    genome : str
        Key of ``CENTROMERES`` (``"chm13"`` or ``"hg38"``).

    Returns
    -------
    DataFrame with the overlapping domains removed.
    """
    if genome not in CENTROMERES:
        raise ValueError(f"Unknown genome '{genome}'. "
                         f"Choose from: {', '.join(available_genomes())}")
    if domains.empty:
        return domains

    drop = np.zeros(len(domains), dtype=bool)
    chrom = domains["Chromosome"].to_numpy()
    start = domains["Start"].to_numpy()
    end = domains["End"].to_numpy()

    for c, cen_start, cen_end in CENTROMERES[genome]:
        on_chrom = chrom == c
        overlaps = ~((end < cen_start) | (start > cen_end))
        drop |= on_chrom & overlaps

    return domains[~drop].reset_index(drop=True)

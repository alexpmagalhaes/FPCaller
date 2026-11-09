#!/usr/bin/env bash
# Run FPCaller on the bundled HCT116 chr19 test data.
# Usage (from the repository root, after `pip install .`):
#   bash examples/run_example.sh

DATA=tests/data
OUT=results/example

fpcaller \
  --fancy $DATA/HCT116_NPM1_Ab_sort_rep1.chr19.w50k_b1k.bedgraph.gz \
  --input $DATA/HCT116_NPM1_Ab_sort_rep1_input.chr19.w50k_b1k.bedgraph.gz \
  --sample HCT116_NPM1_Ab_sort_rep1_chr19 \
  --outdir $OUT \
  --plot-chrom chr19 \
  --write-signal

# Compare with the expected domains shipped with the tests.
diff <(cut -f1-3 $OUT/HCT116_NPM1_Ab_sort_rep1_chr19.bed) tests/data/expected/HCT116_chr19.bed \
  && echo "Output matches tests/data/expected/HCT116_chr19.bed"

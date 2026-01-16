E. coli Read Mapper (C++)

Implements:
- FM-index (BWT + Occ checkpoints) over the reference genome
- Seed-and-extend mapping
  - remove leading/trailing N (and any non-ACGT)
  - internal N splits the read into ACGT-only segments
  - seeds (~22bp) are searched exactly with the FM-index
  - candidate hits are extended with banded Smith-Waterman
- reverse-strand mapping by building a second FM-index over the reverse-complement reference

Build:
  mkdir -p build
  cmake -S . -B build
  cmake --build build -j

Run:
  ./build/ecoli_mapper \
    --ref GCF_000005845.2_ASM584v2_genomic.fna.gz \
    --reads ERR022075_1.fastq.gz \
    --out results/ecoli \
    --seed-len 22 \
    --min-seed 12 \
    --min-score 50 \
    --cov-depth results/ecoli.depth.tsv

Outputs:
- results/ecoli.alignments.tsv  (per-read mapping summary)
- results/ecoli.report.txt      (assignment report)
- results/ecoli.depth.tsv       (optional per-base depth)

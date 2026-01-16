## Mapping algorithms

### Reference indexing (BWT + FM-index)
The reference (and its reversed complement) is indexed with an FM-index over the Burrows–Wheeler Transform (BWT):

- Suffix array is built (doubling + counting/radix sort), then BWT is derived from it.
- FM-index stores `C` (prefix counts) and `Occ` (rank table) to support LF-mapping / backward search.
- The suffix array is used to convert FM hit intervals into reference coordinates.

Seeds are queried via standard FM backward search.

### k-mer seeding (with N handling)
Reads are split into fixed-length k-mers (e.g., `k=20` with stride, step = 10), skipping any seed containing `N`.

Each seed is tracked as `(offset_in_read, seed_string)`.

### Candidate generation + voting
For each seed:
1. FM backward search returns an SA interval `[l, r)`.
2. Each hit yields a candidate start `pos = SA[i] - offset`.
3. Very large hit sets are capped by uniform subsampling to limit runtime in repetitive regions.

All candidate starts are aggregated and **voted** (count duplicates). Candidates are ranked by vote count and only the top set is evaluated.

### Verification + Smith–Waterman fallback
Candidates are checked in two stages:

1. **Exact check (fast path):** compare read vs reference at `pos` with zero mismatches; accept immediately if exact. 
2. **SW fallback:** if not exact, run Smith–Waterman local alignment against a bounded window around `pos` and accept if the alignment is sufficiently long and high-identity.

### Strand handling and multi-mapping
The full procedure runs for both the read and its reverse complement.  
Reads are reported as **unique** if only one high-scoring locus remains, otherwise **multi-mapped**.

## Program output

```
Total reads: 22720100
Mapped reads: 22530184 (99.1641%)
Unique reads: 22018717 (96.9129% of total, 97.7299% of mapped)
Multi-mapped reads: 511467 (2.25117% of total, 2.27014% of mapped)
Unmapped reads: 189916
Total SW calls: 11031391
Reads with at least one verify-accepted alignment: 19387428
Reads with at least one SW-accepted alignment: 3509266
Mean Q (per read average): 36.1585
Mean coverage: 484.9
```

## Summary

As per what was asked in the task:

- **Algorithms used**
  - **BWT + FM-index** for fast exact seed search.
  - **k-mer seeding**.
  - **Candidate voting** across seed hits to rank likely mapping starts.
  - **Exact verification** at candidate start; **Smith–Waterman** local alignment as fallback.

- **Total reads processed:** 22,720,100

- **Mapping rate:** 22,530,184 / 22,720,100 = **99.1641%** mapped  
  (unmapped: **0.8359%**), matched exactly: 19,387,428 / 22,720,100 = **85.3%**

- **Unique / multi-mapped**
  - **Unique:** 22,018,717 (**96.9129%** of total; **97.7299%** of mapped)
  - **Multi-mapped:** 511,467 (**2.2512%** of total; **2.2701%** of mapped)

- **Alignment quality**
  - Mean per-read average base quality: **Q = 36.1585**

- **Genome read coverage:** mean **484.9**

## Interpretation

- **Exact verification accounts for the majority of successful mappings.** 19.39M reads have at least one verify-accepted alignment (~86% of mapped). This indicates that many reads match the reference perfectly at the chosen locus after candidate selection, and the fast path is effective.

- **Smith–Waterman materially increases sensitivity but is the main computational cost driver.** 3.51M reads have at least one SW-accepted alignment (~15.6% of mapped), and there are 11.0M SW calls total. This implies that the pipeline trades compute for sensitivity on reads with mismatches or imperfect seed consensus.

- **Read quality is very high, so errors are not the dominant limitation.** Mean per-read Phred Q ≈ 36 corresponds to an expected error rate around 2.4×10⁻⁴ per base, i.e., well below 0.1%.

- **Coverage is very deep .** Mean coverage ~485 implies the dataset provides dense sampling across the genome.

## How to run
build:
```
mkdir build
cmake -S . -B build
cmake --build build --config Release
```
run:
```
.\build\mapper.exe GCF_000005845.2_ASM584v2_genomic.fna ERR022075_1.fastq
```
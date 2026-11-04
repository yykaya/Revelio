# Changelog

## 0.4.0 — 2026-11-04

- **`revelio call`**: alignment-free 5mC calls from MM/ML-tagged long reads (PacBio HiFi, ONT). Calls are placed on any reference (locus, haplotype, small genome) via unique k-mer anchors and collinear chaining; CpG strands are pooled, CHG/CHH stay stranded. Unaligned BAMs work (`samtools view | revelio call`). Secondary/supplementary records, hard-clipped reads and `MN` mismatches are skipped.
- Stdin is read through a 1 MiB `fread` buffer (libc++ `std::cin` was up to ~20× slower on long lines); stdin read errors now exit non-zero.

## 0.3.0 — 2026-11-04

- **Window profiles**: `--window <bp>` (and `--step <bp>`) write one row per window with GC, CpG observed/expected, site counts per context and strand, and the mean prior per context. Counts equal the per-site rows of the same command.

## 0.2.0 — 2026-11-04

- **Both strands** are scanned by default (`--strand both|+|-`). Minus-strand sites are reported at the + strand coordinate of their G, with context and k-mer read 5'→3' on the minus strand. `--strand +` reproduces the 0.1.x output.
- Weights tables accept per-context intercepts `BIAS_CG`, `BIAS_CHG`, `BIAS_CHH`: `p = sigmoid(BIAS + BIAS_<context> + w[context,kmer])`.
- Worked example (`examples/`, `docs/example_sites.png`).
- ⚠ Tables with `BIAS_<context>` rows need Revelio ≥ 0.2.0. Older binaries read every row starting with `BIAS` as the global bias and give wrong scores without an error.

## 0.1.1 — 2026-11-04

- Input from pipes and stdin works (`-i -`, `-i /dev/stdin`); CRLF line endings are handled.
- `seq_id` is the header up to the first whitespace.
- Invalid options and inputs exit non-zero; `--input`, `--kmer`, `-o` aliases.
- Weights tables are validated on load; an unreadable `--model` is an error instead of a silent fallback to `--heuristic`.
- Header column `pos` renamed to `pos_1based`.
- `make test` runs a regression suite.

## 0.1.0 — 2025-08-19

- Initial release.

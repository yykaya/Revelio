# Changelog

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

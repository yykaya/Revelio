# Revelio 🪄 — Lightweight Methylation-Context Predictor

<table>
  <tr>
    <td width="260" valign="top" style="padding-right:16px;">
      <img src="revelio.png"
           alt="Revelio — Lightweight Methylation-Context Predictor"
           width="240">
    </td>
    <td valign="top">
      <strong>Revelio</strong> (HP charm for revealing hidden things) quickly scans
      FASTA/FASTQ to identify cytosine contexts <strong>CG / CHG / CHH</strong> and assigns a
      <strong>sequence-based prior probability</strong> to each candidate site (either a simple
      heuristic or a user-provided logistic <strong>k-mer</strong> weight table). It’s designed for
      <strong>fast triage</strong> and <strong>candidate region ranking</strong> before heavier
      signal-level methylation callers.
      <br><br>
      <blockquote>⚠️ Revelio uses <strong>sequence context only</strong> (no raw signal). For
      signal-level calling on ONT data, use tools such as Remora/Dorado, Megalodon, or DeepSignal-plant.</blockquote>
    </td>
  </tr>
</table>



**Status:** research prototype (v0.1.1) • **Developer:** Yasin Kaya, MPIPZ (ykaya@mpipz.mpg.de | yyasinkkaya@gmail.com) • **License:** MIT

---

## Features

- **Fast, portable, single binary.** Streams large genomes/reads in one pass.
- **Context-aware priors.** Pick a quick **heuristic** or supply your own **logistic k-mer weights**.
- **FASTA/FASTQ auto-detection.** Handles multi-line FASTA and 4-line FASTQ.
- **Practical TSV output**: `seq_id, pos_1based, context, strand, kmer, prob`.
- **Works from stdin.** Pipe compressed input via `zcat`/`pigz -dc` (`-i -` or `-i /dev/stdin`).

---

## Build

### Prerequisites
- A C++17 compiler (GCC ≥ 8, Clang ≥ 7, or MSVC ≥ 2019).

### Linux / macOS
```bash
g++ -std=c++17 -O3 -march=native -o revelio revelio.cpp
# or
clang++ -std=c++17 -O3 -march=native -o revelio revelio.cpp
```

### Windows (MSVC x64 Native Tools)
```bat
cl /std:c++17 /O2 /EHsc revelio.cpp /Fe:revelio.exe
```

### Optional: Makefile
Drop this into a file named `Makefile`:
```make
CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -march=native -Wall -Wextra -DNDEBUG

revelio: revelio.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f revelio
```

---

## Quick Start

Heuristic priors (no model):
```bash
./revelio -i genome.fa --heuristic > calls.tsv
```

FASTQ with a logistic weights table and 7-mer context:
```bash
./revelio -i reads.fq --model weights.tsv -k 7 --header --out calls.tsv
```

Filter by context (only CG and CHG):
```bash
./revelio -i genome.fa --heuristic --contexts CG,CHG > calls.tsv
```

Stream a gzipped FASTA:
```bash
zcat genome.fa.gz | ./revelio -i - --heuristic > calls.tsv
```

---

## Usage
```text
revelio 0.1.1 — lightweight methylation context predictor
Usage:
  revelio -i <input.fa|fq|-> [--fasta|--fastq] [-k 5]
          [--heuristic | --model weights.tsv]
          [--contexts CG,CHG,CHH] [--header] [--out output.tsv]

Options:
  -i, --input <path>      Input FASTA/FASTQ; '-' or /dev/stdin reads stdin
      --fasta             Force FASTA (auto-detected by default)
      --fastq             Force FASTQ
  -k, --kmer <odd>        Centered k-mer size (default: 5; odd, >= 3)
      --heuristic         Use built-in heuristic priors (no weights table)
      --model <tsv>       Use logistic weights table
      --contexts <list>   Comma-separated subset (CG,CHG,CHH); default: all
      --header            Print header line in TSV output
  -o, --out <path>        Write to file (default: stdout)
  -h, --help              Show help and exit
      --version           Show version and exit
```

**Output columns (TSV)**  
`seq_id  pos_1based  context  strand  kmer  prob`  
`seq_id` is the header text up to the first whitespace. Invalid options or inputs exit with a non-zero status.

Example:
```text
Chr1    123     CG      +       ACGTA   0.903114
Chr1    456     CHH     +       CATTG   0.377540
```

> Note: Positions report the **C on the + strand**. Current version outputs `strand=+` (reverse-strand aggregation is on the roadmap).

---

## Modes

### Heuristic mode (`--heuristic`)
A lightweight prior that scores sites based on **context class** and simple **k-mer checks** (e.g., centered C, valid CG/CHG/CHH, optional flanking patterns). Intended for **fast ranking**, not calibrated probabilities.

### Logistic weights mode (`--model weights.tsv`)
Loads a user-provided **(context, k-mer) → weight** mapping and computes:
```
p = sigmoid( BIAS + weight(context, kmer) )
sigmoid(x) = 1 / (1 + e^(-x))
```
- `-k` must be **odd** (e.g., 5, 7, 9) and the k-mer is centered on the cytosine.
- If a `(context,kmer)` pair is **missing**, its score is `sigmoid(BIAS)` (0.5 when no `BIAS` row is given).
- The table is validated on load: unknown contexts, non-ACGT k-mers, mixed k-mer lengths, or a k-mer length different from `-k` are errors.
- `--heuristic` and `--model` are mutually exclusive; an unreadable `--model` is an error (no silent fallback).

---

## Weights Table Format

A simple **TSV without header**:

```text
# optional global bias row
BIAS    0.0

# context  kmer       weight
CG       ACGTA       2.25
CHG      CAGTC       1.10
CHH      CATTG       0.30
```

Guidelines:
- `context ∈ {CG, CHG, CHH}`
- `kmer` length must equal `-k` and be uppercase DNA (A/C/G/T).
- Use **`BIAS`** to shift all logits up/down (acts like an intercept).

---

## Input Notes

- **FASTA**: multi-line sequences supported; Windows (CRLF) line endings are handled.
- **FASTQ**: 4-line records, qualities ignored (sequence only).
- **Compression**: not handled internally—use shell pipes (`zcat | revelio -i /dev/stdin`).

---

## Output Schema (TSV)

| Column       | Description                                        |
|--------------|----------------------------------------------------|
| `seq_id`     | FASTA/FASTQ record identifier                      |
| `pos_1based` | 1-based coordinate of the **C** on the + strand    |
| `context`    | One of `CG`, `CHG`, `CHH`                          |
| `strand`     | Currently `+` (forward)                            |
| `kmer`       | Centered k-mer (odd `-k`, C is the center)         |
| `prob`       | Heuristic/logistic probability in **[0,1]**        |

---

## Performance Tips

- Compile with `-O3 -march=native` on a recent compiler.
- Prefer smaller odd `-k` (5 or 7) for speed; larger `k` can improve specificity but increases table size and cache pressure.
- Stream from disk or stdin. Records are processed one at a time, so peak memory is about the size of the largest sequence (e.g. one chromosome).

---

## How Revelio fits in your pipeline

- Use Revelio to **pre-rank** CG/CHG/CHH sites or **flag hotspots** before signal-level callers.
- For **ground-truth methylation** on ONT data, switch to **Remora/Dorado**, **Megalodon**, or **DeepSignal-plant** and then reconcile with Revelio’s ranked candidates.

---

## Roadmap

- Reverse-strand aggregation and duplex-aware reporting
- Optional BED/BEDMethyl output
- Native gzip input
- Pluggable k-mer featurizers (e.g., degeneracy, repeats)

---

## License

MIT. See [https://github.com/yykaya/Revelio/blob/main/LICENSE](LICENSE).

---

## Acknowledgements

Thanks to the plant epigenetics community for foundational work on CG/CHG/CHH methylation and to the ONT ecosystem for practical modified-base pipelines that complement Revelio.

#!/usr/bin/env python3
"""Regenerate the synthetic `revelio call` fixture in tests/data/ (deterministic).

The expected table is derived from the simulation itself, not from Revelio, so the
test checks placement against ground truth:

  call_ref.fa       record "loc": U1 + R + U2 + revcomp(R') + U3, where R' is R with
                    3 substitutions (an inverted, 99.4%-identical repeat), and a
                    second unrelated record "loc2"
  call_reads.sam.gz reads from both records and strands, all spanning unique flanks;
                    forward (flag 0), reverse (flag 16, SEQ reverse-complemented) and
                    unaligned (flag 4, SEQ as sequenced) records; some with A<->T
                    sequencing errors; plus records that must be ignored: secondary,
                    supplementary, hard-clipped, and an MN/sequence-length mismatch
  call_expected.tsv `revelio call --header` output a perfect caller gives

Every C of a read carries a call (all-context model). CpGs are methylated in
alternating 50-bp blocks (ML 240 / 15); non-CpG Cs get ML 10. Every fourth read uses
multi-entry tags as jasmine/Dorado write them: A+a. (6mA), C+hm? (5hmC and 5mC,
ML interleaved per position) and T-a. (opposite strand, every third T); only the
5mC values may be used.
"""
import gzip
import os
import random

random.seed(20261104)
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def rand_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


R = rand_seq(500)
Rp = list(R)
for pos in (120, 260, 400):
    Rp[pos] = {"A": "C", "C": "G", "G": "T", "T": "A"}[Rp[pos]]
Rp = "".join(Rp)
refs = {"loc": rand_seq(800) + R + rand_seq(700) + rc(Rp) + rand_seq(800), "loc2": rand_seq(1200)}


def cpg_ml(rec, c):
    """ML for the CpG whose + strand C is at c (0-based)."""
    return 240 if (c // 50) % 2 == 0 else 15


def context(ref, x, strand):
    if strand == "+":
        n1 = ref[x + 1] if x + 1 < len(ref) else "N"
        n2 = ref[x + 2] if x + 2 < len(ref) else "N"
    else:
        n1 = ref[x - 1].translate(COMP) if x >= 1 else "N"
        n2 = ref[x - 2].translate(COMP) if x >= 2 else "N"
    if n1 == "G":
        return "CG"
    if n1 not in "ACT" or n1 == "N":
        return None
    if n2 == "G":
        return "CHG"
    if n2 in "ACT":
        return "CHH"
    return None


expected = {}   # (rec, pos0, strand) -> [ctx, mod, tot, psum]
sam = ["@HD\tVN:1.6\tSO:unsorted"]
n = 0


def emit(rec, start, end, strand, flag, errors, multi=False):
    """Simulate one read of ref[start:end]; record expected calls."""
    global n
    ref = refs[rec]
    read = ref[start:end] if strand == "+" else rc(ref[start:end])
    read = list(read)
    # A<->T errors where both neighbours are A/T, so no C context changes
    cand = [i for i in range(30, len(read) - 30) if read[i] in "AT" and read[i - 1] in "AT" and read[i + 1] in "AT"]
    for i in random.sample(cand, errors):
        read[i] = "T" if read[i] == "A" else "A"
    read = "".join(read)
    ml = []
    for i, b in enumerate(read):
        if b != "C":
            continue
        # plus coordinate of the read base
        x = start + i if strand == "+" else end - 1 - i
        ctx = context(ref, x, strand)
        if ctx == "CG":
            c = x if strand == "+" else x - 1
            v = cpg_ml(rec, c)
            key = (rec, c, ".")
        else:
            v = 10
            key = (rec, x, strand)
        ml.append(v)
        if ctx is not None:
            e = expected.setdefault(key, [ctx, 0, 0, 0.0])
            e[1] += (v + 0.5) / 256 >= 0.5
            e[2] += 1
            e[3] += (v + 0.5) / 256
    if multi:
        n_a = read.count("A")
        n_t = len([i for i, b in enumerate(read) if b == "T"][::3])
        h = [random.randint(0, 40) for _ in ml]
        mm = ("A+a." + ",0" * n_a + ";C+hm?" + ",0" * len(ml) + ";T-a." +
              (",0" + ",2" * (n_t - 1) if n_t else "") + ";")
        vals = ([random.randint(0, 255) for _ in range(n_a)] + [v for pair in zip(h, ml) for v in pair] +
                [random.randint(0, 255) for _ in range(n_t)])
    else:
        mm = "C+m?" + "".join(",0" for _ in ml) + ";"
        vals = ml
    mlt = "ML:B:C," + ",".join(map(str, vals))
    seq = rc(read) if flag & 16 else read          # SAM stores reverse reads as + strand
    n += 1
    sam.append(f"r{n}\t{flag}\t{'*' if flag & 4 else rec}\t{0 if flag & 4 else start + 1}\t{0 if flag & 4 else 60}\t"
               f"{'*' if flag & 4 else str(len(seq)) + 'M'}\t*\t0\t0\t{seq}\t*\tMM:Z:{mm}\t{mlt}")
    return seq, mm, mlt


L = len(refs["loc"])
spans = []
for _ in range(60):                                 # loc: reads spanning unique flanks
    a = random.randint(5, 700)
    b = min(L - 5, a + random.randint(1500, 2600))
    spans.append(("loc", a, b))
for _ in range(15):
    a = random.randint(5, 300)
    spans.append(("loc2", a, a + random.randint(600, 850)))
for k, (rec, a, b) in enumerate(spans):
    strand = "+" if k % 2 == 0 else "-"
    flag = 4 if k % 5 == 4 else (0 if strand == "+" else 16)
    seq, mm, mlt = emit(rec, a, b, strand, flag, errors=2 if k % 3 == 0 else 0, multi=(k % 4 == 1))

# records that must be ignored (not added to `expected`)
junk = rand_seq(900)
sam.append(f"sec1\t256\tloc\t100\t0\t900M\t*\t0\t0\t{junk}\t*\tMM:Z:C+m?,0;\tML:B:C,255")
sam.append(f"sup1\t2048\tloc\t100\t0\t900M\t*\t0\t0\t{junk}\t*\tMM:Z:C+m?,0;\tML:B:C,255")
ref = refs["loc"]
sam.append(f"hc1\t0\tloc\t2001\t60\t50H900M\t*\t0\t0\t{ref[2000:2900]}\t*\tMM:Z:C+m?,0;\tML:B:C,255")
sam.append(f"mn1\t0\tloc\t2001\t60\t900M\t*\t0\t0\t{ref[2000:2900]}\t*\tMM:Z:C+m?,0;\tML:B:C,255\tMN:i:950")

with open(os.path.join(HERE, "call_ref.fa"), "w") as f:
    for name, s in refs.items():
        f.write(f">{name} synthetic\n" + "\n".join(s[i:i + 70] for i in range(0, len(s), 70)) + "\n")
with gzip.open(os.path.join(HERE, "call_reads.sam.gz"), "wt") as f:
    f.write("\n".join(sam) + "\n")
order = {"loc": 0, "loc2": 1}
with open(os.path.join(HERE, "call_expected.tsv"), "w") as f:
    f.write("seq_id\tpos_1based\tcontext\tstrand\tmod\ttot\tpct\tmean_prob\n")
    for (rec, pos, strand), (ctx, mod, tot, ps) in sorted(expected.items(), key=lambda kv: (order[kv[0][0]], kv[0][1], kv[0][2])):
        f.write(f"{rec}\t{pos + 1}\t{ctx}\t{strand}\t{mod}\t{tot}\t{100 * mod / tot:.1f}\t{ps / tot:.4f}\n")
print(f"{len(spans)} reads, {len(expected)} expected sites")

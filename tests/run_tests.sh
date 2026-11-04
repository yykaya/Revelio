#!/usr/bin/env bash
# Revelio regression tests — run via `make test` (or: tests/run_tests.sh ./revelio)
set -u

BIN=${1:-./revelio}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
pass=0; fail=0

ok()   { echo "ok   - $1"; pass=$((pass+1)); }
bad()  { echo "FAIL - $1"; fail=$((fail+1)); }

# expect_rc <name> <expected-rc> <cmd...>
expect_rc() {
    local name=$1 want=$2; shift 2
    "$@" >"$TMP/out" 2>"$TMP/err"; local got=$?
    if [ "$got" -eq "$want" ]; then ok "$name"; else bad "$name (rc=$got, want $want)"; cat "$TMP/err"; fi
}

# expect_out <name> <expected-file> <cmd...>
expect_out() {
    local name=$1 want=$2; shift 2
    if "$@" >"$TMP/out" 2>"$TMP/err" && diff -u "$want" "$TMP/out" >"$TMP/diff"; then
        ok "$name"
    else
        bad "$name"; cat "$TMP/diff" "$TMP/err"
    fi
}

# ---- fixtures ----
printf '>chr1 some description\nAACGTACCAG\nTTTCTGAACCAA\n>chr2\nGGCGAA\n' > "$TMP/t.fa"
printf '>chr1 some description\r\nAACGTACCAG\r\nTTTCTGAACCAA\r\n>chr2\r\nGGCGAA\r\n' > "$TMP/crlf.fa"
printf '@r1 extra\nAACGTACCAG\n+\nIIIIIIIIII\nnoise\n@r2\nTTCTGAA\n+\nIIIIIII\n' > "$TMP/t.fq"
printf 'BIAS\t0.0\nCG\tAACGT\t2.0\n' > "$TMP/w5.tsv"
printf 'CG\tAAACGTA\t2.0\n' > "$TMP/w7.tsv"
printf 'CG\tAACGT\t2.0\nCHX\tAACGT\t1.0\n' > "$TMP/badctx.tsv"
printf 'CG\tAACGT\t2.0\nCG\tAACGTAA\t1.0\n' > "$TMP/mixk.tsv"
printf 'BIAS\t0.5\nBIAS_CG\t1.0\nBIAS_CHH\t-1.0\nCG\tAACGT\t0.5\n' > "$TMP/ctxbias.tsv"
printf 'BIAS_CA\t1.0\n' > "$TMP/badbias.tsv"

cat > "$TMP/exp_heur.tsv" <<'EOF'
seq_id	pos_1based	context	strand	kmer	prob
chr1	3	CG	+	AACGT	0.850000
chr1	4	CG	-	TACGT	0.850000
chr1	7	CHH	+	TACCA	0.250000
chr1	8	CHG	+	ACCAG	0.550000
chr1	10	CHG	-	AACTG	0.550000
chr1	14	CHG	+	TTCTG	0.550000
chr1	16	CHG	-	TTCAG	0.550000
chr1	19	CHH	+	AACCA	0.250000
chr1	20	CHH	+	ACCAA	0.250000
chr2	3	CG	+	GGCGA	0.850000
chr2	4	CG	-	TTCGC	0.850000
EOF
tail -n +2 "$TMP/exp_heur.tsv" > "$TMP/exp_heur_nohdr.tsv"
awk -F'\t' '$4=="+"' "$TMP/exp_heur_nohdr.tsv" > "$TMP/exp_plus.tsv"
awk -F'\t' '$4=="-"' "$TMP/exp_heur_nohdr.tsv" > "$TMP/exp_minus.tsv"

# BIAS 0.5 + BIAS_CG 1.0 (+0.5 for AACGT) / BIAS 0.5 + BIAS_CHH -1.0 / BIAS 0.5 only (CHG)
cat > "$TMP/exp_ctxbias.tsv" <<'EOF'
chr1	3	CG	+	AACGT	0.880797
chr1	4	CG	-	TACGT	0.817574
chr1	7	CHH	+	TACCA	0.377541
chr1	8	CHG	+	ACCAG	0.622459
EOF

cat > "$TMP/exp_model.tsv" <<'EOF'
chr1	3	CG	+	AACGT	0.880797
chr1	4	CG	-	TACGT	0.500000
EOF

cat > "$TMP/exp_fq.tsv" <<'EOF'
r1	3	CG	+	AACGT	0.850000
r1	4	CG	-	TACGT	0.850000
r1	7	CHH	+	TACCA	0.250000
r1	8	CHG	+	ACCAG	0.550000
r2	3	CHG	+	TTCTG	0.550000
r2	5	CHG	-	TTCAG	0.550000
EOF

# ---- output correctness ----
expect_out "FASTA file, heuristic, header, id trimmed" "$TMP/exp_heur.tsv" "$BIN" -i "$TMP/t.fa" --heuristic --header
expect_out "stdin via pipe (-i -)"            "$TMP/exp_heur_nohdr.tsv" sh -c "cat '$TMP/t.fa' | '$BIN' -i - --heuristic"
expect_out "stdin via pipe (-i /dev/stdin)"   "$TMP/exp_heur_nohdr.tsv" sh -c "cat '$TMP/t.fa' | '$BIN' -i /dev/stdin --heuristic"
expect_out "gzip pipe"                        "$TMP/exp_heur_nohdr.tsv" sh -c "gzip -c '$TMP/t.fa' | gzip -dc | '$BIN' -i - --heuristic"
expect_out "CRLF line endings"                "$TMP/exp_heur_nohdr.tsv" "$BIN" -i "$TMP/crlf.fa" --heuristic
expect_out "FASTQ auto-detect, skips noise"   "$TMP/exp_fq.tsv" "$BIN" -i "$TMP/t.fq" --heuristic
expect_out "FASTQ via pipe"                   "$TMP/exp_fq.tsv" sh -c "cat '$TMP/t.fq' | '$BIN' -i - --heuristic"
expect_out "model scoring"                    "$TMP/exp_model.tsv" sh -c "'$BIN' -i '$TMP/t.fa' --model '$TMP/w5.tsv' --contexts CG,CHH | head -2"
expect_out "--strand + (forward only)"        "$TMP/exp_plus.tsv"  "$BIN" -i "$TMP/t.fa" --heuristic --strand +
expect_out "--strand - (reverse only)"        "$TMP/exp_minus.tsv" "$BIN" -i "$TMP/t.fa" --heuristic --strand -
expect_out "per-context BIAS_<ctx> rows"      "$TMP/exp_ctxbias.tsv" sh -c "'$BIN' -i '$TMP/t.fa' --model '$TMP/ctxbias.tsv' | head -4"
expect_out "--kmer/--input/-o aliases"        "$TMP/exp_heur_nohdr.tsv" sh -c "'$BIN' --input '$TMP/t.fa' --kmer 5 --heuristic -o '$TMP/o.tsv' && cat '$TMP/o.tsv'"

# ---- exit codes ----
expect_rc "--help exits 0"                 0 "$BIN" --help
expect_rc "--version exits 0"              0 "$BIN" --version
expect_rc "unknown option exits 1"         1 "$BIN" -i "$TMP/t.fa" --bogus
expect_rc "non-numeric -k exits 1"         1 "$BIN" -i "$TMP/t.fa" -k abc
expect_rc "even -k exits 1"                1 "$BIN" -i "$TMP/t.fa" -k 4
expect_rc "bad --contexts exits 1"         1 "$BIN" -i "$TMP/t.fa" --contexts CG,CHX
expect_rc "missing input exits 2"          2 "$BIN" -i "$TMP/nope.fa"
expect_rc "empty input exits 2"            2 sh -c "printf '' | '$BIN' -i -"
expect_rc "missing model exits 2"          2 "$BIN" -i "$TMP/t.fa" --model "$TMP/nope.tsv"
expect_rc "model k != -k exits 1"          1 "$BIN" -i "$TMP/t.fa" --model "$TMP/w7.tsv"
expect_rc "model bad context exits 2"      2 "$BIN" -i "$TMP/t.fa" --model "$TMP/badctx.tsv"
expect_rc "model mixed k exits 2"          2 "$BIN" -i "$TMP/t.fa" --model "$TMP/mixk.tsv"
expect_rc "bad --strand exits 1"           1 "$BIN" -i "$TMP/t.fa" --strand x
expect_rc "unknown BIAS_ row exits 2"      2 "$BIN" -i "$TMP/t.fa" --model "$TMP/badbias.tsv"
expect_rc "--heuristic + --model exits 1"  1 "$BIN" -i "$TMP/t.fa" --heuristic --model "$TMP/w5.tsv"

echo "# $pass passed, $fail failed"
[ "$fail" -eq 0 ]

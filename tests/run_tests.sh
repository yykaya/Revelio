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

# --window 10 on t.fa: hand-checked counts, GC and CpG o/e
printf 'seq_id\tstart\tend\tn_acgt\tgc\tcpg_oe\tCG_plus\tCG_minus\tCHG_plus\tCHG_minus\tCHH_plus\tCHH_minus\tp_CG\tp_CHG\tp_CHH\n' > "$TMP/exp_win.tsv"
cat >> "$TMP/exp_win.tsv" <<'EOF'
chr1	0	10	10	0.5000	1.6667	1	1	1	1	1	0	0.850000	0.550000	0.250000
chr1	10	20	10	0.4000	0.0000	0	0	1	1	2	0	NA	0.550000	0.250000
chr1	20	22	2	0.0000	NA	0	0	0	0	0	0	NA	NA	NA
chr2	0	6	6	0.6667	2.0000	1	1	0	0	0	0	0.850000	NA	NA
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
# stdin buffer refills: one 3 MB line (crosses the 1 MiB buffer twice) piped vs read from file
awk 'BEGIN{srand(11); printf ">long\n"; for(i=0;i<3000000;i++) printf "%s", substr("ACGT", int(rand()*4)+1, 1); printf "\n"}' > "$TMP/long.fa"
"$BIN" -i "$TMP/long.fa" --heuristic --window 100000 > "$TMP/long_file.tsv"
expect_out "stdin across buffer refills"      "$TMP/long_file.tsv" sh -c "cat '$TMP/long.fa' | '$BIN' -i - --heuristic --window 100000"

expect_out "--window profile"                 "$TMP/exp_win.tsv" "$BIN" -i "$TMP/t.fa" --heuristic --window 10 --header

# window counts must equal the per-site rows of the same command
EX=$(dirname "$0")/../examples/col0_chr1_3631-4630.fa
WT=$(dirname "$0")/../weights.tsv
for args in "--heuristic" "--model $WT --strand -" "--heuristic --contexts CHG,CHH -k 7"; do
    sites=$("$BIN" -i "$EX" $args | awk -F'\t' '{n[$3 $4]++} END{for(k in n) printf "%s=%d ", k, n[k]}' | tr ' ' '\n' | sort | tr '\n' ' ')
    wins=$("$BIN" -i "$EX" $args --window 137 | awk -F'\t' '{a+=$7;b+=$8;c+=$9;d+=$10;e+=$11;f+=$12}
        END{if(a)printf "CG+=%d\n",a; if(b)printf "CG-=%d\n",b; if(c)printf "CHG+=%d\n",c; if(d)printf "CHG-=%d\n",d; if(e)printf "CHH+=%d\n",e; if(f)printf "CHH-=%d\n",f}' | sort | tr '\n' ' ')
    if [ "$sites" = "$wins" ] && [ -n "$sites" ]; then ok "window counts = site rows ($args)"; else bad "window counts = site rows ($args): $sites vs $wins"; fi
done
# overlapping windows: every window = sum of its non-overlapping halves
"$BIN" -i "$EX" --heuristic --window 100 --step 50 | awk -F'\t' '$3-$2==100{print $2, $7+$8+$9+$10+$11+$12}' > "$TMP/w100"
"$BIN" -i "$EX" --heuristic --window 50 | awk -F'\t' '{n[$2]=$7+$8+$9+$10+$11+$12} END{for(s in n) print s, n[s]}' > "$TMP/w50"
if awk 'NR==FNR{h[$1]=$2; next} {if(h[$1]+h[$1+50]!=$2) bad++} END{exit bad>0}' "$TMP/w50" "$TMP/w100"; then ok "--step windows = sum of blocks"; else bad "--step windows = sum of blocks"; fi
# --per-record: one row per record (hand-checked: chr1 6 C, 3 G, 1 CpG -> o/e 22/18),
# identical to a window longer than every record
printf 'chr1\t0\t22\t22\t0.4091\t1.2222\t1\t1\t2\t2\t3\t0\t0.850000\t0.550000\t0.250000\nchr2\t0\t6\t6\t0.6667\t2.0000\t1\t1\t0\t0\t0\t0\t0.850000\tNA\tNA\n' > "$TMP/exp_rec.tsv"
expect_out "--per-record (multi-record)"      "$TMP/exp_rec.tsv" "$BIN" -i "$TMP/t.fa" --heuristic --per-record
"$BIN" -i "$EX" --model "$WT" --window 1000000000 --header > "$TMP/rec_win.tsv"
expect_out "--per-record = window > record"   "$TMP/rec_win.tsv" "$BIN" -i "$EX" --model "$WT" --per-record --header

# ---- revelio call: synthetic reads with known truth (tests/make_call_fixture.py) ----
DATA=$(dirname "$0")/data
gzip -dc "$DATA/call_reads.sam.gz" > "$TMP/reads.sam"
expect_out "call: placement = ground truth"   "$DATA/call_expected.tsv" sh -c "gzip -dc '$DATA/call_reads.sam.gz' | '$BIN' call -r '$DATA/call_ref.fa' --header | cut -f1-8"
"$BIN" call -r "$DATA/call_ref.fa" -i "$TMP/reads.sam" --stats "$TMP/cstats.tsv" > /dev/null
if grep -qx "reads_clipped	2" "$TMP/cstats.tsv" && grep -qx "calls_unplaced	0" "$TMP/cstats.tsv"; then ok "call: clipped/MN records skipped, all calls placed"; else bad "call: stats"; cat "$TMP/cstats.tsv"; fi
"$BIN" call -r "$DATA/call_ref.fa" -i "$TMP/reads.sam" --no-combine > "$TMP/split.tsv"
"$BIN" call -r "$DATA/call_ref.fa" -i "$TMP/reads.sam" > "$TMP/comb.tsv"
if [ "$(awk -F'\t' '$3=="CG"{t+=$6} END{print t}' "$TMP/split.tsv")" = "$(awk -F'\t' '$3=="CG"{t+=$6} END{print t}' "$TMP/comb.tsv")" ] \
   && ! awk -F'\t' '$4=="."' "$TMP/split.tsv" | grep -q .; then ok "call: --no-combine keeps CpG strands apart"; else bad "call: --no-combine"; fi
if awk -F'\t' '$6<20{bad=1} END{exit bad}' <("$BIN" call -r "$DATA/call_ref.fa" -i "$TMP/reads.sam" --min-cov 20); then ok "call: --min-cov filters"; else bad "call: --min-cov"; fi
expect_rc "call: missing -r exits 1"        1 "$BIN" call -i "$TMP/reads.sam"
expect_rc "call: even -k exits 1"           1 "$BIN" call -r "$DATA/call_ref.fa" -k 20 -i "$TMP/reads.sam"
expect_rc "call: non-SAM input exits 2"     2 "$BIN" call -r "$DATA/call_ref.fa" -i "$TMP/t.fa"
expect_rc "call: --help exits 0"            0 "$BIN" call --help

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
expect_rc "--step without --window exits 1"  1 "$BIN" -i "$TMP/t.fa" --step 5
expect_rc "--window not multiple of step"   1 "$BIN" -i "$TMP/t.fa" --window 10 --step 3
expect_rc "--window 0 exits 1"              1 "$BIN" -i "$TMP/t.fa" --window 0
expect_rc "--per-record + --window exits 1" 1 "$BIN" -i "$TMP/t.fa" --per-record --window 10
expect_rc "--heuristic + --model exits 1"  1 "$BIN" -i "$TMP/t.fa" --heuristic --model "$TMP/w5.tsv"

echo "# $pass passed, $fail failed"
[ "$fail" -eq 0 ]

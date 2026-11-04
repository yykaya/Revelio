/*
  Revelio — a lightweight methylation-context predictor for FASTA/FASTQ
  Author: Yasin Kaya | MPIPZ - MPG, Cologne
  Date  : 18.08.2025
  License: MIT

  Features
  - Detect CG / CHG / CHH contexts (H ∈ {A, C, T}) around cytosines on both strands.
  - Extract a centered k-mer (odd k, default k=5).
  - Predict methylation probability via:
      (A) heuristic priors (CG > CHG > CHH), or
      (B) logistic regression weights: prob = sigmoid(bias + w[context,kmer]).
  - Outputs TSV: seq_id  pos_1based  context  strand  kmer  prob
  - FASTA (multi-line) and FASTQ (4-line) support with auto-detection; reads pipes/stdin.
  - Practical CLI: --header, --contexts, --out, --version.
  - --window/--step: per-window profiles (context counts per strand, GC, CpG o/e, mean prior).
  - seq_id is the header up to the first whitespace.

  Notes / Limitations
  - Minus-strand sites are reported at the + strand coordinate of their G, with the k-mer
    and context read 5'->3' on the minus strand.
  - Ambiguity (non-ACGT) breaks context; such windows are skipped.
  - This is a fast *sequence-context prior*. It does NOT analyze raw signal.
    For signal-level methylation (e.g., ONT), see deep learning tools such as DeepSignal-plant.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <sstream>
#include <iomanip>
#include <functional>

namespace rv {

static constexpr const char* kVersion = "0.3.0";

// ---------- small helpers ----------
inline char up(char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }
inline bool is_acgt(char b) { b = up(b); return b=='A'||b=='C'||b=='G'||b=='T'; }
inline bool is_H(char b)    { b = up(b); return b=='A'||b=='C'||b=='T'; }

// Logistic sigmoid for probability mapping
inline double sigmoid(double x) {
    if (x >= 0) {
        const double e = std::exp(-x);
        return 1.0 / (1.0 + e);
    } else {
        const double e = std::exp(x);
        return e / (1.0 + e);
    }
}

// Trim ASCII whitespace
inline std::string trim(const std::string& s) {
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Reverse complement (minus-strand k-mers)
inline char comp(char b) {
    b = up(b);
    switch (b) { case 'A': return 'T'; case 'T': return 'A'; case 'C': return 'G'; case 'G': return 'C'; default: return 'N'; }
}
inline std::string revcomp(const std::string& s) {
    std::string t; t.reserve(s.size());
    for (int i = static_cast<int>(s.size()) - 1; i >= 0; --i) t.push_back(comp(s[i]));
    return t;
}

// ---------- weights model ----------
// TSV (no header):
//   context  kmer  weight
// e.g.
//   BIAS   -0.20      # optional global bias
//   BIAS_CG  1.00     # optional per-context intercepts (BIAS_CG / BIAS_CHG / BIAS_CHH)
//   CG     ACGTA  2.10
//   CHG    CCGTC  0.80
//   CHH    CATTG -0.40
struct Model {
    double bias = 0.0;
    std::unordered_map<std::string, double> ctx_bias; // CG/CHG/CHH -> intercept
    int kmer_len = 0;                          // set from the first k-mer row
    std::unordered_map<std::string, double> w; // key = context + '\t' + kmer

    // Returns false (with a message in err) on unreadable file or malformed rows
    bool load(const std::string& path, std::string& err) {
        std::ifstream in(path);
        if (!in) { err = "cannot open " + path; return false; }
        std::string line;
        int lineno = 0;
        while (std::getline(in, line)) {
            ++lineno;
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;
            const std::string where = path + ":" + std::to_string(lineno) + ": ";

            // BIAS <number>  or  BIAS_<context> <number>
            if (line.rfind("BIAS", 0) == 0) {
                std::istringstream iss(line);
                std::string tag; double wt;
                if (!(iss >> tag >> wt)) { err = where + tag + " needs a number"; return false; }
                for (char& c : tag) c = up(c);
                if (tag == "BIAS") {
                    bias = wt;
                } else if (tag == "BIAS_CG" || tag == "BIAS_CHG" || tag == "BIAS_CHH") {
                    ctx_bias[tag.substr(5)] = wt;
                } else {
                    err = where + "unknown bias row '" + tag + "' (use BIAS, BIAS_CG, BIAS_CHG or BIAS_CHH)";
                    return false;
                }
                continue;
            }

            // context kmer weight
            std::istringstream iss(line);
            std::string ctx, kmer; double wt;
            if (!(iss >> ctx >> kmer >> wt)) { err = where + "expected <context> <kmer> <weight>"; return false; }
            for (char& c : ctx)  c = up(c);
            for (char& c : kmer) c = up(c);
            if (ctx != "CG" && ctx != "CHG" && ctx != "CHH") { err = where + "unknown context '" + ctx + "'"; return false; }
            for (char c : kmer) {
                if (!is_acgt(c)) { err = where + "k-mer '" + kmer + "' has non-ACGT bases"; return false; }
            }
            const int len = static_cast<int>(kmer.size());
            if (kmer_len == 0) kmer_len = len;
            if (len != kmer_len) {
                err = where + "k-mer length " + std::to_string(len) + " differs from earlier rows (" + std::to_string(kmer_len) + ")";
                return false;
            }
            w[ctx + "\t" + kmer] = wt;
        }
        return true;
    }

    // p = sigmoid(BIAS + BIAS_<ctx> + w[ctx,kmer]); missing terms count as 0
    double score(const std::string& ctx, const std::string& kmer) const {
        double s = bias;
        const auto cb = ctx_bias.find(ctx);
        if (cb != ctx_bias.end()) s += cb->second;
        const auto it = w.find(ctx + "\t" + kmer);
        if (it != w.end()) s += it->second;
        return sigmoid(s);
    }
};

// ---------- FASTA/FASTQ reader ----------
enum class FileType { FASTA, FASTQ };

// std::cin on libc++ (macOS) reads stdin almost byte by byte even without stdio sync;
// this buffer fills 1 MiB at a time with fread.
class StdinBuf : public std::streambuf {
public:
    StdinBuf() : buf_(1 << 20) {}
protected:
    int_type underflow() override {
        if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
        const std::size_t n = std::fread(buf_.data(), 1, buf_.size(), stdin);
        if (n == 0) return traits_type::eof();
        setg(buf_.data(), buf_.data(), buf_.data() + n);
        return traits_type::to_int_type(*gptr());
    }
private:
    std::vector<char> buf_;
};

// Buffered stream over stdin, for "-i -"
inline std::istream& stdin_stream() {
    static StdinBuf buf;
    static std::istream in(&buf);
    return in;
}

// Record ID = header text up to the first whitespace (as samtools/bedtools expect)
inline std::string header_id(const std::string& line) {
    const std::string h = trim(line.substr(1));
    return h.substr(0, h.find_first_of(" \t"));
}

// Append sequence characters, skipping whitespace (e.g. '\r' from CRLF files)
inline void append_seq(std::string& seq, const std::string& line) {
    for (char c : line) {
        if (std::isspace(static_cast<unsigned char>(c))) continue;
        seq.push_back(is_acgt(c) ? up(c) : 'N');
    }
}

// Reads forward only (no tellg/seekg), so pipes and /dev/stdin work.
struct Reader {
    std::istream& in;
    FileType type;
    std::string cur_id, cur_seq;
    std::string pending;          // FASTA header already read for the next record
    bool has_pending = false;

    explicit Reader(std::istream& input, FileType t) : in(input), type(t) {}

    // Read the next record (id, seq). Returns false on EOF.
    bool next() {
        cur_id.clear(); cur_seq.clear();
        std::string line;

        if (type == FileType::FASTA) {
            // Find next '>' header
            if (has_pending) {
                line = pending; has_pending = false;
            } else {
                while (std::getline(in, line)) {
                    if (!line.empty() && line[0] == '>') break;
                }
                if (line.empty() || line[0] != '>') return false;
            }
            cur_id = header_id(line);

            // Accumulate sequence lines until next header or EOF
            while (std::getline(in, line)) {
                if (!line.empty() && line[0] == '>') { pending = line; has_pending = true; break; }
                append_seq(cur_seq, line);
            }
            return true;
        } else {
            // FASTQ: @id, seq, +, qual
            do {
                if (!std::getline(in, line)) return false;
            } while (line.empty() || line[0] != '@'); // skip until a header
            cur_id = header_id(line);
            std::string seq, plus, qual;
            if (!std::getline(in, seq))  return false;
            if (!std::getline(in, plus)) return false;
            if (!std::getline(in, qual)) return false;
            cur_seq.reserve(seq.size());
            append_seq(cur_seq, seq);
            return true;
        }
    }
};

// Build a centered k-mer around position p (the cytosine).
// If window exceeds borders or contains ambiguous bases, return empty string.
inline std::string centered_kmer(const std::string& s, std::size_t p, int k) {
    if (k % 2 == 0) return "";
    const int half = k / 2;
    if (static_cast<int>(p) - half < 0) return "";
    if (p + half >= s.size()) return "";
    std::string km; km.reserve(k);
    for (int i = static_cast<int>(p) - half; i <= static_cast<int>(p) + half; ++i) {
        const char b = s[static_cast<std::size_t>(i)];
        if (!is_acgt(b)) return "";
        km.push_back(b);
    }
    return km;
}

struct Options {
    std::string infile;
    std::string mode = "auto";      // "auto"|"fasta"|"fastq"
    std::string model_path = "";
    std::string outfile = "";
    bool heuristic = false;
    bool print_header = false;
    int  k = 5;
    std::string strand = "both";    // "both"|"+"|"-"
    std::int64_t window = 0;        // > 0: one row per window instead of per site
    std::int64_t step = 0;          // window step (default: window)
    std::unordered_set<std::string> ctx_filter; // empty = accept all
};

// Parse comma-separated context list like "CG,CHG,CHH"
inline std::unordered_set<std::string> parse_ctx_list(const std::string& s) {
    std::unordered_set<std::string> S;
    std::istringstream iss(s);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        tok = trim(tok);
        for (char& c : tok) c = up(c);
        if (tok=="CG" || tok=="CHG" || tok=="CHH") {
            S.insert(tok);
        } else {
            std::cerr << "Error: unknown context '" << tok << "'\n";
            return {}; // caller treats empty as an error
        }
    }
    return S;
}

inline void print_help(std::ostream& os) {
    os <<
    "revelio " << kVersion << " — lightweight methylation context predictor\n"
    "Usage:\n"
    "  revelio -i <input.fa|fq|-> [--fasta|--fastq] [-k 5]\n"
    "          [--heuristic | --model weights.tsv]\n"
    "          [--contexts CG,CHG,CHH] [--strand both] [--header] [--out output.tsv]\n"
    "\n"
    "Options:\n"
    "  -i, --input <path>      Input FASTA/FASTQ; '-' or /dev/stdin reads stdin\n"
    "      --fasta             Force FASTA (auto-detected by default)\n"
    "      --fastq             Force FASTQ\n"
    "  -k, --kmer <odd>        Centered k-mer size (default: 5; odd, >= 3)\n"
    "      --heuristic         Use built-in heuristic priors (no weights table)\n"
    "      --model <tsv>       Use logistic weights table\n"
    "      --contexts <list>   Comma-separated subset (CG,CHG,CHH); default: all\n"
    "      --strand <s>        both | + | - (default: both)\n"
    "      --window <bp>       Write one profile row per window instead of one row per site\n"
    "      --step <bp>         Window step (default: --window; --window must be a multiple)\n"
    "      --header            Print header line in TSV output\n"
    "  -o, --out <path>        Write to file (default: stdout)\n"
    "  -h, --help              Show help and exit\n"
    "      --version           Show version and exit\n"
    "\n"
    "Output columns (TSV): seq_id  pos_1based  context  strand  kmer  prob\n"
    "With --window:        seq_id  start  end  n_acgt  gc  cpg_oe  CG_plus  CG_minus\n"
    "                      CHG_plus  CHG_minus  CHH_plus  CHH_minus  p_CG  p_CHG  p_CHH\n"
    "\n"
    "Examples:\n"
    "  revelio -i genome.fa --heuristic > calls.tsv\n"
    "  zcat genome.fa.gz | revelio -i - --heuristic > calls.tsv\n"
    "  revelio -i reads.fq --fastq --model weights.tsv -k 7 --header --out calls.tsv\n"
    "  revelio -i genome.fa --heuristic --window 10000 --header > profile.tsv\n";
}

// Parse a whole-string positive 64-bit integer; false on junk like "abc", "5x" or "0"
inline bool parse_pos_i64(const std::string& s, std::int64_t& v) {
    try {
        std::size_t used = 0;
        v = static_cast<std::int64_t>(std::stoll(s, &used));
        return used == s.size() && v > 0;
    } catch (const std::exception&) {
        return false;
    }
}

// Parse a whole-string integer; false on junk like "abc" or "5x"
inline bool parse_int(const std::string& s, int& v) {
    try {
        std::size_t used = 0;
        v = std::stoi(s, &used);
        return used == s.size();
    } catch (const std::exception&) {
        return false;
    }
}

// Returns -1 to continue, otherwise the process exit code (0 for --help/--version)
inline int parse_cli(int argc, char* argv[], Options& opt) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        const bool takes_value = (a == "-i" || a == "--in" || a == "--input" ||
                                  a == "-k" || a == "--kmer" || a == "--model" ||
                                  a == "--contexts" || a == "--strand" ||
                                  a == "--window" || a == "--step" ||
                                  a == "-o" || a == "--out");
        if (takes_value && i+1 >= argc) { std::cerr << "Error: missing value for " << a << "\n"; return 1; }

        if (a == "-i" || a == "--in" || a == "--input") {
            opt.infile = argv[++i];
        } else if (a == "--fastq") {
            opt.mode = "fastq";
        } else if (a == "--fasta") {
            opt.mode = "fasta";
        } else if (a == "-k" || a == "--kmer") {
            const std::string v = argv[++i];
            if (!parse_int(v, opt.k)) { std::cerr << "Error: -k expects an odd integer, got '" << v << "'\n"; return 1; }
        } else if (a == "--model") {
            opt.model_path = argv[++i];
        } else if (a == "--heuristic") {
            opt.heuristic = true;
        } else if (a == "--header") {
            opt.print_header = true;
        } else if (a == "--contexts") {
            opt.ctx_filter = parse_ctx_list(argv[++i]);
            if (opt.ctx_filter.empty()) {
                std::cerr << "Error: --contexts must list CG, CHG and/or CHH\n";
                return 1;
            }
        } else if (a == "--strand") {
            opt.strand = argv[++i];
            if (opt.strand != "both" && opt.strand != "+" && opt.strand != "-") {
                std::cerr << "Error: --strand must be both, + or -\n";
                return 1;
            }
        } else if (a == "--window" || a == "--step") {
            const std::string v = argv[++i];
            if (!parse_pos_i64(v, a == "--window" ? opt.window : opt.step)) {
                std::cerr << "Error: " << a << " expects a positive integer, got '" << v << "'\n";
                return 1;
            }
        } else if (a == "-o" || a == "--out") {
            opt.outfile = argv[++i];
        } else if (a == "-h" || a == "--help") {
            print_help(std::cout);
            return 0;
        } else if (a == "--version") {
            std::cout << "revelio " << kVersion << "\n";
            return 0;
        } else {
            std::cerr << "Error: unknown option: " << a << " (see --help)\n";
            return 1;
        }
    }
    return -1;
}

// Context (CG/CHG/CHH) of the cytosine at position i on the given strand, or nullptr.
// '+': s[i]=='C', downstream bases are s[i+1], s[i+2].
// '-': s[i]=='G' (a C on the reverse strand); its downstream bases are comp(s[i-1]), comp(s[i-2]).
inline const char* context_at(const std::string& s, std::size_t i, char strand) {
    char n1 = 'N', n2 = 'N';
    if (strand == '+') {
        if (i + 1 < s.size()) n1 = s[i + 1];
        if (i + 2 < s.size()) n2 = s[i + 2];
    } else {
        if (i >= 1) n1 = comp(s[i - 1]);
        if (i >= 2) n2 = comp(s[i - 2]);
    }
    if (n1 == 'G') return "CG";
    if (!is_H(n1)) return nullptr;
    if (n2 == 'G') return "CHG";
    if (is_H(n2))  return "CHH";
    return nullptr;
}

using Predictor = std::function<double(const std::string&, const std::string&)>;

// Call f(i, ctx, strand, kmer, prob) for every reported site of seq.
// Per-site output and window profiles share this, so their counts always agree.
template <class F>
inline void for_each_site(const std::string& seq, int k, const Predictor& predictor,
                          const std::unordered_set<std::string>& ctx_filter,
                          bool scan_plus, bool scan_minus, F&& f)
{
    for (std::size_t i = 0; i < seq.size(); ++i) {
        const char b = seq[i];
        char strand;
        if (b == 'C' && scan_plus)       strand = '+';
        else if (b == 'G' && scan_minus) strand = '-';   // C on the reverse strand
        else continue;

        const char* ctx = context_at(seq, i, strand);
        if (!ctx) continue;

        if (!ctx_filter.empty() && !ctx_filter.count(ctx)) continue;

        std::string km = centered_kmer(seq, i, k);
        if (km.empty()) continue;
        if (strand == '-') km = revcomp(km);   // read 5'->3' on the C's own strand

        f(i, ctx, strand, km, predictor(ctx, km));
    }
}

// Scan one sequence and write hits
inline void scan_sequence(const std::string& id,
                          const std::string& seq,
                          int k,
                          const Predictor& predictor,
                          const std::unordered_set<std::string>& ctx_filter,
                          bool scan_plus,
                          bool scan_minus,
                          std::ostream& out)
{
    for_each_site(seq, k, predictor, ctx_filter, scan_plus, scan_minus,
        [&](std::size_t i, const char* ctx, char strand, const std::string& km, double prob) {
            // seq_id  pos_1based  context  strand  kmer  prob
            out << id << '\t'
                << (i + 1) << '\t'
                << ctx << '\t'
                << strand << '\t'
                << km << '\t'
                << std::fixed << std::setprecision(6) << prob << '\n';
        });
}

// Per-block tallies for --window; a window is the sum of window/step consecutive blocks
struct Block {
    std::int64_t n_acgt = 0, n_c = 0, n_g = 0, n_cpg = 0;
    std::int64_t n[3][2] = {};  // [CG, CHG, CHH][+, -] site counts
    double psum[3] = {};        // summed prior per context

    // sign = +1 adds o, -1 removes it (running window sum)
    void add(const Block& o, int sign = 1) {
        n_acgt += sign * o.n_acgt; n_c += sign * o.n_c; n_g += sign * o.n_g; n_cpg += sign * o.n_cpg;
        for (int c = 0; c < 3; ++c) {
            n[c][0] += sign * o.n[c][0]; n[c][1] += sign * o.n[c][1]; psum[c] += sign * o.psum[c];
        }
    }
};

inline int ctx_index(const char* ctx) {
    return ctx[1] == 'G' ? 0 : (ctx[2] == 'G' ? 1 : 2);   // CG, CHG, CHH
}

inline void print_profile_header(std::ostream& out) {
    out << "seq_id\tstart\tend\tn_acgt\tgc\tcpg_oe\tCG_plus\tCG_minus\tCHG_plus\tCHG_minus"
           "\tCHH_plus\tCHH_minus\tp_CG\tp_CHG\tp_CHH\n";
}

// Profile one sequence in windows of `window` bp moving by `step` bp (BED coordinates)
inline void profile_sequence(const std::string& id,
                             const std::string& seq,
                             int k,
                             const Predictor& predictor,
                             const std::unordered_set<std::string>& ctx_filter,
                             bool scan_plus,
                             bool scan_minus,
                             std::int64_t window,
                             std::int64_t step,
                             std::ostream& out)
{
    const std::int64_t len = static_cast<std::int64_t>(seq.size());
    if (len == 0) return;
    std::vector<Block> blocks(static_cast<std::size_t>((len + step - 1) / step));

    for (std::int64_t i = 0; i < len; ++i) {
        Block& bl = blocks[static_cast<std::size_t>(i / step)];
        const char b = seq[static_cast<std::size_t>(i)];
        if (b == 'N') continue;
        ++bl.n_acgt;
        if (b == 'C') { ++bl.n_c; if (i + 1 < len && seq[static_cast<std::size_t>(i + 1)] == 'G') ++bl.n_cpg; }
        else if (b == 'G') ++bl.n_g;
    }
    for_each_site(seq, k, predictor, ctx_filter, scan_plus, scan_minus,
        [&](std::size_t i, const char* ctx, char strand, const std::string&, double prob) {
            Block& bl = blocks[i / static_cast<std::size_t>(step)];
            const int c = ctx_index(ctx);
            ++bl.n[c][strand == '+' ? 0 : 1];
            bl.psum[c] += prob;
        });

    // running sum over blocks j .. j+per_window-1: O(n_blocks) for any overlap
    const std::size_t per_window = static_cast<std::size_t>(window / step);
    Block w;
    for (std::size_t q = 0; q < per_window && q < blocks.size(); ++q) w.add(blocks[q]);
    for (std::size_t j = 0; j < blocks.size(); ++j) {
        if (j > 0) {
            w.add(blocks[j - 1], -1);
            if (j + per_window - 1 < blocks.size()) w.add(blocks[j + per_window - 1]);
        }
        const std::int64_t start = static_cast<std::int64_t>(j) * step;
        const std::int64_t end = std::min(start + window, len);

        out << id << '\t' << start << '\t' << end << '\t' << w.n_acgt << '\t' << std::fixed << std::setprecision(4);
        if (w.n_acgt > 0) out << static_cast<double>(w.n_c + w.n_g) / static_cast<double>(w.n_acgt);
        else out << "NA";
        out << '\t';
        // CpG observed/expected (Gardiner-Garden & Frommer 1987)
        if (w.n_c > 0 && w.n_g > 0)
            out << static_cast<double>(w.n_cpg) * static_cast<double>(w.n_acgt) / (static_cast<double>(w.n_c) * static_cast<double>(w.n_g));
        else out << "NA";
        for (int c = 0; c < 3; ++c) out << '\t' << w.n[c][0] << '\t' << w.n[c][1];
        out << std::setprecision(6);
        for (int c = 0; c < 3; ++c) {
            const std::int64_t nc = w.n[c][0] + w.n[c][1];
            out << '\t';
            if (nc > 0) out << w.psum[c] / static_cast<double>(nc);
            else out << "NA";
        }
        out << '\n';
    }
}

} // namespace rv

int main(int argc, char* argv[]) {
    using namespace rv;

    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options opt;
    const int cli_rc = parse_cli(argc, argv, opt);
    if (cli_rc >= 0) return cli_rc; // help/version printed, or a CLI error

    if (opt.infile.empty()) {
        std::cerr << "Error: input file required. Use -i <file>\n";
        print_help(std::cerr);
        return 1;
    }
    if (opt.step > 0 && opt.window == 0) {
        std::cerr << "Error: --step needs --window\n";
        return 1;
    }
    if (opt.window > 0) {
        if (opt.step == 0) opt.step = opt.window;
        if (opt.step > opt.window || opt.window % opt.step != 0) {
            std::cerr << "Error: --window must be a multiple of --step (and >= --step)\n";
            return 1;
        }
    }
    if (opt.heuristic && !opt.model_path.empty()) {
        std::cerr << "Error: --heuristic and --model are mutually exclusive\n";
        return 1;
    }
    if (opt.k % 2 == 0 || opt.k < 3) {
        std::cerr << "Error: k must be odd and >= 3 (e.g., -k 5)\n";
        return 1;
    }

    // "-" reads stdin directly; /dev/stdin also works now that the reader never seeks
    std::ifstream fin;
    std::istream* inp = &stdin_stream();
    if (opt.infile != "-") {
        fin.open(opt.infile);
        if (!fin) {
            std::cerr << "Error: cannot open input: " << opt.infile << "\n";
            return 2;
        }
        inp = &fin;
    }
    std::istream& in = *inp;

    // Detect file type if needed
    FileType ft = FileType::FASTA;
    if (opt.mode == "auto") {
        // Peek first non-space char (peek/get only, safe on pipes)
        while (in.peek() != EOF && std::isspace(in.peek())) in.get();
        const int c = in.peek();
        if (c == '@')      ft = FileType::FASTQ;
        else if (c == '>') ft = FileType::FASTA;
        else {
            std::cerr << "Error: input is empty or not FASTA/FASTQ (first record must start with '>' or '@')\n";
            return 2;
        }
    } else {
        ft = (opt.mode == "fastq") ? FileType::FASTQ : FileType::FASTA;
    }

    // Predictor (heuristic or model)
    Model mdl;
    Predictor predictor;

    if (!opt.model_path.empty()) {
        std::string err;
        if (!mdl.load(opt.model_path, err)) {
            std::cerr << "Error: failed to load model: " << err << "\n";
            return 2;
        }
        if (mdl.kmer_len != 0 && mdl.kmer_len != opt.k) {
            std::cerr << "Error: model k-mers have length " << mdl.kmer_len
                      << " but -k is " << opt.k << "; every site would get only BIAS\n";
            return 1;
        }
        if (mdl.w.empty() && mdl.ctx_bias.empty()) {
            std::cerr << "Warning: model has no k-mer or BIAS_<context> rows; every site scores sigmoid(BIAS)\n";
        }
        predictor = [&](const std::string& ctx, const std::string& km){ return mdl.score(ctx, km); };
    }
    if (opt.heuristic || !predictor) {
        // Typical plant prior: CG > CHG > CHH (tune as needed)
        predictor = [&](const std::string& ctx, const std::string& /*km*/){
            if (ctx == "CG")  return 0.85;
            if (ctx == "CHG") return 0.55;
            return 0.25; // CHH
        };
    }

    // Output stream
    std::ostream* out = &std::cout;
    std::ofstream fout;
    if (!opt.outfile.empty()) {
        fout.open(opt.outfile);
        if (!fout) { std::cerr << "Error: cannot open --out: " << opt.outfile << "\n"; return 3; }
        out = &fout;
    }

    if (opt.print_header) {
        if (opt.window > 0) print_profile_header(*out);
        else *out << "seq_id\tpos_1based\tcontext\tstrand\tkmer\tprob\n";
    }

    // Read and scan
    Reader rdr(in, ft);
    while (rdr.next()) {
        if (opt.window > 0)
            profile_sequence(rdr.cur_id, rdr.cur_seq, opt.k, predictor, opt.ctx_filter,
                             opt.strand != "-", opt.strand != "+", opt.window, opt.step, *out);
        else
            scan_sequence(rdr.cur_id, rdr.cur_seq, opt.k, predictor, opt.ctx_filter,
                          opt.strand != "-", opt.strand != "+", *out);
    }

    out->flush();
    if (!*out) {
        std::cerr << "Error: failed writing output (disk full or closed pipe?)\n";
        return 3;
    }
    return 0;
}
/*
  Revelio — a lightweight methylation-context predictor for FASTA/FASTQ
  Author: Yasin Kaya | MPIPZ - MPG, Cologne
  Date  : 18.08.2025
  License: MIT

  Features
  - Detect CG / CHG / CHH contexts (H ∈ {A, C, T}) around cytosines (forward strand).
  - Extract a centered k-mer (odd k, default k=5).
  - Predict methylation probability via:
      (A) heuristic priors (CG > CHG > CHH), or
      (B) logistic regression weights: prob = sigmoid(bias + w[context,kmer]).
  - Outputs TSV: seq_id  pos_1based  context  strand  kmer  prob
  - FASTA (multi-line) and FASTQ (4-line) support with auto-detection.
  - Practical CLI: --header, --contexts, --out, --version.

  Notes / Limitations
  - Forward strand only for simplicity (counts contexts as seen on + strand).
  - Ambiguity (non-ACGT) breaks context; such windows are skipped.
  - This is a fast *sequence-context prior*. It does NOT analyze raw signal.
    For signal-level methylation (e.g., ONT), see deep learning tools such as DeepSignal-plant.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cctype>
#include <cstdint>
#include <sstream>
#include <iomanip>
#include <functional>

namespace rv {

static constexpr const char* kVersion = "0.1.0";

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

// (Prepared for future both-strands mode)
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
//   CG     ACGTA  2.10
//   CHG    CCGTC  0.80
//   CHH    CATTG -0.40
struct Model {
    double bias = 0.0;
    std::unordered_map<std::string, double> w; // key = context + '\t' + kmer

    bool load(const std::string& path) {
        std::ifstream in(path);
        if (!in) return false;
        std::string line;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;

            // BIAS <number>
            if (line.rfind("BIAS", 0) == 0) {
                std::istringstream iss(line);
                std::string tag; double wt;
                if (iss >> tag >> wt) bias = wt;
                continue;
            }

            // context kmer weight
            std::istringstream iss(line);
            std::string ctx, kmer; double wt;
            if (!(iss >> ctx >> kmer >> wt)) continue;
            for (char& c : ctx)  c = up(c);
            for (char& c : kmer) c = up(c);
            w[ctx + "\t" + kmer] = wt;
        }
        return true;
    }

    double score(const std::string& ctx, const std::string& kmer) const {
        const auto it = w.find(ctx + "\t" + kmer);
        double s = bias;
        if (it != w.end()) s += it->second;
        return sigmoid(s);
    }
};

// ---------- FASTA/FASTQ reader ----------
enum class FileType { FASTA, FASTQ };

struct Reader {
    std::istream& in;
    FileType type;
    std::string cur_id, cur_seq;

    explicit Reader(std::istream& input, FileType t) : in(input), type(t) {}

    // Read the next record (id, seq). Returns false on EOF.
    bool next() {
        cur_id.clear(); cur_seq.clear();
        std::string line;

        if (type == FileType::FASTA) {
            // Find next '>' header
            while (std::getline(in, line)) {
                if (!line.empty() && line[0] == '>') { cur_id = trim(line.substr(1)); break; }
            }
            if (cur_id.empty()) return false;

            // Accumulate sequence lines until next header or EOF
            std::streampos pos;
            while (true) {
                pos = in.tellg();
                if (!std::getline(in, line)) break;
                if (!line.empty() && line[0] == '>') { in.seekg(pos); break; }
                for (char c : line) cur_seq.push_back(is_acgt(c) ? up(c) : 'N');
            }
            return true;
        } else {
            // FASTQ: @id, seq, +, qual
            if (!std::getline(in, line)) return false;
            if (line.empty() || line[0] != '@') return next(); // skip until a header
            cur_id = trim(line.substr(1));
            std::string seq, plus, qual;
            if (!std::getline(in, seq))  return false;
            if (!std::getline(in, plus)) return false;
            if (!std::getline(in, qual)) return false;
            cur_seq.reserve(seq.size());
            for (char c : seq) cur_seq.push_back(is_acgt(c) ? up(c) : 'N');
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
    std::unordered_set<std::string> ctx_filter; // empty = accept all
};

// Parse comma-separated context list like "CG,CHG,CHH"
inline std::unordered_set<std::string> parse_ctx_list(const std::string& s) {
    std::unordered_set<std::string> S;
    std::istringstream iss(s);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        for (char& c : tok) c = up(c);
        if (tok=="CG" || tok=="CHG" || tok=="CHH") S.insert(tok);
    }
    return S;
}

inline void print_help() {
    std::cerr <<
    "revelio " << kVersion << " — lightweight methylation context predictor\n"
    "Usage:\n"
    "  revelio -i <input.fa|fq> [--fasta|--fastq] [-k 5]\n"
    "          [--heuristic | --model weights.tsv]\n"
    "          [--contexts CG,CHG,CHH] [--header] [--out output.tsv]\n"
    "\n"
    "Output columns (TSV): seq_id  pos_1based  context  strand  kmer  prob\n"
    "\n"
    "Examples:\n"
    "  revelio -i genome.fa --heuristic > calls.tsv\n"
    "  revelio -i reads.fq --fastq --model weights.tsv -k 7 --header --out calls.tsv\n";
}

inline bool parse_cli(int argc, char* argv[], Options& opt) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-i" || a == "--in") {
            if (i+1 >= argc) { std::cerr << "Missing value for " << a << "\n"; return false; }
            opt.infile = argv[++i];
        } else if (a == "--fastq") {
            opt.mode = "fastq";
        } else if (a == "--fasta") {
            opt.mode = "fasta";
        } else if (a == "-k") {
            if (i+1 >= argc) { std::cerr << "Missing value for -k\n"; return false; }
            opt.k = std::stoi(argv[++i]);
        } else if (a == "--model") {
            if (i+1 >= argc) { std::cerr << "Missing value for --model\n"; return false; }
            opt.model_path = argv[++i];
        } else if (a == "--heuristic") {
            opt.heuristic = true;
        } else if (a == "--header") {
            opt.print_header = true;
        } else if (a == "--contexts") {
            if (i+1 >= argc) { std::cerr << "Missing value for --contexts\n"; return false; }
            opt.ctx_filter = parse_ctx_list(argv[++i]);
            if (opt.ctx_filter.empty()) {
                std::cerr << "Warning: --contexts parsed empty set; accepting all.\n";
            }
        } else if (a == "--out") {
            if (i+1 >= argc) { std::cerr << "Missing value for --out\n"; return false; }
            opt.outfile = argv[++i];
        } else if (a == "-h" || a == "--help") {
            print_help();
            return false; // signal caller to exit(0)
        } else if (a == "--version") {
            std::cerr << "revelio " << kVersion << "\n";
            return false; // signal caller to exit(0)
        } else {
            std::cerr << "Unknown option: " << a << "\n";
            return false;
        }
    }
    return true;
}

// Scan one sequence and write hits
inline void scan_sequence(const std::string& id,
                          const std::string& seq,
                          int k,
                          const std::function<double(const std::string&, const std::string&)>& predictor,
                          const std::unordered_set<std::string>& ctx_filter,
                          std::ostream& out)
{
    for (std::size_t i = 0; i < seq.size(); ++i) {
        const char b = seq[i];
        if (b != 'C') continue;

        std::string ctx;
        if (i + 1 < seq.size() && seq[i + 1] == 'G') {
            ctx = "CG";
        } else if (i + 2 < seq.size() && is_H(seq[i + 1]) && seq[i + 2] == 'G') {
            ctx = "CHG";
        } else if (i + 2 < seq.size() && is_H(seq[i + 1]) && is_H(seq[i + 2])) {
            ctx = "CHH";
        } else {
            continue;
        }

        if (!ctx_filter.empty() && !ctx_filter.count(ctx)) continue;

        const std::string km = centered_kmer(seq, i, k);
        if (km.empty()) continue;

        const double prob = predictor(ctx, km);

        // seq_id  pos_1based  context  strand  kmer  prob
        out << id << '\t'
            << (i + 1) << '\t'
            << ctx << '\t'
            << '+' << '\t'
            << km << '\t'
            << std::fixed << std::setprecision(6) << prob << '\n';
    }
}

} // namespace rv

int main(int argc, char* argv[]) {
    using namespace rv;

    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options opt;
    if (!parse_cli(argc, argv, opt)) return 0; // help/version already printed

    if (opt.infile.empty()) {
        std::cerr << "Error: input file required. Use -i <file>\n";
        print_help();
        return 1;
    }
    if (opt.k % 2 == 0 || opt.k < 3) {
        std::cerr << "Error: k must be odd and >= 3 (e.g., -k 5)\n";
        return 1;
    }

    std::ifstream in(opt.infile);
    if (!in) {
        std::cerr << "Error: cannot open input: " << opt.infile << "\n";
        return 2;
    }

    // Detect file type if needed
    FileType ft = FileType::FASTA;
    if (opt.mode == "auto") {
        // Peek first non-space char
        char c;
        std::streampos pos = in.tellg();
        while (in.get(c)) {
            if (!std::isspace(static_cast<unsigned char>(c))) break;
        }
        if (in && c == '@')      ft = FileType::FASTQ;
        else if (in && c == '>') ft = FileType::FASTA;
        else                     ft = FileType::FASTA; // fallback
        in.clear(); in.seekg(pos);
    } else {
        ft = (opt.mode == "fastq") ? FileType::FASTQ : FileType::FASTA;
    }

    // Predictor (heuristic or model)
    Model mdl;
    std::function<double(const std::string&, const std::string&)> predictor;

    if (!opt.heuristic && !opt.model_path.empty()) {
        if (!mdl.load(opt.model_path)) {
            std::cerr << "Warning: failed to load model; falling back to heuristic.\n";
            opt.heuristic = true;
        } else {
            predictor = [&](const std::string& ctx, const std::string& km){ return mdl.score(ctx, km); };
        }
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
    std::unique_ptr<std::ostream> out_holder;
    std::ostream* out = &std::cout;
    std::ofstream fout;
    if (!opt.outfile.empty()) {
        fout.open(opt.outfile);
        if (!fout) { std::cerr << "Error: cannot open --out: " << opt.outfile << "\n"; return 3; }
        out = &fout;
    }

    if (opt.print_header) {
        *out << "seq_id\tpos\tcontext\tstrand\tkmer\tprob\n";
    }

    // Read and scan
    Reader rdr(in, ft);
    while (rdr.next()) {
        scan_sequence(rdr.cur_id, rdr.cur_seq, opt.k, predictor, opt.ctx_filter, *out);
    }

    return 0;
}
#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

struct FMIndex {
    std::string text, bwt;
    std::vector<int> sa;
    std::array<int,256> C{};
    std::vector<std::array<int,256>> occ;

    static void counting_sort(std::vector<int>& sa, const std::vector<int>& r, int k, int maxv) {
        int n = (int)sa.size();
        std::vector<int> tmp(n), cnt(maxv + 1, 0);
        for (int i = 0; i < n; ++i) cnt[(sa[i] + k < n) ? r[sa[i] + k] : 0]++;
        for (int i = 1; i <= maxv; ++i) cnt[i] += cnt[i - 1];
        for (int i = n - 1; i >= 0; --i) {
            int idx = (sa[i] + k < n) ? r[sa[i] + k] : 0;
            tmp[--cnt[idx]] = sa[i];
        }
        sa.swap(tmp);
    }

    static std::vector<int> build_sa(const std::string& s) {
        int n = (int)s.size();
        std::vector<int> sa(n), r(n), tmp(n);
        for (int i = 0; i < n; ++i) sa[i] = i, r[i] = (unsigned char)s[i] + 1;
        int maxv = 257;
        for (int k = 1; k < n; k <<= 1) {
            counting_sort(sa, r, k, maxv);
            counting_sort(sa, r, 0, maxv);
            tmp[sa[0]] = 1;
            int classes = 1;
            for (int i = 1; i < n; ++i) {
                int a = sa[i], b = sa[i - 1];
                int ra1 = r[a], rb1 = r[b];
                int ra2 = (a + k < n) ? r[a + k] : 0;
                int rb2 = (b + k < n) ? r[b + k] : 0;
                if (ra1 != rb1 || ra2 != rb2) classes++;
                tmp[a] = classes;
            }
            r.swap(tmp);
            maxv = classes;
            if (classes == n) break;
        }
        return sa;
    }

    void build(std::string t) {
        if (t.empty()) return;
        if (t.back() != '$') t.push_back('$');
        text = std::move(t);
        sa = build_sa(text);
        int n = (int)text.size();
        bwt.resize(n);
        for (int i = 0; i < n; ++i) {
            int p = sa[i];
            bwt[i] = (p == 0) ? text[n - 1] : text[p - 1];
        }
        std::array<int,256> freq{};
        for (char c : bwt) freq[(unsigned char)c]++;
        int sum = 0;
        for (int i = 0; i < 256; ++i) C[i] = std::exchange(sum, sum + freq[i]);
        occ.resize(n + 1);
        occ[0].fill(0);
        for (int i = 0; i < n; ++i) {
            occ[i + 1] = occ[i];
            occ[i + 1][(unsigned char)bwt[i]]++;
        }
    }

    std::pair<int,int> search(const std::string& p) const {
        int l = 0, r_ = (int)bwt.size();
        for (int i = (int)p.size() - 1; i >= 0; --i) {
            unsigned char c = (unsigned char)p[i];
            l = C[c] + occ[l][c];
            r_ = C[c] + occ[r_][c];
            if (l >= r_) return {0,0};
        }
        return {l,r_};
    }
};

inline std::string read_fasta(const std::string& path) {
    std::ifstream in(path);
    std::string line, s;
    while (std::getline(in, line))
        if (!line.empty() && line[0] != '>') s += line;
    return s;
}

inline char comp(char c) {
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return 'N';
    }
}

inline std::string revcomp(const std::string& s) {
    std::string r; r.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) r.push_back(comp(*it));
    return r;
}

struct Seed { size_t pos; std::string seq; };
struct Alignment { bool accepted{false}; int score{0}; double identity{0.0}; int ref_start{0}; int ref_end{0}; int aligned_len{0}; };
struct ReadResult { bool mapped{false}; bool unique{false}; bool multi{false}; Alignment best{}; };

inline std::vector<Seed> collect_seeds(const std::string& s, size_t min_len) {
    std::vector<Seed> seeds; size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == 'N') {
            if (i > start && i - start >= min_len) seeds.push_back({start, s.substr(start, i - start)});
            start = i + 1;
        }
    }
    return seeds;
}

inline bool verify_at_start(const std::string& ref, const std::string& read, int pos, int allowed_mismatches) {
    if (pos < 0 || static_cast<size_t>(pos) + read.size() > ref.size()) return false;
    int mismatches = 0;
    for (size_t i = 0; i < read.size(); ++i) {
        if (read[i] != ref[pos + static_cast<int>(i)] && ++mismatches > allowed_mismatches) return false;
    }
    return true;
}

Alignment smith_waterman(const std::string& read, const std::string& ref, int window_start, int window_end) {
    constexpr int MATCH = 2, MISMATCH = -2, GAP = -3;
    size_t n = read.size(), m = static_cast<size_t>(window_end - window_start);
    static thread_local std::vector<int> prev, curr;
    static thread_local std::vector<uint8_t> trace;
    if (prev.size() < m + 1) prev.assign(m + 1, 0); else std::fill(prev.begin(), prev.end(), 0);
    if (curr.size() < m + 1) curr.assign(m + 1, 0); else std::fill(curr.begin(), curr.end(), 0);
    if (trace.size() < (n + 1) * (m + 1)) trace.assign((n + 1) * (m + 1), 0);
    else std::fill(trace.begin(), trace.begin() + (n + 1) * (m + 1), 0);

    int best_score = 0; size_t best_i = 0, best_j = 0;
    for (size_t i = 1; i <= n; ++i) {
        curr[0] = 0;
        for (size_t j = 1; j <= m; ++j) {
            char ref_c = ref[window_start + static_cast<int>(j - 1)];
            int diag = prev[j - 1] + ((read[i - 1] == ref_c) ? MATCH : MISMATCH);
            int up = prev[j] + GAP;
            int left = curr[j - 1] + GAP;
            int score = std::max({0, diag, up, left});
            uint8_t dir = (score == diag && score) ? 1 : (score == up && score) ? 2 : (score == left && score) ? 3 : 0;
            curr[j] = score; trace[i * (m + 1) + j] = dir;
            if (score > best_score) best_score = score, best_i = i, best_j = j;
        }
        prev.swap(curr);
    }

    Alignment res; res.score = best_score; if (!best_score) return res;
    size_t i = best_i, j = best_j, matches = 0; int aligned_len = 0; size_t ref_end = j;
    while (i && j) {
        uint8_t dir = trace[i * (m + 1) + j]; if (!dir) break;
        if (dir == 1) { aligned_len++; if (read[i - 1] == ref[window_start + static_cast<int>(j - 1)]) matches++; --i; --j; }
        else if (dir == 2) { aligned_len++; --i; }
        else { aligned_len++; --j; }
    }
    res.accepted = true; res.aligned_len = aligned_len;
    res.ref_start = window_start + static_cast<int>(j);
    res.ref_end = window_start + static_cast<int>(ref_end);
    res.identity = aligned_len ? static_cast<double>(matches) / aligned_len : 0.0;
    return res;
}

ReadResult map_read(const FMIndex& fm, const std::string& ref, const std::string& read, size_t min_seed_len) {
    constexpr int MAX_CANDIDATES = 200, MARGIN = 30;
    constexpr double ID_THRESH = 0.90, LEN_THRESH = 0.7;
    const int ref_len = static_cast<int>(ref.size());
    std::vector<Alignment> accepted;

    auto try_orientation = [&](const std::string& oriented) {
        auto seeds = collect_seeds(oriented, min_seed_len);
        if (seeds.empty()) return;
        std::sort(seeds.begin(), seeds.end(), [](const Seed& a, const Seed& b) { return a.seq.size() > b.seq.size(); });
        int allowed_mismatches = static_cast<int>(std::floor((1.0 - ID_THRESH) * oriented.size()));
        for (const auto& seed : seeds) {
            auto iv = fm.search(seed.seq);
            int hits = iv.second - iv.first; if (!hits) continue;
            std::vector<int> candidates; candidates.reserve(std::min(hits, MAX_CANDIDATES));
            if (hits <= MAX_CANDIDATES) {
                for (int i = iv.first; i < iv.second && (int)candidates.size() < MAX_CANDIDATES; ++i) {
                    int pos = fm.sa[i] - (int)seed.pos; if (pos < 0 || pos >= ref_len) continue; candidates.push_back(pos);
                }
            } else {
                for (int k = 0; k < MAX_CANDIDATES; ++k) {
                    double frac = (MAX_CANDIDATES == 1) ? 0.0 : (double)k / (double)(MAX_CANDIDATES - 1);
                    int idx = iv.first + (int)(frac * (hits - 1) + 0.5); if (idx >= iv.second) idx = iv.second - 1;
                    int pos = fm.sa[idx] - (int)seed.pos; if (pos < 0 || pos >= ref_len) continue; candidates.push_back(pos);
                }
            }
            if (candidates.empty()) continue;
            for (int pos : candidates) {
                if (verify_at_start(ref, oriented, pos, allowed_mismatches)) {
                    int mismatches = 0;
                    for (size_t i = 0; i < oriented.size(); ++i) if (oriented[i] != ref[pos + (int)i]) mismatches++;
                    int matches = (int)oriented.size() - mismatches;
                    Alignment aln{true, matches * 2 - mismatches * 2, oriented.empty() ? 0.0 : (double)matches / oriented.size(), pos, pos + (int)oriented.size(), (int)oriented.size()};
                    accepted.push_back(aln); continue;
                }
                int ws = std::max(0, pos - MARGIN), we = std::min(ref_len, pos + (int)oriented.size() + MARGIN);
                if (ws >= we) continue;
                Alignment aln = smith_waterman(oriented, ref, ws, we);
                if (aln.aligned_len >= (int)(LEN_THRESH * oriented.size()) && aln.identity >= ID_THRESH) accepted.push_back(aln);
            }
            if (!accepted.empty()) break;
        }
    };

    try_orientation(read);
    try_orientation(revcomp(read));
    ReadResult res;
    if (accepted.empty()) return res;
    std::sort(accepted.begin(), accepted.end(), [](const Alignment& a, const Alignment& b) { return a.score > b.score; });
    res.mapped = true; res.best = accepted.front();
    int best = accepted[0].score, second = (accepted.size() > 1) ? accepted[1].score : std::numeric_limits<int>::min();
    res.unique = (accepted.size() == 1 || best - second >= 10);
    res.multi = !res.unique;
    return res;
}

int main() {
    const size_t MIN_SEED = 20, BATCH = 200000;
    std::cout << "[INFO] Reading reference\n" << std::flush;
    std::string ref = read_fasta("GCF_000005845.2_ASM584v2_genomic.fna");
    size_t ref_len = ref.size();

    std::cout << "[INFO] Building FM-index\n" << std::flush;
    FMIndex fm; fm.build(ref);
    std::cout << "[INFO] FM-index built\n" << std::flush;

    std::ifstream fq("ERR022075_1.fastq");
    if (!fq) { std::cerr << "[ERROR] FASTQ not found (run from project root)\n"; return 1; }

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
    std::cout << "[INFO] OpenMP threads: " << num_threads << "\n" << std::flush;
#else
    std::cout << "[INFO] OpenMP not enabled\n" << std::flush;
#endif

    std::vector<std::vector<uint32_t>> thread_cov(num_threads, std::vector<uint32_t>(ref_len, 0));
    size_t total = 0, mapped = 0, unique = 0, multi = 0;

    std::string id, seq, plus, qual;
    std::vector<std::string> batch; batch.reserve(BATCH);
    std::cout << "[INFO] FASTQ streaming started\n" << std::flush;

    while (true) {
        batch.clear();
        for (size_t i = 0; i < BATCH; ++i) {
            if (!std::getline(fq, id)) break;
            if (!std::getline(fq, seq)) break;
            if (!std::getline(fq, plus)) break;
            if (!std::getline(fq, qual)) break;
            batch.push_back(seq);
        }
        if (batch.empty()) break;
        size_t batch_n = batch.size(); total += batch_n;

        struct ThreadStats { size_t mapped{0}, unique{0}, multi{0}; };
        std::vector<ThreadStats> thread_stats(num_threads);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 256)
#endif
        for (size_t i = 0; i < batch_n; ++i) {
            int tid =
#ifdef _OPENMP
                omp_get_thread_num();
#else
                0;
#endif
            auto res = map_read(fm, ref, batch[i], MIN_SEED);
            if (res.mapped) {
                auto& ts = thread_stats[tid];
                ts.mapped++; res.unique ? ts.unique++ : ts.multi++;
                auto& cov = thread_cov[tid];
                int start = std::max(0, res.best.ref_start), end = std::min((int)ref_len, res.best.ref_end);
                for (int p = start; p < end; ++p) cov[p]++;
            }

            size_t done = (total - batch_n) + i + 1;
            if (done % 500000 == 0) {
                double pct = total ? 100.0 * (double)done / (double)total : 0.0;
                std::cout << "[INFO LOG] Processed: " << done << " (" << pct << " %)\n" << std::flush;
            }
        }
        for (const auto& ts : thread_stats) { mapped += ts.mapped; unique += ts.unique; multi += ts.multi; }
    }

    std::vector<uint32_t> coverage(ref_len, 0);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < (int)ref_len; ++i) {
        uint32_t sum = 0; for (const auto& cov : thread_cov) sum += cov[i]; coverage[i] = sum;
    }

    uint64_t cov_sum = 0; size_t cov_ge1 = 0, cov_ge10 = 0;
    for (auto v : coverage) { cov_sum += v; if (v >= 1) cov_ge1++; if (v >= 10) cov_ge10++; }
    double mean_cov = ref_len ? (double)cov_sum / ref_len : 0.0;
    double cov1_pct = ref_len ? 100.0 * (double)cov_ge1 / ref_len : 0.0;
    double cov10_pct = ref_len ? 100.0 * (double)cov_ge10 / ref_len : 0.0;
    double mapping_rate = total ? 100.0 * (double)mapped / total : 0.0;
    double unique_pct = total ? 100.0 * (double)unique / total : 0.0;
    double multi_pct = total ? 100.0 * (double)multi / total : 0.0;
    size_t unmapped = total - mapped;

    std::cout << "[INFO] Mapping finished\n" << std::flush;
    std::cout << "Total reads: " << total << "\n";
    std::cout << "Mapped reads: " << mapped << " (" << mapping_rate << "%)\n";
    std::cout << "Unique reads: " << unique << " (" << unique_pct << "%)\n";
    std::cout << "Multi-mapped reads: " << multi << " (" << multi_pct << "%)\n";
    std::cout << "Mean coverage: " << mean_cov << "\n";
    std::cout << "Genome covered >=1x: " << cov1_pct << "%\n";
    std::cout << "Genome covered >=10x: " << cov10_pct << "%\n";
    std::cout << "Unmapped reads: " << unmapped << "\n";
}

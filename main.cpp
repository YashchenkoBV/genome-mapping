#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

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

        for (int i = 0; i < n; ++i) {
            int idx = (sa[i] + k < n) ? r[sa[i] + k] : 0;
            cnt[idx]++;
        }
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

        for (int i = 0; i < n; ++i) sa[i] = i;
        for (int i = 0; i < n; ++i) r[i] = (unsigned char)s[i] + 1;

        int maxv = 256 + 1;

        for (int k = 1, round = 0; k < n; k <<= 1, ++round) {
            if ((round % 1) == 0) {
                std::cout << "[INFO] SA round " << round << " (k=" << k << ")\n";
            }

            counting_sort(sa, r, k, maxv);
            counting_sort(sa, r, 0, maxv);

            tmp[sa[0]] = 1;
            int classes = 1;

            for (int i = 1; i < n; ++i) {
                int a = sa[i], b = sa[i - 1];
                int ra1 = r[a];
                int rb1 = r[b];
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

        std::cout << "[INFO] Building suffix array (n=" << text.size() << ")\n";
        sa = build_sa(text);
        std::cout << "[INFO] Suffix array built\n";

        int n = (int)text.size();
        bwt.resize(n);
        for (int i = 0; i < n; ++i) {
            int p = sa[i];
            bwt[i] = (p == 0) ? text[n - 1] : text[p - 1];
        }
        std::cout << "[INFO] BWT built\n";

        std::array<int,256> freq{};
        for (char c : bwt) freq[(unsigned char)c]++;

        int sum = 0;
        for (int i = 0; i < 256; ++i) {
            C[i] = sum;
            sum += freq[i];
        }

        occ.resize(n + 1);
        occ[0].fill(0);
        for (int i = 0; i < n; ++i) {
            occ[i + 1] = occ[i];
            occ[i + 1][(unsigned char)bwt[i]]++;
            if ((i + 1) % 1000000 == 0) {
                std::cout << "[INFO] Occ built: " << (i + 1) << "/" << n << "\n";
            }
        }
        std::cout << "[INFO] FM-index ready\n";
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

std::string read_fasta(const std::string& path) {
    std::ifstream in(path);
    std::string line, s;
    while (std::getline(in, line))
        if (!line.empty() && line[0] != '>')
            s += line;
    return s;
}

inline char comp(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

std::string revcomp(const std::string& s) {
    std::string r;
    r.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) r.push_back(comp(*it));
    return r;
}

struct Seed {
    size_t pos{};
    std::string seq;
};

enum class FailReason { NONE, NO_SEED, TOO_MANY_HITS, IDENTITY_TOO_LOW, ALIGNED_TOO_SHORT };

struct Alignment {
    bool accepted{false};
    int score{0};
    double identity{0.0};
    int ref_start{0};
    int ref_end{0}; // exclusive
    int aligned_len{0};
};

struct ReadResult {
    bool mapped{false};
    bool unique{false};
    bool multi{false};
    Alignment best{};
    FailReason reason{FailReason::NONE};
    size_t sw_attempts{0};
    size_t sw_success{0};
    size_t verify_success{0};
};

std::vector<Seed> collect_seeds(const std::string& s, size_t min_len) {
    std::vector<Seed> seeds;
    size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == 'N') {
            size_t len = i - start;
            if (len >= min_len) seeds.push_back({start, s.substr(start, len)});
            start = i + 1;
        }
    }
    return seeds;
}

bool verify_at_start(const std::string& ref, const std::string& read, int pos, int allowed_mismatches) {
    if (pos < 0) return false;
    if (static_cast<size_t>(pos) + read.size() > ref.size()) return false;
    int mismatches = 0;
    for (size_t i = 0; i < read.size(); ++i) {
        if (read[i] != ref[pos + static_cast<int>(i)]) {
            mismatches++;
            if (mismatches > allowed_mismatches) return false;
        }
    }
    return true;
}

// Local Smith-Waterman alignment; returns best local score and the aligned reference span.
Alignment smith_waterman(const std::string& read, const std::string& ref, int window_start, int window_end) {
    constexpr int MATCH = 2;
    constexpr int MISMATCH = -2;
    constexpr int GAP = -3;

    size_t n = read.size();
    size_t m = static_cast<size_t>(window_end - window_start);

    static thread_local std::vector<int> prev;
    static thread_local std::vector<int> curr;
    static thread_local std::vector<uint8_t> trace;
    if (prev.size() < m + 1) prev.assign(m + 1, 0);
    else std::fill(prev.begin(), prev.end(), 0);
    if (curr.size() < m + 1) curr.assign(m + 1, 0);
    else std::fill(curr.begin(), curr.end(), 0);
    if (trace.size() < (n + 1) * (m + 1)) trace.assign((n + 1) * (m + 1), 0);
    else std::fill(trace.begin(), trace.begin() + (n + 1) * (m + 1), 0);

    int best_score = 0;
    size_t best_i = 0, best_j = 0;

    for (size_t i = 1; i <= n; ++i) {
        curr[0] = 0;
        for (size_t j = 1; j <= m; ++j) {
            char ref_c = ref[window_start + static_cast<int>(j - 1)];
            int diag = prev[j - 1] + ((read[i - 1] == ref_c) ? MATCH : MISMATCH);
            int up = prev[j] + GAP;
            int left = curr[j - 1] + GAP;
            int score = std::max({0, diag, up, left});

            uint8_t dir = 0;
            if (score == diag && score != 0) dir = 1;
            else if (score == up && score != 0) dir = 2;
            else if (score == left && score != 0) dir = 3;

            curr[j] = score;
            trace[i * (m + 1) + j] = dir;

            if (score > best_score) {
                best_score = score;
                best_i = i;
                best_j = j;
            }
        }
        prev.swap(curr);
    }

    Alignment res;
    res.score = best_score;
    if (best_score == 0) return res;

    size_t i = best_i;
    size_t j = best_j;
    size_t matches = 0;
    int aligned_len = 0;
    size_t ref_end = j;

    while (i > 0 && j > 0) {
        uint8_t dir = trace[i * (m + 1) + j];
        if (dir == 0) break;
        if (dir == 1) {
            aligned_len++;
            if (read[i - 1] == ref[window_start + static_cast<int>(j - 1)]) matches++;
            --i; --j;
        } else if (dir == 2) {
            aligned_len++;
            --i;
        } else if (dir == 3) {
            aligned_len++;
            --j;
        }
    }

    size_t ref_start = j;
    res.accepted = true;
    res.aligned_len = aligned_len;
    res.ref_start = window_start + static_cast<int>(ref_start);
    res.ref_end = window_start + static_cast<int>(ref_end);
    res.identity = (aligned_len > 0) ? static_cast<double>(matches) / static_cast<double>(aligned_len) : 0.0;
    return res;
}

ReadResult map_read(const FMIndex& fm, const std::string& ref, const std::string& read, size_t min_seed_len) {
    constexpr int MAX_CANDIDATES = 200;
    constexpr int MARGIN = 30;
    constexpr double ID_THRESH = 0.90;
    constexpr double LEN_THRESH = 0.7;
    const int ref_len = static_cast<int>(ref.size());

    std::vector<Alignment> accepted;
    bool had_seeds = false;
    bool had_nonrepetitive_seed = false;
    bool attempted_alignment = false;
    bool saw_identity_fail = false;
    bool saw_len_fail = false;
    size_t sw_attempts = 0;
    size_t sw_success = 0;
    size_t verify_success = 0;

    // Try one orientation: pick seeds from longest to shortest, generate candidates, extend.
    auto try_orientation = [&](const std::string& oriented) {
        auto seeds = collect_seeds(oriented, min_seed_len);
        if (!seeds.empty()) had_seeds = true;
        if (seeds.empty()) return;
        std::sort(seeds.begin(), seeds.end(), [](const Seed& a, const Seed& b) {
            return a.seq.size() > b.seq.size();
        });
        int allowed_mismatches = static_cast<int>(std::floor((1.0 - ID_THRESH) * oriented.size()));

        for (const auto& seed : seeds) {
            auto iv = fm.search(seed.seq);
            int hits = iv.second - iv.first;
            if (hits == 0) continue;

            std::vector<int> candidates;
            candidates.reserve(std::min(hits, MAX_CANDIDATES));
            if (hits <= MAX_CANDIDATES) {
                for (int i = iv.first; i < iv.second && static_cast<int>(candidates.size()) < MAX_CANDIDATES; ++i) {
                    int pos = fm.sa[i] - static_cast<int>(seed.pos);
                    if (pos < 0 || pos >= ref_len) continue;
                    candidates.push_back(pos);
                }
            } else {
                // Evenly sample across the full interval (inclusive endpoints).
                for (int k = 0; k < MAX_CANDIDATES; ++k) {
                    double frac = (MAX_CANDIDATES == 1) ? 0.0
                        : static_cast<double>(k) / static_cast<double>(MAX_CANDIDATES - 1);
                    int offset = static_cast<int>(frac * (hits - 1) + 0.5);
                    int idx = iv.first + offset;
                    if (idx >= iv.second) idx = iv.second - 1;
                    int pos = fm.sa[idx] - static_cast<int>(seed.pos);
                    if (pos < 0 || pos >= ref_len) continue;
                    candidates.push_back(pos);
                }
            }
            if (!candidates.empty()) had_nonrepetitive_seed = true;
            if (candidates.empty()) continue;

            for (int pos : candidates) {
                if (verify_at_start(ref, oriented, pos, allowed_mismatches)) {
                    int mismatches = 0;
                    for (size_t i = 0; i < oriented.size(); ++i) {
                        if (oriented[i] != ref[pos + static_cast<int>(i)]) mismatches++;
                    }
                    int matches = static_cast<int>(oriented.size()) - mismatches;
                    Alignment aln;
                    aln.accepted = true;
                    aln.aligned_len = static_cast<int>(oriented.size());
                    aln.ref_start = pos;
                    aln.ref_end = pos + static_cast<int>(oriented.size());
                    aln.identity = oriented.empty() ? 0.0 : static_cast<double>(matches) / static_cast<double>(oriented.size());
                    aln.score = matches * 2 - mismatches * 2;
                    accepted.push_back(aln);
                    verify_success++;
                    continue;
                }

                int window_start = std::max(0, pos - MARGIN);
                int window_end = std::min(ref_len, pos + static_cast<int>(oriented.size()) + MARGIN);
                if (window_start >= window_end) continue;
                attempted_alignment = true;
                sw_attempts++;
                Alignment aln = smith_waterman(oriented, ref, window_start, window_end);
                bool len_ok = aln.aligned_len >= static_cast<int>(LEN_THRESH * oriented.size());
                bool id_ok = aln.identity >= ID_THRESH;
                if (len_ok && id_ok) {
                    accepted.push_back(aln);
                    sw_success++;
                } else {
                    if (!len_ok) saw_len_fail = true;
                    if (!id_ok) saw_identity_fail = true;
                }
            }
            if (!accepted.empty()) break; // fallback to next-longest only if no hits
        }
    };

    try_orientation(read);
    try_orientation(revcomp(read));

    ReadResult result;
    result.sw_attempts = sw_attempts;
    result.sw_success = sw_success;
    result.verify_success = verify_success;
    if (accepted.empty()) {
        if (!had_seeds) result.reason = FailReason::NO_SEED;
        else if (!had_nonrepetitive_seed) result.reason = FailReason::TOO_MANY_HITS;
        else if (attempted_alignment && saw_identity_fail) result.reason = FailReason::IDENTITY_TOO_LOW;
        else if (attempted_alignment && saw_len_fail) result.reason = FailReason::ALIGNED_TOO_SHORT;
        return result;
    }

    std::sort(accepted.begin(), accepted.end(), [](const Alignment& a, const Alignment& b) {
        return a.score > b.score;
    });

    result.mapped = true;
    result.best = accepted.front();

    int best_score = accepted[0].score;
    int second_score = (accepted.size() > 1) ? accepted[1].score : std::numeric_limits<int>::min();

    if (accepted.size() == 1 || best_score - second_score >= 10) {
        result.unique = true;
    } else {
        result.multi = true;
    }

    return result;
}

int main() {
    const size_t MIN_SEED = 20;
    const size_t LOG_STEP = 500000;
    const size_t BATCH = 200000;

    std::cout << "111\n";
    std::cout << "[INFO] Reading reference\n" << std::flush;
    std::string ref = read_fasta("GCF_000005845.2_ASM584v2_genomic.fna");
    size_t ref_len = ref.size();

    std::cout << "[INFO] Building FM-index\n" << std::flush;
    FMIndex fm;
    fm.build(ref);
    std::cout << "[INFO] FM-index built\n" << std::flush;

    std::ifstream fq("ERR022075_1.fastq");
    if (!fq) {
        std::cerr << "[ERROR] FASTQ not found (run from project root)\n";
        return 1;
    }

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
    std::cout << "[INFO] OpenMP threads: " << num_threads << "\n" << std::flush;
#else
    std::cout << "[INFO] OpenMP not enabled\n" << std::flush;
#endif

    std::vector<std::vector<uint32_t>> thread_cov(num_threads, std::vector<uint32_t>(ref_len, 0));

    size_t total = 0, mapped = 0, unique = 0, multi = 0;
    size_t no_seed = 0, too_many_hits = 0, identity_low = 0, aligned_short = 0;
    size_t sw_attempts_total = 0, sw_success_total = 0, verify_success_total = 0;

    struct StatAgg {
        double score_sum{0.0};
        int score_min{std::numeric_limits<int>::max()};
        int score_max{std::numeric_limits<int>::min()};
        double id_sum{0.0};
        double id_min{std::numeric_limits<double>::infinity()};
        double id_max{0.0};
        size_t count{0};
    } stats;

    std::string id, seq, plus, qual;
    std::vector<std::string> batch;
    batch.reserve(BATCH);

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

        size_t batch_n = batch.size();
        total += batch_n;

        struct ThreadStats {
            size_t mapped{0};
            size_t unique{0};
            size_t multi{0};
            double score_sum{0.0};
            int score_min{std::numeric_limits<int>::max()};
            int score_max{std::numeric_limits<int>::min()};
            double id_sum{0.0};
            double id_min{std::numeric_limits<double>::infinity()};
            double id_max{0.0};
            size_t count{0};
            size_t no_seed{0};
            size_t too_many_hits{0};
            size_t identity_low{0};
            size_t aligned_short{0};
            size_t sw_attempts{0};
            size_t sw_success{0};
            size_t verify_success{0};
        };

        std::vector<ThreadStats> thread_stats(num_threads);
        size_t start_total = total - batch_n;

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
            {
                auto& ts = thread_stats[tid];
                ts.sw_attempts += res.sw_attempts;
                ts.sw_success += res.sw_success;
                ts.verify_success += res.verify_success;
            }
            if (res.mapped) {
                auto& ts = thread_stats[tid];
                ts.mapped++;
                if (res.unique) ts.unique++;
                else if (res.multi) ts.multi++;

                if (res.best.accepted) {
                    auto& cov = thread_cov[tid];
                    int start = std::max(0, res.best.ref_start);
                    int end = std::min(static_cast<int>(ref_len), res.best.ref_end);
                    for (int p = start; p < end; ++p) cov[p]++;
                }

                ts.score_sum += res.best.score;
                ts.score_min = std::min(ts.score_min, res.best.score);
                ts.score_max = std::max(ts.score_max, res.best.score);
                ts.id_sum += res.best.identity;
                ts.id_min = std::min(ts.id_min, res.best.identity);
                ts.id_max = std::max(ts.id_max, res.best.identity);
                ts.count++;
            } else {
                auto& ts = thread_stats[tid];
                switch (res.reason) {
                    case FailReason::NO_SEED: ts.no_seed++; break;
                    case FailReason::TOO_MANY_HITS: ts.too_many_hits++; break;
                    case FailReason::IDENTITY_TOO_LOW: ts.identity_low++; break;
                    case FailReason::ALIGNED_TOO_SHORT: ts.aligned_short++; break;
                    default: break;
                }
            }

            size_t done = start_total + i + 1;
            if (done % LOG_STEP == 0) {
                double pct_done = total ? (100.0 * static_cast<double>(done) / static_cast<double>(total)) : 0.0;
                std::cout << "[INFO LOG] Processed: " << done << " (" << pct_done << " %)\n" << std::flush;
            }
        }

        for (const auto& ts : thread_stats) {
            mapped += ts.mapped;
            unique += ts.unique;
            multi += ts.multi;
            if (ts.count > 0) {
                stats.count += ts.count;
                stats.score_sum += ts.score_sum;
                stats.score_min = std::min(stats.score_min, ts.score_min);
                stats.score_max = std::max(stats.score_max, ts.score_max);
                stats.id_sum += ts.id_sum;
                stats.id_min = std::min(stats.id_min, ts.id_min);
                stats.id_max = std::max(stats.id_max, ts.id_max);
            }
            no_seed += ts.no_seed;
            too_many_hits += ts.too_many_hits;
            identity_low += ts.identity_low;
            aligned_short += ts.aligned_short;
            sw_attempts_total += ts.sw_attempts;
            sw_success_total += ts.sw_success;
            verify_success_total += ts.verify_success;
        }
    }

    std::vector<uint32_t> coverage(ref_len, 0);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < static_cast<int>(ref_len); ++i) {
        uint32_t sum = 0;
        for (const auto& cov : thread_cov) sum += cov[i];
        coverage[i] = sum;
    }

    uint64_t cov_sum = 0;
    size_t cov_ge1 = 0, cov_ge10 = 0;
    for (auto v : coverage) {
        cov_sum += v;
        if (v >= 1) cov_ge1++;
        if (v >= 10) cov_ge10++;
    }

    double mean_cov = ref_len ? static_cast<double>(cov_sum) / static_cast<double>(ref_len) : 0.0;
    double cov1_pct = ref_len ? 100.0 * static_cast<double>(cov_ge1) / static_cast<double>(ref_len) : 0.0;
    double cov10_pct = ref_len ? 100.0 * static_cast<double>(cov_ge10) / static_cast<double>(ref_len) : 0.0;

    double mapping_rate = total ? 100.0 * static_cast<double>(mapped) / static_cast<double>(total) : 0.0;
    double unique_pct = total ? 100.0 * static_cast<double>(unique) / static_cast<double>(total) : 0.0;
    double multi_pct = total ? 100.0 * static_cast<double>(multi) / static_cast<double>(total) : 0.0;

    double score_avg = stats.count ? stats.score_sum / static_cast<double>(stats.count) : 0.0;
    double id_avg = stats.count ? stats.id_sum / static_cast<double>(stats.count) : 0.0;
    int score_min = (stats.count ? stats.score_min : 0);
    int score_max = (stats.count ? stats.score_max : 0);
    double id_min = (stats.count ? stats.id_min : 0.0);
    double id_max = (stats.count ? stats.id_max : 0.0);
    size_t unmapped = total - mapped;
    auto pct = [&](size_t v) {
        return (unmapped ? 100.0 * static_cast<double>(v) / static_cast<double>(unmapped) : 0.0);
    };
    auto pct_total = [&](size_t v) {
        return (total ? 100.0 * static_cast<double>(v) / static_cast<double>(total) : 0.0);
    };
    auto pct_mapped = [&](size_t v) {
        return (mapped ? 100.0 * static_cast<double>(v) / static_cast<double>(mapped) : 0.0);
    };

    std::cout << "[INFO] Mapping finished\n" << std::flush;
    std::cout << "Total reads: " << total << "\n";
    std::cout << "Mapped reads: " << mapped << " (" << mapping_rate << "%)\n";
    std::cout << "Unique reads: " << unique << " (" << unique_pct << "%)\n";
    std::cout << "Multi-mapped reads: " << multi << " (" << multi_pct << "%)\n";
    std::cout << "Mean coverage: " << mean_cov << "\n";
    std::cout << "Genome covered >=1x: " << cov1_pct << "%\n";
    std::cout << "Genome covered >=10x: " << cov10_pct << "%\n";
    std::cout << "Alignment score (min/avg/max): " << score_min << " / " << score_avg << " / " << score_max << "\n";
    std::cout << "Identity (min/avg/max): " << id_min << " / " << id_avg << " / " << id_max << "\n";
    std::cout << "SW attempts: " << sw_attempts_total << " (" << pct_total(sw_attempts_total) << "% of reads)\n";
    std::cout << "SW accepted: " << sw_success_total << " (" << pct_mapped(sw_success_total) << "% of mapped)\n";
    std::cout << "Verify accepted: " << verify_success_total << " (" << pct_mapped(verify_success_total) << "% of mapped)\n";
    std::cout << "Unmapped reads: " << unmapped << " (100%)\n";
    std::cout << "  no_seed: " << no_seed << " (" << pct(no_seed) << "%)\n";
    std::cout << "  too_many_hits: " << too_many_hits << " (" << pct(too_many_hits) << "%)\n";
    std::cout << "  identity_too_low: " << identity_low << " (" << pct(identity_low) << "%)\n";
    std::cout << "  aligned_too_short: " << aligned_short << " (" << pct(aligned_short) << "%)\n";
}

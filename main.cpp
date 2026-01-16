#include <algorithm>
#include <atomic>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "fm_index.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

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

// Local Smith-Waterman alignment; returns best local score and the aligned reference span.
Alignment smith_waterman(const std::string& read, const std::string& window, int window_start) {
    constexpr int MATCH = 2;
    constexpr int MISMATCH = -2;
    constexpr int GAP = -3;

    size_t n = read.size();
    size_t m = window.size();

    std::vector<int> prev(m + 1, 0), curr(m + 1, 0);
    std::vector<uint8_t> trace((n + 1) * (m + 1), 0);

    int best_score = 0;
    size_t best_i = 0, best_j = 0;

    for (size_t i = 1; i <= n; ++i) {
        curr[0] = 0;
        for (size_t j = 1; j <= m; ++j) {
            int diag = prev[j - 1] + ((read[i - 1] == window[j - 1]) ? MATCH : MISMATCH);
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
            if (read[i - 1] == window[j - 1]) matches++;
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
    const int ref_len = static_cast<int>(ref.size());

    std::vector<Alignment> accepted;
    bool had_seeds = false;
    bool had_nonrepetitive_seed = false;
    bool attempted_alignment = false;
    bool saw_identity_fail = false;
    bool saw_len_fail = false;

    // Try one orientation: pick longest seed with <=MAX_CANDIDATES hits, then extend.
    auto try_orientation = [&](const std::string& oriented) {
        auto seeds = collect_seeds(oriented, min_seed_len);
        if (!seeds.empty()) had_seeds = true;
        if (seeds.empty()) return;
        std::sort(seeds.begin(), seeds.end(), [](const Seed& a, const Seed& b) {
            return a.seq.size() > b.seq.size();
        });

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
                int window_start = std::max(0, pos - MARGIN);
                int window_end = std::min(ref_len, pos + static_cast<int>(oriented.size()) + MARGIN);
                if (window_start >= window_end) continue;
                attempted_alignment = true;
                Alignment aln = smith_waterman(oriented, ref.substr(window_start, window_end - window_start), window_start);
                bool len_ok = aln.aligned_len >= static_cast<int>(0.7 * oriented.size());
                bool id_ok = aln.identity >= 0.90;
                if (len_ok && id_ok) {
                    accepted.push_back(aln);
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
    const size_t LOG_STEP = 100000;
    const size_t BATCH = 200000;

    std::cout << "!!!\n";
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

        std::cout << "[INFO] Loaded reads: " << total << "\n" << std::flush;

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
                std::cout << "[INFO] Reads processed: " << done << "\n" << std::flush;
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
    std::cout << "Unmapped reads: " << unmapped << " (100%)\n";
    std::cout << "  no_seed: " << no_seed << " (" << pct(no_seed) << "%)\n";
    std::cout << "  too_many_hits: " << too_many_hits << " (" << pct(too_many_hits) << "%)\n";
    std::cout << "  identity_too_low: " << identity_low << " (" << pct(identity_low) << "%)\n";
    std::cout << "  aligned_too_short: " << aligned_short << " (" << pct(aligned_short) << "%)\n";
}

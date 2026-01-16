#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// --- Mapping parameters (simple, student-level constants) ---
static const int KMER = 20;
static const int STEP = 10;

static const int MAX_PER_SEED = 80;       // sample this many SA hits per seed when repetitive
static const int MAX_CANDIDATES = 200;    // max candidate starts tested per orientation
static const int MARGIN = 60;             // SW window padding around candidate

static const double VERIFY_ID_THRESH = 0.90; // fast verify mismatch allowance (Hamming)
static const double SW_ID_THRESH = 0.90;     // stricter SW identity threshold
static const double SW_LEN_THRESH = 0.70;    // stricter SW min aligned length fraction

// Low-complexity filters
static const double MAX_BASE_FRAC = 0.80;    // reject if any base >80% of non-N
static const double ENTROPY_THRESH = 0.80;   // reject if Shannon entropy < 0.80 bits

struct FMIndex {
    string text, bwt;
    vector<int> sa;
    vector<int> C;                  // prefix counts, size 256
    vector<vector<int>> occ;        // occ[i][c] counts c in bwt[0..i-1]

    static void counting_sort(vector<int>& sa, const vector<int>& r, int k, int maxv) {
        int n = (int)sa.size();
        vector<int> tmp(n), cnt(maxv + 1, 0);
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

    static vector<int> build_sa(const string& s) {
        int n = (int)s.size();
        vector<int> sa(n), r(n), tmp(n);
        for (int i = 0; i < n; ++i) {
            sa[i] = i;
            r[i] = (int)(unsigned char)s[i] + 1;
        }
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

    void build(string t) {
        if (t.empty()) return;
        if (t.back() != '$') t.push_back('$');
        text = t;
        sa = build_sa(text);

        int n = (int)text.size();
        bwt.assign(n, 0);
        for (int i = 0; i < n; ++i) {
            int p = sa[i];
            bwt[i] = (p == 0) ? text[n - 1] : text[p - 1];
        }

        C.assign(256, 0);
        vector<int> freq(256, 0);
        for (char c : bwt) freq[(int)(unsigned char)c]++;
        int running = 0;
        for (int i = 0; i < 256; ++i) {
            C[i] = running;
            running += freq[i];
        }

        // NOTE: This occ layout is memory-heavy. Works for this assignment but not production-grade.
        occ.assign(n + 1, vector<int>(256, 0));
        for (int i = 0; i < n; ++i) {
            occ[i + 1] = occ[i];
            occ[i + 1][(int)(unsigned char)bwt[i]]++;
        }
    }

    pair<int,int> search(const string& p) const {
        int l = 0, r_ = (int)bwt.size();
        for (int i = (int)p.size() - 1; i >= 0; --i) {
            int c = (int)(unsigned char)p[i];
            l = C[c] + occ[l][c];
            r_ = C[c] + occ[r_][c];
            if (l >= r_) return {0,0};
        }
        return {l,r_};
    }
};

string read_fasta(const string& path) {
    ifstream in(path);
    string line, s;
    while (getline(in, line)) {
        if (!line.empty() && line[0] != '>') s += line;
    }
    return s;
}

char comp(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

string revcomp(const string& s) {
    string r;
    r.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) r.push_back(comp(*it));
    return r;
}

struct Alignment {
    bool accepted = false;
    int score = 0;
    double identity = 0.0;
    int ref_start = 0;
    int ref_end = 0;
    int aligned_len = 0;
};

// Reject very low-complexity reads to prevent spurious local alignments inflating mapping rate.
bool is_low_complexity(const string& s) {
    long long cntA = 0, cntC = 0, cntG = 0, cntT = 0, cntN = 0;
    for (char c : s) {
        switch (c) {
            case 'A': cntA++; break;
            case 'C': cntC++; break;
            case 'G': cntG++; break;
            case 'T': cntT++; break;
            default: cntN++; break;
        }
    }
    long long nonN = (long long)s.size() - cntN;
    if (nonN <= 0) return true;

    long long mx = max(max(cntA, cntC), max(cntG, cntT));
    double max_frac = (double)mx / (double)nonN;
    if (max_frac > MAX_BASE_FRAC) return true;

    // Shannon entropy in bits over A/C/G/T only (ignore N)
    auto Hterm = [&](long long c) -> double {
        if (c <= 0) return 0.0;
        double p = (double)c / (double)nonN;
        return -p * (log(p) / log(2.0));
    };
    double H = Hterm(cntA) + Hterm(cntC) + Hterm(cntG) + Hterm(cntT);
    if (H < ENTROPY_THRESH) return true;

    return false;
}

// k-mer seeding (k=20, step=10). Skip any seed containing 'N'.
// If read shorter than k, fallback to the longest non-N segment.
vector<pair<int,string>> collect_kmer_seeds(const string& s, int k, int step) {
    vector<pair<int,string>> seeds;
    int n = (int)s.size();
    if (n <= 0) return seeds;

    if (n < k) {
        int best_len = 0, best_start = -1;
        int start = 0;
        for (int i = 0; i <= n; ++i) {
            if (i == n || s[i] == 'N') {
                int len = i - start;
                if (len > best_len) { best_len = len; best_start = start; }
                start = i + 1;
            }
        }
        if (best_start >= 0 && best_len > 0) {
            seeds.push_back({best_start, s.substr(best_start, best_len)});
        }
        return seeds;
    }

    auto add_seed = [&](int off) {
        string sub = s.substr(off, k);
        if (sub.find('N') != string::npos) return;
        seeds.push_back({off, sub});
    };

    for (int i = 0; i + k <= n; i += step) add_seed(i);

    // ensure last k-mer is included
    int last = n - k;
    if (last >= 0 && (seeds.empty() || seeds.back().first != last)) add_seed(last);

    return seeds;
}

bool verify_at_start(const string& ref, const string& read, int pos, int allowed_mismatches) {
    if (pos < 0 || pos + (int)read.size() > (int)ref.size()) return false;
    int mismatches = 0;
    for (int i = 0; i < (int)read.size(); ++i) {
        if (read[i] != ref[pos + i]) {
            mismatches++;
            if (mismatches > allowed_mismatches) return false;
        }
    }
    return true;
}

Alignment smith_waterman(const string& read, const string& ref, int ws, int we) {
    const int MATCH = 2, MISMATCH = -2, GAP = -3;
    int n = (int)read.size();
    int m = we - ws;
    vector<int> prev(m + 1, 0), curr(m + 1, 0);
    vector<uint8_t> trace((n + 1) * (m + 1), 0);

    int best_score = 0, best_i = 0, best_j = 0;
    for (int i = 1; i <= n; ++i) {
        curr[0] = 0;
        for (int j = 1; j <= m; ++j) {
            char ref_c = ref[ws + j - 1];
            int diag = prev[j - 1] + ((read[i - 1] == ref_c) ? MATCH : MISMATCH);
            int up = prev[j] + GAP;
            int left = curr[j - 1] + GAP;
            int score = max(0, max(diag, max(up, left)));

            uint8_t dir = 0;
            if (score == diag && score != 0) dir = 1;
            else if (score == up && score != 0) dir = 2;
            else if (score == left && score != 0) dir = 3;

            curr[j] = score;
            trace[i * (m + 1) + j] = dir;

            if (score > best_score) { best_score = score; best_i = i; best_j = j; }
        }
        prev.swap(curr);
    }

    Alignment res; res.score = best_score;
    if (best_score == 0) return res;

    int i = best_i, j = best_j;
    int matches = 0, aligned_len = 0;
    while (i > 0 && j > 0) {
        uint8_t dir = trace[i * (m + 1) + j];
        if (dir == 0) break;
        if (dir == 1) {
            aligned_len++;
            if (read[i - 1] == ref[ws + j - 1]) matches++;
            i--; j--;
        } else if (dir == 2) {
            aligned_len++;
            i--;
        } else {
            aligned_len++;
            j--;
        }
    }

    res.accepted = true;
    res.aligned_len = aligned_len;
    res.ref_start = ws + j;
    res.ref_end = ws + best_j;
    res.identity = aligned_len ? (double)matches / (double)aligned_len : 0.0;
    return res;
}

// Map a single read:
// - low-complexity filter first
// - k-mer seeds -> candidates via vote
// - verify first, then SW if verify fails
bool map_read(
    const FMIndex& fm,
    const string& ref,
    const string& read,
    int& best_start,
    int& best_end,
    bool& unique_hit,
    bool& multi_hit,
    int& sw_calls,
    bool& sw_accepted,
    bool& verify_accepted,
    bool& had_seed_hits,
    bool& filtered_low_complexity
) {
    filtered_low_complexity = false;

    if (is_low_complexity(read)) {
        filtered_low_complexity = true;
        return false;
    }

    int ref_len = (int)ref.size();
    vector<Alignment> accepted;

    sw_calls = 0;
    sw_accepted = false;
    verify_accepted = false;
    had_seed_hits = false;

    auto try_orientation = [&](const string& oriented) {
        if (oriented.empty()) return;

        auto seeds = collect_kmer_seeds(oriented, KMER, STEP);
        if (seeds.empty()) return;

        vector<int> raw;
        raw.reserve((int)seeds.size() * MAX_PER_SEED);

        for (const auto& seed : seeds) {
            auto iv = fm.search(seed.second);
            int hits = iv.second - iv.first;
            if (hits <= 0) continue;

            had_seed_hits = true;

            if (hits <= MAX_PER_SEED) {
                for (int i = iv.first; i < iv.second; ++i) {
                    int pos = fm.sa[i] - seed.first;
                    if (pos < 0 || pos >= ref_len) continue;
                    raw.push_back(pos);
                }
            } else {
                // sample uniformly across the interval
                for (int k = 0; k < MAX_PER_SEED; ++k) {
                    double frac = (MAX_PER_SEED == 1) ? 0.0 : (double)k / (double)(MAX_PER_SEED - 1);
                    int idx = iv.first + (int)(frac * (hits - 1) + 0.5);
                    if (idx >= iv.second) idx = iv.second - 1;
                    int pos = fm.sa[idx] - seed.first;
                    if (pos < 0 || pos >= ref_len) continue;
                    raw.push_back(pos);
                }
            }
        }

        if (raw.empty()) return;

        // vote: same start position seen multiple times gets higher priority
        sort(raw.begin(), raw.end());
        vector<pair<int,int>> voted; // {pos, votes}
        voted.reserve(raw.size());
        for (int i = 0; i < (int)raw.size();) {
            int j = i + 1;
            while (j < (int)raw.size() && raw[j] == raw[i]) j++;
            voted.push_back({raw[i], j - i});
            i = j;
        }

        sort(voted.begin(), voted.end(), [](const pair<int,int>& a, const pair<int,int>& b) {
            if (a.second != b.second) return a.second > b.second;
            return a.first < b.first;
        });

        vector<int> candidates;
        candidates.reserve(min((int)voted.size(), MAX_CANDIDATES));
        for (int i = 0; i < (int)voted.size() && (int)candidates.size() < MAX_CANDIDATES; ++i) {
            candidates.push_back(voted[i].first);
        }

        // requested: cheap dedup step
        sort(candidates.begin(), candidates.end());
        candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());

        int allowed_mismatches = 0;

        for (int pos : candidates) {
            if (verify_at_start(ref, oriented, pos, allowed_mismatches)) {
                int mismatches = 0;
                for (int i = 0; i < (int)oriented.size(); ++i) {
                    if (oriented[i] != ref[pos + i]) mismatches++;
                }
                int matches = (int)oriented.size() - mismatches;

                Alignment aln;
                aln.accepted = true;
                aln.aligned_len = (int)oriented.size();
                aln.ref_start = pos;
                aln.ref_end = pos + (int)oriented.size();
                aln.identity = oriented.empty() ? 0.0 : (double)matches / (double)oriented.size();
                aln.score = matches * 2 - mismatches * 2;

                accepted.push_back(aln);
                verify_accepted = true;
                continue;
            }

            // If verify fails, try SW with stricter thresholds (more credible mapping)
            sw_calls++;
            int ws = max(0, pos - MARGIN);
            int we = min(ref_len, pos + (int)oriented.size() + MARGIN);
            if (ws >= we) continue;

            Alignment aln = smith_waterman(oriented, ref, ws, we);
            if (aln.accepted &&
                aln.aligned_len >= (int)(SW_LEN_THRESH * (double)oriented.size()) &&
                aln.identity >= SW_ID_THRESH
            ) {
                accepted.push_back(aln);
                sw_accepted = true;
            }
        }
    };

    try_orientation(read);
    try_orientation(revcomp(read));

    if (accepted.empty()) return false;

    // Dedup accepted alignments by (ref_start, ref_end): keep best score for same interval.
    sort(accepted.begin(), accepted.end(), [](const Alignment& a, const Alignment& b) {
        if (a.ref_start != b.ref_start) return a.ref_start < b.ref_start;
        if (a.ref_end != b.ref_end) return a.ref_end < b.ref_end;
        return a.score > b.score;
    });
    vector<Alignment> dedup;
    dedup.reserve(accepted.size());
    for (int i = 0; i < (int)accepted.size();) {
        int j = i + 1;
        Alignment best = accepted[i];
        while (j < (int)accepted.size() &&
               accepted[j].ref_start == accepted[i].ref_start &&
               accepted[j].ref_end == accepted[i].ref_end) {
            if (accepted[j].score > best.score) best = accepted[j];
            j++;
        }
        dedup.push_back(best);
        i = j;
    }
    accepted.swap(dedup);

    // Choose best by score
    sort(accepted.begin(), accepted.end(), [](const Alignment& a, const Alignment& b) {
        return a.score > b.score;
    });

    int best = accepted[0].score;
    int second = (accepted.size() > 1) ? accepted[1].score : numeric_limits<int>::min();
    unique_hit = (accepted.size() == 1 || best - second >= 10);
    multi_hit = !unique_hit;

    best_start = accepted.front().ref_start;
    best_end = accepted.front().ref_end;

    return true;
}

int main() {
    const int BATCH = 200000;

    cout << "[INFO] Reading reference\n" << flush;
    string ref = read_fasta("GCF_000005845.2_ASM584v2_genomic.fna");
    int ref_len = (int)ref.size();

    cout << "[INFO] Building FM-index\n" << flush;
    FMIndex fm;
    fm.build(ref);
    cout << "[INFO] FM-index built\n" << flush;

    ifstream fq("ERR022075_1.fastq");
    if (!fq) {
        cerr << "[ERROR] FASTQ not found (run from project root)\n";
        return 1;
    }

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
    cout << "[INFO] OpenMP threads: " << num_threads << "\n" << flush;
#else
    cout << "[INFO] OpenMP not enabled\n" << flush;
#endif

    vector<vector<int>> thread_cov(num_threads, vector<int>(ref_len, 0));

    long long total = 0;
    long long mapped = 0, unique = 0, multi = 0;

    // Stats/counters
    long long reads_with_seed_hits = 0;
    long long reads_that_ran_sw = 0;
    long long total_sw_calls = 0;
    long long reads_verify_accepted = 0;
    long long reads_sw_accepted = 0;
    long long reads_filtered_low_complexity = 0;

    string id, seq, plus, qual;
    vector<string> batch;
    batch.reserve(BATCH);

    cout << "[INFO] FASTQ streaming started\n" << flush;
    cout << "[INFO] Seeding: k=" << KMER << ", step=" << STEP << "\n";
    cout << "[INFO] Verify ID thresh: " << VERIFY_ID_THRESH << "\n";
    cout << "[INFO] SW thresholds: min_identity=" << SW_ID_THRESH
         << ", min_len_frac=" << SW_LEN_THRESH << "\n";
    cout << "[INFO] Low complexity: max_base_frac>" << MAX_BASE_FRAC
         << " OR entropy<" << ENTROPY_THRESH << " bits\n";

    while (true) {
        batch.clear();
        for (int i = 0; i < BATCH; ++i) {
            if (!getline(fq, id)) break;
            if (!getline(fq, seq)) break;
            if (!getline(fq, plus)) break;
            if (!getline(fq, qual)) break;
            batch.push_back(seq);
        }
        if (batch.empty()) break;

        int batch_n = (int)batch.size();
        total += batch_n;

        vector<long long> t_mapped(num_threads, 0), t_unique(num_threads, 0), t_multi(num_threads, 0);
        vector<long long> t_seed_hits(num_threads, 0), t_sw_reads(num_threads, 0), t_sw_calls(num_threads, 0);
        vector<long long> t_verify_acc(num_threads, 0), t_sw_acc(num_threads, 0);
        vector<long long> t_lowcomp(num_threads, 0);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 256)
#endif
        for (int i = 0; i < batch_n; ++i) {
            int tid =
#ifdef _OPENMP
                omp_get_thread_num();
#else
                0;
#endif

            int best_start = 0, best_end = 0;
            bool uniq = false, multi_flag = false;
            int sw_calls = 0;
            bool sw_acc = false, ver_acc = false, seed_hits = false, lowcomp = false;

            bool mapped_flag = map_read(
                fm, ref, batch[i],
                best_start, best_end,
                uniq, multi_flag,
                sw_calls, sw_acc, ver_acc, seed_hits, lowcomp
            );

            if (lowcomp) {
                t_lowcomp[tid]++;
            } else {
                if (seed_hits) t_seed_hits[tid]++;
                if (sw_calls > 0) t_sw_reads[tid]++;
                t_sw_calls[tid] += sw_calls;
                if (ver_acc) t_verify_acc[tid]++;
                if (sw_acc) t_sw_acc[tid]++;
            }

            if (mapped_flag) {
                t_mapped[tid]++;
                if (uniq) t_unique[tid]++; else if (multi_flag) t_multi[tid]++;

                auto& cov = thread_cov[tid];
                int start = max(0, best_start);
                int end = min(ref_len, best_end);
                for (int p = start; p < end; ++p) cov[p]++;
            }

            long long done = (total - batch_n) + i + 1;
            if (done % 500000 == 0) {
                cout << "[INFO LOG] Processed: " << done << " reads\n" << flush;
            }
        }

        for (int t = 0; t < num_threads; ++t) {
            mapped += t_mapped[t];
            unique += t_unique[t];
            multi += t_multi[t];

            reads_with_seed_hits += t_seed_hits[t];
            reads_that_ran_sw += t_sw_reads[t];
            total_sw_calls += t_sw_calls[t];
            reads_verify_accepted += t_verify_acc[t];
            reads_sw_accepted += t_sw_acc[t];
            reads_filtered_low_complexity += t_lowcomp[t];
        }
    }

    // Reduce coverage arrays
    vector<int> coverage(ref_len, 0);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < ref_len; ++i) {
        int sum = 0;
        for (const auto& cov : thread_cov) sum += cov[i];
        coverage[i] = sum;
    }

    long long cov_sum = 0;
    int cov_ge1 = 0, cov_ge10 = 0;
    for (int v : coverage) {
        cov_sum += v;
        if (v >= 1) cov_ge1++;
        if (v >= 10) cov_ge10++;
    }

    double mean_cov = ref_len ? (double)cov_sum / (double)ref_len : 0.0;
    double cov1_pct = ref_len ? 100.0 * (double)cov_ge1 / (double)ref_len : 0.0;
    double cov10_pct = ref_len ? 100.0 * (double)cov_ge10 / (double)ref_len : 0.0;

    long long unmapped = total - mapped;
    double mapping_rate = total ? 100.0 * (double)mapped / (double)total : 0.0;

    // Percent of total (as before)
    double unique_pct_total = total ? 100.0 * (double)unique / (double)total : 0.0;
    double multi_pct_total = total ? 100.0 * (double)multi / (double)total : 0.0;

    // Percent of mapped (more standard)
    double unique_pct_mapped = mapped ? 100.0 * (double)unique / (double)mapped : 0.0;
    double multi_pct_mapped = mapped ? 100.0 * (double)multi / (double)mapped : 0.0;

    double seed_hit_pct = total ? 100.0 * (double)reads_with_seed_hits / (double)total : 0.0;
    double sw_read_pct = total ? 100.0 * (double)reads_that_ran_sw / (double)total : 0.0;

    cout << "[INFO] Mapping finished\n" << flush;
    cout << "Total reads: " << total << "\n";
    cout << "Mapped reads: " << mapped << " (" << mapping_rate << "%)\n";
    cout << "Unique reads: " << unique
         << " (" << unique_pct_total << "% of total, " << unique_pct_mapped << "% of mapped)\n";
    cout << "Multi-mapped reads: " << multi
         << " (" << multi_pct_total << "% of total, " << multi_pct_mapped << "% of mapped)\n";
    cout << "Unmapped reads: " << unmapped << "\n";

    cout << "Low-complexity filtered reads: " << reads_filtered_low_complexity << "\n";
    cout << "Reads with seed hits (not low-comp): " << reads_with_seed_hits << " (" << seed_hit_pct << "%)\n";
    cout << "Reads that ran SW (not low-comp): " << reads_that_ran_sw << " (" << sw_read_pct << "%)\n";
    cout << "Total SW calls: " << total_sw_calls << "\n";
    cout << "Reads with at least one verify-accepted alignment: " << reads_verify_accepted << "\n";
    cout << "Reads with at least one SW-accepted alignment: " << reads_sw_accepted << "\n";

    cout << "Mean coverage: " << mean_cov << "\n";
    cout << "Genome covered >=1x: " << cov1_pct << "%\n";
    cout << "Genome covered >=10x: " << cov10_pct << "%\n";

    return 0;
}

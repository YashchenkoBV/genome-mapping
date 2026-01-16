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

struct Alignment { bool accepted=false; int score=0; double identity=0.0; int ref_start=0; int ref_end=0; int aligned_len=0; };

vector<pair<int,string>> collect_seeds(const string& s, int min_len) {
    vector<pair<int,string>> seeds;
    int start = 0;
    for (int i = 0; i <= (int)s.size(); ++i) {
        if (i == (int)s.size() || s[i] == 'N') {
            int len = i - start;
            if (len >= min_len) seeds.push_back({start, s.substr(start, len)});
            start = i + 1;
        }
    }
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
        if (dir == 1) { aligned_len++; if (read[i - 1] == ref[ws + j - 1]) matches++; i--; j--; }
        else if (dir == 2) { aligned_len++; i--; }
        else { aligned_len++; j--; }
    }
    res.accepted = true;
    res.aligned_len = aligned_len;
    res.ref_start = ws + j;
    res.ref_end = ws + best_j;
    res.identity = aligned_len ? (double)matches / aligned_len : 0.0;
    return res;
}

bool map_read(const FMIndex& fm, const string& ref, const string& read, int min_seed_len, int& best_start, int& best_end, bool& unique_hit, bool& multi_hit) {
    const int MAX_CANDIDATES = 200;
    const int MARGIN = 30;
    const double ID_THRESH = 0.90;
    const double LEN_THRESH = 0.7;
    int ref_len = (int)ref.size();
    vector<Alignment> accepted;

    auto try_orientation = [&](const string& oriented) {
        auto seeds = collect_seeds(oriented, min_seed_len);
        if (seeds.empty()) return;
        sort(seeds.begin(), seeds.end(), [](const pair<int,string>& a, const pair<int,string>& b) { return a.second.size() > b.second.size(); });
        int allowed_mismatches = (int)((1.0 - ID_THRESH) * oriented.size());
        for (const auto& seed : seeds) {
            auto iv = fm.search(seed.second);
            int hits = iv.second - iv.first;
            if (hits == 0) continue;

            vector<int> candidates;
            if (hits <= MAX_CANDIDATES) {
                for (int i = iv.first; i < iv.second && (int)candidates.size() < MAX_CANDIDATES; ++i) {
                    int pos = fm.sa[i] - seed.first;
                    if (pos < 0 || pos >= ref_len) continue;
                    candidates.push_back(pos);
                }
            } else {
                for (int k = 0; k < MAX_CANDIDATES; ++k) {
                    double frac = (MAX_CANDIDATES == 1) ? 0.0 : (double)k / (double)(MAX_CANDIDATES - 1);
                    int idx = iv.first + (int)(frac * (hits - 1) + 0.5);
                    if (idx >= iv.second) idx = iv.second - 1;
                    int pos = fm.sa[idx] - seed.first;
                    if (pos < 0 || pos >= ref_len) continue;
                    candidates.push_back(pos);
                }
            }
            if (candidates.empty()) continue;

            for (int pos : candidates) {
                if (verify_at_start(ref, oriented, pos, allowed_mismatches)) {
                    int mismatches = 0;
                    for (int i = 0; i < (int)oriented.size(); ++i) if (oriented[i] != ref[pos + i]) mismatches++;
                    int matches = (int)oriented.size() - mismatches;
                    Alignment aln;
                    aln.accepted = true;
                    aln.aligned_len = (int)oriented.size();
                    aln.ref_start = pos;
                    aln.ref_end = pos + (int)oriented.size();
                    aln.identity = oriented.empty() ? 0.0 : (double)matches / (double)oriented.size();
                    aln.score = matches * 2 - mismatches * 2;
                    accepted.push_back(aln);
                    continue;
                }
                int ws = max(0, pos - MARGIN);
                int we = min(ref_len, pos + (int)oriented.size() + MARGIN);
                if (ws >= we) continue;
                Alignment aln = smith_waterman(oriented, ref, ws, we);
                if (aln.aligned_len >= (int)(LEN_THRESH * oriented.size()) && aln.identity >= ID_THRESH) {
                    accepted.push_back(aln);
                }
            }
            if (!accepted.empty()) break;
        }
    };

    try_orientation(read);
    try_orientation(revcomp(read));

    if (accepted.empty()) return false;
    sort(accepted.begin(), accepted.end(), [](const Alignment& a, const Alignment& b) { return a.score > b.score; });
    int best = accepted[0].score;
    int second = (accepted.size() > 1) ? accepted[1].score : numeric_limits<int>::min();
    unique_hit = (accepted.size() == 1 || best - second >= 10);
    multi_hit = !unique_hit;
    best_start = accepted.front().ref_start;
    best_end = accepted.front().ref_end;
    return true;
}

int main() {
    const int MIN_SEED = 20;
    const int BATCH = 200000;

    cout << "[INFO] Reading reference\n" << flush;
    string ref = read_fasta("GCF_000005845.2_ASM584v2_genomic.fna");
    int ref_len = (int)ref.size();

    cout << "[INFO] Building FM-index\n" << flush;
    FMIndex fm;
    fm.build(ref);
    cout << "[INFO] FM-index built\n" << flush;

    ifstream fq("ERR022075_1.fastq");
    if (!fq) { cerr << "[ERROR] FASTQ not found (run from project root)\n"; return 1; }

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = omp_get_max_threads();
    cout << "[INFO] OpenMP threads: " << num_threads << "\n" << flush;
#else
    cout << "[INFO] OpenMP not enabled\n" << flush;
#endif

    vector<vector<int>> thread_cov(num_threads, vector<int>(ref_len, 0));
    int total = 0, mapped = 0, unique = 0, multi = 0;

    string id, seq, plus, qual;
    vector<string> batch;
    batch.reserve(BATCH);
    cout << "[INFO] FASTQ streaming started\n" << flush;

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

        vector<int> t_mapped(num_threads, 0), t_unique(num_threads, 0), t_multi(num_threads, 0);

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
            bool mapped_flag = map_read(fm, ref, batch[i], MIN_SEED, best_start, best_end, uniq, multi_flag);
            if (mapped_flag) {
                t_mapped[tid]++;
                if (uniq) t_unique[tid]++; else if (multi_flag) t_multi[tid]++;
                auto& cov = thread_cov[tid];
                int start = max(0, best_start);
                int end = min(ref_len, best_end);
                for (int p = start; p < end; ++p) cov[p]++;
            }

            int done = (total - batch_n) + i + 1;
            if (done % 500000 == 0) {
                double pct = total ? 100.0 * (double)done / (double)total : 0.0;
                cout << "[INFO LOG] Processed: " << done << " (" << pct << " %)\n" << flush;
            }
        }

        for (int t = 0; t < num_threads; ++t) { mapped += t_mapped[t]; unique += t_unique[t]; multi += t_multi[t]; }
    }

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
    double mapping_rate = total ? 100.0 * (double)mapped / (double)total : 0.0;
    double unique_pct = total ? 100.0 * (double)unique / (double)total : 0.0;
    double multi_pct = total ? 100.0 * (double)multi / (double)total : 0.0;
    int unmapped = total - mapped;

    cout << "[INFO] Mapping finished\n" << flush;
    cout << "Total reads: " << total << "\n";
    cout << "Mapped reads: " << mapped << " (" << mapping_rate << "%)\n";
    cout << "Unique reads: " << unique << " (" << unique_pct << "%)\n";
    cout << "Multi-mapped reads: " << multi << " (" << multi_pct << "%)\n";
    cout << "Mean coverage: " << mean_cov << "\n";
    cout << "Genome covered >=1x: " << cov1_pct << "%\n";
    cout << "Genome covered >=10x: " << cov10_pct << "%\n";
    cout << "Unmapped reads: " << unmapped << "\n";
}

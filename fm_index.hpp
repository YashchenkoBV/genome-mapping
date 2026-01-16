#pragma once
#include <vector>
#include <string>
#include <array>
#include <iostream>
#include <algorithm>

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

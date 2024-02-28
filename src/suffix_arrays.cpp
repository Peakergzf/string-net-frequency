#include "suffix_arrays.hpp"

#include <divsufsort.h>

#include <string_view>
#include <iostream>
#include <cassert>
#include <ctime>
#include <stack>
#include <functional>
#include <numeric>
#include <unordered_set>


AugmentedSuffixArray::AugmentedSuffixArray(std::string_view _txt) :
        txt(_txt),
        n((uint32_t) txt.size()),
        sa(std::make_unique<uint32_t[]>(n)),
        isa(n),
        lf(n),
        lcp(n),
        C(n, -1),
        RMQ(n) {
    clock_t t0 = clock();

    assert(divsufsort((sauchar_t *) txt.data(), (saidx_t *) sa.get(), (saidx_t) n) == 0);

    std::cout << "(SA constructed in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n";

    t0 = clock();

    for (uint32_t i = 0; i < n; i++) isa[sa[i]] = i;

    for (size_t i = 0; i < n; i++) if (sa[i] > 0) lf[i] = isa[sa[i] - 1];

    std::cout << "(ISA and LF constructed in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n";

    t0 = clock();

    construct_lcp();

    std::cout << "(LCP constructed in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n";

    t0 = clock();

    construct_c_array();
    RMQ.v = C;
    RMQ.construct();

    std::cout << "(C array and RMQ constructed in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n";
}


// kasai's algorithm
void AugmentedSuffixArray::construct_lcp() {
    // idea: iterate through each suffix T_i in text order,
    // rather than in sorted order (i.e., the order the lcp values are stored)
    // (note that T_i is at index isa[i] in sorted order)
    // at each iteration, we compute the current lcp value l,
    // which is stored at index k := isa[i]+1,
    // also, l is the lcp between T_i and T_j where j := ASA[k],
    // so overall, we want to compute lcp[k] = lcp(T_i, T_j),
    // (as a remainder, i and j are indices of the text,
    //  whereas k is an index of the lcp array)
    // the computation might involve some explicit character comparison,
    // but in text order, from T_i to T_{i+1}, we can reuse the matching results,
    // since T_{i+1} is just T_i without the first character
    uint32_t l = 0;
    for (size_t i = 0; i < n; i++) {
        uint32_t k = isa[i] + 1;
        if (k == n) { // T_{ASA[n]} doesn't exist
            assert(l == 0);
            continue;
        }
        uint32_t j = sa[k];
        // overall there are at most 3n explicit character comparisons:
        //     at most n comparisons ended with mismatch (one per each iteration);
        //     exactly l comparisons ended with match (one per each increment of l),
        //     (l <= 2n because it is always bounded by n,
        //      and it is decremented at most n times)
        while (i + l < n && j + l < n && txt[i + l] == txt[j + l]) l++;
        lcp[k] = l;
        // removing the first character
        if (l > 0) l--;
    }
}


void AugmentedSuffixArray::construct_c_array() {
    std::unordered_map<char, int> last_occ;
    for (int j = 0; j < n; j++) {
        if (sa[j] == 0) continue;
        auto c = txt[sa[j] - 1]; // BWT[j]
        if (last_occ.find(c) != last_occ.end()) {
            C[j] = last_occ.at(c);
        }
        last_occ[c] = j;
    }
}

void AugmentedSuffixArray::_crl_aux(size_t i, size_t j, size_t x, std::vector<size_t> &idx) {
    if (i > j) return;
    auto q = RMQ.query(i, j);
    if (C[q] >= static_cast<int>(x)) return;
    // casting is needed because C[q] could be -1
    idx.push_back(q);
    _crl_aux(i, q - 1, x, idx);
    _crl_aux(q + 1, j, x, idx);
}

std::vector<size_t> AugmentedSuffixArray::crl(size_t i, size_t j) {
    std::vector<size_t> idx;
    _crl_aux(i, j, i, idx);
    return idx;
}

uint32_t AugmentedSuffixArray::single_nf_crl(std::string_view s) {
    auto lb = find_interval_lb(s);
    if (!lb.has_value()) return 0;
    auto l = lb.value();
    auto r = find_interval_rb(s).value();
    if (l == r) return 0;

    uint32_t nf = 0, len = (uint32_t) s.size();
    for (auto i: crl(l, r)) {
        if (lcp[i] <= len && (i == n - 1 || lcp[i + 1] <= len)) {
            auto j = lf[i];
            if (lcp[j] <= len && (j == n - 1 || lcp[j + 1] <= len)) {
                nf++;
            }
        }
    }
    return nf;
}


std::optional<uint32_t> AugmentedSuffixArray::find_interval_lb(std::string_view s) {
    size_t len = s.size();
    int low = -1, high = (int) n - 1;
    while (low < high - 1) {
        int mid = (low + high) / 2;
        auto s_mid = txt.substr(sa[(size_t) mid], len);
        if (s > s_mid) low = mid;
        else high = mid;
    }
    if (txt.substr(sa[(size_t) high], len) == s) return high;
    else return {};
}


std::optional<uint32_t> AugmentedSuffixArray::find_interval_rb(std::string_view s) {
    size_t len = s.size();
    size_t low = 0, high = n;
    while (low < high - 1) {
        size_t mid = (low + high) / 2;
        auto s_mid = txt.substr(sa[mid], len);
        if (s >= s_mid) low = mid;
        else high = mid;
    }
    if (txt.substr(sa[low], len) == s) return low;
    else return {};
}


// define a hash function for std::pair<char, char>, used in std::unordered_map
namespace std {
    template<>
    struct hash<std::pair<char, char>> {
        size_t operator()(const std::pair<char, char> &key) const noexcept {
            return key.first << sizeof(char) | key.second;
        }
    };
}


uint32_t AugmentedSuffixArray::single_nf_hash(std::string_view s) {
    auto lb = find_interval_lb(s);
    if (!lb.has_value()) return 0;
    auto l = lb.value();
    auto r = find_interval_rb(s).value();
    if (l == r) return 0; // s is unique

    std::unordered_map<char, uint32_t> R, L;
    std::unordered_map<std::pair<char, char>, uint32_t> B;
    for (auto i = l; i <= r; i++) {
        if (sa[i] - 1 >= 0 && sa[i] + s.size() <= n - 1) {
            auto y = txt[sa[i] + s.size()];
            auto x = txt[sa[i] - 1];
            R[y]++;
            L[x]++;
            B[{x, y}]++;
        }
    }
    uint32_t nf = 0;
    for (auto [xy, f]: B) {
        auto [x, y] = xy;
        if (L.at(x) == 1 and R.at(y) == 1) {
            assert(f == 1);
            nf++;
        }
    }
    return nf;
}


uint32_t AugmentedSuffixArray::single_nf(std::string_view s) {
    auto lb = find_interval_lb(s);
    if (!lb.has_value()) return 0;
    auto l = lb.value();

    uint32_t nf = 0, len = (uint32_t) s.size(), i = l;
    while (true) {
        if (lcp[i] <= len && (i == n - 1 || lcp[i + 1] <= len)) {
            auto j = lf[i];
            if (lcp[j] <= len && (j == n - 1 || lcp[j + 1] <= len)) {
                nf++;
            }
        }
        if (i == n - 1 || lcp[i + 1] < len) break;
        i++;
    }
    if (i == l) return 0;
    return nf;
}


void AugmentedSuffixArray::all_nf(bool report) {
    clock_t t0 = clock();
    for (size_t i = 0; i < n; i++) {
        auto len = (i == n - 1) ? lcp[i] : std::max(lcp[i + 1], lcp[i]);
        if (len == 0) continue;
        auto j = lf[i];
        if ((j == n - 1 || lcp[j + 1] <= len) && lcp[j] <= len) {
            pos_nf[txt.substr(sa[i], len)]++;
        }
    }
    if (report) {
        std::ofstream out_file("/dev/null");
        for (const auto &[s, nf]: pos_nf) {
            out_file << s.size() << " " << nf << '\n';
        }
    }
    std::cout << "(all-nf " << (report ? "reported" : "computed") << " using SA in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n";
}


void AugmentedSuffixArray::LCPInterval::process(bool report) {
    if (nf == 0) return;
    if (report) {
        std::ofstream out_file("/dev/null");
        out_file << suf_arr.sa[lb] << " " << len << " " << nf << '\n';
    } else {
        auto s = suf_arr.txt.substr(suf_arr.sa[lb], len);
        suf_arr.pos_nf_stack[s] = nf;
    }
}

void AugmentedSuffixArray::LCPInterval::print() {
    auto s = suf_arr.txt.substr(suf_arr.sa[lb], len);
    std::cout << len << "-[" << lb << ", " << rb << "], φ(" << s << ")=" << nf << '\n';
}


// bottom-up traversal of the internal nodes
void AugmentedSuffixArray::all_nf_stack(bool report) {
    clock_t t0 = clock();

    std::stack<LCPInterval> stk;
    // push the root node onto the stack
    // (whenever we push a node, rb is always initialised to 0
    //  because we don't know where the interval ends yet)
    stk.push({0, 0, 0, 0, *this});

    bool for_next = false;

    for (uint32_t i = 1; i < n; i++) {
        std::cout << "i = " << i + 1 << ":\n";
        uint32_t lb = i - 1; // used for push later
        while (lcp[i] < stk.top().len) { // mismatch
            stk.top().rb = i - 1; // we now know where the current streak ends
            auto itv = stk.top();
            stk.pop();
            itv.process(report);
            lb = itv.lb; // update with the parent lb
        }
        if (lcp[i] > stk.top().len) { // can extend the match
            stk.push({lcp[i], lb, 0, 0, *this});
            std::cout << "\tpush: ";
            std::cout << lcp[i] << "-[" << lb << ", -]\n";
            if (for_next) {
                stk.top().nf++;
                for_next = false;
            }
        }
        // check for net occurrence
        auto ell = (i == n - 1) ? lcp[i] : std::max(lcp[i + 1], lcp[i]);
        auto j = lf[i];
        if ((j == n - 1 || lcp[j + 1] <= ell) && lcp[j] <= ell) {
            if (lcp[i] == ell) stk.top().nf++;
            else for_next = true;
        }
        // an invariant at the end of each loop
        assert(lcp[i] == stk.top().len);
    }
    // while there are internal nodes (other than the root) left
    while (stk.size() > 1) {
        stk.top().rb = n - 1;
        std::cout << "out of loop:\n";
        stk.top().process(report);
        stk.pop();
    }
    std::cout << "(all-nf " << (report ? "reported" : "computed") << " using SA (with a stack) in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n";
}


uint32_t AugmentedSuffixArray::freq(std::string_view s) {
    auto lb = find_interval_lb(s);
    if (!lb.has_value()) return 0;
    auto l = lb.value();
    auto r = find_interval_rb(s).value();
    return r - l + 1;
}


std::pair<uint32_t, double> AugmentedSuffixArray::lcp_stats() {
    auto lcp_sum = std::accumulate(lcp.begin(), lcp.end(), (unsigned long long) 0);

    auto _lcp = lcp; // make a copy for partial sorting later
    auto mid_idx = _lcp.size() / 2;
    // the middle element is in sorted position
    std::nth_element(_lcp.begin(), _lcp.begin() + (long) mid_idx, _lcp.end());

    return {_lcp[mid_idx], static_cast<float>(lcp_sum) / static_cast<float>(n - 1)};
}


void AugmentedSuffixArray::get_text_stats(size_t len_ub, const std::string &out_file_name) {
    all_nf(false);
    std::unordered_set<char> alphabet;
    for (auto c: txt) alphabet.insert(c);

    unsigned long long tot_len = 0, tot_nf = 0, tot_nf_len = 0;
    unsigned long long tot_len_ub = 0, tot_nf_len_ub = 0;
    for (auto [s, nf]: pos_nf) {
        tot_len += s.size();
        tot_nf += nf;
        tot_nf_len += nf * s.size();
        if (s.size() <= len_ub) {
            tot_len_ub += s.size();
            tot_nf_len_ub += nf * s.size();
        }
    }
    assert(tot_nf <= n);

    auto [lcp_med, lcp_avg] = lcp_stats();

    auto _n = static_cast<float>(n);

    std::cout << "\ntext stats:\n"
              << " -> |T| = " << txt.size() << ", |Σ| = " << alphabet.size() << "\n"
              << " -> median LCP: " << lcp_med << ", average LCP: " << lcp_avg << "\n"
              << " -> " << pos_nf.size() << " strings with positive NF:\n"
              << "\t -> total string length = " << static_cast<float>(tot_len) / _n << " * n\n"
              << "\t -> total occurrence length = " << static_cast<float>(tot_nf_len) / _n << " * n\n"
              << "\t -> total string length (with u.b.) = " << static_cast<float>(tot_len_ub) / _n << " * n\n"
              << "\t -> total occurrence length (with u.b.)= " << static_cast<float>(tot_nf_len_ub) / _n << " * n\n\n";

    if (!out_file_name.empty()) { // non-empty file name
        std::ofstream out_file(out_file_name);

        // (use map instead of unordered_map because pair doesn't have a default hash function;
        //  use map also to see the result in sorted order)
        std::map<std::pair<size_t, size_t>, size_t> nf_len_dist;

        for (auto [s, nf]: pos_nf) {
            nf_len_dist[std::make_pair(s.size(), nf)]++;
        }

        for (auto [pr, cnt]: nf_len_dist) {
            auto [len, nf] = pr;
            out_file << len << '\t' << nf << '\t' << cnt << '\n';
        }
    }
}


CompressedSuffixArray::CompressedSuffixArray(std::string_view _txt) : n(_txt.size() + 1) {
    std::string txt{_txt};

    clock_t t0 = clock();
    sdsl::construct_im(sa, txt, 1);
    std::cout << "(compressed SA and LF constructed in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n";

    t0 = clock();
    sdsl::construct_im(lcp, txt, 1);
    std::cout << "(compressed LCP constructed in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n";

    assert(sa.size() == n && lcp.size() == n);
    assert(sa[0] == n - 1);

//    std::unordered_map<char, int> last_occ;
//    for (int j = 0; j < n; j++) {
//        auto k = sa[j];
//        if (k == 0) continue;
//        auto c = txt[k - 1]; // BWT[j]
//        if (last_occ.find(c) != last_occ.end()) {
//            C[j] = last_occ.at(c);
//        }
//        last_occ[c] = j;
//    }
}


size_t CompressedSuffixArray::single_nf(std::string_view s) {
    auto [l, r] = sdsl::lex_interval(sa, s.begin(), s.end());
    if (l >= r) return 0;
    size_t nf = 0, len = s.size();
    for (auto i = l; i <= r; i++) {
        if ((i == n - 1 || lcp[i + 1] <= len) && lcp[i] <= len) {
            auto j = sa.lf[i];
            if ((j == n - 1 || lcp[j + 1] <= len) && lcp[j] <= len) {
                nf++;
            }
        }
    }
    return nf;
}

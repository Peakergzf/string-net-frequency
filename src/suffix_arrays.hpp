#pragma once

#include <string_view>
#include <vector>
#include <optional>
#include <unordered_map>
#include <memory>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp.hpp>


// the following implementation of range minimum query is based on
// https://codeforces.com/blog/entry/78931
// author: Bruno Monteiro
template<typename T>
struct rmq {
    std::vector<T> v;
    size_t n;
    static const size_t b = 30;
    std::vector<size_t> mask, t;

    size_t op(size_t x, size_t y) { return v[x] < v[y] ? x : y; }

    size_t msb(size_t x) { return __builtin_clz(1) - __builtin_clz(x); }

    size_t small(size_t r, size_t sz = b) { return r - msb(mask[r] & ((1 << sz) - 1)); }

    explicit rmq(size_t _n) : v(_n), n(_n), mask(n), t(n) {}

    void construct() {
        for (size_t i = 0, at = 0; i < n; mask[i++] = at |= 1) {
            at = (at << 1) & ((1 << b) - 1);
            while (at and op(i, i - msb(at & -at)) == i) at ^= at & -at;
        }
        for (size_t i = 0; i < n / b; i++) t[i] = small(b * i + b - 1);
        for (size_t j = 1; (1 << j) <= n / b; j++)
            for (size_t i = 0; i + (1 << j) <= n / b; i++)
                t[n / b * j + i] = op(t[n / b * (j - 1) + i], t[n / b * (j - 1) + i + (1 << (j - 1))]);
    }

    size_t query(size_t l, size_t r) {
        if (r - l + 1 <= b) return small(r, r - l + 1);
        size_t ans = op(small(l + b - 1), small(r));
        size_t x = l / b + 1, y = r / b - 1;
        if (x <= y) {
            size_t j = msb(y - x + 1);
            ans = op(ans, op(t[n / b * j + x], t[n / b * j + y - (1 << j) + 1]));
        }
        return ans;
    }
};


class AugmentedSuffixArray {
private:
    std::string_view txt;
    uint32_t n;
    // suffix array (a pointer is used for divsufsort (see constructor))
    std::unique_ptr<uint32_t[]> sa;
    // inverse suffix array, LF mapping, and LCP array
    std::vector<uint32_t> isa, lf, lcp;

    // represents an LCP interval `len-[lb, rb]` with net frequency `nf`,
    // (it is only used in `all_nf_stack`)
    class LCPInterval {
    public:
        // each LCP interval corresponds to a string `s` with `len = |s|` and `nf(s) = nf`,
        // moreover, `[lb, rb]` is an interval in SA in which all the suffixes are prefixed by `s`
        // (`s` is also a branching string (i.e., corresponds to an internal node in a suffix tree))
        uint32_t len, lb, rb, nf;

        // a reference to the outer class
        // (this is necessary for the inner class to access the members of the outer class)
        AugmentedSuffixArray &suf_arr;

        // report positive NF strings if report is true, preprocess them otherwise
        void process(bool report);

        void print();
    };

    void construct_lcp();

public:
    // strings with positive net frequency
    // (the first one is computed using `all_nf`, the second using `all_nf_stack`)
    std::unordered_map<std::string_view, uint32_t> pos_nf, pos_nf_stack;

    // constructor
    explicit AugmentedSuffixArray(std::string_view _txt);

    // locate a string and return the suffix array interval
    // (implemented by a binary search)
    std::optional<uint32_t> find_interval_lb(std::string_view s);

    std::optional<uint32_t> find_interval_rb(std::string_view s);

    // single-NF algorithm using a hash table
    uint32_t single_nf_hash(std::string_view s);

    // single-NF algorithm using LF and LCP
    uint32_t single_nf(std::string_view s);

private:
    // C[j] := max { i | i < j and BWT[i] = BWT[j] } if such i exists, -1 otherwise
    std::vector<int> C;
    // the RMQ support on C
    rmq<int> RMQ;

    void construct_c_array();

    // coloured range listing
    // crl(i, j) := the positions of the distinct characters in BWT[i ... j]
    std::vector<size_t> crl(size_t i, size_t j);

    // an auxiliary function used when computing crl(i, j)
    void _crl_aux(size_t i, size_t j, size_t x, std::vector<size_t> &idx);

public:
    // single-NF algorithm using coloured range listing
    uint32_t single_nf_crl(std::string_view s);

    // all-NF using a hash table
    void all_nf(bool report);

    // all-NF using a stack
    void all_nf_stack(bool report);

    // number of occurrences
    uint32_t freq(std::string_view s);

    // compute the median and average of the lcp array
    std::pair<uint32_t, double> lcp_stats();

    // print basic text stats and write net frequency distribution against string length result to a file
    // if `out_file_name` is non-empty (empty by default)
    void get_text_stats(size_t len_ub, const std::string &out_file_name = "");
};


class CompressedSuffixArray {
private:
    sdsl::csa_wt<> sa;
    sdsl::lcp_wt<> lcp;
    size_t n; // the size of SA and LCP (also the size of the text + 1)

public:
    explicit CompressedSuffixArray(std::string_view txt);

    size_t single_nf(std::string_view s);
};

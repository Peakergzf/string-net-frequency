#pragma once

#include <string_view>
#include <vector>
#include <optional>
#include <unordered_map>
#include <memory>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp.hpp>


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

    // all-NF using a hash table
    void all_nf(bool report);
    // all-NF using a stack
    void all_nf_stack(bool report);

    // true frequency
    uint32_t freq(std::string_view s);

    // compute the median and average of the lcp array
    std::pair<uint32_t, double> lcp_stats();
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

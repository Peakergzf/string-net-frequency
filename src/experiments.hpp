#pragma once

#include <utility>
#include "suffix_arrays.hpp"


namespace experiment {

    // ===========================================================================================
    //                                enum and constants
    // ===========================================================================================

    // the net frequency algorithms used in the experiments
    enum class Algo {
        ASA, HSA, CSA, CRL, size
        // the correct size of the enum class relies on the fact that
        // none of the enumerators has an initializer
        // (i.e., they are between 0 and size-1 implicitly)
    };
    // the number of algorithms
    const size_t A = (size_t) Algo::size;

    const std::unordered_map<Algo, std::string> algo_name = {
            {Algo::ASA, "ASA"},
            {Algo::HSA, "HSA"},
            {Algo::CSA, "CSA"},
            {Algo::CRL, "CRL"},
    };

    // query string length lower and upper bounds
    const uint32_t LEN_LB = 5, LEN_UB = 35;

    // ===========================================================================================
    //                              query and result related
    // ===========================================================================================

    // a classed for generating different types of queries from an input text
    class Query {
    private:
        // the input text
        std::string_view txt;

    public:
        explicit Query(std::string_view _txt) : txt(_txt) {}

        // p is the probability of including a query string
        void generate_tokenized_queries(std::string_view file_name, double p);

        void generate_dna_queries(std::string_view file_name, double p);
    };

    // the aggregated result of a single algorithm on multiple query strings,
    // specifically, each run of an algorithm results in a quadruple:
    // {query string length, frequency, net frequency, query time},
    // we record the total query time for each of the three quantities separately
    class Result {
    public:
        // {q}_time[x] := sum of query time for strings of {q} x,
        // where q = {len, freq, nf}
        std::unordered_map<size_t, uint64_t> len_time, freq_time, nf_time;

        void update(size_t len, size_t freq, size_t nf, size_t time);
    };

    // ===========================================================================================
    //                              single-NF experiment
    // ===========================================================================================

    class SingleNF {
    private:
        // the input text
        std::string_view txt;

        // the data structures used
        std::unique_ptr<AugmentedSuffixArray> asa;
        std::unique_ptr<CompressedSuffixArray> csa;

        // the result for each algorithm
        // (query time could be measured in either cpu time or wall time)
        std::vector<Result> results_cpu, results_wall;

        // the following are independent of the algorithms
        // {q}_cnt[x] := number of substrings of {q} x
        // ({q}_time and {q}_cnt together can be used to compute the average query time)
        std::unordered_map<size_t, uint64_t> len_cnt, freq_cnt, nf_cnt;

        // the name of the directory that stores the results
        std::string result_directory;

    public:
        // constructor
        SingleNF(std::string_view _txt, std::string dir_name) :
                txt(_txt),
                asa(std::make_unique<AugmentedSuffixArray>(txt)),
                csa(std::make_unique<CompressedSuffixArray>(txt)),
                results_cpu(A),
                results_wall(A),
                result_directory(std::move(dir_name)) {}

        // run experiment using the queries from the file
        void run_queries(std::string_view query_file_name);

        // run each algorithm on query string s and record the time taken
        void time_algo(std::string_view s);

        // print the total time taken by each algorithm
        void print_total_time();

        void write_results();

        void print_text_stats();
    };

    // ===========================================================================================
    //                             helper functions
    // ===========================================================================================

    // return true if the file is empty, false otherwise
    bool is_empty(std::string_view file_name);

    void process_text(std::string_view in_file_name, std::string_view out_file_name);

    void run();
}

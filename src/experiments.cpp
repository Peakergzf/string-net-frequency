#include "./experiments.hpp"

#include <iostream>
#include <set>
#include <cassert>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <regex>
#include <random>
#include <filesystem>


void experiment::Query::generate_tokenized_queries(std::string_view file_name, double p) {

    assert(is_empty(file_name));
    std::string _file_name {file_name};
    std::ofstream file(_file_name, std::ios_base::app);

    clock_t t0 = clock();
    // convert from string_view to string
    std::string _txt{txt};
    // tokenize the input text based on spaces
    std::istringstream iss(_txt);
    std::vector<std::string> tokens;
    std::copy(
            std::istream_iterator<std::string>(iss),
            std::istream_iterator<std::string>(),
            std::back_inserter(tokens)
    );

    std::random_device rd;
    std::default_random_engine double_gen(rd());
    // uniform distribution on the half-closed interval [0.0, 1.0)
    std::uniform_real_distribution<double> uni_0_1(0.0, 1.0);

    // a query is the concatenation of all the tokens between tokens[i] and tokens[j]
    size_t i = 0, j = i; // start and end index

    // total number queries
    size_t cnt = 0;

    std::string query;
    while (j < tokens.size()) {
        // query becomes too long to add the next token
        if (tokens[i].size() > LEN_UB || query.size() + tokens[j].size() + 1 > LEN_UB) {
            i++;
            j = i;
            query = "";
            continue;
        }
        // no leading space initially (i.e., when j==i)
        if (j > i) query += ' ';
        query += tokens[j];
        j++;
        if (query.size() >= LEN_LB) { // query is long enough
            if (query.front() == txt.front() || query.back() == txt.back()) continue;
            if (uni_0_1(double_gen) > p) continue;
            file << query << '\n';
            cnt++;
        }
    }
    std::cout << "(" << cnt
              << " tokenized queries generated in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n";
}


void experiment::Query::generate_dna_queries(std::string_view file_name, double p) {
    assert(is_empty(file_name));
    std::string _file_name {file_name};
    std::ofstream file(_file_name, std::ios_base::app);

    clock_t t0 = clock();

    std::random_device rd;
    std::default_random_engine double_gen(rd());
    // uniform distribution on the half-closed interval [0.0, 1.0)
    std::uniform_real_distribution<double> uni_0_1(0.0, 1.0);

    size_t cnt = 0;

    for (size_t i = 0; i < txt.size() - LEN_UB; i++) {
        if (uni_0_1(double_gen) > p) continue;
        for (size_t l = LEN_LB; l <= LEN_UB; l++) {
            if (uni_0_1(double_gen) > p) continue;
            file << txt.substr(i, l) << '\n';
            cnt++;
        }
    }

    std::cout << "(" << cnt
              << " DNA queries generated in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n";
}


void experiment::Result::update(size_t len, size_t freq, size_t net, size_t time) {
    len_time[len] += time;
    freq_time[freq] += time;
    nf_time[net] += time;
}


void experiment::SingleNF::run_queries(std::string_view query_file_name) {
    clock_t t0 = clock();
    size_t cnt = 0;

    std::string _query_file_name {query_file_name};
    std::ifstream file(_query_file_name);
    std::string query;
    while (std::getline(file, query)) {
        time_algo(query);
        cnt++;
    }
    std::cout << '(' << cnt << " queries ran in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s)\n ";
}


void experiment::SingleNF::time_algo(std::string_view s) {
    std::set < size_t > nfs;
    size_t nf;
    auto f = asa->freq(s);

    std::ofstream file(result_directory + "all.txt", std::ios_base::app);

    for (size_t a = 0; a < A; a++) {
        // start time t0 (cpu time and wall time, both in nanoseconds)
        struct timespec cpu_t0{};
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_t0);
        auto wall_t0 = std::chrono::high_resolution_clock::now();
        // some notes on timing:
        //
        // 1 second      = 1000 milliseconds (ms)
        // 1 millisecond = 1000 microseconds (\miu s)
        // 1 microsecond = 1000 nanoseconds  (ns)
        //
        // - clock() from time.h only supports precision up to microseconds
        //   (CLOCKS_PER_SEC = 1000000)
        // - clock_gettime() supports nanoseconds but only on linux
        // - timespec_get() measures wall time, not cpu time
        //
        // (ref: https://en.cppreference.com/w/c/chrono/clock,
        //      https://stackoverflow.com/a/12480485/7962088)
        switch ((Algo) a) {
            case Algo::ASA:
                nf = asa->single_nf(s);
                break;
            case Algo::HSA:
                nf = asa->single_nf_hash(s);
                break;
            case Algo::CSA:
                nf = csa->single_nf(s);
                break;
            case Algo::CRL:
                nf = asa->single_nf_crl(s);
                break;
            default:
                assert(false);
        }
        // end time t1
        struct timespec cpu_t1{};
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_t1);
        auto dur_wall = (std::chrono::high_resolution_clock::now() - wall_t0).count();
        // convert cpu_t0 and cpu_t1 to nanoseconds by adding the sec and nsec components
        auto t0ns = (size_t) 1e9 * cpu_t0.tv_sec + cpu_t0.tv_nsec;
        auto t1ns = (size_t) 1e9 * cpu_t1.tv_sec + cpu_t1.tv_nsec;
        auto dur_cpu = t1ns - t0ns;
        // no algorithms are faster a nanosecond
        // assert(dur_wall != 0 && dur_cpu != 0);
        nfs.insert(nf);
        results_cpu[a].update(s.size(), f, nf, dur_cpu);
        results_wall[a].update(s.size(), f, nf, dur_wall);
        file << a << '\t' << s.size() << '\t' << f << '\t' << nf << '\t' << dur_cpu << '\t' << dur_wall << '\n';
    }
    // all the algorithms have the same net frequency value
    assert(nfs.size() == 1);
    len_cnt[s.size()]++;
    freq_cnt[f]++;
    nf_cnt[nf]++;
}


void experiment::SingleNF::print_total_time() {
    std::cout << "\nalgo\tcpu\twall\n";
    for (size_t i = 0; i < A; i++) {
        uint64_t tot_cpu_time = 0, tot_wall_time = 0;
        for (auto [_, time]: results_cpu[i].freq_time) {
            tot_cpu_time += time;
        }
        for (auto [_, time]: results_wall[i].freq_time) {
            tot_wall_time += time;
        }
        std::cout << algo_name.at(Algo(i)) << '\t'
                  << (double) tot_cpu_time / 1e9 << "s\t"
                  << (double) tot_wall_time / 1e9 << "s\n";
    }
}


void experiment::SingleNF::write_results() {
    const std::vector<std::string> Q_str = {"len", "freq", "nf"};
    const std::vector<std::unordered_map<size_t, uint64_t>> Q_cnt = {len_cnt, freq_cnt, nf_cnt};
    for (size_t j = 0; j < Q_str.size(); j++) {
        // for each quantity (len, freq, or nf) has it own output file
        auto name = result_directory + Q_str[j] + ".txt";
        assert(is_empty(name));
        std::ofstream out_file(name, std::ios_base::app);
        for (size_t i = 0; i < A; i++) { // for each algorithm
            std::unordered_map<size_t, uint64_t> q_time_cpu, q_time_wall;
            if (j == 0) {
                q_time_cpu = results_cpu[i].len_time;
                q_time_wall = results_wall[i].len_time;
            }
            if (j == 1) {
                q_time_cpu = results_cpu[i].freq_time;
                q_time_wall = results_wall[i].freq_time;
            }
            if (j == 2) {
                q_time_cpu = results_cpu[i].nf_time;
                q_time_wall = results_wall[i].nf_time;
            }
            for (auto [q, tot_cpu]: q_time_cpu) {
                out_file << i << '\t' << q << '\t'
                         << tot_cpu << '\t' << q_time_wall.at(q) << '\t'
                         << Q_cnt[j].at(q) << '\n';
            }
        }
    }
}

void experiment::SingleNF::print_text_stats() {
    asa->get_text_stats(experiment::LEN_UB, result_directory + "dist.txt");
}


bool experiment::is_empty(std::string_view file_name) {
    std::string _file_name {file_name};
    std::ifstream file(_file_name);
    return file.peek() == std::ifstream::traits_type::eof();
}


void experiment::process_text(std::string_view in_file_name, std::string_view out_file_name) {
    std::string _in_file_name {in_file_name};
    std::ifstream in_file(_in_file_name);
    std::string txt;
    std::getline(in_file, txt);

    clock_t t0 = clock();

    auto before = txt.size();

    // remove all the punctuations
    txt.erase(std::remove_if(txt.begin(), txt.end(), ispunct), txt.end());

    // convert to lowercase
    std::transform(txt.begin(), txt.end(), txt.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    // reduce consecutive spaces to a single space
    std::regex re(R"(\s+)");
    txt = std::regex_replace(txt, re, " ");
    assert(txt.find("  ") == std::string::npos);

    std::cout << "(input text processed in "
              << ((double) (clock() - t0)) / CLOCKS_PER_SEC << "s;\n";

    auto after = txt.size();
    std::cout << " " << before << " - " << after << " = " << before - after << " characters removed)\n";

    assert(is_empty(out_file_name));
    std::string _out_file_name {out_file_name};
    std::ofstream out_file(_out_file_name, std::ios_base::app);
    out_file << txt;
}


void experiment::run() {
    std::cout << std::fixed << std::setprecision(3);

    const std::string input_file_name = "dna.txt";

    std::ifstream in_file("../data/" + input_file_name);
    std::string txt((std::istreambuf_iterator<char>(in_file)),
                    (std::istreambuf_iterator<char>()));

    // use half of the text to generate queries,
    // the other half to build the data structures
    auto txt0 = char(2) + txt.substr(0, txt.size() / 2) + char(3);
    txt = char(2) + txt.substr(txt.size() / 2) + char(3);

    const std::string query_dir_name = "../queries/";
    const std::string query_file_name = query_dir_name + input_file_name;
    const std::string output_dir_name = "../output/";

    std::filesystem::create_directory(query_dir_name);
    std::filesystem::create_directory(output_dir_name);

    Query q{txt0};
    q.generate_dna_queries(query_file_name, 0.01);
//    q.generate_tokenized_queries(query_file_name, 0.01);

    SingleNF instance{txt, output_dir_name};
    instance.print_text_stats();
    instance.run_queries(query_file_name);
    instance.print_total_time();
    instance.write_results();
}

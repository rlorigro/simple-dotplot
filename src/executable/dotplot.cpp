#include "IterativeSummaryStats.hpp"
#include "FastqIterator.hpp"
#include "FastaReader.hpp"
#include "Color.hpp"
#include "CLI11.hpp"

using simple_dotplot::FastqElement;
using simple_dotplot::FastqIterator;

#include <Filesystem.hpp>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <map>

using ghc::filesystem::create_directories;
using ghc::filesystem::path;
using std::sort;
using std::runtime_error;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::string;
using std::vector;
using std::array;
using std::cerr;
using std::cout;
using std::map;

#include "mummer/sparseSA.hpp"


char complement(char c){
    c = char(toupper(c));

    if (c == 'A'){
        return 'T';
    }
    else if (c == 'C'){
        return 'G';
    }
    else if (c == 'G'){
        return 'C';
    }
    else if (c == 'T'){
        return 'A';
    }
    else {
        throw runtime_error("ERROR: non DNA character cannot be complemented: " + string(1,c));
    }
}


void reverse_complement(string& sequence, string& rc_sequence){
    rc_sequence.resize(0);

    for (std::string::reverse_iterator iter=sequence.rbegin(); iter!=sequence.rend(); ++iter){
        rc_sequence += complement(*iter);
    }
}


int64_t find_matches(string& ref_sequence,
                  string& query_sequence,
                  int64_t min_length,
                  vector<mummer::mummer::match_t>& matches
                  ){

    auto matcher = mummer::mummer::sparseSA::create_auto(ref_sequence.c_str(), ref_sequence.size(), 0, true);

    IterativeSummaryStats<double> accumulator;

    matcher.findMEM_each(query_sequence, min_length, false, [&](const mummer::mummer::match_t& match){
        matches.emplace_back(match);

        for (int64_t m=0; m<match.len; m+=50) {
            accumulator.add(double(match.ref)/1000.0);
        }
    });

    return accumulator.get_mean();
}


void plot_matches_as_heatmap(
        const map <int64_t, vector<mummer::mummer::match_t> >& matches,
        const map <int64_t, string>& names,
        const map <int64_t, size_t>& query_sizes,
        size_t ref_size,
        bool mask_diagonal,
        path output_path){

    const size_t size = 1000;
    vector <vector <size_t> > matrix(size, vector<size_t>(size, 0));

    int64_t x_offset = 0;
    int64_t y_offset = 0;

    int64_t query_size_sum = 0;
    for (const auto& item: query_sizes){
        query_size_sum += item.second;
    }

    int64_t total_size = std::max(query_size_sum, int64_t(ref_size));

    if (total_size > int64_t(ref_size)) {
        cerr << "WARNING: rescaling because sum of query sizes (" << total_size
             << ") is greater than ref size (" << ref_size << ")" << '\n';
    }

    for (const auto& item: matches){

        cerr << "Plotting " << names.at(item.first) << " with length " << query_sizes.at(item.first) << '\n';
        cerr << "Matches found: " << item.second.size() << '\n';

        for (const auto& match: item.second){
            // If the user wants to ignore matches on the diagonal, then skip
            if (match.ref == match.query and mask_diagonal){
                continue;
            }

//            cerr << "\n---- " << match.ref << " " << match.query << " " << match.len << " ----\n";

            double m=0.00001;
            while (true) {
                m += double(total_size)/double(size);

                auto x = double(match.ref) + m + double(x_offset);
                auto y = double(match.query) + m + double(y_offset);

                size_t x_bin = round((double(x)/double(total_size))*double(size)) - 1;
                size_t y_bin = round((double(y)/double(total_size))*double(size)) - 1;

//                cerr << m << " " << match.len << " " << x << " " << y << " " << x_bin << " " << y_bin << '\n';

                if (x_bin >= matrix.size()){
                    break;
                }
                else{
                    if (y_bin >= matrix[x_bin].size()){
                        break;
                    }
                }

                if (size_t(match.len) > matrix[x_bin][y_bin]){
                    matrix[x_bin][y_bin] = match.len;
                }

                if (m > match.len){
                    break;
                }
            }
        }

        y_offset += query_sizes.at(item.first);
    }

    cerr << "total length: " << y_offset << '\n';

    size_t max_observed_match_length = 0;
    for (size_t i=0; i<size; i++){
        for (size_t j=0; j<size; j++){
            if (matrix[i][j] > max_observed_match_length){
                max_observed_match_length = matrix[i][j];
            }
        }
    }

    cerr << "Maximum length match observed: " << max_observed_match_length << '\n';

    Viridis color_map;

    ofstream file(output_path);

    if (not (file.good() and file.is_open())){
        throw runtime_error("ERROR: could not write to file: " + output_path.string());
    }


    file << "P3" << '\n';
    file << size << ' ' << size << '\n';
    file << "255" << '\n';
    for (size_t i=0; i<size; i++){
        for (size_t j=0; j<size; j++){
            auto match_density = matrix[i][j];
            auto normalized_density = log(double(match_density) + 1)/log(double(max_observed_match_length) + 1);
            auto color = color_map.get_rgb(normalized_density);

            file << int(color[0]*255) << ' ' << int(color[1]*255) << ' ' << int(color[2]*255) << "    ";
//            cerr << int(match_density) << ' ' << int(color[0]*255) << ' ' << int(color[1]*255) << ' ' << int(color[2]*255) << '\n';
        }
        file << '\n';
    }
}


template <class T, class T2> void fasta_or_fastq_dotplot(
        T& ref_reader,
        T2& query_reader,
        uint32_t min_length,
        bool mask_diagonal,
        path output_directory,
        bool autosort=false){

    auto ref_element = ref_reader.generate_sequence_container();
    auto query_element = query_reader.generate_sequence_container();

    vector<string> ref_sequences;
    vector<string> ref_names;

    vector<string> query_sequences;
    vector<string> query_names;

    cerr << "Loading sequences... " << '\n';
    while (ref_reader.next_element(ref_element)){
        ref_sequences.emplace_back(ref_element.sequence);
        ref_names.emplace_back(ref_element.name);
    }

    while (query_reader.next_element(query_element)){
        query_sequences.emplace_back(query_element.sequence);
        query_names.emplace_back(query_element.name);
    }

    for (size_t i=0; i<ref_sequences.size(); i++){
        cerr << "Constructing suffix array for ref " << ref_names[i] << "..." << '\n';
        map <int64_t, vector <mummer::mummer::match_t> > sorted_query_matches;
        map <int64_t, string> sorted_query_names;
        map <int64_t, size_t> sorted_query_sizes;

        size_t sum_of_query_sizes = 0;

        // Accumulate matches for each query sequence
        for (size_t j=0; j<query_sequences.size(); j++){
            cerr << "Finding matches for " << query_names[j] << " VS " << ref_names[i] << "..." << '\n';

            vector <mummer::mummer::match_t> matches;

            int64_t mean_x_value = find_matches(ref_sequences[i], query_sequences[j], min_length, matches);

            if (matches.empty()){
                cerr << "WARNING: no matches found for query: " << query_names[j] << '\n';
            }

            cerr << "mean x value: " << mean_x_value << '\n';

            int64_t ordinal;
            if (autosort){
                ordinal = mean_x_value;
            }
            else{
                ordinal = j;
            }

            sorted_query_matches.emplace(ordinal, matches);
            sorted_query_names.emplace(ordinal, query_names[j]);
            sorted_query_sizes.emplace(ordinal, query_sequences[j].size());

            sum_of_query_sizes += query_sequences[j].size();
        }

//        string name = "all_VS_" + ref_names[i] + "_dotplot.svg";
//        path plot_path = output_directory / name;
//
//        plot_matches(all_matches, all_names, ref_sequences[i].size(), sum_of_query_sizes, plot_path);

        auto query_prefix = query_reader.file_path.stem();
        string name = query_prefix.string() + "_VS_" + ref_names[i] + "_dotplot.ppm";
        path plot_path = output_directory / name;

        cerr << "Plotting..." << '\n';
        plot_matches_as_heatmap(sorted_query_matches, sorted_query_names, sorted_query_sizes, ref_sequences[i].size(), mask_diagonal, plot_path);
    }
}


void dotplot(path ref_path, path query_path, int64_t min_length, bool mask_diagonal, path output_directory){

    if (exists(output_directory)){
        throw runtime_error("ERROR: output directory already exists");
    }
    else{
        create_directories(ghc::filesystem::absolute(output_directory));
    }

    if (ref_path.extension() == ".fastq" and query_path.extension() == ".fastq"){
        FastqIterator ref_reader(ref_path);
        FastqIterator query_reader(query_path);

        fasta_or_fastq_dotplot(ref_reader, query_reader, min_length, mask_diagonal, output_directory);
    }
    else if (ref_path.extension() == ".fasta" and query_path.extension() == ".fasta"){
        FastaReader ref_reader(ref_path);
        FastaReader query_reader(query_path);

        fasta_or_fastq_dotplot(ref_reader, query_reader, min_length, mask_diagonal, output_directory);
    }
    else if (ref_path.extension() == ".fasta" and query_path.extension() == ".fastq"){
        FastaReader ref_reader(ref_path);
        FastqIterator query_reader(query_path);

        fasta_or_fastq_dotplot(ref_reader, query_reader, min_length, mask_diagonal, output_directory);
    }
    else if (ref_path.extension() == ".fastq" and query_path.extension() == ".fasta"){
        FastqIterator ref_reader(ref_path);
        FastaReader query_reader(query_path);

        fasta_or_fastq_dotplot(ref_reader, query_reader, min_length, mask_diagonal, output_directory);
    }
    else{
        throw runtime_error("ERROR: file extensions do not match '.fastq' or '.fasta'");
    }
}


int main(int argc, char* argv[]){
    path ref_path;
    path query_path;
    int64_t min_length;
    bool mask_diagonal = false;
    path output_directory;

    CLI::App app{"App description"};

    std::string filename = "default";
    app.add_option("-r,--ref", ref_path, "Path to FASTA/Q of ref sequences")
        ->required();

    app.add_option("-q,--query", query_path, "Path to FASTA/Q of query sequences")
        ->required();

    app.add_option("-l,--min_length", min_length, "Minimum length of match to be considered")
        ->required();

    app.add_flag("--mask_diagonal", mask_diagonal, "Whether to zero all diagonal entries. \nUseful for self-dotplots, which have an extreme outlier on the diagonal");

    app.add_option("-o,--output_dir", output_directory, "Path of directory to save output")
        ->required();

    CLI11_PARSE(app, argc, argv);

    dotplot(ref_path, query_path, min_length, mask_diagonal, output_directory);

    return 0;
}


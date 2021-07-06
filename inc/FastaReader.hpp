#ifndef SIMPLE_DOTPLOT_FASTAREADER_HPP
#define SIMPLE_DOTPLOT_FASTAREADER_HPP

#include <Filesystem.hpp>
#include "SequenceElement.hpp"
#include <unordered_map>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using ghc::filesystem::create_directories;
using ghc::filesystem::path;
using std::unordered_map;
using std::pair;
using std::vector;
using std::string;
using std::ifstream;


class FastaIndex {
public:
    uint64_t byte_index;    // Where in the fasta is the start of the sequence
    uint64_t length;        // How long is the sequence

    FastaIndex(uint64_t byte_index, uint64_t length);
    FastaIndex();

    uint64_t size();
};


void get_vector_from_index_map(vector< pair <string,FastaIndex> >& items, unordered_map<string,FastaIndex>& map_object);


class FastaReader{
public:
    /// Attributes ///
    path file_path;
    path index_path;
    char header_symbol;
    string eof_placeholder_name;
    uint64_t line_index;
    bool end_of_file;
    vector <FastaIndex> sequential_read_offsets;
    unordered_map <string, size_t > read_indexes_by_name;

    /// Methods ///
    FastaReader(path file_path);
    FastaReader();

    // Fetch header + sequence assuming read position precedes a header
    bool next_element(SequenceElement& element);

    // Generate index if needed, otherwise load index
    void index();

    // Return a copy of the read indexes
    void get_indexes_mapped_by_name(unordered_map <string, FastaIndex>& indexes);

    // Set this reader's index-related attributes to match the donor FastaReader
    void adopt_index(FastaReader& donor);

    // Fetch a sequence from file, and generate index first if necessary
    void get_sequence(SequenceElement& element, string& sequence_name);
    void get_sequence(SequenceElement& element, string& sequence_name, uint64_t fasta_byte_index);

    // Create a copy of the appropriate sequence container for this reader
    SequenceElement generate_sequence_container();

    // Use htslib to make an .fai
    void build_fasta_index();

private:
    /// Attributes ///
    ifstream fasta_file;

    /// Methods ///
    // Assuming the read position of the ifstream is at a header, parse the header line
    void read_next_header(SequenceElement& element);

    // Assuming the read position of the ifstream is at a sequence, parse the sequence line(s),
    // which precede the next header
    void read_next_sequence(SequenceElement& element);

    // Read fai
    void read_fasta_index();
};


#endif //SIMPLE_DOTPLOT_FASTAREADER_HPP

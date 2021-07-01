#ifndef SIMPLE_DOTPLOT_FASTAWRITER_HPP
#define SIMPLE_DOTPLOT_FASTAWRITER_HPP

#include "FastaReader.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <experimental/filesystem>

using std::string;
using std::to_string;
using std::cout;
using std::ofstream;
using std::vector;
using std::runtime_error;
using std::experimental::filesystem::path;
using std::experimental::filesystem::create_directories;


class FastaWriter{
    /// ATTRIBUTES ///
    path file_path;
    ofstream fasta_file;
    char header_symbol;
    bool file_end;
    uint64_t index;

    /// METHODS ///
public:
    FastaWriter(path file_path);
    template<class T> void write(T& sequence);
    template<class T> void write(vector<T>& sequences);
};


template<class T> void FastaWriter::write(T& sequence){
    this->fasta_file << string(1, this->header_symbol);
    this->fasta_file << sequence.name;
    this->fasta_file << "\n";
    this->fasta_file << sequence.sequence;
    this->fasta_file << "\n";
}


template<class T> void FastaWriter::write(vector<T>& sequences){
    for (auto& element: sequences){
        this->write(element);
    }
}

#endif //SIMPLE_DOTPLOT_FASTAWRITER_HPP

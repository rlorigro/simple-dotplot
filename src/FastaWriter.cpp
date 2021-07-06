#include "FastaWriter.hpp"

#include <string>
#include <stdexcept>
#include <Filesystem.hpp>

using std::string;
using std::to_string;
using std::cout;
using std::ofstream;
using std::vector;
using std::runtime_error;
using ghc::filesystem::create_directories;
using ghc::filesystem::path;


FastaWriter::FastaWriter(path file_path) {
    this->file_path = file_path;
    this->fasta_file = ofstream(file_path);
    this->header_symbol = '>';
    this->index = 0;

    // Check if file is readable or exists
    if (not this->fasta_file.is_open()){
        throw runtime_error("ERROR: file write error: " + string(file_path));
    }
}


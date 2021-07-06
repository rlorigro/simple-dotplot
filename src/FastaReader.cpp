#include "FastaReader.hpp"
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <exception>
#include <Filesystem.hpp>
#include "boost/algorithm/string.hpp"

using std::unordered_map;
using std::runtime_error;
using std::out_of_range;
using std::string;
using std::vector;
using std::to_string;
using std::cout;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::runtime_error;
using ghc::filesystem::path;
using ghc::filesystem::exists;
using ghc::filesystem::absolute;
using boost::trim_left_if;
using boost::trim_right;
using boost::split;


void get_vector_from_index_map(vector< pair <string,FastaIndex> >& items, unordered_map<string,FastaIndex>& map_object){
    for (auto& item: map_object){
    items.emplace_back(item.first, item.second);
    }
}


bool is_caret(char c){
    return (c == '>');
}


bool is_tab(char c){
    return (c == '\t');
}


bool is_comma(char c){
    return (c == ',');
}


FastaIndex::FastaIndex(uint64_t byte_index, uint64_t length){
    this->byte_index = byte_index;
    this->length = length;
}


uint64_t FastaIndex::size(){
    return this->length;
}


FastaIndex::FastaIndex() = default;

FastaReader::FastaReader() = default;

FastaReader::FastaReader(path file_path){
    file_path = absolute(file_path);
    this->file_path = file_path;
    this->fasta_file = ifstream(file_path);
    this->header_symbol = '>';
    this->eof_placeholder_name = ">EOF<";
    this->end_of_file = false;
    this->line_index = 0;

    // Check if file is readable or exists
    if (not this->fasta_file.good()){
        throw runtime_error("ERROR: file read error: " + this->file_path.string());
    }
}


SequenceElement FastaReader::generate_sequence_container(){
    return SequenceElement();
}


// Fetch header which starts at current ifstream position
void FastaReader::read_next_header(SequenceElement& element){
    // Try fetching header line
    getline(this->fasta_file, element.name);
    this->line_index++;

    // Update class status
    if (this->fasta_file.peek() == EOF){
        this->end_of_file = true;
    }

    // Verify that header starts with header symbol '>', and trim the symbol from the line. Also trim trailing space.
    if (element.name[0] == this->header_symbol){
        trim_left_if(element.name, &is_caret);
//        trim_right(element.name);
        element.name = element.name.substr(0,element.name.find_first_of(' '));
    }
    else{
        throw runtime_error("Unrecognized FASTA header character on line: " + to_string(this->line_index));
    }
}


// Fetch sequence which starts at current ifstream position
void FastaReader::read_next_sequence(SequenceElement& element){
    // Make sure the element is empty before starting
    if (element.sequence.size() != 0){
        element = {};
    }

    // Try fetching sequence line (of which there may be multiple)
    while (this->fasta_file.peek() != std::char_traits<char>::to_int_type(this->header_symbol) && !this->end_of_file) {
        // If this is a single line, just place it in the sequence object
        if (element.sequence.length() == 0){
            getline(this->fasta_file, element.sequence);
        }
        // If this is a multi-line FASTA, append the sequence object using a copy of the line (sadly)
        else{
            string buffer;
            getline(this->fasta_file, buffer);
            element.sequence += buffer;
        }
        this->line_index++;

        // Trim any trailing whitespace on the sequence line
        trim_right(element.sequence);

        // Update class status
        if (this->fasta_file.peek() == EOF){
            this->end_of_file = true;
        }
    }
}


// Read the next lines of the file and update the SequenceElement object. If no lines exists, toggle this->file_end
bool FastaReader::next_element(SequenceElement& element){
    if (this->end_of_file) {
        return false;
    }

    element = {};
    this->read_next_header(element);
    this->read_next_sequence(element);

    return true;
}

#include "FastqIterator.hpp"
#include <cmath>

using std::runtime_error;
using std::pow;

namespace simple_dotplot {

double quality_char_to_error_probability(char q){
    return pow(10,double(q - 33) / -10.0);
}


void FastqElement::get_error_probabilities(vector<float>& probabilities) {
    probabilities.resize(0);
    for (auto& q: quality_string) {
        auto p = float(quality_char_to_error_probability(q));
        probabilities.emplace_back(p);
    }
}


FastqIterator::FastqIterator(path file_path):
        file_path(file_path),
        file(file_path)
{
    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: fastq file cannot be opened: " + absolute(file_path).string());
    }
}


size_t FastqIterator::get_line_number() const{
    return line_number;
}


FastqElement FastqIterator::generate_sequence_container(){
    return FastqElement();
}


bool FastqIterator::next_element(FastqElement& element){
    if (file.eof() or file.peek() == EOF){
        return false;
    }

    element = {};

    // Get name (trim delimiter)
    getline(file, element.name);
    element.name = element.name.substr(1,element.name.size()-1);

    if (file.eof()){
        throw runtime_error("Error: FASTQ file terminates early: " + file_path.string());
    }

    // Get sequence
    getline(file, element.sequence);

    if (file.eof()){
        throw runtime_error("Error: FASTQ file terminates early: " + file_path.string());
    }

    // Skip useless 2nd header line
    file.ignore(1, '\n');

    if (file.eof()){
        throw runtime_error("Error: FASTQ file terminates early: " + file_path.string());
    }

    // Get quality string
    getline(file, element.quality_string);

    return true;
}

}
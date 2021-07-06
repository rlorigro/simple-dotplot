#ifndef SIMPLE_DOTPLOT_FASTQITERATOR_HPP
#define SIMPLE_DOTPLOT_FASTQITERATOR_HPP

#include <Filesystem.hpp>
#include <fstream>
#include <string>
#include <vector>

using ghc::filesystem::create_directories;
using ghc::filesystem::path;
using std::ifstream;
using std::to_string;
using std::string;
using std::vector;


namespace simple_dotplot {

double quality_char_to_error_probability(char q);


class FastqElement{
public:
    /// Attributes ///
    string name;
    string sequence;
    string quality_string;

    /// Methods ///
    FastqElement()=default;
    void get_error_probabilities(vector<float>& probabilities);
};


class FastqIterator{
public:
    /// Attributes ///
    const path file_path;

    /// Methods ///
    FastqIterator(path file_path);

    bool next_element(FastqElement& element);
    size_t get_line_number() const;

    FastqElement generate_sequence_container();

private:
    /// Attributes ///
    bool eof;
    ifstream file;
    size_t line_number;

};

}

#endif //SIMPLE_DOTPLOT_FASTQITERATOR_HPP

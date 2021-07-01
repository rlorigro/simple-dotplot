#ifndef SIMPLE_DOTPLOT_PATH_HPP
#define SIMPLE_DOTPLOT_PATH_HPP

#include <string>

using std::string;


class Path {
public:
    string path;

    Path(string p);
    void append(Path b);
    void append(string b);

    operator string() const{
        return path;
    };
};


#endif //SIMPLE_DOTPLOT_PATH_HPP

#ifndef SIMPLE_DOTPLOT_PATH_HPP
#define SIMPLE_DOTPLOT_PATH_HPP

#include <functional>
#include <string>

using std::function;
using std::string;


class Path {
public:
    string path;

    Path(string p);
    void append(Path b);
    void append(string b);
    void for_each_directory(const function<void(Path& dir_path)>& f);

    operator string() const{
        return path;
    };
};


void create_directories(Path p);


void create_directories(string s);


bool directory_exists(const Path& p);


#endif //SIMPLE_DOTPLOT_PATH_HPP

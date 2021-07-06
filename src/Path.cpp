#include "Path.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>

using std::runtime_error;
using std::vector;
using std::string;
using std::cerr;


Path::Path(string p):
    path(p)
{
    for (char& c: path) {
        if (not isalnum(c)) {
            if ((c != '.') and (c != '_') and (c != '/')) {
                throw runtime_error("ERROR: path must not contain disallowed character: " + string(1,c));
            }
        }
    }

    if (path.back() == '/'){
        path.resize(path.size()-1);
    }
}


void Path::append(Path b){
    if (b.path.empty()){
        return;
    }

    if (path.empty()){
        path += b.path;
    }
    else{
        if (b.path[0] == '/'){
            path += b.path.substr(1,b.path.size()-2);
        }
        else {
            path += '/';
            path += b.path;
        }
    }
}


void Path::append(string b){
    Path b_path(b);

    append(b_path);
}


void operator/(Path& a, Path& b){
    a.append(b);
}


void Path::for_each_directory(const function<void(Path& dir_path)>& f){
    vector<string> directories;
    string token;

    for (char& c: path){
        if (c == '/'){
            directories.push_back(token);
            token.resize(0);
        }
        else {
            token += c;
        }
    }
    directories.push_back(token);

    Path output_path("");
    for (auto& item: directories){
        output_path.append(item);
        f(output_path);
    }
}


void create_directories(Path p){
    p.for_each_directory([&](Path& dir_path){
        cerr << dir_path.path << '\n';
    });
    cerr << '\n';
}


void create_directories(string s){
    Path p(s);

    create_directories(p);
}


// ADAPTED FROM https://stackoverflow.com/a/52043954/5095925
/******************************************************************************
 * Checks to see if a directory exists. Note: This method only checks the
 * existence of the full p AND if p leaf is a dir.
 *
 * @return  >0 if dir exists AND is a dir,
 *           0 if dir does not exist OR exists but not a dir,
 *          <0 if an error occurred (errno is also set)
 *****************************************************************************/
bool directory_exists(const Path& p){
    struct stat info;
    int statRC = stat(p.path.c_str(), &info);

    if(statRC != 0){
        if (errno == ENOENT){
            return false; // something along the p does not exist
        }
        if (errno == ENOTDIR){
            return false; // something in p prefix is not a dir
        }

        throw runtime_error("ERROR: unspecified filesystem error in verifying p: " + p.path);
    }

    return (info.st_mode & S_IFDIR);
}

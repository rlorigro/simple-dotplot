#include "Path.hpp"



Path::Path(string p):
    path(p)
{
    if (path.back() == '/'){
        path.resize(path.size()-1);
    }
}


void Path::append(Path b){
    if (b.path.empty()){
        return;
    }

    if (b.path[0] == '/'){
        path += b.path.substr(1,b.path.size()-2);
    }
    else {
        path += '/';
        path += b.path;
    }
}


void Path::append(string b){
    Path b_path(b);

    append(b_path);
}


void operator/(Path& a, Path& b){
    a.append(b);
}


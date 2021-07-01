#include "Path.hpp"
#include <iostream>
#include <string>

using std::string;
using std::cerr;


int main(){
    Path p("/home/ryan/data/test/a/");

    p.append("b");
    p.append("c");

    cerr << string(p) << '\n';

    return 0;
}

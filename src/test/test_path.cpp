#include "Path.hpp"
#include <iostream>
#include <string>

using std::string;
using std::cerr;


int main(){
    Path p("/home/ryan/data/test/a/");


    try {
        string s = "0";

        cerr << "TESTING: " << s << '\n' << std::flush;
        Path good_path(s);
    }
    catch (const std::runtime_error& e){
        cerr << "Error caught:\n\t" <<  e.what() << '\n' << "... end of error\n\n";
    }

    try {
        string s = ".";

        cerr << "TESTING: " << s << '\n' << std::flush;
        Path good_path(s);
    }
    catch (const std::runtime_error& e){
        cerr << "Error caught:\n\t" <<  e.what() << '\n' << "... end of error\n\n";
    }

    try {
        string s = "_";

        cerr << "TESTING: " << s << '\n' << std::flush;
        Path good_path(s);
    }
    catch (const std::runtime_error& e){
        cerr << "Error caught:\n\t" <<  e.what() << '\n' << "... end of error\n\n";
    }

    try {
        string s = ")";

        cerr << "TESTING: " << s << '\n' << std::flush;
        Path good_path(s);
    }
    catch (const std::runtime_error& e){
        cerr << "Error caught:\n\t" <<  e.what() << '\n' << "... end of error\n\n";
    }

    try {
        string s = "*";

        cerr << "TESTING: " << s << '\n' << std::flush;
        Path good_path(s);
    }
    catch (const std::runtime_error& e){
        cerr << "Error caught:\n\t" <<  e.what() << '\n' << "... end of error\n\n";
    }

    p.append("b");
    p.append("c");

    cerr << string(p) << '\n';

    create_directories("/home/ryan/data/test/");

    return 0;
}

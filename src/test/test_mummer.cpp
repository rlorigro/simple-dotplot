#include "mummer/sparseSA.hpp"

#include <string>
#include <vector>
#include <iostream>

using std::string;
using std::vector;
using std::cerr;


int main (){
    string ref = "gattacagattacagattacacgatctgagcatgc";

    auto matcher = mummer::mummer::sparseSA::create_auto(ref.c_str(), ref.size(), 0, true);

    vector<mummer::mummer::match_t> mems;
    vector<mummer::mummer::match_t> mams;
    vector<mummer::mummer::match_t> mums;

    string query = "attacac";

    cerr << "MEMs" << '\n';

    // min_len must be > 1
    matcher.findMEM_each(query, 2, false, [&](const mummer::mummer::match_t& match){
        cerr << match.ref << ' ' << match.len << ' ' << ref.substr(match.ref, match.len) << '\n';

        mems.emplace_back(match);
    });

    cerr << "MAMs" << '\n';

    // min_len must be > 1
    matcher.findMAM_each(query, 2, false, [&](const mummer::mummer::match_t& match){
        cerr << match.ref << ' ' << match.len << ' ' << ref.substr(match.ref, match.len) << '\n';

        mams.emplace_back(match);
    });

    cerr << "MUMs" << '\n';

    // min_len must be > 1
    matcher.findMUM_each(query, 2, false, [&](const mummer::mummer::match_t& match){
        cerr << match.ref << ' ' << match.len << ' ' << ref.substr(match.ref, match.len) << '\n';

        mums.emplace_back(match);
    });

    cerr << "done" << '\n';

    return 0;
}


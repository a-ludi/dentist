#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <iostream>
#include <vector>
#include <exception>

using namespace sdsl;
using namespace std;


#ifndef INCLUDED_FM_INDEX_HPP
#define INCLUDED_FM_INDEX_HPP

int main(int argc, char** argv);
void buildIndex(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex, string referenceFile);
vector<size_t> getRecordStarts(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex);
int hasStdin();
void locateQueries(
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex,
    vector<size_t> recordStarts,
    istream &queries,
    string sourceName,
    string queryBuffer
);
size_t locateQuery(
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex,
    vector<size_t> recordStarts,
    string sourceName,
    size_t queryId,
    string query
);
size_t findSourceId(vector<size_t> recordStarts, size_t hitBegin);


class FmIndexException : public exception {
    const string msg;

public:
    FmIndexException(const string msg) : msg(msg)
    {
    }

    const char* what() const noexcept {
        return msg.c_str();
    }
};

#endif

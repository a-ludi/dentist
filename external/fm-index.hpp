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
void printUsage();
int parseArgs(int argc, char** argv, string *tempDir, bool *includeReverseComplement);
void buildIndex(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex, string referenceFile);
vector<size_t> getRecordStarts(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex, string referenceFile);
int hasStdin();
void locateQueries(
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex,
    vector<size_t> recordStarts,
    istream &queries,
    string sourceName,
    string queryBuffer,
    bool reverseComplement
);
size_t locateQuery(
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex,
    vector<size_t> recordStarts,
    string sourceName,
    size_t queryId,
    string query,
    bool reverseComplement
);
void getReverseComplement(string &revCompBuffer, string query);
size_t findSourceId(vector<size_t> recordStarts, size_t hitBegin);

template<class int_type, class t_int_vec>
bool load_from_plain_array(t_int_vec& v, const std::string& file)
{
    ifstream in(file, std::ios::binary | std::ios::in);

    if (!in)
        return false;

    // get length of file
    in.seekg (0, in.end);
    int fileSize = in.tellg();
    in.seekg (0, in.beg);

    if (fileSize % sizeof(int_type) != 0)
        return false;

    v.resize(fileSize / sizeof(int_type));

    in.read((char*) v.data(), fileSize);

    return (bool) in;
}

bool dirExists(const string& dir);

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

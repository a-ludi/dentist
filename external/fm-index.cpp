#include "fm-index.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <ctime>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <exception>
#include <poll.h>

using namespace sdsl;
using namespace std;



int main(int argc, char** argv)
{
    if (argc <  2)
    {
        cerr << "Usage " << argv[0] << " <in:reference> [<in:queries> ...]" << endl;
        cerr << "    This program constructs a very compact FM-index" << endl;
        cerr << "    of <reference> and locates <queries> if given." << endl;
        cerr << endl;
        cerr << "Positional arguments:" << endl;
        cerr << "    <in:refernce>  Original text file with one record per line." << endl;
        cerr << "    <in:queries>   List of queries (one per line) to locate in <reference>." << endl;
        cerr << "                   Queries given on standard input will be located before" << endl;
        cerr << "                   all others." <<endl;
        cerr << endl;
        cerr << "Output:" << endl;
        cerr << "    Creates an FM-index <reference>.fm9 if not present. Then produces" << endl;
        cerr << "    a list of exact matches for all queries in a TAB-separated format:" << endl;
        cerr << endl;
        cerr << "        refName  refId  refLength  queryId  hitBegin  hitEnd" << endl;
        cerr << endl;
        cerr << "    IDs and coordinates are zero-based. Coordinates are right-open." << endl;
        cerr << endl;
        cerr << "    Note: no output will be produced if no queries are given." << endl;

        return 1;
    }

    try
    {
        string referenceFile = string(argv[1]);
        char** queriesFiles = &argv[2];
        int numQueries = max(0, argc - 2);

        csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fmIndex;
        buildIndex(fmIndex, referenceFile);

        auto recordStarts = getRecordStarts(fmIndex);
        string queryBuffer;
        // Reserve 10 Mb for a single query
        queryBuffer.reserve(10 * (2 << 20));

        if (cin && hasStdin())
            locateQueries(fmIndex, recordStarts, cin, "stdin", queryBuffer);

        for (int i = 0; i < numQueries; ++i)
        {
            string queryFile(queriesFiles[i]);
            ifstream queryData(queryFile);

            if (!queryData)
            {
                cerr << "{\"level\":\"warning\","
                     << "\"info\":\"File does not exist. Skipping.\","
                     << "\"file\":\"" << queryFile << "\"}" << endl;

                continue;
            }

            locateQueries(fmIndex, recordStarts, queryData, queryFile, queryBuffer);
        }
    }
    catch (FmIndexException& e)
    {
        cerr << "error: " << e.what() << endl;

        return 2;
    }
    catch (exception& e)
    {
        cerr << "critical error: " << e.what() << endl;

        return 3;
    }

    return 0;
}

void buildIndex(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex, string referenceFile)
{
    string indexSuffix = ".fm9";
    string indexFile = referenceFile + indexSuffix;

    if (!load_from_file(fmIndex, indexFile)) {
        ifstream referenceData(referenceFile.c_str());

        if (!referenceData)
            throw new FmIndexException("File `" + referenceFile + "` does not exist.");

        cerr << "{\"level\":\"info\","
             << "\"info\":\"Index does not exist. Building it now.\","
             << "\"file\":\"" << indexFile << "\""
             << "}" << endl;

        time_t start = time(NULL);

        construct(fmIndex, referenceFile, 1); // generate index
        store_to_file(fmIndex, indexFile); // save it

        cerr << "{\"level\":\"info\","
             << "\"info\":\"Built index.\","
             << "\"file\":\"" << indexFile << "\","
             << "\"elpasedSecs\":" << difftime(time(NULL), start) << ","
             << "\"sizeMiB\":" << size_in_mega_bytes(fmIndex)
             << "}" << endl;
    }
}

vector<size_t> getRecordStarts(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex)
{
    auto lineEndLocations = locate(fmIndex, "\n");
    sort(lineEndLocations.begin(), lineEndLocations.end());
    vector<size_t> recordStarts(lineEndLocations.size() + 1);

    recordStarts[0] = 0;
    for (int i = 1; i < recordStarts.size(); ++i)
        recordStarts[i] = lineEndLocations[i - 1] + 1;

    return recordStarts;
}

int hasStdin()
{
    int ret;
    struct pollfd pfd[1] = {0};

    pfd[0].fd = STDIN_FILENO;
    pfd[0].events = POLLIN;
    ret = poll(pfd, 1, 0);

    return (ret>0);
}

void locateQueries(
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex,
    vector<size_t> recordStarts,
    istream &queriesData,
    string sourceName,
    string queryBuffer
)
{
    cerr << "{\"level\":\"info\","
         << "\"info\":\"Processing queries.\","
         << "\"source\":\"" << sourceName << "\""
         << "}" << endl;

    time_t start = time(NULL);
    size_t queryId = 0;
    while (queriesData && !queriesData.eof())
    {
        getline(queriesData, queryBuffer);
        if (queryBuffer.size() > 0)
            locateQuery(fmIndex, recordStarts, sourceName, queryId++, queryBuffer);

        if (sourceName == "stdin" && !hasStdin())
            break;
    }

    cerr << "{\"level\":\"info\","
         << "\"info\":\"Finished queries.\","
         << "\"source\":\"" << sourceName << "\","
         << "\"elpasedSecs\":" << difftime(time(NULL), start)
         << "}" << endl;
}

void locateQuery(
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex,
    vector<size_t> recordStarts,
    string sourceName,
    size_t queryId,
    string query
)
{
    auto locations = locate(fmIndex, query);
    sort(locations.begin(), locations.end());

    for (size_t i = 0; i < locations.size(); ++i)
    {
        auto hitBegin = locations[i];
        auto hitEnd = hitBegin + query.size();
        auto sourceId = findSourceId(recordStarts, hitBegin);
        auto sourceBegin = recordStarts[sourceId];
        // NOTE: length of line terminator must be subtracted
        auto sourceLength = recordStarts[sourceId + 1] - recordStarts[sourceId] - 1;

        cout << sourceName << '\t'
             << sourceId << '\t'
             << sourceLength << '\t'
             << queryId << '\t'
             << hitBegin - sourceBegin << '\t'
             << hitEnd - sourceBegin << endl;
    }
}

size_t findSourceId(vector<size_t> recordStarts, size_t hitBegin)
{
    for (size_t i = 0; i < recordStarts.size() - 1; ++i)
        if (recordStarts[i] <= hitBegin && hitBegin < recordStarts[i + 1])
            return i;

    throw new FmIndexException("Invalid hit: cannot associate a record.");
}

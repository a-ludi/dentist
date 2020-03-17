#include "fm-index.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/io.hpp>
#include <ctime>
#include <sys/stat.h>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <exception>
#include <poll.h>

using namespace sdsl;
using namespace std;


char* executable;
sdsl::cache_config sdsl_config;

int main(int argc, char** argv)
{
    executable = argv[0];

    try
    {
        int positionalBegin;
        bool includeReverseComplement = false;

        try
        {
            positionalBegin = parseArgs(argc, argv, &(sdsl_config.dir), &includeReverseComplement);

            if (argc - positionalBegin <  1)
                throw FmIndexException("Missing arguments.");
        }
        catch (const FmIndexException& e)
        {
            cerr << "error: " << e.what() << endl;
            cerr << endl;
            printUsage();

            return 1;
        }

        string referenceFile = string(argv[positionalBegin]);
        char** queriesFiles = &argv[positionalBegin + 1];
        int numQueries = max(0, argc - positionalBegin - 1);

        csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fmIndex(sdsl_config);
        buildIndex(fmIndex, referenceFile);
        auto recordStarts = getRecordStarts(fmIndex, referenceFile);

        string queryBuffer;
        // Reserve 10 Mb for a single query
        queryBuffer.reserve(10 * (2 << 20));

        if (numQueries == 0)
            locateQueries(fmIndex, recordStarts, cin, "stdin", queryBuffer, includeReverseComplement);
        else
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

                locateQueries(fmIndex, recordStarts, queryData, queryFile, queryBuffer, includeReverseComplement);
            }
    }
    catch (const FmIndexException& e)
    {
        cerr << "error: " << e.what() << endl;

        return 2;
    }
    catch (const exception& e)
    {
        cerr << "critical error: " << e.what() << endl;

        return 3;
    }
    catch (...)
    {
        cerr << "critical error: Unexpected exception." << endl;

        return 4;
    }

    return 0;
}

string getDefaultTempDir()
{
    char* tempDir;

    #ifdef _GNU_SOURCE
        tempDir = secure_getenv("TMPDIR");
    #else
        tempDir = getenv("TMPDIR");
    #endif

    if (NULL == tempDir)
        return "/tmp";

    return string(tempDir);
}

void printUsage()
{
    cerr << "Usage " << executable << " [-P<dir>] [-r] <in:reference> [<in:queries> ...]" << endl;
    cerr << "    This program constructs a very compact FM-index" << endl;
    cerr << "    of <reference> and locates <queries> if given." << endl;
    cerr << endl;
    cerr << "Positional arguments:" << endl;
    cerr << "    <in:refernce>  Original text file with one record per line." << endl;
    cerr << "    <in:queries>   List of queries (one per line) to locate in <reference>." << endl;
    cerr << "                   Reads standard input if no <queries> given." <<endl;
    cerr << endl;
    cerr << "Optional arguments:" << endl;
    cerr << "    -P<dir>        Use <dir> as temporary directory (default: " << getDefaultTempDir() << ")." << endl;
    cerr << "    -r             Search the reverse complement of each query as well" << endl;
    cerr << endl;
    cerr << "Output:" << endl;
    cerr << "    Creates an FM-index <reference>.fm9 if not present. Then produces" << endl;
    cerr << "    a list of exact matches for all queries in a TAB-separated format:" << endl;
    cerr << endl;
    cerr << "        refName  refId  refLength  queryId  hitBegin  hitEnd  revComp" << endl;
    cerr << endl;
    cerr << "    IDs and coordinates are zero-based. Coordinates are right-open." << endl;
    cerr << endl;
    cerr << "    Note: no output will be produced if no queries are given." << endl;
}

int parseArgs(int argc, char** argv, string *tempDir, bool *includeReverseComplement)
{
    *tempDir = getDefaultTempDir();

    int i;
    for (i = 1; i < argc && argv[i][0] == '-'; ++i)
        switch (argv[i][1])
        {
            case 'P':
                *tempDir = string(argv[i] + 2);

                if (tempDir->size() == 0)
                    throw FmIndexException("Missing value for -P.");
                break;
            case 'r':
                *includeReverseComplement = true;

                if (*(argv[i] + 2) != '\0')
                    throw FmIndexException("Flag -r takes no value.");
                break;
            default:
                throw FmIndexException(string("Invalid option -") + argv[i][1] + ".");
        }

    if (!dirExists(*tempDir))
        throw FmIndexException(string("Cannot open temporary directory: ") + *tempDir);

    return i;
}

void buildIndex(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex, string referenceFile)
{
    string indexSuffix = ".fm9";
    string indexFile = referenceFile + indexSuffix;

    if (!load_from_file(fmIndex, indexFile))
    {
        ifstream referenceData(referenceFile.c_str());

        if (!referenceData)
            throw FmIndexException("File `" + referenceFile + "` does not exist.");

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

vector<size_t> getRecordStarts(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex, string referenceFile)
{
    string recordIndexSuffix = ".idx";
    string recordIndexFile = referenceFile + recordIndexSuffix;
    vector<size_t> recordStarts;

    if (!load_from_plain_array<size_t>(recordStarts, recordIndexFile))
    {
        cerr << "{\"level\":\"info\","
             << "\"info\":\"Record index does not exist. Building it now.\","
             << "\"file\":\"" << recordIndexFile << "\""
             << "}" << endl;

        time_t start = time(NULL);

        auto lineEndLocations = locate(fmIndex, "\n");
        sort(lineEndLocations.begin(), lineEndLocations.end());
        recordStarts.resize(lineEndLocations.size() + 1);

        recordStarts[0] = 0;
        for (int i = 1; i < recordStarts.size(); ++i)
            recordStarts[i] = lineEndLocations[i - 1] + 1;

        if (!store_to_plain_array<size_t>(recordStarts, recordIndexFile))
            throw FmIndexException("Could not store record index: " + recordIndexFile);

        cerr << "{\"level\":\"info\","
             << "\"info\":\"Built record index.\","
             << "\"file\":\"" << recordIndexFile << "\","
             << "\"numRecords\":" << recordStarts.size() - 1 << ","
             << "\"elpasedSecs\":" << difftime(time(NULL), start) << ","
             << "\"sizeMiB\":" << size_in_mega_bytes(recordStarts)
             << "}" << endl;
    }

    return recordStarts;
}

void locateQueries(
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex,
    vector<size_t> recordStarts,
    istream &queriesData,
    string sourceName,
    string queryBuffer,
    bool includeReverseComplement
)
{
    cerr << "{\"level\":\"info\","
         << "\"info\":\"Processing queries.\","
         << "\"source\":\"" << sourceName << "\""
         << "}" << endl;

    time_t start = time(NULL);
    size_t numHits = 0;
    size_t queryId = 0;
    static string revCompBuffer;
    while (queriesData && !queriesData.eof())
    {
        getline(queriesData, queryBuffer);
        if (queryBuffer.size() > 0)
        {
            numHits += locateQuery(fmIndex, recordStarts, sourceName, queryId, queryBuffer, false);

            if (includeReverseComplement)
            {
                getReverseComplement(revCompBuffer, queryBuffer);
                numHits += locateQuery(fmIndex, recordStarts, sourceName, queryId, revCompBuffer, true);
            }

            ++queryId;
        }
    }

    cerr << "{\"level\":\"info\","
         << "\"info\":\"Finished queries.\","
         << "\"source\":\"" << sourceName << "\","
         << "\"numHits\":" << numHits << ","
         << "\"elpasedSecs\":" << difftime(time(NULL), start)
         << "}" << endl;
}

size_t locateQuery(
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fmIndex,
    vector<size_t> recordStarts,
    string sourceName,
    size_t queryId,
    string query,
    bool reverseComplement
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
             << hitEnd - sourceBegin << '\t'
             << (reverseComplement ? "yes" : "no") << endl;
    }

    return locations.size();
}

char revCompTable[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//      A       C             G                      N                   T
    0, 'T', 0, 'G', 0, 0, 0, 'C', 0, 0, 0, 0, 0, 0, 'N', 0, 0, 0, 0, 0, 'A', 0, 0, 0, 0, 0, 0, 0, 0,

//               a       c             g                      n                   t
    0, 0, 0, 0, 't', 0, 'g', 0, 0, 0, 'c', 0, 0, 0, 0, 0, 0, 'n', 0, 0, 0, 0, 0, 'a', 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0
};

void getReverseComplement(string &revCompBuffer, string query)
{
    auto queryLen = query.size();
    revCompBuffer.resize(queryLen);

    for (int i = 0; i < queryLen; ++i)
        revCompBuffer[queryLen - i - 1] = revCompTable[query[i]];
}

size_t findSourceId(vector<size_t> recordStarts, size_t hitBegin)
{
    for (size_t i = 0; i < recordStarts.size() - 1; ++i)
        if (recordStarts[i] <= hitBegin && hitBegin < recordStarts[i + 1])
            return i;

    throw FmIndexException("Invalid hit: cannot associate a record.");
}

bool dirExists(const string& dir)
{
    struct stat dirStat;

    return 0 == stat(dir.c_str(), &dirStat) && ((dirStat.st_mode & S_IFMT) == S_IFDIR);
}

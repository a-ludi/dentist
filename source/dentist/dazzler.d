/**
    Defines bindings and utilities to/for the dazzler commands.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.dazzler;

import core.memory : GC;
import dentist.common : ReferenceInterval, ReferenceRegion;
import dentist.common.alignments :
    AlignmentChain,
    AlignmentFlag = Flag,
    AlignmentFlags = Flags,
    ChainingOptions,
    chainLocalAlignmentsAlgo = chainLocalAlignments,
    Contig,
    coord_t,
    diff_t,
    FlatLocalAlignment,
    id_t,
    Locus,
    trace_point_t,
    TracePoint,
    TranslatedTracePoint;
import dentist.common.binio : CompressedSequence;
import dentist.common.external : ExternalDependency;
import dentist.util.algorithm : sliceUntil;
import dentist.util.fasta : parseFastaRecord, reverseComplement;
import dentist.util.log;
import dentist.util.math : absdiff, floor, ceil, RoundingMode;
import dentist.util.process : executePipe = pipeLines;
import dentist.util.range : arrayChunks, takeExactly;
import dentist.util.region : convexHull, findTilings, min, sup;
import dentist.util.string :
    EditOp,
    findAlignment,
    longestInputsLength,
    memoryRequired,
    score_t,
    SequenceAlignment;
import dentist.util.tempfile : mkstemp;
import std.algorithm :
    all,
    among,
    cache,
    canFind,
    copy,
    countUntil,
    cumulativeFold,
    endsWith,
    filter,
    find,
    isSorted,
    joiner,
    map,
    max,
    maxElement,
    min,
    minElement,
    remove,
    sort,
    splitter,
    startsWith,
    sum,
    SwapStrategy,
    uniq,
    until;
import std.array :
    appender,
    Appender,
    array,
    replace,
    split,
    uninitializedArray;
import std.conv :
    ConvException,
    to;
import std.exception : enforce;
import std.file : exists, remove;
import std.format : format, formattedRead;
import std.math : isNaN;
import std.meta : AliasSeq;
import std.path :
    absolutePath,
    baseName,
    buildPath,
    dirName,
    extension,
    relativePath,
    stripExtension,
    withExtension;
import std.process : Config, escapeShellCommand, kill, pipeProcess,
    ProcessPipes, Redirect, wait;
import std.range :
    chain,
    chunks,
    drop,
    enumerate,
    generate,
    only,
    repeat,
    slide,
    take,
    takeExactly,
    zip;
import std.range.primitives :
    ElementType,
    empty,
    front,
    isForwardRange,
    isInputRange,
    popFront,
    save,
    walkLength;
import std.stdio : File;
import std.string : lineSplitter, outdent;
import std.traits : isArray, isIntegral, isSomeString, ReturnType, Unqual;
import std.typecons : Flag, No, tuple, Tuple, Yes;
import std.uni : toUpper;
import std.variant : Algebraic;
import vibe.data.json : Json, toJson = serializeToJson;

debug import std.stdio : writeln;


/// File suffixes of hidden .db files.
private enum hiddenDbFileSuffixes = [".bps", ".idx"];

/// File suffixes of hidden .dam files.
private enum hiddenDamFileSuffixes = [".bps", ".hdr", ".idx"];

/// Constant holding the .db file extension.
enum dbFileExtension = ".db";

/// Constant holding the .dam file extension.
enum damFileExtension = ".dam";

/// The Dazzler tools require sequence of a least minSequenceLength base pairs.
enum minSequenceLength = 14;

/// This trace point distance enforces the use of ushort for trace point encoding.
enum forceLargeTracePointType = 126;

/// Minimum allowed value for `-e` option of `daligner`/`damapper`
enum minAverageCorrelationRate = 0.7;

/// Minimum allowed value for `-n` option of `damapper`
version (damapper) enum minBestMatches = 0.7;

enum isOptionsList(T) = isArray!T && isSomeString!(ElementType!T);

/**
    Return a list of hidden files associated to every `.dam`/`.db` file. These
    files contain the actual data used in all the computation. Thus, we
    carefully check for their existence.
*/
auto getHiddenDbFiles(string dbFile)
{
    import std.algorithm : map;

    assert(dbFile.endsWith(dbFileExtension, damFileExtension), "must use with Dazzler DB");
    auto suffixes = dbFile.endsWith(dbFileExtension)
        ? hiddenDbFileSuffixes
        : hiddenDamFileSuffixes;

    return suffixes.map!(suffix => buildPath(dbFile.dirName,
            "." ~ dbFile.baseName.withExtension(suffix).to!string));
}


string stripDbExtension(in string dbFile) pure nothrow @safe
{
    if (dbFile.endsWith(dbFileExtension))
        return dbFile[0 .. $ - dbFileExtension.length];
    else if (dbFile.endsWith(damFileExtension))
        return dbFile[0 .. $ - damFileExtension.length];
    else
        return dbFile[];
}


class DazzlerCommandException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

alias Append = Flag!"append";

/// This funcion is called if writing to an existing DB is attempted without
/// passing `Yes.append` or similar to the called function.
void delegate(string dbFile) handleExistingDb = (dbFile) => removeDB(dbFile);


/// Returns true iff lasFile contains zero parts.
bool lasEmpty(in string lasFile)
{
    auto header = readLasHeader(lasFile);

    return header.numParts == 0;
}

deprecated("use the one argument version")
bool lasEmpty(in string lasFile, in string dbA, in string workdir)
{
    return lasEmpty(lasFile);
}

/// ditto
deprecated("use the one argument version")
bool lasEmpty(in string lasFile, in string dbA, in string dbB, in string workdir)
{
    return lasEmpty(lasFile);
}


/// Returns the number of records in dbFile.
id_t numDbRecords(in string dbFile)
{
    auto dbdumpLines = dbdump(dbFile, []);
    scope (exit) dbdumpLines.destroy();
    auto recordNumberLine = dbdumpLines.find!(l => l.startsWith("+ R")).front;
    id_t numRecords;

    recordNumberLine.formattedRead!"+ R %d"(numRecords);
    // Clean up child process

    return numRecords;
}

/// Returns true iff dbFile is empty.
bool dbEmpty(in string dbFile)
{
    return numDbRecords(dbFile) == 0;
}

/**
    Remove database and hidden files.
*/
@ExternalDependency("DBrm", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
void removeDB(in string dbFile)
{
    executeCommand(only("DBrm", dbFile));
}

/// Build outputDb file by using the given subset of reads in inDbFile.
string dbSubset(Options, R)(in string inDbFile, R readIds, in Options options, Append append = No.append)
        if (isSomeString!(typeof(options.tmpdir)) &&
            isOptionsList!(typeof(options.dbsplitOptions)))
{
    enum outDbNameTemplate = "subset-XXXXXX";

    auto outDbTemplate = buildPath(options.tmpdir, outDbNameTemplate);
    auto outDb = mkstemp(outDbTemplate, inDbFile.extension);

    outDb.file.close();
    remove(outDb.name);

    return dbSubset(outDb.name, inDbFile, readIds, options, append);
}

/**
    Build `outputDb` by using the given subset of reads in `inDbFile`. If no
    `outputDb` is given a temporary file with the same extension as `inDbFile`
    will be created.

    Returns: DB file name
*/
string dbSubset(Options, R)(in string outputDb, in string inDbFile, R readIds, in Options options, Append append = No.append)
        if (isSomeString!(typeof(options.tmpdir)) &&
            isOptionsList!(typeof(options.dbsplitOptions)))
{
    auto _outputDb = outputDb.extension == inDbFile.extension
        ? outputDb
        : outputDb ~ inDbFile.extension;

    buildSubsetDb(inDbFile, _outputDb, readIds, append);
    dbsplit(_outputDb, options.dbsplitOptions);

    return outputDb;
}

/// ditto
AlignmentChain[] getLocalAlignments(Options)(in string dbA, in Options options)
        if (isOptionsList!(typeof(options.dalignerOptions)) &&
            isOptionsList!(typeof(options.ladumpOptions)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    if (!lasFileGenerated(dbA, options.tmpdir))
    {
        dalign(dbA, options.dalignerOptions, options.tmpdir);
    }

    return getGeneratedAlignments(dbA, null, options);
}

AlignmentChain[] getLocalAlignments(Options)(in string dbA, in string dbB, in Options options)
        if (isOptionsList!(typeof(options.dalignerOptions)) &&
            isOptionsList!(typeof(options.ladumpOptions)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    if (!lasFileGenerated(dbA, dbB, options.tmpdir))
    {
        dalign(dbA, dbB, options.dalignerOptions, options.tmpdir);
    }

    return getGeneratedAlignments(dbA, dbB, options);
}

void computeLocalAlignments(Options)(in string[] dbList, in Options options)
        if (isOptionsList!(typeof(options.dalignerOptions)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    dalign(dbList, options.dalignerOptions, options.tmpdir);
}

version (damapper)
{
    AlignmentChain[] getMappings(Options)(in string dbA, in string dbB, in Options options)
            if (isOptionsList!(typeof(options.damapperOptions)) &&
                isOptionsList!(typeof(options.ladumpOptions)) &&
                isSomeString!(typeof(options.tmpdir)))
    {
        if (!lasFileGenerated(dbA, dbB, options.tmpdir))
        {
            damapper(dbA, dbB, options.damapperOptions, options.tmpdir);
        }

        return getGeneratedAlignments(dbA, dbB, options);
    }

    void computeMappings(Options)(in string[] dbList, in Options options)
            if (isOptionsList!(typeof(options.damapperOptions)) &&
                isOptionsList!(typeof(options.ladumpOptions)) &&
                isSomeString!(typeof(options.tmpdir)))
    {
        damapper(dbList, options.damapperOptions, options.tmpdir);
    }
}

private AlignmentChain[] getGeneratedAlignments(Options)(
    in string dbA,
    in string dbB,
    in Options options
)
        if (isOptionsList!(typeof(options.ladumpOptions)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    auto lasFile = getLasFile(dbA, dbB, options.tmpdir);

    return getAlignments(dbA, dbB, lasFile, options);
}


enum AlignmentReaderFlag : uint
{
    none = 0,
    includeTracePoints = 1 << 0,
    sort = 1 << 1,
}


AlignmentChain[] getAlignments(
    in string dbA,
    in string lasFile,
    Flag!"includeTracePoints" includeTracePoints = No.includeTracePoints,
)
{
    return getAlignments(dbA, null, lasFile, includeTracePoints);
}

AlignmentChain[] getAlignments(
    in string dbA,
    in string dbB,
    in string lasFile,
    Flag!"includeTracePoints" includeTracePoints = No.includeTracePoints,
)
{
    auto flags = AlignmentReaderFlag.sort;

    if (includeTracePoints)
        flags |= AlignmentReaderFlag.includeTracePoints;

    return getAlignments(dbA, dbB, lasFile, flags);
}

AlignmentChain[] getAlignments(
    in string dbA,
    in string lasFile,
    AlignmentReaderFlag flags,
)
{
    return getAlignments(dbA, null, lasFile, flags);
}

AlignmentChain[] getAlignments(
    in string dbA,
    in string dbB,
    in string lasFile,
    AlignmentReaderFlag flags,
)
{
    string[] ladumpOptions = [
        LAdumpOptions.coordinates,
        LAdumpOptions.numDiffs,
        LAdumpOptions.lengths,
    ];
    trace_point_t tracePointDistance;

    if (flags & AlignmentReaderFlag.includeTracePoints)
    {
        ladumpOptions ~= LAdumpOptions.tracePoints;
        tracePointDistance = getTracePointDistance(lasFile);
    }

    auto lasdumpReader = readLasDump(ladump(
        lasFile,
        dbA,
        dbB,
        ladumpOptions,
    ), tracePointDistance);
    scope (exit) lasdumpReader.closePipe();
    auto alignmentChains = lasdumpReader.array;

    if (flags & AlignmentReaderFlag.sort)
        alignmentChains.sort!("a < b", SwapStrategy.stable);

    return alignmentChains;
}

deprecated("use version without arguments workdir and tracePointDistance 1")
AlignmentChain[] getAlignments(
    in string dbA,
    in string lasFile,
    in string workdir,
    in trace_point_t tracePointDistance = 0,
)
{
    return getAlignments(
        dbA,
        null,
        lasFile,
        cast(Flag!"includeTracePoints") (tracePointDistance > 0),
    );
}

deprecated("use version without arguments workdir and tracePointDistance 2")
AlignmentChain[] getAlignments(
    in string dbA,
    in string dbB,
    in string lasFile,
    in string workdir,
    in trace_point_t tracePointDistance = 0,
)
{
    return getAlignments(
        dbA,
        dbB,
        lasFile,
        cast(Flag!"includeTracePoints") (tracePointDistance > 0),
    );
}


/// Returns a tuple of contig IDs and first and last begin/end coords.
auto fingerprint(in ref AlignmentChain alignmentChain) pure nothrow
{
    return tuple(
        alignmentChain.contigA.id,
        alignmentChain.contigB.id,
        alignmentChain.first.contigA.begin + 0,
        alignmentChain.first.contigB.begin + 0,
        alignmentChain.last.contigA.end + 0,
        alignmentChain.last.contigB.end + 0,
    );
}


private struct LasDumpLineFormatTuple
{
    char indicator;
    char subIndicator;
    string format;
}

private enum LasDumpLineFormat : LasDumpLineFormatTuple
{
    totalChainPartsCount = LasDumpLineFormatTuple('+', 'P', "+ P %d"),
    totalTracePointsCount = LasDumpLineFormatTuple('+', 'T', "+ T %d"),
    maxChainPartsCountPerPile = LasDumpLineFormatTuple('%', 'P', "%% P %d"),
    maxTracePointsCountPerPile = LasDumpLineFormatTuple('%', 'T', "%% T %d"),
    maxTracePointCount = LasDumpLineFormatTuple('@', 'T', "@ T %d"),
    tracePointDistance = LasDumpLineFormatTuple('X', '\0', "X %d"),
    chainPart = LasDumpLineFormatTuple('P', '\0', "P %d %d %c %c"),
    lengths = LasDumpLineFormatTuple('L', '\0', "L %d %d"),
    coordinates = LasDumpLineFormatTuple('C', '\0', "C %d %d %d %d"),
    numDiffs = LasDumpLineFormatTuple('D', '\0', "D %d"),
    tracePointBegin = LasDumpLineFormatTuple('T', '\0', "T %d"),
    tracePoint = LasDumpLineFormatTuple(' ', '\0', " %d %d"),
}

private enum ChainPartType : char
{
    start = '>',
    continuation = '-',
    alternateStart = '+',
    noChainInFile = '.',
}

private struct LasDumpReader(S) if (isInputRange!S && isSomeString!(ElementType!S))
{
    static alias dstring = immutable(dchar)[];
    static alias LasDump = ReturnType!getDumpLines;
    static alias LocalAlignment = AlignmentChain.LocalAlignment;
    static alias TracePoint = LocalAlignment.TracePoint;

private:
    LasDump lasDump;
    bool _empty;
    AlignmentChain currentAC;
    id_t currentACID;
    Appender!(LocalAlignment[]) localAlignmentsAcc;
    Appender!(TracePoint[]) tracePointsAcc;
    dstring currentDumpLine;
    size_t currentDumpLineNumber;
    dchar currentLineType;
    dchar currentLineSubType;
    size_t maxTracePointCount;
    trace_point_t tracePointDistance;
    debug dchar[] allowedLineTypes;

public:
    this(S lasDump, trace_point_t tracePointDistance)
    {
        this.lasDump = getDumpLines(lasDump);
        this.tracePointDistance = tracePointDistance;
        debug with (LasDumpLineFormat)
        {
            this.allowedLineTypes = [
                totalChainPartsCount.indicator,
                totalTracePointsCount.indicator,
                maxChainPartsCountPerPile.indicator,
                maxTracePointsCountPerPile.indicator,
                maxTracePointCount.indicator,
                tracePointDistance.indicator,
                chainPart.indicator,
            ];
        }
        this.popFront();
    }

    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty LasDumpReader");

        if (lasDump.empty)
        {
            return setEmpty();
        }

        readNextAlignmentChain();
    }

    @property bool empty() const pure nothrow
    {
        return _empty;
    }

    @property AlignmentChain front() pure nothrow
    {
        assert(!empty, "Attempting to fetch the front of an empty LasDumpReader");

        return currentAC;
    }

    static if (__traits(hasMember, lasDump, "destroy"))
        alias closePipe = setEmpty;

    void setEmpty() pure nothrow
    {
        static if (__traits(hasMember, lasDump, "destroy"))
            lasDump.destroy();
        _empty = true;
    }

private:

    static auto getDumpLines(S lasDump)
    {
        return lasDump.enumerate(1).filter!"!a[1].empty";
    }

    void readNextAlignmentChain()
    {
        currentAC = AlignmentChain.init;
        localAlignmentsAcc.clear();
        tracePointsAcc.clear();
        peekDumpLine();

        while (true)
        {
            debug _enforce(allowedLineTypes.canFind(currentLineType), format!"forbidden line type `%c` (allowed: `%(%c%)`)"(currentLineType, allowedLineTypes));

            with (LasDumpLineFormat)
            {
                switch (currentLineType)
                {
                case totalChainPartsCount.indicator:
                    static assert(totalTracePointsCount.indicator == totalChainPartsCount.indicator);

                    switch (currentLineSubType)
                    {
                    case totalChainPartsCount.subIndicator:
                        if (readTotalChainPartsCount() == 0)
                            // do not try to read empty dump
                            return setEmpty();
                        break;
                    case totalTracePointsCount.subIndicator:
                        break; // ignore
                    default:
                        error(format!"unknown line sub-type `%c`"(currentLineSubType));
                    }
                    break;
                case maxChainPartsCountPerPile.indicator:
                    static assert(maxTracePointsCountPerPile.indicator == maxChainPartsCountPerPile.indicator);
                    break; // ignore
                case maxTracePointCount.indicator:
                    _enforce(currentLineSubType == maxTracePointCount.subIndicator, "expected `@ T` line");
                    readMaxTracePointCount();
                    debug disallowCurrentLineType();
                    break;
                case tracePointDistance.indicator:
                    readTracePointDistance();
                    debug disallowCurrentLineType();
                    break;
                case chainPart.indicator:
                    if (readChainPart() == No.wasDumpPartConsumed)
                    {
                        debug allowedLineTypes = [chainPart.indicator];

                        return; // chain completed; stay on current dump line
                    }
                    else
                    {
                        debug allowedLineTypes = [
                            chainPart.indicator,
                            lengths.indicator,
                            coordinates.indicator,
                        ];
                        break;
                    }
                case lengths.indicator:
                    readLengths();
                    debug disallowCurrentLineType();
                    break;
                case coordinates.indicator:
                    readCoordinates();
                    debug
                    {
                        disallowCurrentLineType();
                        allowedLineTypes ~= [
                            numDiffs.indicator,
                            tracePointBegin.indicator,
                        ];
                    }
                    break;
                case numDiffs.indicator:
                    readNumDiffs();
                    debug disallowCurrentLineType();
                    break;
                case tracePointBegin.indicator:
                    readTracePointBegin();
                    debug allowedLineTypes = [tracePoint.indicator];
                    break;
                case tracePoint.indicator:
                    readTracePoint();
                    debug allowedLineTypes = [tracePoint.indicator, chainPart.indicator];
                    break;
                default:
                    error(format!"unknown line type `%c`"(currentLineType));
                }
            }

            if (popDumpLine() == Yes.empty)
            {
                finishCurrentAC();

                return; // EOF reached
            }
        }
    }

    void peekDumpLine()
    {
        auto currentLine = lasDump.front;

        currentDumpLineNumber = currentLine[0];
        currentDumpLine = currentLine[1].array;
        currentLineType = currentDumpLine[0];
        currentLineSubType = currentDumpLine.length >= 3 ? currentDumpLine[2] : '\0';
    }

    Flag!"empty" popDumpLine()
    {
        if (lasDump.empty)
        {
            return Yes.empty;
        }

        lasDump.popFront();
        ++currentDumpLineNumber;

        if (lasDump.empty)
        {
            return Yes.empty;
        }
        else
        {
            peekDumpLine();

            return No.empty;
        }
    }

    id_t readTotalChainPartsCount()
    {
        enum totalChainPartsCountFormat = LasDumpLineFormat.totalChainPartsCount.format;

        id_t totalChainPartsCount;
        currentDumpLine[].formattedRead!totalChainPartsCountFormat(totalChainPartsCount);

        return totalChainPartsCount;
    }

    void readMaxTracePointCount()
    {
        enum maxTracePointCountFormat = LasDumpLineFormat.maxTracePointCount.format;

        currentDumpLine[].formattedRead!maxTracePointCountFormat(maxTracePointCount);

        tracePointsAcc.reserve(maxTracePointCount);
    }

    void readTracePointDistance()
    {
        enum tracePointDistanceFormat = LasDumpLineFormat.tracePointDistance.format;

        currentDumpLine[].formattedRead!tracePointDistanceFormat(tracePointDistance);
    }

    Flag!"wasDumpPartConsumed" readChainPart()
    {
        enum chainPartFormat = LasDumpLineFormat.chainPart.format;
        id_t contigAID;
        id_t contigBID;
        char rawComplement;
        char rawChainPartType;

        currentDumpLine[].formattedRead!chainPartFormat(
            contigAID,
            contigBID,
            rawComplement,
            rawChainPartType,
        );

        auto chainPartType = rawChainPartType.to!ChainPartType;
        AlignmentFlags flags;
        flags.complement = (rawComplement == 'c');
        flags.alternateChain = (chainPartType == ChainPartType.alternateStart);

        auto startingNewChain = currentAC == AlignmentChain.init;

        if (startingNewChain)
        {
            currentAC.contigA.id = contigAID;
            currentAC.contigB.id = contigBID;
            currentAC.flags = flags;
            currentAC.tracePointDistance = tracePointDistance;

            return Yes.wasDumpPartConsumed;
        }
        else if (isChainContinuation(contigAID, contigBID, chainPartType))
        {
            _enforce(
                currentAC.flags.complement == flags.complement,
                "matching both strands in one alignment chain",
            );
            finishCurrentLA();

            return Yes.wasDumpPartConsumed;
        }
        else
        {
            finishCurrentAC();

            return No.wasDumpPartConsumed;
        }

    }

    bool isChainContinuation(in size_t contigAID, in size_t contigBID, in ChainPartType chainPartType)
    {
        return currentAC.contigA.id == contigAID &&
               currentAC.contigB.id == contigBID &&
               chainPartType == ChainPartType.continuation;
    }

    void finishCurrentLA()
    {
        assert(localAlignmentsAcc.data.length > 0);

        localAlignmentsAcc.data[$ - 1].tracePoints = tracePointsAcc.data.dup;
    }

    void finishCurrentLAs()
    {
        if (localAlignmentsAcc.data.length > 0)
        {
            finishCurrentLA();
        }

        currentAC.localAlignments = localAlignmentsAcc.data.dup;
    }

    void finishCurrentAC()
    {
        finishCurrentLAs();
        currentAC.id = currentACID++;

        // ensure alignment chain is valid by triggering invariant of AlignmentChain
        cast(void) currentAC.first;
    }

    void readLengths()
    {
        enum lengthsFormat = LasDumpLineFormat.lengths.format;

        currentDumpLine[].formattedRead!lengthsFormat(
            currentAC.contigA.length,
            currentAC.contigB.length,
        );
    }
    void readCoordinates()
    {
        enum coordinatesFormat = LasDumpLineFormat.coordinates.format;
        LocalAlignment currentLA;

        currentDumpLine[].formattedRead!coordinatesFormat(
            currentLA.contigA.begin,
            currentLA.contigA.end,
            currentLA.contigB.begin,
            currentLA.contigB.end,
        );

        localAlignmentsAcc ~= currentLA;
    }

    void readNumDiffs()
    {
        enum numDiffsFormat = LasDumpLineFormat.numDiffs.format;

        currentDumpLine[].formattedRead!numDiffsFormat(
            localAlignmentsAcc.data[$ - 1].numDiffs,
        );
    }

    void readTracePointBegin()
    {
        enum tracePointBeginFormat = LasDumpLineFormat.tracePointBegin.format;
        size_t numTracePoints;

        currentDumpLine[].formattedRead!tracePointBeginFormat(numTracePoints);

        tracePointsAcc.clear();
        tracePointsAcc.reserve(numTracePoints);
    }

    void readTracePoint()
    {
        enum tracePointFormat = LasDumpLineFormat.tracePoint.format;
        TracePoint currentTP;

        currentDumpLine[].formattedRead!tracePointFormat(
            currentTP.numDiffs,
            currentTP.numBasePairs,
        );

        tracePointsAcc ~= currentTP;
    }

    debug void disallowCurrentLineType() {
        allowedLineTypes = allowedLineTypes
            .filter!(type => type != currentLineType)
            .array;
    }

    void error(in string reason)
    {
        _enforce(false, reason);
    }

    void _enforce(bool condition, lazy string reason)
    {
        enforce!DazzlerCommandException(condition, format!"ill-formatted LAdump output: %s (line %d)"(reason, currentDumpLineNumber));
    }
}

private auto readLasDump(S)(S lasDump, trace_point_t tracePointDistance)
{
    return LasDumpReader!S(lasDump, tracePointDistance);
}

unittest
{
    import std.algorithm : equal;

    enum testLasDump = [
        "+ P 9",
        "% P 9",
        "+ T 12",
        "% T 12",
        "@ T 4",
        "P 1 2 n >",
        "L 13 15",
        "C 3 4 5 6",
        "D 7",
        "P 1 2 n -",
        "L 13 15",
        "C 12 13 14 15",
        "D 16",
        "P 19 20 c +",
        "L 31 33",
        "C 21 22 23 24",
        "D 25",
        "P 19 20 c -",
        "L 31 33",
        "C 30 31 32 33",
        "D 34",
        "T 1",
        "   0 1",
        "P 37 38 n .",
        "L 40 42",
        "C 39 40 41 42",
        "D 43",
        "P 46 47 c .",
        "L 58 60",
        "C 48 49 50 51",
        "D 52",
        "T 3",
        "   2 102",
        "   3 101",
        "   4 104",
        "P 46 47 n .",
        "L 58 60",
        "C 57 58 59 60",
        "D 61",
        "T 3",
        "   3 101",
        "   4 104",
        "   2 102",
        "P 64 65 c >",
        "L 71 72",
        "C 66 67 68 69",
        "D 70",
        "T 4",
        "   6 105",
        "   1 101",
        "   2 100",
        "   3  97",
        "P 55 56 c -",
        "L 80 81",
        "C 75 76 77 78",
        "D 79",
        "T 2",
        "   0   2",
        "   2 102",
        "P 1 3197 c >",
        "C 0 71 12 86",
        "T 1",
        "   3  74",
        "P 1 3197 c -",
        "C 0 8300 0 318",
        "T 3",
        "   6 105",
        "   9 108",
        "   7 105",
    ];

    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias complement = AlignmentFlag.complement;
    alias alternateChain = AlignmentFlag.alternateChain;

    auto alignmentChains = readLasDump(testLasDump, 0).array;
    auto expectedResult = [
        AlignmentChain(
            0,
            Contig(1, 13),
            Contig(2, 15),
            AlignmentFlags(),
            [
                LocalAlignment(
                    Locus(3, 4),
                    Locus(5, 6),
                    7,
                ),
                LocalAlignment(
                    Locus(12, 13),
                    Locus(14, 15),
                    16,
                ),
            ],
        ),
        AlignmentChain(
            1,
            Contig(19, 31),
            Contig(20, 33),
            AlignmentFlags(complement, alternateChain),
            [
                LocalAlignment(
                    Locus(21, 22),
                    Locus(23, 24),
                    25,
                ),
                LocalAlignment(
                    Locus(30, 31),
                    Locus(32, 33),
                    34,
                    [
                        TracePoint(0, 1),
                    ]
                ),
            ],
        ),
        AlignmentChain(
            2,
            Contig(37, 40),
            Contig(38, 42),
            AlignmentFlags(),
            [
                LocalAlignment(
                    Locus(39, 40),
                    Locus(41, 42),
                    43,
                ),
            ],
        ),
        AlignmentChain(
            3,
            Contig(46, 58),
            Contig(47, 60),
            AlignmentFlags(complement),
            [
                LocalAlignment(
                    Locus(48, 49),
                    Locus(50, 51),
                    52,
                    [
                        TracePoint(2, 102),
                        TracePoint(3, 101),
                        TracePoint(4, 104),
                    ],
                ),
            ],
        ),
        AlignmentChain(
            4,
            Contig(46, 58),
            Contig(47, 60),
            AlignmentFlags(),
            [
                LocalAlignment(
                    Locus(57, 58),
                    Locus(59, 60),
                    61,
                    [
                        TracePoint(3, 101),
                        TracePoint(4, 104),
                        TracePoint(2, 102),
                    ],
                ),
            ],
        ),
        AlignmentChain(
            5,
            Contig(64, 71),
            Contig(65, 72),
            AlignmentFlags(complement),
            [
                LocalAlignment(
                    Locus(66, 67),
                    Locus(68, 69),
                    70,
                    [
                        TracePoint(6, 105),
                        TracePoint(1, 101),
                        TracePoint(2, 100),
                        TracePoint(3, 97),
                    ],
                ),
            ],
        ),
        AlignmentChain(
            6,
            Contig(55, 80),
            Contig(56, 81),
            AlignmentFlags(complement),
            [
                LocalAlignment(
                    Locus(75, 76),
                    Locus(77, 78),
                    79,
                    [
                        TracePoint(0, 2),
                        TracePoint(2, 102),
                    ],
                ),
            ],
        ),
        AlignmentChain(
            7,
            Contig(1, 0),
            Contig(3197, 0),
            AlignmentFlags(complement),
            [
                LocalAlignment(
                    Locus(0, 71),
                    Locus(12, 86),
                    0,
                    [
                        TracePoint(3, 74),
                    ],
                ),
                LocalAlignment(
                    Locus(0, 8300),
                    Locus(0, 318),
                    0,
                    [
                        TracePoint(6, 105),
                        TracePoint(9, 108),
                        TracePoint(7, 105),
                    ],
                ),
            ],
        ),
    ];

    assert(alignmentChains == expectedResult);
}


auto getFlatLocalAlignments(
    in string dbA,
    in string lasFile,
    Flag!"includeTracePoints" includeTracePoints = No.includeTracePoints,
)
{
    return getFlatLocalAlignments(dbA, null, lasFile, includeTracePoints);
}

auto getFlatLocalAlignments(
    in string dbA,
    in string dbB,
    in string lasFile,
    Flag!"includeTracePoints" includeTracePoints = No.includeTracePoints,
)
{
    string[] ladumpOptions = [
        LAdumpOptions.coordinates,
        LAdumpOptions.lengths,
    ];
    trace_point_t tracePointDistance;

    if (includeTracePoints)
    {
        ladumpOptions ~= LAdumpOptions.tracePoints;
        tracePointDistance = getTracePointDistance(lasFile);
    }

    return readFlatLasDump(ladump(
        lasFile,
        dbA,
        dbB,
        ladumpOptions,
    ), tracePointDistance);
}


private class FlatLasDumpReader(S) if (isInputRange!S && isSomeString!(ElementType!S))
{
    static alias dstring = immutable(dchar)[];
    static alias LasDump = ReturnType!getDumpLines;
    static alias FlatLocus = FlatLocalAlignment.FlatLocus;

private:

    LasDump lasDump;
    bool _empty;
    FlatLocalAlignment currentLA;
    id_t currentLAId;
    dstring currentDumpLine;
    size_t currentDumpLineNumber;
    dchar currentLineType;
    dchar currentLineSubType;
    size_t numLocalAlignments;
    public size_t maxLocalAlignmentsPerContig;
    public size_t maxTracePointCount;
    size_t numTracePointsLeft;
    TracePoint[] tracePointsAcc;
    public trace_point_t tracePointDistance;

public:

    this(S lasDump, trace_point_t tracePointDistance)
    {
        this.lasDump = getDumpLines(lasDump);
        this.tracePointDistance = tracePointDistance;
        this.popFront();
    }

    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty FlatLasDumpReader");

        if (lasDump.empty)
        {
            return setEmpty();
        }

        readNextFlatLocalAlignment();
    }

    @property bool empty() const pure nothrow
    {
        return _empty;
    }

    @property FlatLocalAlignment front() pure nothrow
    {
        assert(!empty, "Attempting to fetch the front of an empty FlatLasDumpReader");

        return currentLA;
    }

    @property size_t length() const pure nothrow @safe
    {
        return numLocalAlignments;
    }

    static if (__traits(hasMember, lasDump, "destroy"))
        alias closePipe = setEmpty;

    void setEmpty() pure nothrow
    {
        static if (__traits(hasMember, lasDump, "destroy"))
            lasDump.destroy();
        _empty = true;
    }

private:

    static auto getDumpLines(S lasDump)
    {
        return lasDump.enumerate(1).filter!"!a[1].empty";
    }

    void readNextFlatLocalAlignment()
    {
        currentLA = FlatLocalAlignment.init;
        tracePointsAcc.length = 0;
        peekDumpLine();

        while (true)
        {
            with (LasDumpLineFormat)
            {
                switch (currentLineType)
                {
                case totalChainPartsCount.indicator:
                    static assert(totalTracePointsCount.indicator == totalChainPartsCount.indicator);

                    switch (currentLineSubType)
                    {
                    case totalChainPartsCount.subIndicator:
                        numLocalAlignments = readTotalChainPartsCount();

                        if (numLocalAlignments == 0)
                            // do not try to read empty dump
                            return setEmpty();
                        ++numLocalAlignments;
                        break;
                    case totalTracePointsCount.subIndicator:
                        break; // ignore
                    default:
                        error(format!"unknown line sub-type `%c`"(currentLineSubType));
                    }
                    break;
                case maxChainPartsCountPerPile.indicator:
                    switch (currentLineSubType)
                    {
                        case maxChainPartsCountPerPile.subIndicator:
                            readMaxLocalAlignmentsPerContig();
                            break;
                        default:
                            break;
                    }
                    break; // ignore
                case maxTracePointCount.indicator:
                    switch (currentLineSubType)
                    {
                        case maxTracePointCount.subIndicator:
                            readMaxTracePointCount();
                            break;
                        default:
                            break;
                    }
                    break;
                case tracePointDistance.indicator:
                    readTracePointDistance();
                    break;
                case chainPart.indicator:
                    _enforce(
                        numTracePointsLeft == 0,
                        "unexpected end of trace points: next local alignment started before the " ~
                        "expected number of trace points was read",
                    );
                    _enforce(
                        numLocalAlignments > 0,
                        "more local alignments than expected",
                    );
                    if (currentLA != FlatLocalAlignment.init)
                        return finishCurrentLA(); // currentLA finished; stay on current line
                    --numLocalAlignments;
                    readChainPart();
                    break;
                case lengths.indicator:
                    readLengths();
                    break;
                case coordinates.indicator:
                    readCoordinates();
                    break;
                case tracePointBegin.indicator:
                    readTracePointBegin();
                    break;
                case tracePoint.indicator:
                    readTracePoint();
                    break;
                default:
                    error(format!"unknown line type `%c`"(currentLineType));
                }
            }

            if (popDumpLine() == Yes.empty)
            {
                _enforce(
                    numTracePointsLeft == 0,
                    "unexpected end of dump: file ended before the " ~
                    "expected number of trace points was read",
                );

                return finishCurrentLA(); // EOF reached
            }
        }
    }

    void peekDumpLine()
    {
        auto currentLine = lasDump.front;

        currentDumpLineNumber = currentLine[0];
        currentDumpLine = currentLine[1].array;
        currentLineType = currentDumpLine[0];
        currentLineSubType = currentDumpLine.length >= 3 ? currentDumpLine[2] : '\0';
    }

    Flag!"empty" popDumpLine()
    {
        if (lasDump.empty)
        {
            return Yes.empty;
        }

        lasDump.popFront();
        ++currentDumpLineNumber;

        if (lasDump.empty)
        {
            return Yes.empty;
        }
        else
        {
            peekDumpLine();

            return No.empty;
        }
    }

    id_t readTotalChainPartsCount()
    {
        enum totalChainPartsCountFormat = LasDumpLineFormat.totalChainPartsCount.format;

        id_t totalChainPartsCount;
        currentDumpLine[].formattedRead!totalChainPartsCountFormat(totalChainPartsCount);

        return totalChainPartsCount;
    }

    void readMaxLocalAlignmentsPerContig()
    {
        enum maxChainPartsCountPerPile = LasDumpLineFormat.maxChainPartsCountPerPile.format;

        currentDumpLine[].formattedRead!maxChainPartsCountPerPile(maxLocalAlignmentsPerContig);
    }

    void readMaxTracePointCount()
    {
        enum maxTracePointCountFormat = LasDumpLineFormat.maxTracePointCount.format;

        currentDumpLine[].formattedRead!maxTracePointCountFormat(maxTracePointCount);

        tracePointsAcc.reserve(maxTracePointCount);
    }

    void readTracePointDistance()
    {
        enum tracePointDistanceFormat = LasDumpLineFormat.tracePointDistance.format;

        currentDumpLine[].formattedRead!tracePointDistanceFormat(tracePointDistance);
    }

    void readChainPart()
    {
        enum chainPartFormat = LasDumpLineFormat.chainPart.format;
        enum yesComplement = AlignmentFlags(AlignmentFlag.complement);
        enum noComplement = AlignmentFlags();
        id_t contigAID;
        id_t contigBID;
        char rawComplement;
        char rawChainPartType;

        currentDumpLine[].formattedRead!chainPartFormat(
            contigAID,
            contigBID,
            rawComplement,
            rawChainPartType,
        );

        auto flags = rawComplement == 'c' ? yesComplement : noComplement;
        auto chainPartType = rawChainPartType.to!ChainPartType;

        currentLA.id = currentLAId++;
        currentLA.contigA.id = contigAID;
        currentLA.contigB.id = contigBID;
        currentLA.flags = flags;
        currentLA.tracePointDistance = tracePointDistance;
    }

    void readLengths()
    {
        enum lengthsFormat = LasDumpLineFormat.lengths.format;

        currentDumpLine[].formattedRead!lengthsFormat(
            currentLA.contigA.length,
            currentLA.contigB.length,
        );
    }

    void readCoordinates()
    {
        enum coordinatesFormat = LasDumpLineFormat.coordinates.format;

        currentDumpLine[].formattedRead!coordinatesFormat(
            currentLA.contigA.begin,
            currentLA.contigA.end,
            currentLA.contigB.begin,
            currentLA.contigB.end,
        );
    }

    void readTracePointBegin()
    {
        enum tracePointBeginFormat = LasDumpLineFormat.tracePointBegin.format;

        currentDumpLine[].formattedRead!tracePointBeginFormat(numTracePointsLeft);

        tracePointsAcc.length = 0;
    }

    void readTracePoint()
    {
        enum tracePointFormat = LasDumpLineFormat.tracePoint.format;
        TracePoint currentTP;

        currentDumpLine[].formattedRead!tracePointFormat(
            currentTP.numDiffs,
            currentTP.numBasePairs,
        );

        tracePointsAcc ~= currentTP;
        --numTracePointsLeft;
    }

    void finishCurrentLA()
    {
        currentLA.tracePoints = tracePointsAcc;
    }

    void error(in string reason)
    {
        _enforce(false, reason);
    }

    void _enforce(bool condition, lazy string reason)
    {
        enforce!DazzlerCommandException(condition, format!"ill-formatted LAdump output: %s (line %d)"(reason, currentDumpLineNumber));
    }
}

private auto readFlatLasDump(S)(S lasDump, trace_point_t tracePointDistance)
{
    return new FlatLasDumpReader!S(lasDump, tracePointDistance);
}

unittest
{
    import std.algorithm : equal;
    import std.range : tee;

    enum expectedTracePointDistance = 100;
    enum testLasDump = [
        "+ P 11",
        "% P 4",
        "+ T 17",
        "% T 4",
        "@ T 4",
        "X 100",
        "P 1 2 n >",
        "L 13 15",
        "C 3 4 5 6",
        "P 1 2 n -",
        "L 13 15",
        "C 12 13 14 15",
        "P 19 20 c +",
        "L 31 33",
        "C 21 22 23 24",
        "P 19 20 c -",
        "L 31 33",
        "C 30 31 32 33",
        "T 1",
        "   0 1",
        "P 37 38 n .",
        "L 40 42",
        "C 39 40 41 42",
        "P 46 47 c .",
        "L 58 60",
        "C 48 49 50 51",
        "T 3",
        "   2 102",
        "   3 101",
        "   4 104",
        "P 46 47 n .",
        "L 58 60",
        "C 57 58 59 60",
        "T 3",
        "   3 101",
        "   4 104",
        "   2 102",
        "P 64 65 c >",
        "L 71 72",
        "C 66 67 68 69",
        "T 4",
        "   6 105",
        "   1 101",
        "   2 100",
        "   3  97",
        "P 55 56 c -",
        "L 80 81",
        "C 75 76 77 78",
        "T 2",
        "   0   2",
        "   2 102",
        "P 1 3197 c >",
        "C 0 71 12 86",
        "T 1",
        "   3  74",
        "P 1 3197 c -",
        "C 0 8300 0 318",
        "T 3",
        "   6 105",
        "   9 108",
        "   7 105",
    ];

    auto flatLocalAlignments = readFlatLasDump(testLasDump, 0)
        .tee!((ref la) { la.tracePoints = la.tracePoints.dup; })
        .array;
    alias FlatLocus = FlatLocalAlignment.FlatLocus;
    auto expectedResult = [
        FlatLocalAlignment(
            0,
            FlatLocus(1, 13, 3, 4),
            FlatLocus(2, 15, 5, 6),
            AlignmentFlags(),
            expectedTracePointDistance,
        ),
        FlatLocalAlignment(
            1,
            FlatLocus(1, 13, 12, 13),
            FlatLocus(2, 15, 14, 15),
            AlignmentFlags(),
            expectedTracePointDistance,
        ),
        FlatLocalAlignment(
            2,
            FlatLocus(19, 31, 21, 22),
            FlatLocus(20, 33, 23, 24),
            AlignmentFlags(AlignmentFlag.complement),
            expectedTracePointDistance,
        ),
        FlatLocalAlignment(
            3,
            FlatLocus(19, 31, 30, 31),
            FlatLocus(20, 33, 32, 33),
            AlignmentFlags(AlignmentFlag.complement),
            expectedTracePointDistance,
            [
                TracePoint(0, 1),
            ],
        ),
        FlatLocalAlignment(
            4,
            FlatLocus(37, 40, 39, 40),
            FlatLocus(38, 42, 41, 42),
            AlignmentFlags(),
            expectedTracePointDistance,
        ),
        FlatLocalAlignment(
            5,
            FlatLocus(46, 58, 48, 49),
            FlatLocus(47, 60, 50, 51),
            AlignmentFlags(AlignmentFlag.complement),
            expectedTracePointDistance,
            [
                TracePoint(2, 102),
                TracePoint(3, 101),
                TracePoint(4, 104),
            ],
        ),
        FlatLocalAlignment(
            6,
            FlatLocus(46, 58, 57, 58),
            FlatLocus(47, 60, 59, 60),
            AlignmentFlags(),
            expectedTracePointDistance,
            [
                TracePoint(3, 101),
                TracePoint(4, 104),
                TracePoint(2, 102),
            ],
        ),
        FlatLocalAlignment(
            7,
            FlatLocus(64, 71, 66, 67),
            FlatLocus(65, 72, 68, 69),
            AlignmentFlags(AlignmentFlag.complement),
            expectedTracePointDistance,
            [
                TracePoint(6, 105),
                TracePoint(1, 101),
                TracePoint(2, 100),
                TracePoint(3, 97),
            ],
        ),
        FlatLocalAlignment(
            8,
            FlatLocus(55, 80, 75, 76),
            FlatLocus(56, 81, 77, 78),
            AlignmentFlags(AlignmentFlag.complement),
            expectedTracePointDistance,
            [
                TracePoint(0, 2),
                TracePoint(2, 102),
            ],
        ),
        FlatLocalAlignment(
            9,
            FlatLocus(1, 0, 0, 71),
            FlatLocus(3197, 0, 12, 86),
            AlignmentFlags(AlignmentFlag.complement),
            expectedTracePointDistance,
            [
                TracePoint(3, 74),
            ],
        ),
        FlatLocalAlignment(
            10,
            FlatLocus(1, 0, 0, 8300),
            FlatLocus(3197, 0, 0, 318),
            AlignmentFlags(AlignmentFlag.complement),
            expectedTracePointDistance,
            [
                TracePoint(6, 105),
                TracePoint(9, 108),
                TracePoint(7, 105),
            ],
        ),
    ];

    assert(flatLocalAlignments == expectedResult);
}


struct AlignmentHeader
{
    size_t numAlignments;
    size_t numLocalAlignments;
    size_t maxLocalAlignments;
    size_t numTracePoints;
    size_t maxTracePoints;
    size_t tracePointDistance;

    static AlignmentHeader inferFrom(R)(R alignmentChains) if (isInputRange!R)
    {
        AlignmentHeader headerData;

        if (alignmentChains.empty)
            headerData.tracePointDistance = 100;
        else
            headerData.tracePointDistance = alignmentChains.front.tracePointDistance;

        foreach (alignmentChain; alignmentChains)
        {
            ++headerData.numAlignments;
            headerData.numLocalAlignments += alignmentChain.localAlignments.length;
            headerData.maxLocalAlignments = max(
                headerData.maxLocalAlignments,
                alignmentChain.localAlignments.length,
            );

            foreach (localAlignment; alignmentChain.localAlignments)
            {
                headerData.numTracePoints += localAlignment.tracePoints.length;
                headerData.maxTracePoints = max(
                    headerData.maxTracePoints,
                    localAlignment.tracePoints.length,
                );
            }
        }

        return headerData;
    }
}

AlignmentHeader writeAlignments(R)(const string lasFile, R alignmentChains)
    if (isForwardRange!R)
{
    auto headerData = AlignmentHeader.inferFrom(alignmentChains.save);

    return lasFile.writeAlignments(alignmentChains, headerData);
}

@ExternalDependency("dumpLA", null, "https://github.com/thegenemyers/DALIGNER")
AlignmentHeader writeAlignments(R)(
    const string lasFile,
    R alignmentChains,
    AlignmentHeader headerData,
) if (isInputRange!R)
{
    auto dumpCommand = ["dumpLA", lasFile];
    auto converter = pipeProcess(dumpCommand, Redirect.stdin);
    auto writer = converter.stdin;

    logJsonDiagnostic(
        "action", "execute",
        "type", "pipe",
        "command", dumpCommand.toJson,
        "state", "pre",
    );

    writer.writefln!"@ T %d"(2 * headerData.maxTracePoints);
    writer.writefln!"X %d"(headerData.tracePointDistance);

    foreach (alignmentChain; alignmentChains)
        writer.dumpAlignment(alignmentChain);

    writer.close();
    wait(converter.pid());

    return headerData;
}

private auto dumpAlignment(File writer, const AlignmentChain alignmentChain)
{
    if (alignmentChain.flags.disabled || alignmentChain.localAlignments.length == 0)
        return;


    foreach (i, localAlignment; alignmentChain.localAlignments)
    {
        auto chainPartType = i == 0
            ? (alignmentChain.flags.alternateChain
                ? ChainPartType.alternateStart
                : ChainPartType.start)
            : ChainPartType.continuation;

        writer.writefln!"P %d %d %c %c"(
            alignmentChain.contigA.id,
            alignmentChain.contigB.id,
            alignmentChain.flags.complement ? 'c' : 'n',
            chainPartType,
        );
        writer.writefln!"C %d %d %d %d"(
            localAlignment.contigA.begin,
            localAlignment.contigA.end,
            localAlignment.contigB.begin,
            localAlignment.contigB.end,
        );

        size_t totalDiffs;
        writer.writefln!"T %d"(localAlignment.tracePoints.length);
        foreach (tracePoint; localAlignment.tracePoints)
        {
            totalDiffs += tracePoint.numDiffs;
            writer.writefln!"  %d %d"(tracePoint.numDiffs, tracePoint.numBasePairs);
        }

        writer.writefln!"D %d"(totalDiffs);
    }
}


/// Returns the trace point distance in lasFile.
trace_point_t getTracePointDistance(in string lasFile)
{
    return readLasHeader(lasFile).tracePointSpacing;
}

/// ditto
trace_point_t getTracePointDistance(in string dbA, in string lasFile)
{
    return getTracePointDistance(lasFile);
}

/// ditto
trace_point_t getTracePointDistance(in string dbA, in string dbB, in string lasFile)
{
    return getTracePointDistance(lasFile);
}


auto getExactAlignment(
    in string dbA,
    in string dbB,
    in AlignmentChain ac,
    in size_t memoryLimit = 2^^20,
)
{
    return getExactAlignment(
        dbA,
        dbB,
        ac,
        ac.first.contigA.begin,
        ac.last.contigA.end,
        memoryLimit,
    );
}

auto getExactAlignment(
    in string dbA,
    in string dbB,
    in AlignmentChain ac,
    in coord_t beginA,
    in coord_t endA,
    in size_t memoryLimit = 2^^20,
)
{
    assert(ac.tracePointDistance > 0, "trace points required for getExactAlignment");
    assert(
        ac.first.contigA.begin <= beginA &&
        beginA < endA &&
        endA <= ac.last.contigA.end
    );

    // Translate input coords to tracepoints to get exact alignments
    auto begin = ac.translateTracePoint(beginA, RoundingMode.floor);
    auto end = ac.translateTracePoint(endA, RoundingMode.ceil);

    // Fetch relevant sequences from DBs
    auto aSequence = getFastaSequence(dbA, ac.contigA.id);
    auto bSequence = getFastaSequence(dbB, ac.contigB.id);

    // Slice sequences to translated coordinates
    aSequence = aSequence[begin.contigA .. end.contigA];
    if (ac.flags.complement)
        bSequence = reverseComplement(bSequence[$ - end.contigB .. $ - begin.contigB]);
    else
        bSequence = bSequence[begin.contigB .. end.contigB];

    auto paddedAlignment = getPaddedAlignment(
        ac,
        begin,
        end,
        aSequence,
        bSequence,
        memoryLimit,
    );

    assert(begin.contigA <= beginA);

    return paddedAlignment[(beginA - begin.contigA) .. min(endA - begin.contigA, $)];
}

auto getPaddedAlignment(S, TranslatedTracePoint)(
    in AlignmentChain ac,
    in TranslatedTracePoint begin,
    in TranslatedTracePoint end,
    S aSequence,
    S bSequence,
    in size_t memoryLimit = 2^^20,
)
{
    static struct AlignmentPadder
    {
        static enum indelPenalty = 1;

        alias LocalAlignment = AlignmentChain.LocalAlignment;

        private const AlignmentChain ac;
        private const TranslatedTracePoint begin;
        private const TranslatedTracePoint end;
        private S aSequence;
        private S bSequence;
        private const size_t memoryLimit;
        private SequenceAlignment!S _paddedAlignment;

        this(
            in AlignmentChain ac,
            in TranslatedTracePoint begin,
            in TranslatedTracePoint end,
            S aSequence,
            S bSequence,
            in size_t memoryLimit = 2^^20,
        )
        {
            this.ac = ac;
            this.begin = begin;
            this.end = end;
            this.aSequence = aSequence;
            this.bSequence = bSequence;
            this.memoryLimit = memoryLimit;
            this._paddedAlignment = typeof(this._paddedAlignment)(
                0,
                [],
                aSequence[0 .. 0],
                bSequence[0 .. 0],
                indelPenalty,
                No.freeShift,
            );
            this._paddedAlignment.editPath.reserve(aSequence.length + bSequence.length);
        }

        @property auto paddedAlignment()
        {
            if (_paddedAlignment.editPath.length == 0)
                computePaddedAlignment();

            return _paddedAlignment;
        }

    private:

        coord_t aSeqPos;
        coord_t bSeqPos;

        void computePaddedAlignment()
        {
            aSeqPos = begin.contigA + 0;
            bSeqPos = begin.contigB + 0;
            auto cleanedLocalAlignments = cleanUpLocalAlignments(ac.localAlignments);
            auto coveringLocalAlignments = cleanedLocalAlignments
                .find!(la => la.contigA.begin <= begin.contigA)
                .sliceUntil!(la => end.contigA < la.contigA.begin);
            assert(!coveringLocalAlignments.empty);

            foreach (i, localAlignment; coveringLocalAlignments.enumerate)
            {
                if (i > 0)
                {
                    if (
                        aSeqPos <= localAlignment.contigA.begin &&
                        bSeqPos <= localAlignment.contigB.begin
                    )
                        addGapAlignment(localAlignment);
                    else
                        resolveOverlappingAlignments(
                            coveringLocalAlignments[i - 1],
                            localAlignment,
                        );
                }

                auto skip = skipTracePointsToASeqPos(localAlignment);
                foreach (tracePoint; localAlignment.tracePoints[skip .. $])
                {
                    if (aSeqPos >= end.contigA)
                        return;

                    addTracePointAlignment(localAlignment, tracePoint);
                }
            }
        }

        // NOTE: damapper may produce false/bad chains that include a
        //       complete stack of local alignments that have nearly
        //       100% overlap. In most cases they decrease performance
        //       drastically and in some cases they can even cause
        //       errors in this procedure. So let's ignore them at
        //       this point.
        const(LocalAlignment[]) cleanUpLocalAlignments(in LocalAlignment[] localAlignments)
        {
            ReferenceInterval toInterval(in LocalAlignment* la) pure
            {
                return ReferenceInterval(
                    0,
                    la.contigA.begin,
                    la.contigA.end,
                );
            }

            assert(localAlignments.length > 0, "localAlignments must not be empty");
            auto coveredInterval = ReferenceInterval(
                0,
                localAlignments[0].contigA.begin,
                localAlignments[$ - 1].contigA.end,
            );
            auto combinations = findTilings!toInterval(
                localAlignments.map!((ref la) => &la).array,
                longestInputsLength(memoryLimit),
            );
            assert(combinations.length > 0, "no viable tilings found");

            alias Combination = typeof(combinations[0]);
            auto combinationScore(in Combination combination)
            {
                long coveredBasePairs = combination.region.size;
                long score =
                     + combination.region.size
                     - combination.totalOverlap
                     - combination.elements.map!"a.numDiffs".sum;

                return score;
            }

            return combinations
                .filter!(combination => convexHull(combination.region) == coveredInterval)
                .maxElement!combinationScore
                .elements
                .map!"*a".array;
        }

        size_t skipTracePointsToASeqPos(in LocalAlignment localAlignment)
        {
            return localAlignment.tracePointsUpTo!"contigA"(
                aSeqPos,
                ac.tracePointDistance,
                RoundingMode.floor,
            );
        }

        void addTracePointAlignment(
            in AlignmentChain.LocalAlignment la,
            in AlignmentChain.LocalAlignment.TracePoint tracePoint,
        )
        {
            auto aSeqBegin = aSeqPos - begin.contigA;
            auto aSeqEnd = nextTracePoint(la, aSeqPos) - begin.contigA;
            auto bSeqBegin = bSeqPos - begin.contigB;
            auto bSeqEnd = bSeqBegin + tracePoint.numBasePairs;

            auto tracePointAlignment = findAlignment(
                aSequence[aSeqBegin .. aSeqEnd],
                bSequence[bSeqBegin .. bSeqEnd],
                indelPenalty,
                No.freeShift,
                memoryLimit,
            );

            appendPartialAlignment(tracePointAlignment, aSeqEnd, bSeqEnd, No.freeShift);
            aSeqPos = nextTracePoint(la, aSeqPos);
            bSeqPos += tracePoint.numBasePairs;
        }

        void addGapAlignment(in AlignmentChain.LocalAlignment afterGap)
        {
            auto aSeqBegin = aSeqPos - begin.contigA;
            auto aSeqEnd = afterGap.contigA.begin - begin.contigA;
            auto bSeqBegin = bSeqPos - begin.contigB;
            auto bSeqEnd = afterGap.contigB.begin - begin.contigB;
            auto aSubsequence = aSequence[aSeqBegin .. aSeqEnd];
            auto bSubsequence = bSequence[bSeqBegin .. bSeqEnd];

            if (memoryRequired(aSubsequence, bSubsequence) <= memoryLimit)
            {
                GC.collect();
                auto gapAlignment = findAlignment(
                    aSubsequence,
                    bSubsequence,
                    indelPenalty,
                    Yes.freeShift,
                    memoryLimit,
                );

                appendPartialAlignment(gapAlignment, aSeqEnd, bSeqEnd, Yes.freeShift);
            }
            else
            {
                auto numIndels = cast(score_t) absdiff(aSubsequence.length, bSubsequence.length);
                auto numSubst = cast(score_t) min(aSubsequence.length, bSubsequence.length);
                auto indelOp = aSubsequence.length < bSubsequence.length
                    ? EditOp.insertion
                    : EditOp.deletetion;

                logJsonWarn(
                    "info", "faking too long alignment gap",
                    "gapSize", aSubsequence.length,
                    "acFingerprint", [fingerprint(ac).expand].toJson,
                );
                auto fakeAlignment = typeof(_paddedAlignment)(
                    0,
                    chain(indelOp.repeat(numIndels), EditOp.substitution.repeat(numSubst)).array,
                    aSubsequence,
                    bSubsequence,
                    indelPenalty,
                    Yes.freeShift,
                );
                fakeAlignment.score = fakeAlignment.computeScore();
                assert(fakeAlignment.isValid());
                appendPartialAlignment(fakeAlignment, aSeqEnd, bSeqEnd, Yes.freeShift);
            }

            aSeqPos = afterGap.contigA.begin;
            bSeqPos = afterGap.contigB.begin;
        }

        void resolveOverlappingAlignments(
            in AlignmentChain.LocalAlignment lhs,
            in AlignmentChain.LocalAlignment rhs,
        )
        {
            auto overlap = getNonOverlappingCoordinates(lhs, rhs);

            assert(aSeqPos >= overlap.begin.contigA);
            assert(bSeqPos >= overlap.begin.contigB);

            auto aSeqBegin = overlap.begin.contigA - begin.contigA;
            auto aSeqEnd = overlap.end.contigA - begin.contigA;
            // Crop alignment according to aSeqPos
            _paddedAlignment = _paddedAlignment[0 .. aSeqBegin];
            // Use resulting position on contig B as it may differ from
            // trace points due to previous overlap resolution
            auto bSeqBegin = _paddedAlignment.query.length;
            auto bSeqEnd = overlap.end.contigB - begin.contigB;

            GC.collect();
            auto overlapAlignment = findAlignment(
                aSequence[aSeqBegin .. aSeqEnd],
                bSequence[bSeqBegin .. bSeqEnd],
                indelPenalty,
                Yes.freeShift,
                memoryLimit,
            );
            appendPartialAlignment(overlapAlignment, aSeqEnd, bSeqEnd, Yes.freeShift);
            aSeqPos = overlap.end.contigA;
            bSeqPos = overlap.end.contigB;
        }

        auto getNonOverlappingCoordinates(
            in AlignmentChain.LocalAlignment lhs,
            in AlignmentChain.LocalAlignment rhs,
        )
        {
            alias ResolvedCoords = Tuple!(
                TranslatedTracePoint, "begin",
                TranslatedTracePoint, "end",
            );

            alias resolveOn(string contig) = () => ResolvedCoords(
                lhs.translateTracePoint!contig(
                    mixin(`rhs.` ~ contig ~ `.begin`),
                    ac.tracePointDistance,
                    RoundingMode.floor,
                ),
                rhs.translateTracePoint!contig(
                    mixin(`lhs.` ~ contig ~ `.end`),
                    ac.tracePointDistance,
                    RoundingMode.ceil,
                ),
            );

            if (
                lhs.contigA.end > rhs.contigA.begin &&
                lhs.contigB.end > rhs.contigB.begin
            )
            {
                auto resolvedOnContigA = resolveOn!"contigA"();
                auto resolvedOnContigB = resolveOn!"contigB"();

                return ResolvedCoords(
                    minElement!"a.contigA"(only(
                        resolvedOnContigA.begin,
                        resolvedOnContigB.begin,
                    )),
                    maxElement!"a.contigA"(only(
                        resolvedOnContigA.end,
                        resolvedOnContigB.end,
                    )),
                );
            }
            else if (lhs.contigA.end > rhs.contigA.begin)
                return resolveOn!"contigA"();
            else if (lhs.contigB.end > rhs.contigB.begin)
                return resolveOn!"contigB"();
            else
                assert(0, "unreachable");
        }

        void appendPartialAlignment(PartialAlignment, FreeShift)(
            PartialAlignment partialAlignment,
            in size_t aSeqEnd,
            in size_t bSeqEnd,
            in FreeShift freeShift,
        )
        {
            if (freeShift)
            {
                auto firstSubstitution = partialAlignment.editPath.countUntil(EditOp.substitution);

                _paddedAlignment.score += indelPenalty * (firstSubstitution < 0
                    ? partialAlignment.editPath.length
                    : firstSubstitution);
            }

            _paddedAlignment.score += partialAlignment.score;
            _paddedAlignment.editPath ~= partialAlignment.editPath;
            _paddedAlignment.reference = aSequence[0 .. aSeqEnd];
            _paddedAlignment.query = bSequence[0 .. bSeqEnd];
            assert(
                _paddedAlignment.isValid(),
                format
                    !"exact alignment invalid after stitching: %1$d %2$d [%3$d, %5$d) [%4$d, %6$d)"
                    (fingerprint(ac).expand),
            );
        }

        coord_t nextTracePoint(
            in AlignmentChain.LocalAlignment la,
            coord_t contigAPos,
        )
        {
            alias isAlignedTracePoint = (pos) => pos % ac.tracePointDistance == 0;

            return min(
                la.contigA.end,
                isAlignedTracePoint(contigAPos)
                    ? contigAPos + ac.tracePointDistance
                    : ceil(contigAPos, ac.tracePointDistance),
            );
        }
    }

    return AlignmentPadder(ac, begin, end, aSequence, bSequence, memoryLimit).paddedAlignment;
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    enum complement = AlignmentFlag.complement;
    enum tracePointDistance = 100;

    auto ac = AlignmentChain(
        0,
        Contig(1, 1092),
        Contig(12407, 12767),
        AlignmentFlags(complement),
        [
            LocalAlignment(
                Locus(0, 439),
                Locus(6604, 7076),
                67,
                [
                    TracePoint(16, 105),
                    TracePoint(17, 106),
                    TracePoint(11, 108),
                    TracePoint(19, 111),
                    TracePoint( 4,  42),
                ],
            ),
            LocalAlignment(
                Locus(400, 439),
                Locus(7034, 7076),
                4,
                [
                    TracePoint( 4,  42),
                ],
            ),
            LocalAlignment(
                Locus(590, 798),
                Locus(7263, 7475),
                41,
                [
                    TracePoint( 1,  10),
                    TracePoint(20, 107),
                    TracePoint(20,  95),
                ],
            ),
        ],
        tracePointDistance,
    );
    auto begin = ac.translateTracePoint(400, RoundingMode.floor);
    auto end = ac.translateTracePoint(600, RoundingMode.ceil);
    auto aSequence = "agtgctggccccagtcaagggcatgtacctgggttgcaggttccccagcc" ~
                     "cccgtaggggtgtgtgtgtggaaggcaaccaattgatgtgtctcctttat" ~
                     "gttgatgtttctctttctctctccctcctccctgtcttccactctctcta" ~
                     "gaaatcaatgggaaaaatatcctcaagtgaagattaaaaaaaaaaaaagg";
    auto bSequence = "agatgatggccccagtgcaagggcatgtacctggtgttgcagcggttcca" ~
                     "gccagcacccgtacggggcgtgcgtgtggaaagggcaaccaatttgatgt" ~
                     "ccctacgttttcatgttgatgttttctcttttctctctctctctcccttc" ~
                     "cttcccacttttcttcccagcttctctctagaaaacaagggaaaagtgtc" ~
                     "ctttcagtgtaggatttttaaaaaaagccaaaaaaaggg";

    auto exactAlignment = getPaddedAlignment(ac, begin, end, aSequence, bSequence);

    assert(exactAlignment.toString(50) ==
        "ag-tgctggccccagt-caagggcatgtacctgg-gttgcag--gttcc-\n" ~
        "|| ||*|||||||||| ||||||||||||||||| |||||||  ||||| \n" ~
        "agatgatggccccagtgcaagggcatgtacctggtgttgcagcggttcca\n" ~
        "\n" ~
        "-ccagcccccgta-ggggtgtgtgtgtggaa-gg-caaccaatt-gatgt\n" ~
        " |||||*|||||| ||||*|||*|||||||| || ||||||||| |||||\n" ~
        "gccagcacccgtacggggcgtgcgtgtggaaagggcaaccaatttgatgt\n" ~
        "\n" ~
        "gtct-ccttt--atgttgatgttt-ctcttt-ctctctc-c-ctcc-t-c\n" ~
        "**|| |*|||  |||||||||||| |||||| ||||||| | |||| | |\n" ~
        "ccctacgttttcatgttgatgttttctcttttctctctctctctcccttc\n" ~
        "\n" ~
        "c--c----tgt-cttcc-a-ct-ctctctagaaatcaatgggaaaaatat\n" ~
        "|  |    |*| ||||| | || |||||||||||*||| |||||||*|*|\n" ~
        "cttcccacttttcttcccagcttctctctagaaaacaa-gggaaaagtgt\n" ~
        "\n" ~
        "cct--caagtgaag-att---aaaaa-----aaaaaaaagg\n" ~
        "|||  || |||*|| |||   |||||     |||||||*||\n" ~
        "cctttca-gtgtaggatttttaaaaaaagccaaaaaaaggg");
}

unittest
{
    import std.range : repeat;

    alias LocalAlignment = AlignmentChain.LocalAlignment;
    enum complement = AlignmentFlag.complement;

    // TODO reduce test case to two consecutive "overlap" scenarios
    auto ac = AlignmentChain(
        574,
        Contig(5, 33046024),
        Contig(344, 73824),
        AlignmentFlags(complement),
        [
            LocalAlignment(
                Locus(8372381, 8419325),
                Locus(0, 46966),
                388,
                array(chain(
                    only(TracePoint(3, 19), TracePoint(4, 102), TracePoint(5, 101), TracePoint(3, 99), TracePoint(4, 100), TracePoint(5, 98), TracePoint(3, 99), TracePoint(2, 100), TracePoint(3, 99), TracePoint(1, 101), TracePoint(5, 101), TracePoint(5, 103), TracePoint(8, 102), TracePoint(8, 102), TracePoint(6, 101), TracePoint(3, 101), TracePoint(3, 99), TracePoint(7, 103), TracePoint(5, 99), TracePoint(4, 99), TracePoint(2, 102), TracePoint(5, 101), TracePoint(6, 96), TracePoint(3, 99), TracePoint(4, 100), TracePoint(7, 101), TracePoint(1, 101), TracePoint(4, 98), TracePoint(3, 101), TracePoint(4, 102), TracePoint(3, 101), TracePoint(4, 104), TracePoint(7, 105), TracePoint(2, 102), TracePoint(8, 102), TracePoint(3, 99), TracePoint(10, 104), TracePoint(5, 99), TracePoint(5, 95), TracePoint(3, 101), TracePoint(3, 101), TracePoint(4, 99), TracePoint(3, 99), TracePoint(4, 100), TracePoint(7, 98), TracePoint(8, 105), TracePoint(4, 98), TracePoint(6, 98), TracePoint(8, 100), TracePoint(4, 100), TracePoint(8, 101), TracePoint(2, 98), TracePoint(5, 99), TracePoint(3, 101), TracePoint(3, 98), TracePoint(1, 101)),
                    TracePoint(0, 100).repeat(17),
                    only(TracePoint(8, 104), TracePoint(2, 102), TracePoint(2, 100), TracePoint(3, 101), TracePoint(2, 98), TracePoint(9, 97), TracePoint(3, 99), TracePoint(4, 98), TracePoint(2, 98), TracePoint(2, 99), TracePoint(4, 101), TracePoint(3, 99), TracePoint(5, 101), TracePoint(4, 99), TracePoint(5, 100), TracePoint(4, 98), TracePoint(1, 99), TracePoint(1, 99), TracePoint(3, 99), TracePoint(3, 101), TracePoint(5, 101), TracePoint(4, 100), TracePoint(4, 102), TracePoint(1, 99), TracePoint(4, 104)),
                    TracePoint(0, 100).repeat(138),
                    only(TracePoint(3, 102), TracePoint(5, 99), TracePoint(4, 100), TracePoint(5, 101), TracePoint(5, 97), TracePoint(1, 101)),
                    TracePoint(0, 100).repeat(163),
                    only(TracePoint(7, 103)),
                    TracePoint(0, 100).repeat(47),
                    only(TracePoint(3, 102), TracePoint(5, 99), TracePoint(1, 99), TracePoint(4, 100), TracePoint(4, 104), TracePoint(4, 99)),
                    TracePoint(0, 100).repeat(11),
                    only(TracePoint(0, 25)),
                )),
            ),
            LocalAlignment(
                Locus(8419319, 8419520),
                Locus(46907, 47119),
                27,
                [
                    TracePoint(10, 87),
                    TracePoint(16, 104),
                    TracePoint(1, 21),
                ],
            ),
            LocalAlignment(
                Locus(8419355, 8446270),
                Locus(46907, 73824),
                31,
                array(chain(
                    only(TracePoint(6, 47), TracePoint(24, 99), TracePoint(1, 101)),
                    TracePoint(0, 100).repeat(266),
                    only(TracePoint(0, 70)),
                )),
            ),
        ],
        100,
    );
    auto begin = ac.translateTracePoint(ac.first.contigA.begin, RoundingMode.floor);
    auto end = ac.translateTracePoint(ac.last.contigA.end, RoundingMode.ceil);
    auto aSequence = "aaagaaggtgtagagattagccaaagaacatatatgcaaagcccacgaacccagaaaacagtgtggtgaaggccagagaggggtgggtagggactgggtggaagtgggcaaaggggggaggaatgggggacatgtgtaatagtgttaacaataaaaaacccaaaaccaaaaccaagaaggaatattaccactattcactattaatgattaagacatcacatgtgactcatcaggaaaaaacagtttacagtttcatatcctaagtggagtacgtacatttaatagttaaacaagttggctttgggatatggaaattataagggagaataaatatttaagcccatgatgcaatttaatttattttttctgtgacttactgaactaagcatgatcagtgagactatgacttagagaatcaaatgtctcctttgagccttgctgaggcaagatttatatttctcagatgtctacacttggaagacatgatcagaaataaaaatgtctgtgaaacgccgccgattgacttgcacggtttgcagaacaggatccagttttaagtgctgtcttttgtgttcctgcagtagctgtgagggagggaactctggtccatgcaaggaaaaggagcagatagctccctgttccaacaagaacagtttgtgcaattaagtttttcacagccaccttgcatgatagatagtaggctcatttcacggaacttttatcaacgtgtaaatccaatcagcactccaacatcctgtttattagaaggtttctgaaatcatagctccattgtgtgtgtgtgtgcgcacatgtggatgtgggttgtcattaatgtgcctgcttcccactgcagcccctgcttgataaagaatgtgaagtatgagggtggagcatgactgaatttttctgagttgagaagcccggtcagccctccctggcatggctgagctaagcgtgtgtggtcccatttgtgtgtcacagttttggtcctttggcaatcaggctgagagtattttgtgattctgcatcaagtttatctgtttaatattgtcacctgcattttaagggcctccttttgcattaaatagctgtgtcaagcaatcaactggatttctttgagatgtgccaagtgtataacactgaactaccaatgtgaccatagaaatccagatacttgtacatcataaagttgtgtattataataaaatggacaatattgataaatttaaatgtgacatgtaccatttgccccatggttttactttattgaggataaattattactgataaatttgttttattattgtgcaatttgttgttgtttatcctcactggattttttttcattgatatttttagaaagagtggaagggaggggtggagggagagagaaaaatatcgatgtgagagagacacatcgattggttgcctgcccaggtatgttcccttgactgggaatcaaacctgagacccttcagtgcacaggctgatgctctaaccgctgagacacactggccagagcatcttatgcaattttcaaaatactttgtggtgagttgtgagagagcagacaagaacaggaattaaaatcagggacgtagacccagaaggatcatgcagtcaaatgtttggactgttaacaaataaactatacttgtaagtatgaaataattttaaatcttatcatgttatgttttcaaattttggtggtagagtggagttttagtcaatgttgctactgaattatttacaaataatatgaaaataccaggtgataactacaatggaatgtgtgacactgtagctgaatagtataacacaccactcagatttcaaggtcacctgtgagtcatcttttggggacaatttccattaggattttaatttttgtcatctgagtcttaggaacccaggcaatggcttcatgggtttgcttatggaaactgtctggcttgggccatctgagagaaatttttatttgtaatgcttggtttaggtctcttgattaaaaaaaaataatttcattttatttcattcgttgggataatttttgtattgaaccttcccctctcacaagtttaaagtataatttacatccagtgaaatgcacacatcctaaatgtatacaaatatataagtactagtggcccggtgcacaaaattcgtgcatggagcggggggttcccttagcccagtctgcaccctctccaattggggaccccttgggggatgtcccacagggatttggactaaaccggcagtcagacatccctcttgcaatcggagactgctggctcctaaccactcacctgccttcctgcataatcgtccctaactgcctctgcctgcctgattgcccctaacccctctgactgcctgcctgatcacccctaactgccctcccctgccggcctgatctcgccccctacagccatcccctgcagcctgatctcgcccccaactgccctctcctgccagcctgattgcccctaactgcctctacctgcctgataatccctaactgccctcctctgctgtcctgatctgcccccaactgctctcccctgcaggcctgatcttgcccccaactgccctcccctgctggccaatttggttctgattggtcagtttctatgccagtcagcatcaaaagctctgcctcctaggcagccattggctcctcactgttcacccagatttggttctgattggttagtttctatgccagtcaacatctctgggcctatctccaggcctgatcagagaggtgggactgataagcagccccctcacggaggccttgagagaaagaggcgtggctgctggtgaagcccggggagagagagagagaggcaacagttgatcaacagctgccatggaggctgcagatcagaccgtgcctctctctctgggcctgatccgcagctccctcggcagtcagtgctgggttgctgtgcccgcaccggcggtagtcagtgctgggtcgccacggcagcccagcactgactgcaggactgctggtgttcagttgagccttcggtaggtcgttatggatcctgggtttttatatattaggatacttacaccatatgtatatgtactcaaatatccatgtgactaagacctggacaaaatgtaaaatgtcaccgtcagctgagaaagttcccctgtgcccctttcaatcaataagcccacctttctcagaagataaccactattctgacttctatcatcacaacttgcttgtttattaacctcacgtaaatgtgattatacattagctacccttttctgtctggcttcatttgttcaacataatgttattgagactcatccatgttggagaatgtatcagtaggtacacctggtattatcaatcttttacattttagccatatagtatatgtatatatatgtaaacattcttaattggtatttttttatggtttaaagtattacatatgtctcccttttccctgattgacccgccccccaaccactcccacccctgaggaaaagtccccaccgccccagtgtctgtgtccattggctatgctaatatgcatgcataatgcatacaagtcctttagttgatctctaaccaccgcccatgccacttgtccctggatctatttttgttcacagtttttgttgatcattatatcccacatcttagtgagatcatgtgatatttatctttctccaactggcttatttcccttagcataatgctctccagttccatccatgctgtttcatccatcccccagcccagcttgcaccctctccaatccagggccccttgggggatgtccgactgccggtttaggcccgatccccagaataaggcctaaaccggcagttggacatccctctcacaatctgagactgctggcttctaactgctcacctgcctgcctgcctgatcacccccaactgccctcctgctggcctggttgcccccaactgccccccgctggcctggttgccccatgcagcctgctgttcagtggtttggtcgtccctctctaacccccctgcaggcctggtcaccccacgcagcctgctgttcagtcatttggtcatccttcactaaatccctgcctgcctggtcgccttatgcagcctgctgttcagtcgtttgattgtccctcactaatccccctgccggccttgtcgccccatgcagcctgctgttcggtcatccagcccgttgttttggttgtgacagcccctggctttttatatattaggaatagtattccattgtgtaggtgtaccacagttttttaatccactcatctgctgatgggcacttaggctgtttccaaatcttagctattgtaaattgtgctgctatgaacataggggtgcatatatcctttctgattggtgtttctagtttcttgggatgtagtcctagaagtgggattactgggtcaaatgggaattccatttttaattttttgaggaaattccatactgtcttccacagtggctgcaccagtctgcattcccaccagcagtgcatgaaggttcctttttctccacatccttgccagtacttgtcctttgttgatttgttgatgatagccattcttacaagtgtgagatggtatctcattgttgttttgatttgcatctctcagaggattagtgactttgagcatgttttcatatgtctcttggccttctctatgtcctcttttgaaaagtgtctatttaggtccattgcccattttttgattggagtgtttatcttccttttgttaagttgtatgagttctctgtaaattttggagattaaacccttatctgaaatagcattggcaaatgttctcccatgcagtgggctttcttgttgttttgttgatgatttcttttgctgtgcagaagctttttattttgatgtagtcctatttgtctatgttctccttattttccattgccctagaacagccgtgggcaaactacggcccgccggccggatccggcccatttgaaatgaataaaactaaaaaaaacccaaaaaacagaccatacccttttatgtaatgatgtttactttgaatttatattagttcacacaaacactccatccatgcttttgttccggccttccagtccagtttaagaacccattgtggccctcgagtcaaaacgtttgcccacccctgccctagaagctgtatcagtaaagatattgctgtgacatatgcctgatattttgctgcctatgtagtcttttaagatttttatggtttcccatcttaaatttaagtcctttaaccattttgagtttgtttttgtgtatggtataggttggtgatctagtttcatttttttttttggtatgtatttgaccaaatttcccagtaccatttattgaagagactgtcttaactccattgtatgctcttgcttcctttgtcaaatattaattgagcataatggcttgggtcaatttctgggttctctgttctgttccggtgcaaaagtatatttgtatacaaaagtatgaaggttaatttcatttgtcaacttgactgggccatgagtgctcagaaatttggttaagcataattctggatgtgtctgtaatggtgtttattctctgcatgagttcagcatttgaatcgatagactgagtaaagcaaataaccttctttaatgtaggcaagccccatcgaatccactggaggactggatagaataaaaggcaggccaagaaagcattctttctctctgactgcctatctttaagctgggacatccatcttttcccgcctttggacttagacttgaactagaatttacactagtggctctcttcgttctcaggcctttaggctcaggctggaactataccctcagctctccttgtgtccctaccttttgactgtagatcttaaacttctcagcctccattattgtgtgagctagttccttataataaaaccaaacaactaatctctctctctctctctctctctctctctctctctctctctctctctctctttctctcaatatatatatgtacacacacacacacacacacacacacacacaaattcattctttctttattataaagtgatctaggcttactggttatttctgcctccttttcagatgttatcttatctacccatttccctgatgttgtaccatgaaacatacaagacaatggctaacaaagagtctttcttgcggggaaggaaggaggcacaagctaatacatattgtctctgctgcaaataattttaaaaataacatcaagatgctttctgataacagttaacaatcactttttatgaagctacaaacactcagctctccacattatgagactctccaagagtagtgttttctaacagatattttctgaagtcatttcctctagaggacacaatgagtcaccttaaaatggctcaaaactcattccctcctctccctcgccatgaaatgacttacaacatttcattgttacttcaactttagggacaacctttcttagtattattactgtgtctggtatccaattaaacaataaggagtaatgtagaaagaattgcccataaaatgtatttatttttctgaattgaaaattaaaaataagagtttttgaatgtgctcagtgagagtctatacctactttgcttccaatatttgagttcctctataggactgtagggacacaactagatgctgatctgatgcccagagcagagtggcagggtggagggacgcatgtgctgtctatgtcactggcaggtcagctttcctctccctgtctccagttgcaagagtgtggacaggcttgtatgatccagcgcacacctggatactcacatctaaagaaagagaccagagctttcctgtttgtttcactcatcaagtcagagtcagagcttcctgttccctatatactttggctgtcactgacaggtctaaagtttcatggtgatataatcctatctaataaaagagaaacatgcaaattgaccataaccccaacacccttcccaggatatgcaggatatcagcaggatatgcaaattaactgccaaccaagatggcggccagcagccacgcagctgaagcgaacaggaggcttgcttgctccagtgatggaggaagccaacattcccggcctgccatgaccagcctctgagctacactctaagcaactatgttgcaattatagaagctaaacaagctccagagacctgctttcagccagcttaggcctcagagcttagagctgcagtgatggcaacagagtttcaattatagaagccaaacagacccagaacctgctttcagcccccagagctggagcctgctctcagctccagtgacagctatacaaggtaaataaatcccagaataaaaaaaaaaagaaaaaaaaggagaggctgggagcttaagttgccctccagcctgaaaatggctctcaggccctcacccagactggccaggcaccccagtggggacccccacctgaagggggtgtgaccagatgcaaacagccatcatcccctcatgcaggctggccaggcacccaagcgggacccccaccctgatctgggacactttcagggcaaaccagccagcccccacccgtgcaccaggcctctatcctatactactggcagggcgggacagctcgagctgccatccgggccccatgtgcagccctgggtagtggggcataagtcctgcccccagccgggatccctgtgtgccacctgatcccccagctgggatccccgtgctgcctaagccctgcactgctaagccccgcccccagccgggaccctatgactattggagagcaggaaagtgccgctgggccctgggtcactgtgttgatctggccatcggtaagtagaaaagtgcagcctgtaggcagcccccagctgcacctgatctgggagcaggaaagcaagcctggcctgtgggcagccccaagctgggcctgaaacgggagcaggaaagcctggcctgcgggcagcccccagataggcctgcttcccagtcaccatggggatccgggaacaggagaacaagcccggcctgcgggcagctccaagctgggcctgctccctggttgccatggagacctgggagcaggaaagcccagcctgtgggcagctcccagctgggcctgctccccggacgccatggcgatccaggagcaggaaagtctggcttgcaagcacccccatcctggccatagctcaagggagagtacaacatgcagtggaggcccggagcctaagagatcaacacagaggggctcctgctccaagcctctatggtgtggctgggggcaggcccccagcagacacacacacacacacacacacacacacacacacacacagtgcctgctctgagcctccaatgtggctggggacagacccccagtggggacacacacacacacacacacacacacacacacacacacacacggcaattgctttgagcctacgtggcgtgccgggagccggtccatgcttgctgtttcaagggacctggcatatatggcatactgttcttaatatgtttgctcaccttcttggcactatgtgttttaaccaaggtctctgagaaaggttgtttccccaggtagggatttttcactgaagttagggagggaataaaacctcttaactaagtgccaggcgggtaattaatcactttaactatgaacaatcatgcttaagctacataatctttactccctggaatggagataagaaacgccctaacctttggaatagagattgataggattggaatcaactggtataaatacagatgtaacaacacagaacttaggagacagaacttagaacacagaactaagaagacagaactaagaacacagaacctacacagaacctagagacagaagaactttgctggagagaacatggcaaaagatcctggactgaacctgactacagaaattggcaagagaacctgactagaacctggtgactgagcctggctggagaacctggacacaacctggctggagaacctagtgagggaacatggctacagaacctggctggagaacctggagaacctagcaagagaacatggacacagaacctggctggagatcctaaccagaacctcactggagatcctgaccagaacgtggctagagaaccttgctatgatgatcatctgaatgccctctccgtgtcattccttcttcgctgactccgtccacacctttggggacccctggacccgctggggttggaccccggcagtggcgcggatgggggcaggccctagcagacacacacacacacacacacacacacacacacacacacacacgggacgcctgctccgagcctctgctgccagactgtggccatcagtaggacatccactgagggctcccggactgtgagagggacaggctaggctgagggaaccccacccccccagtgcacaaattttgtgcaccaggcctctagtcctatataataaaaggctaatatgcaaatcaacccaatggtagaaggactggttgctatgatgcgcactggccaccagggggcagacactcaatgcaggagctgcccccctggtggtcagtgtgctcccacaggggaagcgccgctcagccagaagccgggctcatggcagtctgacatcccctgagggctcccagactgcaagagggtgcaggccaagctttgggaacccctccctgaatgcacgaattttgtgcactgggcctctcgttttatatgaaaaggaaacatggtgtgtaccccctgaaaagccttcatgtcacaaaggaaaaatgccgtgaccaaaaaggcaaatttgttccatgcaagttgtctactgaaatcatttcctatttccaaatgtcttcattaaaacaaaacagataataaatcacttacacatctcccaacaatcctctcaatttaaaaaaaaaaaaaaatagagccctggctgatgttgctcaatggttggagcgtcggcctgcatgttgaagagtctcaggttatatttccaatctagggcacgtaagttcccttgggttgcaaatttgatccttggcctcggttagggtgcacacggagggaaccaatcgattgtctctcttacatccatatttctctctgtccctctctcaccccctccctcccactatctctaaaaatcaatgggaaaaatatccttgggtgaggattaagaacaacaacaaaaaacagtaaaaattcagtcacatgattttccattttgaatgaattaatctaagtggatgggaaacccaaatgcagcagagtccaggaggccattttcagttagctccattcaaaatgattaggaacactttttagaaaagtcttctttttgagtaatgccattggatatagaaaagctgctttcccaaaatgaaaaatcattggtaaaaataactgacatttatcttccagcatgatttgatttttaaaaaaatatattttattgattttttacagagtggaaggaagaggaatagagagttagaaacatcctgcacaccccctattggggatgtgctcacaactaaggtacatgcccttgactggaattgaacctgggacccttcagtccacagaccaacactctatccactgagccaaaccggttaccaccccagcatgattttaaatagcgttcctatatgtcttaatttattaagaaaagcatcaataccaaaatatttcattctatatcataaagtagaatgtacttttaataaaaatcggaggaaaaatttcagatgtttaaagaaaatgtagtttgacataataatttggagctctatccttgaggtgaaaaatgaattttactaaaagcaatataatttttttctggaaattaagcaaaggattattttttttttgtggtctcttagccaacggacaagggggaaaagaaaagaactggtctattctttttagaatcagagagaaggcaattgaaagttccctgtgaaagtaattttttttcccctcaagggctcttgagggttcacactgaaatgtcaatgaccttgtttgttcagagtttgtaacagccattttcaactgccctctgtgcttcgataaatatgttccagtcattgtaggttctgtgtctaaaacctatcaacatcattccctaagcatcaaggggatgctttttatgacaaaacatgctttatgcaggactctcactaaacccatgtgattccgagaacatccatcccaaagtgcagcgttaagaaattacacaatataaaattatttggaggcagttctggtgaattattttaatagtcctaaggcttggaaattccagacaatctcaacaaaccatacaacggacaactggcaaagatggtctatattttttaagaaaagccaagttattatttcaccatcttttcttaggataggtttcaaaattctgatcaactctagaccaccaggtagcccaacttagaacagtttaccagaagggccttacttatttctgctacttagtatagaaaatttttatttgagtgaagtaagcataaagagactttgactgatttgggattttcctcgtcaacacatattaaaaactaaaattaaattatatatattttaatttattttctctccatcataatcataaatgatagctactatttatggaattcttctattcctaggctttgcaataatactgttccatgtattttatcacttaatatttacactgatgtgtgaggtgttccacaatgatcttttcttatcactgaagaaattgaaaattagagaggttaagtaatttatccaatgatcttgctcaagaactccaagacacacacacacacacacacacacacacacacacacacacacacactcatatatattaaccactctgtcatattgcttagtatgtttaactaggggcccagtgcacgaattcatgcaccttgaaaggaactgtgggctgcgaggctgcagtgggcacaggggcgggtctcggctcatcctccatgcccttgcctggactctcccactgcacccccagtcccttgtctgccagcagccctgctcccgccaccaccactccagcgcactgacagtgctggccctgctcgtacccactgaaggtgcggagcaattgggtcctggcaccagcaatgggtatgagtagggctggcgccatcatcgggtgcaagcagcagttgctgccctgttcgcccctcaggagcaggtggaggtggagaagccctcaggggcaattgggattggcagccacagcttgcacccactgatggcaccgagagattggggctggcgccgggtgccagcagtgggtgcgaatggggccagcactggcagcaggtgcaagcaccaggcaggaccatggtacacaagagcaaagaattttcagtaaccaccagatgcttgcctcaatgacagcaaccggcgccctgccttggtctggtgcccccactcacctgctccaccatcccgtcgcagcggatgcctgccatgttttgcgcacgccccctggtggtcagtgcatgtcatagcaactggttgttcggttgttcggctgttctaccatttggtctatttgcatattagccttttattatataggatactgtcgtggaccttttatgattattctacatcagcaggctataaataggctgcagaggatatggatgtagaactgtaaaaaccagaaacacaattctctttttccttctttttcccagacctcccttcttcctctcttttttcttctatgccacaaggtctgccctcctatggtttgaaatataattgagaagacaacatatgtaattctcaaaagtaaaacatagaaaagtttgataacattgagtaaagatggagaggaagtaggtaatcatatttttttaatttctgaaatattccaatagccaaataaatgaatcattctcaaaaataggaaaataaagttttaataggaaaagagaatggagtaactgggggaaatcacaatgtgcctttcttatagtttttaattcctttaagtttattgaatgattcagattttctattcctttttatgtcagtttgagtaagttgtattttttcctagaaacttgtccaatttacctaaattttcaaatttaggaacataaagttgttaataataaatatcattaaatgattttataataataaaaaaaatatctgtaggctctgtagtgagatatccctttccatttctggttcattttcttgaggtacaggacattgtatgtgaaaaaaaagttacaattatttgaggcctagaattatcttcctccaaggaagacttacatttgcctttggtaggtggttagggccccagctatcttagatcaccttcatctaattatagggactgaacccagctggcatggctcagtggttgagcatcaacctatgagccaggaggtcacagttaaattccccattagggcacacaccgggtttcaggcttcacctgcagtaggggtcttgcaggaggcagctgatcaatgagtctcgtcactaatgtttctatctctctcatttcctctctgaaatcaataaaaaatattttttaattatagggactgtgatgatttaaactttagtgcagctccaagttttttccagttactccttactcataggtggagccctttgggatcccatctagatgtgtgaacatttaccagagttcactcctgaacagaactttgattttttgctaaattatatatatatatatattttatttactttaaccatctagccctatgggtctgccaaaacctctgccgaaattctaccttaatcatctcttccagaaccatcaaccaatcctaggattcacctctctggtttccttttcttctttatttcataattctttactatcttggtggctataagatttaaatttttttatcttgctttgttaatttgtcttagcagaggtagtttatgttatgtgcatgttttgtgtatttttaaaaaaatattcctgcaaagcaaaggttcttcacactttgttttaaatattttaccattgatctgtgatttcttggggtcattttgcatatgttcatatttttccccataatatggaagttcattagggctatgacaagtgctgagatggctgggatcacccaattggtccaccagctggtattagaatcaacagtaatactaggagcttctgaagacttggttatcatttgtctgtcctctggacgaagctcccgagtgctatatgttttggacagttctctggcatctgtagagttccaggtgtcctcaaagtttttggtagcatcacctccagcttgttccctcaggacttctgccccaccaggatgctcctccaaaaatctggtcaaatcatgcaccttttggtgcaggatcggtcaggtgctcctgctgtggctgtgcttctggagcccttccagggtgtatgtctgatggccttgggaacgaggcaagcctctggccctgcagacaccagtgagctgggcagaggtaatagctgcaattttaaagcccactttcagaggaatgaacacaatgtaccagctactgagtagcatcatttttttttttttttttttgagaagtgaggttgtaaacttagaacaggataaaacatcagtgagagagaaacatcatcaattggcagtttcctgcatgcccccgactggggactgagccagcaacctggcatgtgccctgactgggaattgaacctgaatctgagacctcttggttcatgggttgacactcaaccactgatccatagcggcaagactgtttcttaacatttttgagcatcattctgagatcaaggatgttgggccgatgcaaattatgatcatcatggtaaagcagtgaaaattcttgattaccccgctcataacataaattcactcagtggttactcatacttgtttagaaggttccgctatagcccatgtattttggtaatattagcttcagattgttttggaaacttattcgggaagatctacctttcttgttcttttatatctctggatctaaatttcatgagtgtacttctctgatttggtctttgcaaataaaagcaggttgatcacttgcatttacctgtatagaaatctggtttggccttgtgatctccaggtaggcccaggagtatctttaactgtcacctttatcttcattaaatggtaataggtttttgaaaattgctcccaaaggagtcctactttgattgctggaatacctgaagttgtagttgacttctctgtgttaaataagaagggattgtgaaaaatagtgtaaaatgaagaaatttgaatgcaaatataataatttactttggggagtgctcacgatttatttttttttctttaggaactgatcattttttagtttctagtccactgcttctgactttacatcataatcccaattcacaaacacataggggtgcacccaccaccaccacacatttatatcctttccctcccccgctttgtctcactaaagtgaagtatggtagctatggttcttgtctgcctatcatgtctttccttaatccctgctttctagaaatagaaatcctctttctcataaatttattccatgtggtccaagtggacagtcttctcagaatcgttagcaatggggacaggaagcataaagctactagtgcccatcttggcaaggctgtctgacaataaagccaatgctcaggtaagagggcatgggctaaaaggagagtgagaatgcctgggtgatatagtctttttcagcttaaatttttttttttttttttttttttactatttttgtaagtactcaaatagtcataaaatcattctatgcacataaactttctcctttttaaaaatatattttttattgatttcagagaggaagggagagagagagaagtagaaacatcaatgatgagagagaatcattgactggctgcctcctgcacaccccctactgggaattgaaccgtgaccttctagatcagcgttcaaccacggagccatgccaggtgagctcttttgtagctttaaaatgaggggttggtccacagtcatgtcagtggctaccatacactgaaatttcattttatgctaaactaggggctcagtgcacaaattcatgcaccttgaaagaaactgtgggccgtgaggctaccataggcataggggcaggtctcagcccatcctccatgcccccgcccagcccctcccaccacagccccctggtcccctgtctgccagcagtcccactcctgccgctcccatgcactgacgatgccagcccctttcatacctgctgaaggtacggagtgattgaggctgatgctagcagcgggtgtgagcagcagctgctgtcccgatcacccccaggagcagggggaggtggagaagccctcaggggcgatcagggccagcagccacacttgcacccattgatggcattgagtgattgggactggtgcagagcatcggcagtaggtgtgagcggtggttccagcactggctgaaggtgcgagcagtggctccggtgctggtagcaggtggaagctctgggcaagaccacagtgcacaggagcaaagaattttcagtaaccacaggaggctcgacccgatgacaatgaatggcacctcaccttggtctggtgtcccctgctcacctgctccaccatcctgccatggccaacacccaccatgttccgcacatgccccctggtggtcaacgcatgtcacagtgactggttgtttggttgtttggttgtttcaccgttcagtctatttgcatattagccttttattatataggatattttccatgcattaatttcatttaatattaatgacaactctatgaggcatatgtttttatctttctcccactcctataactaagaaaactgaataccactgagaaaagcaatgtagtgtaagctagtaagagaagttatattcatcattttaagtttgtttctatcagtgtctggcttaatggtgagatatcaaaaataattgatcaagtaaatagataaagcacgatctctttgacttcaaagtccatgtccataaccaccatcaataacttaaaacattattaacttgagcaactgtctaaactcactgagtggtatctaggtcaataagtctatgtttctaacaaatttgaaatgattaaagatagaagcaaacatgtacttccctcaccattcctaaacacaaaactgaatctacccctggttctaatacaagaataggtgaaataggctttggggtttggcctgagacgctgaccaccctctggaacaaaaaaaaaagcattcaaaagccctagccggtttggctccttggatagagtgctggcctgcagactgaagggtccttggtttgatactggtcaggggcacatgcccgggttgcgggctcgattcccagtagggggcgtgcaggaggcagccaatccatgattctctctcatcattgatgtttctattatctctccctctcccttcctctctgaaatcaataacaataaaaaaaaatgcattcaaaggtctaaaaaaaaccccacaaaaaacaaaacccacagatttacaagcaaacttttgcaatcttgccaactcgaaactgaacacagagagagaaacatcaatgtgagagagacacaccaactggttgcttccagcatgcatcagagtgggggcagggagcaaacctcaaacccaggtacacgcccttgaccaggactcccagagatactatttaagaaagcattctagatgcagtgcagaactttatggcaagaaaactcagtgcttatttagtttctcaaatagagggatttaggttgtttgtgaggtccagtttctgtgcttcaaaaaataacacgtgacttagagaggttttgagttcatatcaaaacatggatcatgatcctaatactcagaaccagaccatttctgtatttataccacataaacagattctaagaactctatccccaacacttaaatcaccggtctatggaattgaaatgaaaattagaatcaggtgaaggtttatgtttaaggaacatctgatccaaaacaaaatgaaacttaaggcttaaattctgttcttcctgagattctagtttctccctgtcattgcaattggccaaaactctttgtgtgaagaaccagctctttggaaagcttttaggaaataactttcagagactgtgttctgctaacatgccatgaaccatagctgcaggaactatgaccaggtctccaagtcattttatgcctcctacttggtgtgggaattctatgatcaaatcccttagcgttagtttggttgccaaccctggtcccttccaatgaaattgtggtaattagagtatgttaaacagtccagtggcaccaactgtgtagtttgtttggtgaaacattcagggcttgtactgacaagtccagagaggtgttaagcataagatcaagaggacaagtaggaactggagaagaaatctattgaaagctaattttttttcaggcactttgtctggcttttcatatatgtgtatatatatcctatataataaaagcctaggtggcattgtgcaacgtcctcacatgatatcatcacaagatggtttccacgacattgtcacaagatggccaccacaagatggccgccataggatggctagcaagggagggaagttgtgggcaatcaggtcggcaggggagggcagctgggagcaatcaggcctgcaggcgagggaagttgcgggggaccaggcctgcaggggagggcagttgggagtgattgggccagcaggggagggtagttggaggtgatcggccggcaggggagcagttaggcatcaatcaagctggcaggagaagttaggaggtgatcaggctggcgggcagaagctgttaggggcaatcagacaggcaggcaggtgaatggttaggagctagcagttctggattgtgagagggatgtccgactgccagttgagggattcacatggattaggcctaaaccggcagtcagacattttccacgggggtcccagattggagagggtgcaggctgagctgagggacaccccccctcctgtgcacgaattttgtgcactgggctatacatatatatatatatatatatatatatatatatatatatatatatatatgtctcatctatttattcttccttaaatgtcctagaaacaactactcctgcacaattctggatgctaacaccatcagaaagatatttgtgatatttgtgcctggctggtgtggctcagtgattgagcattgacccatgaaccaggaggttacagttcaattcccagtcagggcgggctcaaaccccagtaaggggcatgcagaaggcagtgatcgagggttctctctcatcattgatgtttctctttctccctctcccttcctctctgaaatcaataaaaacatattttttaaaaaataaataaaacagagtgttttctggtttttgttgttgttttgttttgttttttgttttttttaagatacttgtcctcatgtgcctcagtcaggagacttctaagatcatgaattggctgctaatacctagcccctctctctggtctaacccatagactcagatcatctcatttcctcctggggacacaaacactccataaaggcagctcttctagggctctgagcaggaccatcctgacattgttttcactcctcttccagcctgcaagcaaatgaaaaacagctccatctcacaggagggctggggaaccttgataaagaacagatctccaattgccttttagactgaaataattggacctgttagatatttttgaaagatttatttactgcatttcaacaattaaggtgaattttataggcccatctctcgtgggtatatttgttttaatcatgaattaaagactataccctggctgatggctcagttggatggagtgtcatcccgaacaccaaaaggttgtgggtttgattcctggtcagggcacatacccaggatgcatgtttaatcccccgtcagggtgtgtctgggaggtaaccaatccatgtttctctctcacatcgatgtttctttctctctctgtgtcccttcctctgtctctaaaaatcaattaaaaacatatcctaaggtggagattaaaaaaaaaaaaaaagattagcgcacatttcccccctgcacaaacatcatgacccataaggttaacgcatgaaatcctttcatcacgaactaggattacactctattccaaatcaagccacgtcaacatagatgatgccagaccaccaaacagtgttttctcaggctctaaaagggcagagtaactgtcacccgtgggccatacgtgtggtccggtgggaagagtaagctttggaggatggaatccgagctcacgggctgggaggtccaggtccatcacttaatctcccgggaacccggtgctcagtgtcaggaggctcttgactcacaacagagtatagtagacaacaaagagagactatacattttcactggcttgaaaattcctagaagacagggcctgcgtgtttagcttgtagcaacaatgtttatcacttatctaataaatcctatattaacatgttcatcagaaaaaaacacctaaatgtatcgacaatcaagtatacctgttaagattctccttgcaaaagccacatagggaagacaaaacaaaactcaattcacattttccaggggacattgttcagtgtcattgaaacaaagcaacatctgttttcctagaagttactatatatgactaaatattctaattagtatttctcttaggccataaatattttaccataattttaaaaacttgataaactaatgtaaaagatggctctttgaaaaaaaattaataaaacagacaaacatcacagcaagattgattaggaaaaagacaggaaaaaacacaatgctagtatataattgcattagaaataaaaaaagaggctacaaatactgtgcatatttataaaaataacttgagacaacaaatataaaactaaaatttttctaacatatataccttttgtttttaagaagtgaagcaaggcagtgctgtggtgaagctggtatagttcatactgtgagggcaggttgtttaccctctttaagcctttgctcatttgtaaaatggaggtggtagtgtcctcacttcttatttacaggattattatgacatcagtggaattaatactgtacatgcaaaagaaatttaatctaaggaagagatgagagagagagagagagagagagagagagagagagagagagagaaagaaagagagagagaagacactagtagatattgagaataaatatgaatgaagagtctatgaatcaggctgacttggaaggaaagaatggaaatctgagtagaaccaatgatgactgagaagtatagcttttcccccagaaggtctaccactaacttcaccctcccagacaccagcaggctcatgtcgtaaggagtcccggtctcttagtttccttgactggggcagggagtccttaactcacttggacaggagaagggaagcatgaaaatattccttgtgttcccagtagtttcttgggcatattaagttagaaagacattattaaacatagtgttgggtgctagaatatttctggtgatgcttaagggacatggtaaataaaacctgagagtgaattctggcttgcagctgattatgggatggtattagatgcctggcatcttaggggttaaaaagtctttggtaggccagccagaaagtttaataagcatgtaaaacacattccaaacaaatgactgctggaacataaccctgtaagttggggactgcctgtttttataacctttgccaaagttactaaaattagcagcataatcgaaataaaaaaattatatattcttatgactgaacagcatggcagcatcttcaaacccagagctgaattttaaagtcaaaattattttttgggcattttctagaaccaaagagaatgatttggcaagccatcccctctggcaatatgaatggagaaattttgcttttatctttataattaaagtagctctgtgcagtttcattcccagtgtgaacgggggcacaggtggtatctgctgatgttacaaaacagagaaatagcaaatatgttaggtaaaagacacaaacattctgtttctcacccatctccttgatttctcaaaccaaaatgctaacagttcagaggatcattttgtatatgtaatggaggaaagcactggaatgtttagcttagaaaagacttcttataccagggaggaggaaagtcattaaatagaagtgacgttctttctgatttggtctatggcagtggttctcaaccttctggccctttaaatacagttcctcatgttgtgacccaaccataaaattattttcgttgctacttcataactgtaatgttgctactgttatgaatcataatgtaaatgcctgatatgcaggatggtcttaggcgacccctgtgaaagggtcgttcgactgccaaaggggttgcgacccacaggttgagaaccgctggtctgtggggatttagtggacacatctaagaccaatggtggtagattttctccagtgtaggaaaggttaacagaattcttcaacaacgagatagtaagtttccttcctttattcaaaggacgactatcatgaatcaagtacaggagattcttacactgggtagaaagttaaagcagatggtctacttttccaatactaaagttctgtaagttcatgaaaatcataataactagtccctctttctttcccaagcaatatgatatcctcaaaattgccaatgaggaaattcaagtggatgctaagaccccagttggtcatacaattgtgcaatcaagttttaagtgtttctgatggagccagccaacctgagactgcccatctccttactcatggagatgtagttttggaaaccaaaatctgtagtaggcaaaaactggtccaagccatgtggctttccaacagggtatgattataaaagaggccagaattattgtagctcagaaggggaaatgtcagcggaatttattatccttataaaacatattctttgatccttccttctatcgtattgcacgcctgtttcttggactgactagtggtcagtgaggtcaagataatgttgtatagatcacattgcctttattttggctttagagagactaagagccaagattagtctctttttagaatagttattatttaggagatgtaagtttttcattgttggaagccatttggagtgctctaataggtggcagttgcctctgcctgatctttgtaagtgatagcaccagggcgttgcagtgtgctcaggaagaaagtcacaatggaagtctatttgtcttctttacaaacccaccacacttgtaattattttttaaatatatttttattaatttcagagaggaagggggggggagagagagagagagagagagagagagagagagagagatcaatgatgagagagaatcattgatctgctgcctcctgcacaccccctactggggattgagcccacaacccaggcatgtgcccccaaatctgggacccttcagtccacaggctgacactctatctatggagccaaaccagctagggccacacttgtaatttttttttttttttacatattttattgatttttttacagagaggaagggagagagatagagagttagaaacatcgatgagagagaaacatcgatcagctgcctcctgcacatctcctactggggatgtagcccgcaacccaggtacatgcccttgaccggaatcgaacctgggacctttcagtccgcaggccgacgctctatccactgagccaaaccggtttcggctaattttttttatcaacattggccgtccatctggaatgcaatttcatgagagtaaggccaaattatcttgtgcagcattccatttccattagctagcacagtgccagcactcaaaaaatatccaattgaaacaataaatctagcaaggaaatgatgacccgctccctccaatgctctttatatcttccacatcaattcctctcccccttacacacatggggcaggatagtcctaagaataaaattatattggtttgtatattgtctacttattccataaacacctatgcataatatttcatttgatctttgtaataatacttcaaggtagttgacagatattattattattgctattcccatctttttttccatatacaaaactgcaggatttaaatgaccatttgagggttaaagtgctagctagaggcaatgctgggtcttgatctccaggagacaaagcagtgtgtgtagcaatagcaggactaaagcctgtggaattccaactccatctgttatctatgtaacctcggaaaaatcctttgcatttgcttagctttagttttttaaaaatatgtctttattgatttcagagaggaagggagaggggagagagagagagaaacatcaatgataagagagaatcattgctcagcagactcctgcatgcccccaactggggattgagcccatagtgctcaaccactgagtaaccatggccaagcagctttagtttttttaaaaatatatatttttaaatggggctaaagaggtaaaaagggaagggcagaggcaatatttatcttgttcacctcatggaattattatgagatggggataacaaaggtcgtagggtgtaactgtgagatggtgtcactctaagaagaataatacttcccaacagagatgaatagatggggctttggggcccgagggtccgagcaccatgcctttctagcggtatgacctcaactgagttcctgaacctttagtttcctcatctttaaaatgggaacaataaaacctgcctcataaggctgacatgtggattaaatgagattaagtacagaaaatagaagtaactaactgtatgataatagtggtagcagaataatttaaatagtaattttctctccagttgaagaaactgaggtttagaagcttactagatttggtttcagggtacccagagagcagagcagagacttggtagaccagtgtgaagtcttgaatttgcttcggaagaacagatctggtgagtggaagctggacgagaagaaacttttttcttttcctctccatcccaccagctccacccacctcagtccacaaggtctatttctttaggatttgggaaaacaggagccactggagggccaaggagaaaccagtacaagagcgatttccctcctttattgctgtcagaatcccaagaaatcgtagacaagagattggatgtttagtggtgccattgcttaaagcccaggttttccccaattactgggacaacatttccttccctggaggacctggtctgaatttctggcccccttaggcaactagaagtaaacctttgagataaaatcctcttgtatgttggtgctccttctgtgccgggagccggtccatgcttgctgtttcaagggacctggcatatatggcatactgttcttaatatgtttgctcaccttcttggcactatgtgttttaaccaaggtctctgagaaagggtgtttccccaggtagggattttcccctgaagttagggagggaataaaacctcttaactaagtgccaggcgggtaattaatcactttaactatgaacagtcatgcttaagctacataatctttactccctggaatggagataagaaacaccctaacctttgaaatagagattgacaggattgaatcaactggtataaatacagatgtaactagacagcaagacacagaacttagaacacagaacttagaataagacagaacttagaacacagaactaagaagacagaaccaagagagacagaaccaggaggcacagaacctacacagaacgttctctagagacagaagaacctcactggagagaacatggcaaaagatcctggactgaacctgactacagaaattggcaagagaacctgactagaacctggtgactgaacctggctggagaacctagcgagggaacatggctacagaacctggctggagaacctggagaacctagcaagagaacatggacacagaacctggctggagatcctagccagaacttcgctggagatccagaccagaacttggctggagatcctggatgggctgctgatcagctgaacgctgtctctgtgtctttccttcttcgccaactccgtccacacctttggggacccctggacctgctggggttggaccccggcacttctgctggttggctgaaagcttactatgacagagcaacgtgactcagaaaattctgaaataaacaaatagtagtcatttcctgggtcataggagaaagtttccaagaataaagaagaactgtagtcatgtgtggcatttgtattgcacaataattcatttaaaaatagttcaaatctcccaggaaagggagtagacactggggaagcaaagagcatccaaacagtcaaatagctgggccacttaactgctctgagcttctgattattctgtaagttgaggattctccgtcagaattgttttgagaattaaataaaataagataatgtaagctattcatacattctaaaccactatacaaatagtagcaattatttaaatcttttttcttctttccatgttcattgccctttacctctgttacagactttaacatattgcatgagagaaaagtttacgtccagtttcctgactggaatgtgattgactcaagaagaaaggctattacactttttatacaacccacaacccaaaacactgaccagcacatacaggtcagtgctttaataagtaccttttccctttaatttttttgaaaaccagaactaagtgatcctttctagtctaacatttatttccccaaattccctggaaattccacatagtgacttctgagatttatgagttgattgctttatcactctctgggcttgttggtgtaaatgaatttaaagtagcttctcgccctgaccaggtggctcagttggttggagcatcatcccatgcaccagaaggttgtaggttcaactcctagtgagggcacatacttgggttgcgggtttgatccctagtcagggtgtgtattagaggttgccaatcagtgtttctctctctcccccttcctctctttctaaaaatcattaaaactatgtcctcggggaataattaaaaaaacaaaacacataaaaaaaaaataaagcagcttctccaccaatatgtcagctgattttgatctgcaatgatctgttcaaagtgtaaagtttatgattagaccttccttaggctcctgtttgtcaaccgttttcatctgctaacttcacagtatcccctggcagtcaactgacccttcaattctccttctccactctgaacactacttcaaaaggtccccttctgtgggtgcctcccactaataatttatcttcaatttccagactttttttttttttttacttctccttcaccagtatctgacattggaatctggaaagagtggcagcgttaataataaaatgaagcagaatccttggacctgtcatcattaattgctttcaaagtctggcctgtctttatgagatttgtttctttatttcctctttgaaaacaaacctgggttattgaacatctgggctggaatttctgctgaagaaaaacattcagaggctcttgaagtctgtccaagctctcaagtctttttcgtaagtggatagaaccacttgcttgagaaaagaatcccccccccacccccttcccacttcagtgtttttttgtcatgtgtgctggtgagcaactttcattgtaggtatgagttacatatggggatagagagtgataggctaaaggatggacttttgtactgtattgagtgatttcatgggaagaggagtttataggtcttttaattttgaatcataagatcagatgagaatatatagttttttaaaaagccatgtatggatttgcaaattactcgtgattaagagagaaaaatcaggacatgttcaatctgtgttttaatcttcttattcatcacagtcactatttacttgtataaaacgtttaacaattattcctgttgcaggccactctcttaaagaaaactgccaccagttttgtgtttctagcatcagtgtatttaaaatctaacccagtaaaagtcacaacaggctgagcttcagggatagtttgttggccacatcttacatgcacctgtgttccactctggactcaaatgaaaggcctcccaatagggcatggaagccatctttccagggagccatagaaggggccagatgacacaaaggaagtgcagggagtaacttggatcaatgagaaacaggaggtgggtgagttccttctcatgaactgccgcaaagcacggccttcacagagcctgtccacaggtaccccgtgtggccgtgcagaggcacttgccacattacctgctgggtttctcctgtctctttgggaaggaagcacggcccccccccccccaccccccgccttacttctcttgttcttcactctcacattagcaagttaactctgttttccagggaacccaggataagatggctattggttaaaaagatctactagccagcacactgaaggttcacttcaatcctggataattttgtgaatctgctatgtgcaaggcactaggttaggtatcctaagggaagcaaaggtgaatgaataccattcaatccctgaccttaggtgccaatagtctagaacaaccgtgggcaaactacggcccacgggccagatccggcccgtttgaaatgaataaaattaaaaaaaaaaaaaaaagaccgtacccttttatgtaatgatgtttactttgaatttatattagttcacacaaacactccatccatgcttttgttccggccctccggtccagtttaagaacctattgtggccctcgagtcaaaaagtttgcccacccctggtctagaagataaaacaattgcaaagtcagcatcactctgaaagaggaattcaggactgttaggctcccctgagagcaatagaagccaatgatctgttacacagggaattgaggacctcatccacagtacattttataattattactctggctgctgtgtttaataagcacagaagaaactccatttatagtggctcaggtgggaaacaatgatggtggtgtgaactttgatagcaataacaagaatggaatcgagatagattgaagacatgtttagattggtaatgattatggatttgatattaggaatgagagagaagacagtatcaaggataactctatgctgtttgttttgaggatctgggaaaatgtgaactattcactgagatccctacacagaaagagttttggctaagaagataatgaggtcacttttagaaatgctaaatgttaggaatctgtaagctatctttataggcagatggataatatagacctagtgctcagaagagacctagcctaaagatctgaaatggaaattcatgagagtctcatttgatgtgtgaatgggagttgacaaggtggctggtggaagtgtggcctgaaaagagaacctaagaaacatcaagtaagtattaggtagaggaaaagaagcttgtaaagaaaaatccaaggaagtgcaaactgaaaagataaaaaaaaaaaagaaaggatgtgacagctctggaaatgtgccactctgatctccctccaagaaacagctttctgtttaactgtcaagaagacatttgctaccagcctctatctgtagcacatttgggatctgtcacatgtttgaaccaaggccctaatcttcctggtaagcccacagctaatggatgagcatggcaggggcactagggcctggctgttatgcccaatgtggggcccctctcacaagccgtcttcactaagagttccccactgggttggccaagactttgtccaatttgcagagaatgtgaggctgtatctaccaagtcctgcttcctctctctttatctttcacagatgttacctcctgcaaaccacatgcattcccaactctgtctcagcatttgctccagtataactgaacttggttaatgggtcatttgatttggtatgtagatgttttttcaatatgtccacctaaaatgaacatactgaatttggctgcatgtggccttgaaagacagcagcaagagaaaatcttcctaatgggaggatttcatggaaaagttgcctgagataagactacatatggaatcgtgagcaaaaggattgcccagctagtcaggggcctgggaagaaaaagttaataattttgcagagaaaaaaatctgggctagatgcttgtggatgaatatatgagagtgggtatgaaatgtgaagatcttcgtatcatatattaatgaccaccaaagcttggccaccacagaagaagcattgaatagccaagacaaaattgattagttgacagtagccatcttttgtcattggccgctccagaattagcatgattggcacattaacagaatgaacattacatgtcaatagcatggactcccctacttatcaaggctgatctactactactatcattgagagtccaaccagtttgcaagagagatcaatgtaacactattacttgaggggaccacatggccacttggtggcaagatgaccactttgagccttttcttactggaagggccagtgatttattgttatagtagtaggcacttattcttggcagcctagagcctcagctagcaccactatctgggagatttcaaagtatctgatctaggcatacagtctcatataacataccatctgatcatgagactcacttcacagagaaggagatagggaggtaagcccatggccatatcaagtgccacacaatccagaagcaattggttttacagagcaataaagtcacctgccaaaggcacatggaagtgctaacttgtaggtaatacttctgccatgataggataccatcatgcagattacagaatgatttacaaactatttggatgttgtgtcctgttagcaagaatatataggtctaggaaccaagggacagaagcagaagtggtcccagttgccattatggtcaatgacccaccacggaattctgtgcttcccatccccacaattctaggctctgcagggaactacaagctatgattcctacctggacattttggattccttatgttcagggacataggatggcagtaattgatgggagtaactgattctgatcaccaggaagaggtagtagacatttgttatataatgatggcagggggaaatatgtatagaattcaggaatccacctggccatgtcatagtacttgctcaattatgactaaaatagtcaattctagcaaccctagtgtgagaagagtatggttaccacaggcctagatctctgaggaatgagatcttggatcatatcagataagtcatcaaggccagaagggatgagagctggggaatttagaatagtcagtggaggatagagatggtaagtaccaattgtggtcctgagaccaactgctgctattgaaaatgaagatactagacacactataaaaagaaagcattgcatccccaagtcaccctacttttcttggtcatgagagacatactttctggaaagtgccaagtttcttatgaataaagggaattcataggtacgttagaatgtctattatctcacagaactgaatcagacagaacctgaacggatcctttcatgcttcatagagcacgtttcagagcttctggatttgctgttttagatcctggggatggaggctggttgtgaaatggcaacttggttatgctacctaagaagaccagccaagcattaagtagggacatctggagagtttggaaatactaacaggcttttcagccaaaatctgacattgaacagataaaatccttcttacagagttatggaccctattgggacagaaatacccaatagggtccagtcagaatttcggccactttacactacaaacagtaactccccagggataagtatctctggcagaagacaagaaatcctgtaaatctttaaaatgtcccttttgcatggaggacatgtgtttgctcttcttactaaggctgcacagccctcaagcctaacagctgctgaatggtccctgcgtgagagcttgacctctacaaactgatgtaatgcccacactgccgggtctgggtctcagacccactactgcattagtgataatagttaacttgtgatcaaagctttattctggtttaccaggtatttcatgtagagctgctagcattatctcgcgtgatcctcaaaaccactctgtgacgtagctgcagcagatctcattagtttgttcatcttgtgagttgagcaaatggtctcccacatagtgagttggtggagagctagaatttgaatgcacatcttcctagaaacccaagcccagagtgggctggaaagaaatagctgagaagtaagctctggtaaacagaaggggaagactactctatttctctcattgaaccttttttaaaaaaattcatagcaaattttctttattagcatgttgattaaattttttaatggagaaacacaataagacatagaagccaatgttattatgagaaaacccattacaaaaaggatactgtttgagttctgcagcttccctagaacagagaaaagtcctttgcataagtattggctctttaaagggtaaagagttaacggacaatgactcatgcattgctatgaagtcacctctttctagaagtgttgacgaagattgccatgggatgcaaaagggataagcggttggtgttcctcagccaactggctagaaacacttcactggccggggactcttaaaacctgtgacatggggaagcaatttcatggttctaagcctctgagaataagcatgggcattcttgctgttcagatggttagatataaggtgtaaggcatttagacatggagtgggactgctctctactccagtgtgtcagaaatgaatctgcctgaaacatcctgaaagtacctaaagcttggcctcttgacctgtgttttcccctcaggccacttgtctggctgcccacaccctttccatccctacagtgtgcatgaaacatgaaaacactgctaagtatgtaaaaaaagccaggcacagaggtcacatgtgatatgattccatttatataaactgtccagaatagacaaattcatagagacagaagcagattagtggttgcttgggcaatgggcagagaggggatggggagtgactggagagtcaagtgcagggtcttaccttccaccctggagcttttctgtcactatctgggaactctttcttcatttatacacaagaaacatggaaggttggtgcaattaacaccctgacccaggggactagaagctgatggatgtaatttccctttcctccctcctagagcaagggtcctgagatgcatagcatcggaattcacagggagtcccatgggatcagcaactatttgcttcgaagcagcagccgactctccctccctgccaatgtcattcccccatcctaccgttctgctcccccaaatcacacaccccagtaagaattgtcacacaggctctgatttctggagaacccaggctaagacattaactcttgatttcactcagactccacacagtctgttaggaaatcctgttggctttgccttcatcataaaatccaagattcaaccatatttcacctccttccccgctgttaagctggtccaagccattgtcatcacttactgattaacaaacacagttgccagggtgatcctaataaaatgtaaatcagttcatgtgattttttttgaactaaatatcctattaagcttcttatcccatttagaacatatacatatttaatgtcttatagtgtctacaaagtacatcatctgtcttataacttctctgatctcattttccataaccctctctggttcagccacatttctctctttattacttcttaaatatacacatttctcttcttgttacttctaaaatatacatgcatgctcctcccatgtactcctacctcgaggtctttgcatttattcagaatttctctcccaagattgatattcaaatgtcacctgctcagtgaggccttctcggaccatttaaaatgcaacttccagcattccttctccctctcccttaatttttctccattatacttatcaccatctgacattctatgtgtattattattattagtattattattatctggctatttccccctgctagaatgtcagatcaatgaggtcatgtattagatatcaccaatgcctaaaatagtacctgtgtatataggtgatcaataaatatttattgaaaaaaaagaagaatggtaggaaggaagaaaagctgcaataaaatgccacctgtccagaggtaaggtttctggtcatctcaaatccttggatggccacgtagccagaaataactgagaagagtgttgagaaggaaacagagcagtcacaaagcccactcaacactgttatgttcattgttttgtacagatcatttccctgaggtcaggtttctcctttcccgttcacaaaagatggaaataacaaaatacggagtgacagtagctggcagtaagaagagacgctttcaatccctgctctggagtaattacaataacactcaggagtgggtaaaattaagctaatgcctccaagatgtaagcccaaggcttcacatatcatataaatacagggatcaggtttttcaaagaaaaactgggacacagggaactcttcccaggtgcaagaggggtgtggaaaattgagttaagatgatgatgaaatatttgagtttctggactatcctgagggtaaaatattttaaatagggctggagcagaaaaatccagggcatatttacatgtcactacatgtaagagaacaacatttataagatggtttaaaaatcacgactgctgctttgtcattttaagaaaaaaatttttgaaacaccctgaatgttcttaaggattcaaccaaagtcaagcaaagtaagtagagtcttatacctaaatggaaccctttaaaaatccaatttgcagctttgaacacatgggaccatttcccccatacatttccttcaactgtgttttaattgcaccaatgaaaggcaaaatggtaattaagcatcccacccagcgtcagtgaaaatggaaacagggccagtatggagcatttaatgtcagatacgaggcatcagaaacttgtaaaatctgtttgtaaactactttaaaatacatttagtgaacttctgcatgaagagcctatttaagggctaagtgcaagacatcacgaggtaatacagtgtttcagagaaaaaatatcacatccgaaggtcaggcaaatttttatttatttaaaaatatgttttagtgttctttatcctatataataaaagcctaatatgctaagtgtctggtcatccagtttgccattcaaccaatcaaaacgtaatatgctagtgatatgctaaggctgctcaaccgcttgctatgacatgcactgaccaccaggggagcagactctttgactggtaggttagctggctgctggggtctggcggatcaggactgagagagatgggaccgacacgccctggagccctcccagggtccctccctggctggccaacctcccgcatccctccctgaccccgattgtgcactggttgggtccctcagcctggcctgcactgtcttgcaatcccagctgagaggacctcccaccccgagtgcacaaatttcatgcactgggcctctagtaataaaaggactgagtcacatttaataccaaatatacagttattaaacaaaaatgccaaattctcccaaaacacaatttcttgacttttttgggggtggggttaagggcacctctggaaatctgagcaaggttatagattcagcagtttacttgtattcaaaatttcatttcatcacagattcccaggaggactgaatggacccccacgaagcccatgaacctccaggtaattaatcttcatctgaagggtaataaacatactttgtacaacaaggaaattttaaaataatcaaatagatcaattcaagttaacaatttctcagctctggagggaaaaaactttttcaattgattgataaagcagataaaattcactcacttaataatatgttttatgtgtattactatgctatcagagctgaagtaagaaaatagatgccattttgctgcctgtaaagtcttgtatgatttttatggtttcccatcttacatgaaaatcctttagccattttgagtttgtttttgtgtatagtgtaagttggtcatctagtttctttcttttcttttttttttttgcatgtatctgaccaaatttcccaacaccatttattgaagagactgtcttgactccattgtatgctcttgcctctcttgtcaaatattaattgagcataataatggcttgggtaaatttctgggttctctgttctgttccactggtctgtatgtctgttcttgtgccagtaccaggcagttttgagaaccgtggcttggtaatatagcttgatatctggttttgtgatccctccaactttgtttttctttctcaggattactgtggctatttggtgtctttttttattccagatgaatttttggagagtttgttctagacctttgaaatattttagtggggattgcattgaatctgtagcaattcctttaccgacacagctcctagggcaatggagactaaggagaatataaacaaataggactacatcaaaataaaaagcttttgcacagcaaaacaaactatcaacaaaacaacaagaaagcccactgcatgggagaacatatttgccaatgctatttcagataagggtttaatctccaaaatttacaaagaattcatacaacttaacaaaaggaattcatataacttaacaaaaggaagataaacaatccaatcaaaaaatgggcaaaggacctaaatagacactttttgaaagaggacatacagaagggaatccactgtgacttttcaaaactttgtggggcagtagagttctattgccccacaaaattttcaaaaatcagtttgctactgttgtagatggttgtagtcattatttaaattagcataccactttacccatacttttaaatgttctcttctctatctcagaataacttaattcttcagaatagtcccaccgccttcctccccctaattgtttttaaggtcagattctccttgctcctgtctctatccagttcttagcctttgtggctcagtggattccctcttgagcatcacagacatcagaaagcacctaactaaaccctctgttttacccagaaggaaacttgggttcttctaaatccatgatctcacctaatttcaccatcagtgagaggggggtaaagagaatgcaggtttcctgactcctcctccatcgctctttctctgatgttaccctttctcctcatcatttggtggaattgacaggatgaatccaaattcacttccactgtcagagcctcatttctcttcaatgacttggaattctgccacctttaaagaagtccccgaagaagtgatctctctggccttctttatcctcacagcgaataatgagactagcatatttatgttcggtcattcactcacttaaacttttttttgaacatcaaagtagttctaaggttgaatgacagggtttccgtgttagaggagctcgtgtcctacttggtaggcatgcaaagaaataagttacagaaaagcatgccatgtttaaggactaggatagaggttgtgcacagaatgttggggggggggggggagctcacccagtctaggggtcatagatcaagagaggcttctggaaagaacggacccctaatctgtatctttaaaaaatatatgcttatattgatttcagagaaaggaagggggagggagagagagatagaaagatcaatgatgagagaatcatccatcagcagcttcctgcactcccctcagagcacagaacctgtactgcaacattaatcaccaaacaacattcaccactatttaatgtccctatcagctacagtcctgagaaaatcctgaagtgcaacttcatcattgcttttcatgtatctacattaataacttgtattatcttacctgaactaggtttaaccagcctaggttgttttactatttttttttttctttgtgctccatttttaaccacccaagaaagcctttttcctggcttttgtctgtaattatcctttcaatggtcattgatttaaatgcctcttaaggaataagcatacaatagaacttgttgacccttgtttttaatttggcaacttttaaaaataaattttatctattacagttgacatacatattacatgaggttcaggtgtacagcatagtgattagacatttatataacttactagaagcccggtgcacaaaatgcatgcatgggtagggtccctaggccaggccagtgatcaggaccaatctgtggggcgactggtggggtgatcgggggggcccccgctggcacccatcttggccggcctggtgccgctcaccagctagccccatcccctgatcgccccatcccctgatcaccccgtcgctgcagccagttgcctccctctgcagggtgacccacaggatgatcggggggatcccccactggcacctgcctaggctggcttgggacctgcaggctgggggcagctcctgcgttgagcgtctgcccccctgtggtcagtgcacatcatagcaaccggtcgtcccgctggtcattctgctgtttggtcaatttgcatattaggcttttattatataggatgaagtgataacaccccaataaggctaataaccatctgacaatttttggcatcatttaatatgcagataaaatgcattttatgaaattttacctctttattcctctcccatgcctatggcaaaggaaacttctggatccactatgttagtgggcagttaacctttgagaaagtaactcaggtacttcatacaaagttaccttcatggaaaatataaagggtatattttaaaaaagaatacctgactaaattttttcaactatttgatgaaccagcaagaccagctaagataatttcccattttggaagactaaaaggagaaatacattttttttgcaaatttttattacccaatatcaaatagaaaagtcagttcttgaattgctactaatcttgttcatacctgaggtaaaaggaagatttagcctgggcgaggagcagggaatatgcagaatagctcatcttgtaccaaagttttgttaccgtttcctgtgttattaaataactatgaggtcagcaatcccaaatttatgggattgagtttttattcattgcgcgcgcgtgcgcgcgcacacacacacacacacacacacacacacacacacacacaccctttggcctgttctcaagttctttttgtgcacttccgtccacctttccggctgcccttggcctctgggacctaccctctcctctcccaggtactcacagagtctgtgtcccgggcggtaggactgggatgcgatggatatccttgcaggtgactctgaagttgtcctcctggtggcactcacagggcggagaagcacaccctttcccccctaagctcctgggcaggccgagaagcagcgccagctgcagtaggggcatcggcctcattttccaggggtccgaggccatttcttcttcctccactatcgcctcattctcagccctctgtgctccgagctccccgcgcggaaggagcgagggaggagtgggagaaagaggactcctggtggttgtgacttcagcaccaaagaggaggctggagagggcaaagaggctctcaagtgttctgtttccctaggcaactctccacgcaaacaagcccattggctttactttcctagccctgagccacgcagcccgctctgtctgttccctgacacttaacaccagtccctgctgtccccaccctgctttccagccagcctcgcgcgttccttccattagctctttccagagcctagttcagaagcgccccccccccccgccgccctttccttgtcgttaagtccagtgggccaccttttcccgtgggtagattccgaaaagcgctgcctgtcaggacccctctcccccgctcccccttctcccccaaggctcgcacgcactgcagagaaaacttcccgcaaaatgcaccgtccagggcagacgcgggtcaagatccttctggcctaaggtgcctctgaggatgccagcctgaggtgcggagaagccggggcgcaatacctctggcggaccctcactgctccccaacttcagcgcaaacccgaccgctttgcacagaaattgtattgaggagtcacctggaggcgtcagctcacattgcagcaggcggacttggctctaggcgggccgggaactaggctctggggccatctgctgggcagaaccaggaatggcaggaacaggtctggtgctctggagccggccttccgacctggcttgaagcctaaaggtaaattctggactgccaggaatgttattgcatcatatagctccgctttccccaaaatccccccaaagatgggaagaaaaataacacttaagtggcacaccatcatctacttaccattcttcacacagtagccaaaaataataatttaaaatataaatcagattgtgccatttccttcttaaaagccttcagttactaattcttcttgggatgaaaagaacatcaattcttcacgtactgtgtaattccatttatatgactttctgggacagaatggggggaactatagtgatggaaatcagaacacgatttaggggagactgattcaaaaaagataggagagcactcaccagggttatagaaatgttctatatcttaattagactggaggtaagagaatgcataatttttttcaaagttcattggactgtatgtttaaaatttgtatatttcaatggttttagggctttattcataattaacccaaactggaaacaacccaagtgttcctcaactaggcaatgggtaaacaaacagtgatacatccatacagtggaatattatccaatataatcatcctatataataaaagggtaatgtgcaaattgaccctaatggcaaaatgaccagaatgaccgctggaccagtcactatgaggtgcactgaccacctcttggtccctttccctggccagcaggctctgatggcctgatggccgtggcaggggtggggctggcgagtgggcaggaacagccaacctctcggtcccttcccctggccgtcaggcaccgatcattcaatggcaaacagggaaccgggggttggtggcggcgggggcagggccggctgtgggcagctggggaagatggccttgatcgcaggccaggacaaggaaccgtatgtatgcacgaatttcatgcaccgggcctctagtttttcaataaaaaaacaaaacaaaaaacagccctggtcagtgtggctgagttggttggagcatcattcagtacactgaaaggtggatttgattccaagtcagggcacatatctatgttgtgggttcaatccctggttggggcacatgtaggaggcaactgattgatgtttctctctcacatcaacgtttctctctctctctctctctctctctctctctctctctatatatatatatatatatatatatatataaaatctcactccccccttcctctctctaaaatcaattaaaaaaatacataatcttgcttttactttgaacacttttgttaacatttcttgtggtacagatctgttgggtatgcattattccatcttttttgtgtctgaacaagtattaattttgccttcatttttcaaggatagcattactaggcacagagttatagattgaaagtatttttttcctttcagtactttaaaaatattgctccactatcttattatttccaaggacaaatgttctgtcattattttcttttgcaacatgcttttttctccctcagaatgcttttaagactagttatcgctggttggaagtattttgattatcatatgctttgatgtactttcttcatatttcttctgcagtattttgaggttcttgaatttttgcatttatgtttgtttttaaaaatgttttattgattttagagagagaggaagggaggtggagagagagagagaaacatcaatgagagggaaatatccatcagcttcctcctgcacaacctctattggggactgaatccacaacctgggtatgtgccctgactgggaatcgagcccactatcttttggtgtactggacgatgttccaaccaactcaactatactgctagggtgaatttctgtttttcgtcacttttggaaaaatttctgccattatttctccatttatttttatgtccccttatgtctctctactctttccaattattacacgtatattaggctaattaaagttgttctacaaaaaagctcatttatgctttttttgtttttgttttttgtttcattgatgcttttttaaattacaaattaattttctctctttttcattttgcatagtttcttttttatatcttcaaattcaatattctttttttttttttttctgaaatagttcatctgccactgattccatcagtatattttttatgtcaggtgaaatgttcatgtctagaagttcaatttgggcctttttcatatttcccctgtctctatgtaacttttgacatctaaaacacttttgatgtctgtgttagactcctagagctactgtaacaaagtaccacatacttgatggctgaaagcaacagaaatgaattcactcacagttctggtggctggaagtccaaaattaagctgttggcagggtcatggctccctctgaaggttctaagggcgattccttccctgcctcttcccagatcctgatggttgccagcagtccttggtgcccttggttcgcatcaccccagtctctcctctgtcatcacacagtattctctcttgtgtctctgctctaaattccctctttccatgaaagcatcagtcactggattagggcccatattaattcagcatgacctcaatttaacttgcatatctgcaaacaccctgtttccaaataaagtcacaatccaggtattagaacttcaacttatcttttcaggggacacaattctacccctaacaatgtccctgtttgataactctaatatctttatcagttctgggtaagtttcaattggttgatttttctccttaaaatgagttgcatttttcttcttctttgcatgactggtaatctttatttggattgtagacattgtaaatttcacctctttgggtactggattatttgttttcatagaaatatttttgagatttgttatcagacgtgggaaaattacttggaaatagtttgatatttttgggtcttgctttcaaggcttgtgaggtgggatgggaatagcatttaatttagcaccaattactcctccctcctgaggcaagactcttctgagcactctacataatactacgtgaattaagaggtgtctccatccggctggtggaaataatcatcattcctggccctgtgtgagtgacaggtactttctcctccccagtcttggtagatgcttccttctccagcctccgtttcttcacgtattcatctgaatattccaagaagacgttccctagttctttagagttttctctgtatgtagctctcttctctttgctactttgtccatctacactggttttcacaaactcccatctctgtcttgtcaactcaaagcatccacggagctctgcctgggttccttctccttgagccatggcctggaaattctctcaaggtaggaagttgtggggatcacaaggctcatgttgtttctttcctgtctgttaggaaatcactggcctttatataatgtctattgtctgaaaaagcattgtttcatatgttttgtatgtttcctttttcttttttaaaaatatttttttgttgatttcagagaggaagggagagggagagagatagaaacatcaatgatgagagaatccttgattggctgtctcctgcatgctccccattagggatcaagcccaaaacccaggcatgtgccctgattggaattgaaccttgacttcctggttcataggtcaactggctcaaccatgctggctagtggtgtctttttactttttaaaaattctgggtgcagggaaagtctgctctctgcttctatgctattccatctttgctgggaacataagtggtatgccaaaaatctgtaacattttaatcaagacttaggcaacatgaaattggctcagtggtggtcataatgtaggtaaaatcttccatgcttacctgtttcatttcaggaaagtggaccaaatcttcagggtaaccatctctatatataaaaagccagtaactgaaacactgtaactaccagaacgaccagccgctatgacgcccactgtggccagagaaccggcctgatctgggggtggggccagctgaccaatctcctgcagcccctccccctggccaaccccacccctgatttccccctctaccctaatagggggcagggccagccagccaactgcatgaagccgctccccctggctggcctcacccccaatgggtccccaccaccccaatcggcctgaggcccctccccccagccagctctgtccccgatcacaccccacacaccaatcgggggtggggcttgcagaccaaactcccatggcccctcccctggccaatgacccctgattggcccccccccccccccccccaatcgggggcagggcctgctggccaatcacccacggtccctccccccagctggcccgaccctgatctgcctgattggggcagactggctggccaacctcctgccatctcctcttccctacaggcctggacccaatctgcccctattggggctgggctggctggaccccacctgtacacaaattcatgcactgggcctctaatattttaataaccttctcaacaaattgctatgataatcagtaaggcacttttctttccttagcctaactttgaaactggggccaataaaagcatgggatgagaagatcaaccagtttaagtaagagagaagtgtttacattctgatccattaaccacaaaaatgcttttaaatataacagatgaatgagaatggcttaagaacattacattttccttcttgtagaagaattgttagcatggtccagatatcagaaagtagtatttagcatggcattactgtctagtgttttatccaccccccccccccccccacacacacaaaggtacctattaaagtggagcacaggctatctgagttagctataataaaaatattttatgattttggatttgtgcattaatagttattcattttgacttaaatgcaaacacattgctgtgttttacatttgtgatgctcttgacagtgacatatacctataattaataagtgtgtggctgaaaagaatattatgggatcaatggcacacatgaatagctgaagatgattccttgagtagttggcaataccatctttgcaacatcttgtcttctattttaactttcccagattttcagttagcttgaaaagttgatgcaggaagagttaacaatcattaatgatcattttgtgctttggttcagagagaattttgggtagttccaattaggatgttggagatttctctttaaccttacccatagcatcccacttcggagctctgtcttccttcttatcccttataattgtctggtgatattacaaggaggctcgaccatctagtggctagatggaacaatgtggaatggaacagatggttactacttgtgctggacaaccactgggtaaaatgcaggctgtcttctcttcacgccaatgcgcggacagcgaggtccagactatgaatagtacttgcttgagtagaagctgagcagtctatggacatgggctttcctccaccaccttcaggtctgcattgcacccccagccccgaaatctccaactctccattctgggactcttttctctgtttcttatcacattcccactcaggaatcccttttctcttggcactgctcttcttcatattttcttgactctgagcaacatgagggcagaaactaggttgttattattctcattttcattacatttcaagtgtctaatataatgtctggatcctagaaagattcagcatacatttgctgaattaaatgaatgagtacaattgttttttttagtacaattgttttgaaaaattaaaatgtaaacaaataaaggtactctaatgcttctgcatgtgcatggggcaatttaatgactgttcctatttatgttaacatttgtgtccggtgcagtgtattctcttttcttttcttctttcttttcttttctttttcttttcttttcttttcttttcttttcttttcttcttttctttttttcttttcttttcttttcttttcttttcttttctttttctttcttttttctttcttttctttccttttcctttccttttcccttccttcccttcctttcctttcctctttcttttcttccttcccttcccttcccttcccttcccttcccttcccttcccttccttcccttcccttcccttcccttcccttcccttcccttcccttcccttcccttcccttttctttttttcttttctttcttaatatggaacgcttcatgcatttgcatgtcatcctcgagcaggggccatgctaatcttctctgtatcattccaattgtagtatatgtgctgccgaagcgagcatgagtgtattcttgaatacagctgtgagcaagcaacgtgaagacagatgcaaatacataggtatgtccttagaataatcatcttatgttcaaaccggttataaaataacatctacctcacagtattggccagccaggaaaaatattttaatgcttctgatatatcaaatatgatataatagatcataattattaatgtgtgtgtgtgagagggagggagagagagagagagagaggcagagagagagagagagagagagagacagagagagagagaagggaggatggaaattttaacaggttaggtccctgtgtttgaaaggcttagaatctaatgtaagtaccaacagggaatggataattacagcacacaatgattaggatagggaatccatgtcgttaaagtgatcttgctttatctgtgtaggaaggcttcacagataaagtaaaacttgactttagccttgaggaggaaagaaagcatcaacagcagaaggaaaggaccatcttcgtgacgaatttggattttcacttcttcttagcttttccaaagctcagttccttgttgaagtttttaaagttactttttcgcaattgaaaccgacttctctgagcggactatgtctctttcaggttccacaataagcggagagatgacaccctggcagcagttttgaagtgacattctcaaattccacttttaatttgtcatcagtgagcttactttttccaggcctaaaaggtgtgactagataatcttctgtttccttccggatctaagaatttgtgaaccacaagacattgggcactggtattattacacgaaggggattttaacatttctctgctggaaattacacctgcctgttccccaaggaataacaaaaatatagtctatgttagtattgcaactactttgattttttttttaagtaaaggaaatatcttagactttcctgattaaattacagtctgggccaaggcagctctctctctagctggtgtttgtcttcaatatacgcttctggccctttaaaaaaacaaaacagggaaagagggtgcgtcctagggcctcaggtgccgacgaatcctgcttcttccgggaaagtcaccgaagacaatcgtttggcatcctgggacttgcagtttcttttaggggcagagtcccttagttaactcattgctttagaggacgatgcggactacaactcccagagggccaagcggcaccagtagggtcctgtactcttgccgcggctgggcggctagttcaagttgccatagctgccagtctgggcacagctcagttcgttgcgccaccgcgtcgggtcagttgcaccttctcggtcacaggtggccgtgggagccgggtggggcctcggcgatgatccggcattaaggacctgggctctcaatctcaggtaattgcccgccaccgcccgtggccctttctttgccacttccctgtcctgtagacgttttaagtcccacttctccttgtcccaggggtggggacactggccgcgcggagatgtggggtggggccgcggcttggtagtgtgagcgccagggtggggagagacggggtgcccggcctgcgaccccctccgcatgggcatccacagcctggtctgtgccaagctattcagggaatttgaacttcctggtcctccttttgcgagcggcatgtcctgagctgtcaaggagatcggactcattttggaggtgactcgcctggcctgggatagagctgggatagcgctggccgttttgtacctgttctctagtcacagacgctgtatcgctgtatcggcctctcctgtgccaccatgggaattgccccaacgttagctgcccgtgatatgagttatttcgtaggaacatctgacattgaattccagctgcaactgtgatcggtttgatggggcccagtggcttccgcaggctcctcgtaggaccgcactgtaacgttgttccaagtatgaagcaaagtatctagcacagggcctgacacattggaagtactcagtatttgttatttatatttactaattatgcagaacagtgccgtcaagtagaaatattatgtgagacacacacacacacacacacacacacacacacacacacacacacacgtttaaattttccagtatccacattgtaccagtaaaaaggaaactggtgatgttaattttaatagtttatttagcacaatttatttatatacaaaatgttttcatttcaacatgtaatcagtatggaaacttattgattgaggcattttacagccttaaatgtatatatgtttttaatgtttctaagtgtccaaaatgtgcattgtatatttacagcacctgtgagctcaaactgggcagatttcaagagctcagtagccacgtatggagtgactactgtattgcacagtgtagattcaaggtgaaatctttagcagcttcatattggcacagagaaattgggatgtgtctgtgacccagcccatggggtgaaggttatgctctgttagaagagtgaagtatgcataaaacctgagtgctctcagcgtgtccgctctgccattgcaagacctgtgttgtttggagctctccacagtttattttgattaacttttccaggatagttttttactttttcttgaagtggtttatagatgtagtctcaggtcttcttatgtgacaactactagttctctggatacttcctgtgatttcctacagactgctttaaaagcagtcacaaagaattgaaaatatttctgtccacagactttcttacaaccaaagacattttatggttacttatttaattgcagtactaaatgtgtgtggtactgggtagcttacgaaatatatcgtttaggcgtactagcatatatagttcgttaactggcttagaaattcaactgggtcattaaattcagaagtatcatgttaaattaattatttgacttagtttaccctctgatagacatgctgaaagagagagtccttttaggtggacttgggctacaattgcatgatctggccagccctagggatgaccaggcatgggcagaagtaacaagtgaggagtatgttggttgggcacttggatgtggagagtttttttgagttcagtggggcttttgcatgacccaaggagcattttctctgttctctgatatgccatgaactgtatatagtgagggtgagcttattttcatcttttatttctttaacaaatttaaagttgtaggtgtgttttaatccttgtttcatttatctgtaggtggtttggggaaaatatagaggcaaccaaagattatctacttgactaggccttttgccctgaacctcatgaagaaatgataggaggcagacatatgtgcctaagaagagcactgagctcagtggagagcaacacggcgatttgggggtgtgctttgtaagtatgctttctcccggtgtttctaaacacatttttgcatatttatagtggaaaagagtatgtcaacaaataagcttttgagaggcaagaagtccctttctaacagcgtaaaatttaaatagaatggttgtgcgtgagagttataaaaattcaagggaacagactagataccttaacaaaaacaaaaatgatgcaattgtttccattagacaaacatttgagaacctattttgccccaggcacagtaatagaaaccttctatgtgacatatatatataatatttttattacactattctttaaaaaaaatatttttattgattttagagaggaagggagagggagagagagagaaacatcaatgatgagagagaatcatggattagttgcctcctgcacaccccctactggggattgagcccacaacctggtcatgtgccctaaccgggaatcgaactgtgacctcttggttcatagtagaggcttaatcgctgagccatgctggctggacattattatgctactagagttccgatacatgaaatttgtacaagagtaggcctttcttcccttggctgccggcaccaggacctgggcttccttcgcagagggaggccctgtccacctgccacagcccaccagtcggtggccttggcctccctctttggggtgatgtggagccccctattgatcgccccgcaggccggtggcctgggtcttcctccacggggcgataatggggcgatggccgggcctccgaccaatcgcatcgcacctaccttggctggcctggcgccagtgggtgtcatagcatggtcccggatggttgtccggatggtcattccgctgttcagacgtttggtcgattttcatattatgcttttattatataggattctttgttttataaatgagacagttgaagctcagatggtacatgtgatttgctggggattctggagcacaaagagtggttttcattcccgggtgttctttctacaaatctcaatcctgcccagtatgaagaggtgccgttgctgtctaagattgccctcagtgggatttgattgacagtcagcatttttagattaagagtaccatataagcaaaatttactgttatttatcttttatttatttattcagaaaaatatatttttattatagatttcagagaggaagggagatgggctatttaatctctttaaatattcaaagggttctattatcttagtttttccaagatataatggttcttacattattattgaggtaaagagttataattcgtgattactgtaaaatgaacatattcagcatgaattgtcattgataatgaaggatgtagtattaagaactcaagttttcagatgtagaggtaaaactctggacactttttagcattttcctcttgcattttgataagggacacttaaacctatataggagaccctttagggccttatgaaaacactgtccaagatctatgctgcacaaagttcttagaattaaagaattttagagctggaagagccagataaggtcatataatcacaaattctcagatgtaatcttgacatatttgggaataatataaaaataatctagagtagattcagtggtatccatgtagtatcactgagctgattctaatatgaaatgtttcaatatattcatgttattctgtagtgtgtatttaattgttgcaataatgtagtattacatgtcagtacataaaattctataatgtctcttttaactgctacataatattccaaatgactatatttaatgtatcccttatcagtagatttgtaattctttcttattttcttttttactattttaattctgcaattaaactctttctacatatacacatgaaatgtatacatacatcttaaacatcttttgttttcctaggatactatacagggtggggcaaaagtaggtttatagttgctcgtatggaaaacaatacaataattaataaataataatataagaataaactctgttgtgcatacttacaactataaacctacttttgccctaccctgtacatagaattagactatctgggtcaaaacatacgcatattatttattttttacttgtcttgccagattgttctccagaaaagctatctatttcagttaatattcttacagcaaggtgtgcccttgtcaatcctggaggtgagagctattttttttaaatttgtacaaaacagataaacaagcatatcatattattattttaattttcaattttttgattaggaatataaagaacattctatacttgttagggccatttgtattctgtgaattgcctgtttatataccctttgcatcttttagtactggtttgttcttttaaaatcattttttcagagagatgatgaatcctttggtttcgatataatgcgaatatttttccgtttgtccttctttccccctaaccttgcttctagtatctttttatcaataattgtaaatgctcctaacaggctacagggtccttatcagatattcactcctccctgaaattcacgtaacccttaacattgatcagtacaatgctgttgtgtgggtgtttgggagatggacatatatgtgtttatatgtatccatttagctgaagaacatagaatctacaatgtgcacagggtttgccatgacacaggtgctcaggatgaatggataagtgaagaagagagacattgtggaagaagtgaacgatgattgttattttaaagaaagaagagataagcagttttaaaaaattatttgtttatttttaagtacagcttacattcagtattattttgtatgagtttcaggtgtatagcatagtgactagaccaggggtggacaaactttttgactcgagggccacaatgggttcttaaactgtaccagagtgccggaacaaaagcatggatggagtgtttgtgtgaactaatataaattcaaagtaaacatcattacataaaagggtacggtctttttttttattttattttagttttattcatttcaaactggccggatctggcccgcgggctgtagtttgcccacggctggactagacaatcatatactttacaatgttccccctgatatttctagtactcacctgtcaccaatcatagttattacaatattttggactatatttcctatgctatactttatattcccttaactgttttgtgtaactagcaatctgtacttttcaatcccttcacctttttcacccagtctcccagccccctctcctctggcaaccattagtctgttctctgtgtctatgagtctgtttctgctttgtttgcttgtttgtattgttctttagattccatgtataagtgaaatcatatggtatttgtctttctctggctgactatttcacttaccataataccctctaggggagataaagtaattttatagattgtagaaaaaaggaaaaggctccatggtgagagactgaaggctttaaattggcaagaaatgggtttgcttatgagcagcaattggtacattgaaagcatactgagaatggacctagatgataaagaacctaagtaatttggagtgtagtgctttatagtacagacaaataataattgaggagtgcagatagctgtagatgacaataatgtatttaaatgatttgcaataattcaaacaagcttaagtgtgattcatttgttttaacgcatttattcaaatggatttaaagtgtttatctttgttttctctataaaggacatattatcaagtgttcatatttttactattttgcagtttaatccatattgtttggagatgtaagataattaatctaggccatttctacaagaatcattttacccaacaattatttaattagtatttactaattgcaaggtcctgtactttgaaattcatctctggattaatcagacaggataataagaaaatgtaattgtaaattcatttaaaagtgttccaactacatgccaggaaatatgctaggtgttgaaatacaaccatggattacataaatagtctctggtaagatttccagtctagtggggagagagacaagtaaatagccagttagcaatctctgggcagtgttgagttggtatccactatggatacttggtgggctaggggtggggggagtgaaagaagactttctaggggaagtttcctttgagctgagacttgaaagaaaagtagaagttagcaaagggaaaacgaagtgaaaagaatgttccaagtggtgaatccagcataggtttgatgacaaaatgaacatgacatattctgtgacagccctggccagatggctcagttgattgcagcatcatcccttacaccaaaaaagttttgggtttgattccccatcagggcacatacctaggttgtggcttgatccctggtctggacacctgtgggggcaactgattgatgtttctctctcacattggtgtttctctctctctcccttcctctctctcagaccaataaaacatatccttaggtgagcataaagaaaaatttaaaaaaattctgtggcacaaataattgtatttcatgcagatattttgacatttcataaatgttttgctataattatagagctatacaaatgcagaggaaggaaagatcacttctggctagaggtggacagcttcattgagattgtgttatttgagctatgcctggacaagtgagtgtttaaataggtggcaatactagatagcagaaggcattccacccagacctaatgacattaatgaggggtggcagtgggggaagaatagttgttgtaatagctgaatgtaggccttaagttttctgtggttcttttattaagacaggtggcttgccgagggccttgaattcctgattcaagaacctagattcagtcaggtggccagtggagttggaatgtataggaaaagactcaaaagtcaaagaaaccattttggtggctagaaatgagttgagaaaaagtcttaaatttatgagtagtgcacgtggcagtgagcacggttagatggtcagaaatctggaccttggagaaagattcaagactcagctctgctacttagtagctgtgtagccttgagtacctttcttagtgcctttaagcctcttaggtatatatcatgacagtgtactgaccttgtagggttattgtcaaacaacaatctgatattgaaacattagttattcagcattattcaacttgagtttcttgaaggtttttgaaattgaggttcttaacttcagtcctcatccatttactcccaataattcaatttcattataattggatttctatctaatccagtatctttaaaaattttttatgtacatttaaaaacattccaggaaacagttcatggggaaatcaaagttaagaacttcaattagctttttttcttaaaatttagaagtctatacattttgaaaagaatatgccaattattctgtttatccactttgaaaattaagaacaagaaataaaattgaggatattatgttaagcaaaataaagtcagccagaaaacctaagaaccatatgatttcactcatatatgggatataaaactgaaactcatagtttccagaggagagggggtggaggattagtaaaggataaatgggaccaaatatatggtcacagaaaatggtttgactttggatggtgggcacacagtgcaatatatagatcatgtatcatagaaatgtacccttgaaacctatataatcttattaaccaatgccaccccaataaatttaataatgataaaaaagaaatacaaaactcttgaaatgtatttcagtcctttattatagaagagttaaaaggaaagatctaaatatataaagtaagtgtggggcaaggaagaaatatgggaattgagaagagtggaggaaaaagagaaatagaaaatacgtgtacaccaaaaagaaaatacaatttaaaaaataaaatttaagattttatgatgtttctgaaaaaaaaaaagaaaaaaaaagaatgcttttctgcctgaagaggtaaatattcttaacaagagaaagggttatttaaaagacattttaaaatagttgtagaaagttggatttataggatatttacatattatcttgtattcattatagaaatctctagttttaagtagataagcacatatggatgacatactctggccttaaatgcatttcttgggggagtgagggtgataagaaatatgatagggatggacaaagaaaatatttcagtgattgggcaatgtatttagaatttacacagttctgctcctagaaatgtagtggaaaagcactggcttagagacacaagatatgggttctagttctgttactgtaactttttatcaaagcaatggtgggcaagtattttattcaaagctttgtttttcttagcagtaaaataaagagttgcaattatacccactttgtaggattgaatggaagtcaaaagagaaagtgctaaagaaaatgctttgtaagctataaagcaacatagaaatttgaggtataaatagcagtcagaaatatttccctttcttaaaatttttagtggtgttttttcactaggaaagtcaagaggagtttctgtttctactttaaaaaataggaggaaaaaaaccaagtattatgctttataatatatttaagaaacagcttgtttcttatgaaaatcttaatgttatggggagtagattaaatagccccaattgctatttatgattgcctagtttttctctctcaataattgctagtctgtcacttgtcatatcaattgcttatgcttattttgtgcatcctaatttgaatttgccctccccactccctcctcaccaccagaagacagtctgaaaacaacttgttttctcaaggagagtagtctttagaggaaacttcaggaaattaaattttttctgaaccaataaaatgatttattagtaaactacatattgtcaacagatgatcattggtgccaattcaacgccaacaggaatccccaatttggggtgtggtgaagaaggaaaatttattcagtgcaaaacagctcaattgtgataaagcatagccgtgaaaacagcaagccgatgaggtggacctggggccaggcccctctctgctagagatctctgttccttgctggtctgttcagctatgtgctggtctgctcttctccattctgtgtttggtgcctgtcttctctggcccatgttggtctgttcagctctgctccagtgagacagagatatctttggtgttttatctcagccagggatcagtctttggtggaaaggaaaggtgcactctcagagccagaggggagctggcttatatagacagaagcccttgcctgtgaccttgattggtctgtcctcatgcaaacaaggactccagatccttgtagttggatttgtccaaaaggcactttcctgattggtcagaaaagagctgctctgattggtcagtgaagatgttgatggggataaagctgtaagctctgtttggatgggggaagcctcagtcccattggttgaaataggattccagcacttcctttataagggctggctcagcagacagaaaaacagtgcaggcaggcagttcagtgcaggcttctgtgtgaagtgcagtttgtgtggaaatggctgctagtctctgtttttaaaatttgagcccagttagccaccaggagaccttcttggtaggtatcttctttcgcagtgttcacaatatatataccacttagcttttgctgcataacaaataaaacatacaaggcgtccgcatatataaaagcctaagcgactgatgactgaatgactggaacgaccagttgaccagttgctatgacgtgtgctgaccaccagggggcagacgcttaatgcaggagctgccccctggtggtcagtgcgctcccacagccaacctcctgtggctggtcggcagcttcagcgtccggcgtccattccttccgtgaccgtcgcctcctctcccaaccctgctgggacctttggcccagggctgggatggggaaccagggtgggtggtggtggatgtgacccctttcccagggggctgagggccggttcccttgtagaggctggaggccaacctcccataccccatccctcctcccagctggcaggccccaattggcctgcagtccctctctttgtccccccccccccccagttgtcggcccagccctaatcactggccaggcctagggaccccacctgtgcacaaattcatgcaccaggtttctagtacaacaataactttatcgcttgtacagttagtcagctaagctatctgcttattatgactggactcacttatttatctgggttcagctggctgatcacctggtctgctgacctaagctggaagcctctggtgcaactatgtattcccacaccttccaacaggataacaggatatcctgggcatgtcttcttggagatgggataggtgcagggaggcaagttccagtgtgtaatgccttttcgtgactccacttgtgccatgttagctagaatcccattgaccaaagcaggtcacatagctgagcccagaggaaggggcctgggcaggtaactcgacctctgagagcagtgcacgtagggatgggtaaagaattggagctaaatcatgaaaggtactacattttttatctgaatttagtttttgatatatatataatagtaggaaaatgttctaagtagatggattaaaatcaagcccttttgtttgtagtggtattttgttctggttttggagttagtaggcatggctgatatcaaagctgttacattaactctcaagcatttaaccctgactaattcttcccaatcaaattggtcatttttaactagggaaactcctttaaaaaggtaaacttcctttgtcttggttctttattcaaaatagggttgcaaactttgctagtaagtgggaccaatcggttaacacactgggtgaagcagctagctagtgggaacttaggccaacagaaaggtacatgcccaactccaaggcattcatggtttatgtccatgttgccatgggaaaattcatcccagtgttgccagatctctgggttttttaagagaagctggaaattaaaatttttatatgaaatttagttgttttaaaattagtgttttaaaataacatagtgagtgttattttcaaatacatgttagtcagctgggtttggtctatagtcctccagttttcatcctctgatttaagcaaagtggtgtcagaaactcgggattctagctccatctttgtccaggctgatctcttgcccagcctggggcagtttggacaatcagtagaagttgtctgtacccccagaggactttattcgtatctgctattttaagaatattgtagtatatatttcattttcaatcggttagaaaatgcaaaagtccaattccagggtgtaactacaagtaatagatccatctcgagtttgctagcatttctctgtaattttcttgagtttttgttgttatttttaatcttgactttttaaaagctttttacttaattttgtaaactgcctcatttattttagaagtgtgtgtgggtgtatgaatctaaaatatctgtgtcattctgcttgtgagacagggaggaggttattaagtaactttcttttaaaactgcattgatctttaacaagtggcacatatgtaattatctgagcaacaactaatttttttctctgtgcctgaataaggagtcattttgtatgctgcttgtccttttccagaaatgcccctgctccctgcttggtgggcactggtccctgtgactctttccctggactcaggagctggagcgaggggagcaaagtgcatctgggagctgagcctgacccgctccacccataggatagagggtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtataagagagagagagggaggaaagggagggcagcaacaggaataagtggaacatctatactaataaaagggtaatatgctaattagaccagacatcttctggatgaccgtctagacgtccttccggacacagccactgtggtggggccgaggcagaggcagttaggggcaatcaggcgaaaaggggggagtggttagggggatcaggcaggcaggggagagcagttggagtgatcgggcaggcaggcagagaggttaggggcgactgggcagttggacatcccccgagggatcctggattggagagggtgcaggccaagctgagggaccccctcccccccccccgccgtgcacgaatatcatgcaccgggcctcttgtacaatataaaatcaaatctagtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtataataattaaaaacattaattctctaagtcagaagtctcaaagtgggatgattttacccagtcacccctccctaccccagggacatttagcagtactggagatgggtttaaatccattccctgatggaagaatgcaaggtgctttgcttatttgtaccttctctagttcttctgaaaaagaagttttctaactggatcccctatccttgtctcatctatttgtcaaacttatgctcctccagatgggtgggaggatcctgaccccactgctccggcgttttgttgtagtacttagctccctgtagagagcacttttacagctctttatctcccatcaaggactgtcacttaatttacttactttcatctttaatttctaatgcagatataattgagtaaatgattgcacaaatatgttcaactttctctctctctctccttttaatcacctgaggatattttttttattgattttagagagacagaaagggggagagagagagggagaaacatctatgtcagagagaaacattgattagttgcctctcagtggatctctgatcgggatcaaaccctcaaccttggtatgtgtcctgactaggaattgaaccggaaacctgcaggggatgatgctctaaccaactgagtcacccagccaggattcatgtatctctttgatttttcccctatgtcttatgttatctattgacttcacagggaattggctgttgattttgatcttatatgaacttaccaaaagggtttgaggagattaagttacagaggagtattgatatcaagaagccaggcatgactgagaattgtttgtacatgataacttttcctggctagttcactttgtctcttgaaggaggaggcaaagctctcaagacctgcagctgttgaaacagtgagaaatgccttctaaaatttgaactttgataatgattaagcttagaatacactccccacatcctgtgatccccttaaatatcaacatcatttccttgaatgtgttctttgaggacagttgccttggacagggatggctgaggcacttggaggccacagaagtaatcacctacttaacatgcactagtcagcttttgtaaaatgtcacccttgatgctgggtgcttatcttttatggaacagcttgcaaattctaattgttttaaactggctctgaggaaaagagacttgagaatcaataaattacagctctcatacatacaaatatgcatcattataatgaggacaatctgatggctttttaggtcaaataacagaagaggtggacttgagtttgtacagactcactgccagtcctcttattacagggccaaatgcacttttactccccctgctcacagtgttcaagtttgcacatgggcatttgtttttgatggactttttcctcagtgggagcgaattgaggagcttcattactcactctggtggcggtatgtcatctttataaggctcgttaaggcagaagcaacaatatttatgcctttaaatactgatgtgttgggaacagaaaaataaagtatttttagagttctaccatctggtttcagacccagagctttctctaagactatggagactttcacccaaggttgcctgtcctagtcctctctagtcttggagaaagctcaatgtctgaaaagcagacagatcaatgctattgcttcctgatcaatgagaagttttgttgtactagatggtggtagtaagagaatattttctgctttttctccagttgaaggaaacgcttttctctcctttaatcagaatttgggtgattcagttttttaaaattgttctttaagagttcttcaatcgttcttcaaatttagcaattcttttccccccaacagatgtcatgttgaaacttcttgtgataatgagaatggatttttaatccaactgaactttttccttttgaatatttcaaaagaatgatgaagattttccttgagcatacatagatgtttttctaccagtcttaataattcttctttgggcttcatattgtaattgtggttaatataattcatgaaccatttaaatttgtttttcaaagcagaatttgattgcatactttatgtgagtacacattaatgtgttttcacatagtccttcaacttttaagttgagccttggccaatgtggctcagttggttgagctttgtctcattcacctagaggtcagggaagggcagttgtgggcgatcaggccggcaggggagggcagttgtaggtgatggggccggcaggggagggcagttgtggccaatcaggccagcaggggtgggcagttgtgggcaattgggctggcaggggagggcagttgcgggtgatcaggccggcaggggagggcagtagtgggtgatcgggccagcaggggagcagttaggcgtcaatcaggtcagcaggggaatggttaggggcgatcaggctggcaggcagaagcggttaggggcaatcaggcaggcaggcaggcaggcgagtggttaggagacagcggtccctgtgggatcaggcctaaactggcagtcggatagcccccaaggggtcccagattggagagggtgcatgctgggctgagggacacccccccccacccctgtccacaaattttgtgcaccaggcctctagttggtatataagaaggccatagatttctgggtgttaattttgtatcctgctacattgccgaattcatttattaagtctagtagttttttggtggagttctttagggttttcaatgtacattcccatgtcatctgcaaataatgacagttttacttcttcttttccaatttggatgccttttatttcttcttcttgtgtgattcctatggctagcatttccagcactatgttgaacaggagtggtgaaagtggacatgtctgtcttgttcctgttcttaggggaaatggttttagtttttgcccattgagtatgatttggctgtaggtttgtcatataaggcttttattatgttgaggtatcatccctctagtcccactttgctgagagttttttttttttttttttataaaaaaagggtgttggattttgtcaaatgccttttctgtatcaattgatatgatcatgtgatttttgtttctcaatttgtttatgtgatgtatcacatttattgatttgcatatattgtaccagccttgcatctctggaataaatctcacttggtcatggtgcatgatctttctaatgtactgctggatccgatttgctagaattttgttgaggattttagcatctatgttcattggggatattggtctgtaattctcttttttttgtggtgtctttatctggttttggaattggtgtaatgctggcttcatggagcttggaagtgtgccttcctcttgaatttttggaatagtctgagaattgaactgtgacctcctggttcataggtcgatgctcaaccactgagtcacatcggccgggcaggagttgctcttttaaagcagaaaaaacaaaggagggagaattacgtgaggaagtgacttatattttcaaattgtatgcaattattatgataaaccaaaggactaatgtaataaatactgtaaaagatcaaaagtaaatagctaaattaatctcttattcccattctctttaccatttgcttgaccatgatctgaattcatttttaaggatttgtgtacatcaatggcagaatcctccagtgattcagaccacttccgctgtcatgaccggttgagtcgatgggctgccaggtcaatacacagggatattagaaatcgtcctacagtcgatgtcaccaagaaggtcaacactatcacaagtactttacaggttagtttaacactggatgctttaaatcaatattaagtaataaatagttaagctgagtaatagtttcagtcactgtattagtttttttgcgcccataacaaattaccacaaaattaatgtcttaaataacagaaatttatttttttaaaattctgtaggtcagaatcggacactaatctcattgggcagtaatgagggtgttggtagggctgcatttctttctggagaccctagggaaatacctctcctcttgcctttttctagcttctaggctgcccatattttctgtcttgtggccccgttcatctccaaagccatcaggggtgggttgagtccttctcatagcatatcagtctgactacctctcctgcttctctcttccacttttaaggacatttgtgatcacattgagtgcacccagtaatccaggatattgctcatctcagcatccttaatttaatcacatctgcacaatcccttttgacatgtaaagttccaggaattagggtttgaacatctttggaagcccttattctgcttaccacagttagtaaggatatttgaaatatatgataggttgtggaaagtgttatatttccaaataaaatattataggaatgaaaacattagggattaccatcctttagcttcttgtggacatgtgccttggagaaattgctcatattttggcttgggtaaatatatcaatcagatagattgccaactagatgatatacatttagcttctattgaatttgaaagtgttaagtcagagggggaagaactcaggacttgaattaaactgagagttactgattattcatccaagccaaactatgattattcatggatgaataatagtgagagtccattttgatttggaggctaactttgacagtgtttggagatttgagtgagggcagttttcattcttccataaaaggggcttgggaaatgaggtgggggaaggtactacagtcacaatgggaaaaagaaaactcagaaatagaattttttcctcctacttattactcccttacctttgttttacctggataagcatggattttggctttttattgctctgggttaatatgcttagattaaaaaaaggctaatttctagtggagaggccataatttgatttaatttttctttttaaacttaatatggcaatgaattttggagtacaaattcttgtactttcccaccacagagacacattagtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtattaaattatcacctaaggagtgaggctatttttttgcattggttttaaagagagtggaaggaaggggggttgggtagaaagagaaacatcattgtgacagagagacatcaattggttcctcccacattcaccctgattgggggtgggaattgttatctgtgacccagctacctgcccttgacctggatttgaactcaagacatttcggtgggcaggccagcgctctaaccattgacccacaccagccggggcagtaagtatatatttttgaaagaaatacctccattaatctctcaggtgaaacctagcatatttaggtagaaaccatttttgttatgattagagctaaatttgttgtttatagctaaaggttaaaggatacttttcccccctcctaattcagttctcagataaagcaataattgaaataacatgaatcccaatgatttgcttctcttgtagataatctactagattagttcttttagcattctttttcatgaagggatgtgaattagttgaggtgtacaagttgataatgagggttttttccagtttatggttttgttgctgacttacggaaccttggataaggatttaaattaggtactttagatttcacacaagtaaattgcagatggattttagtttctctctcagttgtgttaagatttatttgatatttatgaaatatagatatctgacttaactatgatataaattgtacctttgcaaatggttgacttagctgtttgaaatactactatgactgattgggaaaaacatcaactaataaagtccattttctttctcaaggacaccagtcgaaacctgcgacaagtagaccagatgcttggacagtatcgagactatagtaacggacaggcgggtgctatagaacacgtaagatgtttgcatgtttgtattttctcttaccttcttgggaatgcagagctgagtattagggatgtgagctaaaggagataataattgttaacagctttaacttttaaactgctccactattttgtagtataaatttaggttagctcatcagttagattgccaactagatgatatacatttagcttctattgaatttgaaagtgttaagtcagaggggaaagaactcaggatatgaattaaactgagagttactaactttgggctctcatggggaggtaaaaattttaaggacagtaaaaaaattttaaagtaatccatgctcattgtttaaataataaaactccgatctcatagcccttttgagtcccattccaaaaaggtaacagttttgatagtttcttatatattttttcaaacattttctatgcatgcataaacttacatgtgattgctttataaattgtaatcatatgtattgttatgcagactttttcatttagtaatatatgtattattttccatatgacacattcttttgaaattacattatactgcatagtataggtgtattaagtttttttgtttttttttttttttggtattgtaaataatgcaccaatgaccatcatggttctcatgtttttctgtaatcatatgtgattcatgggaacagagttgctgggttgaaggatatgtatacttacaatttgctagatatttcaaaattgccctcccaaggaagtgccagtttgtacttacttcacatgggggtgtaagaaggaaacatatccttattctggagtaaaaggaatatttttttattaatattttgttaatgtaatctatattattttctcttttgaaatttgttaattcatctgcaaactcttggaaatccaattaaatttatttctccatgttttcaatatgagtttacactttaataagactagaattatatttgtactataatgattgccatatatatacatatcctatataataaaagcctaatatgctaagtgcctggtggtctgttcaaccaatcaaagcataatatactaatgatatgctaagtctgctcaaccgctcactgtgacgtgcactgaccaccagggggcagatagtcgactggtcgaccagtcactatgacatgcactgaccaccagggggcagacgctgaatgcaggagctgctccctggtggtcagtgcgctcccacaggaggagctccgctcagccagaagccaggcttatgactggcaagtgcatcgggggtggtgggagcctggcagttggcagttggacatctcccaagggctctcagactgctagagggcacagggcaggctgagggacacccccccgttcccccccgagtgcacgaatttcatgcaccaggcctctagtgtgtgatattcaatatatatatactgaattgtcacataacttcctagttataattgaagatttggtattttctttttataaaacagaagaatagttagcattatttgcttgtaattggaaactaattgcttgtgcacagggaaatgtatcttctgtcttaaggtttcccactatcatttagaaaataaacttggatttttgcagctaaaagattttcatgatttatgatgtatttgacagatggagctatagattttgtatatttgggtagctcattaaatgctgtagattaaactactttattaataatgggagttaatactctgccttttaattttagtgaacaattaaaaacactaatagcagagaacaatcctactttgtcatagtgaatattagatgcttgcttatagtagcaaattttctttgaatttgcaattgggaaatgaaattaagaattagagaaacttctgtgcaataagcatttaaaacacaaagtatgttctactatggcagtaagcctttgtggtaggccttgagttgaataatgctgaactctgggggctttatatatgggaagtaaatacattcaacaggtacttcatttccatctgagttgtaagtacatagtatgcaaagccactttataaaaaataagacacaggagagcaagtgtctttctttggcttcagactttttcttagtttggtctagagtctcctcctgaatttctccaaattcttgtcctgttaaagcagttggttctcagtgattttccctgccaacaaatcagaggagtaatgtcactcattgcctgacattatactcacttgatagagcagtgatcttcatcgccccacctataagttacctggggagcatctgccataatctattggaagatgtgcaactcaggaagagtatttctcacaataagtgattctgatgggctcctcaatgaaaattgttttcttagcattaaaatatttactgttatggaatattttacacctaaaaatgagtatatgtaatgtgctaaggatcatgtgatttatctgaatctttctatttatattctttttttaaaaacatatttttattgattttagagaggaagggagtggaagggaaagatagaaatatcaatgatgagagagaatcattgactggttgcctcctgcacacttcctactggggatcgagcctgcaatccaggcacgtgccctgactgggaaccgaaacgtgacctcctggttcattggtcgatgctcaaccactgagcaacaccagttgggctctatttatattcttgaagtccttgaagtacatggtagtattacccattttttcttttcttttaatcaatgtgtgcattgtttagggaggaaccctggaaattgatgtttttatttttgaggctatgctacaacttatacaagtcctcttgctgtgaaacttgttataaaaaaaaaaaaaatcccgttacagtttaatctgttttaaccagtactcaactgttaagaattttctcgctgctaaatatctattgcctcattgaaatgagcaattcttttttaaaatataagggatttcctagaaaatgagatcaactcttttgctttagggaagacaagaatttttattattaacttagattttataaatgattctattgctaagcctcaggacatggatgttttatatatgtatttaacatagttatttgattggggtttaattttgtgaacataatttgccagttggtgagagggaatataaataatcgatggggtgttttcctcctgaggtgacctcacccaacaagagaactatgagaatcattgctgttggttttgtgaggtcttagggtcaagagtcatgttgggtcatatttttgactacttttccttgaaatgtacataatttcattgttgaccccttgaatctttctaacatagatttaattctttttttttttttaaaaaatatatattttattgatttttttatagagaggaagggagagggatagagagttagaaacattgatgagagagaaacatcaatcagctgcctcctgcacactccctactggggatgtgcctgcaaccaaggtacatgcccacaaccaaggtacatgcccttgaccggaatcgaacctgggacctttcagtctgtaggccgacgctgtccactgagccaagccggttagggccatagatttaattctttgcttataatatgtgaagctgagcagtgttagagggtagcagggcatctgtagtgttttgtagagaagcagcacctgccttttcacagtttccctgggaatttctctgtaagcttagcatatcccaggtgcttttaaactatcaactcatttaattctctcaccaatcctgttaggtagatactattattttcctcatttacaagataaaaaactgggacacagaatggttgctaacttgcccatagtctcttagctagtcagtaacagtctacgatgttattctaaaacctagtgcctggattcatttacctcatttacctcatggttgctcacatttttgtgatgttgtgaatctatagttctgaaaattgggggtttgaactcaaagcttaaagggttatataattcacattgataaactagattttggagattaccttaactgattgaaaaacagaagcccattgtcacacaagtcattagtgttgtagctggtgctgtgtcagaaattagagttagatggaacaaatacagaaacagcattgggaagaaatctggtgctgagctaggactggcatggacagattgtcagggttggcagttttggctggcacttaggcaacatggtggattctgcagcttccctacaccatcattgacttccccttttctcccgagtaacaaccattgatttttcactattccagccagtgagccatatttagttggatacaatgacatagggtctttgggctccttgactttcataaccttaaagaagggtcattatttctcattctgtaaatttcttgggaggtctgagaacatttgggtcaaagaggacctgtggggaaggttttagttaaaagcttattttaaacatttgttatatagtgcaattgattagcttggtgtttagaatatttctagaactctgtgttcctattcataagacaatagggactgtttttgagttaatgctataaaaatggagtgccctcaagaaagaatattcaacataaaaagttaatgctttctgcagtggaagtgtagctgatgatgctagagtagttcatattgggtataaaagcaaaatacagttcccttgattcagatcagctgatttattcatttttaaggtaaagatgggtgaaaatatagtaatagcaaattaagccttagaggcttgaattaataggtaatttttctcctaattaagtgttctgtttaaaatttcttttgggagtaataggatgttaccatctggtacatgaactttggaagaatttaaatgttaaataataaactccaggtttttgttctatttggtgggtttcagttgtcaaatatactgagtgataccctaaggttaatggtttagagccaggcacaagccacttgtaatctctcaaagctgaaggcatgtgcccctttctctctggccaggttttactagttgtcctgaagttgagaaccagtcacagtcttttaataactttaatgctagatagtaagagtctatgagtaagaacttgtaatcgaacaattaaagtcttctggcttctttttgtaaatacatttcaaatcttagaaattaggttcttagaattctgtttaaaattttattaccaagtctcacatgaaacaactggataagagaatattttttacagacttttgctttgagtctagaggaaaaatctagattttggaatcagacttccaattgaatattgatttttctgctagaaaccatttagttactttgggcaagctactaaaccacctctatgttcatgagtaaatcgtctgcccacccaaaatacagtaatgtactattaagggtaccttccattttacatagtagtaatgcttcagtttttttcggtctaatatcaaatgagcagaggtaaaagaatagatttaacatgtaaatgaagcatttatagcagattgtgaaaatgtcattgcatgttttgtgatgtgggacagagagaaatggaatcatcacttggtctgtggtgatgttgtcatgctagactctaagaatctgcacttaaaaaaatttttttttttaatttcagagagagaaaaagagagggagatagagatagaaacatcgatgagagagaaatatcgatcgactgcctcctgcttggcccctgctggggatcgagcccgcagactgggcatgtgccctaatcaggaattgaccttctggttcataggttgacattcacccactgagccacactggctgggcaagaatctgtacctttggaatcagtaccacagtctcctgtggagcattgttttcatttcagcatagtagcatactgaattaatgttaaattgagctttattagtgttgttttgttttatttaacttcaagttttagtttgcttttgcacctagttaggagttataatcctaaaacatttgaagaaacttttattattttacattttgtctcttcatctctttgattttgtgtttcaggagtgatcaggttagtttcatgcatcatgcttagattttgtgcagtgtctgctgtgttattgctactttaagtgtatattttaattctgttacatttttagttttgtattgtaggtggaggacatgtttttgccttgatattaaaaaaatctttatcctgcttatttacttttgcatttcattttttggtctttatcagcctgttcatctgttcaatggcttttgatataccttttaaagtttaaaatgaagtttcaaatagtttatcttctgaaatgattaatttctcttgttctgcaatatttttccctatgttccatggtggctttccccccatttttgttaacgggttgtgtgagtgtgtgtgtaagtttattttattcatcagcttttagttacagctcaggatttgttagcaaagaagcccacgtaaaaaggatgctggcatactctgcctctagcatcagtagataagattattttggatctactctcttactctgaacaactagaaagg";
    auto bSequence = reverseComplement("cctttctagttgttcagagtaagagagtagatccaaaataatcttatctactgatgctagaggcagagtatgccagcatcctttttacgtgggcttctttgctaacaaatcctgagctgtaactaaaagctgatgaataaaataaacttacacacacactcacacaacccgttaacaaaaatggggggaaagccaccatggaacatagggaaaaatattgcagaacaagagaaattaatcatttcagaagataaactatttgaaacttcattttaaactttaaaaggtatatcaaaagccattgaacagatgaacaggctgataaagaccaaaaaatgaaatgcaaaagtaaataagcaggataaagatttttttaatatcaaggcaaaaacatgtcctccacctacaatacaaaactaaaaatgtaacagaattaaaatatacacttaaagtagcaataacacagcagacactgcacaaaatctaagcatgatgcatgaaactaacctgatcactcctgaaacacaaaatcaaagagatgaagagacaaaatgtaaaataataaaagtttcttcaaatgttttaggattataactcctaactaggtgcaaaagcaaactaaaacttgaagttaaataaaacaaaacaacactaataaagctcaatttaacattaattcagtatgctactatgctgaaatgaaaacaatgctccacaggagactgtggtactgattccaaaggtacagattcttgcccagccagtgtggctcagtgggtgaatgtcaacctatgaaccagaaggtcaattcctgattagggcacatgcccagtctgcgggctcgatccccagcaggggccaagcaggaggcagtcgatcgatatttctctctcatcgatgtttctatctctatctccctctctttttctctctctgaaattaaaaaaaaaaatttttttaagtgcagattcttagagtctagcatgacaacatcaccacagaccaagtgatgattccatttctctctgtcccacatcacaaaacatgcaatgacattttcacaatctgctataaatgcttcatttacatgttaaatctattcttttacctctgctcatttgatattagaccgaaaaaaactgaagcattactactatgtaaaatggaaggtacccttaatagtacattactgtattttgggtgggcagacgatttactcatgaacatagaggtggtttagtagcttgcccaaagtaactaaatggtttctagcagaaaaatcaatattcaattggaagtctgattccaaaatctagatttttcctctagactcaaagcaaaagtctgtaaaaaatattctcttatccagttgtttcatgtgagacttggtaataaaattttaaacagaattctaagaacctaatttctaagatttgaaatgtatttacaaaaagaagccagaagactttaattgttcgattacaagttcttactcatagactcttactatctagcattaaagttattaaaagactgtgactggttctcaacttcaggacaactagtaaaacctggccagagagaaaggggcacatgccttcagctttgagagattacaagtggcttgtgcctggctctaaaccattaaccttagggtatcactcagtatatttgacaactgaaacccaccaaatagaacaaaaacctggagtttattatttaacatttaaattcttccaaagttcatgtaccagatggtaacatcctattactcccaaaagaaattttaaacagaacacttaattaggagaaaaattacctattaattcaagcctctaaggcttaatttgctattactatattttcacccatctttaccttaaaaatgaataaatcagctgatctgaatcaagggaactgtattttgcttttatacccaatatgaactactctagcatcatcagctacacttccactgcagaaagcattaactttttatgttgaatattctttcttgagggcactccatttttatagcattaactcaaaaacagtccctattgtcttatgaataggaacacagagttctagaaatattctaaacaccaagctaatcaattgcactatataacaaatgtttaaaataagcttttaactaaaaccttccccacaggtcctctttgacccaaatgttctcagacctcccaagaaatttacagaatgagaaataatgacccttctttaaggttatgaaagtcaaggagcccaaagaccctatgtcattgtatccaactaaatatggctcactggctggaatagtgaaaaatcaatggttgttactcgggagaaaaggggaagtcaatgatggtgtagggaagctgcagaatccaccatgttgcctaagtgccagccaaaactgccaaccctgacaatctgtccatgccagtcctagctcagcaccagatttcttcccaatgctgtttctgtatttgttccatctaactctaatttctgacacagcaccagctacaacactaatgacttgtgtgacaatgggcttctgtttttcaatcagttaaggtaatctccaaaatctagtttatcaatgtgaattatataaccctttaagctttgagttcaaacccccaattttcagaactatagattcacaacatcacaaaaatgtgagcaaccatgaggtaaatgaggtaaatgaatccaggcactaggttttagaataacatcgtagactgttactgactagctaagagactatgggcaagttagcaaccattctgtgtcccagttttttatcttgtaaatgaggaaaataatagtatctacctaacaggattggtgagagaattaaatgagttgatagtttaaaagcacctgggatatgctaagcttacagagaaattcccagggaaactgtgaaaaggcaggtgctgcttctctacaaaacactacagatgccctgctaccctctaacactgctcagcttcacatattataagcaaagaattaaatctatggccctaaccggcttggctcagtggacagcgtcggcctacagactgaaaggtcccaggttcgattccggtcaagggcatgtaccttggttgtgggcatgtaccttggttgcaggcacatccccagtagggagtgtgcaggaggcagctgattgatgtttctctctcatcaatgtttctaactctctatccctctcccttcctctctataaaaaaatcaataaaatatatattttttaaaaaaaaaaaaagaattaaatctatgttagaaagattcaaggggtcaacaatgaaattatgtacatttcaaggaaaagtagtcaaaaatatgacccaacatgactcttgaccctaagacctcacaaaaccaacagcaatgattctcatagttctcttgttgggtgaggtcacctcaggaggaaaacaccccatcgattatttatattccctctcaccaactggcaaattatgttcacaaaattaaaccccaatcaaataactatgttaaatacatatataaaacatccatgtcctgaggcttagcaatagaatcatttataaaatctaagttaataataaaaattcttgtcttccctaaagcaaaagagttgatctcattttctaggaaatcccttatattttaaaaaagaattgctcatttcaatgaggcaatagatatttagcagcgagaaaattcttaacagttgagtactggttaaaacagattaaactgtaacgggattttttttttttttataacaagtttcacagcaagaggacttgtataagttgtagcatagcctcaaaaataaaaacatcaatttccagggttcctccctaaacaatgcacacattgattaaaagaaaagaaaaaatgggtaatactaccatgtacttcaaggacttcaagaatataaatagagcccaactggtgttgctcagtggttgagcatcgaccaatgaaccaggaggtcacgtttcggttcccagtcagggcacgtgcctggattgcaggctcgatccccagtaggaagtgtgcaggaggcaaccagtcaatgattctctctcatcattgatatttctatctttcccttccactcccttcctctctaaaatcaataaaaatatgtttttaaaaaaagaatataaatagaaagattcagataaatcacatgatccttagcacattacatatactcatttttaggtgtaaaatattccataacagtaaatattttaatgctaagaaaacaattttcattgaggagcccatcagaatcacttattgtgagaaatactcttcctgagttgcacatcttccaatagattatggcagatgctccccaggtaacttataggtggggcgatgaagatcactgctctatcaagtgagtataatgtcaggcaatgagtgacattactcctctgatttgttggcagggaaaatcactgagaaccaactgctttaacaggacaagaatttggagaaattcaggaggagactctagaccaaactaagaaaaagtctgaagccaaagaaagacacttgctctcctgtgtcttattttttataaagtggctttgcatactatgtacttacaactcagatggaaatgaagtacctgttgaatgtatttacttcccatatataaagcccccagagttcagcattattcaactcaaggcctaccacaaaggcttactgccatagtagaacatactttgtgttttaaatgcttattgcacagaagtttctctaattcttaatttcatttcccaattgcaaattcaaagaaaatttgctactataagcaagcatctaatattcactatgacaaagtaggattgttctctgctattagtgtttttaattgttcactaaaattaaaaggcagagtattaactcccattattaataaagtagtttaatctacagcatttaatgagctacccaaatatacaaaatctatagctccatctgtcaaatacatcataaatcatgaaaatcttttagctgcaaaaatccaagtttattttctaaatgatagtgggaaaccttaagacagaagatacatttccctgtgcacaagcaattagtttccaattacaagcaaataatgctaactattcttctgttttataaaaagaaaataccaaatcttcaattataactaggaagttatgtgacaattcagtatatatatattgaatatcacacactagaggcctggtgcatgaaattcgtgcactcgggggggaacgggggggtgtccctcagcctgccctgtgccctctagcagtctgagagcccttgggagatgtccaactgccaactgccaggctcccaccacccccgatgcacttgccagtcataagcctggcttctggctgagcggagctcctcctgtgggagcgcactgaccaccagggagcagctcctgcattcagcgtctgccccctggtggtcagtgcatgtcatagtgactggtcgaccagtcgactatctgccccctggtggtcagtgcacgtcacagtgagcggttgagcagacttagcatatcattagtatattatgctttgattggttgaacagaccaccaggcacttagcatattaggcttttattatataggatatgtatatatatggcaatcattatagtacaaatataattctagtcttattaaagtgtaaactcatattgaaaacatggagaaataaatttaattggatttccaagagtttgcagatgaattaacaaatttcaaaagagaaaataatatagattacattaacaaaatattaataaaaaaatattccttttactccagaataaggatatgtttccttcttacacccccatgtgaagtaagtacaaactggcacttccttgggagggcaattttgaaatatctagcaaattgtaagtatacatatccttcaacccagcaactctgttcccatgaatcacatatgattacagaaaaacatgagaaccatgatggtcattggtgcattatttacaataccaaaaaaaaaaaaaaacaaaaaaacttaatacacctatactatgcagtataatgtaatttcaaaagaatgtgtcatatggaaaataatacatatattactaaatgaaaaagtctgcataacaatacatatgattacaatttataaagcaatcacatgtaagtttatgcatgcatagaaaatgtttgaaaaaatatataagaaactatcaaaactgttacctttttggaatgggactcaaaagggctatgagatcggagttttattatttaaacaatgagcatggattactttaaaatttttttactgtccttaaaatttttacctccccatgagagcccaaagttagtaactctcagtttaattcatatcctgagttctttcccctctgacttaacactttcaaattcaatagaagctaaatgtatatcatctagttggcaatctaactgatgagctaacctaaatttatactacaaaatagtggagcagtttaaaagttaaagctgttaacaattattatctcctttagctcacatccctaatactcagctctgcattcccaagaaggtaagagaaaatacaaacatgcaaacatcttacgtgttctatagcacccgcctgtccgttactatagtctcgatactgtccaagcatctggtctacttgtcgcaggtttcgactggtgtccttgagaaagaaaatggactttattagttgatgtttttcccaatcagtcatagtagtatttcaaacagctaagtcaaccatttgcaaaggtacaatttatatcatagttaagtcagatatctatatttcataaatatcaaataaatcttaacacaactgagagagaaactaaaatccatctgcaatttacttgtgtgaaatctaaagtacctaatttaaatccttatccaaggttccgtaagtcagcaacaaaaccataaactggaaaaaaccctcattatcaacttgtacacctcaactaattcacatcccttcatgaaaaagaatgctaaaagaactaatctagtagattatctacaagagaagcaaatcattgggattcatgttatttcaattattgctttatctgagaactgaattaggaggggggaaaagtatcctttaacctttagctataaacaacaaatttagctctaatcataacaaaaatggtttctacctaaatatgctaggtttcacctgagagattaatggaggtatttctttcaaaaatatatacttactgccccggctggtgtgggtcaatggttagagcgctggcctgcccaccgaaatgtcttgagttcaaatccaggtcaagggcaggtagctgggtcacagataacaattcccacccccaatcagggtgaatgtgggaggaaccaattgatgtctctctgtcacaatgatgtttctctttctacccaaccccccttccttccactctctttaaaaccaatgcaaaaaaatagcctcactccttaggtgataatttaatacacacacacacacacacacacacacacacactaatgtgtctctgtggtgggaaagtacaagaatttgtactccaaaattcattgccatattaagtttaaaaagaaaaattaaatcaaattatggcctctccactagaaattagcctttttttaatctaagcatattaacccagagcaataaaaagccaaaatccatgcttatccaggtaaaacaaaggtaagggagtaataagtaggaggaaaaaattctatttctgagttttctttttcccattgtgactgtagtaccttcccccacctcatttcccaagccccttttatggaagaatgaaaactgccctcactcaaatctccaaacactgtcaaagttagcctccaaatcaaaatggactctcactattattcatccatgaataatcatagtttggcttggatgaataatcagtaactctcagtttaattcaagtcctgagttcttccccctctgacttaacactttcaaattcaatagaagctaaatgtatatcatctagttggcaatctatctgattgatatatttacccaagccaaaatatgagcaatttctccaaggcacatgtccacaagaagctaaaggatggtaatccctaatgttttcattcctataatattttatttggaaatataacactttccacaacctatcatatatttcaaatatccttactaactgtggtaagcagaataagggcttccaaagatgttcaaaccctaattcctggaactttacatgtcaaaagggattgtgcagatgtgattaaattaaggatgctgagatgagcaatatcctggattactgggtgcactcaatgtgatcacaaatgtccttaaaagtggaagagagaagcaggagaggtagtcagactgatatgctatgagaaggactcaacccacccctgatggctttggagatgaacggggccacaagacagaaaatatgggcagcctagaagctagaaaaaggcaagaggagaggtatttccctagggtctccagaaagaaatgcagccctaccaacaccctcattactgcccaatgagattagtgtccgattctgacctacagaattttaaaaaaataaatttctgttatttaagacattaattttgtggtaatttgttatgggcgcaaaaaaactaatacagtgactgaaactattactcagcttaactatttattacttaatattgatttaaagcatccagtgttaaactaacctgtaaagtacttgtgatagtgttgaccttcttggtgacatcgactgtaggacgatttctaatatccctgtgtattgacctggcagcccatcgactcaaccggtcatgacagcggaagtggtctgaatcactggaggattctgccattgatgtacacaaatccttaaaaatgaattcagatcatggtcaagcaaatggtaaagagaatgggaataagagattaatttagctatttacttttgatcttttacagtatttattacattagtcctttggtttatcataataattgcatacaatttgaaaatataagtcacttcctcacgtaattctccctcctttgttttttctgctttaaaagagcaactcctgcccggccgatgtgactcagtggttgagcatcgacctatgaaccaggaggtcacagttcaattctcagactattccaaaaattcaagaggaaggcacacttccaagctccatgaagccagcattacaccaattccaaaaccagataaagacaccacaaaaaaaagagaattacagaccaatatccccaatgaacatagatgctaaaatcctcaacaaaattctagcaaatcggatccagcagtacattagaaagatcatgcaccatgaccaagtgagatttattccagagatgcaaggctggtacaatatatgcaaatcaataaatgtgatacatcacataaacaaattgagaaacaaaaatcacatgatcatatcaattgatacagaaaaggcatttgacaaaatccaacaccctttttttataaaaaaaaaaaaaaaaaactctcagcaaagtgggactagagggatgatacctcaacataataaaagccttatatgacaaacctacagccaaatcatactcaatgggcaaaaactaaaaccatttcccctaagaacaggaacaagacagacatgtccactttcaccactcctgttcaacatagtgctggaaatgctagccataggaatcacacaagaagaagaaataaaaggcatccaaattggaaaagaagaagtaaaactgtcattatttgcagatgacatgggaatgtacattgaaaaccctaaagaactccaccaaaaaactactagacttaataaatgaattcggcaatgtagcaggatacaaaattaacacccagaaatctatggccttcttatataccaactagaggcctggtgcacaaaatttgtggacaggggtgggggggggtgtccctcagcccagcatgcaccctctccaatctgggaccccttgggggctatccgactgccagtttaggcctgatcccacagggaccgctgtctcctaaccactcgcctgcctgcctgcctgcctgattgcccctaaccgcttctgcctgccagcctgatcgcccctaaccattcccctgctgacctgattgacgcctaactgctcccctgctggcccgatcacccactactgccctcccctgccggcctgatcacccgcaactgccctcccctgccagcccaattgcccacaactgcccacccctgctggcctgattggccacaactgccctcccctgccggccccatcacctacaactgccctcccctgccggcctgatcgcccacaactgcccttccctgacctctaggtgaatgagacaaagctcaaccaactgagccacattggccaaggctcaacttaaaagttgaaggactatgtgaaaacacattaatgtgtactcacataaagtatgcaatcaaattctgctttgaaaaacaaatttaaatggttcatgaattatattaaccacaattacaatatgaagcccaaagaagaattattaagactggtagaaaaacatctatgtatgctcaaggaaaatcttcatcattcttttgaaatattcaaaaggaaaaagttcagttggattaaaaatccattctcattatcacaagaagtttcaacatgacatctgttggggggaaaagaattgctaaatttgaagaacgattgaagaactcttaaagaacaattttaaaaaactgaatcacccaaattctgattaaaggagagaaaagcgtttccttcaactggagaaaaagcagaaaatattctcttactaccaccatctagtacaacaaaacttctcattgatcaggaagcaatagcattgatctgtctgcttttcagacattgagctttctccaagactagagaggactaggacaggcaaccttgggtgaaagtctccatagtcttagagaaagctctgggtctgaaaccagatggtagaactctaaaaatactttatttttctgttcccaacacatcagtatttaaaggcataaatattgttgcttctgccttaacgagccttataaagatgacataccgccaccagagtgagtaatgaagctcctcaattcgctcccactgaggaaaaagtccatcaaaaacaaatgcccatgtgcaaacttgaacactgtgagcagggggagtaaaagtgcatttggccctgtaataagaggactggcagtgagtctgtacaaactcaagtccacctcttctgttatttgacctaaaaagccatcagattgtcctcattataatgatgcatatttgtatgtatgagagctgtaatttattgattctcaagtctcttttcctcagagccagtttaaaacaattagaatttgcaagctgttccataaaagataagcacccagcatcaagggtgacattttacaaaagctgactagtgcatgttaagtaggtgattacttctgtggcctccaagtgcctcagccatccctgtccaaggcaactgtcctcaaagaacacattcaaggaaatgatgttgatatttaaggggatcacaggatgtggggagtgtattctaagcttaatcattatcaaagttcaaattttagaaggcatttctcactgtttcaacagctgcaggtcttgagagctttgcctcctccttcaagagacaaagtgaactagccaggaaaagttatcatgtacaaacaattctcagtcatgcctggcttcttgatatcaatactcctctgtaacttaatctcctcaaacccttttggtaagttcatataagatcaaaatcaacagccaattccctgtgaagtcaatagataacataagacataggggaaaaatcaaagagatacatgaatcctggctgggtgactcagttggttagagcatcatcccctgcaggtttccggttcaattcctagtcaggacacataccaaggttgagggtttgatcccgatcagagatccactgagaggcaactaatcaatgtttctctctgacatagatgtttctccctctctctctccccctttctgtctctctaaaatcaataaaaaaaatatcctcaggtgattaaaaggagagagagagagaaagttgaacatatttgtgcaatcatttactcaattatatctgcattagaaattaaagatgaaagtaagtaaattaagtgacagtccttgatgggagataaagagctgtaaaagtgctctctacagggagctaagtactacaacaaaacgccggagcagtggggtcaggatcctcccacccatctggaggagcataagtttgacaaatagatgagacaaggataggggatccagttagaaaacttctttttcagaagaactagagaaggtacaaataagcaaagcaccttgcattcttccatcagggaatggatttaaacccatctccagtactgctaaatgtccctggggtagggaggggtgactgggtaaaatcatcccactttgagacttctgacttagagaattaatgtttttaattattatacacacacacacacacacacacacacacacacacacacacacactagatttgattttatattgtacaagaggcccggtgcatgatattcgtgcacggcgggggggggggagggggtccctcagcttggcctgcaccctctccaatccaggatccctcgggggatgtccaactgcccagtcgcccctaacctctctgcctgcctgcccgatcactccaactgctctcccctgcctgcctgatccccctaaccactccccccttttcgcctgattgcccctaactgcctctgcctcggccccaccacagtggctgtgtccggaaggacgtctagacggtcatccagaagatgtctggtctaattagcatattacccttttattagtatagatgttccacttattcctgttgctgccctccctttcctccctctctctctcttatacacacacacacacacacacacacacacacaccctctatcctatgggtggagcgggtcaggctcagctcccagatgcactttgctcccctcgctccagctcctgagtccagggaaagagtcacagggaccagtgcccaccaagcagggagcaggggcatttctggaaaaggacaagcagcatacaaaatgactccttattcaggcacagagaaaaaaattagttgttgctcagataattacatatgtgccacttgttaaagatcaatgcagttttaaaagaaagttacttaataacctcctccctgtctcacaagcagaatgacacagatattttagattcatacacccacacacacttctaaaataaatgaggcagtttacaaaattaagtaaaaagcttttaaaaagtcaagattaaaaataacaacaaaaactcaagaaaattacagagaaatgctagcaaactcgagatggatctattacttgtagttacaccctggaattggacttttgcattttctaaccgattgaaaatgaaatatatactacaatattcttaaaatagcagatacgaataaagtcctctgggggtacagacaacttctactgattgtccaaactgccccaggctgggcaagagatcagcctggacaaagatggagctagaatcccgagtttctgacaccactttgcttaaatcagaggatgaaaactggaggactatagaccaaacccagctgactaacatgtatttgaaaataacactcactatgttattttaaaacactaattttaaaacaactaaatttcatataaaaattttaatttccagcttctcttaaaaaacccagagatctggcaacactgggatgaattttcccatggcaacatggacataaaccatgaatgccttggagttgggcatgtacctttctgttggcctaagttcccactagctagctgcttcacccagtgtgttaaccgattggtcccacttactagcaaagtttgcaaccctattttgaataaagaaccaagacaaaggaagtttacctttttaaaggagtttccctagttaaaaatgaccaatttgattgggaagaattagtcagggttaaatgcttgagagttaatgtaacagctttgatatcagccatgcctactaactccaaaaccagaacaaaataccactacaaacaaaagggcttgattttaatccatctacttagaacattttcctactattatatatatatcaaaaactaaattcagataaaaaatgtagtacctttcatgatttagctccaattctttacccatccctacgtgcactgctctcagaggtcgagttacctgcccaggccccttcctctgggctcagctatgtgacctgctttggtcaatgggattctagctaacatggcacaagtggagtcacgaaaaggcattacacactggaacttgcctccctgcacctatcccatctccaagaagacatgcccaggatatcctgttatcctgttggaaggtgtgggaatacatagttgcaccagaggcttccagcttaggtcagcagaccaggtgatcagccagctgaacccagataaataagtgagtccagtcataataagcagatagcttagctgactaactgtacaagcgataaagttattgttgtactagaaacctggtgcatgaatttgtgcacaggtggggtccctaggcctggccagtgattagggctgggccgacaactggggggggggggggacaaagagagggactgcaggccaattggggcctgccagctgggaggagggatggggtatgggaggttggcctccagcctctacaagggaaccggccctcagccccctgggaaaggggtcacatccaccaccacccaccctggttccccatcccagccctgggccaaaggtcccagcagggttgggagaggaggcgacggtcacggaaggaatggacgccggacgctgaagctgccgaccagccacaggaggttggctgtgggagcgcactgaccaccagggggcagctcctgcattaagcgtctgccccctggtggtcagcacacgtcatagcaactggtcaactggtcgttccagtcattcagtcatcagtcgcttaggcttttatatatgcggacgccttgtatgttttatttgttatgcagcaaaagctaagtggtatatatattgtgaacactgcgaaagaagatacctaccaagaaggtctcctggtggctaactgggctcaaattttaaaaacagagactagcagccatttccacacaaactgcacttcacacagaagcctgcactgaactgcctgcctgcactgtttttctgtctgctgagccagcccttataaaggaagtgctggaatcctatttcaaccaatgggactgaggcttcccccatccaaacagagcttacagctttatccccatcaacatcttcactgaccaatcagagcagctcttttctgaccaatcaggaaagtgccttttggacaaatccaactacaaggatctggagtccttgtttgcatgaggacagaccaatcaaggtcacaggcaagggcttctgtctatataagccagctcccctctggctctgagagtgcacctttcctttccaccaaagactgatccctggctgagataaaacaccaaagatatctctgtctcactggagcagagctgaacagaccaacatgggccagagaagacaggcaccaaacacagaatggagaagagcagaccagcacatagctgaacagaccagcaaggaacagagatctctagcagagaggggcctggccccaggtccacctcatcggcttgctgttttcacggctatgctttatcacaattgagctgttttgcactgaataaattttccttcttcaccacaccccaaattggggattcctgttggcgttgaattggcaccaatgatcatctgttgacaatatgtagtttactaataaatcattttattggttcagaaaaaatttaatttcctgaagtttcctctaaagactactctccttgagaaaacaagttgttttcagactgtcttctggtggtgaggagggagtggggagggcaaattcaaattaggatgcacaaaataagcataagcaattgatatgacaagtgacagactagcaattattgagagagaaaaactaggcaatcataaatagcaattggggctatttaatctactccccataacattaagattttcataagaaacaagctgtttcttaaatatattataaagcataatacttggtttttttcctcctattttttaaagtagaaacagaaactcctcttgactttcctagtgaaaaaacaccactaaaaattttaagaaagggaaatatttctgactgctatttatacctcaaatttctatgttgctttatagcttacaaagcattttctttagcactttctcttttgacttccattcaatcctacaaagtgggtataattgcaactctttattttactgctaagaaaaacaaagctttgaataaaatacttgcccaccattgctttgataaaaagttacagtaacagaactagaacccatatcttgtgtctctaagccagtgcttttccactacatttctaggagcagaactgtgtaaattctaaatacattgcccaatcactgaaatattttctttgtccatccctatcatatttcttatcaccctcactcccccaagaaatgcatttaaggccagagtatgtcatccatatgtgcttatctacttaaaactagagatttctataatgaatacaagataatatgtaaatatcctataaatccaactttctacaactattttaaaatgtcttttaaataaccctttctcttgttaagaatatttacctcttcaggcagaaaagcattctttttttttctttttttttttcagaaacatcataaaatcttaaattttattttttaaattgtattttctttttggtgtacacgtattttctatttctctttttcctccactcttctcaattcccatatttcttccttgccccacacttactttatatatttagatctttccttttaactcttctataataaaggactgaaatacatttcaagagttttgtatttcttttttatcattattaaatttattggggtggcattggttaataagattatataggtttcaagggtacatttctatgatacatgatctatatattgcactgtgtgcccaccatccaaagtcaaaccattttctgtgaccatatatttggtcccatttatcctttactaatcctccaccccctctcctctggaaactatgagtttcagttttatatcccatatatgagtgaaatcatatggttcttaggttttctggctgactttattttgcttaacataatatcctcaattttatttcttgttcttaattttcaaagtggataaacagaataattggcatattcttttcaaaatgtatagacttctaaattttaagaaaaaaagctaattgaagttcttaactttgatttccccatgaactgtttcctggaatgtttttaaatgtacataaaaaatttttaaagatactggattagatagaaatccaattataatgaaattgaattattgggagtaaatggatgaggactgaagttaagaacctcaatttcaaaaaccttcaagaaactcaagttgaataatgctgaataactaatgtttcaatatcagattgttgtttgacaataaccctacaaggtcagtacactgtcatgatatatacctaagaggcttaaaggcactaagaaaggtactcaaggctacacagctactaagtagcagagctgagtcttgaatctttctccaaggtccagatttctgaccatctaaccgtgctcactgccacgtgcactactcataaatttaagactttttctcaactcatttctagccaccaaaatggtttctttgacttttgagtcttttcctatacattccaactccactggccacctgactgaatctaggttcttgaatcaggaattcaaggccctcggcaagccacctgtcttaataaaagaaccacagaaaacttaaggcctacattcagctattacaacaactattcttcccccactgccacccctcattaatgtcattaggtctgggtggaatgccttctgctatctagtattgccacctatttaaacactcacttgtccaggcatagctcaaataacacaatctcaatgaagctgtccacctctagccagaagtgatctttccttcctctgcatttgtatagctctataattatagcaaaacatttatgaaatgtcaaaatatctgcatgaaatacaattatttgtgccacagaatttttttaaatttttctttatgctcacctaaggatatgttttattggtctgagagagaggaagggagagagagagaaacaccaatgtgagagagaaacatcaatcagttgcccccacaggtgtccagaccagggatcaagccacaacctaggtatgtgccctgatggggaatcaaacccaaaacttttttggtgtaagggatgatgctgcaatcaactgagccatctggccagggctgtcacagaatatgtcatgttcattttgtcatcaaacctatgctggattcaccacttggaacattcttttcacttcgttttccctttgctaacttctacttttctttcaagtctcagctcaaaggaaacttcccctagaaagtcttctttcactccccccacccctagcccaccaagtatccatagtggataccaactcaacactgcccagagattgctaactggctatttacttgtctctctccccactagactggaaatcttaccagagactatttatgtaatccatggttgtatttcaacacctagcatatttcctggcatgtagttggaacacttttaaatgaatttacaattacattttcttattatcctgtctgattaatccagagatgaatttcaaagtacaggaccttgcaattagtaaatactaattaaataattgttgggtaaaatgattcttgtagaaatggcctagattaattatcttacatctccaaacaatatggattaaactgcaaaatagtaaaaatatgaacacttgataatatgtcctttatagagaaaacaaagataaacactttaaatccatttgaataaatgcgttaaaacaaatgaatcacacttaagcttgtttgaattattgcaaatcatttaaatacattattgtcatctacagctatctgcactcctcaattattatttgtctgtactataaagcactacactccaaattacttaggttctttatcatctaggtccattctcagtatgctttcaatgtaccaattgctgctcataagcaaacccatttcttgccaatttaaagccttcagtctctcaccatggagccttttccttttttctacaatctataaaattactttatctcccctagagggtattatggtaagtgaaatagtcagccagagaaagacaaataccatatgatttcacttatacatggaatctaaagaacaatacaaacaagcaaacaaagcagaaacagactcatagacacagagaacagactaatggttgccagaggagagggggctgggagactgggtgaaaaaggtgaagggattgaaaagtacagattgctagttacacaaaacagttaagggaatataaagtatagcataggaaatatagtccaaaatattgtaataactatgattggtgacaggtgagtactagaaatatcagggggaacattgtaaagtatatgattgtctagtccagccgtgggcaaactacagcccgcgggccagatccggccagtttgaaatgaataaaactaaaataaaataaaaaaaaagaccgtacccttttatgtaatgatgtttactttgaatttatattagttcacacaaacactccatccatgcttttgttccggcactctggtacagtttaagaacccattgtggccctcgagtcaaaaagtttgtccacccctggtctagtcactatgctatacacctgaaactcatacaaaataatactgaatgtaagctgtacttaaaaataaacaaataattttttaaaactgcttatctcttctttctttaaaataacaatcatcgttcacttcttccacaatgtctctcttcttcacttatccattcatcctgagcacctgtgtcatggcaaaccctgtgcacattgtagattctatgttcttcagctaaatggatacatataaacacatatatgtccatctcccaaacacccacacaacagcattgtactgatcaatgttaagggttacgtgaatttcagggaggagtgaatatctgataaggaccctgtagcctgttaggagcatttacaattattgataaaaagatactagaagcaaggttagggggaaagaaggacaaacggaaaaatattcgcattatatcgaaaccaaaggattcatcatctctctgaaaaaatgattttaaaagaacaaaccagtactaaaagatgcaaagggtatataaacaggcaattcacagaatacaaatggccctaacaagtatagaatgttctttatattcctaatcaaaaaattgaaaattaaaataataatatgatatgcttgtttatctgttttgtacaaatttaaaaaaaatagctctcacctccaggattgacaagggcacaccttgctgtaagaatattaactgaaatagatagcttttctggagaacaatctggcaagacaagtaaaaaataaataatatgcgtatgttttgacccagatagtctaattctatgtacagggtagggcaaaagtaggtttatagttgtaagtatgcacaacagagtttattcttatattattatttattaattattgtattgttttccatacgagcaactataaacctacttttgccccaccctgtatagtatcctaggaaaacaaaagatgtttaagatgtatgtatacatttcatgtgtatatgtagaaagagtttaattgcagaattaaaatagtaaaaaagaaaataagaaagaattacaaatctactgataagggatacattaaatatagtcatttggaatattatgtagcagttaaaagagacattatagaattttatgtactgacatgtaatactacattattgcaacaattaaatacacactacagaataacatgaatatattgaaacatttcatattagaatcagctcagtgatactacatggataccactgaatctactctagattatttttatattattcccaaatatgtcaagattacatctgagaatttgtgattatatgaccttatctggctcttccagctctaaaattctttaattctaagaactttgtgcagcatagatcttggacagtgttttcataaggccctaaagggtctcctatataggtttaagtgtcccttatcaaaatgcaagaggaaaatgctaaaaagtgtccagagttttacctctacatctgaaaacttgagttcttaatactacatccttcattatcaatgacaattcatgctgaatatgttcattttacagtaatcacgaattataactctttacctcaataataatgtaagaaccattatatcttggaaaaactaagataatagaaccctttgaatatttaaagagattaaatagcccatctcccttcctctctgaaatctataataaaaatatatttttctgaataaataaataaaagataaataacagtaaattttgcttatatggtactcttaatctaaaaatgctgactgtcaatcaaatcccactgagggcaatcttagacagcaacggcacctcttcatactgggcaggattgagatttgtagaaagaacacccgggaatgaaaaccactctttgtgctccagaatccccagcaaatcacatgtaccatctgagcttcaactgtctcatttataaaacaaagaatcctatataataaaagcataatatgaaaatcgaccaaacgtctgaacagcggaatgaccatccggacaaccatccgggaccatgctatgacacccactggcgccaggccagccaaggtaggtgcgatgcgattggtcggaggcccggccatcgccccattatcgccccgtggaggaagacccaggccaccggcctgcggggcgatcaatagggggctccacatcaccccaaagagggaggccaaggccaccgactggtgggctgtggcaggtggacagggcctccctctgcgaaggaagcccaggtcctggtgccggcagccaagggaagaaaggcctactcttgtacaaatttcatgtatcggaactctagtagcataataatgtccagccagcatggctcagcgattaagcctctactatgaaccaagaggtcacagttcgattcccggttagggcacatgaccaggttgtgggctcaatccccagtagggggtgtgcaggaggcaactaatccatgattctctctcatcattgatgtttctctctctctccctctcccttcctctctaaaatcaataaaaatattttttttaaagaatagtgtaataaaaatattatatatatatgtcacatagaaggtttctattactgtgcctggggcaaaataggttctcaaatgtttgtctaatggaaacaattgcatcatttttgtttttgttaaggtatctagtctgttcccttgaatttttataactctcacgcacaaccattctatttaaattttacgctgttagaaagggacttcttgcctctcaaaagcttatttgttgacatactcttttccactataaatatgcaaaaatgtgtttagaaacaccgggagaaagcatacttacaaagcacacccccaaatcgccgtgttgctctccactgagctcagtgctcttcttaggcacatatgtctgcctcctatcatttcttcatgaggttcagggcaaaaggcctagtcaagtagataatctttggttgcctctatattttccccaaaccacctacagataaatgaaacaaggattaaaacacacctacaactttaaatttgttaaagaaataaaagatgaaaataagctcaccctcactatatacagttcatggcatatcagagaacagagaaaatgctccttgggtcatgcaaaagccccactgaactcaaaaaaactctccacatccaagtgcccaaccaacatactcctcacttgttacttctgcccatgcctggtcatccctagggctggccagatcatgcaattgtagcccaagtccacctaaaaggactctctctttcagcatgtctatcagagggtaaactaagtcaaataattaatttaacatgatacttctgaatttaatgacccagttgaatttctaagccagttaacgaactatatatgctagtacgcctaaacgatatatttcgtaagctacccagtaccacacacatttagtactgcaattaaataagtaaccataaaatgtctttggttgtaagaaagtctgtggacagaaatattttcaattctttgtgactgcttttaaagcagtctgtaggaaatcacaggaagtatccagagaactagtagttgtcacataagaagacctgagactacatctataaaccacttcaagaaaaagtaaaaaactatcctggaaaagttaatcaaaataaactgtggagagctccaaacaacacaggtcttgcaatggcagagcggacacgctgagagcactcaggttttatgcatacttcactcttctaacagagcataaccttcaccccatgggctgggtcacagacacatcccaatttctctgtgccaatatgaagctgctaaagatttcaccttgaatctacactgtgcaatacagtagtcactccatacgtggctactgagctcttgaaatctgcccagtttgagctcacaggtgctgtaaatatacaatgcacattttggacacttagaaacattaaaaacatatatacatttaaggctgtaaaatgcctcaatcaataagtttccatactgattacatgttgaaatgaaaacattttgtatataaataaattgtgctaaataaactattaaaattaacatcaccagtttcctttttactggtacaatgtggatactggaaaatttaaacgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtctcacataatatttctacttgacggcactgttctgcataattagtaaatataaataacaaatactgagtacttccaatgtgtcaggccctgtgctagatactttgcttcatacttggaacaacgttacagtgcggtcctacgaggagcctgcggaagccactgggccccatcaaaccgatcacagttgcagctggaattcaatgtcagatgttcctacgaaataactcatatcacgggcagctaacgttggggcaattcccatggtggcacaggagaggccgatacagcgatacagcgtctgtgactagagaacaggtacaaaacggccagcgctatcccagctctatcccaggccaggcgagtcacctccaaaatgagtccgatctccttgacagctcaggacatgccgctcgcaaaaggaggaccaggaagttcaaattccctgaatagcttggcacagaccaggctgtggatgcccatgcggagggggtcgcaggccgggcaccccgtctctccccaccctggcgctcacactaccaagccgcggccccaccccacatctccgcgcggccagtgtccccacccctgggacaaggagaagtgggacttaaaacgtctacaggacagggaagtggcaaagaaagggccacgggcggtggcgggcaattacctgagattgagagcccaggtccttaatgccggatcatcgccgaggccccacccggctcccacggccacctgtgaccgagaaggtgcaactgacccgacgcggtggcgcaacgaactgagctgtgcccagactggcagctatggcaacttgaactagccgcccagccgcggcaagagtacaggaccctactggtgccgcttggccctctgggagttgtagtccgcatcgtcctctaaagcaatgagttaactaagggactctgcccctaaaagaaactgcaagtcccaggatgccaaacgattgtcttcggtgactttcccggaagaagcaggattcgtcggcacctgaggccctaggacgcaccctctttccctgttttgtttttttaaagggccagaagcgtatattgaagacaaacaccagctagagagagagctgccttggcccagactgtaatttaatcaggaaagtctaagatatttcctttacttaaaaaaaaaatcaaagtagttgcaatactaacatagactatatttttgttattccttggggaacaggcaggtgtaatttccagcagagaaatgttaaaatccccttcgtgtaataataccagtgcccaatgtcttgtggttcacaaattcttagatccggaaggaaacagaagattatctagtcacaccttttaggcctggaaaaagtaagctcactgatgacaaattaaaagtggaatttgagaatgtcacttcaaaactgctgccagggtgtcatctctccgcttattgtggaacctgaaagagacatagtccgctcagagaagtcggtttcaattgcgaaaaagtaactttaaaaacttcaacaaggaactgagctttggaaaagctaagaagaagtgaaaatccaaattcgtcacgaagatggtcctttccttctgctgttgatgctttctttcctcctcaaggctaaagtcaagttttactttatctgtgaagccttcctacacagataaagcaagatcactttaacgacatggattccctatcctaatcattgtgtgctgtaattatccattccctgttggtacttacattagattctaagcctttcaaacacagggacctaacctgttaaaatttccatcctcccttctctctctctctgtctctctctctctctctctctctctgcctctctctctctctctctccctccctctcacacacacacattaataattatgatctattatatcatatttgatatatcagaagcattaaaatatttttcctggctggccaatactgtgaggtagatgttattttataaccggtttgaacataagatgattattctaaggacatacctatgtatttgcatctgtcttcacgttgcttgctcacagctgtattcaagaatacactcatgctcgcttcggcagcacatatactacaattggaatgatacagagaagattagcatggcccctgctcgaggatgacatgcaaatgcatgaagcgttccatattaagaaagaaaagaaaaaaagaaaagggaagggaagggaagggaagggaagggaagggaagggaaggggaagggaagggaaggaagggaagggaggagagaaaagaaagaggaaaggaaaggaagggaaggaagggaaaaggaaggaaaaggaagaaaagaaagaaaaaagaaagaaaaagaaaagaaaagaaaagaaaagaaaaagaaaagaaaagaaagaagaaaagaaaagagaatacactgcaccggacacaaatgttaacataaataggaacagtcattaaattgccccatgcacatgcagaagcattagagtacctttatttgtttacattttaatttttcaaaacaattgtactaaaaaaaacaattgtactcattcatttaattcagcaaatgtatgctgaatctttctaggatccagacattatattagacacttgaaatgtaatgaaaatgagaataataacaacctagtttctgccctcatgttgctcagagtcaagaaaatatgaagaagagcagtgccaagagaaaagggattcctgagtgggaatgtgataagaaacagagaaaagagtcccagaatggagagttggagatttcggggctgggggtgcaatgcagacctgaaggtggtggaggaaagcccatgtccatagactgctcagcttctactcaagcaagtactattcatagtctggacctcgctgtccgcgcattggcgtgaagagaagacagcctgcattttacccagtggttgtccagcacaagtagtaaccatctgttccattccacattgttccatctagccactagatggtcgagcctccttgtaatatcaccagacaattataagggataagaaggaagacagagctccgaagtgggatgctatgggtaaggttaaagagaaatctccaacatcctaattggaactacccaaaattctctctgaaccaaagcacaaaatgatcattaatgattgttaactcttcctgcatcaacttttcaagctaactgaaaatctgggaaagttaaaatagaagacaagatgttgcaaagatggtattgccaactactcaaggaatcatcttcagctattcatgtgtgccattgatcccataatattcttttcagccacacacttattaattataggtatatgtcactgtcaagagcatcacaaatgtaaaacacagcaatgtgtttgcatttaagtcaaaatgaataactattaatgcacaaatccaaaatcataaaatatttttattatagctaactcagatagcctgtgctccactttaataggtacctttgtgtgtggtcggggggggggggggtggataaacactagacagtaatgccatgctaaatactactttctgatatctggaccatgctaacaattcttctaacaagaaggaaaatgtaatgttcttaagccatttctcattcattctgtttatatttaaaagcatttttgtggttaatggatcagaatgtaaacacttcttctcttacttaaactggttgatctctcatcccatgcttttattggccccagtttcaaaagttaggctaaggaaagaaaagtgcctactgattatcatagcaatttgttgagaaggttattaaatattagaggcccagtgcatgaatttgtgtacaggtggggtccagccagcccagccccaataggggcagattggtccaggcctgtagggaagaggagatggcaggaggtttggccagccagttcgccccaatcagcagatcagggtcgggccagctggggggagggaccgtgggtgattggccagcaggccctgcccccgattgggggggggggggggggttgggcaatcaggggtcattggccaggggaggggccatgggagtttggtctgcaagccccacccccgattggtgtgtggggtgtgatcggggacagagctggctggggggaggggcctcaggccgattggggtggtggggacccattgggggtgaggccagccagggggagcggcttcatgcagttggctggctggccctgccccctattagggtagagggggaaatcaggggtggggttggccagggggaggggctgcaggagattggtcagctggccccacccccagatcaggccggttctctggccacagtgggcgtcatagcggctggtcgttctggtagttacagtgtttcagttactggctttttatatatagagatggttaccctgaagatttggtccactttcctgaaatgaaacaggtaagcatggaagattttacctacattatgaccaccactgagccaatttcatgttgcctaagtcttgattaaaatgttacagatttttggcataccacttatgttcccagcaaagatggaatagcatagaagcagagagcagactttccctgcacccagaatttttaaaaagtaaaaagacaccactagccagcatggttgagccagttgacctatgaaccaggaagtcaaggttcaattccaatcagggcacatgcctgggttttgggcttgatccctaatggggagcatgcaggagacagccaatcaaggattctctcatcattgatgtttctatctctctccctctcccttcctctctgaaatcaacaaaaaaatatttttaaaaaagaaaaaggaaacatacaaaacatatgaaacaatgctttttcagacaatagacattatataaaggccagtgatttcctaacagacaggaaagaaacaacatgagccttgtgatccccacaacttcctaccttgagagaatttccaggccatggctcaaggagaaggaacccaggcagagctccgtggatgctttgagttgacaagacagagatgggagtttgtgaaaaccagtgtagatggacaaagtagcaaagagaagagagctacatacagagaaaactctaaagaactagggaacgtcttcttggaatattcagatgaatacgtgaagaaacggaggctggagaaggaagcatctaccaagactggggaggagaaagtacctgtcactcacacagggccaggaatgatgattatttccaccagccggatggagacacctcttaattcacgtagtattatgtagagtgctcagaagagtcttgcctcaggagggaggagtaattggtgctaaattaaatgctattcccatcccacctcacaagccttgaaagcaagacccaaaaatatcaaactatttccaagtaattttcccacgtctgataacaaatctcaaaaatatttctatgaaaacaaataatccagtacccaaagaggtgaaatttacaatgtctacaatccaaataaagattaccagtcatgcaaagaagaagaaaaatgcaactcattttaaggagaaaaatcaaccaattgaaacttacccagaactgataaagatattagagttatcaaacagggacattgttaggggtagaattgtgtcccctgaaaagataagttgaagttctaatacctggattgtgactttatttggaaacagggtgtttgcagatatgcaagttaaattgaggtcatgctgaattaatatgggccctaatccagtgactgatgctttcatggaaagagggaatttagagcagagacacaagagagaatactgtgtgatgacagaggagagactggggtgatgcgaaccaagggcaccaaggactgctggcaaccatcaggatctgggaagaggcagggaaggaatcgcccttagaaccttcagagggagccatgaccctgccaacagcttaattttggacttccagccaccagaactgtgagtgaattcatttctgttgctttcagccatcaagtatgtggtactttgttacagtagctctaggagtctaacacagacatcaaaagtgttttagatgtcaaaagttacatagagacaggggaaatatgaaaaaggcccaaattgaacttctagacatgaacatttcacctgacataaaaaatatactgatggaatcagtggcagatgaactatttcagaaaaaaaaaaaaaaagaatattgaatttgaagatataaaaaagaaactatgcaaaatgaaaaagagagaaaattaatttgtaatttaaaaaagcatcaatgaaacaaaaaacaaaaacaaaaaaagcataaatgagcttttttgtagaacaactttaattagcctaatatacgtgtaataattggaaagagtagagagacataaggggacataaaaataaatggagaaataatggcagaaatttttccaaaagtgacgaaaaacagaaattcaccctagcagtatagttgagttggttggaacatcgtccagtacaccaaaagatagtgggctcgattcccagtcagggcacatacccaggttgtggattcagtccccaatagaggttgtgcaggaggaagctgatggatatttccctctcattgatgtttctctctctctctccacctcccttcctctctctctaaaatcaataaaacatttttaaaaacaaacataaatgcaaaaattcaagaacctcaaaatactgcagaagaaatatgaagaaagtacatcaaagcatatgataatcaaaatacttccaaccagcgataactagtcttaaaagcattctgagggagaaaaaagcatgttgcaaaagaaaataatgacagaacatttgtccttggaaataataagatagtggagcaatatttttaaagtactgaaaggaaaaaaatactttcaatctataactctgtgcctagtaatgctatccttgaaaaatgaaggcaaaattaatacttgttcagacacaaaaaagatggaataatgcatacccaacagatctgtaccacaagaaatgttaacaaaagtgttcaaagtaaaagcaagattatgtatttttttaattgattttagagagaggaaggggggagtgagattttatatatatatatatatatatatatatatagagagagagagagagagagagagagagagagagagaaacgttgatgtgagagagaaacatcaatcagttgcctcctacatgtgccccaaccagggattgaacccacaacatagatatgtgccctgacttggaatcaaatccacctttcagtgtactgaatgatgctccaaccaactcagccacactgaccagggctgttttttgttttgtttttttattgaaaaactagaggcccggtgcatgaaattcgtgcatacatacggttccttgtcctggcctgcgatcaaggccatcttccccagctgcccacagccggccctgcccccgccgccaccaacccccggttccctgtttgccattgaatgatcggtgcctgacggccaggggaagggaccgagaggttggctgttcctgcccactcgccagccccacccctgccacggccatcaggccatcagagcctgctggccagggaaagggaccaagaggtggtcagtgcacctcatagtgactggtccagcggtcattctggtcattttgccattagggtcaatttgcacattacccttttattatataggatgattatattggataatattccactgtatggatgtatcactgtttgtttacccattgcctagttgaggaacacttgggttgtttccagtttgggttaattatgaataaagccctaaaaccattgaaatatacaaattttaaacatacagtccaatgaactttgaaaaaaattatgcattctcttacctccagtctaattaagatatagaacatttctataaccctggtgagtgctctcctatcttttttgaatcagtctcccctaaatcgtgttctgatttccatcactatagttccccccattctgtcccagaaagtcatataaatggaattacacagtacgtgaagaattgatgttcttttcatcccaagaagaattagtaactgaaggcttttaagaaggaaatggcacaatctgatttatattttaaattattatttttggctactgtgtgaagaatggtaagtagatgatggtgtgccacttaagtgttatttttcttcccatctttggggggattttggggaaagcggagctatatgatgcaataacattcctggcagtccagaatttacctttaggcttcaagccaggtcggaaggccggctccagagcaccagacctgttcctgccattcctggttctgcccagcagatggccccagagcctagttcccggcccgcctagagccaagtccgcctgctgcaatgtgagctgacgcctccaggtgactcctcaatacaatttctgtgcaaagcggtcgggtttgcgctgaagttggggagcagtgagggtccgccagaggtattgcgccccggcttctccgcacctcaggctggcatcctcagaggcaccttaggccagaaggatcttgacccgcgtctgccctggacggtgcattttgcgggaagttttctctgcagtgcgtgcgagccttgggggagaagggggagcgggggagaggggtcctgacaggcagcgcttttcggaatctacccacgggaaaaggtggcccactggacttaacgacaagggaagggcggcgggggggggggcgctctccggaacatatggctctggaaagagctaatggaaggaacgcgcgaggctggctggaaagcagggtggggacagcagggactggtgttaagtgtcagggaacagacagagcgggctgcgtggctcagggctaggaaagtaaagccaatgggcttgtttgcgtggagagttgcctagggaaacagaacacttgagagcctctttgccctctccagcctcctctttggtgctgaagtcacaaccaccaggagtcctctttctcccactcctccctcgctccttccgcgcggggagctcggagcacagagggctgagaatgaggcgatagtggaggaagaagaaatggcctcggacccctggaaaatgaggccgatgcccctactgcagctggcgctgcttctcggcctgcccaggagcttaggggggaaagggtgtgcttctccgccctgtgagtgccaccaggaggacaacttcagagtcacctgcaaggatatccatcgcatcccagtcctaccgcccgggacacagactctgtgagtacctgggagaggagagggtaggtcccagaggccaagggcagccggaaaggtggacggaagtgcacaaaaagaacttgagaacaggccaaagggtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgcgcgcgcacgcgcgcgcaatgaataaaaactcaatcccataaatttgggattgctgacctcatagttatttaataacacaggaaacggtaacaaaactttggtacaagatgagctattctgcatattccctgctcctcgcccaggctaaatcttccttttacctcaggtatgaacaagattagtagcaattcaagaactgacttttctatttgatattgggtaataaaaatttgcaaaaaaaatgtatttctccttttagtcttccaaaatgggaaattatcttagctggtcttgctggttcatcaaatagttgaaaaaatttagtcaggtattcttttttaaaatataccctttatattttccatgaaggtaactttgtatgaagtacctgagttactttctcaaaggttaactgcccactaacatagtggatccagaagtttcctttgccataggcatgggagaggaataaagaggtaaaatttcataaaatgcattttatctgcatattaaatgatgccaaaaattgtcagatggttattagccttattggggtgttatcacttcatcctatataataaaagcctaatatgcaaattgaccaaacagcagaatgaccagcgggacgaccggttgctatgatgtgcactgaccacaggggggcagacgctcaacgcaggagctgcccccagcctgcaggtcccaagccagcctaggcaggtgccagtgggggatccccccgatcatcctgtgggtcaccctgcagagggaggcaactggctgcagcgacggggtgatcaggggatggggcgatcaggggatggggctagctggtgagcggcaccaggccggccaagatgggtgccagcgggggcccccccgatcaccccaccagtcgccccacagattggtcctgatcactggcctggcctagggaccctacccatgcatgcattttgtgcaccgggcttctagtaagttatataaatgtctaatcactatgctgtacacctgaacctcatgtaatatgtatgtcaactgtaatagataaaatttatttttaaaagttgccaaattaaaaacaagggtcaacaagttctattgtatgcttattccttaagaggcatttaaatcaatgaccattgaaaggataattacagacaaaagccaggaaaaaggctttcttgggtggttaaaaatggagcacaaagaaaaaaaaaaatagtaaaacaacctaggctggttaaacctagttcaggtaagataatacaagttattaatgtagatacatgaaaagcaatgatgaagttgcacttcaggattttctcaggactgtagctgatagggacattaaatagtggtgaatgttgtttggtgattaatgttgcagtacaggttctgtgctctgaggggagtgcaggaagctgctgatggatgattctctcatcattgatctttctatctctctctccctcccccttcctttctctgaaatcaatataagcatatattttttaaagatacagattaggggtccgttctttccagaagcctctcttgatctatgacccctagactgggtgagctccccccccccccccaacattctgtgcacaacctctatcctagtccttaaacatggcatgcttttctgtaacttatttctttgcatgcctaccaagtaggacacgagctcctctaacacggaaaccctgtcattcaaccttagaactactttgatgttcaaaaaaaagtttaagtgagtgaatgaccgaacataaatatgctagtctcattattcgctgtgaggataaagaaggccagagagatcacttcttcggggacttctttaaaggtggcagaattccaagtcattgaagagaaatgaggctctgacagtggaagtgaatttggattcatcctgtcaattccaccaaatgatgaggagaaagggtaacatcagagaaagagcgatggaggaggagtcaggaaacctgcattctctttacccccctctcactgatggtgaaattaggtgagatcatggatttagaagaacccaagtttccttctgggtaaaacagagggtttagttaggtgctttctgatgtctgtgatgctcaagagggaatccactgagccacaaaggctaagaactggatagagacaggagcaaggagaatctgaccttaaaaacaattagggggaggaaggcggtgggactattctgaagaattaagttattctgagatagagaagagaacatttaaaagtatgggtaaagtggtatgctaatttaaataatgactacaaccatctacaacagtagcaaactgatttttgaaaattttgtggggcaatagaactctactgccccacaaagttttgaaaagtcacagtggattcccttctgtatgtcctctttcaaaaagtgtctatttaggtcctttgcccattttttgattggattgtttatcttccttttgttaagttatatgaattccttttgttaagttgtatgaattctttgtaaattttggagattaaacccttatctgaaatagcattggcaaatatgttctcccatgcagtgggctttcttgttgttttgttgatagtttgttttgctgtgcaaaagctttttattttgatgtagtcctatttgtttatattctccttagtctccattgccctaggagctgtgtcggtaaaggaattgctacagattcaatgcaatccccactaaaatatttcaaaggtctagaacaaactctccaaaaattcatctggaataaaaaaagacaccaaatagccacagtaatcctgagaaagaaaaacaaagttggagggatcacaaaaccagatatcaagctatattaccaagccacggttctcaaaactgcctggtactggcacaagaacagacatacagaccagtggaacagaacagagaacccagaaatttacccaagccattattatgctcaattaatatttgacaagagaggcaagagcatacaatggagtcaagacagtctcttcaataaatggtgttgggaaatttggtcagatacatgcaaaaaaaaaaaagaaaagaaagaaactagatgaccaacttacactatacacaaaaacaaactcaaaatggctaaaggattttcatgtaagatgggaaaccataaaaatcatacaagactttacaggcagcaaaatggcatctattttcttacttcagctctgatagcatagtaatacacataaaacatattattaagtgagtgaattttatctgctttatcaatcaattgaaaaagttttttccctccagagctgagaaattgttaacttgaattgatctatttgattattttaaaatttccttgttgtacaaagtatgtttattacccttcagatgaagattaattacctggaggttcatgggcttcgtgggggtccattcagtcctcctgggaatctgtgatgaaatgaaattttgaatacaagtaaactgctgaatctataaccttgctcagatttccagaggtgcccttaaccccacccccaaaaaagtcaagaaattgtgttttgggagaatttggcatttttgtttaataactgtatatttggtattaaatgtgactcagtccttttattactagaggcccagtgcatgaaatttgtgcactcggggtgggaggtcctctcagctgggattgcaagacagtgcaggccaggctgagggacccaaccagtgcacaatcggggtcagggagggatgcgggaggttggccagccagggagggaccctgggagggctccagggcgtgtcggtcccatctctctcagtcctgatccgccagaccccagcagccagctaacctaccagtcaaagagtctgctcccctggtggtcagtgcatgtcatagcaagcggttgagcagccttagcatatcactagcatattacgttttgattggttgaatggcaaactggatgaccagacacttagcatattaggcttttattatataggataaagaacactaaaacatatttttaaataaataaaaatttgcctgaccttcggatgtgatattttttctctgaaacactgtattacctcgtgatgtcttgcacttagcccttaaataggctcttcatgcagaagttcactaaatgtattttaaagtagtttacaaacagattttacaagtttctgatgcctcgtatctgacattaaatgctccatactggccctgtttccattttcactgacgctgggtgggatgcttaattaccattttgcctttcattggtgcaattaaaacacagttgaaggaaatgtatgggggaaatggtcccatgtgttcaaagctgcaaattggatttttaaagggttccatttaggtataagactctacttactttgcttgactttggttgaatccttaagaacattcagggtgtttcaaaaatttttttcttaaaatgacaaagcagcagtcgtgatttttaaaccatcttataaatgttgttctcttacatgtagtgacatgtaaatatgccctggatttttctgctccagccctatttaaaatattttaccctcaggatagtccagaaactcaaatatttcatcatcatcttaactcaattttccacacccctcttgcacctgggaagagttccctgtgtcccagtttttctttgaaaaacctgatccctgtatttatatgatatgtgaagccttgggcttacatcttggaggcattagcttaattttacccactcctgagtgttattgtaattactccagagcagggattgaaagcgtctcttcttactgccagctactgtcactccgtattttgttatttccatcttttgtgaacgggaaaggagaaacctgacctcagggaaatgatctgtacaaaacaatgaacataacagtgttgagtgggctttgtgactgctctgtttccttctcaacactcttctcagttatttctggctacgtggccatccaaggatttgagatgaccagaaaccttacctctggacaggtggcattttattgcagcttttcttccttcctaccattcttcttttttttcaataaatatttattgatcacctatatacacaggtactattttaggcattggtgatatctaatacatgacctcattgatctgacattctagcagggggaaatagccagataataataatactaataataataatacacatagaatgtcagatggtgataagtataatggagaaaaattaagggagagggagaaggaatgctggaagttgcattttaaatggtccgagaaggcctcactgagcaggtgacatttgaatatcaatcttgggagagaaattctgaataaatgcaaagacctcgaggtaggagtacatgggaggagcatgcatgtatattttagaagtaacaagaagagaaatgtgtatatttaagaagtaataaagagagaaatgtggctgaaccagagagggttatggaaaatgagatcagagaagttataagacagatgatgtactttgtagacactataagacattaaatatgtatatgttctaaatgggataagaagcttaataggatatttagttcaaaaaaaatcacatgaactgatttacattttattaggatcaccctggcaactgtgtttgttaatcagtaagtgatgacaatggcttggaccagcttaacagcggggaaggaggtgaaatatggttgaatcttggattttatgatgaaggcaaagccaacaggatttcctaacagactgtgtggagtctgagtgaaatcaagagttaatgtcttagcctgggttctccagaaatcagagcctgtgtgacaattcttactggggtgtgtgatttgggggagcagaacggtaggatgggggaatgacattggcagggagggagagtcggctgctgcttcgaagcaaatagttgctgatcccatgggactccctgtgaattccgatgctatgcatctcaggacccttgctctaggagggaggaaagggaaattacatccatcagcttctagtcccctgggtcagggtgttaattgcaccaaccttccatgtttcttgtgtataaatgaagaaagagttcccagatagtgacagaaaagctccagggtggaaggtaagaccctgcacttgactctccagtcactccccatcccctctctgcccattgcccaagcaaccactaatctgcttctgtctctatgaatttgtctattctggacagtttatataaatggaatcatatcacatgtgacctctgtgcctggctttttttacatacttagcagtgttttcatgtttcatgcacactgtagggatggaaagggtgtgggcagccagacaagtggcctgaggggaaaacacaggtcaagaggccaagctttaggtactttcaggatgtttcaggcagattcatttctgacacactggagtagagagcagtcccactccatgtctaaatgccttacaccttatatctaaccatctgaacagcaagaatgcccatgcttattctcagaggcttagaaccatgaaattgcttccccatgtcacaggttttaagagtccccggccagtgaagtgtttctagccagttggctgaggaacaccaaccgcttatcccttttgcatcccatggcaatcttcgtcaacacttctagaaagaggtgacttcatagcaatgcatgagtcattgtccgttaactctttaccctttaaagagccaatacttatgcaaaggacttttctctgttctagggaagctgcagaactcaaacagtatcctttttgtaatgggttttctcataataacattggcttctatgtcttattgtgtttctccattaaaaaatttaatcaacatgctaataaagaaaatttgctatgaatttttttaaaaaaggttcaatgagagaaatagagtagtcttccccttctgtttaccagagcttacttctcagctatttctttccagcccactctgggcttgggtttctaggaagatgtgcattcaaattctagctctccaccaactcactatgtgggagaccatttgctcaactcacaagatgaacaaactaatgagatctgctgcagctacgtcacagagtggttttgaggatcacgcgagataatgctagcagctctacatgaaatacctggtaaaccagaataaagctttgatcacaagttaactattatcactaatgcagtagtgggtctgagacccagacccggcagtgtgggcattacatcagtttgtagaggtcaagctctcacgcagggaccattcagcagctgttaggcttgagggctgtgcagccttagtaagaagagcaaacacatgtcctccatgcaaaagggacattttaaagatttacaggatttcttgtcttctgccagagatacttatccctggggagttactgtttgtagtgtaaagtggccgaaattctgactggaccctattgggtatttctgtcccaatagggtccataactctgtaagaaggattttatctgttcaatgtcagattttggctgaaaagcctgttagtatttccaaactctccagatgtccctacttaatgcttggctggtcttcttaggtagcataaccaagttgccatttcacaaccagcctccatccccaggatctaaaacagcaaatccagaagctctgaaacgtgctctatgaagcatgaaaggatccgttcaggttctgtctgattcagttctgtgagataatagacattctaacgtacctatgaattccctttattcataagaaacttggcactttccagaaagtatgtctctcatgaccaagaaaagtagggtgacttggggatgcaatgctttctttttatagtgtgtctagtatcttcattttcaatagcagcagttggtctcaggaccacaattggtacttaccatctctatcctccactgactattctaaattccccagctctcatcccttctggccttgatgacttatctgatatgatccaagatctcattcctcagagatctaggcctgtggtaaccatactcttctcacactagggttgctagaattgactattttagtcataattgagcaagtactatgacatggccaggtggattcctgaattctatacatatttccccctgccatcattatataacaaatgtctactacctcttcctggtgatcagaatcagttactcccatcaattactgccatcctatgtccctgaacataaggaatccaaaatgtccaggtaggaatcatagcttgtagttccctgcagagcctagaattgtggggatgggaagcacagaattccgtggtgggtcattgaccataatggcaactgggaccacttctgcttctgtcccttggttcctagacctatatattcttgctaacaggacacaacatccaaatagtttgtaaatcattctgtaatctgcatgatggtatcctatcatggcagaagtattacctacaagttagcacttccatgtgcctttggcaggtgactttattgctctgtaaaaccaattgcttctggattgtgtggcacttgatatggccatgggcttacctccctatctccttctctgtgaagtgagtctcatgatcagatggtatgttatatgagactgtatgcctagatcagatactttgaaatctcccagatagtggtgctagctgaggctctaggctgccaagaataagtgcctactactataacaataaatcactggcccttccagtaagaaaaggctcaaagtggtcatcttgccaccaagtggccatgtggtcccctcaagtaatagtgttacattgatctctcttgcaaactggttggactctcaatgatagtagtagtagatcagccttgataagtaggggagtccatgctattgacatgtaatgttcattctgttaatgtgccaatcatgctaattctggagcggccaatgacaaaagatggctactgtcaactaatcaattttgtcttggctattcaatgcttcttctgtggtggccaagctttggtggtcattaatatatgatacgaagatcttcacatttcatacccactctcatatattcatccacaagcatctagcccagatttttttctctgcaaaattattaactttttcttcccaggcccctgactagctgggcaatccttttgctcacgattccatatgtagtcttatctcaggcaacttttccatgaaatcctcccattaggaagattttctcttgctgctgtctttcaaggccacatgcagccaaattcagtatgttcattttaggtggacatattgaaaaaacatctacataccaaatcaaatgacccattaaccaagttcagttatactggagcaaatgctgagacagagttgggaatgcatgtggtttgcaggaggtaacatctgtgaaagataaagagagaggaagcaggacttggtagatacagcctcacattctctgcaaattggacaaagtcttggccaacccagtggggaactcttagtgaagacggcttgtgagaggggccccacattgggcataacagccaggccctagtgcccctgccatgctcatccattagctgtgggcttaccaggaagattagggccttggttcaaacatgtgacagatcccaaatgtgctacagatagaggctggtagcaaatgtcttcttgacagttaaacagaaagctgtttcttggagggagatcagagtggcacatttccagagctgtcacatcctttcttttttttttttatcttttcagtttgcacttccttggatttttctttacaagcttcttttcctctacctaatacttacttgatgtttcttaggttctcttttcaggccacacttccaccagccaccttgtcaactcccattcacacatcaaatgagactctcatgaatttccatttcagatctttaggctaggtctcttctgagcactaggtctatattatccatctgcctataaagatagcttacagattcctaacatttagcatttctaaaagtgacctcattatcttcttagccaaaactctttctgtgtagggatctcagtgaatagttcacattttcccagatcctcaaaacaaacagcatagagttatccttgatactgtcttctctctcattcctaatatcaaatccataatcattaccaatctaaacatgtcttcaatctatctcgattccattcttgttattgctatcaaagttcacaccaccatcattgtttcccacctgagccactataaatggagtttcttctgtgcttattaaacacagcagccagagtaataattataaaatgtactgtggatgaggtcctcaattccctgtgtaacagatcattggcttctattgctctcaggggagcctaacagtcctgaattcctctttcagagtgatgctgactttgcaattgttttatcttctagaccaggggtgggcaaactttttgactcgagggccacaataggttcttaaactggaccggagggccggaacaaaagcatggatggagtgtttgtgtgaactaatataaattcaaagtaaacatcattacataaaagggtacggtcttttttttttttttttaattttattcatttcaaacgggccggatctggcccgtgggccgtagtttgcccacggttgttctagactattggcacctaaggtcagggattgaatggtattcattcacctttgcttcccttaggatacctaacctagtgccttgcacatagcagattcacaaaattatccaggattgaagtgaaccttcagtgtgctggctagtagatctttttaaccaatagccatcttatcctgggttccctggaaaacagagttaacttgctaatgtgagagtgaagaacaagagaagtaaggcggggggtgggggggggggggccgtgcttccttcccaaagagacaggagaaacccagcaggtaatgtggcaagtgcctctgcacggccacacggggtacctgtggacaggctctgtgaaggccgtgctttgcggcagttcatgagaaggaactcacccacctcctgtttctcattgatccaagttactccctgcacttcctttgtgtcatctggccccttctatggctccctggaaagatggcttccatgccctattgggaggcctttcatttgagtccagagtggaacacaggtgcatgtaagatgtggccaacaaactatccctgaagctcagcctgttgtgacttttactgggttagattttaaatacactgatgctagaaacacaaaactggtggcagttttctttaagagagtggcctgcaacaggaataattgttaaacgttttatacaagtaaatagtgactgtgatgaataagaagattaaaacacagattgaacatgtcctgatttttctctcttaatcacgagtaatttgcaaatccatacatggctttttaaaaaactatatattctcatctgatcttatgattcaaaattaaaagacctataaactcctcttcccatgaaatcactcaatacagtacaaaagtccatcctttagcctatcactctctatccccatatgtaactcatacctacaatgaaagttgctcaccagcacacatgacaaaaaaacactgaagtgggaagggggtgggggggggattcttttctcaagcaagtggttctatccacttacgaaaaagacttgagagcttggacagacttcaagagcctctgaatgtttttcttcagcagaaattccagcccagatgttcaataacccaggtttgttttcaaagaggaaataaagaaacaaatctcataaagacaggccagactttgaaagcaattaatgatgacaggtccaaggattctgcttcattttattattaacgctgccactctttccagattccaatgtcagatactggtgaaggagaagtaaaaaaaaaaaaaaagtctggaaattgaagataaattattagtgggaggcacccacagaaggggaccttttgaagtagtgttcagagtggagaaggagaattgaagggtcagttgactgccaggggatactgtgaagttagcagatgaaaacggttgacaaacaggagcctaaggaaggtctaatcataaactttacactttgaacagatcattgcagatcaaaatcagctgacatattggtggagaagctgctttattttttttttatgtgttttgtttttttaattattccccgaggacatagttttaatgatttttagaaagagaggaagggggagagagagaaacactgattggcaacctctaatacacaccctgactagggatcaaacccgcaacccaagtatgtgccctcactaggagttgaacctacaaccttctggtgcatgggatgatgctccaaccaactgagccacctggtcagggcgagaagctactttaaattcatttacaccaacaagcccagagagtgataaagcaatcaactcataaatctcagaagtcactatgtggaatttccagggaatttggggaaataaatgttagactagaaaggatcacttagttctggttttcaaaaaaattaaagggaaaaggtacttattaaagcactgacctgtatgtgctggtcagtgttttgggttgtgggttgtataaaaagtgtaatagcctttcttcttgagtcaatcacattccagtcaggaaactggacgtaaacttttctctcatgcaatatgttaaagtctgtaacagaggtaaagggcaatgaacatggaaagaagaaaaaagatttaaataattgctactatttgtatagtggtttagaatgtatgaatagcttacattatcttattttatttaattctcaaaacaattctgacggagaatcctcaacttacagaataatcagaagctcagagcagttaagtggcccagctatttgactgtttggatgctctttgcttccccagtgtctactccctttcctgggagatttgaactatttttaaatgaattattgtgcaatacaaatgccacacatgactacagttcttctttattcttggaaactttctcctatgacccaggaaatgactactatttgtttatttcagaattttctgagtcacgttgctctgtcatagtaagctttcagccaaccagcagaagtgccggggtccaaccccagcaggtccaggggtccccaaaggtgtggacggagttggcgaagaaggaaagacacagagacagcgttcagctgatcagcagcccatccaggatctccagccaagttctggtctggatctccagcgaagttctggctaggatctccagccaggttctgtgtccatgttctcttgctaggttctccaggttctccagccaggttctgtagccatgttccctcgctaggttctccagccaggttcagtcaccaggttctagtcaggttctcttgccaatttctgtagtcaggttcagtccaggatcttttgccatgttctctccagtgaggttcttctgtctctagagaacgttctgtgtaggttctgtgcctcctggttctgtctctcttggttctgtcttcttagttctgtgttctaagttctgtcttattctaagttctgtgttctaagttctgtgtcttgctgtctagttacatctgtatttataccagttgattcaatcctgtcaatctctatttcaaaggttagggtgtttcttatctccattccagggagtaaagattatgtagcttaagcatgactgttcatagttaaagtgattaattacccgcctggcacttagttaagaggttttattccctccctaacttcaggggaaaatccctacctggggaaacaccctttctcagagaccttggttaaaacacatagtgccaagaaggtgagcaaacatattaagaacagtatgccatatatgccaggtcccttgaaacagcaagcatggaccggctcccggcacagaaggagcaccaacatacaagaggattttatctcaaaggtttacttctagttgcctaagggggccagaaattcagaccaggtcctccagggaaggaaatgttgtcccagtaattggggaaaacctgggctttaagcaatggcaccactaaacatccaatctcttgtctacgatttcttgggattctgacagcaataaaggagggaaatcgctcttgtactggtttctccttggccctccagtggctcctgttttcccaaatcctaaagaaatagaccttgtggactgaggtgggtggagctggtgggatggagaggaaaagaaaaaagtttcttctcgtccagcttccactcaccagatctgttcttccgaagcaaattcaagacttcacactggtctaccaagtctctgctctgctctctgggtaccctgaaaccaaatctagtaagcttctaaacctcagtttcttcaactggagagaaaattactatttaaattattctgctaccactattatcatacagttagttacttctattttctgtacttaatctcatttaatccacatgtcagccttatgaggcaggttttattgttcccattttaaagatgaggaaactaaaggttcaggaactcagttgaggtcataccgctagaaaggcatggtgctcggaccctcgggccccaaagccccatctattcatctctgttgggaagtattattcttcttagagtgacaccatctcacagttacaccctacgacctttgttatccccatctcataataattccatgaggtgaacaagataaatattgcctctgcccttccctttttacctctttagccccatttaaaaatatatatttttaaaaaaactaaagctgcttggccatggttactcagtggttgagcactatgggctcaatccccagttgggggcatgcaggagtctgctgagcaatgattctctcttatcattgatgtttctctctctctctcccctctcccttcctctctgaaatcaataaagacatatttttaaaaaactaaagctaagcaaatgcaaaggatttttccgaggttacatagataacagatggagttggaattccacaggctttagtcctgctattgctacacacactgctttgtctcctggagatcaagacccagcattgcctctagctagcactttaaccctcaaatggtcatttaaatcctgcagttttgtatatggaaaaaaagatgggaatagcaataataataatatctgtcaactaccttgaagtattattacaaagatcaaatgaaatattatgcataggtgtttatggaataagtagacaatatacaaaccaatataattttattcttaggactatcctgccccatgtgtgtaagggggagaggaattgatgtggaagatataaagagcattggagggagcgggtcatcatttccttgctagatttattgtttcaattggatattttttgagtgctggcactgtgctagctaatggaaatggaatgctgcacaagataatttggccttactctcatgaaattgcattccagatggacggccaatgttgataaaaaaaattagccgaaaccggtttggctcagtggatagagcgtcggcctgcgtgactgaaaggtcccaggtcgatccggtcaagggctgtacctgggtgcgggctacatccccagtaggagatgtgcaggagggcagctgatcgatgtttctctctcatcgatgtttctactctctatcttcttcccttcctctctgtaaaaaaaatcaataaaatatgtaaaaaaaaaaaaaaattacaagtgtggccctagcttggtttggctccatagatagagttgtcagcctgtgactgaagggtcccagatttgggggcacatcctgggttgtgggcttcaatccccagtagggggtgtgcaggagagcagatcatgattctctctcatcattgatctctctctctctctctctctctctctctctctctctctccccccccttcctctctgaaattaaataaaaataaatttaaaaaataattacaagtgtggtgaggttctgtaaagaagacaaatagacttccattgtgactttcttcctgagcacactgcaacgccctggtgctatcacttacaaagatcaggcagaggcaactgccacctattagagcactccaaatggcttccaacaatgaaaaacttacatctcctaaataataactattctaaaaagagactaatcttggctcttagtctctctaaagccaaaataaaggcaatgtgatctatacaacattatcttgacctcactgaccactagtcagtccaagaaacaggcgtgcaatacgatagaaggaaggatcaaagaatatgttttataaggataataaattccgctgacatttccccttctgagctacaataattctggcctcttttataatcataccctgttggaaagccacatggcttggaccagtttttgcctactacagattttggtttccaaaactacatctccatgagtaaggagatgggcagtctcaggttggctggctccatcagaaacacttaaaacttgattgcacaattgtatgaccaactggggtcttagcatccacttgaatttcctcattggcaattttgaggatatcatattgcttgggaaagaaagagggactagttattatgattttcatgaacttacagaactttagtattggaaaagtagaccatctgctttaactttctacccagtgtaagaatctcctgtacttgattcatgatagtcgtcctttgaataaaggaaggaaacttactatctcgttgttgaagaattctgttaacctttcctacactggagaaaatctaccaccattggtcttagatgtgtccactaaatccccacagaccagcggttctcaacctgtgggtcgcaacccctttggcagtcgaacgaccctttcacaggggtcgcctaagaccatcctgcatatcaggcatttacattatgattcataacagtagcaacattacagttatgaagtagcaacgaaaataattttatggttgggtcacaacatgaggaactgtatttaaagggccagaaggttgagaaccactgccatagaccaaatcagaaagaacgtcacttctatttaatgactttcctcctccctggtataagaagtcttttctaagctaaacattccagtgctttcctccattacatatacaaaatgatcctctgaactgttagcattttggtttgagaaatcaaggagatgggtgagaaacagaatgtttgtgtcttttacctaacatatttgctatttctctgttttgtaacatcagcagataccacctgtgcccccgttcacactgggaatgaaactgcacagagctactttaattataaagataaaagcaaaatttctccattcatattgccagaggggatggcttgccaaatcattctctttggttctagaaaatgcccaaaaaataattttgactttaaaattcagctctgggtttgaagatgctgccatgctgttcagtcataagaatatataatttttttatttcgattatgctgctaattttagtaactttggcaaaggttataaaaacaggcagtccccaacttacagggttatgttccagcagtcatttgtttggaatgtgttttacatgcttattaaactttctggctggcctaccaaagactttttaacccctaagatgccaggcatctaataccatcccataatcagctgcaagccagaattcactctcaggttttatttaccatgtcccttaagcatcaccagaaatattctagcacccaacactatgtttaataatgtctttctaacttaatatgcccaagaaactactgggaacacaaggaatattttcatgcttcccttctcctgtccaagtgagttaaggactccctgccccagtcaaggaaactaagagaccgggactccttacgacatgagcctgctggtgtctgggagggtgaagttagtggtagaccttctgggggaaaagctatacttctcagtcatcattggttctactcagatttccattctttccttccaagtcagcctgattcatagactcttcattcatatttattctcaatatctactagtgtcttctctctctctttctttctctctctctctctctctctctctctctctctctctctctctcatctcttccttagattaaatttcttttgcatgtacagtattaattccactgatgtcataataatcctgtaaataagaagtgaggacactaccacctccattttacaaatgagcaaaggcttaaagagggtaaacaacctgccctcacagtatgaactataccagcttcaccacagcactgccttgcttcacttcttaaaaacaaaaggtatatatgttagaaaaattttagttttatatttgttgtctcaagttatttttataaatatgcacagtatttgtagcctctttttttatttctaatgcaattatatactagcattgtgttttttcctgtctttttcctaatcaatcttgctgtgatgtttgtctgttttattaattttttttcaaagagccatcttttacattagtttatcaagtttttaaaattatggtaaaatatttatggcctaagagaaatactaattagaatatttagtcatatatagtaacttctaggaaaacagatgttgctttgtttcaatgacactgaacaatgtcccctggaaaatgtgaattgagttttgttttgtcttccctatgtggcttttgcaaggagaatcttaacaggtatacttgattgtcgatacatttaggtgtttttttctgatgaacatgttaatataggatttattagataagtgataaacattgttgctacaagctaaacacgcaggccctgtcttctaggaattttcaagccagtgaaaatgtatagtctctctttgttgtctactatactctgttgtgagtcaagagcctcctgacactgagcaccgggttcccgggagattaagtgatggacctggacctcccagcccgtgagctcggattccatcctccaaagcttactcttcccaccggaccacacgtatggcccacgggtgacagttactctgcccttttagagcctgagaaaacactgtttggtggtctggcatcatctatgttgacgtggcttgatttggaatagagtgtaatcctagttcgtgatgaaaggatttcatgcgttaaccttatgggtcatgatgtttgtgcaggggggaaatgtgcgctaatctttttttttttttttaatctccaccttaggatatgtttttaattgatttttagagacagaggaagggacacagagagagaaagaaacatcgatgtgagagagaaacatggattggttacctcccagacacaccctgacgggggattaaacatgcatcctgggtatgtgccctgaccaggaatcaaacccacaaccttttggtgttcgggatgacactccatccaactgagccatcagccagggtatagtctttaattcatgattaaaacaaatatacccacgagagatgggcctataaaattcaccttaattgttgaaatgcagtaaataaatctttcaaaaatatctaacaggtccaattatttcagtctaaaaggcaattggagatctgttctttatcaaggttccccagccctcctgtgagatggagctgtttttcatttgcttgcaggctggaagaggagtgaaaacaatgtcaggatggtcctgctcagagccctagaagagctgcctttatggagtgtttgtgtccccaggaggaaatgagatgatctgagtctatgggttagaccagagagaggggctaggtattagcagccaattcatgatcttagaagtctcctgactgaggcacatgaggacaagtatcttaaaaaaaacaaaaaacaaaacaaaacaacaacaaaaaccagaaaacactctgttttatttattttttaaaaaatatgtttttattgatttcagagaggaagggagagggagaaagagaaacatcaatgatgagagagaaccctcgatcactgccttctgcatgccccttactggggtttgagcccgccctgactgggaattgaactgtaacctcctggttcatgggtcaatgctcaatcactgagccacaccagccaggcacaaatatcacaaatatctttctgatggtgttagcatccagaattgtgcaggagtagttgtttctaggacatttaaggaagaataaatagatgagacatatatatatatatatatatatatatatatatatatatatatatatatgtatagcccagtgcacaaaattcgtgcacaggagggggggtgtccctcagctcagcctgcaccctctccaatctgggacccccgtggaaaatgtctgactgccggtttaggcctaatccatgtgaatccctcaactggcagtcggacatccctctcacaatccagaactgctagctcctaaccattcacctgcctgcctgtctgattgcccctaacagcttctgcccgccagcctgatcacctcctaacttctcctgccagcttgattgatgcctaactgctcccctgccggccgatcacctccaactaccctcccctgctggcccaatcactcccaactgccctcccctgcaggcctggtcccccgcaacttccctcgcctgcaggcctgattgctcccagctgccctcccctgccgacctgattgcccacaacttccctcccttgctagccatcctatggcggccatcttgtggtggccatcttgtgacaatgtcgtggaaaccatcttgtgatgatatcatgtgaggacgttgcacaatgccacctaggcttttattatataggatatatatacacatatatgaaaagccagacaaagtgcctgaaaaaaaattagctttcaatagatttcttctccagttcctacttgtcctcttgatcttatgcttaacacctctctggacttgtcagtacaagccctgaatgtttcaccaaacaaactacacagttggtgccactggactgtttaacatactctaattaccacaatttcattggaagggaccagggttggcaaccaaactaacgctaagggatttgatcatagaattcccacaccaagtaggaggcataaaatgacttggagacctggtcatagttcctgcagctatggttcatggcatgttagcagaacacagtctctgaaagttatttcctaaaagctttccaaagagctggttcttcacacaaagagttttggccaattgcaatgacagggagaaactagaatctcaggaagaacagaatttaagccttaagtttcattttgttttggatcagatgttccttaaacataaaccttcacctgattctaattttcatttcaattccatagaccggtgatttaagtgttggggatagagttcttagaatctgtttatgtggtataaatacagaaatggtctggttctgagtattaggatcatgatccatgttttgatatgaactcaaaacctctctaagtcacgtgttattttttgaagcacagaaactggacctcacaaacaacctaaatccctctatttgagaaactaaataagcactgagttttcttgccataaagttctgcactgcatctagaatgctttcttaaatagtatctctgggagtcctggtcaagggcgtgtacctgggtttgaggtttgctccctgcccccactctgatgcatgctggaagcaaccagttggtgtgtctctctcacattgatgtttctctctctgtgttcagtttcgagttggcaagattgcaaaagtttgcttgtaaatctgtgggttttgttttttgtggggttttttttagacctttgaatgcatttttttttattgttattgatttcagagaggaagggagagggagagataatagaaacatcaatgatgagagagaatcatggattggctgcctcctgcacgccccctactgggaatcgagcccgcaacccgggcatgtgcccctgaccagtatcaaaccaaggacccttcagtctgcaggccagcactctatccaaggagccaaaccggctagggcttttgaatgcttttttttttgttccagagggtggtcagcgtctcaggccaaaccccaaagcctatttcacctattcttgtattagaaccaggggtagattcagttttgtgtttaggaatggtgagggaagtacatgtttgcttctatctttaatcatttcaaatttgttagaaacatagacttattgacctagataccactcagtgagtttagacagttgctcaagttaataatgttttaagttattgatggtggttatggacatggactttgaagtcaaagagatcgtgctttatctatttacttgatcaattatttttgatatctcaccattaagccagacactgatagaaacaaacttaaaatgatgaatataacttctcttactagcttacactacattgcttttctcagtggtattcagttttcttagttataggagtgggagaaagataaaaacatatgcctcatagagttgtcattaatattaaatgaaattaatgcatggaaaatatcctatataataaaaggctaatatgcaaatagactgaacggtgaaacaaccaaacaaccaaacaaccagtcactgtgacatgcgttgaccaccagggggcatgtgcggaacatggtgggtgttggccatggcaggatggtggagcaggtgagcaggggacaccagaccaaggtgaggtgccattcattgtcatcgggtcgagcctcctgtggttactgaaaattctttgctcctgtgcactgtggtcttgcccagagcttccacctgctaccagcaccggagccactgctcgcaccttcagccagtgctggaaccaccgctcacacctactgccgatgctctgcaccagtcccaatcactcaatgccatcaatgggtgcaagtgtggctgctggccctgatcgcccctgagggcttctccacctccccctgctcctgggggtgatcgggacagcagctgctgctcacacccgctgctagcatcagcctcaatcactccgtaccttcagcaggtatgaaaggggctggcatcgtcagtgcatgggagcggcaggagtgggactgctggcagacaggggaccagggggctgtggtgggaggggctgggcgggggcatggaggatgggctgagacctgcccctatgcctatggtagcctcacggcccacagtttctttcaaggtgcatgaatttgtgcactgagcccctagtttagcataaaatgaaatttcagtgtatggtagccactgacatgactgtggaccaacccctcattttaaagctacaaaagagctcacctggcatggctccgtggttgaacgctgatctagaaggtcacggttcaattcccagtagggggtgtgcaggaggcagccagtcaatgattctctctcatcattgatgtttctacttctctctctctcccttcctctctgaaatcaataaaaaatatatttttaaaaaggagaaagtttatgtgcatagaatgattttatgactatttgagtacttacaaaaatagtaaaaaaaaaaaaaaaaaaaaaaatttaagctgaaaaagactatatcacccaggcattctcactctccttttagcccatgccctcttacctgagcattggctttattgtcagacagccttgccaagatgggcactagtagctttatgcttcctgtccccattgctaacgattctgagaagactgtccacttggaccacatggaataaatttatgagaaagaggatttctatttctagaaagcagggattaaggaaagacatgataggcagacaagaaccatagctaccatacttcactttagtgagacaaagcgggggagggaaaggatataaatgtgtggtggtggtgggtgcacccctatgtgtttgtgaattgggattatgatgtaaagtcagaagcagtggactagaaactaaaaaatgatcagttcctaaagaaaaaaaaataaatcgtgagcactccccaaagtaaattattatatttgcattcaaatttcttcattttacactatttttcacaatcccttcttatttaacacagagaagtcaactacaacttcaggtattccagcaatcaaagtaggactcctttgggagcaattttcaaaaacctattaccatttaatgaagataaaggtgacagttaaagatactcctgggcctacctggagatcacaaggccaaaccagatttctatacaggtaaatgcaagtgatcaacctgcttttatttgcaaagaccaaatcagagaagtacactcatgaaatttagatccagagatataaaagaacaagaaaggtagatcttcccgaataagtttccaaaacaatctgaagctaatattaccaaaatacatgggctatagcggaaccttctaaacaagtatgagtaaccactgagtgaatttatgttatgagcggggtaatcaagaattttcactgctttaccatgatgatcataatttgcatcggcccaacatccttgatctcagaatgatgctcaaaaatgttaagaaacagtcttgccgctatggatcagtggttgagtgtcaacccatgaaccaagaggtctcagattcaggttcaattcccagtcagggcacatgccaggttgctggctcagtccccagtcgggggcatgcaggaaactgccaattgatgatgtttctctctcactgatgttttatcctgttctaagtttacaacctcacttctcaaaaaaaaaaaaaaaaaaatgatgctactcagtagctggtacattgtgttcattcctctgaaagtgggctttaaaattgcagctattacctctgcccagctcactggtgtctgcagggccagaggcttgcctcgttcccaaggccatcagacatacaccctggaagggctccagaagcacagccacagcaggagcacctgaccgatcctgcaccaaaaggtgcatgatttgaccagatttttggaggagcatcctggtggggcagaagtcctgagggaacaagctggaggtgatgctaccaaaaactttgaggacacctggaactctacagatgccagagaactgtccaaaacatatagcactcgggagcttcgtccagaggacagacaaatgataaccaagtcttcagaagctcctagtattactgttgattctaataccagctggtggaccaattgggtgatcccagccatctcagcacttgtcatagccctaatgaacttccatattatggggaaaaatatgaacatatgcaaaatgaccccaagaaatcacagatcaatggtaaaatatttaaaacaaagtgtgaagaacctttgctttgcaggaatatttttttaaaaatacacaaaacatgcacataacataaactacctctgctaagacaaattaacaaagcaagataaaaaaatttaaatcttatagccaccaagatagtaaagaattatgaaataaagaagaaaaggaaaccagagaggtgaatcctaggattggttgatggttctggaagagatgattaaggtagaatttcggcagaggttttggcagacccatagggctagatggttaaagtaaataaaatatatatatatatataatttagcaaaaaatcaaagttctgttcaggagtgaactctggtaaatgttcacacatctagatgggatcccaaagggctccacctatgagtaaggagtaactggaaaaaacttggagctgcactaaagtttaaatcatcacagtccctataattaaaaaatattttttattgatttcagagaggaaatgagagagatagaaacattagtgacgagactcattgatcagctgcctcctgcaagacccctactgcaggtgaagcctgaaacccggtgtgtgccctaatggggaatttaactgtgacctcctggctcataggttgatgctcaaccactgagccatgccagctgggttcagtccctataattagatgaaggtgatctaagatagctggggccctaaccacctaccaaaggcaaatgtaagtcttccttggaggaagataattctaggcctcaaataattgtaacttttttttcacatacaatgtcctgtacctcaagaaaatgaaccagaaatggaaagggatatctcactacagagcctacagatattttttttattattataaaatcatttaatgatatttattattaacaactttatgttcctaaatttgaaaatttaggtaaattggacaagtttctaggaaaaaatacaacttactcaaactgacataaaaaggaatagaaaatctgaatcattcaataaacttaaaggaattaaaaactataagaaaggcacattgtgatttcccccagttactccattctcttttcctattaaaactttattttcctatttttgagaatgattcatttatttggctattggaatatttcagaaattaaaaaaatatgattacctacttcctctccatctttactcaatgttatcaaacttttctatgttttacttttgagaattacatatgttgtcttctcaattatatttcaaaccataggagggcagaccttgtggcatagaagaaaaaagagaggaagaagggaggtctgggaaaaagaaggaaaaagagaattgtgtttctggtttttacagttctacatccatatcctctgcagcctatttatagcctgctgatgtagaataatcataaaaggtccacgacagtatcctatataataaaaggctaatatgcaaatagaccaaatggtagaacagccgaacaaccgaacaaccagttgctatgacatgcactgaccaccagggggcgtgcgcaaaacatggcaggcatccgctgcgacgggatggtggagcaggtgagtgggggcaccagaccaaggcagggcgccggttgctgtcattgaggcaagcatctggtggttactgaaaattctttgctcttgtgtaccatggtcctgcctggtgcttgcacctgctgccagtgctggccccattcgcacccactgctggcacccggcgccagccccaatctctcggtgccatcagtgggtgcaagctgtggctgccaatcccaattgcccctgagggcttctccacctccacctgctcctgaggggcgaacagggcagcaactgctgcttgcacccgatgatggcgccagccctactcatacccattgctggtgccaggacccaattgctccgcaccttcagtgggtacgagcagggccagcactgtcagtgcgctggagtggtggtggcgggagcagggctgctggcagacaagggactgggggtgcagtgggagagtccaggcaagggcatggaggatgagccgagacccgcccctgtgcccactgcagcctcgcagcccacagttcctttcaaggtgcatgaattcgtgcactgggcccctagttaaacatactaagcaatatgacagagtggttaatatatatgagtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtcttggagttcttgagcaagatcattggataaattacttaacctctctaattttcaatttcttcagtgataagaaaagatcattgtggaacacctcacacatcagtgtaaatattaagtgataaaatacatggaacagtattattgcaaagcctaggaatagaagaattccataaatagtagctatcatttatgattatgatggagagaaaataaattaaaatatatataatttaattttagtttttaatatgtgttgacgaggaaaatcccaaatcagtcaaagtctctttatgcttacttcactcaaataaaaattttctatactaagtagcagaaataagtaaggcccttctggtaaactgttctaagttgggctacctggtggtctagagttgatcagaattttgaaacctatcctaagaaaagatggtgaaataataacttggcttttcttaaaaaatatagaccatctttgccagttgtccgttgtatggtttgttgagattgtctggaatttccaagccttaggactattaaaataattcaccagaactgcctccaaataattttatattgtgtaatttcttaacgctgcactttgggatggatgttctcggaatcacatgggtttagtgagagtcctgcataaagcatgttttgtcataaaaagcatccccttgatgcttagggaatgatgttgataggttttagacacagaacctacaatgactggaacatatttatcgaagcacagagggcagttgaaaatggctgttacaaactctgaacaaacaaggtcattgacatttcagtgtgaaccctcaagagcccttgaggggaaaaaaaattactttcacagggaactttcaattgccttctctctgattctaaaaagaatagaccagttcttttcttttcccccttgtccgttggctaagagaccacaaaaaaaaaataatcctttgcttaatttccagaaaaaaattatattgcttttagtaaaattcatttttcacctcaaggatagagctccaaattattatgtcaaactacattttctttaaacatctgaaatttttcctccgatttttattaaaagtacattctactttatgatatagaatgaaatattttggtattgatgcttttcttaataaattaagacatataggaacgctatttaaaatcatgctggggtggtaaccggtttggctcagtggatagagtgttggtctgtggactgaagggtcccaggttcaattccagtcaagggcatgtaccttagttgtgagcacatccccaatagggggtgtgcaggatgtttctaactctctattcctcttccttccactctgtaaaaaatcaataaaatatatttttttaaaaatcaaatcatgctggaagataaatgtcagttatttttaccaatgatttttcattttgggaaagcagcttttctatatccaatggcattactcaaaaagaagacttttctaaaaagtgttcctaatcattttgaatggagctaactgaaaatggcctcctggactctgctgcatttgggtttcccatccacttagattaattcattcaaaatggaaaatcatgtgactgaatttttactgttttttgttgttgttcttaatcctcacccaaggatatttttcccattgatttttagagatagtgggagggagggggtgagagagggacagagagaaatatggatgtaagagagacaatcgattggttccctccgtgtgcaccctaaccgaggccaaggatcaaatttgcaacccaagggaacttacgtgccctagattggaaatataacctgagactcttcaacatgcaggccgacgctccaaccattgagcaacatcagccagggctctatttttttttttttttaaattgagaggattgttgggagatgtgtaagtgatttattatctgttttgttttaatgaagacatttggaaataggaaatgatttcagtagacaacttgcatggaacaaatttgcctttttggtcacggcatttttcctttgtgacatgaaggcttttcagggggtacacaccatgtttccttttcatataaaacgagaggcccagtgcacaaaattcgtgcattcagggaggggttcccaaagcttggcctgcaccctcttgcagtctgggagccctcaggggatgtcagacatgccaatgagcccggcttctgtggctgagcggcgcttcccctgtgggagcacactgaccaccaggggggcagctcctgcattgagtgtctgccccctggtggccagtgcgcatcatagcaaccagtccttctaccattgggttgattgcatattagccttttattatataggactagaggcctggtgcacaaaatttgtgcactctgggggggtgggttccctcagcctagcctgtccctctccacagtccgggagccctcagtggatgtcctactgatggccacagtctggagcagaggctcggagcaggcgtcccgtgtgtgtgtgtgtgtgtgtgtgttgtgtgtgtgttgtgtgtctgctaggggcctgcccccatccgcgcactgccggggtccaaccccagcggggtcaggggtcccccaaaggtgtggacggagtcagcgaagaaggaatgacacgggagagggcattcagatgatcatcatagcaaggttctctagccacgttctggtcaggatctccagtgaggtcttggttaggatctccagccaggttctgtgtccatggttctcttgctaggttctccaggttctcagccagttctgtagccatgttccctcactaggttcttccagccaggttgtgtccaggttctccagccaggctcagtcaccaggttctagtcaggttctcttgccaatttctgtagtcaggttcagtccaggatctttgccatgttctctccagcaaagttcttctgtctctaggttctgtgtaggttctgtgttcttagttctgtcttcttagttctgtgtctaagttctgtctcctaagttctgtgttgttacatctgtatttataccagttgattccaatcctatcaatctctattccaaaggttagggcgtttctatctccattccagggagtaaagattatgtagcttagcatgttgttcattagtatagttgattaattacccgcctggcacttagttaagaggttttattccctccctaacttcagtgaaaatcccctacctggggaaacaacctttctcagagaccttggttaaaacacatagtccaagaaggtgacaaacatattaaggaacagtatgccatatattccaggtcccttgaaacagcaagcatggaccggctcccggcacgccacgtaggctcaaagcaatattgccgtgtgttgtgtgtgtgtgttgtgtgtgtgtgtgtgtgtgtcccccactgggggtctgtccccagccacattggaggctcagagcaggcactgtgtgtgtgtgtgtgtgtgtgtgtgttgtgtgtgtgttctgctgggggcctgcccccagcacaccatagaggcttggagcaggagcccctctgtgttgatctcttaggctccgggcccccaccttgcatgttgtacctcccttgagctatggccaggatgggggtgcttgcaagccagactttcctgctcctggatcgccatggcgtccggggagcaggcccagctgggagctgcccacaggctgggctttcctgctcccaggtctccgggcaaccagggagcagcccagcttggagctgcccgcaggccgggcttgttctcctgttcccggatccccatggtgactgggaagcaggcctatcgggggctgcccgccaggccaggctttcctgctcccgtttcaggcccagcttggggctgcccacaggccaggctgctttccctcccagatcaggtgcagctgggggctgcctacaggctgcacttttctacttacgatggccagattcaacacagtgacccagggcccagcggcactttcctgctctccaatagtcataggtcccgggcgtggggcggggcttagcagtgcaggcttaggcagcacgggatccccagctggggatcaggtggcacacagggatccggcgggggcaggacttatgcccactacccagggctgcacatgggcccggatggcagctcgagctgtcccgccctgccagtagtataggatagaggcctggtgcacgggtgggggctggctggtttgccctgaaagtgtcccagatcagggtggggggtcccgcttgggtggcctggccagctgcatgaggggatgatggctgtttgcatctggtcacacccccttcaggtgggggtccccactggggtgcctggccagtctgggtgagggcctgagagccattttcaggctggagggcaacttaagctcccacgctctccttttttttcttttttttttttattctgggattttatttaccttgtatagctgtcactggagctgagagcaggctccagctctgggggctgaaagcaggttctgggcttttttggcttctataattgaaactctgttgcccatcactgcagctctaagctctgagggcctaagctggctgaaagcaagtctcgtggaggcttgtttagcttctataattgcaacatagttgcttagagtgtagctcagaggctggtcatggcaggccgggaatgttggcttcctccatcactggagcaagcaagcctcctgttcgcttcagctgcgtggctgctggccgccatcttggttggcagttaatttgcatatcctgctgatatcctgcatatcctgggaagggtgttggggttatggtcaatttgcatgtttctcttttattagataggattatatcaccatgaaactttagacctgtcagtgacagccaaagtatatagggaacaggaagctctgactctgacttgatgagtgaaacaaacaggaaagctctggtctctttctttagatgtgagtatccaggtgtgcgctggatcatacaagcctgtccacactcttgcaactggagacagggagaggaaagctgacctgccagtgacatagacagcacatgcgtccctccaccctgccactctgctctgggcatcagatcagcatctagttgtgtccctacagtcctatagaggaactcaaatattggaagcaaagtaggtatagactctcactgagcacattcaaaaactcttatttttaattttcaattcagaaaaataaatacattttatgggcaattctttctacattactccttattgtttaattggataccagacacagtaataatactaagaaaggttgtccctaaagttgaagtaacaatgaaatgttgtaagtcatttcatggcgagggagaggagggaatgagttttgagccattttaaggtgactcattgtgtcctctagaggaaatgacttcagaaaatatctgttagaaaacactactcttggagagtctcataatgtggagagctgagtgtttgtagcttcataaaaagtgattgttaactgttatcagaaagcatcttgatgttatttttaaaattatttgcagcagagacaatatgtattagcttgtgcctccttccttccccgcaagaaagactctttgttagccattgtcttgtatgtttcatggtacaacatcagggaaatgggtagataagataacatctgaaaaggaggcagaaataaccagtaagcctagatcactttataataaagaaagaatgaatttgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtacatatatatattgagagaaagagagagagagagagagagagagagagagagagagagagagagagagagattagttgtttggttttattataaggaactagctcacacaataatggaggctgagaagtttaagatctacagtcaaaaggtagggacacaaggagagctgagggtatagttccagcctgagcctaaaggcctgagaacgaagagagccactagtgtaaattctagttcaagtctaagtccaaaggcgggaaaagatggatgtcccagcttaaagataggcagtcagagagaaagaatgctttcttggcctgccttttattctatccagtcctccagtggattcgatggggcttgcctacattaaagaaggttatttgctttactcagtctatcgattcaaatgctgaactcatgcagagaataaacaccattacagacacatccagaattatgcttaaccaaatttctgagcactcatggcccagtcaagttgacaaatgaaattaaccttcatacttttgtatacaaatatacttttgcaccggaacagaacagagaacccagaaattgacccaagccattatgctcaattaatatttgacaaaaggaagcaagagcatacaatggagtaagacagtctcttcaataatggtactgggaaatttggtcaaatccataccaaaaaaaaaaatgaaactagatcaccaacctataccatacacaaaaacaaactcaaaatggttaaaggacttaaatttaagatgggaaaccaaaaaatcttaaaagactacataggcagcaaaaaatatcaggcatatgtcacagcaatatctttactgatacagcttctagggacgggtgggcaaacgttttgaatcgagggccacaatgggttcttaaactggacgtgaaggccggaacaaaagcatggatggagtgtttgtgtgaactaatataaattcaaagtaaacatcatacataaaagggtatggtctgtttttgggttttttttagttttatcatttcaaatgggccggatcgcccggcgggccgtagtttgcccaggctgttctaagggcaatgggaaaataagggagaacatagaccaaataggactacatcaaaataaaaagctttctgcacagcaaaagaaatcatcaacaaaacaacaagaagcccactggcatgggagaacatttgccaatgctattcagataagggtttatcctccaaaattttacagagaacttcatacaactacaaaaggaagataaacactcaatcaaaaaatgggcaatggacctaaattagacacttttaaaagaggacatagaaaggccaagagacaatgaaacatgctcaaagtcactaatcctctgaggagatgcaaatcaaaaacaacaatgagataccactcacactgtaagaatggctatcatcaacaaataacaaaggacaagtactggcaaggatggtggagaaaaaggaaccttcatgcactgctggtgggaatgcagactgggtgcagccaactgttggaagacgtatggaatttcctccaaaaattaaaaatggaattcccatttgacccagttaattcccacttctagagactaccatcccaagaactagaacaccaatcagaaaggatatatgccaccctatgttcatagcagcacaatttacaatagctagatttgaaaacagcctaagtgccccatcagcagatgagtggattaaaaaactgtggtacacctaacacaatggaatactattcctaatatataaaaagccaggggctgtaaaccaaaacaaacgggctggatgaccgaacagcaggctgcatggggcgacaaggccggcagggggattagtgagggacaatcaacgactgaacagcagctgcataaggcgaccaggcaggcagggattttagtgaaggatgaccaaatgactgaacagcagctgcgtggggtgacccgcctgcaggggggttagagagggacgaccaaaccactgaacagcaggctgcatggggcaaccaggccagcggggggcagttgggggcaacacggccagcaggagggcagttgggggtgattcaggcagcaggcaggtgagcagttagaaggccagcagtctcagattgtgagaggggatgtccaactgccggtttaggccttattctggggatcgggcctaaaccggcagtcgacacccccaaggggccctgattggagagggtgcaagctggctgggggatggatgaaacagcatggtggaactggagagcattatgctaagggaaataagccagttggagaaagataatatatcacatgatactcactagatgtgggatataatgatcaacaaaactgtgaacaaaatagatccagggacaagtggcatgggcggtgggttagaagatcaacctaaagacttgtatgcatatatatgcatattagcattagccaatggacacagacactggggggcggtggggacttttcctcaggggttgggagtggttgggggcgggtcaatcagggaaaagggagacatatgtaatactttaaaccataaaaaataccaattaagaatgtttacatatatatacattatctatagtgctaaaatgtaaaagattgatataccaggtgtacctacttgatacatttctccaacatggatggagtctcaataacattatgtttgaacaaatgaagccagacagaaaagggtagctaatgtataatcacatttacgtgaggttaataaacaagcaagttggtgatgatagaagtcagaatagtggtgttatccttctgagaaaggtgggcttattggattgaaaggggcacaggggaactttctccagctgacggtgacaattttacatttttccaggtcttagtcacatggatatttgagtaacatataacatatggtgtaagtatcctaatattataaaaacccagggatccataacgacctaccgaaggctcaactgaacaccagcagtcctgcagtcagtgctggggctgccgtggcgaccagcactgactaccgccggtgcgggcacagccaacccagcactgactgccgagggagctgcggatcaggcccagagagagaggcacgtctgatctgcaggcctccatggcagctgttgatcaactgttgcctcttctctctctctccccgggcttcaccagcagcccacgctctttctctcaaggcctccgtgagggggctgcttatcagtcccaccttctctgatcaggcctgggagataggcccagagatgttgactggcatagaaactaaccaatcagaccaaatctgggtgacagtaggagccaatggctgcctaggaggcagagcttttgatgctgactggcatagaaacttgaccaatcagaaccaaattggccagcaggggagggcagttgggggcaaagatcaggcctgcaggggagagcagttgggggcagatcaggacagcagaggagggcagttagggattatcaggcaggtagaggcagtagggcaattcagctggcaaggagagggcagtttgggggcgagatcaggctgcaggggatggctggtagggggcgagatcaggccggcaggggagggcaggttaggggtgatcaggcaggcagtcagaggggttaggggcaatcaggcaggcagaggagttagggaccgatatgcaggaaggcaggtgagtggttaggagccagcagtctccgattgcaaggagggatgtctgactgccggtttagtccaaatccctgtgggacatccccaagggtccccaaatggagagggtgcagactgggctaagggaacccccgctccatgcggattttgtgcaccgggccactagtacttatatattgtatacatttaggatgtgtgcatttcactggatgtaaattatactttaaacttgtgaggaggggaggttcaatacaaaaattattcccaacgaatgaaataaatgaaatttatttttttttaatcaagagacctaaaccaagcattacaaataaaaatttctctcagatgggcccaagccagacagtttccataagcaaacccatgaaagccattgcctgggttcctaagactcagatgacaaaatttaaaatcctatggaattgtccccaaaaagatgactcacaggtgaccttgaaatctgagtggtgtgtatactatttcagctacagtgtcacacattccattgtagttatcacctggtatttcataattatttgtaaataattcagtagcaacattgactaaactccactctaccaccaaaattgaaaacataacattgatagattttaaaattatttcatacttacaagtattagtttatttgttaacaggtccaaacattttgactgcatggatccttctgggtctacgtcctgattttaattcctgttctgtctgctctctcacaactcaccacaaagtattttgaaaattgcataagatgctctggccagtgtgtctcagcgggttagagatcagcctgtgcactggaagggtctcaggtttgattcccagtcaagggaacatacctgggcaggcaaccaatcgatgtgtcttctctcacatcgatattttttctctccccctccaccctcccttccactctttctaaaaaatatcaatgaaaaaaatccagtgaggataaacaacaacaaattgcaccaataataaaaccaaattatcagtaataatttattcctaataaagataaaaccatggggccaaatggtacatgtcacatttaatttatcaatattgtcatttattataatacacaactttatgatgtcaaagtatctggatttctatggtcacattgggtagttcagtgttatacacctttgggcacatctcaaaagaatccagttgattcgcttgacacaggctatttaatgcaaaaggaggcccttaaaatgcaggtgacaatattaaacagataaacttgatgcagaaatcacaaatacctctcagcctgattgcccaaagaccaaaactgtgacacacaaatgggaccacacacgcttagctcagccatgccaggggagggctgaccgggcttctcaactcagaaaaattcagtcatgctccaccctcatacttcacatttctttatcaagcaggggctgcagtgggaagcaggcacattaatgacaacccacattccacatgtgcgcacacacaccacaatggagctatgatttcagaaaccttctaataaacaggatgttggagtgctgattgatttacacgttgataaaagttccgttgaaatgagcctactatctatcatgcaagtggctgtgaaaaacttaattgcacaaactgttcttgttggaacagggagctatctgctccttttccttgcatggaccagagttcctccctcacagctactgcaggacacaaaaagacagcacttaaaactggatcctgttctgcaaaccgtgcaagtcaatccggcgccgttcacagacatttttatttctgatcatgtctccaagtgtagacatctgagaaatataatcttgcctcagcaaggctcaaaggagacatttgatctctaagtcatagtctcactgatcatgcttagttagtaagtcacagaaaaaataaattaaattgcatcatgggcttaaaatattttattctcccttataatttccatatcccaaagccaacttgtttaactattaaatgtacgtactccacttaggaatgaaaactgtaaactgttttttcctgatgagcacatgtgatgtctttatcattaatagtgaataggtggtaatattccttcttggttttgtttttgggttttttatgttaacactattacacatgtcccccattccttccccccttttgcccacttccacccagtccctacccacccctctctggccttcaccaccactgttttctgggttcgtgggcttgcatatattgttctttggctaatcctctacacttcctt");
    auto exactAlignment = getPaddedAlignment(ac, begin, end, aSequence, bSequence);

    cast(void) exactAlignment[47139 .. 73889];
}

/**
    Get the designated records of `dbFile`.

    Throws: DazzlerCommandException if recordNumber is not in dbFile
*/
auto getDbRecords(in string dbFile, in DBdumpOptions[] dbdumpOptions = [])
{
    enum size_t[] allRecords = [];

    return getDbRecords(dbFile, allRecords, dbdumpOptions);
}

/// ditto
auto getDbRecords(Range)(
    in string dbFile,
    Range recordNumbers,
    in DBdumpOptions[] dbdumpOptions = [],
)
        if (isForwardRange!Range && is(ElementType!Range : size_t))
{
    return readDbDump(dbdump(dbFile, recordNumbers, cast(string[]) dbdumpOptions));
}


private struct DbDumpLineFormatTuple
{
    char indicator;
    char subIndicator;
    string format;
}

private enum DbDumpLineFormat : DbDumpLineFormatTuple
{
    totalReadNumberCount = DbDumpLineFormatTuple('+', 'R', "+ R %d"),
    totalMCount = DbDumpLineFormatTuple('+', 'M', "+ M %d"),
    totalHeaderCount = DbDumpLineFormatTuple('+', 'H', "+ H %d"),
    maxHeaderLength = DbDumpLineFormatTuple('@', 'H', "@ H %d"),
    totalSequenceCount = DbDumpLineFormatTuple('+', 'S', "+ S %d"),
    maxSequenceLength = DbDumpLineFormatTuple('@', 'S', "@ S %d"),
    totalIntrinsicQualityVector = DbDumpLineFormatTuple('+', 'I', "+ I %d"),
    maxIntrinsicQualityVector = DbDumpLineFormatTuple('@', 'I', "@ I %d"),
    readNumber = DbDumpLineFormatTuple('R', '\0', "R %d"),
    header = DbDumpLineFormatTuple('H', '\0', "H %d %s"),
    pbLocation = DbDumpLineFormatTuple('L', '\0', "L %d %d %d"),
    pbQuality = DbDumpLineFormatTuple('Q', '\0', "Q %f"),
    sequence = DbDumpLineFormatTuple('S', '\0', "S %d %s"),
    arrowSNR = DbDumpLineFormatTuple('N', '\0', "S %d %d %d %d"),
    arrowPulseWidth = DbDumpLineFormatTuple('A', '\0', "A %d %s"),
    intrinsicQualityVector = DbDumpLineFormatTuple('I', '\0', "I %d %s"),
    quivaDeletionValues = DbDumpLineFormatTuple('d', '\0', "d %d %s"),
    quivaDeletionString = DbDumpLineFormatTuple('c', '\0', "c %d %s"),
    quivaInsertionValues = DbDumpLineFormatTuple('i', '\0', "i %d %s"),
    quivaMergeValues = DbDumpLineFormatTuple('m', '\0', "m %d %s"),
    quivaSubstitutionValues = DbDumpLineFormatTuple('s', '\0', "s %d %s"),
    repeatProfileVector = DbDumpLineFormatTuple('P', '\0', "P %d %s"),
    maskTrack = DbDumpLineFormatTuple('T', '\0', "T%d %d %(%d %d%)"),
}

struct DbRecord
{
    static struct PacBioReadInfo
    {
        id_t well;
        coord_t pulseStart;
        alias begin = pulseStart;
        coord_t pulseEnd;
        alias end = pulseEnd;
        float readQuality;


        @property coord_t pulseLength() const pure nothrow @safe
        {
            return pulseEnd - pulseStart;
        }

        alias length = pulseLength;
    }

    id_t readNumber;
    alias contigId = readNumber;
    string header;
    PacBioReadInfo pacBioReadInfo;
    alias location = pacBioReadInfo;
    string sequence;
    byte[] intrinsicQualityVector;
    alias intrinsicQVs = intrinsicQualityVector;

    enum maxQV = 50;

    static byte fromQVChar(const char qv)
    {
        if ('a' <= qv && qv <= 'z')
            return cast(byte) (qv - 'a');
        else if ('A' <= qv && qv <= 'Y')
            return cast(byte) (26 + qv - 'A');
        else
            return -1;
    }


    static char toQVChar(const byte qv)
    {
        if (qv <= 25)
            return cast(char) ('a' + qv);
        else if (qv <= 50)
            return cast(char) ('A' + qv - 26);
        else
            return '\xFF';
    }
}

private struct DbDumpReader(S) if (isInputRange!S && isSomeString!(ElementType!S))
{
    static alias dstring = immutable(dchar)[];
    static alias DbDump = ReturnType!getDumpLines;

private:
    DbDump dbDump;
    bool _empty;
    id_t numReads;
    DbRecord currentRecord;
    dstring currentDumpLine;
    size_t currentDumpLineNumber;
    dchar currentLineType;
    dchar currentLineSubType;
    debug dchar[] allowedLineTypes;

public:
    coord_t totalSequence;

    this(S dbDump)
    {
        this.dbDump = getDumpLines(dbDump);
        debug with (DbDumpLineFormat)
        {
            this.allowedLineTypes = [
                totalReadNumberCount.indicator,
                totalMCount.indicator,
                totalHeaderCount.indicator,
                maxHeaderLength.indicator,
                totalSequenceCount.indicator,
                maxSequenceLength.indicator,
                totalIntrinsicQualityVector.indicator,
                maxIntrinsicQualityVector.indicator,
            ];
        }
        this.popFront();
    }

    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty DbDumpReader");

        if (dbDump.empty)
        {
            return setEmpty();
        }

        readNextDbRecord();
    }

    @property bool empty() const pure nothrow
    {
        return _empty;
    }

    @property DbRecord front() pure nothrow
    {
        assert(!empty, "Attempting to fetch the front of an empty DbDumpReader");

        return currentRecord;
    }

    @property size_t length() pure nothrow
    {
        return numReads;
    }

    static if (__traits(hasMember, dbDump, "destroy"))
        alias closePipe = setEmpty;

    void setEmpty() pure nothrow
    {
        static if (__traits(hasMember, dbDump, "destroy"))
            dbDump.destroy();
        _empty = true;
    }

private:

    static auto getDumpLines(S dbDump)
    {
        return dbDump.enumerate(1).filter!"!a[1].empty";
    }

    void readNextDbRecord()
    {
        currentRecord = DbRecord.init;
        peekDumpLine();

        while (true)
        {
            debug _enforce(allowedLineTypes.canFind(currentLineType), format!"forbidden line type `%c` (allowed: `%(%c%)`)"(currentLineType, allowedLineTypes));

            with (DbDumpLineFormat)
            {
                switch (currentLineType)
                {
                case totalReadNumberCount.indicator:
                    static assert(totalMCount.indicator == totalReadNumberCount.indicator);
                    static assert(totalHeaderCount.indicator == totalReadNumberCount.indicator);
                    static assert(totalSequenceCount.indicator == totalReadNumberCount.indicator);
                    static assert(totalIntrinsicQualityVector.indicator == totalReadNumberCount.indicator);

                    switch (currentLineSubType)
                    {
                    case totalReadNumberCount.subIndicator:
                        read!totalReadNumberCount();

                        debug allowedLineTypes ~= [
                            readNumber.indicator,
                            header.indicator,
                            pbLocation.indicator,
                            pbQuality.indicator,
                            sequence.indicator,
                            intrinsicQualityVector.indicator,
                        ];

                        // do not try to read empty dump
                        if (numReads == 0)
                            return setEmpty();
                        break;
                    case totalMCount.subIndicator:
                        assert(readTotal!(totalMCount, size_t) == 0, "unexpected totalMCount != 0");
                        break;
                    case totalHeaderCount.subIndicator:
                        break; // ignore
                    case totalSequenceCount.subIndicator:
                        read!totalSequenceCount();
                        break;
                    case totalIntrinsicQualityVector.subIndicator:
                        break; // ignore
                    default:
                        error(format!"unknown line sub-type `%c`"(currentLineSubType));
                    }
                    break;
                case maxHeaderLength.indicator:
                    static assert(maxSequenceLength.indicator == maxHeaderLength.indicator);
                    static assert(maxIntrinsicQualityVector.indicator == maxHeaderLength.indicator);
                    break; // ignore both
                static foreach (lineTypeFormat; [readNumber, header, pbLocation, pbQuality, sequence, intrinsicQualityVector])
                {
                    case lineTypeFormat.indicator:
                        static if (lineTypeFormat == pbLocation)
                        {
                            auto pbReadInfo = currentRecord.pacBioReadInfo;
                            auto wasLineTypeDone = pbReadInfo.well != pbReadInfo.well.init ||
                                                   pbReadInfo.pulseStart != pbReadInfo.pulseStart.init ||
                                                   pbReadInfo.pulseEnd != pbReadInfo.pulseEnd.init;
                        }
                        else static if (lineTypeFormat == pbQuality)
                        {
                            auto pbReadInfo = currentRecord.pacBioReadInfo;
                            auto wasLineTypeDone = !isNaN(pbReadInfo.readQuality);
                        }
                        else
                        {
                            auto wasLineTypeDone = mixin(
                                "currentRecord." ~ lineTypeFormat.to!string ~ " != " ~
                                "currentRecord." ~ lineTypeFormat.to!string ~ ".init"
                            );
                        }

                        if (wasLineTypeDone)
                        {
                            debug allowedLineTypes = [
                                readNumber.indicator,
                                header.indicator,
                                pbLocation.indicator,
                                pbQuality.indicator,
                                sequence.indicator,
                                intrinsicQualityVector.indicator,
                            ];
                            return; // DB record completed; stay on current dump line
                        }

                        read!lineTypeFormat();
                        goto break_;
                }
                default:
                    error(format!"unknown line type `%c`"(currentLineType));
                }

                break_: // pseudo-break for use with static foreach
            }

            if (popDumpLine() == Yes.empty)
                return; // EOF reached
        }
    }

    void peekDumpLine()
    {
        auto currentLine = dbDump.front;

        currentDumpLineNumber = currentLine[0];
        currentDumpLine = currentLine[1].array;
        currentLineType = currentDumpLine[0];
        currentLineSubType = currentDumpLine.length >= 3 ? currentDumpLine[2] : '\0';
    }

    Flag!"empty" popDumpLine()
    {
        if (dbDump.empty)
        {
            return Yes.empty;
        }

        dbDump.popFront();
        ++currentDumpLineNumber;

        if (dbDump.empty)
        {
            return Yes.empty;
        }
        else
        {
            peekDumpLine();

            return No.empty;
        }
    }

    Int readTotal(DbDumpLineFormat lineTypeFormat, Int)()
    {
        Int total;
        currentDumpLine[].formattedRead!(lineTypeFormat.format)(total);

        return total;
    }

    void read(DbDumpLineFormat lineTypeFormat)()
        if (lineTypeFormat == DbDumpLineFormat.totalReadNumberCount)
    {
        currentDumpLine[].formattedRead!(lineTypeFormat.format)(numReads);
    }

    void read(DbDumpLineFormat lineTypeFormat)()
        if (lineTypeFormat == DbDumpLineFormat.totalSequenceCount)
    {
        currentDumpLine[].formattedRead!(lineTypeFormat.format)(totalSequence);
    }

    void read(DbDumpLineFormat lineTypeFormat)()
        if (lineTypeFormat == DbDumpLineFormat.header)
    {
        coord_t headerLength;

        currentDumpLine[].formattedRead!(lineTypeFormat.format)(
            headerLength,
            currentRecord.header,
        );

        assert(currentRecord.header.length == headerLength, "mismatched header length");
    }

    void read(DbDumpLineFormat lineTypeFormat)()
        if (lineTypeFormat == DbDumpLineFormat.readNumber)
    {
        currentDumpLine[].formattedRead!(lineTypeFormat.format)(currentRecord.readNumber);
    }

    void read(DbDumpLineFormat lineTypeFormat)()
        if (lineTypeFormat == DbDumpLineFormat.pbLocation)
    {
        currentDumpLine[].formattedRead!(lineTypeFormat.format)(
            currentRecord.pacBioReadInfo.well,
            currentRecord.pacBioReadInfo.pulseStart,
            currentRecord.pacBioReadInfo.pulseEnd,
        );
    }

    void read(DbDumpLineFormat lineTypeFormat)()
        if (lineTypeFormat == DbDumpLineFormat.pbQuality)
    {
        currentDumpLine[].formattedRead!(lineTypeFormat.format)(
            currentRecord.pacBioReadInfo.readQuality,
        );
    }

    void read(DbDumpLineFormat lineTypeFormat)()
        if (lineTypeFormat == DbDumpLineFormat.sequence)
    {
        coord_t sequenceLength;

        currentDumpLine[].formattedRead!(lineTypeFormat.format)(
            sequenceLength,
            currentRecord.sequence,
        );

        assert(currentRecord.sequence.length == sequenceLength, "mismatched sequence length");
    }

    void read(DbDumpLineFormat lineTypeFormat)()
        if (lineTypeFormat == DbDumpLineFormat.intrinsicQualityVector)
    {
        coord_t vectorLength;
        char[] qvChars;

        currentDumpLine[].formattedRead!(lineTypeFormat.format)(
            vectorLength,
            qvChars,
        );

        assert(qvChars.length == vectorLength, "mismatched QVs length");

        currentRecord.intrinsicQualityVector = cast(byte[]) qvChars;
        foreach (ref qv; currentRecord.intrinsicQualityVector)
            qv = DbRecord.fromQVChar(qv);
    }

    debug void disallowCurrentLineType() {
        allowedLineTypes = allowedLineTypes
            .filter!(type => type != currentLineType)
            .array;
    }

    void error(in string reason)
    {
        _enforce(false, reason);
    }

    void _enforce(bool condition, lazy string reason)
    {
        enforce!DazzlerCommandException(condition, format!"ill-formatted DBdump output: %s (line %d)"(reason, currentDumpLineNumber));
    }
}

private auto readDbDump(S)(S lasDump)
{
    return DbDumpReader!S(lasDump);
}

unittest
{
    enum testDbDump = q"EOF
        + R 5
        + M 0
        + H 15
        @ H 3
        + S 281
        @ S 63
        R 1
        H 3 Sim
        L 1 0 62
        Q 0.851
        S 62 ctaaattaacacttgtgatgaaccagtgaggaaggaggctggctaaacaatgtgaacggttc
        I 1 q
        R 2
        H 3 Sim
        L 2 0 63
        Q 0.852
        S 63 cctaactaaaccttctgaaactacagcgcaagatcagagggggtttgaaggtcatattattat
        I 1 l
        R 3
        H 3 Sim
        L 3 0 62
        Q 0.853
        S 62 aaccgatgagaaatccatatatctgggagctagagacaccaagaaaaagataccagccaaaa
        I 1 m
        R 4
        H 3 Sim
        L 4 0 62
        Q 0.854
        S 62 ttttgttcatcaaatgcaggccataaatccaatttagccactggctttcacgtaaccgttca
        I 1 S
        R 5
        H 3 Sim
        L 5 0 32
        Q 0.855
        S 32 gtgtctgctgttttttttcttttagtggacat
        I 1 j
EOF".outdent;

    import std.algorithm : equal;
    import std.string : lineSplitter;

    auto dbDump = readDbDump(testDbDump.lineSplitter);
    auto expectedResult = [
        DbRecord(
            1,
            "Sim",
            DbRecord.PacBioReadInfo(
                1,
                0,
                62,
                0.851,
            ),
            "ctaaattaacacttgtgatgaaccagtgaggaaggaggctggctaaacaatgtgaacggttc",
            [16]
        ),
        DbRecord(
            2,
            "Sim",
            DbRecord.PacBioReadInfo(
                2,
                0,
                63,
                0.852,
            ),
            "cctaactaaaccttctgaaactacagcgcaagatcagagggggtttgaaggtcatattattat",
            [11]
        ),
        DbRecord(
            3,
            "Sim",
            DbRecord.PacBioReadInfo(
                3,
                0,
                62,
                0.853,
            ),
            "aaccgatgagaaatccatatatctgggagctagagacaccaagaaaaagataccagccaaaa",
            [12]
        ),
        DbRecord(
            4,
            "Sim",
            DbRecord.PacBioReadInfo(
                4,
                0,
                62,
                0.854,
            ),
            "ttttgttcatcaaatgcaggccataaatccaatttagccactggctttcacgtaaccgttca",
            [44]
        ),
        DbRecord(
            5,
            "Sim",
            DbRecord.PacBioReadInfo(
                5,
                0,
                32,
                0.855,
            ),
            "gtgtctgctgttttttttcttttagtggacat",
            [9]
        ),
    ];

    assert(dbDump.length == expectedResult.length);
    assert(dbDump.equal(expectedResult));
}

/**
    Get the FASTA sequences of the designated records.

    Throws: DazzlerCommandException if recordNumber is not in dbFile
*/
auto getFastaSequences(in string dbFile)
{
    enum size_t[] allRecords = [];

    return getFastaSequences(dbFile, allRecords);
}

/// ditto
auto getFastaSequences(Range)(in string dbFile, Range recordNumbers)
        if (isForwardRange!Range && is(ElementType!Range : size_t))
{
    string[] dbdumpOptions = [DBdumpOptions.sequenceString];
    auto numRecords = recordNumbers.save.walkLength;

    if (numRecords == 0)
        numRecords = getNumContigs(dbFile);

    auto sequences = readSequences(dbdump(dbFile, recordNumbers, dbdumpOptions));
    size_t numFoundSequences;

    string countedSequences()
    {
        enforce!DazzlerCommandException(
            !sequences.empty || numFoundSequences >= numRecords,
            "cannot read sequence: dump too short"
        );

        if (sequences.empty)
        {
            assert(numFoundSequences == numRecords, "unexpected excessive sequence in dump");
            return null;
        }

        ++numFoundSequences;
        auto currentSequence = sequences.front;
        sequences.popFront();

        return currentSequence;
    }

    return generate!countedSequences.takeExactly(numRecords);
}

/**
    Get the FASTA sequence of the designated record with prefetching to reduce `fork`s.

    Throws: DazzlerCommandException if recordNumber is not in dbFile
*/
string getFastaSequence(in string dbFile, id_t recordNumber, in id_t cacheSize = 1024)
{
    // FIXME the cache size should limit the number of `char`s retrieved, ie. control the memory
    // requirements of this function
    static uint _dbIdx;
    static id_t[2] _firstRecord;
    static string[2] _dbFile;
    static id_t[2] _numRecords;
    static string[][2] _cache;

    if (!dbFile.among(_dbFile[0], _dbFile[1]))
    {
        // Select least recently used DB cache
        _dbIdx = 1 - _dbIdx;
        _firstRecord[_dbIdx] = 0;
        _dbFile[_dbIdx] = dbFile;
        _numRecords[_dbIdx] = cast(id_t) getNumContigs(dbFile);
        _cache[_dbIdx].length = 0;
    }

    _dbIdx = _dbFile[0] == dbFile ? 0 : 1;
    assert(_dbFile[_dbIdx] == dbFile);

    if (
        recordNumber >= _firstRecord[_dbIdx] + _cache[_dbIdx].length ||
        recordNumber < _firstRecord[_dbIdx]
    )
    {
        enum string[] dbdumpOptions = [DBdumpOptions.sequenceString];
        _cache[_dbIdx].length = cacheSize;
        auto dbdumpLines = dbdump(
            dbFile,
            recordNumber,
            min(recordNumber + cacheSize - 1, _numRecords[_dbIdx]),
            dbdumpOptions,
        );
        scope (exit) dbdumpLines.destroy();
        auto bufferRest = readSequences(dbdumpLines).copy(_cache[_dbIdx]);
        _cache[_dbIdx] = _cache[_dbIdx][0 .. $ - bufferRest.length];
        _firstRecord[_dbIdx] = recordNumber;
        enforce!DazzlerCommandException(_cache[_dbIdx].length > 0, "cannot read sequence: empty dump");
    }


    return _cache[_dbIdx][recordNumber - _firstRecord[_dbIdx]];
}

private auto readSequences(R)(R dbdump)
{
    enum baseLetters = AliasSeq!('A', 'C', 'G', 'N', 'T', 'a', 'c', 'g', 'n', 't');

    return dbdump
        .filter!(dumpLine => dumpLine[0] == 'S')
        .map!(dumpLine => dumpLine.find!(among!baseLetters));
}

unittest
{
    import std.algorithm : equal;

    enum testDbDump = q"EOF
        + R 2
        + M 0
        + S 100
        @ S 100
        S 58 tgtgatatcggtacagtaaaccacagttgggtttaaggagggacgatcaacgaacacc
        S 42 atgccaactactttgaacgcgccgcaaggcacaggtgcgcct
EOF".outdent;

    size_t[] recordIds = [];
    auto fastaSequences = readSequences(testDbDump.lineSplitter);
    assert(fastaSequences.equal([
        "tgtgatatcggtacagtaaaccacagttgggtttaaggagggacgatcaacgaacacc",
        "atgccaactactttgaacgcgccgcaaggcacaggtgcgcct",
    ]));
}

/**
    Get the designated set of records in FASTA format. If recordNumbers is
    empty the whole DB will be converted.
*/
auto getFastaEntries(Options, Range)(in string dbFile, Range recordNumbers, in Options options)
        if (isIntegral!(typeof(options.fastaLineWidth)) &&
            isSomeString!(typeof(options.tmpdir)) &&
            isInputRange!Range && is(ElementType!Range : size_t))
{
    string[] dbdumpOptions = [
        DBdumpOptions.readNumber,
        DBdumpOptions.originalHeader,
        DBdumpOptions.sequenceString,
    ];

    return readDbDumpForFastaEntries(
        dbdump(dbFile, recordNumbers, dbdumpOptions),
        recordNumbers,
        options.fastaLineWidth,
    );
}

private auto readDbDumpForFastaEntries(S, Range)(S dbDump, Range recordNumbers, in size_t lineLength)
        if (isInputRange!S && isSomeString!(ElementType!S)
            && isInputRange!Range && is(ElementType!Range : size_t))
{
    import std.algorithm : count, filter, sort;
    import std.array : appender;
    import std.range : chunks, drop;

    enum lineSeparator = '\n';
    enum subrecordSeparator = ';';
    enum recordFormat = "R %d;H %d %s;L %d %d %d;S %d %s";
    enum numRecordLines = recordFormat.count(subrecordSeparator) + 1;

    /// Build chunks of numRecordLines lines.
    alias byRecordSplitter = dbDump => dbDump.drop(6).arrayChunks(numRecordLines);
    /// Parse chunks of numRecordLines lines into FASTA format.
    alias parseRecord = recordLines => {
        size_t recordNumber;
        size_t headerLineLength;
        string headerLine;
        size_t locationWell;
        size_t locationPulseStart;
        size_t locationPulseEnd;
        size_t sequenceLength;
        string sequence;

        auto joinedLines = recordLines.joiner(only(subrecordSeparator)).array;

        int numMatches = joinedLines
            .formattedRead!recordFormat(
                recordNumber,
                headerLineLength,
                headerLine,
                locationWell,
                locationPulseStart,
                locationPulseEnd,
                sequenceLength,
                sequence,
            );
        assert(numMatches == 8, format!"%d matches in chunk: `%s`"(numMatches, joinedLines.array));

        bool isSkipping = recordNumbers.length > 0 && !recordNumbers.canFind(recordNumber);

        debug logJsonDebug(
            "isSkipping", isSkipping,
            "wantedRecordNumbers", recordNumbers.toJson,
            "recordNumber", recordNumber,
            "headerLine", headerLine,
        );

        // skip unwanted records
        if (isSkipping)
            return null;

        auto fastaData = appender!string;
        fastaData.reserve(headerLine.length + sequence.length + sequence.length / lineLength + 1);

        fastaData ~= headerLine ~ lineSeparator;
        fastaData ~= sequence.chunks(lineLength).joiner(only(lineSeparator));

        return fastaData.data;
    };

    return byRecordSplitter(dbDump).map!parseRecord
        .map!"a()"
        .cache
        .filter!"a !is null";
}

unittest
{
    enum testDbDump = q"EOF
        + R 4
        + M 0
        + H 104539
        @ H 26
        + S 38
        @ S 19574
        R 1
        H 22 >Sim/1/0_14 RQ=0.975
        L 0 0 14
        S 14 ggcccaggcagccc
        R 2
        H 22 >Sim/2/0_9 RQ=0.975
        L 0 0 9
        S 9 cacattgtg
        R 3
        H 23 >Sim/3/0_11 RQ=0.975
        L 0 0 11
        S 11 gagtgcagtgg
        R 4
        H 23 >Sim/4/0_4 RQ=0.975
        L 0 0 4
        S 4 gagc
        R 5
        H 24 >Sim/5/0_60 RQ=0.975
        L 0 0 60
        S 60 gagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagc
EOF".outdent;

    {
        size_t[] recordIds = [];
        auto fastaEntries = readDbDumpForFastaEntries(testDbDump.lineSplitter, recordIds, 50).array;
        assert(fastaEntries == [
            ">Sim/1/0_14 RQ=0.975\nggcccaggcagccc",
            ">Sim/2/0_9 RQ=0.975\ncacattgtg",
            ">Sim/3/0_11 RQ=0.975\ngagtgcagtgg",
            ">Sim/4/0_4 RQ=0.975\ngagc",
            ">Sim/5/0_60 RQ=0.975\ngagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcga\ngcgagcgagc",
        ], fastaEntries.to!string);
    }
    {
        size_t[] recordIds = [1, 3];
        auto fastaEntries = readDbDumpForFastaEntries(testDbDump.lineSplitter, recordIds, 50).array;
        assert(fastaEntries == [
            ">Sim/1/0_14 RQ=0.975\nggcccaggcagccc",
            ">Sim/3/0_11 RQ=0.975\ngagtgcagtgg",
        ], fastaEntries.to!string);
    }
}

/**
    Build `outputDb` with the given set of FASTA records. If no `outputDb`
    is given a temporary `.dam` file will be created.

    Returns: DB file name
*/
string buildDamFile(Range)(Range fastaRecords, in string tmpdir, in string[] dbsplitOptions = [], Append append = No.append)
        if (isInputRange!Range && isSomeString!(ElementType!Range))
{
    enum tempDbNameTemplate = "auxiliary-XXXXXX";

    auto tempDbTemplate = buildPath(tmpdir, tempDbNameTemplate);
    auto tempDb = mkstemp(tempDbTemplate, damFileExtension);

    tempDb.file.close();
    remove(tempDb.name);

    return buildDamFile(tempDb.name, fastaRecords, dbsplitOptions);
}

/// ditto
string buildDamFile(Range)(string outputDb, Range fastaRecords, in string[] dbsplitOptions = [], Append append = No.append)
        if (isInputRange!Range && isSomeString!(ElementType!Range))
{
    assert(outputDb.endsWith(damFileExtension), "outputDb must end with " ~ damFileExtension);

    fasta2dam(outputDb, fastaRecords, append);
    dbsplit(outputDb, dbsplitOptions);

    return outputDb;
}

unittest
{
    import dentist.util.tempfile : mkdtemp;
    import std.file : rmdirRecurse, isFile;

    auto fastaRecords = [
        ">Sim/1/0_14 RQ=0.975\nggcccacccaggcagccc",
        ">Sim/3/0_11 RQ=0.975\ngagtgcgtgcagtgg",
    ];

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    scope (exit)
        rmdirRecurse(tmpDir);

    {
        string dbName = buildDamFile(fastaRecords[], tmpDir);

        assert(dbName.isFile);
        foreach (hiddenDbFile; getHiddenDbFiles(dbName))
            assert(hiddenDbFile.isFile);
    }
    {
        string wantedDbName = buildPath(tmpDir, "unit-test.dam");
        string dbName = buildDamFile(wantedDbName, fastaRecords[]);

        assert(dbName == wantedDbName);
        assert(dbName.isFile);
        foreach (hiddenDbFile; getHiddenDbFiles(dbName))
            assert(hiddenDbFile.isFile);
    }
}

/**
    Build `outputDb` with the given set of FASTA records. If no `outputDb`
    is given a temporary `.db` file will be created.

    Returns: DB file name
*/
string buildDbFile(Range)(Range fastaRecords, in string tmpdir, in string[] dbsplitOptions = [], Append append = No.append)
        if (isInputRange!Range && isSomeString!(ElementType!Range))
{
    enum tempDbNameTemplate = "auxiliary-XXXXXX";

    auto tempDbTemplate = buildPath(tmpdir, tempDbNameTemplate);
    auto tempDb = mkstemp(tempDbTemplate, dbFileExtension);

    tempDb.file.close();
    remove(tempDb.name);

    return buildDbFile(tempDb.name, fastaRecords, dbsplitOptions);
}

/// ditto
string buildDbFile(Range)(string outputDb, Range fastaRecords, in string[] dbsplitOptions = [], Append append = No.append)
        if (isInputRange!Range && isSomeString!(ElementType!Range))
{
    assert(outputDb.endsWith(dbFileExtension), "outputDb must end with " ~ dbFileExtension);

    fasta2db(outputDb, fastaRecords, append);
    dbsplit(outputDb, dbsplitOptions);

    return outputDb;
}

unittest
{
    import dentist.util.tempfile : mkdtemp;
    import std.file : rmdirRecurse, isFile;

    auto fastaRecords = [
        ">Sim/1/0_14 RQ=0.975\nggcccacccaggcagccc",
        ">Sim/3/0_11 RQ=0.975\ngagtgcgtgcagtgg",
    ];

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    scope (exit)
        rmdirRecurse(tmpDir);

    {
        string dbName = buildDamFile(fastaRecords[], tmpDir);

        assert(dbName.isFile);
        foreach (hiddenDbFile; getHiddenDbFiles(dbName))
            assert(hiddenDbFile.isFile);
    }
    {
        string wantedDbName = buildPath(tmpDir, "unit-test.dam");
        string dbName = buildDamFile(wantedDbName, fastaRecords[]);

        assert(dbName == wantedDbName);
        assert(dbName.isFile);
        foreach (hiddenDbFile; getHiddenDbFiles(dbName))
            assert(hiddenDbFile.isFile);
    }
}



enum id_t minQVCoverage = 4;


void computeQVs(in string dbFile, in string lasFile, id_t coverage = 0)
{
    computeQVs(dbFile, lasFile, [], coverage);
}


void computeQVs(in string dbFile, in string lasFile, in string[] masks, id_t coverage = 0)
{
    dascover(dbFile, lasFile, masks);
    dasqv(dbFile, lasFile, coverage);
}


/// Options for `DBdump`.
enum DBdustOptions : string
{
    // DUST algorithm window size.
    windowSize = "-w",
    // DUST algorithm threshold.
    dustThreshold = "-t",
    // Record only low-complexity intervals >= this size.
    minIntervalSize = "-m",
    // Take into account base composition bias.
    baseCompositionBias = "-b",
}

enum dbdustMaskName = "dust";

/// Run DBdust on dbFile.
@ExternalDependency("DBdust", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
void dbdust(in string dbFile, in string[] dbdustOptions)
{
    executeCommand(chain(only("DBdust"), dbdustOptions, only(dbFile.stripDbExtension)));
}

/**
    Align DB(s) to each other using `daligner`.

    Returns: path to las-file.
*/
string getDalignment(in string dbFile, in string[] dalignerOptions, in string outdir)
{
    dalign(dbFile, dalignerOptions, outdir);
    auto lasFile = getLasFile(dbFile, outdir);

    return lasFile;
}

/// ditto
string getDalignment(in string referenceDb, in string queryDb, in string[] dalignerOptions, in string outdir)
{
    dalign(referenceDb, queryDb, dalignerOptions, outdir);
    auto lasFile = getLasFile(referenceDb, queryDb, outdir);

    return lasFile;
}

/**
    Map dbFile.

    Returns: path to las-file.
*/
version(damapper) string getDamapping(
    in string refDb,
    in string queryDb,
    in string[] damapperOptions,
    in string outdir,
)
{
    damapper(refDb, queryDb, damapperOptions, outdir);
    auto lasFile = getLasFile(refDb, queryDb, outdir);

    return lasFile;
}


string chainLocalAlignments(
    in string dbFile,
    in string lasFile,
    in ChainingOptions options,
)
{
    string chainedLasFile = lasFile.stripExtension.to!string ~ "-chained.las";
    auto flatLocalAlignments = getFlatLocalAlignments(
        dbFile,
        lasFile,
        Yes.includeTracePoints,
    );

    AlignmentHeader headerData;
    headerData.maxTracePoints = flatLocalAlignments.maxTracePointCount;
    headerData.tracePointDistance = flatLocalAlignments.tracePointDistance;

    auto chainedAlignments = chainLocalAlignmentsAlgo(
        flatLocalAlignments,
        options,
    );

    chainedLasFile.writeAlignments(chainedAlignments, headerData);

    return chainedLasFile;
}


string filterPileUpAlignments(
    in string dbFile,
    in string lasFile,
    in coord_t properAlignmentAllowance,
)
{
    auto alignments = getAlignments(dbFile, lasFile, Yes.includeTracePoints);

    filterPileUpAlignments(alignments, properAlignmentAllowance);

    string filteredLasFile = lasFile.stripExtension.to!string ~ "-filtered.las";
    writeAlignments(filteredLasFile, alignments);

    return filteredLasFile;
}

void filterPileUpAlignments(
    ref AlignmentChain[] alignments,
    in coord_t properAlignmentAllowance,
)
{
    /// An alignment in a pile up is valid iff it is proper and the begin/end
    /// of both reads match.
    static bool isValidPileUpAlignment(const ref AlignmentChain ac, const coord_t allowance)
    {
        alias isLeftAnchored = () =>
            ac.beginsWith!"contigA"(allowance) && ac.beginsWith!"contigB"(allowance);
        alias isLeftProper = () =>
            ac.beginsWith!"contigA"(allowance) || ac.beginsWith!"contigB"(allowance);
        alias isRightAnchored = () =>
            ac.endsWith!"contigA"(allowance) && ac.endsWith!"contigB"(allowance);
        alias isRightProper = () =>
            ac.endsWith!"contigA"(allowance) || ac.endsWith!"contigB"(allowance);

        return ac.contigA.id != ac.contigB.id && (
            (isLeftAnchored() && isRightProper()) ||
            (isRightAnchored() && isLeftProper())
        );
    }

    foreach (ref alignment; alignments)
        alignment.disableIf(!isValidPileUpAlignment(alignment, properAlignmentAllowance));
}


/**
    Self-dalign dbFile and build consensus using daccord.

    Returns: filename of consensus DB.
*/
string getConsensus(Options)(in string dbFile, in size_t readId, in Options options)
        if (isOptionsList!(typeof(options.daccordOptions)) &&
            isOptionsList!(typeof(options.dalignerOptions)) &&
            isOptionsList!(typeof(options.dbsplitOptions)) &&
            is(typeof(options.properAlignmentAllowance) == const(coord_t)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    static struct ModifiedOptions
    {
        string[] daccordOptions;
        string[] dalignerOptions;
        string[] dbsplitOptions;
        string tmpdir;
        coord_t properAlignmentAllowance;
    }

    auto readIdx = readId - 1;
    auto consensusDb = getConsensus(dbFile, const(ModifiedOptions)(
        options.daccordOptions ~ format!"%s%d,%d"(cast(string) DaccordOptions.readInterval, readIdx, readIdx),
        options.dalignerOptions,
        options.dbsplitOptions,
        options.tmpdir,
        options.properAlignmentAllowance,
    ));

    if (consensusDb is null)
    {
        throw new Exception("empty consensus");
    }

    return consensusDb;
}

/// ditto
string getConsensus(Options)(in string dbFile, in Options options)
        if (isOptionsList!(typeof(options.daccordOptions)) &&
            isOptionsList!(typeof(options.dalignerOptions)) &&
            isOptionsList!(typeof(options.dbsplitOptions)) &&
            is(typeof(options.properAlignmentAllowance) == const(coord_t)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    dalign(dbFile, options.dalignerOptions, options.tmpdir);
    auto lasFile = getLasFile(dbFile, options.tmpdir);
    enforce!DazzlerCommandException(!lasEmpty(lasFile), "empty pre-consensus alignment");

    auto filteredLasFile = filterPileUpAlignments(dbFile, lasFile, options.properAlignmentAllowance);

    return getConsensus(dbFile, filteredLasFile, options);
}

/// ditto
string getConsensus(Options)(in string dbFile, in string filteredLasFile, in size_t readId, in Options options)
        if (isOptionsList!(typeof(options.daccordOptions)) &&
            isOptionsList!(typeof(options.dbsplitOptions)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    static struct ModifiedOptions
    {
        string[] daccordOptions;
        string[] dbsplitOptions;
        string tmpdir;
    }

    auto readIdx = readId - 1;
    auto consensusDb = getConsensus(dbFile, filteredLasFile, const(ModifiedOptions)(
        options.daccordOptions ~ format!"%s%d,%d"(cast(string) DaccordOptions.readInterval, readIdx, readIdx),
        options.dbsplitOptions,
        options.tmpdir,
    ));

    if (consensusDb is null)
    {
        throw new Exception("empty consensus");
    }

    return consensusDb;
}

/// ditto
string getConsensus(Options)(in string dbFile, in string filteredLasFile, in Options options)
        if (isOptionsList!(typeof(options.daccordOptions)) &&
            isOptionsList!(typeof(options.dbsplitOptions)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    computeIntrinsticQualityValuesForConsensus(dbFile, filteredLasFile);

    enforce!DazzlerCommandException(!lasEmpty(filteredLasFile), "empty pre-consensus alignment");

    computeErrorProfile(dbFile, filteredLasFile, options);
    auto consensusDb = daccord(dbFile, filteredLasFile, options.daccordOptions);
    dbsplit(consensusDb, options.dbsplitOptions);

    return consensusDb;
}

private void computeIntrinsticQualityValuesForConsensus(in string dbFile, in string lasFile)
{
    auto readDepth = getNumContigs(dbFile);

    computeIntrinsicQV(dbFile, lasFile, readDepth);
}

private void computeErrorProfile(Options)(in string dbFile, in string lasFile, in Options options)
        if (isOptionsList!(typeof(options.daccordOptions)))
{
    auto eProfOptions = options
        .daccordOptions
        .filter!(option => !option.startsWith(
            cast(string) DaccordOptions.produceFullSequences,
            cast(string) DaccordOptions.readsPart,
            cast(string) DaccordOptions.errorProfileFileName,
        ))
        .chain(only(DaccordOptions.computeErrorProfileOnly))
        .array;

    // Produce error profile
    silentDaccord(dbFile, lasFile, eProfOptions);
}


unittest
{
    import dentist.util.tempfile : mkdtemp;
    import std.file : rmdirRecurse, isFile;

    auto fastaRecords = [
        ">Sim/1/0_1050 RQ=0.975\nattTgaggcatcagccactgcacccagccttgtgccctttctgagagccgggaagatgctcccaggagccctcg\nggaggcttccctccggtcgtcgtggccagaattgtctcctgctcgtgtggagtcggtggcgggccaggcgaatg\nggagctaccggggctgccgctttggactgctcggcatttgccccatggggctgcacaggggcccaggctggctg\nagaatgtccctgggtccaggaggcagacggaggtacagcccagcagccaggaggtgttcaggatgttccccagt\ncagcacccgtggaggggagggaggaggcagggtgggcgaggaaggtccaacagtggacggcctgcccacaagag\nagctctgagctgggagctggcagagttgctgcaagtgggtgtgggccaggactgactgggcctgtgcacctgcc\ntggatgcatcagtggtcgtggtgctgcccgggaagggcgtgaagctccctgcagccaaggatcctggaggtgca\ngacatcacccagcccaccggacaacagcctgccctacttcgaggagctctgggcagcccagccccatgtccccc\ntcacgccccaccccacactgacaaaaagaccacaggattccaacagtccaaccagggggaggccgttgaattcg\nggggacaaccagaaacgcctgaaacagagataaagagactgatatggaaaagactgggctggcatggtggctcc\ncaactgggatcccagtgcttgtgagaggccgaggcgggaggatcacttgagcccagaagttcaagaccagcgtg\nggcaacatagtgagaccccgtctcttttaaaaatccttttttaattaggcaggcataggtagttgcgtgcctgc\nttttcccagctgctagggaggtagaggcaggagaatcacgggagtttcgaagtccaaggtcacagtgagctgtg\nattgcaccactgcactccagcctgggcaacatggcaagaccccatctctaaaagaaagaaacaagaagacatgg\nagagaaatatccaa",
        ">Sim/2/0_1050 RQ=0.975\nattagagCcatcagccactgcacccagccttgtgccctttctgagagccgggaagatgctcccaggagccctcg\nggaggcttccctccggtcgtcgtggccagaattgtctcctgctcgtgtggagtcggtggcgggccaggcgaatg\nggagctaccggggctgccgctttggactgctcggcatttgccccatggggctgcacaggggcccaggctggctg\nagaatgtccctgggtccaggaggcagacggaggtacagcccagcagccaggaggtgttcaggatgttccccagt\ncagcacccgtggaggggagggaggaggcagggtgggcgaggaaggtccaacagtggacggcctgcccacaagag\nagctctgagctgggagctggcagagttgctgcaagtgggtgtgggccaggactgactgggcctgtgcacctgcc\ntggatgcatcagtggtcgtggtgctgcccgggaagggcgtgaagctccctgcagccaaggatcctggaggtgca\ngacatcacccagcccaccggacaacagcctgccctacttcgaggagctctgggcagcccagccccatgtccccc\ntcacgccccaccccacactgacaaaaagaccacaggattccaacagtccaaccagggggaggccgttgaattcg\nggggacaaccagaaacgcctgaaacagagataaagagactgatatggaaaagactgggctggcatggtggctcc\ncaactgggatcccagtgcttgtgagaggccgaggcgggaggatcacttgagcccagaagttcaagaccagcgtg\nggcaacatagtgagaccccgtctcttttaaaaatccttttttaattaggcaggcataggtagttgcgtgcctgc\nttttcccagctgctagggaggtagaggcaggagaatcacgggagtttcgaagtccaaggtcacagtgagctgtg\nattgcaccactgcactccagcctgggcaacatggcaagaccccatctctaaaagaaagaaacaagaagacatgg\nagagaaatatccaa",
        ">Sim/3/0_1050 RQ=0.975\nattagaggcatcagccactgcacccagccttgtgccctttctgagagccgggaagatgctcccaggagccctcg\nggaggcttccctccggtcgtcgtggccagaattgtctcctgctcgtgtggagtcggtggcgggccaggcgaatg\nggagctaccggggctgccgctttggactgctcggcatttgccccatggggctgcacaggggcccaggctggctg\nagaatgtccctgggtccaggaggcagacggaggtacagcccagcagccaggaggtgttcaggatgttccccagt\ncagcacccgtggaggggagggaggaggcagggtgggcgaggaaggtccaacagtggacggcctgcccacaagag\nagctctgagctgggagctggcagagttgctgcaagtgggtgtgggccaggactgactgggcctgtgcacctgcc\ntggatgcatcagtggtcgtggtgctgcccgggaagggcgtgaagctccctgcagccaaggatcctggaggtgca\ngacatcacccagcccaccggacaacagcctgccctacttcgaggagctctgggcagcccagccccatgtccccc\ntcacgccccaccccacactgacaaaaagaccacaggattccaacagtccaaccagggggaggccgttgaattcg\nggggacaaccagaaacgcctgaaacagagataaagagactgatatggaaaagactgggctggcatggtggctcc\ncaactgggatcccagtgcttgtgagaggccgaggcgggaggatcacttgagcccagaagttcaagaccagcgtg\nggcaacatagtgagaccccgtctcttttaaaaatccttttttaattaggcaggcataggtagttgcgtgcctgc\nttttcccagctgctagggaggtagaggcaggagaatcacgggagtttcgaagtccaaggtcacagtgagctgtg\nattgcaccactgcactccagcctgggcaacatggcaagaccccatctctaaaagaaagaaacaagaagacatgg\nagagaaatatccaa",
    ];

    struct Options
    {
        string[] dbsplitOptions;
        string[] dalignerOptions;
        string[] daccordOptions;
        size_t fastaLineWidth;
        string tmpdir;
        coord_t properAlignmentAllowance;
    }

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    auto options = Options(
        [],
        [DalignerOptions.minAlignmentLength ~ "15"],
        [],
        74,
        tmpDir,
        100,
    );
    scope (exit)
        rmdirRecurse(tmpDir);

    string dbName = buildDamFile(fastaRecords[], tmpDir);
    string consensusDb = getConsensus(dbName, options);
    assert(consensusDb !is null);
    auto consensusFasta = getFastaEntries(consensusDb, cast(size_t[])[], options);
    auto expectedSequence = fastaRecords[$ - 1].lineSplitter.drop(1).joiner.array;
    auto consensusSequence = consensusFasta.front.lineSplitter.drop(1).joiner.array;

    assert(expectedSequence == consensusSequence,
            format!"expected %s but got %s"(expectedSequence, consensusSequence));
}

string getLasFile(in string dbA, in string baseDirectory)
{
    return getLasFile(dbA, null, baseDirectory);
}

string getLasFile(in string dbA, in string dbB, in string baseDirectory)
{
    alias dbName = dbFile => dbFile.baseName.stripExtension;

    enum fileTemplate = "%s/%s.%s.las";
    auto dbAName = dbName(dbA);
    auto dbBName = dbB is null ? dbAName : dbName(dbB);

    return format!fileTemplate(baseDirectory, dbAName, dbBName);
}

bool lasFileGenerated(in string dbA, in string baseDirectory)
{
    return lasFileGenerated(dbA, null, baseDirectory);
}

bool lasFileGenerated(in string dbA, in string dbB, in string baseDirectory)
{
    return getLasFile(dbA, dbB, baseDirectory).exists;
}

id_t getNumBlocks(in string damFile)
{
    // see also in dazzler's DB.h:394
    //     #define DB_NBLOCK "blocks = %9d\n"  //  number of blocks
    enum blockNumFormat = "blocks = %d";
    enum blockNumFormatStart = blockNumFormat[0 .. 6];
    id_t numBlocks;
    auto matchingLines = File(damFile.stripBlock).byLine.filter!(
            line => line.startsWith(blockNumFormatStart));

    if (matchingLines.empty)
    {
        auto errorMessage = format!"could not read the block count in `%s`"(damFile.stripBlock);
        throw new DazzlerCommandException(errorMessage);
    }

    auto matchingLine = matchingLines.front;

    if (formattedRead!blockNumFormat(matchingLine, numBlocks) != 1)
    {
        auto errorMessage = format!"could not read the block count in `%s`"(damFile.stripBlock);
        throw new DazzlerCommandException(errorMessage);
    }

    return numBlocks;
}

coord_t getContigCutoff(in string dbFile)
{
    // see also in dazzler's DB.h:394
    //     #define DB_PARAMS "size = %10lld cutoff = %9d all = %1d\n"
    enum paramsFormat = "size = %d cutoff = %d all = %d";
    enum paramsFormatStart = paramsFormat[0 .. 6];
    coord_t contigCutoff;
    size_t dummy;
    auto matchingLines = File(dbFile.stripBlock).byLine.filter!(
            line => line.startsWith(paramsFormatStart));

    if (matchingLines.empty)
    {
        enum msg = "could not read the contig cutoff in `%s`; " ~
                   "try using DBsplit to fix";
        auto errorMessage = format!msg(dbFile.stripBlock);
        throw new DazzlerCommandException(errorMessage);
    }

    auto matchingLine = matchingLines.front;

    if (matchingLine.formattedRead!paramsFormat(dummy, contigCutoff, dummy) != 3)
    {
        auto errorMessage = format!"could not read the contig cutoff in `%s`"(dbFile.stripBlock);
        throw new DazzlerCommandException(errorMessage);
    }

    return contigCutoff;
}

id_t getNumContigs(in string damFile)
{
    return getNumContigs(damFile, No.untrimmedDb);
}

id_t getNumContigs(in string damFile, Flag!"untrimmedDb" untrimmedDb = No.untrimmedDb)
{
    enum contigNumFormat = "+ R %d";
    enum contigNumFormatStart = contigNumFormat[0 .. 4];
    id_t numContigs;
    id_t[] empty;
    string[] options = untrimmedDb ? [DBdumpOptions.untrimmedDatabase] : [];
    auto dbdumpLines = dbdump(damFile, empty, options);
    scope (exit) dbdumpLines.destroy();
    auto matchingLine = dbdumpLines
        .filter!(line => line.startsWith(contigNumFormatStart))
        .front;

    if (!matchingLine)
    {
        auto errorMessage = format!"could not read the contig count in `%s`"(damFile);
        throw new DazzlerCommandException(errorMessage);
    }

    if (formattedRead!contigNumFormat(matchingLine, numContigs) != 1)
    {
        auto errorMessage = format!"could not read the contig count in `%s`"(damFile);
        throw new DazzlerCommandException(errorMessage);
    }

    return numContigs;
}

auto getScaffoldStructure(in string damFile)
{
    enum string[] dbshowOptions = [DBshowOptions.noSequence];

    auto rawScaffoldInfo = dbshow(damFile, dbshowOptions);

    return ScaffoldStructureReader(rawScaffoldInfo);
}

alias ScaffoldSegment = Algebraic!(ContigSegment, GapSegment);

struct ContigSegment
{
    size_t globalContigId;
    size_t scaffoldId;
    size_t contigId;
    size_t begin;
    size_t end;
    string header;

    invariant
    {
        assert(begin < end);
    }

    @property size_t length() const pure nothrow
    {
        return end - begin;
    }
}

struct GapSegment
{
    size_t beginGlobalContigId;
    size_t endGlobalContigId;
    size_t scaffoldId;
    size_t beginContigId;
    size_t endContigId;
    size_t begin;
    size_t end;

    invariant
    {
        assert(begin < end);
    }

    @property size_t length() const pure nothrow
    {
        return end - begin;
    }
}

private struct ScaffoldStructureReader
{
    static enum scaffoldInfoLineSeparator = " :: ";
    static enum scaffoldInfoLineFormat = "Contig %d[%d,%d]";
    alias RawScaffoldInfo = typeof("".lineSplitter);

    private RawScaffoldInfo rawScaffoldInfo;
    private ContigSegment lastContigPart;
    private ScaffoldSegment currentPart;
    private bool _empty;

    this(string rawScaffoldInfo)
    {
        this.rawScaffoldInfo = rawScaffoldInfo.lineSplitter;
        // Force the first element to be a contigPart.
        this.currentPart = GapSegment();
        // Make `scaffoldId`s start at 0.
        this.lastContigPart.scaffoldId = -1UL;

        if (!empty)
        {
            popFront();
        }
    }

    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty ScaffoldStructureReader");

        if (rawScaffoldInfo.empty)
        {
            _empty = true;

            return;
        }

        auto nextContigPart = ContigSegment(
            lastContigPart.globalContigId + 1,
            lastContigPart.scaffoldId,
        );

        auto currentLine = rawScaffoldInfo.front;
        auto parts = rawScaffoldInfo.front.split(scaffoldInfoLineSeparator);

        assert(parts.length == 2, "illegally formatted line from `DBshow -n`");
        nextContigPart.header = parts[0];
        parts[1].formattedRead!scaffoldInfoLineFormat(
            nextContigPart.contigId,
            nextContigPart.begin,
            nextContigPart.end,
        );
        auto hasHeaderChanged = lastContigPart.header != nextContigPart.header;

        if (hasHeaderChanged)
            ++nextContigPart.scaffoldId;

        if (currentPart.peek!GapSegment !is null
                || lastContigPart.scaffoldId != nextContigPart.scaffoldId)
        {
            lastContigPart = nextContigPart;
            currentPart = nextContigPart;
            rawScaffoldInfo.popFront();
        }
        else
        {
            currentPart = GapSegment(
                lastContigPart.globalContigId,
                nextContigPart.globalContigId,
                lastContigPart.scaffoldId,
                lastContigPart.contigId,
                nextContigPart.contigId,
                lastContigPart.end,
                nextContigPart.begin,
            );
        }
    }

    @property ScaffoldSegment front() const
    {
        assert(!empty, "Attempting to fetch the front of an empty ScaffoldStructureReader");
        return currentPart;
    }

    @property bool empty() const
    {
        return _empty;
    }

    ScaffoldStructureReader save()
    {
        ScaffoldStructureReader copy;

        copy.rawScaffoldInfo = this.rawScaffoldInfo.save;
        copy.lastContigPart = this.lastContigPart;
        copy.currentPart = this.currentPart;
        copy._empty = this._empty;

        return copy;
    }
}

unittest
{
    auto exampleDump = q"EOS
>reference_mod/1/0_837550 RQ=0.850 :: Contig 0[0,8300]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 1[12400,20750]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 2[29200,154900]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 3[159900,169900]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 4[174900,200650]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 5[203650,216400]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 6[218900,235150]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 7[238750,260150]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 8[263650,837500]
>reference_mod/2/0_1450 RQ=0.850 :: Contig 0[0,1450]
EOS";

    auto reader = ScaffoldStructureReader(exampleDump);
    auto scaffoldStructure = reader.array;

    assert(scaffoldStructure == [
        ScaffoldSegment(ContigSegment(
            1, 0, 0, 0, 8300,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldSegment(GapSegment(1, 2, 0, 0, 1, 8300, 12400)),
        ScaffoldSegment(ContigSegment(
            2, 0, 1, 12400, 20750,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldSegment(GapSegment(2, 3, 0, 1, 2, 20750, 29200)),
        ScaffoldSegment(ContigSegment(
            3, 0, 2, 29200, 154900,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldSegment(GapSegment(3, 4, 0, 2, 3, 154900, 159900)),
        ScaffoldSegment(ContigSegment(
            4, 0, 3, 159900, 169900,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldSegment(GapSegment(4, 5, 0, 3, 4, 169900, 174900)),
        ScaffoldSegment(ContigSegment(
            5, 0, 4, 174900, 200650,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldSegment(GapSegment(5, 6, 0, 4, 5, 200650, 203650)),
        ScaffoldSegment(ContigSegment(
            6, 0, 5, 203650, 216400,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldSegment(GapSegment(6, 7, 0, 5, 6, 216400, 218900)),
        ScaffoldSegment(ContigSegment(
            7, 0, 6, 218900, 235150,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldSegment(GapSegment(7, 8, 0, 6, 7, 235150, 238750)),
        ScaffoldSegment(ContigSegment(
            8, 0, 7, 238750, 260150,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldSegment(GapSegment(8, 9, 0, 7, 8, 260150, 263650)),
        ScaffoldSegment(ContigSegment(
            9, 0, 8, 263650, 837500,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldSegment(ContigSegment(
            10, 1, 0, 0, 1450,
            ">reference_mod/2/0_1450 RQ=0.850",
        )),
    ]);
}

/**
    Get the hidden files comprising the designated mask.
*/
auto getMaskFiles(in string dbFile, in string maskName, Flag!"allowBlock" allowBlock = No.allowBlock)
{
    enum errorDotNotAllowed = "mask name must not contain dots `.`";

    enforce!DazzlerCommandException(!maskName.canFind("/"), "mask name must not contain slashes `/`");

    auto maskNameParts = maskName.split(".");

    enforce!DazzlerCommandException(maskNameParts.length <= 2, errorDotNotAllowed);

    id_t block;
    if (maskNameParts.length == 2)
    {
        try
        {
            block = maskNameParts[0].to!id_t;
        }
        catch (ConvException e)
        {
            throw new DazzlerCommandException(allowBlock
                ? "ill-formed block mask name: block part must be integral"
                : errorDotNotAllowed
            );
        }

        enforce!DazzlerCommandException(block > 0, allowBlock
            ? "ill-formed block mask name: block part must be positive"
            : errorDotNotAllowed
        );

        enforce!DazzlerCommandException(allowBlock, "block mask not allowed");
    }

    auto destinationDir = dbFile.dirName;
    auto dbName = dbFile.baseName.stripExtension;
    auto maskHeader = format!"%s/.%s.%s.anno"(destinationDir, dbName, maskName);
    auto maskData = format!"%s/.%s.%s.data"(destinationDir, dbName, maskName);

    return tuple!("header", "data")(maskHeader, maskData);
}

/// Thrown on failure while reading a Dazzler mask.
///
/// See_Also: `readMask`
class MaskReaderException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

private
{
    alias MaskHeaderEntry = int;
    alias MaskDataPointer = long;
    alias MaskDataEntry = int;
}

/**
    Read the `Region`s of a Dazzler mask for `dbFile`.

    Throws: MaskReaderException
    See_Also: `writeMask`, `getMaskFiles`
*/
Region[] readMask(Region)(in string dbFile, in string maskName)
{
    alias _enforce = enforce!MaskReaderException;

    auto maskFileNames = getMaskFiles(dbFile, maskName, Yes.allowBlock);
    auto maskHeader = readMaskHeader(maskFileNames.header);
    auto maskData = getBinaryFile!MaskDataEntry(maskFileNames.data);

    auto maskRegions = appender!(Region[]);
    alias RegionContigId = typeof(maskRegions.data[0].tag);
    alias RegionBegin = typeof(maskRegions.data[0].begin);
    alias RegionEnd = typeof(maskRegions.data[0].end);
    auto numReads = getNumContigs(dbFile, No.untrimmedDb).to!int;
    id_t[] trimmedDbTranslateTable;

    size_t currentContig = 1;

    if (dbFile.endsWith(damFileExtension) && maskHeader.numReads > numReads)
    {
        logJsonWarn(
            "info", "reading mask for untrimmed DB",
            "dbFile", dbFile,
            "maskName", maskName,
        );
        numReads = getNumContigs(dbFile, Yes.untrimmedDb).to!int;
        trimmedDbTranslateTable = getTrimmedDbTranslateTable(dbFile);
    }

    if (maskHeader.numReads != numReads)
        logJsonWarn(
            "info", "mask does not match DB: number of reads does not match",
            "dbFile", dbFile,
            "maskName", maskName,
        );
    _enforce(maskHeader.size == 0, "corrupted mask: expected 0");
    _enforce(maskHeader.dataPointers.length == maskHeader.numReads + 1,
             "corrupted mask: unexpected number of data pointers");

    foreach (dataPtrRange; maskHeader.dataPointers[].slide!(No.withPartial)(2))
    {
        auto dataPtrs = dataPtrRange.map!(ptr => ptr / MaskDataEntry.sizeof)
            .takeExactly!2;
        _enforce(0 <= dataPtrs[0] && dataPtrs[0] <= dataPtrs[1]
                && dataPtrs[1] <= maskData.length, "corrupted mask: data pointer out of bounds");
        _enforce(dataPtrs[0] % 2 == 0 && dataPtrs[1] % 2 == 0,
                "corrupted mask: non-sense data pointers");

        foreach (interval; maskData[dataPtrs[0] .. dataPtrs[1]].chunks(2))
        {
            enforce!MaskReaderException(interval.length == 2 && 0 <= interval[0]
                    && interval[0] <= interval[1], "corrupted mask: invalid interval");

            Region newRegion;
            newRegion.tag = currentContig.to!RegionContigId;
            newRegion.begin = interval[0].to!RegionBegin;
            newRegion.end = interval[1].to!RegionEnd;

            if (trimmedDbTranslateTable.length > 0)
                newRegion.tag = trimmedDbTranslateTable[newRegion.tag - 1].to!RegionContigId;

            if (newRegion.tag < id_t.max)
                maskRegions ~= newRegion;
        }

        ++currentContig;
    }

    return maskRegions.data;
}

private auto readMaskHeader(in string fileName)
{
    auto headerFile = File(fileName, "rb");

    return readMaskHeader(headerFile);
}

private auto readMaskHeader(File headerFile, Flag!"readPointers" readPointers = Yes.readPointers)
{
    import core.stdc.stdio : SEEK_CUR;

    MaskHeaderEntry[2] headerBuffer;

    enforce!MaskReaderException(headerFile.rawRead(headerBuffer).length == headerBuffer.length,
            format!"error while reading mask header `%s`: file too short"(headerFile.name));

    auto numPointers = headerBuffer[0] + 1;
    auto pointerBuffer = readPointers
        ? uninitializedArray!(MaskDataPointer[])(numPointers)
        : [];

    if (readPointers)
        enforce!MaskReaderException(headerFile.rawRead(pointerBuffer).length == numPointers,
                format!"error while reading mask header `%s`: file too short"(headerFile.name));
    else
        headerFile.seek(MaskDataPointer.sizeof * numPointers, SEEK_CUR);

    return tuple!("numReads", "size", "dataPointers")(headerBuffer[0],
            headerBuffer[1], pointerBuffer);
}

private void skipMaskHeader(File headerFile)
{
    cast(void) readMaskHeader(headerFile, No.readPointers);
}

private T[] getBinaryFile(T)(in string fileName)
{
    auto file = File(fileName, "rb");
    auto bufferLength = file.size / T.sizeof;

    if (bufferLength == 0)
        return [];

    auto dataBuffer = file.rawRead(uninitializedArray!(T[])(bufferLength));

    enforce!MaskReaderException(dataBuffer.length == bufferLength,
            format!"error while reading binary file `%s`: expected %d bytes of data but read only %d"(
                fileName, bufferLength * T.sizeof, dataBuffer.length * T.sizeof));

    return dataBuffer;
}

private id_t[] getTrimmedDbTranslateTable(in string dbFile)
{
    assert(dbFile.endsWith(damFileExtension), "only implemented for DAM files");

    auto cutoff = getContigCutoff(dbFile);
    auto untrimmedSize = getNumContigs(dbFile, Yes.untrimmedDb);
    auto dbRecords = getDbRecords(dbFile, [
        DBdumpOptions.readNumber,
        DBdumpOptions.originalHeader,
        DBdumpOptions.untrimmedDatabase,
    ]);

    auto trimmedDbTranslateTable = uninitializedArray!(id_t[])(untrimmedSize);
    trimmedDbTranslateTable[] = id_t.max;
    id_t untrimmedReadNumber;
    foreach (dbRecord; dbRecords)
        if (dbRecord.location.length >= cutoff)
            trimmedDbTranslateTable[dbRecord.contigId - 1] = (++untrimmedReadNumber);

    return trimmedDbTranslateTable;
}

/**
    Write the list of regions to a Dazzler mask for `dbFile`.

    See_Also: `readMask`, `getMaskFiles`
*/
void writeMask(Regions)(in string dbFile, in string maskName, Regions regions)
    if (isInputRange!Regions)
{
    alias MaskRegion = Tuple!(
        MaskHeaderEntry, "tag",
        MaskDataEntry, "begin",
        MaskDataEntry, "end",
    );

    auto maskFileNames = getMaskFiles(dbFile, maskName, Yes.allowBlock);
    auto maskHeader = File(maskFileNames.header, "wb");
    auto maskData = File(maskFileNames.data, "wb");

    auto maskRegions = regions
        .map!(region => MaskRegion(
            region.tag.to!MaskHeaderEntry,
            region.begin.to!MaskDataEntry,
            region.end.to!MaskDataEntry,
        ))
        .array;
    maskRegions.sort();

    auto numReads = getNumContigs(dbFile, No.untrimmedDb).to!MaskHeaderEntry;
    MaskHeaderEntry size = 0; // Mark the DAZZ_TRACK as a mask (see DAZZ_DB/DB.c:1183)
    MaskHeaderEntry currentContig = 1;
    MaskDataPointer dataPointer = 0;

    maskHeader.rawWrite([numReads, size]);
    maskHeader.rawWrite([dataPointer]);
    foreach (maskRegion; maskRegions)
    {
        assert(maskRegion.tag >= currentContig);

        while (maskRegion.tag > currentContig)
        {
            maskHeader.rawWrite([dataPointer]);
            ++currentContig;
        }

        if (maskRegion.tag == currentContig)
        {
            maskData.rawWrite([maskRegion.begin, maskRegion.end]);
            dataPointer += typeof(maskRegion.begin).sizeof + typeof(maskRegion.end).sizeof;
        }
    }

    foreach (emptyContig; currentContig .. numReads + 1)
    {
        maskHeader.rawWrite([dataPointer]);
    }
}


enum AccumMode : int
{
    exact = 0,
    sum = 1,
}


struct DazzExtra(T) if (is(T == long) || is(T == double))
{
    static if (is(T == long))
        enum int vtype = 0;
    else static if (is(T == double))
        enum int vtype = 1;
    else
        static assert(0);

    string name;
    T[] data;
    alias data this;
    AccumMode accumMode;
}


DazzExtra!T dazzExtra(T)(string name, T[] data, AccumMode accumMode = AccumMode.init)
    if (is(T == long) || is(T == double))
{
    return DazzExtra!T(name, data, accumMode);
}


/// Thrown on failure while reading a Dazzler mask.
///
/// See_Also: `readMask`
class DazzExtraNotFound : Exception
{
    import std.exception : basicExceptionCtors;

    mixin basicExceptionCtors;
}


/**
    Read an extra from Dazzler mask for `dbFile`.

    Returns:  fully populated DazzExtra!T.
    Throws:   MaskReaderException on read errors
              DazzExtraNotFound if no extra with given name exists
    See_Also: `writeDazzExtra`, `getMaskFiles`
*/
DazzExtra!T readDazzExtra(T)(in string dbFile, in string maskName, string extraName)
    if (is(T == long) || is(T == double))
{
    import core.stdc.stdio : SEEK_CUR;

    auto maskFileNames = getMaskFiles(dbFile, maskName, Yes.allowBlock);
    auto maskHeaderFile = File(maskFileNames.header, "rb");

    skipMaskHeader(maskHeaderFile);

    size_t extraIndex;
    int vtype;
    int dataLength;
    AccumMode accumMode;
    auto currentExtraName = new char[1024];

    alias _enforce = (cond, lazy err) => enforce!MaskReaderException(
        cond,
        format!"error while reading mask extra `%s.%s` (current extra index %d): %s"(
            maskHeaderFile.name,
            extraName,
            extraIndex,
            err,
        ),
    );
    alias _safeRead = (data, lazy what) => _enforce(
        maskHeaderFile.rawRead(data).length == data.length,
        format!"unexpected end of file while reading %s"(what),
    );

    while (!maskHeaderFile.eof)
    {
        int[4] header;

        _safeRead(header[], "header data");

        vtype = header[0];
        dataLength = header[1];
        accumMode = header[2].to!AccumMode;
        auto nameLength = header[3];

        currentExtraName.reserve(nameLength);
        currentExtraName.length = nameLength;

        if (currentExtraName.length > 0)
            _safeRead(currentExtraName, "name");

        if (currentExtraName == extraName)
            break;

        maskHeaderFile.seek(T.sizeof * dataLength, SEEK_CUR);
        ++extraIndex;
    }

    enforce!DazzExtraNotFound(currentExtraName == extraName);
    _enforce(
        DazzExtra!T.vtype == vtype,
        format!"vtype does not match: expected %s but got %s"(
            DazzExtra!T.vtype.to!string,
            vtype.to!string,
        ),
    );

    auto data = uninitializedArray!(T[])(dataLength);
    _safeRead(data, "data");

    return dazzExtra(
        extraName,
        data,
        accumMode,
    );
}


/**
    Write given extras to an existing Dazzler mask for `dbFile`.

    See_Also: `readMask`, `getMaskFiles`
*/
void writeDazzExtra(T)(in string dbFile, in string maskName, DazzExtra!T extra)
    if (is(T == long) || is(T == double))
{
    auto maskFileNames = getMaskFiles(dbFile, maskName, Yes.allowBlock);
    auto maskHeader = File(maskFileNames.header, "ab");

    int[4] header = [
        extra.vtype,
        extra.length.to!int,
        cast(int) extra.accumMode,
        extra.name.length.to!int,
    ];
    maskHeader.rawWrite(header);
    maskHeader.rawWrite(extra.name);
    maskHeader.rawWrite(extra.data);
}


unittest
{
    import dentist.util.tempfile : mkdtemp;
    import std.file : rmdirRecurse;
    import std.range : iota;

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    scope (exit)
        rmdirRecurse(tmpDir);

    auto dbFile = tmpDir ~ "/test.db";
    enum maskName = "test-mask";

    buildDbFile(dbFile, [
        ">Sim/1/0_14 RQ=0.975\nggcccacccaggcagccc",
        ">Sim/3/0_11 RQ=0.975\ngagtgcgtgcagtgg",
    ]);

    alias getData(T) = () => iota(1024)
        .map!(i => cast(T) i)
        .array;

    writeMask(dbFile, maskName, ReferenceInterval[].init);
    auto intExtra = dazzExtra("int-extra", getData!long(), AccumMode.sum);
    writeDazzExtra(dbFile, maskName, intExtra);
    auto floatExtra = dazzExtra("float-extra", getData!double());
    writeDazzExtra(dbFile, maskName, floatExtra);

    auto recoveredIntExtra = readDazzExtra!long(dbFile, maskName, intExtra.name);
    auto recoveredFloatExtra = readDazzExtra!double(dbFile, maskName, floatExtra.name);

    assert(recoveredIntExtra == intExtra);
    assert(recoveredFloatExtra == floatExtra);
}


/// Specifies the type of action `withOption` should take.
enum OptionModifier : ubyte
{
    /// Replace existing option or add as new option.
    replaceOrAdd,
    /// Make sure option is present, e.i. option name and value must match.
    ensurePresent,
    /// Make sure option is present, e.i. an option with the given name exists.
    /// Does not modify existing options. This does nothing if no option value
    /// is given.
    defaultValue,
    /// Append new option.
    add,
    /// Remove existing options.
    remove,
    /// Replace existing options. Does not add new option.
    replace,
}


/// Lazily modifies `dazzlerOptions` as specified.
auto withOption(R)(R dazzlerOptions, string optionName, OptionModifier mod)
    if (isInputRange!R && is(ElementType!R == string))
{
    return dazzlerOptions.withOption(optionName, null, mod);
}

auto withOption(R)(
    R dazzlerOptions,
    string optionName,
    string value,
    OptionModifier mod,
) if (isInputRange!R && is(ElementType!R == string))
{
    return WithOptionsImpl!R(
        dazzlerOptions,
        optionName,
        value,
        mod,
    );
}

/// `replaceOrAdd`
unittest
{
    import std.algorithm : equal;

    enum dazzlerOptions = [
        "-C",
        "-n.7",
        "-e0.841500",
        "-M25",
        "-T8",
        "-P/tmp",
        "-mdust",
        "-mdentist-self",
        "-mtan",
    ];

    auto modifiedDazzlerOptions = dazzlerOptions
        .withOption("-T", "16", OptionModifier.replaceOrAdd)
        .withOption("-k", "20", OptionModifier.replaceOrAdd);

    assert(equal(modifiedDazzlerOptions, [
        "-C",
        "-n.7",
        "-e0.841500",
        "-M25",
        "-T16",
        "-P/tmp",
        "-mdust",
        "-mdentist-self",
        "-mtan",
        "-k20",
    ]));
}

/// `ensurePresent`
unittest
{
    import std.algorithm : equal;

    enum dazzlerOptions = [
        "-C",
        "-n.7",
        "-e0.841500",
        "-M25",
        "-T8",
        "-P/tmp",
        "-mdust",
        "-mdentist-self",
        "-mtan",
    ];

    auto modifiedDazzlerOptions = dazzlerOptions
        .withOption("-C", OptionModifier.ensurePresent)
        .withOption("-e", "0.841500", OptionModifier.ensurePresent)
        .withOption("-m", "dust", OptionModifier.ensurePresent)
        .withOption("-m", "dentist-reads", OptionModifier.ensurePresent);

    assert(equal(modifiedDazzlerOptions, [
        "-C",
        "-n.7",
        "-e0.841500",
        "-M25",
        "-T8",
        "-P/tmp",
        "-mdust",
        "-mdentist-self",
        "-mtan",
        "-mdentist-reads",
    ]));
}

/// `defaultValue`
unittest
{
    import std.algorithm : equal;

    enum dazzlerOptions = [
        "-C",
        "-n.7",
        "-e0.841500",
        "-M25",
        "-T8",
        "-P/tmp",
        "-mdust",
        "-mdentist-self",
        "-mtan",
    ];

    auto modifiedDazzlerOptions = dazzlerOptions
        .withOption("-e", "0.7", OptionModifier.defaultValue)
        .withOption("-s", "126", OptionModifier.defaultValue)
        .withOption("-P", null, OptionModifier.defaultValue);

    assert(equal(modifiedDazzlerOptions, [
        "-C",
        "-n.7",
        "-e0.841500",
        "-M25",
        "-T8",
        "-P/tmp",
        "-mdust",
        "-mdentist-self",
        "-mtan",
        "-s126",
    ]));
}

/// `add` and `remove`
unittest
{
    import std.algorithm : equal;

    enum dazzlerOptions = [
        "-C",
        "-n.7",
        "-e0.841500",
        "-M25",
        "-T8",
        "-P/tmp",
        "-mdust",
        "-mdentist-self",
        "-mtan",
    ];

    auto modifiedDazzlerOptions = dazzlerOptions
        .withOption("-T", OptionModifier.remove)
        .withOption("-M", "128", OptionModifier.remove)
        .withOption("-P", "/tmp", OptionModifier.remove)
        .withOption("-m", "dentist-reads", OptionModifier.add);

    assert(equal(modifiedDazzlerOptions, [
        "-C",
        "-n.7",
        "-e0.841500",
        "-M25",
        "-mdust",
        "-mdentist-self",
        "-mtan",
        "-mdentist-reads",
    ]));
}

// empty input range
unittest
{
    import std.algorithm : equal;

    enum dazzlerOptions = new string[0];

    assert(equal(
        dazzlerOptions.withOption("-I", OptionModifier.ensurePresent),
        ["-I"],
    ));
    assert(equal(
        dazzlerOptions.withOption("-T", "8", OptionModifier.replaceOrAdd),
        ["-T8", ],
    ));
    assert(equal(
        dazzlerOptions.withOption("-m", "dentist-reads", OptionModifier.add),
        ["-mdentist-reads", ],
    ));
}


struct WithOptionsImpl(R) if (isInputRange!R && is(ElementType!R == string))
{
    private R options;
    private string optionName;
    private string value;
    private OptionModifier mod;
    private string currentOption;
    private bool optionPresent;


    this(R options, string optionName, string value, OptionModifier mod)
    {
        assert(
            optionName.length > 0 && optionName[0] == '-',
            "optionName must start with a dash `-`",
        );

        this.options = options;
        this.optionName = optionName;
        this.value = value;
        this.mod = mod;

        if (!empty)
            popFront();
    }


    this(this)
    {
        static if (isForwardRange!R)
            options = options.save;
    }


    @property string newOption() const pure nothrow @safe
    {
        return optionName ~ value;
    }


    @property bool empty() const
    {
        return options.empty &&
               currentOption is null &&
               (!isAdditiveModifier || optionPresent);
    }


    @property string front() const pure nothrow @safe
    {
        assert(!empty, "Attempting to fetch the front of an empty WithOptionsImpl");

        return currentOption;
    }


    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty WithOptionsImpl");

        if (!options.empty)
        {
            auto option = options.front;

            final switch (mod)
            {
                case OptionModifier.replaceOrAdd:
                    if (option.startsWith(optionName))
                    {
                        optionPresent = true;

                        return emitOption(newOption);
                    }
                    else
                    {
                        return emitOption(option);
                    }
                case OptionModifier.remove:
                    while (shouldRemove(option))
                    {
                        options.popFront();
                        option = options.front;
                    }

                    return emitOption(option);
                case OptionModifier.replace:
                    if (option.startsWith(optionName))
                        return emitOption(newOption);
                    else
                        return emitOption(option);
                case OptionModifier.ensurePresent:
                    optionPresent |= (option == newOption);

                    return emitOption(option);
                case OptionModifier.defaultValue:
                    optionPresent |= value is null || option.startsWith(optionName);

                    return emitOption(option);
                case OptionModifier.add:
                    return emitOption(option);
            }
        }

        if (!isAdditiveModifier || optionPresent)
        {
            return emitOption(null);
        }
        else
        {
            optionPresent = true;

            return emitOption(newOption);
        }
    }


    private void emitOption(string option)
    {
        if (!options.empty)
            options.popFront();

        this.currentOption = option;
    }


    private bool shouldRemove(string option) const pure nothrow @safe
    {
        return value is null
            ? option.startsWith(optionName)
            : option == newOption;
    }


    private bool isAdditiveModifier() const pure nothrow @safe
    {
        return 0 < mod.among(
            OptionModifier.replaceOrAdd,
            OptionModifier.ensurePresent,
            OptionModifier.defaultValue,
            OptionModifier.add,
        );
    }


    static if (isForwardRange!R)
        @property auto save()
        {
            // Postblit constructor does the actual work

            return this;
        }
}


/// Options for `daccord`.
enum DaccordOptions : string
{
    /// number of threads (default 4)
    numberOfThreads = "-t",
    /// window size (default 40)
    windowSize = "-w",
    /// advance size (default 10)
    advanceSize = "-a",
    /// max depth (default 18446744073709551615)
    maxDepth = "-d",
    /// produce full sequences (default 0)
    produceFullSequences = "-f",
    /// verbosity (default 18446744073709551615)
    verbosity = "-V",
    /// read interval (default 0,18446744073709551615)
    readInterval = "-I",
    /// reads part (default 0,1)
    readsPart = "-J",
    /// error profile file name (default input.las.eprof)
    errorProfileFileName = "-E",
    /// minimum window coverage (default 3)
    minWindowCoverage = "-m",
    /// maximum window error (default 18446744073709551615)
    maxWindowError = "-e",
    /// minimum length of output (default 0)
    minLengthOfOutput = "-l",
    /// minimum k-mer filter frequency (default 0)
    minKMerFilterFrequency = "--minfilterfreq",
    /// maximum k-mer filter frequency (default 2)
    maxKMerFilterFrequency = "--maxfilterfreq",
    /// temporary file prefix (default daccord_ozelot_4500_1529654843)
    temporaryFilePrefix = "-T",
    /// maximum number of alignments considered per read (default 5000)
    maxAlignmentsPerRead = "-D",
    /// maximum number of alignments considered per read (default 0)
    maxAlignmentsPerReadVard = "--vard",
    /// compute error profile only (default disable)
    computeErrorProfileOnly = "--eprofonly",
    /// compute error distribution estimate (default disable)
    computeErrorDistributionEstimate = "--deepprofileonly",
    /// kmer size (default 8)
    kmerSize = "-k",
}

/// Options for `daligner`.
enum DalignerOptions : string
{
    verbose = "-v",
    /// If the -b option is set, then the daligner assumes the data has a
    /// strong compositional bias (e.g. >65% AT rich).
    strongCompositionalBias = "-b",
    /// If the -A option is set (“A” for “asymmetric”) then just overlaps
    /// where the a-read is in block X and the b-read is in block Y are
    /// reported, and if X = Y then it further reports only those overlaps
    /// where the a-read index is less than the b-read index.
    asymmetric = "-A",
    /// If the -I option is set (“I” for “identity”) then when X = Y, overlaps
    /// between different portions of the same read will also be found and
    /// reported.
    identity = "-I",
    /// Search code looks for a pair of diagonal bands of width 2^^w
    /// (default 26 = 64) that contain a collection of exact matching k-mers
    /// (default 14) between the two reads, such that the total number of
    /// bases covered by the k-mer hits is h (default 35).
    kMerSize = "-k",
    /// ditto
    bandWidth = "-w",
    /// ditto
    hitBaseCoverage = "-h",
    /// Modimer percentage (take % of the k-mers)
    modimerPercentage = "-%",
    /// Suppresses the use of any k-mer that occurs more than t times in
    /// either the subject or target block.
    maxKmerOccurence = "-t",
    /// Let the program automatically select a value of t that meets a given
    /// memory usage limit specified (in Gb) by the -M parameter.
    maxKmerMemory = "-M",
    /// Bridge consecutive aligned segments into one if possible
    bridge = "-B",
    tempDir = "-P",
    /// Searching for local alignments involving at least -l base pairs
    /// (default 1000) or more, that have an average correlation rate of
    /// -e (default 70%).
    minAlignmentLength = "-l",
    /// ditto
    averageCorrelationRate = "-e",
    /// The local alignments found will be output in a sparse encoding where
    /// a trace point on the alignment is recorded every -s base pairs of
    /// the a-read (default 100bp).
    tracePointDistance = "-s",
    /// By setting the -H parameter to say N, one alters daligner so that it
    /// only reports overlaps where the a-read is over N base-pairs long.
    minAReadLength = "-H",
    /// The program runs with 4 threads by default, but this may be set to
    /// any power of 2 with the -T option.
    numThreads = "-T",
    /// If there are one or more interval tracks specified with the -m option
    /// (m for mask), then the reads of the DB or DB’s to which the track
    /// applies are soft masked with the union of the intervals of all the
    /// interval tracks that apply, that is any k-mers that contain any bases
    /// in any of the masked intervals are ignored for the purposes of seeding
    /// a match.
    masks = "-m",
    /// sort .las by A-read,A-position pairs for map usecase;
    /// off => sort .las by A,B-read pairs for overlap piles
    sortMap = "-a",
}


/// Options for `datander`.
enum DatanderOptions : string
{
    verbose = "-v",
    /// Search code looks for a pair of diagonal bands of width 2^^w
    /// (default 26 = 64) that contain a collection of exact matching k-mers
    /// (default 12) between the two reads, such that the total number of
    /// bases covered by the k-mer hits is h (default 35).
    kMerSize = "-k",
    /// ditto
    bandWidth = "-w",
    /// ditto
    hitBaseCoverage = "-h",
    tempDir = "-P",
    /// Searching for local alignments involving at least -l base pairs
    /// (default 1000) or more, that have an average correlation rate of
    /// -e (default 70%).
    minAlignmentLength = "-l",
    /// ditto
    averageCorrelationRate = "-e",
    /// The local alignments found will be output in a sparse encoding where
    /// a trace point on the alignment is recorded every -s base pairs of
    /// the a-read (default 100bp).
    tracePointDistance = "-s",
    /// The program runs with 4 threads by default, but this may be set to
    /// any power of 2 with the -T option.
    numThreads = "-T",
}

/// Options for `damapper`.
version (damapper) enum DamapperOptions : string
{
    verbose = "-v",
    /// If the -b option is set, then the daligner assumes the data has a
    /// strong compositional bias (e.g. >65% AT rich).
    strongCompositionalBias = "-b",
    /// Search code looks for a pair of diagonal bands of width 2^^w
    /// (default 26 = 64) that contain a collection of exact matching k-mers
    /// (default 14) between the two reads, such that the total number of
    /// bases covered by the k-mer hits is h (default 35).
    kMerSize = "-k",
    /// Suppresses the use of any k-mer that occurs more than t times in
    /// either the subject or target block.
    maxKmerOccurence = "-t",
    /// Let the program automatically select a value of t that meets a given
    /// memory usage limit specified (in Gb) by the -M parameter.
    maxKmerMemory = "-M",
    /// ditto
    averageCorrelationRate = "-e",
    /// The local alignments found will be output in a sparse encoding where
    /// a trace point on the alignment is recorded every -s base pairs of
    /// the a-read (default 100bp).
    tracePointDistance = "-s",
    /// The program runs with 4 threads by default, but this may be set to
    /// any power of 2 with the -T option.
    numThreads = "-T",
    tempDir = "-P",
    /// If there are one or more interval tracks specified with the -m option
    /// (m for mask), then the reads of the DB or DB’s to which the track
    /// applies are soft masked with the union of the intervals of all the
    /// interval tracks that apply, that is any k-mers that contain any bases
    /// in any of the masked intervals are ignored for the purposes of seeding
    /// a match.
    masks = "-m",
    /// If the -n option is given then all chains that are within the given
    /// fraction of the best are also reported, e.g. -n.95 reports all
    /// matches within 95% of the top match.
    bestMatches = "-n",
    /// The -p option requests that damapper produce a repeat profile track
    /// for each read.
    repeatProfileTrack = "-p",
    /// The parameter -z asks that LAs are sorted in pile order as opposed to
    /// map order (see the -a option of daligner for which this is the
    /// negation).
    sortPileOrder = "-z",
    /// If the -C option is set, then damapper also outputs a file Y.X.las
    /// for a given block pair that contains all the same matches as in
    /// X.Y.las but where the A-read is a contig of the reference and the
    /// B-read is a mapped read. And if the -N options is set, then the file
    /// Y.X.las is not produced.
    symmetric = "-C",
    /// ditto
    oneDirection = "-N",
}

/// Options for `DBdump`.
enum DBdumpOptions : string
{
    readNumber = "-r",
    originalHeader = "-h",
    sequenceString = "-s",
    sNROfACGTChannels = "-a",
    intrinsicQualityVector = "-i",
    quivaValues = "-q",
    repeatProfileVector = "-p",
    masks = "-m",
    untrimmedDatabase = "-u",
    upperCase = "-U",
}

/// Options for `DBshow`.
enum DBshowOptions : string
{
    untrimmedDatabase = "-u",
    showQuiva = "-q",
    showArrowPulseSequence = "-a",
    noSequence = "-n",
    masks = "-m",
    produceQuivaFile = "-Q",
    produceArrowFile = "-A",
    upperCase = "-U",
    fastaLineWidth = "-w",
}

/// Options for `fasta2DAM` and `fasta2DB`.
enum Fasta2DazzlerOptions : string
{
    verbose = "-v",
    /// Import files listed 1/line in given file.
    fromFile = "-f",
    /// Import data from stdin, use optional name as data source.
    fromStdin = "-i",
}

/// Options for `LAdump`.
enum LAdumpOptions : string
{
    coordinates = "-c",
    numDiffs = "-d",
    tracePoints = "-t",
    lengths = "-l",
    properOverlapsOnly = "-o",
}

/// Options for `computeintrinsicqv`.
enum ComputeIntrinsicQVOptions : string
{
    /// Read depth aka. read coverage. (mandatory)
    readDepth = "-d",
}

/// Options for `DBsplit`.
enum DbSplitOptions : string
{
    /// Target size of blocks (in Mbp).
    blockSize = "-s",
    /// Trimmed DB has reads >= this threshold.
    minReadLength = "-x",
    /// Trimmed DB contains all reads from a well (not just longest).
    allReads = "-a",
    /// Force the split to occur even if already split.
    force = "-f",
    /// Set primary read for a well to be the longest.
    onlyLongest = "-l",
    /// Set primary read for a well to be the median.
    onlyMedian = "-m",
}

private
{
    auto readLasHeader(in string lasFile)
    {
        long numParts;
        int tracePointSpacing;

        auto las = File(lasFile, "rb");

        enforce!DazzlerCommandException(
            las.rawRead((&numParts)[0 .. 1]).length == 1,
            "corrupted las header: could not read number of local alignments",
        );
        enforce!DazzlerCommandException(
            las.rawRead((&tracePointSpacing)[0 .. 1]).length == 1,
            "corrupted las header: could not read trace point spacing",
        );

        return tuple!(
            "numParts",
            "tracePointSpacing",
        )(
            numParts.to!size_t,
            tracePointSpacing.to!trace_point_t,
        );
    }

    unittest
    {
        import dentist.util.tempfile : mkdtemp;
        import std.file : rmdirRecurse;

        alias LocalAlignment = AlignmentChain.LocalAlignment;
        enum complement = AlignmentFlag.complement;

        auto tmpDir = mkdtemp("./.unittest-XXXXXX");
        scope (exit)
            rmdirRecurse(tmpDir);

        auto lasFile = buildPath(tmpDir, "test.las");
        enum tracePointSpacing = 1337;

        lasFile.writeAlignments([
            AlignmentChain(
                0,
                Contig(1, 13),
                Contig(2, 15),
                AlignmentFlags(),
                [
                    LocalAlignment(
                        Locus(3, 4),
                        Locus(5, 6),
                        7,
                    ),
                    LocalAlignment(
                        Locus(12, 13),
                        Locus(14, 15),
                        16,
                    ),
                ],
                tracePointSpacing,
            ),
            AlignmentChain(
                1,
                Contig(19, 31),
                Contig(20, 33),
                AlignmentFlags(complement),
                [
                    LocalAlignment(
                        Locus(21, 22),
                        Locus(23, 24),
                        25,
                    ),
                    LocalAlignment(
                        Locus(30, 31),
                        Locus(32, 33),
                        34,
                        [
                            TracePoint(0, 1),
                        ]
                    ),
                ],
                tracePointSpacing,
            ),
        ]);

        auto h = readLasHeader(lasFile);

        assert(h.numParts == 4);
        assert(h.tracePointSpacing == tracePointSpacing);
    }


    void ensureWritableDb(string dbFile, Append append)
    {
        if (!append && dbFile.exists)
            handleExistingDb(dbFile);
    }

    void dalign(in string refDam, in string[] dalignerOpts, in string outdir)
    {
        dalign([refDam], dalignerOpts, outdir);
    }

    void dalign(in string refDam, in string readsDam, in string[] dalignerOpts, in string outdir)
    {
        dalign([refDam, readsDam], dalignerOpts, outdir);
    }

    @ExternalDependency("daligner", "DALIGNER", "https://github.com/thegenemyers/DALIGNER")
    void dalign(in string[] dbList, in string[] dalignerOpts, in string outdir)
    {
        assert(dbList.length >= 1);
        auto isSelfAlignment = dbList.length == 1;
        auto inputFiles = isSelfAlignment ? [dbList[0], dbList[0]] : dbList;
        const(string[]) absInputFiles = inputFiles.map!absolutePath.array;

        executeCommand(chain(only("daligner"), dalignerOpts, absInputFiles), outdir);
    }

    @ExternalDependency("DAScover", "DASCRUBBER", "https://github.com/thegenemyers/DASCRUBBER")
    void dascover(in string dbFile, in string lasFile, in string[] masks)
    {
        auto maskArgs = masks.map!"`-m` ~ a";

        executeCommand(chain(only("DAScover", "-v"), maskArgs, only(dbFile[], lasFile[])), null);
    }

    @ExternalDependency("DASqv", "DASCRUBBER", "https://github.com/thegenemyers/DASCRUBBER")
    void dasqv(in string dbFile, in string lasFile, in id_t coverage)
    {
        auto coverageArg = coverage > 0 ? format!"-c%d"(coverage) : null;

        executeCommand(only("DASqv", "-v", coverageArg, dbFile, lasFile), null);
    }

    version (damapper)
    void damapper(in string refDam, in string readsDam, in string[] damapperOpts, in string outdir)
    {
        damapper([refDam, readsDam], damapperOpts, outdir);
    }

    version (damapper)
    @ExternalDependency("damapper", "DAMAPPER", "https://github.com/thegenemyers/DAMAPPER")
    void damapper(in string[] dbList, in string[] damapperOpts, in string outdir)
    {
        const(string[]) absDbList = dbList.map!absolutePath.array;

        executeCommand(chain(only("damapper", DamapperOptions.symmetric),
                damapperOpts, absDbList), outdir);
    }

    @ExternalDependency("computeintrinsicqv", "daccord", "https://gitlab.com/german.tischler/daccord")
    void computeIntrinsicQV(in string dbFile, in string lasFile, in size_t readDepth)
    {
        executeCommand(chain(
            only("computeintrinsicqv"),
            only(
                ComputeIntrinsicQVOptions.readDepth ~ readDepth.to!string,
                dbFile.stripBlock,
                lasFile[],
            ),
        ));
    }

    @ExternalDependency("daccord", "daccord", "https://gitlab.com/german.tischler/daccord")
    @ExternalDependency("fasta2DAM", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    string daccord(in string dbFile, in string lasFile, in string[] daccordOpts)
    {
        alias esc = escapeShellCommand;

        auto readIntervalOptFinder = daccordOpts.find!(opt => opt.startsWith(cast(string) DaccordOptions.readInterval));
        string daccordedDb;

        if (readIntervalOptFinder.empty)
        {
            daccordedDb = dbFile.stripExtension.to!string ~ "-daccord.dam";
        }
        else
        {
            auto readInterval = readIntervalOptFinder
                .front
                .replace(",", "-");

            daccordedDb = dbFile.stripExtension.to!string ~ "-daccord" ~ readInterval ~ ".dam";
        }

        ensureWritableDb(daccordedDb, No.append);

        executeShell(chain(
            only("daccord"),
            only(esc(daccordOpts)),
            only(esc(lasFile)),
            only(esc(dbFile.stripBlock)),
            only("|"),
            only("fasta2DAM", Fasta2DazzlerOptions.fromStdin),
            only(esc(daccordedDb)),
        ));

        return daccordedDb;
    }

    @ExternalDependency("daccord", "daccord", "https://gitlab.com/german.tischler/daccord")
    void silentDaccord(in string dbFile, in string lasFile, in string[] daccordOpts)
    {
        executeCommand(chain(
            only("daccord"),
            daccordOpts,
            only(lasFile),
            only(dbFile.stripBlock),
        ));
    }

    @ExternalDependency("DBshow", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    @ExternalDependency("fasta2DAM", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    @ExternalDependency("fasta2DB", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    void buildSubsetDb(R)(in string inDbFile, in string outDbFile, R readIds, Append append)
    {
        ensureWritableDb(outDbFile, append);

        alias esc = escapeShellCommand;
        auto escapedReadIds = readIds
            .map!(to!size_t)
            .map!(to!string)
            .map!esc;
        auto fastaToDbCommand = "fasta2" ~ inDbFile.extension[1 .. $].toUpper;

        executeShell(chain(
            only("DBshow"),
            only(esc(inDbFile)),
            escapedReadIds,
            only("|"),
            only(fastaToDbCommand, Fasta2DazzlerOptions.fromStdin),
            only(esc(outDbFile)),
        ));
    }

    @ExternalDependency("fasta2DAM", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    void fasta2dam(Range)(in string outFile, Range fastaRecords, Append append = No.append)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.algorithm : each, joiner;
        import std.process : Config, pipeProcess, Redirect, wait;
        import std.range : chunks;

        ensureWritableDb(outFile, append);

        enum writeChunkSize = 1024 * 1024;
        auto command = ["fasta2DAM", Fasta2DazzlerOptions.fromStdin, outFile];

        if (shouldLog(LogLevel.diagnostic))
        {
            static if (isForwardRange!Range)
                auto input = fastaRecords
                    .save
                    .map!(record => record[0 .. min(1024, $)].toJson)
                    .array[0 .. min(1024, $)]
                    .toJson;
            else
                auto input = toJson(null);

            logJsonDiagnostic(
                "action", "execute",
                "type", "pipe",
                "command", command.map!Json.array,
                "input", input,
                "state", "pre",
            );
        }

        auto process = pipeProcess(
            ["fasta2DAM", Fasta2DazzlerOptions.fromStdin, outFile],
            Redirect.stdin,
            null, // env
            Config.none,
        );
        fastaRecords
            .filter!(fastaRecord => parseFastaRecord(fastaRecord).length >= minSequenceLength)
            .joiner(only('\n'))
            .chain("\n")
            .chunks(writeChunkSize)
            .each!(chunk => process.stdin.write(chunk.array));
        process.stdin.close();
        auto exitStatus = wait(process.pid);
        if (exitStatus != 0)
        {
            throw new DazzlerCommandException(
                    format!"command `fasta2dam` failed with exit code %d"(exitStatus));
        }

        return;
    }

    @ExternalDependency("fasta2DAM", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    void fasta2dam(in string inFile, in string outFile, Append append = No.append)
    {
        ensureWritableDb(outFile, append);

        executeCommand(only("fasta2DAM", outFile, inFile));
    }

    @ExternalDependency("fasta2DB", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    void fasta2db(Range)(in string outFile, Range fastaRecords, Append append = No.append)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.algorithm : each, joiner;
        import std.process : Config, pipeProcess, Redirect, wait;
        import std.range : chunks;

        ensureWritableDb(outFile, append);

        enum writeChunkSize = 1024 * 1024;
        auto command = ["fasta2DB", Fasta2DazzlerOptions.fromStdin, outFile];

        if (shouldLog(LogLevel.diagnostic))
        {
            static if (isForwardRange!Range)
                auto input = fastaRecords
                    .save
                    .map!(record => record[0 .. min(1024, $)].toJson)
                    .array[0 .. min(1024, $)]
                    .toJson;
            else
                auto input = toJson(null);

            logJsonDiagnostic(
                "action", "execute",
                "type", "pipe",
                "command", command.map!Json.array,
                "input", input,
                "state", "pre",
            );
        }

        auto process = pipeProcess(
            ["fasta2DB", Fasta2DazzlerOptions.fromStdin, outFile],
            Redirect.stdin,
            null, // env
            Config.none,
        );
        fastaRecords
            .filter!(fastaRecord => parseFastaRecord(fastaRecord).length >= minSequenceLength)
            .joiner(only('\n'))
            .chain("\n")
            .chunks(writeChunkSize)
            .each!(chunk => process.stdin.write(chunk.array));
        process.stdin.close();
        auto exitStatus = wait(process.pid);
        if (exitStatus != 0)
        {
            throw new DazzlerCommandException(
                    format!"command `fasta2db` failed with exit code %d"(exitStatus));
        }

        return;
    }

    @ExternalDependency("fasta2DB", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    void fasta2db(in string inFile, in string outFile, Append append = No.append)
    {
        ensureWritableDb(outFile, append);

        executeCommand(only("fasta2DB", outFile, inFile));
    }

    @ExternalDependency("DBsplit", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    void dbsplit(in string dbFile, in string[] dbsplitOptions)
    {
        executeCommand(chain(only("DBsplit"), dbsplitOptions, only(dbFile.stripBlock)));
    }

    auto ladump(in string lasFile, in string dbA, in string dbB, in string[] ladumpOpts)
    {
        return ladump(lasFile, dbA, dbB, [], ladumpOpts);
    }

    @ExternalDependency("LAdump", "DALIGNER", "https://github.com/thegenemyers/DALIGNER")
    auto ladump(
        in string lasFile,
        in string dbA,
        in string dbB,
        in id_t[] readIds,
        in string[] ladumpOpts,
    )
    {
        return executePipe(chain(
            only("LAdump"),
            ladumpOpts,
            only(
                dbA.stripBlock,
                dbB.stripBlock,
                lasFile[],
            ),
            readIds.map!(to!string),
        ));
    }

    @ExternalDependency("DBdump", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    auto dbdump(in string dbFile, in string[] dbdumpOptions)
    {
        return executePipe(chain(
            only("DBdump"),
            dbdumpOptions,
            only(dbFile),
        ));
    }

    @ExternalDependency("DBdump", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    auto dbdump(Range)(in string dbFile, Range recordNumbers,
            in string[] dbdumpOptions)
        if (
            isForwardRange!Range &&
            (
                (is(Range.length) && Range.length == 0) ||
                is(ElementType!Range : size_t)
            )
        )
    {
        version (assert)
        {
            auto numRecords = numDbRecords(dbFile);

            assert(
                recordNumbers.save.all!(n => 0 < n && n <= numRecords),
                "illegal records for `DBdump`",
            );
        }

        return executePipe(chain(
            only("DBdump"),
            dbdumpOptions,
            only(dbFile),
            recordNumbers.map!(to!string)
        ));
    }

    @ExternalDependency("DBdump", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    auto dbdump(in string dbFile, id_t firstRecord, id_t lastRecord, in string[] dbdumpOptions)
    {
        version (assert)
        {
            auto numRecords = numDbRecords(dbFile);

            assert(
                0 < firstRecord && firstRecord <= lastRecord && lastRecord <= numRecords,
                "illegal record range for `DBdump`",
            );
        }

        return executePipe(chain(
            only("DBdump"),
            dbdumpOptions,
            only(
                dbFile,
                format!"%d-%d"(firstRecord, lastRecord),
            ),
        ));
    }

    @ExternalDependency("DBshow", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    string dbshow(in string dbFile, in string contigId)
    {
        return executeCommand(only("DBshow", dbFile, contigId));
    }

    @ExternalDependency("DBshow", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
    string dbshow(in string dbFile, in string[] dbshowOptions)
    {
        return executeCommand(chain(only("DBshow"), dbshowOptions, only(dbFile)));
    }

    string executeCommand(Range)(Range command, in string workdir = null)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, execute;

        string output = command.executeWrapper!("command",
                sCmd => execute(sCmd, null, // env
                    Config.none, size_t.max, workdir));
        return output;
    }

    void executeShell(Range)(Range command, in string workdir = null)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.algorithm : joiner;
        import std.process : Config, executeShell;

        string output = command.executeWrapper!("shell",
                sCmd => executeShell(sCmd.joiner(" ").array.to!string, null, // env
                    Config.none, size_t.max, workdir));
    }

    void executeScript(Range)(Range command, in string workdir = null)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, executeShell;

        string output = command.executeWrapper!("script",
                sCmd => executeShell(sCmd.buildScriptLine, null, // env
                    Config.none, size_t.max, workdir));
    }

    string executeWrapper(string type, alias execCall, Range)(Range command)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.array : array;
        import std.algorithm : filter;
        import std.string : lineSplitter;

        auto sanitizedCommand = command.filter!"a != null".array;

        logJsonDiagnostic(
            "action", "execute",
            "type", type,
            "command", sanitizedCommand.map!Json.array,
            "state", "pre",
        );
        auto result = execCall(sanitizedCommand);
        logJsonDiagnostic(
            "action", "execute",
            "type", type,
            "command", sanitizedCommand.map!Json.array,
            "output", result
                .output[0 .. min(1024, $)]
                .lineSplitter
                .map!Json
                .array,
            "exitStatus", result.status,
            "state", "post",
        );
        if (result.status > 0)
        {
            throw new DazzlerCommandException(
                    format("process %s returned with non-zero exit code %d: %s",
                    sanitizedCommand[0], result.status, result.output));
        }

        return result.output;
    }

    string buildScriptLine(in string[] command)
    {
        return escapeShellCommand(command) ~ " | sh -sve";
    }

    string stripBlock(in string fileName)
    {
        import std.regex : ctRegex, replaceFirst;

        enum blockNumRegex = ctRegex!(`\.[1-9][0-9]*\.(dam|db)$`);

        return fileName.replaceFirst(blockNumRegex, `.$1`);
    }

    unittest
    {
        assert("foo_bar.1.dam".stripBlock == "foo_bar.dam");
        assert("foo_bar.1024.db".stripBlock == "foo_bar.db");
        assert("foo_bar.dam".stripBlock == "foo_bar.dam");
        assert("foo_bar.db".stripBlock == "foo_bar.db");
    }
}

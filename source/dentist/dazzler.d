/**
    Defines bindings and utilities to/for the dazzler commands.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.dazzler;

import dentist.common.alignments : AlignmentChain, coord_t, diff_t, id_t, trace_point_t;
import dentist.util.fasta : parseFastaRecord;
import dentist.util.log;
import dentist.util.range : arrayChunks, takeExactly;
import dentist.util.tempfile : mkstemp;
import std.algorithm :
    all,
    among,
    cache,
    canFind,
    copy,
    countUntil,
    endsWith,
    filter,
    find,
    isSorted,
    joiner,
    map,
    min,
    sort,
    splitter,
    startsWith,
    SwapStrategy,
    uniq;
import std.array : appender, Appender, array, uninitializedArray;
import std.conv : to;
import std.exception : enforce;
import std.file : exists, remove;
import std.format : format, formattedRead;
import std.meta : AliasSeq;
import std.path : absolutePath, baseName, buildPath, dirName, relativePath,
    stripExtension, withExtension;
import std.process : Config, escapeShellCommand, kill, pipeProcess,
    ProcessPipes, Redirect, wait;
import std.range : chain, chunks, drop, enumerate, generate, only, slide, takeExactly;
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
import std.variant : Algebraic;
import vibe.data.json : Json, toJson = serializeToJson;

debug import std.stdio : writeln;

/// File suffixes of hidden .db files.
private immutable hiddenDbFileSuffixes = [".bps", ".idx"];

/// File suffixes of hidden .dam files.
private immutable hiddenDamFileSuffixes = [".bps", ".hdr", ".idx"];

/// Constant holding the .db file extension.
immutable dbFileExtension = ".db";

/// Constant holding the .dam file extension.
immutable damFileExtension = ".dam";

/// The Dazzler tools require sequence of a least minSequenceLength base pairs.
immutable minSequenceLength = 14;

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

class DazzlerCommandException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

/// Determines how data should be provided in the working directory.
enum ProvideMethod
{
    copy,
    symlink,
}
enum provideMethods = __traits(allMembers, ProvideMethod);

/**
    Provide dbFile in `workdir`.

    Returns: Path of the dbFile in `workdir`.
*/
string provideDamFileInWorkdir(in string dbFile, ProvideMethod provideMethod, in string workdir)
{
    foreach (hiddenDbFile; getHiddenDbFiles(dbFile))
    {
        provideFileInWorkdir(hiddenDbFile, provideMethod, workdir);
    }

    return provideFileInWorkdir(dbFile, provideMethod, workdir);
}

/**
    Provide lasFile in `workdir`.

    Returns: Path of the lasFile in `workdir`.
*/
string provideLasFileInWorkdir(in string lasFile, ProvideMethod provideMethod, in string workdir)
{
    return provideFileInWorkdir(lasFile, provideMethod, workdir);
}

/// Provide file in `workdir`.
private string provideFileInWorkdir(in string file, ProvideMethod provideMethod, in string workdir)
{
    import std.file : copy, symlink;
    import std.path : absolutePath;

    auto fileInWorkdir = buildPath(workdir, file.baseName);

    final switch (provideMethod)
    {
    case ProvideMethod.copy:
        copy(file.absolutePath, fileInWorkdir);
        break;
    case ProvideMethod.symlink:
        symlink(file.absolutePath, fileInWorkdir);
        break;
    }

    return fileInWorkdir;
}

/// Returns true iff lasFile contains zero parts.
bool lasEmpty(in string lasFile, in string dbA, in string workdir)
{
    return lasEmpty(lasFile, dbA, null, workdir);
}

/// ditto
bool lasEmpty(in string lasFile, in string dbA, in string dbB, in string workdir)
{
    auto dumpHeader = ladump(lasFile, dbA, dbB, [], workdir);

    if (dumpHeader.empty)
    {
        return true;
    }

    size_t numParts;

    dumpHeader.front.formattedRead!"+ P %d"(numParts);

    return numParts == 0;
}

/// Build a new .dam file by using the given subset of reads in inDbFile.
string dbSubset(Options, R)(in string inDbFile, R readIds, in Options options)
        if (isSomeString!(typeof(options.workdir)) &&
            isOptionsList!(typeof(options.dbsplitOptions)))
{
    immutable outDbNameTemplate = "subset-XXXXXX";

    auto outDbTemplate = buildPath(options.workdir, outDbNameTemplate);
    auto outDb = mkstemp(outDbTemplate, damFileExtension);

    outDb.file.close();
    remove(outDb.name);
    buildSubsetDb(inDbFile, outDb.name, readIds, options.workdir);
    dbsplit(outDb.name, options.dbsplitOptions, options.workdir);

    return outDb.name;
}

AlignmentChain[] getLocalAlignments(Options)(in string dbA, in Options options)
        if (isOptionsList!(typeof(options.dalignerOptions)) &&
            isOptionsList!(typeof(options.ladumpOptions)) &&
            isSomeString!(typeof(options.workdir)))
{
    if (!lasFileGenerated(dbA, options.workdir))
    {
        dalign(dbA, options.dalignerOptions, options.workdir);
    }

    return getGeneratedAlignments(dbA, null, options);
}

AlignmentChain[] getLocalAlignments(Options)(in string dbA, in string dbB, in Options options)
        if (isOptionsList!(typeof(options.dalignerOptions)) &&
            isOptionsList!(typeof(options.ladumpOptions)) &&
            isSomeString!(typeof(options.workdir)))
{
    if (!lasFileGenerated(dbA, dbB, options.workdir))
    {
        dalign(dbA, dbB, options.dalignerOptions, options.workdir);
    }

    return getGeneratedAlignments(dbA, dbB, options);
}

void computeLocalAlignments(Options)(in string[] dbList, in Options options)
        if (isOptionsList!(typeof(options.dalignerOptions)) &&
            isSomeString!(typeof(options.workdir)))
{
    dalign(dbList, options.dalignerOptions, options.workdir);
}

AlignmentChain[] getMappings(Options)(in string dbA, in string dbB, in Options options)
        if (isOptionsList!(typeof(options.damapperOptions)) &&
            isOptionsList!(typeof(options.ladumpOptions)) &&
            isSomeString!(typeof(options.workdir)))
{
    if (!lasFileGenerated(dbA, dbB, options.workdir))
    {
        damapper(dbA, dbB, options.damapperOptions, options.workdir);
    }

    return getGeneratedAlignments(dbA, dbB, options);
}

void computeMappings(Options)(in string[] dbList, in Options options)
        if (isOptionsList!(typeof(options.damapperOptions)) &&
            isOptionsList!(typeof(options.ladumpOptions)) &&
            isSomeString!(typeof(options.workdir)))
{
    damapper(dbList, options.damapperOptions, options.workdir);
}

private AlignmentChain[] getGeneratedAlignments(Options)(
    in string dbA,
    in string dbB,
    in Options options
)
        if (isOptionsList!(typeof(options.ladumpOptions)) &&
            isSomeString!(typeof(options.workdir)))
{
    auto lasFile = getLasFile(dbA, dbB, options.workdir);

    return getAlignments(dbA, dbB, lasFile, options);
}

AlignmentChain[] getAlignments(in string dbA, in string lasFile, in string workdir)
{
    return getAlignments(dbA, null, lasFile, workdir);
}

AlignmentChain[] getAlignments(
    in string dbA,
    in string dbB,
    in string lasFile,
    in string workdir
)
{
    static enum ladumpOptions = [
        LAdumpOptions.coordinates,
        LAdumpOptions.numDiffs,
        LAdumpOptions.lengths,
    ];

    auto alignmentChains = readLasDump(ladump(
        lasFile,
        dbA,
        dbB,
        ladumpOptions,
        workdir,
    )).array;
    alignmentChains.sort!("a < b", SwapStrategy.stable);

    return alignmentChains;
}

void attachTracePoints(
    AlignmentChain*[] alignmentChains,
    in string dbA,
    in string dbB,
    in string lasFile,
    in trace_point_t tracePointDistance,
    in string workdir
)
{
    static enum ladumpOptions = [
        LAdumpOptions.coordinates,
        LAdumpOptions.tracePoints,
    ];

    // NOTE: dump only for matching A reads; better would be for matching
    //       B reads but that is not possible with `LAdump`.
    auto aReadIds = alignmentChains.map!"a.contigA.id".uniq.array;
    auto acsWithTracePoints = readLasDump(ladump(
        lasFile,
        dbA,
        dbB,
        aReadIds,
        ladumpOptions,
        workdir,
    )).array;
    acsWithTracePoints.sort!"a < b";
    assert(isSorted!"*a < *b"(alignmentChains), "alignmentChains must be sorted");

    auto numAttached = alignmentChains.attachTracePoints(acsWithTracePoints, tracePointDistance);
    assert(numAttached == alignmentChains.length,
            "missing trace point lists for some alignment chains");

    debug logJsonDebug("acsWithTracePoints", acsWithTracePoints.toJson);
    version (assert)
    {
        // ensure all alignment chains are valid...
        foreach (alignmentChain; alignmentChains)
        {
            // trigger invariant of AlignmentChain
            cast(void) alignmentChain.first;
        }
    }
}

auto fingerprint(AlignmentChain* alignmentChain) pure nothrow
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

private id_t attachTracePoints(AlignmentChain*[] alignmentChains,
        ref AlignmentChain[] acsWithTracePoints, in trace_point_t tracePointDistance) pure
{
    assert(tracePointDistance > 0);

    id_t numAlignmentChainsAffected = 0;
    id_t numLoops = 0;
    id_t i = 0;
    id_t j = 0;

    while (i < alignmentChains.length && j < acsWithTracePoints.length
            && numLoops < alignmentChains.length + acsWithTracePoints.length)
    {
        auto alignmentChain = alignmentChains[i];
        auto acWithTracePoints = &acsWithTracePoints[j];
        auto acFingerprint = alignmentChain.fingerprint;
        auto tpFingerprint = acWithTracePoints.fingerprint;

        debug logJsonDebug(
            "acFingerprint", acFingerprint.toJson,
            "tpFingerprint", tpFingerprint.toJson,
        );

        if (acFingerprint == tpFingerprint)
        {
            foreach (k, ref localAlignment; alignmentChain.localAlignments)
            {
                localAlignment.tracePoints = acWithTracePoints.localAlignments[k].tracePoints;
            }
            alignmentChain.tracePointDistance = tracePointDistance;
            numAlignmentChainsAffected = i + 1;
            ++i;
        }
        else if (tpFingerprint < acFingerprint)
        {
            ++j;
        }
        else
        {
            assert(tpFingerprint > acFingerprint);
            throw new Exception(
                    format!"missing trace point data for alignment chain: %s"(*alignmentChain));
        }

        ++numLoops;
    }

    return numAlignmentChainsAffected;
}

unittest
{
    with (AlignmentChain) with (LocalAlignment)
    {
        auto alignmentChains = [
            new AlignmentChain(
                1,
                Contig(1, 1337),
                Contig(1, 307),
                emptyFlags,
                [
                    LocalAlignment(Locus(0, 1337), Locus(0, 307), 1),
                ],
            ),
            new AlignmentChain(
                2,
                Contig(1, 1338),
                Contig(1, 309),
                emptyFlags,
                [
                    LocalAlignment(Locus(0, 1338), Locus(0, 309), 2),
                ],
            ),
            new AlignmentChain(
                3,
                Contig(1, 1340),
                Contig(3, 508),
                emptyFlags,
                [
                    LocalAlignment(Locus(0, 1339), Locus(0, 403), 3),
                    LocalAlignment(Locus(0, 1340), Locus(404, 508), 4),
                ],
            ),
        ];
        auto acsWithTracePoints = [
            AlignmentChain(
                1,
                Contig(1, 1337),
                Contig(1, 307),
                emptyFlags,
                [
                    LocalAlignment(Locus(0, 1337), Locus(0, 307), 0, [
                        TracePoint(1, 102),
                        TracePoint(2, 101),
                        TracePoint(3, 104),
                    ]),
                ],
                101
            ),
            AlignmentChain(
                2,
                Contig(1, 1338),
                Contig(1, 309),
                emptyFlags,
                [
                    LocalAlignment(Locus(0, 1338), Locus(0, 309), 0, [
                        TracePoint(4, 103),
                        TracePoint(5, 104),
                        TracePoint(6, 102),
                    ]),
                ],
                101
            ),
            AlignmentChain(
                3,
                Contig(1, 1340),
                Contig(3, 508),
                emptyFlags,
                [
                    LocalAlignment(Locus(0, 1339), Locus(0, 403), 0, [
                        TracePoint(7, 105),
                        TracePoint(8, 101),
                        TracePoint(9, 100),
                        TracePoint(10, 97),
                    ]),
                    LocalAlignment(Locus(0, 1340), Locus(404, 508), 0, [
                        TracePoint(11, 2),
                        TracePoint(12, 102),
                    ]),
                ],
                101
            ),
        ];
        auto expectedAlignmentChains = [
            AlignmentChain(
                1,
                Contig(1, 1337),
                Contig(1, 307),
                emptyFlags,
                [
                    LocalAlignment(Locus(0, 1337), Locus(0, 307), 1, [
                        TracePoint(1, 102),
                        TracePoint(2, 101),
                        TracePoint(3, 104),
                    ]),
                ],
                101
            ),
            AlignmentChain(
                2,
                Contig(1, 1338),
                Contig(1, 309),
                emptyFlags,
                [
                    LocalAlignment(Locus(0, 1338), Locus(0, 309), 2, [
                        TracePoint(4, 103),
                        TracePoint(5, 104),
                        TracePoint(6, 102),
                    ]),
                ],
                101
            ),
            AlignmentChain(
                3,
                Contig(1, 1340),
                Contig(3, 508),
                emptyFlags,
                [
                    LocalAlignment(Locus(0, 1339), Locus(0, 403), 3, [
                        TracePoint(7, 105),
                        TracePoint(8, 101),
                        TracePoint(9, 100),
                        TracePoint(10, 97),
                    ]),
                    LocalAlignment(Locus(0, 1340), Locus(404, 508), 4, [
                        TracePoint(11, 2),
                        TracePoint(12, 102),
                    ]),
                ],
                101
            ),
        ];

        assert(alignmentChains.attachTracePoints(acsWithTracePoints, 101) == 3);
        assert(alignmentChains.map!"*a".array == expectedAlignmentChains);
    }
}

private struct LasDumpLineFormatTuple
{
    char indicator;
    char subIndicator;
    string format;
}

private enum LasDumpLineFormat
{
    totalChainPartsCount = LasDumpLineFormatTuple('+', 'P', "+ P %d"),
    totalTracePointsCount = LasDumpLineFormatTuple('+', 'T', "+ T %d"),
    maxChainPartsCountPerPile = LasDumpLineFormatTuple('%', 'P', "% P %d"),
    maxTracePointsCountPerPile = LasDumpLineFormatTuple('%', 'T', "% T %d"),
    maxTracePointCount = LasDumpLineFormatTuple('@', 'T', "@ T %d"),
    chainPart = LasDumpLineFormatTuple('P', '\0', "P %d %d %c %c"),
    lengths = LasDumpLineFormatTuple('L', '\0', "L %d %d"),
    coordinates = LasDumpLineFormatTuple('C', '\0', "C %d %d %d %d"),
    numDiffs = LasDumpLineFormatTuple('D', '\0', "D %d"),
    tracePointBegin = LasDumpLineFormatTuple('T', '\0', "T %d"),
    tracePoint = LasDumpLineFormatTuple(' ', '\0', " %d %d"),
}

private enum ChainPartType
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
    debug dchar[] allowedLineTypes;

public:
    this(S lasDump)
    {
        this.lasDump = getDumpLines(lasDump);
        debug with (LasDumpLineFormat)
        {
            this.allowedLineTypes = [
                totalChainPartsCount.indicator,
                totalTracePointsCount.indicator,
                maxChainPartsCountPerPile.indicator,
                maxTracePointsCountPerPile.indicator,
                maxTracePointCount.indicator,
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

private:

    static auto getDumpLines(S lasDump)
    {
        return lasDump.enumerate(1).filter!"!a[1].empty";
    }

    void setEmpty() pure nothrow
    {
        _empty = true;
    }

    void readNextAlignmentChain()
    {
        currentAC = AlignmentChain.init;
        localAlignmentsAcc.clear();
        tracePointsAcc.clear();
        peekDumpLine();

        while (true)
        {
            debug _enforce(allowedLineTypes.canFind(currentLineType), format!"forbidden line type `%c`"(currentLineType));

            with (LasDumpLineFormat)
            {
                switch (currentLineType)
                {
                case totalChainPartsCount.indicator: goto case;
                case maxChainPartsCountPerPile.indicator:
                    break; // ignore
                case maxTracePointCount.indicator:
                    _enforce(currentLineSubType == maxTracePointCount.subIndicator, "expected `@ T` line");
                    readMaxTracePointCount();
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

    void readMaxTracePointCount()
    {
        immutable maxTracePointCountFormat = LasDumpLineFormat.maxTracePointCount.format;

        currentDumpLine[].formattedRead!maxTracePointCountFormat(maxTracePointCount);

        tracePointsAcc.reserve(maxTracePointCount);
    }

    Flag!"wasDumpPartConsumed" readChainPart()
    {
        immutable chainPartFormat = LasDumpLineFormat.chainPart.format;
        immutable yesComplement = AlignmentChain.Flags(AlignmentChain.Flag.complement);
        immutable noComplement = AlignmentChain.emptyFlags;
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
        auto startingNewChain = currentAC == AlignmentChain.init;

        if (startingNewChain)
        {
            currentAC.contigA.id = contigAID;
            currentAC.contigB.id = contigBID;
            currentAC.flags = flags;

            return Yes.wasDumpPartConsumed;
        }
        else if (isChainContinuation(contigAID, contigBID, chainPartType))
        {
            _enforce(currentAC.flags == flags, "matching both strands in one alignment chain");
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
    }

    void readLengths()
    {
        immutable lengthsFormat = LasDumpLineFormat.lengths.format;

        currentDumpLine[].formattedRead!lengthsFormat(
            currentAC.contigA.length,
            currentAC.contigB.length,
        );
    }
    void readCoordinates()
    {
        immutable coordinatesFormat = LasDumpLineFormat.coordinates.format;
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
        immutable numDiffsFormat = LasDumpLineFormat.numDiffs.format;

        currentDumpLine[].formattedRead!numDiffsFormat(
            localAlignmentsAcc.data[$ - 1].numDiffs,
        );
    }

    void readTracePointBegin()
    {
        immutable tracePointBeginFormat = LasDumpLineFormat.tracePointBegin.format;
        size_t numTracePoints;

        currentDumpLine[].formattedRead!tracePointBeginFormat(numTracePoints);

        tracePointsAcc.clear();
        tracePointsAcc.reserve(numTracePoints);
    }

    void readTracePoint()
    {
        immutable tracePointFormat = LasDumpLineFormat.tracePoint.format;
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

private auto readLasDump(S)(S lasDump)
{
    return LasDumpReader!S(lasDump);
}

unittest
{
    import std.algorithm : equal;

    immutable testLasDump = [
        "+ P 9",
        "% P 9",
        "+ T 12",
        "% T 12",
        "@ T 4",
        "P 1 2 n >",
        "L 8 9",
        "C 3 4 5 6",
        "D 7",
        "P 1 2 n -",
        "L 8 9",
        "C 12 13 14 15",
        "D 16",
        "P 19 20 c +",
        "L 26 27",
        "C 21 22 23 24",
        "D 25",
        "P 19 20 c -",
        "L 26 27",
        "C 30 31 32 33",
        "D 34",
        "T 1",
        "   0 1",
        "P 37 38 n .",
        "L 35 36",
        "C 39 40 41 42",
        "D 43",
        "P 46 47 c .",
        "L 53 54",
        "C 48 49 50 51",
        "D 52",
        "T 3",
        "   2 102",
        "   3 101",
        "   4 104",
        "P 46 47 n .",
        "L 53 54",
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

    with (AlignmentChain) with (Flag) with (LocalAlignment)
    {
        auto alignmentChains = readLasDump(testLasDump).array;
        auto expectedResult = [
            AlignmentChain(
                0,
                Contig(1, 8),
                Contig(2, 9),
                emptyFlags,
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
                Contig(19, 26),
                Contig(20, 27),
                Flags(complement),
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
                Contig(37, 35),
                Contig(38, 36),
                emptyFlags,
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
                Contig(46, 53),
                Contig(47, 54),
                Flags(complement),
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
                Contig(46, 53),
                Contig(47, 54),
                emptyFlags,
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
                Flags(complement),
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
                Flags(complement),
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
                Flags(complement),
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
}

trace_point_t getTracePointDistance(in string[] dazzlerOptions = []) pure
{
    immutable defaultTracePointDistance = 100;

    foreach (option; dazzlerOptions)
    {
        if (option.startsWith(cast(const(string)) DalignerOptions.tracePointDistance))
        {
            return option[2 .. $].to!trace_point_t;
        }
    }

    return defaultTracePointDistance;
}

unittest
{
    with (DalignerOptions)
    {
        assert(getTracePointDistance([]) == 100);
        assert(getTracePointDistance([identity]) == 100);
        assert(getTracePointDistance([
            tracePointDistance ~ "42",
            averageCorrelationRate ~ ".8",
            identity,
        ]) == 42);
        assert(getTracePointDistance([
            identity,
            tracePointDistance ~ "42",
            averageCorrelationRate ~ ".8",
        ]) == 42);
        assert(getTracePointDistance([
            averageCorrelationRate ~ ".8",
            identity,
            tracePointDistance ~ "42",
        ]) == 42);
    }
}

/**
    Get the FASTA sequences of the designated records.

    Throws: DazzlerCommandException if recordNumber is not in dbFile
*/
auto getFastaSequences(Range)(in string dbFile, Range recordNumbers, in string workdir)
        if (isForwardRange!Range && is(ElementType!Range : size_t))
{
    string[] dbdumpOptions = [DBdumpOptions.sequenceString];
    auto numRecords = recordNumbers.save.walkLength;
    auto sequences = readSequences(dbdump(dbFile, recordNumbers, dbdumpOptions, workdir));
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
auto getFastaSequence(in string dbFile, id_t recordNumber, in string workdir, in id_t cacheSize = 1024)
{
    // FIXME the cache size should limit the number of `char`s retrieved, ie. control the memory
    // requirements of this function
    static id_t _firstRecord;
    static string _dbFile;
    static id_t _numRecords;
    static string[] _cache;

    if (dbFile != _dbFile)
    {
        _dbFile = dbFile;
        _numRecords = cast(id_t) getNumContigs(dbFile, workdir);
        _cache.length = 0;
    }

    if (recordNumber >= _firstRecord + _cache.length || recordNumber < _firstRecord)
    {
        enum string[] dbdumpOptions = [DBdumpOptions.sequenceString];
        _cache.length = cacheSize;
        auto bufferRest = readSequences(dbdump(
            dbFile,
            recordNumber,
            min(recordNumber + cacheSize - 1, _numRecords),
            dbdumpOptions,
            workdir,
        )).copy(_cache);
        _cache = _cache[0 .. $ - bufferRest.length];
        _firstRecord = recordNumber;
        enforce!DazzlerCommandException(_cache.length > 0, "cannot read sequence: empty dump");
    }


    return _cache[recordNumber - _firstRecord];
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

    immutable testDbDump = q"EOF
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
            isSomeString!(typeof(options.workdir)) &&
            isInputRange!Range && is(ElementType!Range : size_t))
{
    string[] dbdumpOptions = [
        DBdumpOptions.readNumber,
        DBdumpOptions.originalHeader,
        DBdumpOptions.sequenceString,
    ];

    return readDbDump(dbdump(dbFile, recordNumbers, dbdumpOptions,
            options.workdir), recordNumbers, options.fastaLineWidth);
}

private auto readDbDump(S, Range)(S dbDump, Range recordNumbers, in size_t lineLength)
        if (isInputRange!S && isSomeString!(ElementType!S)
            && isInputRange!Range && is(ElementType!Range : size_t))
{
    import std.algorithm : count, filter, sort;
    import std.array : appender;
    import std.range : chunks, drop;

    immutable lineSeparator = '\n';
    immutable subrecordSeparator = ';';
    immutable recordFormat = "R %d;H %d %s;L %d %d %d;S %d %s";
    immutable numRecordLines = recordFormat.count(subrecordSeparator) + 1;

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
    immutable testDbDump = q"EOF
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
        auto fastaEntries = readDbDump(testDbDump.lineSplitter, recordIds, 50).array;
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
        auto fastaEntries = readDbDump(testDbDump.lineSplitter, recordIds, 50).array;
        assert(fastaEntries == [
            ">Sim/1/0_14 RQ=0.975\nggcccaggcagccc",
            ">Sim/3/0_11 RQ=0.975\ngagtgcagtgg",
        ], fastaEntries.to!string);
    }
}

/// Build a .dam file with the given set of FASTA records.
string buildDamFile(Range)(Range fastaRecords, in string workdir, in string[] dbsplitOptions = [])
        if (isInputRange!Range && isSomeString!(ElementType!Range))
{
    immutable tempDbNameTemplate = "auxiliary-XXXXXX";

    auto tempDbTemplate = buildPath(workdir, tempDbNameTemplate);
    auto tempDb = mkstemp(tempDbTemplate, damFileExtension);

    tempDb.file.close();
    remove(tempDb.name);
    fasta2dam(tempDb.name, fastaRecords, workdir);
    dbsplit(tempDb.name, dbsplitOptions, workdir);

    return tempDb.name;
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

    string dbName = buildDamFile(fastaRecords[], tmpDir);

    assert(dbName.isFile);
    foreach (hiddenDbFile; getHiddenDbFiles(dbName))
    {
        assert(hiddenDbFile.isFile);
    }
}

/**
    Self-dalign dbFile and build consensus using daccord.

    Returns: list of consensus DBs.
*/
string getConsensus(Options)(in string dbFile, in size_t readId, in Options options)
        if (isOptionsList!(typeof(options.daccordOptions)) &&
            isOptionsList!(typeof(options.dalignerOptions)) &&
            isOptionsList!(typeof(options.dbsplitOptions)) &&
            isSomeString!(typeof(options.workdir)))
{
    static struct ModifiedOptions
    {
        string[] daccordOptions;
        string[] dalignerOptions;
        string[] dbsplitOptions;
        string workdir;
    }

    auto readIdx = readId - 1;
    auto consensusDb = getConsensus(dbFile, const(ModifiedOptions)(
        options.daccordOptions ~ format!"%s%d,%d"(cast(string) DaccordOptions.readInterval, readIdx, readIdx),
        options.dalignerOptions,
        options.dbsplitOptions,
        options.workdir,
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
            isSomeString!(typeof(options.workdir)))
{
    dalign(dbFile, options.dalignerOptions, options.workdir);
    computeErrorProfile(dbFile, options);

    auto lasFile = getLasFile(dbFile, options.workdir);

    if (lasEmpty(lasFile, dbFile, null, options.workdir))
    {
        return null;
    }

    auto consensusDb = daccord(dbFile, lasFile, options.daccordOptions, options.workdir);
    dbsplit(consensusDb, options.dbsplitOptions, options.workdir);

    return consensusDb;
}

private void computeErrorProfile(Options)(in string dbFile, in Options options)
        if (isOptionsList!(typeof(options.daccordOptions)) &&
            isSomeString!(typeof(options.workdir)))
{
    auto eProfOptions = options
        .daccordOptions
        .filter!(option => !option.startsWith(
            cast(string) DaccordOptions.produceFullSequences,
            //cast(string) DaccordOptions.readInterval,
            cast(string) DaccordOptions.readsPart,
            cast(string) DaccordOptions.errorProfileFileName,
        ))
        .chain(only(DaccordOptions.computeErrorProfileOnly))
        .array;
    auto lasFile = getLasFile(dbFile, options.workdir);
    enforce!DazzlerCommandException(!lasEmpty(lasFile, dbFile, null,
            options.workdir), "empty pre-consensus alignment");

    // Produce error profile
    silentDaccord(dbFile, lasFile, eProfOptions, options.workdir);
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
        string workdir;
    }

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    auto options = Options(
        [],
        [DalignerOptions.minAlignmentLength ~ "15"],
        [],
        74,
        tmpDir,
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

    immutable fileTemplate = "%s/%s.%s.las";
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
    immutable blockNumFormat = "blocks = %d";
    immutable blockNumFormatStart = blockNumFormat[0 .. 6];
    id_t numBlocks;
    auto matchingLine = File(damFile.stripBlock).byLine.filter!(
            line => line.startsWith(blockNumFormatStart)).front;

    if (!matchingLine)
    {
        auto errorMessage = format!"could not read the block count in `%s`"(damFile.stripBlock);
        throw new DazzlerCommandException(errorMessage);
    }

    if (formattedRead!blockNumFormat(matchingLine, numBlocks) != 1)
    {
        auto errorMessage = format!"could not read the block count in `%s`"(damFile.stripBlock);
        throw new DazzlerCommandException(errorMessage);
    }

    return numBlocks;
}

id_t getNumContigs(in string damFile, in string workdir)
{
    immutable contigNumFormat = "+ R %d";
    immutable contigNumFormatStart = contigNumFormat[0 .. 4];
    id_t numContigs;
    id_t[] empty;
    auto matchingLine = dbdump(damFile, empty, [], workdir)
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

auto getScaffoldStructure(Options)(in string damFile, in Options options)
        if (isSomeString!(typeof(options.workdir)))
{
    immutable string[] dbshowOptions = [DBshowOptions.noSequence];

    auto rawScaffoldInfo = dbshow(damFile, dbshowOptions, options.workdir);

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
    static immutable scaffoldInfoLineFormat = "%s:: Contig %d[%d,%d]";
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

        auto nextContigPart = ContigSegment(lastContigPart.globalContigId + 1);

        rawScaffoldInfo.front.formattedRead!scaffoldInfoLineFormat(
            nextContigPart.header,
            nextContigPart.contigId,
            nextContigPart.begin,
            nextContigPart.end,
        );
        if (lastContigPart.header != nextContigPart.header[0 .. $ - 1])
        {
            nextContigPart.scaffoldId = lastContigPart.scaffoldId + 1;
        }

        if (currentPart.peek!GapSegment !is null
                || lastContigPart.scaffoldId != nextContigPart.scaffoldId)
        {
            assert(nextContigPart.header[$ - 1] == ' ');
            // Remove the trailing space
            nextContigPart.header = nextContigPart.header[0 .. $ - 1];
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
auto getMaskFiles(in string dbFile, in string maskDestination)
{
    auto destinationDir = maskDestination.dirName;
    auto maskName = maskDestination.baseName;
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
Region[] readMask(Region)(in string dbFile, in string maskDestination, in string workdir)
{
    alias _enforce = enforce!MaskReaderException;

    auto maskFileNames = getMaskFiles(dbFile, maskDestination);
    auto maskHeader = readMaskHeader(maskFileNames.header);
    auto maskData = getBinaryFile!MaskDataEntry(maskFileNames.data);

    auto maskRegions = appender!(Region[]);
    alias RegionContigId = typeof(maskRegions.data[0].tag);
    alias RegionBegin = typeof(maskRegions.data[0].begin);
    alias RegionEnd = typeof(maskRegions.data[0].end);
    auto numReads = getNumContigs(dbFile, workdir).to!int;

    size_t currentContig = 1;

    _enforce(maskHeader.numReads == numReads, "mask does not match DB");
    _enforce(maskHeader.size == 0, "corrupted mask: expected 0");
    _enforce(maskHeader.dataPointers.length == numReads + 1,
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

            maskRegions ~= newRegion;
        }

        ++currentContig;
    }

    return maskRegions.data;
}

private auto readMaskHeader(in string fileName)
{
    auto headerFile = File(fileName, "rb");
    MaskHeaderEntry[2] headerBuffer;
    auto numPointers = (headerFile.size - headerBuffer.sizeof) / MaskDataPointer.sizeof;
    auto pointerBuffer = uninitializedArray!(MaskDataPointer[])(numPointers);

    enforce!MaskReaderException(headerFile.rawRead(headerBuffer).length == headerBuffer.length,
            format!"error while reading mask header `%s`: file too short"(fileName));
    enforce!MaskReaderException(headerFile.rawRead(pointerBuffer).length == numPointers,
            format!"error while reading mask header `%s`: file too short"(fileName));

    return tuple!("numReads", "size", "dataPointers")(headerBuffer[0],
            headerBuffer[1], pointerBuffer);
}

private T[] getBinaryFile(T)(in string fileName)
{
    auto file = File(fileName, "rb");
    auto bufferLength = file.size / T.sizeof;
    auto dataBuffer = file.rawRead(uninitializedArray!(T[])(bufferLength));

    enforce!MaskReaderException(dataBuffer.length == bufferLength,
            format!"error while reading binary file `%s`: expected %d bytes of data but read only %d"(
                fileName, bufferLength * T.sizeof, dataBuffer.length * T.sizeof));

    return dataBuffer;
}

/**
    Write the list of regions to a Dazzler mask for `dbFile`.

    See_Also: `readMask`, `getMaskFiles`
*/
void writeMask(Region)(
    in string dbFile,
    in string maskDestination,
    in Region[] regions,
    in string workdir,
)
{
    alias MaskRegion = Tuple!(
        MaskHeaderEntry, "tag",
        MaskDataEntry, "begin",
        MaskDataEntry, "end",
    );

    if (regions.length == 0)
    {
        logJsonDiagnostic(
            "notice", "skipping empty mask",
            "dbFile", dbFile,
            "maskDestination", maskDestination,
        );

        return;
    }

    auto maskFileNames = getMaskFiles(dbFile, maskDestination);
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

    auto numReads = getNumContigs(dbFile, workdir).to!MaskHeaderEntry;
    MaskHeaderEntry size = 0; // this seems to be zero always (see DAMASKER/TANmask.c:422)
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
    /// Suppresses the use of any k-mer that occurs more than t times in
    /// either the subject or target block.
    maxKmerOccurence = "-t",
    /// Let the program automatically select a value of t that meets a given
    /// memory usage limit specified (in Gb) by the -M parameter.
    maxKmerMemory = "-M",
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
}

/// Options for `damapper`.
enum DamapperOptions : string
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
enum DBdumpOptions
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
enum DBshowOptions
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
enum Fasta2DazzlerOptions
{
    verbose = "-v",
    /// Import files listed 1/line in given file.
    fromFile = "-f",
    /// Import data from stdin, use optional name as data source.
    fromStdin = "-i",
}

/// Options for `LAdump`.
enum LAdumpOptions
{
    coordinates = "-c",
    numDiffs = "-d",
    tracePoints = "-t",
    lengths = "-l",
    properOverlapsOnly = "-o",
}

private
{

    void dalign(in string refDam, in string[] dalignerOpts, in string workdir)
    {
        dalign([refDam], dalignerOpts, workdir);
    }

    void dalign(in string refDam, in string readsDam, in string[] dalignerOpts, in string workdir)
    {
        dalign([refDam, readsDam], dalignerOpts, workdir);
    }

    void dalign(in string[] dbList, in string[] dalignerOpts, in string workdir)
    {
        assert(dbList.length >= 1);
        auto isSelfAlignment = dbList.length == 1;
        auto additionalOptions = only(isSelfAlignment ? DalignerOptions.identity : null);
        auto inputFiles = isSelfAlignment ? [dbList[0], dbList[0]] : dbList;
        const(string[]) inputFilesRelativeToWorkDir = inputFiles.map!(
                f => f.relativeToWorkdir(workdir)).array;

        executeCommand(chain(only("daligner"), additionalOptions, dalignerOpts,
                inputFilesRelativeToWorkDir), workdir);
    }

    void damapper(in string refDam, in string readsDam, in string[] damapperOpts, in string workdir)
    {
        damapper([refDam, readsDam], damapperOpts, workdir);
    }

    void damapper(in string[] dbList, in string[] damapperOpts, in string workdir)
    {
        const(string[]) dbListRelativeToWorkDir = dbList.map!(
                f => f.relativeToWorkdir(workdir)).array;

        executeCommand(chain(only("damapper", DamapperOptions.symmetric),
                damapperOpts, dbListRelativeToWorkDir), workdir);
    }

    string daccord(in string dbFile, in string lasFile, in string[] daccordOpts, in string workdir)
    {
        alias esc = escapeShellCommand;
        string daccordedDb = dbFile.stripExtension.to!string ~ "-daccord.dam";

        executeShell(chain(
            only("daccord"),
            only(esc(daccordOpts)),
            only(esc(lasFile.relativeToWorkdir(workdir))),
            only(esc(dbFile.stripBlock.relativeToWorkdir(workdir))),
            only("|"),
            only("fasta2DAM", Fasta2DazzlerOptions.fromStdin),
            only(esc(daccordedDb.relativeToWorkdir(workdir))),
        ), workdir);

        return daccordedDb;
    }

    void silentDaccord(in string dbFile, in string lasFile, in string[] daccordOpts,
            in string workdir)
    {
        executeCommand(chain(
            only("daccord"),
            daccordOpts,
            only(lasFile.relativeToWorkdir(workdir)),
            only(dbFile.stripBlock.relativeToWorkdir(workdir)),
        ), workdir);
    }

    void buildSubsetDb(R)(in string inDbFile, in string outDbFile, R readIds, in string workdir)
    {
        alias esc = escapeShellCommand;
        auto escapedReadIds = readIds
            .map!(to!size_t)
            .map!(to!string)
            .map!esc;

        executeShell(chain(
            only("DBshow"),
            only(esc(inDbFile.relativeToWorkdir(workdir))),
            escapedReadIds,
            only("|"),
            only("fasta2DAM", Fasta2DazzlerOptions.fromStdin),
            only(esc(outDbFile.relativeToWorkdir(workdir))),
        ), workdir);
    }

    void fasta2dam(Range)(in string outFile, Range fastaRecords, in string workdir)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.algorithm : each, joiner;
        import std.process : Config, pipeProcess, Redirect, wait;
        import std.range : chunks;

        immutable writeChunkSize = 1024 * 1024;
        auto outFileArg = outFile.relativeToWorkdir(workdir);
        auto command = ["fasta2DAM", Fasta2DazzlerOptions.fromStdin, outFileArg];

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
            ["fasta2DAM", Fasta2DazzlerOptions.fromStdin, outFileArg],
            Redirect.stdin,
            null, // env
            Config.none,
            workdir
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

    void fasta2dam(in string inFile, in string outFile, in string workdir)
    {
        executeCommand(only("fasta2DAM", outFile.relativeToWorkdir(workdir), inFile), workdir);
    }

    void dbsplit(in string dbFile, in string[] dbsplitOptions, in string workdir)
    {
        executeCommand(chain(only("DBsplit"), dbsplitOptions,
                only(dbFile.stripBlock.relativeToWorkdir(workdir))), workdir);
    }

    auto ladump(in string lasFile, in string dbA, in string dbB,
            in string[] ladumpOpts, in string workdir)
    {
        return ladump(lasFile, dbA, dbB, [], ladumpOpts, workdir);
    }

    auto ladump(in string lasFile, in string dbA, in string dbB, in id_t[] readIds,
            in string[] ladumpOpts, in string workdir)
    {
        return executePipe(chain(
            only("LAdump"),
            ladumpOpts,
            only(
                dbA.stripBlock.relativeToWorkdir(workdir),
                dbB.stripBlock.relativeToWorkdir(workdir),
                lasFile.relativeToWorkdir(workdir)
            ),
            readIds.map!(to!string),
        ), workdir);
    }

    auto dbdump(Range)(in string dbFile, Range recordNumbers,
            in string[] dbdumpOptions, in string workdir)
            if (isForwardRange!Range && is(ElementType!Range : size_t))
    {
        return executePipe(chain(
            only("DBdump"),
            dbdumpOptions,
            only(dbFile.relativeToWorkdir(workdir)),
            recordNumbers.map!(to!string)
        ), workdir);
    }

    auto dbdump(
        in string dbFile,
        id_t firstRecord,
        id_t lastRecord,
        in string[] dbdumpOptions,
        in string workdir,
    )
    {
        return executePipe(chain(
            only("DBdump"),
            dbdumpOptions,
            only(
                dbFile.relativeToWorkdir(workdir),
                format!"%d-%d"(firstRecord, lastRecord),
            ),
        ), workdir);
    }

    string dbshow(in string dbFile, in string contigId, in string workdir)
    {
        return executeCommand(only("DBshow", dbFile.relativeToWorkdir(workdir), contigId), workdir);
    }

    string dbshow(in string dbFile, in string[] dbshowOptions, in string workdir)
    {
        return executeCommand(chain(only("DBshow"), dbshowOptions,
                only(dbFile.relativeToWorkdir(workdir))), workdir);
    }

    auto executePipe(Range)(Range command, in string workdir = null)
            if (isInputRange!Range && isSomeString!(ElementType!Range))
    {
        static struct LinesPipe
        {
            static immutable lineTerminator = "\n";

            private const string[] command;
            private const string workdir;
            private ProcessPipes process;
            private string currentLine;

            ~this()
            {
                if (!(process.pid is null))
                    releaseProcess();
            }

            void releaseProcess()
            {
                if (!process.stdout.isOpen)
                    return;

                process.stdout.close();

                version (Posix)
                {
                    import core.sys.posix.signal : SIGKILL;

                    process.pid.kill(SIGKILL);
                }
                else
                {
                    static assert(0, "Only intended for use on POSIX compliant OS.");
                }
                process.pid.wait();
            }

            private void ensureInitialized()
            {
                if (!(process.pid is null))
                    return;

                logJsonDiagnostic(
                    "action", "execute",
                    "type", "pipe",
                    "command", command.toJson,
                    "state", "pre",
                );
                process = pipeProcess(command, Redirect.stdout, null, Config.none, workdir);

                if (!empty)
                    popFront();
            }

            void popFront()
            {
                ensureInitialized();
                assert(!empty, "Attempting to popFront an empty LinesPipe");
                currentLine = process.stdout.readln();

                if (currentLine.empty)
                {
                    currentLine = null;
                    releaseProcess();
                }

                if (currentLine.endsWith(lineTerminator))
                    currentLine = currentLine[0 .. $ - lineTerminator.length];
            }

            @property string front()
            {
                ensureInitialized();
                assert(!empty, "Attempting to fetch the front of an empty LinesPipe");

                return currentLine;
            }

            @property bool empty()
            {
                ensureInitialized();

                if (!process.stdout.isOpen || process.stdout.eof)
                {
                    releaseProcess();

                    return true;
                }
                else
                {
                    return false;
                }
            }
        }

        auto sanitizedCommand = command.filter!"a != null".array;

        return new LinesPipe(sanitizedCommand, workdir);
    }

    unittest
    {
        import std.algorithm : equal;
        import std.range : take;

        auto cheers = executePipe(only("yes", "Cheers!"), ".");
        assert(cheers.take(5).equal([
            "Cheers!",
            "Cheers!",
            "Cheers!",
            "Cheers!",
            "Cheers!",
        ]));

        auto helloWorld = executePipe(only("echo", "Hello World!"), ".");
        assert(helloWorld.equal(["Hello World!"]));
    }

    string executeCommand(Range)(in Range command, in string workdir = null)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, execute;

        string output = command.executeWrapper!("command",
                sCmd => execute(sCmd, null, // env
                    Config.none, size_t.max, workdir));
        return output;
    }

    void executeShell(Range)(in Range command, in string workdir = null)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.algorithm : joiner;
        import std.process : Config, executeShell;

        string output = command.executeWrapper!("shell",
                sCmd => executeShell(sCmd.joiner(" ").array.to!string, null, // env
                    Config.none, size_t.max, workdir));
    }

    void executeScript(Range)(in Range command, in string workdir = null)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, executeShell;

        string output = command.executeWrapper!("script",
                sCmd => executeShell(sCmd.buildScriptLine, null, // env
                    Config.none, size_t.max, workdir));
    }

    string executeWrapper(string type, alias execCall, Range)(in Range command)
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

        immutable blockNumRegex = ctRegex!(`\.[1-9][0-9]*\.(dam|db)$`);

        return fileName.replaceFirst(blockNumRegex, `.$1`);
    }

    unittest
    {
        assert("foo_bar.1.dam".stripBlock == "foo_bar.dam");
        assert("foo_bar.1024.db".stripBlock == "foo_bar.db");
        assert("foo_bar.dam".stripBlock == "foo_bar.dam");
        assert("foo_bar.db".stripBlock == "foo_bar.db");
    }

    string relativeToWorkdir(in string fileName, in string workdir)
    {
        return relativePath(absolutePath(fileName), absolutePath(workdir));
    }
}

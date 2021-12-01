/**
    Defines bindings to and utilities for the Dazzler tool suite.

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
    iota,
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
    hasLength,
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

version (unittest)
{
    version = dumpLA;
}


/// File suffixes of hidden .db files.
private enum hiddenDbFileSuffixes = [".bps", ".idx"];

/// File suffixes of hidden .dam files.
private enum hiddenDamFileSuffixes = [".bps", ".hdr", ".idx"];

/// Constant holding the .db file extension.
enum dbFileExtension = ".db";

/// Constant holding the .dam file extension.
enum damFileExtension = ".dam";

/// The Dazzler tools require sequence of
/// at least `minSequenceLength` base pairs.
enum minSequenceLength = 14;

/// This trace point spacing enforces the use
/// of `ushort` for trace point encoding.
enum forceLargeTracePointType = 126;

/// Minimum allowed value for `-e` option of `daligner`/`damapper`
enum minAverageCorrelationRate = 0.7;

/// Minimum allowed value for `-n` option of `damapper`
enum minBestMatches = 0.7;

/// True if `T` is some array of strings.
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


/// Strip extension of `dbFile` if its is `dbFileExtension`
/// or `damFileExtension`. Otherwise return `dbFile` untouched.
string stripDbExtension(in string dbFile) pure nothrow @safe
{
    if (dbFile.endsWith(dbFileExtension))
        return dbFile[0 .. $ - dbFileExtension.length];
    else if (dbFile.endsWith(damFileExtension))
        return dbFile[0 .. $ - damFileExtension.length];
    else
        return dbFile[];
}


/// Signal error from the Dazzler bindings.
class DazzlerCommandException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}


/// Flag used to signal if DB write operations should append to an existing
/// DB or call `handleExistingDb` instead (default).
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


/// Returns the (trimmed) number of records in `dbFile`.
id_t numDbRecords(in string dbFile)
{
    auto dbdumpLines = dbdump(dbFile, []);
    scope (exit)
        // Clean up child process
        dbdumpLines.destroy();

    auto recordNumberLine = dbdumpLines.find!(l => l.startsWith("+ R")).front;
    id_t numRecords;
    recordNumberLine.formattedRead!"+ R %d"(numRecords);

    return numRecords;
}


/// Returns true iff `dbFile` is empty.
bool dbEmpty(in string dbFile)
{
    return numDbRecords(dbFile) == 0;
}


/// Remove database and hidden files.
@ExternalDependency("DBrm", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
void removeDB(in string dbFile)
{
    executeCommand(only("DBrm", dbFile));
}


/// Build DB file by using the given subset of reads in `inDbFile`.
/// Writes to `outputDB` if given; otherwise a name with the same extenion as
/// `inDbFile` under `options.tmpdir` is safely generated. The resulting DB
/// is implictly split using `options.dbsplitOptions`.
///
/// Note that the returned DB path always has the same extension as `inDbFile`.
/// If `outputDb` is given and does not have the expected extension it will be
/// append.
///
/// Returns: path to resulting DB
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

/// ditto
string dbSubset(Options, R)(in string outputDb, in string inDbFile, R readIds, in Options options, Append append = No.append)
        if (isOptionsList!(typeof(options.dbsplitOptions)))
{
    auto _outputDb = outputDb.extension == inDbFile.extension
        ? outputDb
        : outputDb ~ inDbFile.extension;

    buildSubsetDb(inDbFile, _outputDb, readIds, append);
    dbsplit(_outputDb, options.dbsplitOptions);

    return outputDb;
}


/// Merge given las files.
@ExternalDependency("LAmerge", "DALIGNER", "https://github.com/thegenemyers/DALIGNER")
void LAmerge(R)(in string mergedLas, R lasFiles)
{
    executeCommand(chain(only("LAmerge", mergedLas), lasFiles));
}


/// Compute local alignments of given DBs using `daligner`. Uses `dbB = dbA`
/// if `dbB` is omitted.
///
/// The resulting LAS file is placed in `options.tmpdir`. Its name can be
/// retrieved by calling `getLasFile(dbA, dbB, options.tmpdir)`. `daligner`
/// will not be executed if a file with that name already exists to avoid
/// redundant computations.
///
/// Returns: alignment data
/// See_also: `getAlignments`, `getLasFile`
AlignmentChain[] getLocalAlignments(Options)(in string dbA, in Options options)
        if (isOptionsList!(typeof(options.dalignerOptions)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    if (!lasFileGenerated(dbA, options.tmpdir))
        dalign(dbA, options.dalignerOptions, options.tmpdir);

    return getGeneratedAlignments(dbA, null, options);
}

/// ditto
AlignmentChain[] getLocalAlignments(Options)(in string dbA, in string dbB, in Options options)
        if (isOptionsList!(typeof(options.dalignerOptions)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    if (!lasFileGenerated(dbA, dbB, options.tmpdir))
        dalign(dbA, dbB, options.dalignerOptions, options.tmpdir);

    return getGeneratedAlignments(dbA, dbB, options);
}


/// Compute local alignments of given DBs using `damapper`.
///
/// The resulting LAS file is placed in `options.tmpdir`. Its name can be
/// retrieved by calling `getLasFile(dbA, dbB, options.tmpdir)`. `damapper`
/// will not be executed if a file with that name already exists to avoid
/// redundant computations.
///
/// Returns: alignment data
/// See_also: `getAlignments`, `getLasFile`
AlignmentChain[] getMappings(Options)(in string dbA, in string dbB, in Options options)
        if (isOptionsList!(typeof(options.damapperOptions)) &&
            isSomeString!(typeof(options.tmpdir)))
{
    if (!lasFileGenerated(dbA, dbB, options.tmpdir))
        damapper(dbA, dbB, options.damapperOptions, options.tmpdir);

    return getGeneratedAlignments(dbA, dbB, options);
}


private AlignmentChain[] getGeneratedAlignments(Options)(
    in string dbA,
    in string dbB,
    in Options options
) if (isSomeString!(typeof(options.tmpdir)))
{
    auto lasFile = getLasFile(dbA, dbB, options.tmpdir);

    return getAlignments(dbA, dbB, lasFile, options);
}


/// Control operation of `LocalAlignmentReader`.
enum AlignmentReaderFlag : uint
{
    /// Perform only minimal set of operations.
    none = 0,

    /// Read trace points and include them in the result.
    includeTracePoints = 1 << 0,

    /// Sort the alignments after reading.
    sort = 1 << 1,
}


/// Read all alignment chains from `lasFile`.
///
/// Params:
///     dbA         = A-read DB
///     dbB         = optional B-read DB; if omitted `dbB = dbA` is used
///     lasFile     = file with alignment
///     includeTracePoints = include trace points in output; this signature
///         implies `AlignmentReaderFlag.sort`. Skipping trace points is
///         faster and requires less memory.
///     flags       = specify additional operations (see `AlignmentReaderFlag`)
AlignmentChain[] getAlignments(
    in string dbA,
    in string lasFile,
    Flag!"includeTracePoints" includeTracePoints = No.includeTracePoints,
)
{
    return getAlignments(dbA, null, lasFile, includeTracePoints);
}

/// ditto
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

/// ditto
AlignmentChain[] getAlignments(
    in string dbA,
    in string lasFile,
    AlignmentReaderFlag flags,
)
{
    return getAlignments(dbA, null, lasFile, flags);
}

/// ditto
AlignmentChain[] getAlignments(
    in string dbA,
    in string dbB,
    in string lasFile,
    AlignmentReaderFlag flags,
)
{
    auto alignmentHeader = AlignmentHeader.inferFrom(lasFile);

    if (
        // there is nothing to do
        alignmentHeader.numAlignments == 0 ||
        // should be unreachable but is kept to prevent assertion failure in
        // LocalAlignmentReader
        alignmentHeader.numTracePoints == 0
    )
        return [];

    auto localAlignmentReader = new LocalAlignmentReader(
        lasFile,
        dbA,
        dbB,
        flags & AlignmentReaderFlag.includeTracePoints
            ? BufferMode.preallocated
            : BufferMode.skip,
        flags & AlignmentReaderFlag.includeTracePoints
            ? uninitializedArray!(TracePoint[])(alignmentHeader.numTracePoints)
            : [],
    );
    auto packer = localAlignmentReader.alignmentChainPacker(
        BufferMode.preallocated,
        uninitializedArray!(AlignmentChain.LocalAlignment[])(alignmentHeader.numLocalAlignments),
    );
    auto alignmentChainsBuffer = uninitializedArray!(AlignmentChain[])(
        alignmentHeader.numAlignments,
    );
    if (alignmentChainsBuffer.length > 0)
    {
        auto bufferRest = packer.copy(alignmentChainsBuffer);
        alignmentChainsBuffer.length -= bufferRest.length;
    }

    if (flags & AlignmentReaderFlag.sort)
        alignmentChainsBuffer.sort!("a < b", SwapStrategy.stable);

    return alignmentChainsBuffer;
}

unittest
{
    import dentist.util.tempfile : mkdtemp;
    import std.file : rmdirRecurse;
    import std.algorithm : equal;

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    scope (exit)
        rmdirRecurse(tmpDir);

    auto lasFile = buildPath(tmpDir, "test.las");
    dumpLA(lasFile, testLasDump);

    auto alignmentChains = getAlignments(
        null,
        null,
        lasFile,
        AlignmentReaderFlag.includeTracePoints,
    );
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    auto expectedResult = [
        AlignmentChain(
            0,
            Contig(1, 0),
            Contig(2, 0),
            AlignmentFlags(),
            [
                LocalAlignment(
                    Locus(3, 4),
                    Locus(5, 6),
                    0,
                    [TracePoint(0, 1)],
                ),
                LocalAlignment(
                    Locus(12, 13),
                    Locus(14, 15),
                    0,
                    [TracePoint(0, 1)],
                ),
            ],
            expectedTracePointDistance,
        ),
        AlignmentChain(
            1,
            Contig(19, 0),
            Contig(20, 0),
            AlignmentFlags(AlignmentFlag.complement, AlignmentFlag.alternateChain),
            [
                LocalAlignment(
                    Locus(21, 22),
                    Locus(23, 24),
                    0,
                    [TracePoint(0, 1)],
                ),
                LocalAlignment(
                    Locus(30, 31),
                    Locus(32, 33),
                    0,
                    [TracePoint(0, 1)],
                ),
            ],
            expectedTracePointDistance,
        ),
        AlignmentChain(
            2,
            Contig(37, 0),
            Contig(38, 0),
            AlignmentFlags(AlignmentFlag.unchained),
            [
                LocalAlignment(
                    Locus(39, 40),
                    Locus(41, 42),
                    0,
                    [TracePoint(0, 1)],
                ),
            ],
            expectedTracePointDistance,
        ),
        AlignmentChain(
            3,
            Contig(46, 0),
            Contig(47, 0),
            AlignmentFlags(AlignmentFlag.complement, AlignmentFlag.unchained),
            [
                LocalAlignment(
                    Locus(48, 49),
                    Locus(50, 51),
                    9,
                    [
                        TracePoint(2, 102),
                        TracePoint(3, 101),
                        TracePoint(4, 104),
                    ],
                ),
            ],
            expectedTracePointDistance,
        ),
        AlignmentChain(
            4,
            Contig(46, 0),
            Contig(47, 0),
            AlignmentFlags(AlignmentFlag.unchained),
            [
                LocalAlignment(
                    Locus(57, 58),
                    Locus(59, 60),
                    9,
                    [
                        TracePoint(3, 101),
                        TracePoint(4, 104),
                        TracePoint(2, 102),
                    ],
                ),
            ],
            expectedTracePointDistance,
        ),
        AlignmentChain(
            5,
            Contig(64, 0),
            Contig(65, 0),
            AlignmentFlags(AlignmentFlag.complement),
            [
                LocalAlignment(
                    Locus(66, 67),
                    Locus(68, 69),
                    12,
                    [
                        TracePoint(6, 105),
                        TracePoint(1, 101),
                        TracePoint(2, 100),
                        TracePoint(3, 97),
                    ],
                ),
                LocalAlignment(
                    Locus(75, 76),
                    Locus(77, 78),
                    2,
                    [
                        TracePoint(0, 2),
                        TracePoint(2, 102),
                    ],
                ),
            ],
            expectedTracePointDistance,
        ),
        AlignmentChain(
            6,
            Contig(1, 0),
            Contig(3197, 0),
            AlignmentFlags(AlignmentFlag.complement),
            [
                LocalAlignment(
                    Locus(0, 71),
                    Locus(12, 86),
                    3,
                    [
                        TracePoint(3, 74),
                    ],
                ),
                LocalAlignment(
                    Locus(0, 8300),
                    Locus(0, 318),
                    22,
                    [
                        TracePoint(6, 105),
                        TracePoint(9, 108),
                        TracePoint(7, 105),
                    ],
                ),
            ],
            expectedTracePointDistance,
        ),
    ];

    assert(alignmentChains == expectedResult);
}


/// Transforms a range `R` of `FlatLocalAlignment`s into a range of
/// `AlignmentChain`s.
///
/// See_also: `alignmentChainPacker` for more details.
struct AlignmentChainPacker(R)
{
    private alias LocalAlignment = AlignmentChain.LocalAlignment;

    private R alignments;
    private BufferMode bufferMode;
    private LocalAlignment[] localAlignmentBuffer;
    private size_t numBufferedLocalAlignments;
    private Appender!(LocalAlignment[]) localAlignmentAcc;
    private AlignmentChain currentChain;
    // used to check if the array is being reused
    private TracePoint* lastTracePointLocation;


    protected this(
        R alignments,
        BufferMode bufferMode,
        LocalAlignment[] localAlignmentBuffer = [],
    )
    {
        this.alignments = alignments;
        this.bufferMode = bufferMode;
        this.localAlignmentBuffer = localAlignmentBuffer;

        if (this.alignments.empty)
        {
            setEmpty();
        }
        else
        {
            popFront();
            // make sure ids start at zero
            --currentChain.id;
        }
    }


    /// Input range interface.
    @property bool empty() const pure nothrow @safe
    {
        return currentChain.id == id_t.max;
    }

    /// ditto
    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty AlignmentChainPacker");

        if (alignments.empty)
            return setEmpty();

        ++currentChain.id;
        currentChain.contigA = currentFLA.contigA.contig;
        currentChain.contigB = currentFLA.contigB.contig;
        currentChain.flags = currentFLA.flags & ~AlignmentFlags(AlignmentFlag.chainContinuation);
        if (bufferMode != BufferMode.skip && currentFLA.tracePoints.length > 0)
            currentChain.tracePointDistance = currentFLA.tracePointDistance;
        else
            currentChain.tracePointDistance = 0;
        bufferCurrentLocalAlignment();
        const chainStartFlags = currentFLA.flags;
        alignments.popFront();

        if (!chainStartFlags.unchained && !chainStartFlags.chainContinuation)
        {
            while (!alignments.empty && currentFLA.flags.chainContinuation)
            {
                version (assert)
                    assertMatchingAlignmentHead();
                bufferCurrentLocalAlignment();
                alignments.popFront();
            }
        }
        else
        {
            enforce!DazzlerCommandException(chainStartFlags.unchained, "chain is missing a start");
        }

        currentChain.localAlignments = finishLocalAlignmentBuffer();
    }


    /// ditto
    @property AlignmentChain front() pure nothrow @safe
    {
        assert(!empty, "Attempting to fetch the front of an empty AlignmentChainPacker");

        return currentChain;
    }


    private void setEmpty() pure nothrow @safe
    {
        currentChain.id = id_t.max;
    }


    version (assert) private void assertMatchingAlignmentHead() pure nothrow @safe
    {
        assert(currentChain.contigA == currentFLA.contigA.contig);
        assert(currentChain.contigB == currentFLA.contigB.contig);
        enum ignoredFlags = AlignmentFlags(
            AlignmentFlag.alternateChain,
            AlignmentFlag.chainContinuation,
        );
        assert((currentChain.flags & ~ignoredFlags) == (currentFLA.flags & ~ignoredFlags));
        assert(
            // we don't collect tracePoints -or-
            currentChain.tracePointDistance == 0 ||
            // the tracePointDistance must not change
            currentChain.tracePointDistance == currentFLA.tracePointDistance
        );
    }


    private void bufferCurrentLocalAlignment() pure nothrow @safe
    {
        final switch (bufferMode)
        {
            case BufferMode.overwrite: goto case;
            case BufferMode.preallocated:
                localAlignmentBuffer[numBufferedLocalAlignments++] = makeCurrentLocalAlignment();
                break;
            case BufferMode.dynamic:
                localAlignmentAcc ~= makeCurrentLocalAlignment();
                break;
            case BufferMode.skip:
                break;
        }
    }


    private LocalAlignment[] finishLocalAlignmentBuffer() pure nothrow @safe
    {
        typeof(return) localAlignments;

        final switch (bufferMode)
        {
            case BufferMode.overwrite:
                localAlignments = localAlignmentBuffer[0 .. numBufferedLocalAlignments];
                numBufferedLocalAlignments = 0;
                break;
            case BufferMode.preallocated:
                localAlignments = localAlignmentBuffer[0 .. numBufferedLocalAlignments];
                localAlignmentBuffer = localAlignmentBuffer[numBufferedLocalAlignments .. $];
                numBufferedLocalAlignments = 0;
                break;
            case BufferMode.dynamic:
                localAlignments = localAlignmentAcc.data;
                localAlignmentAcc = appender!(LocalAlignment[]);
                break;
            case BufferMode.skip:
                break;
        }

        return localAlignments;
    }


    private AlignmentChain.LocalAlignment makeCurrentLocalAlignment() pure nothrow @safe
    {
        if (currentChain.tracePointDistance > 0)
        {
            // we actually collect tracePoints
            assert(
                lastTracePointLocation is null || lastTracePointLocation != &currentFLA.tracePoints[0],
                "tracePoints buffer was reused; trace points will be invalid",
            );
            lastTracePointLocation = &currentFLA.tracePoints[0];
        }

        return AlignmentChain.LocalAlignment(
            currentFLA.contigA.locus,
            currentFLA.contigB.locus,
            currentFLA.tracePoints.map!"a.numDiffs".sum,
            currentFLA.tracePoints,
        );
    }


    private @property auto currentFLA() pure nothrow @safe
    {
        return alignments.front;
    }
}


/// Transforms a range `R` of `FlatLocalAlignment`s into a range of
/// `AlignmentChain`s.
///
/// Params:
///     localAlignments = input range of `FlatLocalAlignment`s
///     bufferMode      = selects how `LocalAlignment`s are buffered
///     localAlignmentBuffer = buffer for `LocalAlignment`s. If provided by
///         the caller it must be adequately sized with respect to
///         `bufferMode`; otherwise an adequate buffer is allocated.
/// See_also: `BufferMode`, `getFlatLocalAlignments`
auto alignmentChainPacker(R)(
    R localAlignments,
    BufferMode bufferMode = BufferMode.skip,
    AlignmentChain.LocalAlignment[] localAlignmentBuffer = [],
)
    if (isInputRange!R && is(ElementType!R == FlatLocalAlignment))
{
    return AlignmentChainPacker!R(localAlignments, bufferMode, localAlignmentBuffer);
}

// test empty input
unittest
{
    FlatLocalAlignment[] emptyAlignments;

    assert(emptyAlignments.alignmentChainPacker().empty);
}


/// Returns a tuple of contig IDs and first and last begin/end coords.
protected auto fingerprint(in ref AlignmentChain alignmentChain) pure nothrow
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


/// Lazily read individual local alignments from `lasFile`.
///
/// Returns: input range of `FlatLocalAlignment` with defined length.
/// Params:
///     dbA         = A-read DB used to insert contig lengths if present
///     dbB         = optional B-read DB; if omitted but `dbA` is given
///         `dbB = dbA` is used
///     lasFile     = file with alignment
///     bufferMode  = select buffer strategy for trace points
///     tracePointBuffer = buffer for `TracePoint`s. If provided by the caller
///         it must be adequately sized with respect to `bufferMode`;
///         otherwise an adequate buffer is allocated.
/// See_also: `BufferMode`, `FlatLocalAlignment`, `LocalAlignmentReader`
auto getFlatLocalAlignments(
    in string lasFile,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    return getFlatLocalAlignments(null, null, lasFile, bufferMode, tracePointBuffer);
}

/// ditto
auto getFlatLocalAlignments(
    in string dbA,
    in string lasFile,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    return getFlatLocalAlignments(dbA, null, lasFile, bufferMode, tracePointBuffer);
}

/// ditto
auto getFlatLocalAlignments(
    in string dbA,
    in string dbB,
    in string lasFile,
    BufferMode bufferMode = BufferMode.skip,
    TracePoint[] tracePointBuffer = [],
)
{
    if (
        bufferMode.among(BufferMode.overwrite, BufferMode.preallocated) &&
        tracePointBuffer.length == 0
    )
    {
        auto alignmentHeader = AlignmentHeader.inferFrom(lasFile);

        if (bufferMode == BufferMode.overwrite)
            tracePointBuffer = uninitializedArray!(typeof(tracePointBuffer))(
                alignmentHeader.maxTracePoints,
            );
        else if (bufferMode == BufferMode.preallocated)
            tracePointBuffer = uninitializedArray!(typeof(tracePointBuffer))(
                alignmentHeader.numTracePoints,
            );
    }

    return new LocalAlignmentReader(
        lasFile,
        dbA,
        dbB,
        bufferMode,
        tracePointBuffer,
    );
}

version (unittest)
{
    enum expectedTracePointDistance = 100;
    enum testLasDump = [
        "+ P 11",
        "% P 2",
        "+ T 42",
        "% T 12",
        "@ T 8",
        "X 100",
        "P 1 2 n >",
        "C 3 4 5 6",
        "T 1",
        "   0 1",
        "P 1 2 n -",
        "C 12 13 14 15",
        "T 1",
        "   0 1",
        "P 19 20 c +",
        "C 21 22 23 24",
        "T 1",
        "   0 1",
        "P 19 20 c -",
        "C 30 31 32 33",
        "T 1",
        "   0 1",
        "P 37 38 n .",
        "C 39 40 41 42",
        "T 1",
        "   0 1",
        "P 46 47 c .",
        "C 48 49 50 51",
        "T 3",
        "   2 102",
        "   3 101",
        "   4 104",
        "P 46 47 n .",
        "C 57 58 59 60",
        "T 3",
        "   3 101",
        "   4 104",
        "   2 102",
        "P 64 65 c >",
        "C 66 67 68 69",
        "T 4",
        "   6 105",
        "   1 101",
        "   2 100",
        "   3  97",
        "P 64 65 c -",
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
}

unittest
{
    import dentist.util.tempfile : mkdtemp;
    import std.file : rmdirRecurse;
    import std.algorithm : equal;
    import std.range : tee;

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    scope (exit)
        rmdirRecurse(tmpDir);

    auto lasFile = buildPath(tmpDir, "test.las");
    dumpLA(lasFile, testLasDump);

    auto flatLocalAlignments = getFlatLocalAlignments(lasFile, BufferMode.preallocated).array;
    alias FlatLocus = FlatLocalAlignment.FlatLocus;
    auto expectedResult = [
        FlatLocalAlignment(
            0,
            FlatLocus(1, 0, 3, 4),
            FlatLocus(2, 0, 5, 6),
            AlignmentFlags(),
            expectedTracePointDistance,
            [
                TracePoint(0, 1),
            ],
        ),
        FlatLocalAlignment(
            1,
            FlatLocus(1, 0, 12, 13),
            FlatLocus(2, 0, 14, 15),
            AlignmentFlags(AlignmentFlag.chainContinuation),
            expectedTracePointDistance,
            [
                TracePoint(0, 1),
            ],
        ),
        FlatLocalAlignment(
            2,
            FlatLocus(19, 0, 21, 22),
            FlatLocus(20, 0, 23, 24),
            AlignmentFlags(AlignmentFlag.complement, AlignmentFlag.alternateChain),
            expectedTracePointDistance,
            [
                TracePoint(0, 1),
            ],
        ),
        FlatLocalAlignment(
            3,
            FlatLocus(19, 0, 30, 31),
            FlatLocus(20, 0, 32, 33),
            AlignmentFlags(AlignmentFlag.complement, AlignmentFlag.chainContinuation),
            expectedTracePointDistance,
            [
                TracePoint(0, 1),
            ],
        ),
        FlatLocalAlignment(
            4,
            FlatLocus(37, 0, 39, 40),
            FlatLocus(38, 0, 41, 42),
            AlignmentFlags(AlignmentFlag.unchained),
            expectedTracePointDistance,
            [
                TracePoint(0, 1),
            ],
        ),
        FlatLocalAlignment(
            5,
            FlatLocus(46, 0, 48, 49),
            FlatLocus(47, 0, 50, 51),
            AlignmentFlags(AlignmentFlag.complement, AlignmentFlag.unchained),
            expectedTracePointDistance,
            [
                TracePoint(2, 102),
                TracePoint(3, 101),
                TracePoint(4, 104),
            ],
        ),
        FlatLocalAlignment(
            6,
            FlatLocus(46, 0, 57, 58),
            FlatLocus(47, 0, 59, 60),
            AlignmentFlags(AlignmentFlag.unchained),
            expectedTracePointDistance,
            [
                TracePoint(3, 101),
                TracePoint(4, 104),
                TracePoint(2, 102),
            ],
        ),
        FlatLocalAlignment(
            7,
            FlatLocus(64, 0, 66, 67),
            FlatLocus(65, 0, 68, 69),
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
            FlatLocus(64, 0, 75, 76),
            FlatLocus(65, 0, 77, 78),
            AlignmentFlags(AlignmentFlag.complement, AlignmentFlag.chainContinuation),
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
            AlignmentFlags(AlignmentFlag.complement, AlignmentFlag.chainContinuation),
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


/// Meta information about a set of alignments, e.g. a LAS file.
///
/// Usually this is created by reading a LAS file as efficiently as possible
/// to improve performance of subsequent operations on the LAS file, e.g.
/// by pre-allocating memory.
struct AlignmentHeader
{
    /// Sum of number of alignment chains and number of unchained alignments
    size_t numAlignments;
    /// Number of alignment chains (i.e. the unchained flag is unset)
    size_t numAlignmentChains;
    /// Total number of local alignments disregarding chaining
    size_t numLocalAlignments;
    /// Maximum number of local alignments per chain
    size_t maxLocalAlignments;
    /// Maximum total number of local alignments per contig
    size_t maxLocalAlignmentsPerContig;
    /// Total number of trace points
    size_t numTracePoints;
    /// Maximum number of trace points per local alignment
    size_t maxTracePoints;
    // TODO size_t maxTracePointsPerContig;
    /// Trace point distance
    size_t tracePointDistance;


    /// Infer header data from a range of alignments with a single pass
    /// across the data.
    ///
    /// The range signatures accepts currently only a range of
    /// `AlignmentChain`s.
    ///
    /// The `lasFile` signature walks through `lasFile` using a
    /// `LocalAlignmentReader` with `BufferMode.skip`.
    static AlignmentHeader inferFrom(R)(R alignmentChains) if (isInputRange!R)
    {
        AlignmentHeader headerData;

        headerData.tracePointDistance = inferTracePointDistanceFrom(alignmentChains);

        id_t lastContig;
        id_t numLocalAlignmentsSinceLastContig;
        foreach (alignmentChain; alignmentChains)
        {
            if (lastContig != alignmentChain.contigA.id)
            {
                headerData.maxLocalAlignmentsPerContig = max(
                    headerData.maxLocalAlignmentsPerContig,
                    numLocalAlignmentsSinceLastContig,
                );
                numLocalAlignmentsSinceLastContig = 0;
            }

            ++headerData.numAlignments;
            if (!alignmentChain.flags.unchained)
                ++headerData.numAlignmentChains;
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

            numLocalAlignmentsSinceLastContig += alignmentChain.localAlignments.length;
            lastContig = alignmentChain.contigA.id;
        }

        return headerData;
    }


    /// ditto
    static AlignmentHeader inferFrom(string lasFile)
    {
        AlignmentHeader headerData;

        auto lasScanner = new LocalAlignmentReader(lasFile, BufferMode.skip);

        size_t currentChainLength;
        id_t lastContig;
        id_t numLocalAlignmentsSinceLastContig;
        headerData.tracePointDistance = lasScanner.tracePointDistance;
        foreach (flatLocalAlignment; lasScanner)
        {
            if (lastContig != flatLocalAlignment.contigA.id)
            {
                headerData.maxLocalAlignmentsPerContig = max(
                    headerData.maxLocalAlignmentsPerContig,
                    numLocalAlignmentsSinceLastContig,
                );
                numLocalAlignmentsSinceLastContig = 0;
            }

            if (flatLocalAlignment.flags.unchained)
            {
                currentChainLength = 1;
                ++headerData.numAlignments;
            }
            else if (flatLocalAlignment.flags.chainContinuation)
            {
                ++currentChainLength;
            }
            else
            {
                currentChainLength = 1;
                ++headerData.numAlignmentChains;
                ++headerData.numAlignments;
            }

            ++headerData.numLocalAlignments;
            ++numLocalAlignmentsSinceLastContig;

            headerData.numTracePoints += lasScanner.currentNumTracePoints;
            headerData.maxTracePoints = max(
                headerData.maxTracePoints,
                lasScanner.currentNumTracePoints,
            );
            headerData.maxLocalAlignments = max(
                headerData.maxLocalAlignments,
                currentChainLength,
            );

            lastContig = flatLocalAlignment.contigA.id;
        }

        return headerData;
    }


    /// Infer trace points spacing from `alignments`.
    ///
    /// This do NOT advance `alignments`.
    ///
    /// Returns: `tracePointDistance` of the first element of `alignments` or
    ///     `100` if `alignments` is empty.
    static size_t inferTracePointDistanceFrom(R)(R alignments)
        if (isInputRange!R && (
            is(const(ElementType!R) == const(AlignmentChain)) ||
            is(const(ElementType!R) == const(FlatLocalAlignment))
        ))
    {
        if (alignments.empty)
            return 100;

        return alignments.front.tracePointDistance;
    }
}


version (unittest)
{
    private AlignmentChain[] getTestAlignmentChains(trace_point_t tracePointSpacing)
    {
        alias LocalAlignment = AlignmentChain.LocalAlignment;
        alias TracePoint = LocalAlignment.TracePoint;
        enum complement = AlignmentFlag.complement;
        enum alternateChain = AlignmentFlag.alternateChain;

        return [
            AlignmentChain(
                0,
                Contig(1, 15),
                Contig(2, 16),
                AlignmentFlags(),
                [
                    LocalAlignment(
                        Locus(3, 4),
                        Locus(5, 6),
                        7,
                        [TracePoint(7, 1)],
                    ),
                    LocalAlignment(
                        Locus(12, 13),
                        Locus(14, 15),
                        16,
                        [TracePoint(16, 1)],
                    ),
                ],
                tracePointSpacing,
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
                        [TracePoint(25, 1)],
                    ),
                    LocalAlignment(
                        Locus(30, 31),
                        Locus(32, 33),
                        0,
                        [TracePoint(0, 1)],
                    ),
                ],
                tracePointSpacing,
            ),
        ];
    }

    private string[] getTestFastaRecords()
    {
        auto fakeSeq = (coord_t length) => iota(length).map!"'a'".array.to!string;
        auto fakeRead = (id_t id, coord_t length) => format!">faked/%d/0_%d RQ=0.85\n%s"(
            id,
            id,
            fakeSeq(length),
        );

        return [
            fakeRead(1, 15),
            fakeRead(2, 16),
            fakeRead(3, 42),
            fakeRead(4, 42),
            fakeRead(5, 42),
            fakeRead(6, 42),
            fakeRead(7, 42),
            fakeRead(8, 42),
            fakeRead(9, 42),
            fakeRead(10, 42),
            fakeRead(11, 42),
            fakeRead(12, 42),
            fakeRead(13, 42),
            fakeRead(14, 42),
            fakeRead(15, 42),
            fakeRead(16, 42),
            fakeRead(17, 42),
            fakeRead(18, 42),
            fakeRead(19, 31),
            fakeRead(20, 33),
            fakeRead(21, 42),
            fakeRead(22, 42),
            fakeRead(23, 42),
            fakeRead(24, 42),
            fakeRead(25, 42),
            fakeRead(26, 42),
        ];
    }
}


/// Controls how to manage the buffer when reading data into memory.
enum BufferMode : ubyte
{
    /// Keep a single buffer and keep overwriting it with every new record.
    overwrite,
    /// Allocate a new buffer for every record. Use a `std.array.Appender`
    /// if the number of records is unknown.
    dynamic,
    /// Write all records to a continuous stretch of memory pre-allocated by
    /// the caller. The memory may be uninitialized since it will always be
    /// written to before any read occurs.
    preallocated,
    /// Do not use memory at all instead just skip over the records.
    skip,
}


/// Read local alignments from a LAS file.
///
/// See_also: `getFlatLocalAlignments` for more details
class LocalAlignmentReader
{
    private File las;
    private coord_t[] aLengths;
    private coord_t[] bLengths;
    private BufferMode bufferMode;
    private size_t _numLocalAlignments;
    private trace_point_t _tracePointDistance;
    private bool isLargeTraceType;

    private size_t numLocalAlignmentsLeft;
    private FlatLocalAlignment currentLA;
    private DazzlerOverlap overlapHead;
    private TracePoint[] tracePointBuffer;
    private TracePoint[] fullTracePointBuffer;


    /// Total number of local alignments in the LAS file. This is not affected
    /// by range operation, e.g. `popFront`.
    ///
    /// See_also: `length`
    @property size_t numLocalAlignments() const pure nothrow @safe @nogc
    {
        return _numLocalAlignments;
    }

    private @property void numLocalAlignments(size_t newValue) pure nothrow @safe @nogc
    {
        _numLocalAlignments = newValue;
    }


    /// Trace point spacing of the LAS file.
    @property trace_point_t tracePointDistance() const pure nothrow @safe @nogc
    {
        return _tracePointDistance;
    }

    private @property void tracePointDistance(trace_point_t newValue) pure nothrow @safe @nogc
    {
        _tracePointDistance = newValue;
    }


    protected this(const string lasFile, BufferMode bufferMode, TracePoint[] tracePointBuffer = [])
    {
        this(
            lasFile,
            cast(string) null,
            cast(string) null,
            bufferMode,
            tracePointBuffer
        );
    }


    protected this(
        const string lasFile,
        string dbA,
        string dbB,
        BufferMode bufferMode,
        TracePoint[] tracePointBuffer = [],
    )
    {
        this(
            lasFile,
            contigLengths(dbA),
            dbB !is null && dbA != dbB
                ? contigLengths(dbB)
                : contigLengths(dbA),
            bufferMode,
            tracePointBuffer
        );
    }


    private this(
        const string lasFile,
        coord_t[] aLengths,
        coord_t[] bLengths,
        BufferMode bufferMode,
        TracePoint[] tracePointBuffer,
    )
    {
        this.las = File(lasFile, "rb");
        this.bufferMode = bufferMode;
        this.tracePointBuffer = tracePointBuffer;
        if (bufferMode.among(BufferMode.dynamic, BufferMode.skip))
            assert(
                tracePointBuffer.length == 0,
                "preallocated buffer was supplied but buffer mode does not one: "~
                "please remove the buffer or adjust the mode according to your needs."
            );
        else
            assert(
                tracePointBuffer.length > 0,
                "buffer mode requires a preallocated buffer but none was supplied: "~
                "please supply a preallocated buffer or adjust the mode according to your needs."
            );
        this.fullTracePointBuffer = tracePointBuffer;
        this.aLengths = aLengths;
        this.bLengths = bLengths;

        initialize();
    }


    private coord_t[] contigLengths(string db)
    {
        if (db is null)
            return [];

        return getDbRecords(db, [DBdumpOptions.originalHeader])
            .map!"a.location.length"
            .array;
    }


    private void initialize()
    {
        readHeader();
        this.numLocalAlignmentsLeft = this.numLocalAlignments;

        if (!empty)
        {
            ++numLocalAlignmentsLeft;
            popFront();
        }
    }


    /// Reset the range. This can be used to walk over a LAS file multiple
    /// times without reopening the `std.stdio.File` or recreating the
    /// buffers.
    ///
    /// This requires a seek-able file as it seeks to the begin of the
    /// underlying LAS file and adjusts the range state accordingly.
    void reset()
    {
        las.rewind();
        tracePointBuffer = fullTracePointBuffer;
        initialize();
    }


    /// Range interface.
    @property bool empty() const pure nothrow @safe
    {
        return numLocalAlignmentsLeft == 0;
    }


    /// ditto
    @property size_t length() const pure nothrow @safe
    {
        return numLocalAlignmentsLeft;
    }


    /// ditto
    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty LocalAlignmentReader");

        --numLocalAlignmentsLeft;
        if (!empty)
            readLocalAlignment();
    }


    /// ditto
    @property FlatLocalAlignment front() pure nothrow @safe
    {
        assert(!empty, "Attempting to fetch the front of an empty LocalAlignmentReader");

        return currentLA;
    }


    /// Returns whether reading or skipping trace points.
    @property bool skipTracePoints() const pure nothrow @safe
    {
        return bufferMode == BufferMode.skip;
    }


    /// Returns the number of trace points in the current local alignment
    /// even if `skipTracePoints` is true.
    ///
    /// This is required because the number of trace points is unavailable in
    /// the `front` element if `skipTracePoints` is true.
    @property size_t currentNumTracePoints() const pure nothrow @safe
    {
        return overlapHead.path.tlen / 2;
    }

protected:

    void readHeader()
    {
        readNumLocalAlignments();
        readTracePointDistance();
    }


    void readNumLocalAlignments()
    {
        long numLocalAlignments;
        auto nlaBuffer = las.rawRead((&numLocalAlignments)[0 .. 1]);
        unexpectedEOF!"numLocalAlignments"(nlaBuffer, (&numLocalAlignments)[0 .. 1]);
        this.numLocalAlignments = numLocalAlignments.to!size_t;
    }


    void readTracePointDistance()
    {
        int tracePointDistance;
        auto tpdBuffer = las.rawRead((&tracePointDistance)[0 .. 1]);
        unexpectedEOF!"tracePointDistance"(tpdBuffer, (&tracePointDistance)[0 .. 1]);
        this.tracePointDistance = tracePointDistance.to!trace_point_t;
        this.isLargeTraceType = DazzlerOverlap.isLargeTraceType(tracePointDistance);
    }


    void readLocalAlignment()
    {
        readOverlapHead();
        fillInOverlapHead();
        if (aLengths.length > 0)
            fillInContigLengths();
        auto traceLength = getTraceVectorLength();

        if (skipTracePoints)
            skipTraceVector(traceLength);
        else
            readTraceVector(traceLength);

        if (!skipTracePoints)
        {
            if (!isLargeTraceType)
                rewriteTracePointBufferToLargeTracePointType();

            auto numTracePoints = overlapHead.path.tlen / 2;
            currentLA.tracePoints = tracePointBuffer[0 .. numTracePoints];
            if (bufferMode == BufferMode.preallocated)
                tracePointBuffer = tracePointBuffer[numTracePoints .. $];
        }
    }


    void readOverlapHead()
    {
        auto overlapBytes = (cast(void*) &overlapHead)[
            typeof(DazzlerOverlap.path.trace).sizeof ..
            DazzlerOverlap.sizeof
        ];
        auto buffer = las.rawRead(overlapBytes);
        unexpectedEOF!"overlapHead"(buffer, overlapBytes);
    }


    void fillInOverlapHead() pure nothrow @safe
    {
        currentLA.id = numLocalAlignments - numLocalAlignmentsLeft;
        currentLA.contigA.id = overlapHead.aread + 1;
        currentLA.contigA.begin = overlapHead.path.abpos;
        currentLA.contigA.end = overlapHead.path.aepos;
        currentLA.contigB.id = overlapHead.bread + 1;
        currentLA.contigB.begin = overlapHead.path.bbpos;
        currentLA.contigB.end = overlapHead.path.bepos;

        currentLA.flags = AlignmentFlags();
        if (overlapHead.flags & DazzlerOverlap.Flag.disabled)
            currentLA.flags |= AlignmentFlag.disabled;
        if (overlapHead.flags & DazzlerOverlap.Flag.complement)
            currentLA.flags |= AlignmentFlag.complement;
        if (
            (overlapHead.flags & DazzlerOverlap.Flag.chainStart) &&
            !(overlapHead.flags & DazzlerOverlap.Flag.bestChain)
        )
            currentLA.flags |= AlignmentFlag.alternateChain;
        if (overlapHead.flags & DazzlerOverlap.Flag.chainContinuation)
            currentLA.flags |= AlignmentFlag.chainContinuation;
        if (!(overlapHead.flags & (
            DazzlerOverlap.Flag.chainStart |
            DazzlerOverlap.Flag.bestChain |
            DazzlerOverlap.Flag.chainContinuation
        )))
            currentLA.flags |= AlignmentFlag.unchained;

        currentLA.tracePointDistance = tracePointDistance;
    }


    void fillInContigLengths() pure @safe
    {
        if (aLengths.length > 0)
        {
            enforce!DazzlerCommandException(
                currentLA.contigA.id - 1 < aLengths.length,
                format!"contigA.id out of bounds: %d >= %d"(currentLA.contigA.id - 1, aLengths.length),
            );
            currentLA.contigA.length = aLengths[currentLA.contigA.id - 1];
        }

        if (bLengths.length)
        {
            enforce!DazzlerCommandException(
                currentLA.contigB.id - 1 < bLengths.length,
                format!"contigB.id out of bounds: %d >= %d"(currentLA.contigB.id - 1, bLengths.length),
            );
            currentLA.contigB.length = bLengths[currentLA.contigB.id - 1];
        }
    }


    size_t getTraceVectorLength() const pure @safe
    {
        auto traceVectorLength = overlapHead.path.tlen * (isLargeTraceType
            ? DazzlerOverlap.largeTraceType.sizeof
            : DazzlerOverlap.smallTraceType.sizeof
        );
        enforce!DazzlerCommandException(
            traceVectorLength % 2 == 0,
            "illegal value for tlen: must be multiple of 2",
        );

        return traceVectorLength;
    }


    void readTraceVector(size_t traceLength)
    {
        if (bufferMode == BufferMode.dynamic)
            tracePointBuffer = uninitializedArray!(TracePoint[])(overlapHead.path.tlen / 2);
        auto rawBuffer = getRawTracePointBuffer(traceLength);
        auto tpBuffer = las.rawRead(rawBuffer);
        unexpectedEOF!"tracePoints"(tpBuffer, rawBuffer);
    }


    void skipTraceVector(size_t traceLength)
    {
        import core.stdc.stdio : SEEK_CUR;

        las.seek(traceLength, SEEK_CUR);
    }


    void rewriteTracePointBufferToLargeTracePointType()
    {
        static assert(is(DazzlerOverlap.largeTraceType == typeof(tracePointBuffer[0].numDiffs)));

        auto traceVectorLength = overlapHead.path.tlen;
        // make sure enough memory is allocated
        auto smallTraceVector = getRawTracePointBuffer(traceVectorLength);
        static assert(is(DazzlerOverlap.smallTraceType == typeof(smallTraceVector[0])));
        // select appropriate portion of tracePointBuffer
        auto largeTracePoints = tracePointBuffer.ptr[0 .. traceVectorLength / 2];
        // convert numbers in reverse order as to avoid overwriting
        foreach_reverse (i, ref largeTracePoint; largeTracePoints)
        {
            largeTracePoint = TracePoint(
                smallTraceVector[2 * i],
                smallTraceVector[2 * i + 1],
            );
        }
    }


    @property ubyte[] getRawTracePointBuffer(size_t traceLength) pure nothrow
    {
        // prepare bytes buffer
        auto bufferBytes = tracePointBuffer.length * TracePoint.sizeof;
        auto rawBuffer = (cast(ubyte*) tracePointBuffer.ptr)[0 .. bufferBytes];
        // reduce size (triggers bounds check)
        rawBuffer = rawBuffer[0 .. traceLength];

        return rawBuffer;
    }


    void unexpectedEOF(string what, T)(const T[] gotBuffer, const T[] expBuffer)
    {
        enum msgFormat = "error reading LAS file `%s`: unexpected end of file; expected " ~ what;

        enforce!DazzlerCommandException(
            gotBuffer.length == expBuffer.length,
            format!msgFormat(las.name),
        );
    }
}

static assert(isInputRange!LocalAlignmentReader);
static assert(hasLength!LocalAlignmentReader);
static assert(is(ElementType!LocalAlignmentReader == FlatLocalAlignment));

unittest
{
    import dentist.util.tempfile : mkdtemp;
    import std.file : rmdirRecurse;
    import std.algorithm : equal;

    alias LocalAlignment = AlignmentChain.LocalAlignment;
    enum complement = AlignmentFlag.complement;

    static foreach (tracePointDistance; [100, 200])
    {{
        auto tmpDir = mkdtemp("./.unittest-XXXXXX");
        scope (exit)
            rmdirRecurse(tmpDir);

        auto dbFile = buildPath(tmpDir, "test.db");
        buildDbFile(dbFile, getTestFastaRecords());

        auto lasFile = buildPath(tmpDir, "test.las");
        auto alignmentChains = getTestAlignmentChains(100);
        auto computedFlatLocalAlignments = alignmentChains
            .map!"a.toFlatLocalAlignments()"
            .joiner
            .enumerate
            .map!((enumFla) {
                auto fla = enumFla.value;
                fla.id = enumFla.index.to!id_t;
                fla.contigA.length = 0;
                fla.contigB.length = 0;

                return fla;
            });

        lasFile.writeAlignments(alignmentChains);
        auto recoveredFlatLocalAlignments = new LocalAlignmentReader(
            lasFile,
            BufferMode.dynamic,
        );

        assert(equal(computedFlatLocalAlignments, recoveredFlatLocalAlignments));
    }}
}


/// Write alignments to `lasFile`. This method takes care of the differences
/// between alignment flags in DENTIST and Dazzler code.
///
/// `lasFile` must point to a seek-able file because the LAS header is updated
/// after writing all data to avoid traversing the data twice.
void writeAlignments(R)(const string lasFile, R flatLocalAlignments)
    if (isInputRange!R && is(const(ElementType!R) == const(FlatLocalAlignment)))
{
    auto las = File(lasFile, "wb");

    long numLocalAlignments = 0; // will be overwritten at the end
    auto tracePointDistance = AlignmentHeader
        .inferTracePointDistanceFrom(flatLocalAlignments)
        .to!int;
    las.rawWrite([numLocalAlignments]);
    las.rawWrite([tracePointDistance]);

    foreach (flatLocalAlignment; flatLocalAlignments)
    {
        las.writeFlatLocalAlignment(flatLocalAlignment, tracePointDistance);
        ++numLocalAlignments;
    }

    las.rewind();
    las.rawWrite([numLocalAlignments]);

    las.close();
}

/// ditto
void writeAlignments(R)(const string lasFile, R alignmentChains)
    if (isInputRange!R && is(const(ElementType!R) == const(AlignmentChain)))
{
    auto las = File(lasFile, "wb");

    long numLocalAlignments = 0; // will be overwritten at the end
    auto tracePointDistance = AlignmentHeader
        .inferTracePointDistanceFrom(alignmentChains)
        .to!int;
    las.rawWrite([numLocalAlignments]);
    las.rawWrite([tracePointDistance]);

    foreach (alignmentChain; alignmentChains)
    {
        las.writeAlignmentChain(alignmentChain, tracePointDistance);
        numLocalAlignments += alignmentChain.localAlignments.length;
    }

    las.rewind();
    las.rawWrite([numLocalAlignments]);

    las.close();
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

    auto dbFile = buildPath(tmpDir, "test.db");
    buildDbFile(dbFile, getTestFastaRecords());

    auto lasFile = buildPath(tmpDir, "test.las");
    const alignmentChains = getTestAlignmentChains(100);

    lasFile.writeAlignments(alignmentChains);
    auto recoveredAlignmentChains = getAlignments(dbFile, lasFile, Yes.includeTracePoints);

    assert(alignmentChains == recoveredAlignmentChains);
}


/// Struct `Overlap` from `dalign.h` used for LAS file I/O.
private struct DazzlerOverlap
{
    /// Flags as defined in `dalign.h`
    static enum Flag : uint
    {
        complement = 0x1,
        chainStart = 0x4,
        chainContinuation = 0x8,
        bestChain = 0x10,
        disabled = 0x20,
    }

    /// Struct `Path` from `dalign.h`
    static struct Path
    {
        void* trace;
        int tlen;
        int diffs;
        int abpos;
        int bbpos;
        int aepos;
        int bepos;
    }

    Path path;
    uint flags;
    int aread;
    int bread;


    /// Utility function to determine the type of trace points.
    static bool isLargeTraceType(const int tracePointDistance) pure nothrow @safe
    {
        // Macro definition from `dalign.h`
        enum TRACE_XOVR = 125u;

        return tracePointDistance > TRACE_XOVR;
    }


    /// Aliases for trace point types.
    alias largeTraceType = ushort;
    /// ditto
    alias smallTraceType = ubyte;
}


/// Write a single `AlignmentChain` to `las` taking care of the appropriate
/// flags for chaining.
private auto writeAlignmentChain(
    File las,
    const AlignmentChain alignmentChain,
    const int tracePointDistance,
)
{
    if (alignmentChain.localAlignments.length == 0)
        return;

    DazzlerOverlap dazzlerOverlap;

    // set read IDs
    dazzlerOverlap.aread = alignmentChain.contigA.id - 1;
    dazzlerOverlap.bread = alignmentChain.contigB.id - 1;
    // set initial flags
    if (alignmentChain.flags.disabled)
        dazzlerOverlap.flags |= DazzlerOverlap.Flag.disabled;
    if (alignmentChain.flags.complement)
        dazzlerOverlap.flags |= DazzlerOverlap.Flag.complement;

    // store initial flags
    const initialFlags = dazzlerOverlap.flags;

    foreach (i, localAlignment; alignmentChain.localAlignments)
    {
        // set flags for this local alignment
        dazzlerOverlap.flags = initialFlags;
        if (i == 0)
        {
            dazzlerOverlap.flags |= DazzlerOverlap.Flag.chainStart;
            if (!alignmentChain.flags.alternateChain)
                dazzlerOverlap.flags |= DazzlerOverlap.Flag.bestChain;
        }
        else
        {
            dazzlerOverlap.flags |= DazzlerOverlap.Flag.chainContinuation;
        }

        // set coordinates
        dazzlerOverlap.path.abpos = localAlignment.contigA.begin;
        dazzlerOverlap.path.aepos = localAlignment.contigA.end;
        dazzlerOverlap.path.bbpos = localAlignment.contigB.begin;
        dazzlerOverlap.path.bepos = localAlignment.contigB.end;

        writeDazzlerOverlap(las, dazzlerOverlap, localAlignment.tracePoints, tracePointDistance);
    }
}


/// Write a single `AlignmentChain` to `las` taking care of translating flags
/// from DENTIST to Dazzler.
private auto writeFlatLocalAlignment(
    File las,
    const FlatLocalAlignment flatLocalAlignment,
    const int tracePointDistance,
)
{
    DazzlerOverlap dazzlerOverlap;

    // set read IDs
    dazzlerOverlap.aread = flatLocalAlignment.contigA.id - 1;
    dazzlerOverlap.bread = flatLocalAlignment.contigB.id - 1;
    // set initial flags
    if (flatLocalAlignment.flags.disabled)
        dazzlerOverlap.flags |= DazzlerOverlap.Flag.disabled;
    if (flatLocalAlignment.flags.complement)
        dazzlerOverlap.flags |= DazzlerOverlap.Flag.complement;
    if (flatLocalAlignment.flags.chainContinuation)
    {
        dazzlerOverlap.flags |= DazzlerOverlap.Flag.chainContinuation;
    }
    else if (!flatLocalAlignment.flags.unchained)
    {
        dazzlerOverlap.flags |= DazzlerOverlap.Flag.chainStart;

        if (!flatLocalAlignment.flags.alternateChain)
            dazzlerOverlap.flags |= DazzlerOverlap.Flag.bestChain;
    }

    // set coordinates
    dazzlerOverlap.path.abpos = flatLocalAlignment.contigA.begin;
    dazzlerOverlap.path.aepos = flatLocalAlignment.contigA.end;
    dazzlerOverlap.path.bbpos = flatLocalAlignment.contigB.begin;
    dazzlerOverlap.path.bepos = flatLocalAlignment.contigB.end;

    writeDazzlerOverlap(las, dazzlerOverlap, flatLocalAlignment.tracePoints, tracePointDistance);
}


/// Write `dazzlerOverlap` to `las` compressing the trace vector if
/// `tracePointDistance` is sufficiently small.
///
/// See_also: `DazzlerOverlap.isLargeTraceType`
private void writeDazzlerOverlap(
    File las,
    ref DazzlerOverlap dazzlerOverlap,
    const TracePoint[] tracePoints,
    const int tracePointDistance,
)
{
    // select appropriate type for trace
    auto isLargeTraceType = DazzlerOverlap.isLargeTraceType(tracePointDistance);

    // set trace vector length
    dazzlerOverlap.path.tlen = to!int(2*tracePoints.length);
    // set diffs
    dazzlerOverlap.path.diffs = tracePoints.map!"a.numDiffs".sum;

    // write overlap "header"
    auto overlapBytes = (cast(void*) &dazzlerOverlap)[
        typeof(dazzlerOverlap.path.trace).sizeof ..
        DazzlerOverlap.sizeof
    ];
    las.rawWrite(overlapBytes);

    // write trace vector
    if (tracePoints.length > 0)
    {
        if (isLargeTraceType)
            las.rawWrite(tracePoints);
        else
            // rewrite trace to use `ubyte`s
            las.rawWrite(
                tracePoints
                    .map!(tp => [
                        tp.numDiffs.to!(DazzlerOverlap.smallTraceType),
                        tp.numBasePairs.to!(DazzlerOverlap.smallTraceType),
                    ])
                    .joiner
                    .takeExactly(2*tracePoints.length)
                    .array
            );
    }
}


/// Returns the trace point distance in `lasFile`.
trace_point_t getTracePointDistance(in string lasFile)
{
    return readLasHeader(lasFile).tracePointSpacing;
}


/// $(RED [Experimental]) Reconstruct the alignment matrix for `ac`
/// force-filling gaps in the chain. Restrict reconstruction to
/// `[beginA, endA)` on contig A if given.
///
/// The implementation is rather inefficient and should be avoided.
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

/// ditto
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


private auto getPaddedAlignment(S, TranslatedTracePoint)(
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
    Get the designated records of `dbFile`. Returns all records unless
    `recordNumbers` is given.

    Returns: lazy range of designated `DbRecord`s.
    Throws: `DazzlerCommandException` if `recordNumber` is not in `dbFile`
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


/// Captures information about a single entry in a Dazzler DB.
struct DbRecord
{
    static struct PacBioReadInfo
    {
        /// ID of the well where the read occurred or zero-based contig index
        /// in the containing scaffold.
        id_t well;
        /// ditto
        alias contigIdx = well;

        /// Start of the high-quality region for reads or begin coordinate
        /// on the scaffold (including gaps) for contigs.
        coord_t pulseStart;
        alias begin = pulseStart;

        /// End of the high-quality region for reads or end coordinate
        /// on the scaffold (including gaps) for contigs.
        coord_t pulseEnd;
        alias end = pulseEnd;

        /// Read quality in [0, 1] for reads. Does not apply to contigs.
        float readQuality;


        /// Length of the high-quality region for reads or length of the
        /// contig for contigs.
        @property coord_t pulseLength() const pure nothrow @safe
        {
            return pulseEnd - pulseStart;
        }

        /// ditto
        alias length = pulseLength;
    }

    /// One-based read/contig number.
    id_t readNumber;
    /// ditto
    alias contigId = readNumber;

    /// FASTA header
    string header;

    /// PacBio read information for reads and assembly location for contigs.
    PacBioReadInfo pacBioReadInfo;
    /// ditto
    alias location = pacBioReadInfo;

    /// Base pair sequence.
    string sequence;

    /// Intrinsic QVs.
    byte[] intrinsicQualityVector;
    alias intrinsicQVs = intrinsicQualityVector;

    /// Maximum allowed QV value.
    enum maxQV = 50;


    /// Get numeric QV from character `qv`.
    static byte fromQVChar(const char qv)
    {
        if ('a' <= qv && qv <= 'z')
            return cast(byte) (qv - 'a');
        else if ('A' <= qv && qv <= 'Y')
            return cast(byte) (26 + qv - 'A');
        else
            return -1;
    }


    /// Get QV character from numeric `qv` value.
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
    Get the base pair sequences of the designated records. This does NOT
    include the header.

    Returns sequences from all records if `recordNumbers` is empty.

    Throws: `DazzlerCommandException` if `recordNumber` is not in `dbFile`
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
    Get the base pair sequence of the designated record.

    This methods caches and pre-fetches sequences in order to improve
    performance and reduce the number of forks. It uses two caches selected
    by `dbFile` which each hold up to `cacheSize` many read sequence at once.
    If `recordNumber` is not in the pre-fetched segment the next `cacheSize`
    many sequences starting with `recordNumber` are fetched.

    Params:
        dbFile = path to DB
        recordNumber = one-based number of the desired record
        cacheSize = maximum number of reads in the cache
    Throws: DazzlerCommandException if recordNumber is not in dbFile
*/
string getFastaSequence(in string dbFile, id_t recordNumber, in id_t cacheSize = 1024)
{
    // FIXME the cache size should limit the number of `char`s retrieved,
    // ie. control the memory requirements of this function
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


/// Lazily extract sequence from `S` lines in `dbdump`.
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
    Get the designated set of records in FASTA format. If `recordNumbers` is
    empty the whole DB will be converted.

    Returns: lazy range of strings
*/
auto getFastaEntries(Range)(in string dbFile, Range recordNumbers, in size_t fastaLineWidth)
if (isInputRange!Range && is(ElementType!Range : size_t))
{
    string[] dbdumpOptions = [
        DBdumpOptions.readNumber,
        DBdumpOptions.originalHeader,
        DBdumpOptions.sequenceString,
    ];

    return readDbDumpForFastaEntries(
        dbdump(dbFile, recordNumbers, dbdumpOptions),
        recordNumbers,
        fastaLineWidth,
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
    alias parseRecord = (recordLines) {
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

    return byRecordSplitter(dbDump)
        .map!parseRecord
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
    is given a temporary `.dam` file in `tmpdir` will be created.

    Returns: path to DB
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
    is given a temporary `.db` file in `tmpdir` will be created.

    Returns: path to DB
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



/// Minimum alignment coverage required to compute intrinsic quality values
/// (QVs).
enum id_t minQVCoverage = 4;


/// Compute intrinsic quality values (QVs) from the alignments in `lasFile`.
///
/// Params:
///     dbFile = path to DB
///     lasFile = alignments of entries in `dbFile`
///     masks = list of soft masks used with `DAScover`
///     coverage = coverage estimate passed to `DASqv`; estimate from
///         `DAScover` will be used unless a positive `coverage` is provided.
void computeQVs(in string dbFile, in string lasFile, id_t coverage = 0)
{
    computeQVs(dbFile, lasFile, [], coverage);
}

/// ditto
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


/// `DBdust` always produces a mask with this name.
enum dbdustMaskName = "dust";


/// Run `DBdust` on `dbFile`.
@ExternalDependency("DBdust", "DAZZ_DB", "https://github.com/thegenemyers/DAZZ_DB")
void dbdust(in string dbFile, in string[] dbdustOptions)
{
    executeCommand(chain(only("DBdust"), dbdustOptions, only(dbFile.stripDbExtension)));
}


/**
    Align DB(s) to each other using `daligner` executed in `outdir`.

    The paths `dbFile`, `referenceDb`, `queryDb` are relative to the current
    working directory – NOT `outdir`.

    Returns: path to LAS file.
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
    Align `queryDb` against `refDb` using `damapper` executed in `outdir`.

    The paths `refDb`, `queryDb` are relative to the current working
    directory – NOT `outdir`.

    Returns: path to LAS file.
*/
string getDamapping(
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


/// Filter local alignments in `lasFile` producing a new LAS file.
///
/// The `lasFile` is read using `getFlatLocalAlignments` and `pred`
/// is applied to every `FlatLocalAlignment`. The alignment is kept if `pred`
/// evaluates to a truthy value and discarded otherwise.
///
/// Read/contig lengths are filled in if `dbFile` is provided.
///
/// Bugs: LAs files referring to two different DBs cannot makes use of the
///     fill-in feature for contig lengths.
string filterLocalAlignments(alias pred)(in string lasFile)
{
    return filterLocalAlignments!pred(null, lasFile);
}

/// ditto
string filterLocalAlignments(alias pred)(in string dbFile, in string lasFile)
{
    string filteredLasFile = lasFile.stripExtension.to!string ~ "-filtered.las";
    auto flatLocalAlignments = getFlatLocalAlignments(
        dbFile,
        lasFile,
        BufferMode.overwrite,
    );

    filteredLasFile.writeAlignments(flatLocalAlignments.filter!pred);

    return filteredLasFile;
}

unittest
{
    import dentist.util.tempfile : mkdtemp;
    import std.file : rmdirRecurse;
    import std.algorithm : equal;
    import std.range : tee;

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    scope (exit)
        rmdirRecurse(tmpDir);

    auto lasFile = buildPath(tmpDir, "test.las");
    dumpLA(lasFile, testLasDump);

    auto filteredLas = filterLocalAlignments!"a.id % 2 == 0"(lasFile);

    auto flatLocalAlignments = getFlatLocalAlignments(filteredLas, BufferMode.preallocated).array;
    alias FlatLocus = FlatLocalAlignment.FlatLocus;
    auto expectedResult = [
        FlatLocalAlignment(
            0,
            FlatLocus(1, 0, 3, 4),
            FlatLocus(2, 0, 5, 6),
            AlignmentFlags(),
            expectedTracePointDistance,
            [
                TracePoint(0, 1),
            ],
        ),
        FlatLocalAlignment(
            1,
            FlatLocus(19, 0, 21, 22),
            FlatLocus(20, 0, 23, 24),
            AlignmentFlags(AlignmentFlag.complement, AlignmentFlag.alternateChain),
            expectedTracePointDistance,
            [
                TracePoint(0, 1),
            ],
        ),
        FlatLocalAlignment(
            2,
            FlatLocus(37, 0, 39, 40),
            FlatLocus(38, 0, 41, 42),
            AlignmentFlags(AlignmentFlag.unchained),
            expectedTracePointDistance,
            [
                TracePoint(0, 1),
            ],
        ),
        FlatLocalAlignment(
            3,
            FlatLocus(46, 0, 57, 58),
            FlatLocus(47, 0, 59, 60),
            AlignmentFlags(AlignmentFlag.unchained),
            expectedTracePointDistance,
            [
                TracePoint(3, 101),
                TracePoint(4, 104),
                TracePoint(2, 102),
            ],
        ),
        FlatLocalAlignment(
            4,
            FlatLocus(64, 0, 75, 76),
            FlatLocus(65, 0, 77, 78),
            AlignmentFlags(AlignmentFlag.complement, AlignmentFlag.chainContinuation),
            expectedTracePointDistance,
            [
                TracePoint(0, 2),
                TracePoint(2, 102),
            ],
        ),
        FlatLocalAlignment(
            5,
            FlatLocus(1, 0, 0, 8300),
            FlatLocus(3197, 0, 0, 318),
            AlignmentFlags(AlignmentFlag.complement, AlignmentFlag.chainContinuation),
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


/// Chain local alignments in `lasFile` and write back to a LAS file.
///
/// Returns: path to LAS file
/// See_also: `dentist.common.alignments.chaining`
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
        BufferMode.preallocated,
    );

    auto chainedAlignments = chainLocalAlignmentsAlgo(
        flatLocalAlignments,
        options,
    );

    chainedLasFile.writeAlignments(chainedAlignments);

    return chainedLasFile;
}


/// Disables alignments that should not appear in a pile up.
///
/// If using the `lasFile` signature then this reads `lasFile` with contig
/// lengths from `dbFile`, applies the filter and writes the result to the
/// returned LAS file. The resulting LAS file will have the same amount of
/// data as the input LAS file.
///
/// If using the array signature then the disabled flag of the array entries
/// is modified according to the filter.
///
/// Returns: path to LAS file or nothing.
/// Params:
///     dbFile = DB of involved reads
///     lasFile = file of alignments
///     alignments = array of alignments
///     properAlignmentAllowance = allowance when determining if an alignment
///         ends or starts at the tips of a contig.
///     forceFlat = set alignment flags to make it look unchained and sort
///         the unchained alignments. This is required if theisValidPileUpAlignment alignments are
///         chained and sorted as chains. $(I Warning: this changes the passed
///         array!)
/// See_also: `isValidPileUpAlignment`
string filterPileUpAlignments(
    in string dbFile,
    in string lasFile,
    in coord_t properAlignmentAllowance,
    Flag!"forceFlat" forceFlat = No.forceFlat,
)
{
    auto alignments = getFlatLocalAlignments(dbFile, lasFile, BufferMode.preallocated)
        .array;

    filterPileUpAlignments(alignments, properAlignmentAllowance, forceFlat);

    string filteredLasFile = lasFile.stripExtension.to!string ~ "-filtered.las";
    writeAlignments(filteredLasFile, alignments);

    return filteredLasFile;
}

/// ditto
void filterPileUpAlignments(
    ref AlignmentChain[] alignments,
    in coord_t properAlignmentAllowance,
)
{
    foreach (ref alignment; alignments)
        alignment.disableIf(!isValidPileUpAlignment(alignment, properAlignmentAllowance));
}


/// ditto
void filterPileUpAlignments(
    ref FlatLocalAlignment[] alignments,
    in coord_t properAlignmentAllowance,
    Flag!"forceFlat" forceFlat = No.forceFlat,
)
{
    foreach (ref alignment; alignments)
    {
        alignment.disableIf(!isValidPileUpAlignment(alignment, properAlignmentAllowance));

        if (forceFlat && ! alignment.flags.disabled)
        {
            alignment.flags.alternateChain = false;
            alignment.flags.chainContinuation = false;
            alignment.flags.unchained = true;
        }
    }

    if (forceFlat)
        alignments.sort;
}


/// An alignment in a pile up is valid iff it is proper and the begin/end
/// of both reads match.
///
/// Params:
///     alignment = alignment to investigate
///     allowance = allowance when determining if an alignment
///         ends or starts at the tips of a contig.
/// See_also: `dentist.common.alignments.base.AlignmentChain.beginsWith`,
///     `dentist.common.alignments.base.AlignmentChain.endsWith`,
///     `dentist.common.alignments.base.FlatLocalAlignment.FlatLocus.beginsWithin`,
///     `dentist.common.alignments.base.FlatLocalAlignment.FlatLocus.endsWithin`,
bool isValidPileUpAlignment(const ref AlignmentChain alignment, const coord_t allowance)
{
    alias isLeftAnchored = () =>
        alignment.beginsWith!"contigA"(allowance) && alignment.beginsWith!"contigB"(allowance);
    alias isLeftProper = () =>
        alignment.beginsWith!"contigA"(allowance) || alignment.beginsWith!"contigB"(allowance);
    alias isRightAnchored = () =>
        alignment.endsWith!"contigA"(allowance) && alignment.endsWith!"contigB"(allowance);
    alias isRightProper = () =>
        alignment.endsWith!"contigA"(allowance) || alignment.endsWith!"contigB"(allowance);

    return alignment.contigA.id != alignment.contigB.id && (
        (isLeftAnchored() && isRightProper()) ||
        (isRightAnchored() && isLeftProper())
    );
}

/// ditto
bool isValidPileUpAlignment(const ref FlatLocalAlignment alignment, const coord_t allowance)
{
    alias isLeftAnchored = () =>
        alignment.contigA.beginsWithin(allowance) && alignment.contigB.beginsWithin(allowance);
    alias isLeftProper = () =>
        alignment.contigA.beginsWithin(allowance) || alignment.contigB.beginsWithin(allowance);
    alias isRightAnchored = () =>
        alignment.contigA.endsWithin(allowance) && alignment.contigB.endsWithin(allowance);
    alias isRightProper = () =>
        alignment.contigA.endsWithin(allowance) || alignment.contigB.endsWithin(allowance);

    return alignment.contigA.id != alignment.contigB.id && (
        (isLeftAnchored() && isRightProper()) ||
        (isRightAnchored() && isLeftProper())
    );
}


/**
    Build consensus using `daccord`.

    If no `filteredLasFile` is provided the reads will be self-aligned using
    `daligner` and filtered with `filterPileUpAlignments`.

    If `readId` is provided the option `"-I{readId},{readId}"` will be appended
    to a copy of `options.daccordOptions`. This means the user may not provide
    `readId` and limit the operation themselves or run `daccord` on all reads.

    Params:
        dbFile          = DB of reads
        filteredLasFile = only "true" local alignments, i.e. without
            repeat-induced alignments –  as good as possible.
        readId          = the reference read used in the consensus procedure
        options         = control the various aspects of the process
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
    auto consensusFasta = getFastaEntries(consensusDb, cast(size_t[])[], options.fastaLineWidth);
    auto expectedSequence = fastaRecords[$ - 1].lineSplitter.drop(1).joiner.array;
    auto consensusSequence = consensusFasta.front.lineSplitter.drop(1).joiner.array;

    assert(expectedSequence == consensusSequence,
            format!"expected %s but got %s"(expectedSequence, consensusSequence));
}


// required by `daccord`
private void computeIntrinsticQualityValuesForConsensus(in string dbFile, in string lasFile)
{
    auto readDepth = getNumContigs(dbFile);

    computeIntrinsicQV(dbFile, lasFile, readDepth);
}


// required by `daccord`
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


/// Generate the LAS file name as `daligner`/`damapper` run in `baseDirectory`
/// does.
///
/// Uses `dbB = dbA` if no `dbB` is provided.
///
/// Params:
///     dbA           = A-read DB
///     dbB           = B-read DB
///     baseDirectory = directory where LAS file is/should be located
string getLasFile(in string dbA, in string baseDirectory)
{
    return getLasFile(dbA, null, baseDirectory);
}

/// ditto
string getLasFile(in string dbA, in string dbB, in string baseDirectory)
{
    alias dbName = dbFile => dbFile.baseName.stripDbExtension;

    enum fileTemplate = "%s/%s.%s.las";
    auto dbAName = dbName(dbA);
    auto dbBName = dbB is null ? dbAName : dbName(dbB);

    return format!fileTemplate(baseDirectory, dbAName, dbBName);
}


/// Check if a LAS file was generated by `daligner`/`damapper` run in
/// `baseDirectory`.
///
/// Uses `dbB = dbA` if no `dbB` is provided.
///
/// Params:
///     dbA           = A-read DB
///     dbB           = B-read DB
///     baseDirectory = directory where LAS file is/should be located
/// See_also: `getLasFile`
bool lasFileGenerated(in string dbA, in string baseDirectory)
{
    return lasFileGenerated(dbA, null, baseDirectory);
}

/// ditto
bool lasFileGenerated(in string dbA, in string dbB, in string baseDirectory)
{
    return getLasFile(dbA, dbB, baseDirectory).exists;
}


/// Return the number of blocks in `dbFile`.
///
/// Throws: `DazzlerCommandException` if the DB is not split or the block
///     count could not be read.
id_t getNumBlocks(in string dbFile)
{
    // see also in Dazzler's DB.h:394
    //     #define DB_NBLOCK "blocks = %9d\n"  //  number of blocks
    enum blockNumFormat = "blocks = %d";
    enum blockNumFormatPrefix = blockNumFormat[0 .. 6];
    id_t numBlocks;
    auto matchingLines = File(dbFile.stripBlock).byLine.filter!(
            line => line.startsWith(blockNumFormatPrefix));

    if (matchingLines.empty)
    {
        auto errorMessage = format!("could not read the block count in `%s`; " ~
            "maybe DBsplit is required?")(dbFile.stripBlock);
        throw new DazzlerCommandException(errorMessage);
    }

    auto matchingLine = matchingLines.front;

    if (formattedRead!blockNumFormat(matchingLine, numBlocks) != 1)
    {
        auto errorMessage = format!("could not read the block count in `%s`: " ~
            "DB header seems to be garbage")(dbFile.stripBlock);
        throw new DazzlerCommandException(errorMessage);
    }

    return numBlocks;
}


/// Return the of block size of `dbFile`.
///
/// Throws: `DazzlerCommandException` if the DB is not split or the block
///     size could not be read.
coord_t getBlockSize(in string dbFile)
{
    // see also in dazzler's DB.h:434
    //     #define DB_PARAMS "size = %11lld cutoff = %9d all = %1d\n"
    enum blockSizeFormat = "size = %d";
    enum blockSizeFormatStart = blockSizeFormat[0 .. 4];
    coord_t blockSize;
    auto matchingLines = File(dbFile.stripBlock).byLine.filter!(
            line => line.startsWith(blockSizeFormatStart));

    if (matchingLines.empty)
    {
        auto errorMessage = format!(
            "could not read the block size in `%s`; try using DBsplit to fix"
        )(dbFile.stripBlock);
        throw new DazzlerCommandException(errorMessage);
    }

    auto matchingLine = matchingLines.front;

    if (formattedRead!blockSizeFormat(matchingLine, blockSize) != 1)
    {
        auto errorMessage = format!"could not read the block count in `%s`"(dbFile.stripBlock);
        throw new DazzlerCommandException(errorMessage);
    }

    return blockSize;
}


/// Return the contig/read length cutoff of `dbFile`.
///
/// Throws: `DazzlerCommandException` if the DB is not split or the cutoff
///     could not be read.
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


/// Get the number of reads/contig in `dbFile`.
///
/// Params:
///     dbFile      = path to DB
///     untrimmedDb = return count for trimmed/untrimmed DB
id_t getNumContigs(in string dbFile)
{
    return getNumContigs(dbFile, No.untrimmedDb);
}

/// ditto
id_t getNumContigs(in string dbFile, Flag!"untrimmedDb" untrimmedDb = No.untrimmedDb)
{
    enum contigNumFormat = "+ R %d";
    enum contigNumFormatStart = contigNumFormat[0 .. 4];
    id_t numContigs;
    id_t[] empty;
    string[] options = untrimmedDb ? [DBdumpOptions.untrimmedDatabase] : [];
    auto dbdumpLines = dbdump(dbFile, empty, options);
    scope (exit) dbdumpLines.destroy();
    auto matchingLine = dbdumpLines
        .filter!(line => line.startsWith(contigNumFormatStart))
        .front;

    if (!matchingLine)
    {
        auto errorMessage = format!"could not read the contig count in `%s`"(dbFile);
        throw new DazzlerCommandException(errorMessage);
    }

    if (formattedRead!contigNumFormat(matchingLine, numContigs) != 1)
    {
        auto errorMessage = format!"could not read the contig count in `%s`"(dbFile);
        throw new DazzlerCommandException(errorMessage);
    }

    return numContigs;
}


/// Returns a lazy range of `ScaffoldSegment`s describing the scaffold
/// structure.
///
/// Returns: range of alternating `ContigSegment`s and `GapSegment`s packed in
///     `ScaffoldSegment`s.
auto getScaffoldStructure(in string dbFile)
{
    enum string[] dbshowOptions = [DBshowOptions.noSequence];

    auto rawScaffoldInfo = dbshow(dbFile, dbshowOptions);

    return ScaffoldStructureReader(rawScaffoldInfo);
}


/// Either `ContigSegment` or `GapSegment`. Used to create a mixed range of
/// the two base types.
alias ScaffoldSegment = Algebraic!(ContigSegment, GapSegment);


/// Describes the location of a contig inside its scaffold.
struct ContigSegment
{
    /// One-based contig ID in the Dazzler DB.
    size_t globalContigId;
    /// Zero-based ID of the scaffold.
    size_t scaffoldId;
    /// Zero-based ID of the contig within the scaffold.
    size_t contigId;
    /// Begin coordinate within the scaffold.
    size_t begin;
    /// End coordinate within the scaffold.
    size_t end;
    /// FASTA header of the scaffold.
    string header;

    invariant
    {
        assert(begin < end);
    }

    /// Length of the contig.
    @property size_t length() const pure nothrow
    {
        return end - begin;
    }
}

struct GapSegment
{
    /// One-based contig ID in the Dazzler DB of the contig preceding the gap.
    size_t beginGlobalContigId;
    /// One-based contig ID in the Dazzler DB of the contig following the gap.
    size_t endGlobalContigId;
    /// Zero-based ID of the scaffold.
    size_t scaffoldId;
    /// Zero-based ID of the contig preceding the gap within the scaffold.
    size_t beginContigId;
    /// Zero-based ID of the contig following the gap within the scaffold.
    size_t endContigId;
    /// Begin coordinate within the scaffold.
    size_t begin;
    /// End coordinate within the scaffold.
    size_t end;

    invariant
    {
        assert(begin < end);
    }


    /// Length of the gap.
    @property size_t length() const pure nothrow
    {
        return end - begin;
    }
}


/// Parses the output of `DBshow` to generate `ScaffoldSegment`s.
///
/// Bugs: occurrences of `" :: "` in FASTA headers cause an exception. This
///     can be fixed by relying on `DBdump` instead.
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

    Returns: tuple with fields `header` and `data`
    Params:
        dbFile = path to DB
        maskName = name of the mask. This must not contain dots (`"."`) unless
            `allowBlock` is true in which case a single dot followed by the
            block ID is allowed. Slashes (`"/`) are never allowed.
        allowBlock = whether or not `maskName` contains a block ID
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
    Read the `Interval`s of a Dazzler mask for `dbFile`.

    Params:
        Interval = target type for intervals the is default-constructible
            and has three fields `tag`, `begin` and `end`.
        dbFile   = path to DB
        maskName = name of the mask
    Throws: `MaskReaderException`
    See_Also: `writeMask`, `getMaskFiles`
*/
Interval[] readMask(Interval)(in string dbFile, in string maskName)
{
    alias _enforce = enforce!MaskReaderException;

    auto maskFileNames = getMaskFiles(dbFile, maskName, Yes.allowBlock);
    auto maskHeader = readMaskHeader(maskFileNames.header);
    auto maskData = getBinaryFile!MaskDataEntry(maskFileNames.data);

    auto maskIntervals = appender!(Interval[]);
    alias IntervalContigId = typeof(maskIntervals.data[0].tag);
    alias IntervalBegin = typeof(maskIntervals.data[0].begin);
    alias IntervalEnd = typeof(maskIntervals.data[0].end);
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

            Interval newInterval;
            newInterval.tag = currentContig.to!IntervalContigId;
            newInterval.begin = interval[0].to!IntervalBegin;
            newInterval.end = interval[1].to!IntervalEnd;

            if (trimmedDbTranslateTable.length > 0)
                newInterval.tag = trimmedDbTranslateTable[newInterval.tag - 1].to!IntervalContigId;

            if (newInterval.tag < id_t.max)
                maskIntervals ~= newInterval;
        }

        ++currentContig;
    }

    return maskIntervals.data;
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


/// Read the contents of `fileName` into a `T[]`.
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


/// Construct an translation table from untrimmed to trimmed DB IDs.
///
/// Example:
/// ---
/// auto tr = getTrimmedDbTranslateTable(dbFile);
///
/// auto untrimmedId = 42;
/// auto trimmedId = tr[untrimmedId - 1];
/// auto isContainedInTrimmedDb = trimmedId < id_t.max;
/// ---
///
/// Returns: an array of IDs in the trimmed DB. Excluded contigs/reads are
///     marked by a value of `id_t.max`.
/// Bugs: this disregards the `-alm` options to `DBsplit`, i.e. only the
///     length cutoff `-x` is applied.
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
    Write the list of intervals to a Dazzler mask for `dbFile`.

    Params:
        dbFile   = path to DB
        maskName = name of the mask
        intervals = range of intervals the have three fields `tag`, `begin`
            and `end`
    See_Also: `readMask`, `getMaskFiles`
*/
void writeMask(Intervals)(in string dbFile, in string maskName, Intervals intervals)
    if (isInputRange!Intervals)
{
    alias MaskInterval = Tuple!(
        MaskHeaderEntry, "tag",
        MaskDataEntry, "begin",
        MaskDataEntry, "end",
    );

    auto maskFileNames = getMaskFiles(dbFile, maskName, Yes.allowBlock);
    auto maskHeader = File(maskFileNames.header, "wb");
    auto maskData = File(maskFileNames.data, "wb");

    auto maskIntervals = intervals
        .map!(interval => MaskInterval(
            interval.tag.to!MaskHeaderEntry,
            interval.begin.to!MaskDataEntry,
            interval.end.to!MaskDataEntry,
        ))
        .array;
    maskIntervals.sort();

    auto numReads = getNumContigs(dbFile, No.untrimmedDb).to!MaskHeaderEntry;
    MaskHeaderEntry size = 0; // Mark the DAZZ_TRACK as a mask (see DAZZ_DB/DB.c:1183)
    MaskHeaderEntry currentContig = 1;
    MaskDataPointer dataPointer = 0;

    maskHeader.rawWrite([numReads, size]);
    maskHeader.rawWrite([dataPointer]);
    foreach (maskInterval; maskIntervals)
    {
        assert(maskInterval.tag >= currentContig);

        while (maskInterval.tag > currentContig)
        {
            maskHeader.rawWrite([dataPointer]);
            ++currentContig;
        }

        if (maskInterval.tag == currentContig)
        {
            maskData.rawWrite([maskInterval.begin, maskInterval.end]);
            dataPointer += typeof(maskInterval.begin).sizeof + typeof(maskInterval.end).sizeof;
        }
    }

    foreach (emptyContig; currentContig .. numReads + 1)
    {
        maskHeader.rawWrite([dataPointer]);
    }
}


/// Accumulation mode for merging extras of block tracks.
///
/// See_also: `DazzExtra`
enum AccumMode : int
{
    /// Extra contents must match exactly for all blocks.
    exact = 0,
    /// Extra contents are summed like vectors across blocks.
    sum = 1,
}


/// Represents a track extra, i.e. additional data stored alongside a DB track.
///
/// Use `dazzExtra` to construct objects.
///
/// See_also: `dazzExtra`, `readDazzExtra`, `writeDazzExtra`
struct DazzExtra(T) if (is(T == long) || is(T == double))
{
    /// Dynamic value type representation.
    static if (is(T == long))
        enum int vtype = 0;
    else static if (is(T == double))
        enum int vtype = 1;
    else
        static assert(0);

    /// Extra name.
    string name;
    /// Extra data. Alias this allows accessing data directly in this
    /// `DazzExtra`, e.g. `extra[0 .. 10]` returns a slice of the first 10
    /// elements.
    T[] data;
    alias data this;

    /// Defines how data from a different block-level extras is combined
    /// in `Catrack`.
    AccumMode accumMode;
}


/// Represents a track extra, i.e. additional data stored alongside a DB track.
///
/// See_also: `DazzExtra`, `readDazzExtra`, `writeDazzExtra`
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
    Read extra `extraName` of track `maskName` with elements of type `T` from
    `dbFile`.

    Returns: fully populated `DazzExtra!T`.
    Throws: $(UL
        $(LI `MaskReaderException` on read errors)
        $(LI `DazzExtraNotFound` if no extra with given name exists)
    )
    See_Also: `writeDazzExtra`, `DazzExtra`
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
        data.length == 0 || maskHeaderFile.rawRead(data).length == data.length,
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
    Write extra `extraName` of track `maskName` with elements of type `T` from
    `dbFile`.

    Use `dazzExtra` to construct `DazzExtra` objects.

    Note: the extra will be appended to the list of extras regardless of any
        existing extras. You may use `readDazzExtra` to check for existing
        extras.
    See_Also: `dazzExtra`, `readDazzExtra`
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

/// ditto
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

/// `OptionModifier.replaceOrAdd`
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

/// `OptionModifier.ensurePresent`
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

/// `OptionModifier.defaultValue`
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

/// `OptionModifier.add` and `OptionModifier.remove`
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


private struct WithOptionsImpl(R) if (isInputRange!R && is(ElementType!R == string))
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
version (LAdump)
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

    void damapper(in string refDam, in string readsDam, in string[] damapperOpts, in string outdir)
    {
        damapper([refDam, readsDam], damapperOpts, outdir);
    }

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

    version (LAdump)
    auto ladump(in string lasFile, in string dbA, in string dbB, in string[] ladumpOpts)
    {
        return ladump(lasFile, dbA, dbB, [], ladumpOpts);
    }

    version (LAdump)
    @ExternalDependency("LAdump", "DALIGNER", "https://github.com/thegenemyers/DALIGNER")
    auto ladump(
        in string lasFile,
        in string dbA,
        in string dbB,
        in id_t[] readIds,
        in string[] ladumpOpts,
    )
    {
        return executePipe!(Yes.isBuffered)(chain(
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

    version (dumpLA)
    @ExternalDependency("dumpLA", "DALIGNER", "https://github.com/thegenemyers/DALIGNER")
    void dumpLA(in string outLas, in string[] dump)
    {
        auto pipes = pipeProcess(["dumpLA", outLas], Redirect.stdin);

        foreach (dumpLine; dump)
            pipes.stdin.writeln(dumpLine);

        pipes.stdin.flush();
        pipes.stdin.close();
        wait(pipes.pid);
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


version (NoAppMain) version (ReadLasTest)
{
    void main(string[] args)
    {
        import std.datetime.stopwatch;
        import std.stdio;

        auto db = args[1];
        auto las = args[2];

        StopWatch timer;

        timer.start();
        cast(void) getFlatLocalAlignments(db, las, Yes.includeTracePoints).array;
        timer.stop();

        writefln!"getFlatLocalAlignments took %fs"(timer.peek.total!"usecs" / 1e6);

        timer.start();
        cast(void) getAlignments(db, las, AlignmentReaderFlag.includeTracePoints);
        timer.stop();

        writefln!"getAlignments took %fs"(timer.peek.total!"usecs" / 1e6);
    }
}

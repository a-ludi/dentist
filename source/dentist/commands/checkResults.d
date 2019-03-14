/**
    This is the `checkResults` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.checkResults;

import dentist.common : isTesting;

static if (isTesting):

import dentist.commandline : OptionsFor;
import dentist.common :
    dentistEnforce,
    ReferenceInterval,
    ReferenceRegion,
    ReferencePoint;
import dentist.common.alignments :
    AlignmentChain,
    coord_t,
    diff_t,
    id_t;
import dentist.common.commands : TestingCommand;
import dentist.dazzler :
    ContigSegment,
    GapSegment,
    getContigCutoff,
    getFastaSequence,
    getScaffoldStructure,
    readMask,
    ScaffoldSegment;
import dentist.util.algorithm :
    first,
    orderLexicographically,
    sliceBy;
import dentist.util.fasta : getFastaLength;
import dentist.util.log;
import dentist.util.math :
    median,
    N,
    NaturalNumberSet;
import dentist.util.process : pipeLines;
import dentist.util.range : tupleMap;
import std.algorithm :
    all,
    copy,
    filter,
    find,
    group,
    joiner,
    map,
    max,
    maxElement,
    minElement,
    sort,
    startsWith,
    sum,
    until;
import std.array : array;
import std.ascii :
    newline,
    toLower;
import std.conv : to;
import std.format :
    format,
    formattedRead;
import std.math :
    ceil,
    log10;
import std.parallelism : parallel;
import std.range :
    chain,
    enumerate,
    only,
    repeat,
    slide,
    StoppingPolicy,
    zip;
import std.range.primitives;
import std.regex :
    ctRegex,
    replaceAll;
import std.stdio :
    File,
    toFile,
    writeln;
import std.string :
    join,
    tr;
import std.typecons :
    No,
    Tuple;
import vibe.data.json :
    Json,
    toJson = serializeToJson,
    toJsonString = serializeToPrettyJson;


/// Options for the `collectPileUps` command.
alias Options = OptionsFor!(TestingCommand.checkResults);

/// Execute the `checkResults` command with `options`.
void execute(in Options options)
{
    auto analyzer = ResultAnalyzer(options);
    auto stats = analyzer.collect();

    if (options.useJson)
        writeln(stats.toJsonString());
    else
        writeln(stats.toTabular());
}

private struct ResultAnalyzer
{
    static enum identityLevels = Stats.identityLevels;

    static struct InsertionMapping
    {
        ReferenceInterval trueAssembly;
        ReferencePoint resultBegin;
        ReferencePoint resultEnd;
        id_t leftRefContig;
        id_t rightRefContig;
    }

    static struct MappingAnalysisResult
    {
        size_t[][identityLevels.length] correctGapsPerIdentityLevel;
        size_t totalDiffsInFilledGaps;
    }

    const(Options) options;
    protected const(ScaffoldSegment)[] trueAssemblyScaffoldStructure;
    protected const(ScaffoldSegment)[] resultScaffoldStructure;
    protected coord_t referenceOffset;
    protected ReferenceRegion mappedRegionsMask;
    protected ReferenceRegion referenceGaps;
    protected RawMummerAlignment[] contigAlignments;
    protected InsertionMapping[] insertionMappings;
    protected size_t[][identityLevels.length] correctGapsPerIdentityLevel;
    protected size_t totalDiffsInFilledGaps;

    Stats collect()
    {
        mixin(traceExecution);

        init();

        Stats stats;

        stats.numBpsExpected = getNumBpsExpected();
        stats.numBpsKnown = getNumBpsKnown();
        stats.numBpsResult = getNumBpsResult();
        stats.numBpsInGaps = getNumBpsInGaps();
        stats.numDiffsInGaps = totalDiffsInFilledGaps;
        stats.numTranslocatedGaps = getNumTranslocatedGaps();
        stats.numContigsExpected = getNumContigsExpected();
        stats.numMappedContigs = getNumMappedContigs();
        stats.numCorrectGaps = getNumCorrectGaps();
        stats.numClosedGaps = getNumClosedGaps();
        stats.numPartlyClosedGaps = getNumPartlyClosedGaps();
        stats.maximumN50 = getMaximumN50();
        stats.inputN50 = getInputN50();
        stats.resultN50 = getResultN50();
        stats.gapMedian = getGapMedian();
        stats.closedGapMedian = getClosedGapMedian();
        stats.minClosedGap = getExtremumClosedGap!minElement();
        stats.maxClosedGap = getExtremumClosedGap!maxElement();
        if (options.bucketSize > 0)
        {
            stats.correctGapLengthHistograms = getCorrectGapLengthHistograms();
            stats.closedGapLengthHistogram = getClosedGapLengthHistogram();
            stats.gapLengthHistogram = getGapLengthHistogram();
        }

        return stats;
    }

    void init()
    {
        trueAssemblyScaffoldStructure = getScaffoldStructure(options.trueAssemblyDb).array;
        resultScaffoldStructure = getScaffoldStructure(options.resultDb).array;
        auto contigCutoff = getContigCutoff(options.refDb);
        mappedRegionsMask = ReferenceRegion(readMask!ReferenceInterval(
            options.trueAssemblyDb,
            options.mappedRegionsMask,
            options.workdir,
        ).filter!(interval => interval.size >= contigCutoff).array);
        referenceOffset = cast(coord_t) mappedRegionsMask.intervals[0].begin;

        contigAlignments = findReferenceContigs();
        insertionMappings = getInsertionMappings(contigAlignments);

        logJsonDiagnostic(
            "referenceOffset", referenceOffset,
            "numContigAlignments", contigAlignments.length.toJson,
            "contigAlignments", shouldLog(LogLevel.debug_)
                ? contigAlignments.toJson
                : toJson(null),
            "numInsertionMappings", insertionMappings.length.toJson,
            "insertionMappings", shouldLog(LogLevel.debug_)
                ? insertionMappings.toJson
                : toJson(null),
        );

        referenceGaps = getReferenceGaps();

        debug logJsonDebug("referenceGaps", referenceGaps.toJson);

        auto analysisResult = analyzeMappingQualities(
            insertionMappings,
            options.gapDetailsTabular !is null
                ? File(options.gapDetailsTabular, "w")
                : File(),
        );
        correctGapsPerIdentityLevel = analysisResult.correctGapsPerIdentityLevel;
        totalDiffsInFilledGaps = analysisResult.totalDiffsInFilledGaps;
    }

    RawMummerAlignment[] findReferenceContigs()
    {
        auto duplicateContigIds = NaturalNumberSet.create(
            findPerfectAlignments(options.refDb, null, options.workdir)
                .map!(RawMummerAlignment.fromString)
                .map!(selfAlignment => cast(size_t) selfAlignment.refContigId)
                .array
        );
        logJsonDiagnostic(
            "numDuplicateContigIds", duplicateContigIds.size,
            "duplicateContigIds", shouldLog(LogLevel.debug_)
                ? duplicateContigIds.elements.array.toJson
                : toJson(null),
        );

        auto perfectContigAlignments = findPerfectAlignments(options.refDb, options.resultDb, options.workdir)
            .map!(RawMummerAlignment.fromString)
            // Ignore contigs that have exact copies in `refDb`
            .filter!(rawAlignment => rawAlignment.refContigId !in duplicateContigIds)
            .array;

        auto bufferRest = perfectContigAlignments
            .sort!((lhs, rhs) => lhs.refContigId < rhs.refContigId)
            .group!((lhs, rhs) => lhs.refContigId == rhs.refContigId)
            .filter!(group => group[1] == 1)
            .map!(group => group[0])
            .copy(perfectContigAlignments);
        perfectContigAlignments.length -= bufferRest.length;

        alias resultOrder = orderLexicographically!(RawMummerAlignment,
            rawAlignment => rawAlignment.resultContigId,
            rawAlignment => rawAlignment.resultBegin,
            rawAlignment => rawAlignment.refContigId,
        );

        return perfectContigAlignments
            .sort!resultOrder
            .release;
    }

    auto getInsertionMappings(in RawMummerAlignment[] contigAlignments)
    {
        static bool isGapClosed(in RawMummerAlignment lhs, in RawMummerAlignment rhs)
        {
            return lhs.resultContigId == rhs.resultContigId;
        }

        bool isGapPartlyClosed(in RawMummerAlignment lhs, in RawMummerAlignment rhs)
        {
            return lhs.resultContigId + 1 == rhs.resultContigId &&
                   lhs.resultEnd < lhs.resultContigLength &&
                   1 < rhs.resultBegin &&
                   resultGapSize(lhs.resultContigId) > 0;
        }

        bool isConsecutiveContigs(in RawMummerAlignment lhs, in RawMummerAlignment rhs)
        {
            return lhs.refContigId + 1 == rhs.refContigId &&
                   trueAssemblyContig(lhs).contigId == trueAssemblyContig(rhs).contigId;
        }

        InsertionMapping getInsertionMapping(in RawMummerAlignment lhs, in RawMummerAlignment rhs)
        {
            auto trueAssemblyInterval = ReferenceInterval(
                trueAssemblyContig(lhs).contigId,
                trueAssemblyContig(lhs).end,
                trueAssemblyContig(rhs).begin,
            );
            ReferencePoint resultBegin;
            ReferencePoint resultEnd;

            if (
                (isGapClosed(lhs, rhs) && lhs.resultEnd < rhs.resultBegin) ||
                isGapPartlyClosed(lhs, rhs)
            )
            {
                resultBegin = ReferencePoint(
                    lhs.resultContigId,
                    lhs.resultEnd,
                );
                resultEnd = ReferencePoint(
                    rhs.resultContigId,
                    rhs.resultBegin - 1,
                );
            }

            return InsertionMapping(
                trueAssemblyInterval,
                resultBegin,
                resultEnd,
                lhs.refContigId,
                rhs.refContigId,
            );
        }

        return contigAlignments
            .slide!(No.withPartial)(2)
            .filter!(mappedRegionPair => isGapClosed(mappedRegionPair[0], mappedRegionPair[1]) ||
                                         isGapPartlyClosed(mappedRegionPair[0], mappedRegionPair[1]))
            .filter!(mappedRegionPair => isConsecutiveContigs(mappedRegionPair[0], mappedRegionPair[1]))
            .map!(mappedRegionPair => getInsertionMapping(mappedRegionPair[0], mappedRegionPair[1]))
            .array;
    }

    ReferenceInterval trueAssemblyContig(in RawMummerAlignment mapping)
    {
        return mappedRegionsMask.intervals[mapping.refContigId - 1];
    }

    coord_t resultGapSize(in id_t contigId)
    {
        auto findGap = resultScaffoldStructure
            .filter!(gapPart => gapPart.peek!GapSegment !is null)
            .map!(gapPart => gapPart.get!GapSegment)
            .find!(gapPart => gapPart.beginGlobalContigId == contigId);

        if (findGap.empty)
            return 0;
        else
            return cast(coord_t) findGap.front.length;
    }

    ReferenceRegion getReferenceGaps()
    {
        mixin(traceExecution);

        alias scaffoldRegion = (contigPart) => ReferenceRegion(ReferenceInterval(
            contigPart.globalContigId,
            0,
            contigPart.length,
        ));
        alias mappedGapsMask = (contigPart) => scaffoldRegion(contigPart) - mappedRegionsMask;

        return ReferenceRegion(trueAssemblyScaffoldStructure
            .filter!(contigPart => contigPart.peek!ContigSegment !is null)
            .map!(contigPart => contigPart.get!ContigSegment)
            .map!(contigPart => cast(ReferenceInterval[]) mappedGapsMask(contigPart).intervals)
            .joiner
            .filter!(gap => isInnerGap(gap))
            .array
        );
    }

    coord_t getResultContigBegin(id_t contigId)
    {
        return cast(coord_t) (resultScaffoldStructure
            .filter!(contigPart => contigPart.peek!ContigSegment !is null)
            .map!(contigPart => contigPart.peek!ContigSegment)
            .find!(contigPart => contigPart.globalContigId == contigId)
            .map!(contigPart => contigPart.begin)
            .front + 0);
    }

    size_t getNumBpsExpected()
    {
        mixin(traceExecution);

        return mappedRegionsMask
            .intervals
            .sliceBy!"a.contigId == b.contigId"
            .map!(trueScaffoldSlice => trueScaffoldSlice[$ - 1].end - trueScaffoldSlice[0].begin)
            .sum;
    }

    size_t getNumBpsKnown()
    {
        mixin(traceExecution);

        return mappedRegionsMask.size;
    }

    size_t getNumBpsResult()
    {
        mixin(traceExecution);

        return resultScaffoldStructure
            .filter!(contigPart => contigPart.peek!ContigSegment !is null)
            .map!(contigPart => contigPart.peek!ContigSegment.length)
            .sum;
    }

    size_t getNumBpsInGaps()
    {
        mixin(traceExecution);

        return insertionMappings
            .map!(mapping => mapping.trueAssembly.size)
            .sum;
    }

    size_t getNumTranslocatedGaps()
    {
        mixin(traceExecution);

        return referenceGaps.intervals.length;
    }

    size_t getNumContigsExpected()
    {
        mixin(traceExecution);

        return mappedRegionsMask.intervals.length;
    }

    size_t getNumMappedContigs()
    {
        return contigAlignments.length;
    }

    size_t getNumCorrectGaps()
    {
        return correctGapsPerIdentityLevel[0].length;
    }

    size_t getNumClosedGaps()
    {
        mixin(traceExecution);

        return insertionMappings
            .filter!(mapping => mapping.resultBegin.contigId == mapping.resultEnd.contigId)
            .walkLength;
    }

    size_t getNumPartlyClosedGaps()
    {
        mixin(traceExecution);

        return insertionMappings
            .filter!(mapping => mapping.resultBegin.contigId < mapping.resultEnd.contigId)
            .walkLength;
    }

    MappingAnalysisResult analyzeMappingQualities(in InsertionMapping[] insertionMappings, File detailsTabular = File())
    {
        mixin(traceExecution);

        if (detailsTabular.isOpen)
            detailsTabular.writefln!"%10s %10s %10s %10s %10s %s %s"(
                "contigId",
                "begin",
                "end",
                "numDiffs",
                "percentIdentity",
                "sequenceAlignment",
                "gapId",
            );

        size_t[][identityLevels.length] identicalLengthsPerLevel;
        size_t totalDiffs;
        foreach (identicalLengths; identicalLengthsPerLevel)
            identicalLengths.reserve(insertionMappings.length);

        foreach (insertionMapping; parallel(insertionMappings))
        {
            auto insertionAlignment = computeInsertionQuality(insertionMapping);

            if (detailsTabular.isOpen)
                synchronized detailsTabular.writefln!"%10d %10d %10d %10d %.8f %s %s"(
                    insertionMapping.trueAssembly.contigId,
                    insertionMapping.trueAssembly.begin,
                    insertionMapping.trueAssembly.end,
                    insertionAlignment.numDiffs,
                    insertionAlignment.percentIdentity,
                    insertionAlignment.alignmentString.replaceAll(ctRegex!("\n", `g`), `\n`),
                    getGapId(insertionMapping),
                );

            foreach (i, identityLevel; identityLevels)
                if (identityLevel < insertionAlignment.percentIdentity)
                    synchronized
                    {
                        identicalLengthsPerLevel[i] ~= insertionMapping.trueAssembly.size;
                        totalDiffs += insertionAlignment.numDiffs;
                    }
        }

        return MappingAnalysisResult(identicalLengthsPerLevel, totalDiffs);
    }

    StretcherAlignment computeInsertionQuality(in InsertionMapping insertionMapping)
    {
        mixin(traceExecution);

        auto gapId = getGapId(insertionMapping);
        auto trueSequenceFile = trueAssemblySubseqFasta(insertionMapping.trueAssembly, gapId);
        auto insertedSequenceFile = resultSubseqFasta(
            insertionMapping.resultBegin,
            insertionMapping.resultEnd,
            gapId,
        );

        return stretcher(trueSequenceFile, insertedSequenceFile);
    }

    static string getGapId(in InsertionMapping insertionMapping)
    {
        return format!"%d-%d@%d-%d"(
            insertionMapping.resultBegin.contigId,
            insertionMapping.resultEnd.contigId,
            insertionMapping.leftRefContig,
            insertionMapping.rightRefContig,
        );
    }

    enum subseqFastaFileFormat = "%s/insertion-%s.%s.fasta";
    enum subseqFastaHeaderFormat = ">%s-%s [contig-%d@%d, contig-%d@%d)";

    string trueAssemblySubseqFasta(in ReferenceInterval interval, in string gapId)
    {
        enum origin = "true";
        auto fastaFile = format!subseqFastaFileFormat(options.workdir, gapId, origin);
        auto fastaHeader = format!subseqFastaHeaderFormat(
            origin,
            gapId,
            interval.contigId,
            interval.begin,
            interval.contigId,
            interval.end,
        );
        auto subSequence = options.trueAssemblyDb.getFastaSequence(
            cast(id_t) interval.contigId,
            options.workdir,
        )[interval.begin .. interval.end];

        chain(
            fastaHeader,
            newline,
            subSequence,
            newline,
        ).map!(to!(immutable(char))).toFile(fastaFile);

        return fastaFile;
    }

    string resultSubseqFasta(in ReferencePoint begin, in ReferencePoint end, in string gapId)
    {
        enum origin = "inserted";
        auto fastaFile = format!subseqFastaFileFormat(options.workdir, gapId, origin);
        auto fastaHeader = format!subseqFastaHeaderFormat(
            origin,
            gapId,
            begin.contigId,
            begin.value,
            end.contigId,
            end.value,
        );
        auto dbFile = origin == "true"
            ? options.trueAssemblyDb
            : options.resultDb;
        string subSequence;

        if (begin.contigId == end.contigId)
        {
            subSequence = options.resultDb.getFastaSequence(
                cast(id_t) begin.contigId,
                options.workdir,
            )[begin.value .. end.value];
        }
        else
        {
            auto leftFlank = options.resultDb.getFastaSequence(
                cast(id_t) begin.contigId,
                options.workdir,
            )[begin.value .. $];
            auto rightFlank = options.resultDb.getFastaSequence(
                cast(id_t) end.contigId,
                options.workdir,
            )[0 .. end.value];
            auto gapSize = resultGapSize(cast(id_t) begin.contigId);

            enum char unkownBase = 'n';
            subSequence = chain(
                leftFlank,
                unkownBase.repeat(gapSize),
                rightFlank,
            ).map!(to!(immutable(char))).array;
        }

        chain(
            fastaHeader,
            newline,
            subSequence,
            newline,
        ).map!(to!(immutable(char))).toFile(fastaFile);

        return fastaFile;
    }

    bool isInnerGap(in ReferenceInterval gap)
    {
        return 0 < gap.begin && gap.end < trueScaffoldLength(cast(id_t) gap.contigId);
    }

    coord_t trueScaffoldLength(in id_t contigId)
    {
        return cast(coord_t) trueAssemblyScaffoldStructure
            .filter!(contigPart => contigPart.peek!ContigSegment !is null)
            .map!(contigPart => contigPart.get!ContigSegment)
            .find!(contigPart => contigPart.globalContigId == contigId)
            .first
            .length;
    }

    size_t getMaximumN50()
    {
        mixin(traceExecution);

        return trueAssemblyScaffoldStructure
            .filter!(contigPart => contigPart.peek!ContigSegment !is null)
            .map!(contigPart => contigPart.peek!ContigSegment.length)
            .array
            .N!50(getNumBpsExpected);
    }

    size_t getInputN50()
    {
        mixin(traceExecution);

        return mappedRegionsMask
            .intervals
            .map!(mappedInterval => mappedInterval.size)
            .array
            .N!50(getNumBpsExpected);
    }

    size_t getResultN50()
    {
        mixin(traceExecution);

        return resultScaffoldStructure
            .filter!(contigPart => contigPart.peek!ContigSegment !is null)
            .map!(contigPart => contigPart.peek!ContigSegment.length)
            .array
            .N!50(getNumBpsExpected);
    }

    size_t getGapMedian()
    {
        mixin(traceExecution);

        auto values = referenceGaps
            .intervals
            .map!(gap => gap.size)
            .array;

        if (values.length == 0)
            return size_t.max;
        else
            return median(values);
    }

    size_t getClosedGapMedian()
    {
        mixin(traceExecution);

        auto values = insertionMappings
            .map!(mapping => mapping.trueAssembly.size)
            .array;

        if (values.length == 0)
            return size_t.max;
        else
            return median(values);
    }

    size_t getExtremumClosedGap(alias extremeElement)()
    {
        mixin(traceExecution);

        auto values = insertionMappings.map!(mapping => mapping.trueAssembly.size);

        if (values.length == 0)
            return size_t.max;
        else
            return extremeElement(values);
    }

    Histogram!coord_t[identityLevels.length] getCorrectGapLengthHistograms()
    {
        mixin(traceExecution);

        typeof(return) correctGapLengthHistograms;

        foreach (i; 0 .. identityLevels.length)
        {
            correctGapLengthHistograms[i] = histogram(
                options.bucketSize,
                correctGapsPerIdentityLevel[i],
            );
            // NOTE: for some reason bucketSize is not correctly set; force it
            correctGapLengthHistograms[i].bucketSize = options.bucketSize;
        }

        return correctGapLengthHistograms;
    }

    Histogram!coord_t getClosedGapLengthHistogram()
    {
        mixin(traceExecution);

        return histogram(
            options.bucketSize,
            insertionMappings
                .map!(insertionMapping => cast(coord_t) insertionMapping.trueAssembly.size)
                .array,
        );
    }

    Histogram!coord_t getGapLengthHistogram()
    {
        mixin(traceExecution);

        return histogram(
            options.bucketSize,
            referenceGaps
                .intervals
                .map!(gap => cast(coord_t) gap.size)
                .array,
        );
    }
}

private struct Histogram(value_t)
{
    alias Bucket = Tuple!(
        value_t, "begin",
        value_t, "end",
        size_t, "count",
    );

    value_t bucketSize;
    size_t[] histogram;
    alias histogram this;

    this(in value_t bucketSize, in value_t[] values)
    {
        this.bucketSize = bucketSize + 0;

        if (values.length == 0)
            return; // empty histogram

        auto largestValue = maxElement(values);
        this.length = largestValue / bucketSize + 1;

        foreach (value; values)
            ++this[value / bucketSize];
    }

    auto buckets() const pure nothrow
    {
        return histogram
            .enumerate
            .map!(rawBucket => Bucket(
                bucketSize * cast(value_t) rawBucket.index,
                bucketSize * cast(value_t) (rawBucket.index + 1),
                rawBucket.value + 0,
            ));
    }
}

auto histogram(value_t)(in value_t bucketSize, in value_t[] values)
{
    return Histogram!value_t(bucketSize, values);
}

private struct Stats
{
    static enum identityLevels = [.999, .99, .95, .90];

    size_t numBpsExpected;
    size_t numBpsKnown;
    size_t numBpsResult;
    size_t numBpsInGaps;
    size_t numDiffsInGaps;
    size_t numTranslocatedGaps;
    size_t numCorrectGaps;
    size_t numContigsExpected;
    size_t numMappedContigs;
    size_t numClosedGaps;
    size_t numPartlyClosedGaps;
    size_t maximumN50;
    size_t inputN50;
    size_t resultN50;
    size_t gapMedian;
    size_t closedGapMedian;
    size_t minClosedGap;
    size_t maxClosedGap;
    Histogram!coord_t[identityLevels.length] correctGapLengthHistograms;
    Histogram!coord_t closedGapLengthHistogram;
    Histogram!coord_t gapLengthHistogram;

    @property double averageInsertionError() const pure
    {
        return numDiffsInGaps.to!double / numBpsInGaps.to!double;
    }

    string toJsonString() const
    {
        return [
            "numBpsExpected": numBpsExpected.toJson,
            "numBpsKnown": numBpsKnown.toJson,
            "numBpsResult": numBpsResult.toJson,
            "numBpsInGaps": numBpsInGaps.toJson,
            "numDiffsInGaps": numDiffsInGaps.toJson,
            "averageInsertionError": averageInsertionError.toJson,
            "numTranslocatedGaps": numTranslocatedGaps.toJson,
            "numCorrectGaps": numCorrectGaps.toJson,
            "numContigsExpected": numContigsExpected.toJson,
            "numMappedContigs": numMappedContigs.toJson,
            "numClosedGaps": numClosedGaps.toJson,
            "numPartlyClosedGaps": numPartlyClosedGaps.toJson,
            "maximumN50": maximumN50.toJson,
            "inputN50": inputN50.toJson,
            "resultN50": resultN50.toJson,
            "gapMedian": gapMedian.toJson,
            "closedGapMedian": closedGapMedian.toJson,
            "minClosedGap": minClosedGap.toJson,
            "maxClosedGap": maxClosedGap.toJson,
            "gapLengthHistogram": histsToJson(
                correctGapLengthHistograms[0],
                correctGapLengthHistograms[1],
                correctGapLengthHistograms[2],
                correctGapLengthHistograms[3],
                closedGapLengthHistogram,
                gapLengthHistogram,
            ),
        ].toJsonString;
    }

    string toTabular() const
    {
        auto columnWidth = columnWidth();
        auto rows =[
            format!"numContigsExpected:    %*d"(columnWidth, numContigsExpected),
            format!"numMappedContigs:      %*d"(columnWidth, numMappedContigs),

            format!"numTranslocatedGaps:   %*d"(columnWidth, numTranslocatedGaps),
            format!"numClosedGaps:         %*d"(columnWidth, numClosedGaps),
            format!"numPartlyClosedGaps:   %*d"(columnWidth, numPartlyClosedGaps),
            format!"numCorrectGaps:        %*d"(columnWidth, numCorrectGaps),

            format!"numBpsExpected:        %*d"(columnWidth, numBpsExpected),
            format!"numBpsKnown:           %*d"(columnWidth, numBpsKnown),
            format!"numBpsResult:          %*d"(columnWidth, numBpsResult),
            format!"numBpsInGaps:          %*d"(columnWidth, numBpsInGaps),
            format!"numDiffsInGaps:        %*d"(columnWidth, numDiffsInGaps),
            format!"averageInsertionError: %*.3e"(columnWidth, averageInsertionError),

            format!"maximumN50:            %*d"(columnWidth, maximumN50),
            format!"inputN50:              %*d"(columnWidth, inputN50),
            format!"resultN50:             %*d"(columnWidth, resultN50),
            format!"gapMedian:             %*d"(columnWidth, gapMedian),
            format!"closedGapMedian:       %*d"(columnWidth, closedGapMedian),
            format!"minClosedGap:          %*s"(columnWidth, minClosedGap),
            format!"maxClosedGap:          %*s"(columnWidth, maxClosedGap),
            gapLengthHistogram.length > 0
                ? format!"gapLengthHistogram:\n    \n%s"(histsToString(
                    correctGapLengthHistograms[0],
                    correctGapLengthHistograms[1],
                    correctGapLengthHistograms[2],
                    correctGapLengthHistograms[3],
                    closedGapLengthHistogram,
                    gapLengthHistogram,
                ))
                : null,
        ];

        return format!"%-(%s\n%)"(rows.filter!"a !is null");
    }

    const string histsToString(Hists...)(in Hists hists)
    {
        assert(hists.length >= 1);
        assert([hists].all!(h => h.bucketSize == hists[0].bucketSize));
        auto limitsWidth = numWidth(hists[0].bucketSize * [hists].map!"a.length".maxElement);
        auto countWidths = [hists]
            .map!(h => numWidth(max(1, h.histogram.length == 0 ? 1 : maxElement(h.histogram))))
            .array;

        auto bucketsList = tupleMap!(h => h.buckets.array)(hists);
        auto rows = zip(StoppingPolicy.longest, bucketsList.expand)
            .map!(buckets => format!"    %*d%s"(
                limitsWidth,
                max(tupleMap!"a.end"(buckets.expand).expand),
                zip(countWidths, [buckets.expand])
                    .map!(pair => format!" %*d"(pair[0], pair[1].count))
                    .join,
            ))
            .array;

        return format!"%-(%s\n%)"(rows);
    }

    const Json histsToJson(Hists...)(in Hists hists)
    {
        assert(hists.length >= 1);
        assert([hists].all!(h => h.bucketSize == hists[0].bucketSize));

        auto bucketsList = tupleMap!(h => h.buckets.array)(hists);
        auto rows = zip(StoppingPolicy.longest, bucketsList.expand)
            .map!(buckets => [
                "limit": max(tupleMap!"a.end"(buckets.expand).expand).toJson,
                "counts": [buckets.expand]
                    .map!"a.count"
                    .array
                    .toJson,
            ])
            .array;

        return rows.toJson;
    }

    protected size_t columnWidth() const
    {
        auto numWidth = numWidth(max(
            numBpsExpected,
            numBpsKnown,
            numBpsResult,
            numBpsInGaps,
            numDiffsInGaps,
            numTranslocatedGaps,
            numCorrectGaps,
            numContigsExpected,
            numMappedContigs,
            numClosedGaps,
            numPartlyClosedGaps,
            maximumN50,
            inputN50,
            resultN50,
            gapMedian,
            closedGapMedian,
            minClosedGap,
            maxClosedGap,
        ));
        auto strWidth = only(
            format!"%.3e"(averageInsertionError).length,
        ).maxElement;

        return max(numWidth, strWidth);
    }

    static protected size_t numWidth(Int)(Int i) nothrow
    {
        if (i == 0)
            return 1;

        auto numWidth = cast(size_t) ceil(log10(i));

        return numWidth;
    }
}

struct RawMummerAlignment
{
    coord_t resultBegin;
    coord_t resultEnd;
    coord_t refBegin;
    coord_t refEnd;
    double sequenceIdentity;
    coord_t resultContigLength;
    coord_t refContigLength;
    id_t resultContigId;
    id_t refContigId;

    static RawMummerAlignment fromString(in string dumpLine) pure
    {
        enum dumpLineFormat = "%d %d %d %d %d %d %f %d %d result-%d ref-%d";
        RawMummerAlignment rawAlignment;
        coord_t dummy;

        dumpLine[].formattedRead!dumpLineFormat(
            rawAlignment.resultBegin,
            rawAlignment.resultEnd,
            rawAlignment.refBegin,
            rawAlignment.refEnd,
            dummy,
            dummy,
            rawAlignment.sequenceIdentity,
            rawAlignment.resultContigLength,
            rawAlignment.refContigLength,
            rawAlignment.resultContigId,
            rawAlignment.refContigId,
        );

        return rawAlignment;
    }

    T to(T : AlignmentChain)() const pure
    {
        import std.conv : to;

        T ac;
        T.LocalAlignment la;

        ac.contigA.id = refContigId;
        ac.contigA.length = refContigId;
        ac.contigB.id = resultContigId;
        ac.contigB.length = resultContigId;

        la.contigA.begin = refBegin - 1;
        la.contigA.end = refEnd;
        la.contigB.begin = resultBegin - 1;
        la.contigB.end = resultEnd;
        la.numDiffs = to!diff_t(la.contigB.length * sequenceIdentity / 100.0);

        ac.localAlignments = [la];

        // ensure alignment chain is valid...
        version (assert)
            // trigger invariant of AlignmentChain
            cast(void) ac.first;

        return ac;
    }
}

enum findPerfectAlignmentsScript = q"{
    #!/bin/bash

    function main()
    {
        while getopts "kT:" OPTION; do
            case "$OPTION" in
                k)
                    KEEP_TEMP=1
                    ;;
                T)
                    TEMP_DIR="$OPTARG"
                    ;;
                *)
                    echo "$(basename "$0"): unkown option $OPTION" >&2
                    echo >&2
                    usage
                    ;;
            esac
        done
        shift $(($OPTIND - 1))

        REFERENCE="$(realpath "$1")"
        shift
        TARGETS=($*)
        for (( i = 0; i < ${#TARGETS[*]}; i++ )); do
            TARGETS[i]="$(realpath "${TARGETS[i]}")"
        done

        WORKDIR="$(mktemp -d find-perfect-alignments.XXXXXX)"
        cd "$WORKDIR"

        REFERENCE_CONTIGS="$(assert_contigs_fasta "$REFERENCE" 'ref')"
        for TARGET in ${TARGETS[*]}; do
            TARGET_CONTIGS="$(assert_contigs_fasta "$TARGET" 'result')"
            PREFIX="$(db_name "$TARGET")"

            nucmer \
                --threads=8 \
                --forward \
                --prefix=$PREFIX \
                "$TARGET_CONTIGS" \
                "$REFERENCE_CONTIGS"
            delta-filter -i 100 -q "$PREFIX.delta" > "$PREFIX.F.delta"

            if [[ $REFERENCE == $TARGET ]]; then
                FILTER='{ if ($3 < $4 && ($5 == $9 || $6 == $9) && $10 != $11) print $0 }'
            else
                FILTER='{ if ($3 < $4 && $5 == $6 && $6 == $9) print $0 }'
            fi

            show-coords -l -T -q "$PREFIX.F.delta" | \
                tail -n+5 | \
                awk -F'\t' "$FILTER"
        done
    }

    function assert_contigs_fasta()
    {
        DB="$1"
        CONTIGS_FASTA="$(db_name "$1").contigs.fasta"

        contigs_fasta "$DB" "$2" > "$CONTIGS_FASTA"

        echo "$CONTIGS_FASTA"
    }

    function contigs_fasta()
    {
        contig_dump "$1" | awk -F'\t' '{ printf ">'"$2"'-%09d %s\n%s\n", $1, substr($2, 2), $3 }'
    }

    function contig_dump()
    {
        DBdump -rhs "$1" | sed -nE '/^R/ {s/^R\s+//;h;n; s/^H[^>]+//;H;n;n; s/^S[^actg]+//;H;x;s/\n/\t/g; p}'
    }

    function db_name()
    {
        basename "$1" ".dam"
    }

    function mktemp()
    {
        if ! [[ -v __TEMP_FILES ]]; then
            declare -a __TEMP_FILES
        fi

        declare -a OPTS=()
        if [[ -v TEMP_DIR ]]; then
            OPTS+=("--tmpdir=$TEMP_DIR")
        else
            OPTS+=("--tmpdir")
        fi

        local TMPFILE="$(/usr/bin/env mktemp "${OPTS[@]}" "$@")"

        local IFS=$'\n'
        __TEMP_FILES+=($TMPFILE)

        echo -n "$TMPFILE"
    }

    function __cleanup()
    {
        if [[ ! -v KEEP_TEMP && -v __TEMP_FILES ]]; then
            local IFS=$'\n'

            for TMPFILE in ${__TEMP_FILES[*]}; do
                rm -rf "$TMPFILE"
            done

            unset __TEMP_FILES
        fi
    }

    trap __cleanup exit err

    function alternate_ifs()
    {
        OLD_IFS="$IFS"
        IFS="$1"
    }

    function reset_ifs()
    {
        IFS="$OLD_IFS"
        unset OLD_IFS
    }

    if ! [[ $- =~ i ]]; then
        # Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
        set -euo pipefail
        IFS=$'\n'

        main "$@"
    fi
}";

auto findPerfectAlignments(in string refDb, in string queryDb = null, in string tmpDir = null)
{
    auto findPerfectAlignmentsCmd = only(
        "bash",
        "-c",
        findPerfectAlignmentsScript,
        "--",
        "-k",
        tmpDir !is null ? "-T" ~ tmpDir : null,
        refDb,
        queryDb !is null ? queryDb : refDb,
    );

    return findPerfectAlignmentsCmd.pipeLines();
}

struct StretcherAlignment
{
    diff_t numIdentical;
    coord_t length;
    string alignmentString;

    @property diff_t numDiffs() const pure nothrow
    {
        return length - numIdentical;
    }

    @property double percentIdentity() const pure
    {
        return numIdentical.to!double / length.to!double;
    }
}

StretcherAlignment stretcher(in string refFasta, in string queryFasta)
{
    if (refFasta.getFastaLength() == 0)
        return StretcherAlignment();
    if (queryFasta.getFastaLength() == 0)
        return StretcherAlignment(0, cast(coord_t) refFasta.getFastaLength());

    auto stretcherCmd = only(
        "stretcher",
        "--auto",
        "--stdout",
        "--aformat=pair",
        "--awidth=" ~ uint.max.to!string,
        refFasta,
        queryFasta,
    );

    StretcherAlignment alignment;
    auto stretcherResult = stretcherCmd.pipeLines();

    stretcherResult.stretcherReadLength(alignment);
    stretcherResult.stretcherReadIdentity(alignment);
    stretcherResult.stretcherReadAlignmentString(alignment);

    return alignment;
}

private void stretcherReadLength(StretcherResult)(ref StretcherResult stretcherResult, ref StretcherAlignment alignment)
{
    enum lengthLineFormat = "# Length: %d";
    enum lengthLinePrefix = lengthLineFormat.until(':', No.openRight);

    auto findLengthLine = stretcherResult.find!(line => line.startsWith(lengthLinePrefix));

    dentistEnforce(!findLengthLine.empty, "Missing length line in stretcher output");

    findLengthLine.front.formattedRead!lengthLineFormat(alignment.length);
}

private void stretcherReadIdentity(StretcherResult)(ref StretcherResult stretcherResult, ref StretcherAlignment alignment)
{
    enum identityLineFormat = "# Identity: %d/%d";
    enum identityLinePrefix = identityLineFormat.until(':', No.openRight);
    auto findIdentityLine = stretcherResult.find!(line => line.startsWith(identityLinePrefix));

    dentistEnforce(!findIdentityLine.empty, "Missing identity line in stretcher output");

    findIdentityLine.front.formattedRead!identityLineFormat(
        alignment.numIdentical,
        alignment.length,
    );
}

private void stretcherReadAlignmentString(StretcherResult)(ref StretcherResult stretcherResult, ref StretcherAlignment alignment)
{
    enum commentChar = '#';

    auto alignmentLines = stretcherResult
        .filter!(line => line.length > 0 && line[0] != commentChar)
        .array;

    dentistEnforce(alignmentLines.length > 0, "Missing alignment in stretcher output");
    dentistEnforce(alignmentLines.length == 3, "Missing alignment in stretcher output");

    enum sequenceRegex = ctRegex!`^\s*(?P<prefix>(?:true|inserted)-[0-9@-]+\s+\d+\s+)(?P<seq>[ACTGN-]+).*$`;
    alias captureSequence = (line) => line
        .replaceAll(sequenceRegex, `$2`)
        .map!toLower;
    alias sequencePrefix = (line) => line
        .replaceAll(sequenceRegex, `$1`);
    alias captureAlignment = (seqLine, alignmentLine) =>
        alignmentLine[sequencePrefix(seqLine).length .. $].tr(" ", "-");

    alignment.alignmentString = chain(
        captureSequence(alignmentLines[0]),
        newline,
        captureAlignment(alignmentLines[0], alignmentLines[1]),
        newline,
        captureSequence(alignmentLines[2]),
    ).map!(to!(immutable(char))).to!string;
}

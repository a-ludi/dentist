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
import dentist.common.external : ExternalDependency;
import dentist.dazzler :
    ContigSegment,
    DBdumpOptions,
    GapSegment,
    getContigCutoff,
    getFastaSequence,
    getDbRecords,
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
import dentist.util.suffixtree : SuffixTree;
import std.algorithm :
    all,
    among,
    copy,
    countUntil,
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
import std.array :
    array,
    minimallyInitializedArray,
    split;
import std.ascii :
    newline,
    toLower;
import std.conv : to;
import std.file :
    exists,
    getSize;
import std.format :
    format,
    formattedRead;
import std.math :
    ceil,
    log10;
import std.parallelism : parallel;
import std.path :
    baseName,
    buildPath,
    setExtension;
import std.process :
    escapeShellFileName,
    executeShell;
import std.range :
    assumeSorted,
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
    splitLines,
    tr;
import std.traits : Unqual;
import std.typecons :
    No,
    tuple,
    Tuple;
import vibe.data.json :
    Json,
    toJson = serializeToJson,
    toJsonCompressed = serializeToJsonString,
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

alias GenomeIndex = SuffixTree!"acgt";

struct ContigMapping
{
    ReferenceInterval reference;
    coord_t referenceContigLength;
    id_t queryContigId;
}

static alias queryOrder = orderLexicographically!(const(ContigMapping),
    mapping => mapping.queryContigId,
);
static alias queryEquiv = (const(ContigMapping) a, const(ContigMapping) b) =>
        a.queryContigId == b.queryContigId;


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

    static enum GapState
    {
        /// Unsure what happened; probably due to an assembly error.
        unkown,
        /// Flanking contigs could be found but not in the expected configuration.
        broken,
        /// Flanking contigs could be found and no sequence was inserted.
        unclosed,
        /// Flanking contigs could be found and sequence was inserted but a gap remains.
        partiallyClosed,
        /// Flanking contigs could be found and contiguous sequence was inserted.
        closed,
        /// This is not a gap – ignore it.
        ignored,
    }

    static struct GapSummary
    {
        id_t lhsContigId;
        GapState state;
        coord_t gapLength;
        StretcherAlignment alignment;

        @property id_t rhsContigId() const
        {
            return lhsContigId + 1;
        }
    }

    const(Options) options;
    protected const(ScaffoldSegment)[] trueAssemblyScaffoldStructure;
    protected const(ScaffoldSegment)[] resultScaffoldStructure;
    protected coord_t referenceOffset;
    protected ReferenceRegion mappedRegionsMask;
    protected ReferenceRegion referenceGaps;
    protected ContigMapping[] contigAlignments;
    protected GapSummary[] gapSummaries;
    protected size_t[][identityLevels.length] correctGapsPerIdentityLevel;

    Stats collect()
    {
        mixin(traceExecution);

        init();

        Stats stats;

        stats.numBpsExpected = getNumBpsExpected();
        stats.numBpsKnown = getNumBpsKnown();
        stats.numBpsResult = getNumBpsResult();
        stats.numBpsInGaps = getNumBpsInGaps();
        stats.numDiffsInGaps = getTotalDiffsInClosedGaps();
        stats.numTranslocatedGaps = getNumTranslocatedGaps();
        stats.numContigsExpected = getNumContigsExpected();
        stats.numMappedContigs = getNumMappedContigs();
        stats.numCorrectGaps = getNumCorrectGaps();
        stats.numClosedGaps = getNumClosedGaps();
        stats.numPartiallyClosedGaps = getNumPartiallyClosedGaps();
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
        referenceGaps = getReferenceGaps();
        gapSummaries = analyzeGaps();
        correctGapsPerIdentityLevel = makeIdentityLevelStats();

        if (options.gapDetailsJson !is null)
            writeGapDetailsJson();

        logJsonDebug(
            "referenceOffset", referenceOffset,
            "numContigAlignments", contigAlignments.length.toJson,
            "contigAlignments", shouldLog(LogLevel.debug_)
                ? contigAlignments.toJson
                : toJson(null),
            "referenceGaps", referenceGaps.toJson,
        );
    }

    ContigMapping[] findReferenceContigs()
    {

        auto duplicateContigIds = NaturalNumberSet.create(
            findPerfectAlignments(options.refDb)
                .map!(selfAlignment => cast(size_t) selfAlignment.reference.contigId)
                .array
        );
        logJsonDiagnostic(
            "numDuplicateContigIds", duplicateContigIds.size,
            "duplicateContigIds", shouldLog(LogLevel.debug_)
                ? duplicateContigIds.elements.array.toJson
                : toJson(null),
        );

        auto perfectContigAlignments = findPerfectAlignments(
            options.resultDb,
            options.refDb,
        )
            // Ignore contigs that have exact copies in `refDb`
            .filter!(contigAlignment => contigAlignment.reference.contigId !in duplicateContigIds)
            .array;

        auto bufferRest = perfectContigAlignments
            .sort!queryOrder
            .group!queryEquiv
            .filter!(group => group[1] == 1)
            .map!(group => group[0])
            .copy(perfectContigAlignments);
        perfectContigAlignments.length -= bufferRest.length;

        return perfectContigAlignments
            .sort!queryOrder
            .release;
    }

    @ExternalDependency("fm-index", "DENTIST", "README.testing")
    ContigMapping[] findPerfectAlignments(in string refDb, in string queryDb = null)
    {
        mixin(traceExecution);

        alias getContigIds = (dbFile) => getDbRecords(dbFile, [DBdumpOptions.readNumber])
            .map!"a.readNumber"
            .array;

        auto isSelfAlignment = queryDb is null;
        auto referenceContigIds = getContigIds(refDb);
        auto refSequenceList = makeIndexedSequenceList(refDb);

        logJsonDiagnostic("indexSize", getSize(refSequenceList ~ ".fm9"));

        auto queryContigIds = isSelfAlignment
            ? referenceContigIds
            : getContigIds(queryDb);
        auto querySequenceList = isSelfAlignment
            ? refSequenceList
            : getSequenceListFile(queryDb);

        assert(querySequenceList.exists, "file should have been created by a prior self-alignment");

        enum resultFieldSeparator = '\t';
        alias FmIndexResult = Tuple!(
            id_t, "refId",
            coord_t, "refLength",
            id_t, "queryId",
            coord_t, "begin",
            coord_t, "end",
        );

        return pipeLines(["fm-index", refSequenceList, querySequenceList])
            .map!(resultLine => resultLine.split(resultFieldSeparator))
            .map!(resultFields => FmIndexResult(
                resultFields[1].to!id_t,
                resultFields[2].to!coord_t,
                resultFields[3].to!id_t,
                resultFields[4].to!coord_t,
                resultFields[5].to!coord_t,
            ))
            .filter!(findResult => !isSelfAlignment || findResult.refId != findResult.queryId)
            .map!(findResult => ContigMapping(
                ReferenceInterval(
                    referenceContigIds[findResult.refId],
                    findResult.begin,
                    findResult.end,
                ),
                findResult.refLength,
                queryContigIds[findResult.queryId],
            ))
            .array;
    }

    @ExternalDependency("DBdump", null, "https://github.com/thegenemyers/DAZZ_DB")
    @ExternalDependency("fm-index", "DENTIST", "README.testing")
    @ExternalDependency("sed", "POSIX")
    string makeIndexedSequenceList(in string dbFile)
    {
        mixin(traceExecution);

        enum commandTemplate = `DBdump -s %s | sed -nE 's/^S [0-9]+ //p' > %s && fm-index %2$s`;

        auto sequenceListFile = getSequenceListFile(dbFile);

        auto commandResult = executeShell(format!commandTemplate(
            escapeShellFileName(dbFile),
            escapeShellFileName(sequenceListFile),
        ));

        dentistEnforce(
            commandResult.status == 0,
            format!"failed to make sequence list `%s`: %s"(sequenceListFile, commandResult.output),
        );

        return sequenceListFile;
    }

    string getSequenceListFile(in string dbFile)
    {
        return buildPath(options.workdir, baseName(dbFile).setExtension(".seq"));
    }

    GapSummary[] analyzeGaps()
    {
        mixin(traceExecution);

        auto contigAlignments = this.contigAlignments.assumeSorted!queryOrder;
        const batchBegin = options.referenceContigBatch[0];
        auto gapSummaries = new GapSummary[options.referenceContigBatchSize];

        ContigMapping needleForContig(in id_t contigId) const
        {
            typeof(return) needleAlignment;

            needleAlignment.queryContigId = contigId;

            return needleAlignment;
        }

        foreach (i, ref gapSummary; parallel(gapSummaries))
        {
            id_t lhsContigId = cast(id_t) (batchBegin + i + 1);
            id_t rhsContigId = lhsContigId + 1;

            gapSummary.lhsContigId = lhsContigId;

            if (!isGap(lhsContigId, rhsContigId))
            {
                gapSummary.state = GapState.ignored;
                continue;
            }

            auto lhsContigAlignments = contigAlignments.equalRange(needleForContig(lhsContigId));
            auto rhsContigAlignments = contigAlignments.equalRange(needleForContig(rhsContigId));

            if (lhsContigAlignments.length != 1 || rhsContigAlignments.length != 1)
            {
                logJsonDiagnostic(
                    "info", "broken gap",
                    "lhs", [
                        "contigId": lhsContigId.toJson,
                        "numContigAlignments": lhsContigAlignments.length.toJson,
                        "contigAlignments": shouldLog(LogLevel.debug_)
                            ? lhsContigAlignments.toJson
                            : null.toJson,
                    ].toJson,
                    "rhs", [
                        "contigId": rhsContigId.toJson,
                        "numContigAlignments": rhsContigAlignments.length.toJson,
                        "contigAlignments": shouldLog(LogLevel.debug_)
                            ? rhsContigAlignments.toJson
                            : null.toJson,
                    ].toJson,
                );
                // One contig alignment either does not exist or is ambiguous:
                // the gap state is unkown (initial value)
                gapSummary.state = GapState.unkown;
                continue;
            }

            auto lhsContigAlignment = lhsContigAlignments[0];
            auto rhsContigAlignment = rhsContigAlignments[0];

            gapSummary.state = getGapState(lhsContigAlignment, rhsContigAlignment);
            gapSummary.gapLength = inputGapSize(lhsContigId);

            if (gapSummary.state.among(GapState.partiallyClosed, GapState.closed))
                gapSummary.alignment = computeInsertionAlignment(getInsertionMapping(
                    lhsContigAlignment,
                    rhsContigAlignment,
                ));
        }

        return gapSummaries;
    }

    private GapState getGapState(in ContigMapping lhs, in ContigMapping rhs) const
    {
        if (!isGap(lhs, rhs))
            return GapState.ignored;
        else if (isGapClosed(lhs, rhs))
            return GapState.closed;
        else if (isGapPartiallyClosed(lhs, rhs))
            return GapState.partiallyClosed;
        else if (isGapUnclosed(lhs, rhs))
            return GapState.unclosed;
        else
            return GapState.broken;
    }

    private bool isGap(in ContigMapping lhs, in ContigMapping rhs) const
    {
        return isGap(lhs.queryContigId, rhs.queryContigId);
    }

    private bool isGap(in id_t lhsContigId, in id_t rhsContigId) const
    {
        auto lhsTrueContigId = mappedIntervalOf(lhsContigId).contigId;
        alias rhsTrueContigId = () => mappedIntervalOf(rhsContigId).contigId;

        return lhsTrueContigId > 0 && lhsTrueContigId == rhsTrueContigId();
    }

    private bool isGapClosed(in ContigMapping lhs, in ContigMapping rhs) const
    {
        return lhs.reference.contigId == rhs.reference.contigId &&
               lhs.reference.end <= rhs.reference.begin;
    }

    private bool isGapPartiallyClosed(in ContigMapping lhs, in ContigMapping rhs) const
    {
        return lhs.reference.contigId + 1 == rhs.reference.contigId &&
               (
                   lhs.reference.end < lhs.referenceContigLength ||
                   0 < rhs.reference.begin
               ) &&
               resultGapSize(cast(id_t) lhs.reference.contigId) > 0;
    }

    private bool isGapUnclosed(in ContigMapping lhs, in ContigMapping rhs) const
    {
        return lhs.reference.contigId + 1 == rhs.reference.contigId &&
               lhs.reference.end == lhs.referenceContigLength &&
               0 == rhs.reference.begin &&
               resultGapSize(cast(id_t) lhs.reference.contigId) > 0;
    }

    InsertionMapping getInsertionMapping(in ContigMapping lhs, in ContigMapping rhs)
    {
        auto trueAssemblyInterval = ReferenceInterval(
            mappedIntervalOf(lhs).contigId,
            mappedIntervalOf(lhs).end,
            mappedIntervalOf(rhs).begin,
        );
        auto resultBegin = ReferencePoint(
            lhs.reference.contigId,
            lhs.reference.end,
        );
        auto resultEnd = ReferencePoint(
            rhs.reference.contigId,
            rhs.reference.begin,
        );

        return InsertionMapping(
            trueAssemblyInterval,
            resultBegin,
            resultEnd,
            cast(id_t) lhs.reference.contigId,
            cast(id_t) rhs.reference.contigId,
        );
    }

    ReferenceInterval mappedIntervalOf(in ContigMapping mapping) const
    {
        return mappedIntervalOf(mapping.queryContigId);
    }

    ReferenceInterval mappedIntervalOf(in id_t contigId) const
    {
        if (contigId > mappedRegionsMask.intervals.length)
            return ReferenceInterval();

        return mappedRegionsMask.intervals[contigId - 1];
    }

    coord_t resultGapSize(in id_t contigId) const
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

    coord_t inputGapSize(in id_t contigId) const
    {
        auto lhsMappedContig = mappedRegionsMask.intervals[contigId - 1];
        auto rhsMappedContig = mappedRegionsMask.intervals[contigId];

        assert(lhsMappedContig.contigId == rhsMappedContig.contigId);
        assert(rhsMappedContig.begin >= lhsMappedContig.end);

        return cast(coord_t) (rhsMappedContig.begin - lhsMappedContig.end);
    }

    size_t[][identityLevels.length] makeIdentityLevelStats()
    {
        typeof(return) correctGapsPerIdentityLevel;
        auto lengthsBuffer = minimallyInitializedArray!(size_t[])(
            getNumClosedGaps() * identityLevels.length,
        );

        foreach (identityLevel, minIdentity; identityLevels)
        {
            auto bufferRest = gapSummaries
                .filter!(gapSummary => gapSummary.state == GapState.closed)
                .filter!(gapSummary => minIdentity <= gapSummary.alignment.percentIdentity)
                .map!(gapSummary => gapSummary.gapLength)
                .copy(lengthsBuffer);

            correctGapsPerIdentityLevel[identityLevel] = lengthsBuffer;
            correctGapsPerIdentityLevel[identityLevel].length -= bufferRest.length;
            lengthsBuffer = bufferRest;
        }

        return correctGapsPerIdentityLevel;
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

        return getGapInfos!"gapLength".sum;
    }

    size_t getTotalDiffsInClosedGaps()
    {
        mixin(traceExecution);

        return getGapInfos!"alignment.numDiffs"(GapState.closed).sum;
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
        return contigAlignments
            .group!queryEquiv
            .filter!(group => group[1] == 1)
            .walkLength;
    }

    size_t getNumCorrectGaps()
    {
        return correctGapsPerIdentityLevel[0].length;
    }

    size_t getNumClosedGaps()
    {
        mixin(traceExecution);

        return gapSummaries
            .filter!(gapSummary => gapSummary.state == GapState.closed)
            .walkLength;
    }

    size_t getNumPartiallyClosedGaps()
    {
        mixin(traceExecution);

        return gapSummaries
            .filter!(gapSummary => gapSummary.state == GapState.partiallyClosed)
            .walkLength;
    }

    void writeGapDetailsJson()
    {
        mixin(traceExecution);

        assert(options.gapDetailsJson !is null);

        auto gapDetailsFile = File(options.gapDetailsJson, "w");

        foreach (gapSummary; gapSummaries)
        {
            if (gapSummary.state != GapState.ignored)
                gapDetailsFile.writeln([
                    "leftContig": gapSummary.lhsContigId.toJson,
                    "rightContig": gapSummary.rhsContigId.toJson,
                    "state": gapSummary.state.to!string.toJson,
                    "gapLength": gapSummary.gapLength.toJson,
                    "details": gapSummary.state.among(GapState.closed, GapState.partiallyClosed)
                        ? (alignment => [
                            "length": alignment.length.toJson,
                            "numIdentical": alignment.numIdentical.toJson,
                            "numDiffs": alignment.numDiffs.toJson,
                            "percentIdentity": alignment.percentIdentity.toJson,
                            "alignment": alignment.alignmentString.splitLines.toJson,
                        ].toJson)(gapSummary.alignment)
                        : toJson(null)
                ].toJsonCompressed);
        }
    }

    StretcherAlignment computeInsertionAlignment(in InsertionMapping insertionMapping)
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

        auto gapSizes = getGapInfos!"gapLength".array;

        if (gapSizes.length == 0)
            return size_t.max;
        else
            return median(gapSizes);
    }

    size_t getClosedGapMedian()
    {
        mixin(traceExecution);

        auto closedGapSizes = getGapInfos!"gapLength"(GapState.closed).array;

        if (closedGapSizes.length == 0)
            return size_t.max;
        else
            return median(closedGapSizes);
    }

    size_t getExtremumClosedGap(alias extremeElement)()
    {
        mixin(traceExecution);

        auto closeGapLengths = getGapInfos!"gapLength"(GapState.closed);

        if (closeGapLengths.empty)
            return size_t.max;
        else
            return extremeElement(closeGapLengths);
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
            getGapInfos!"gapLength"(GapState.closed).array,
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

    auto getGapInfos(string what)(in GapState wantedState) const
    {
        return gapSummaries
            .filter!(gapSummary => gapSummary.state == wantedState)
            .map!(gapSummary => mixin("gapSummary." ~ what))
            .map!(info => cast(Unqual!(typeof(info))) info);
    }

    auto getGapInfos(string what)() const
    {
        return gapSummaries
            .filter!(gapSummary => gapSummary.state != GapState.ignored)
            .map!(gapSummary => mixin("gapSummary." ~ what))
            .map!(info => cast(Unqual!(typeof(info))) info);
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
    static enum identityLevels = [
        1.0,
        0.999,
        0.99,
        0.95,
        0.90,
        0.70,
    ];

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
    size_t numPartiallyClosedGaps;
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
            "numPartiallyClosedGaps": numPartiallyClosedGaps.toJson,
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
            format!"numContigsExpected:     %*d"(columnWidth, numContigsExpected),
            format!"numMappedContigs:       %*d"(columnWidth, numMappedContigs),

            format!"numTranslocatedGaps:    %*d"(columnWidth, numTranslocatedGaps),
            format!"numClosedGaps:          %*d"(columnWidth, numClosedGaps),
            format!"numPartiallyClosedGaps: %*d"(columnWidth, numPartiallyClosedGaps),
            format!"numCorrectGaps:         %*d"(columnWidth, numCorrectGaps),

            format!"numBpsExpected:         %*d"(columnWidth, numBpsExpected),
            format!"numBpsKnown:            %*d"(columnWidth, numBpsKnown),
            format!"numBpsResult:           %*d"(columnWidth, numBpsResult),
            format!"numBpsInGaps:           %*d"(columnWidth, numBpsInGaps),
            format!"numDiffsInGaps:         %*d"(columnWidth, numDiffsInGaps),
            format!"averageInsertionError:  %*.3e"(columnWidth, averageInsertionError),

            format!"maximumN50:             %*d"(columnWidth, maximumN50),
            format!"inputN50:               %*d"(columnWidth, inputN50),
            format!"resultN50:              %*d"(columnWidth, resultN50),
            format!"gapMedian:              %*d"(columnWidth, gapMedian),
            format!"closedGapMedian:        %*d"(columnWidth, closedGapMedian),
            format!"minClosedGap:           %*s"(columnWidth, minClosedGap),
            format!"maxClosedGap:           %*s"(columnWidth, maxClosedGap),
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
            numPartiallyClosedGaps,
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

@ExternalDependency("stretcher", "EMBOSS >=6.0.0", "http://emboss.sourceforge.net/apps/")
StretcherAlignment stretcher(in string refFasta, in string queryFasta)
{
    alias onlyNs = (fasta) => fasta.all!(base => base.among('n', 'N'));

    if (refFasta.getFastaLength() == 0)
        return StretcherAlignment();
    if (queryFasta.getFastaLength() == 0 || onlyNs(queryFasta))
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

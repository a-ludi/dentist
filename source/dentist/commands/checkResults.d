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

import core.atomic : atomicOp;
import dentist.commandline : TestingCommand, OptionsFor;
import dentist.common :
    ReferenceInterval,
    ReferenceRegion,
    to;
import dentist.common.alignments :
    AlignmentChain,
    coord_t,
    id_t;
import dentist.dazzler :
    ContigSegment,
    fingerprint,
    GapSegment,
    getAlignments,
    getExactAlignment,
    getScaffoldStructure,
    readMask,
    ScaffoldSegment;
import dentist.util.algorithm :
    first,
    last,
    sliceBy;
import dentist.util.log;
import dentist.util.math :
    absdiff,
    ceildiv,
    longestIncreasingSubsequence,
    mean,
    median,
    N;
import dentist.util.range : tupleMap;
import std.algorithm :
    all,
    among,
    chunkBy,
    copy,
    count,
    filter,
    find,
    fold,
    joiner,
    map,
    max,
    maxElement,
    min,
    sum;
import std.array : array;
import std.format : format;
import std.math :
    ceil,
    log10;
import std.parallelism : parallel;
import std.range :
    assumeSorted,
    enumerate,
    StoppingPolicy,
    zip;
import std.range.primitives;
import std.stdio : File, writeln;
import std.string : join;
import std.typecons : Tuple;
import vibe.data.json : Json, toJson = serializeToJson, toJsonString = serializeToPrettyJson;


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
    alias ExactAlignment = typeof(getExactAlignment("", "", resultAlignment[0], ""));
    alias Fingerprint = typeof(fingerprint(&resultAlignment[0]));
    static enum identityLevels = Stats.identityLevels;

    const(Options) options;
    protected const(ScaffoldSegment)[] trueAssemblyScaffoldStructure;
    protected const(ScaffoldSegment)[] resultScaffoldStructure;
    protected coord_t referenceOffset;
    protected AlignmentChain[] resultAlignment;
    protected ExactAlignment[Fingerprint] exactResultAlignments;
    protected ReferenceRegion mappedRegionsMask;
    protected ReferenceRegion referenceGaps;
    protected ReferenceRegion reconstructedRegions;
    protected ReferenceRegion[] reconstructedGaps;
    protected size_t[][identityLevels.length] correctContigsPerIdentityLevel;
    protected size_t[][identityLevels.length] correctGapsPerIdentityLevel;

    Stats collect()
    {
        mixin(traceExecution);

        init();

        if (options.alignmentTabular !is null)
            writeAlignmentTabular(
                File(options.alignmentTabular, "w"),
                resultAlignment,
            );

        Stats stats;

        stats.numBpsExpected = getNumBpsExpected();
        stats.numBpsKnown = getNumBpsKnown();
        stats.numBpsResult = getNumBpsResult();
        stats.numTranslocatedGaps = getNumTranslocatedGaps();
        stats.numBpsCorrect = getNumBpsCorrect();
        stats.numContigsExpected = getNumContigsExpected();
        stats.numCorrectContigs = getNumCorrectContigs();
        stats.numCorrectGaps = getNumCorrectGaps();
        stats.numClosedGaps = getNumClosedGaps();
        stats.numPartlyClosedGaps = getNumPartlyClosedGaps();
        stats.numUnaddressedGaps = getNumUnaddressedGaps();
        stats.numBpsInClosedGaps = getNumBpsInClosedGaps();
        stats.numBpsInGaps = getNumBpsInGaps();
        stats.inputN50 = getInputN50();
        stats.resultN50 = getResultN50();
        stats.gapMedian = getGapMedian();
        stats.closedGapMedian = getClosedGapMedian();
        stats.maxClosedGap = getMaxClosedGap();
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
        resultAlignment = getAlignments(
            options.trueAssemblyDb,
            options.resultDb,
            options.resultsAlignmentFile,
            options.workdir,
            options.tracePointDistance,
        );
        calculateExactResultAlignments();
        mappedRegionsMask = ReferenceRegion(readMask!ReferenceInterval(
            options.trueAssemblyDb,
            options.mappedRegionsMask,
            options.workdir,
        ));
        referenceOffset = cast(coord_t) mappedRegionsMask.intervals[0].begin;
        debug logJsonDebug("referenceOffset", referenceOffset);
        referenceGaps = getReferenceGaps();
        reconstructedRegions = getReconstructedRegions();
        reconstructedGaps = getReconstructedGaps();
        correctContigsPerIdentityLevel = getCorrectRegions(mappedRegionsMask.intervals);
        correctGapsPerIdentityLevel = getCorrectRegions(referenceGaps
                .intervals
                .filter!(gap => isInnerGap(gap))
                .array);
    }

    void calculateExactResultAlignments()
    {
        mixin(traceExecution);

        foreach (alignmentChain; parallel(resultAlignment))
        {
            auto acFingerprint = fingerprint(&alignmentChain);
            auto exactAlignment = getExactAlignment(
                options.trueAssemblyDb,
                options.resultDb,
                alignmentChain,
                options.workdir,
                10*2^^30,
            );

            synchronized
                exactResultAlignments[acFingerprint] = exactAlignment;
        }
    }

    ReferenceRegion getReferenceGaps()
    {
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
            .array
        );
    }

    ReferenceRegion getReconstructedRegions()
    {
        return ReferenceRegion(resultAlignment
            .map!(ac => ac
                .localAlignments
                .map!(la => ReferenceInterval(
                    ac.contigA.id,
                    la.contigA.begin,
                    la.contigA.end,
                )))
            .joiner
            .array);
    }

    ReferenceRegion[] getReconstructedGaps()
    {
        return referenceGaps
            .intervals
            .map!(gap => reconstructedRegions & gap)
            .array;
    }

    void writeAlignmentTabular(File tabularFile, AlignmentChain[] alignmentChains)
    {
        coord_t contigBId;
        coord_t contigBOffset;
        foreach (alignmentChain; alignmentChains)
        {
            if (contigBId != alignmentChain.contigB.id)
            {
                contigBId = alignmentChain.contigB.id;
                contigBOffset = getResultContigBegin(contigBId);
            }

            alias percentSimilarity = (la) => 1.0 - cast(double) la.numDiffs /
                                                    cast(double) (la.contigA.end - la.contigA.begin);

            foreach (i, localAlignment; alignmentChain.localAlignments)
            {
                tabularFile.writefln!"%12d %12d %12d %12d %12f %d"(
                    alignmentChain.contigA.id,
                    alignmentChain.contigB.id,
                    localAlignment.contigA.begin,
                    alignmentChain.flags.complement
                        ? alignmentChain.contigB.length - (contigBOffset + localAlignment.contigB.end)
                        : contigBOffset + localAlignment.contigB.begin,
                    0.0,
                    alignmentChain.flags.complement ? -1 : 1,
                );
                tabularFile.writefln!"%12d %12d %12d %12d %12f %d"(
                    alignmentChain.contigA.id,
                    alignmentChain.contigB.id,
                    localAlignment.contigA.end,
                    alignmentChain.flags.complement
                        ? alignmentChain.contigB.length - (contigBOffset + localAlignment.contigB.begin)
                        : contigBOffset + localAlignment.contigB.end,
                    percentSimilarity(localAlignment),
                    alignmentChain.flags.complement ? -1 : 1,
                );
            }
            tabularFile.writeln();
        }
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

    size_t getNumBpsCorrect()
    {
        mixin(traceExecution);

        return resultAlignment
            .map!(ac => ac
                .localAlignments
                .map!(la => la.contigA.end - la.contigA.begin - la.numDiffs))
            .joiner
            .sum;
    }

    size_t getNumTranslocatedGaps()
    {
        mixin(traceExecution);

        version (assert)
        {
            foreach (pair; zip(referenceGaps.intervals, reconstructedGaps))
                assert(
                    isGapClosed(pair) ^ isGapPartlyClosed(pair) ^ isGapUnaddressed(pair),
                    format!"%s-%s-%s"(
                        isGapClosed(pair),
                        isGapPartlyClosed(pair),
                        isGapUnaddressed(pair),
                    ),
                );
        }

        return referenceGaps.intervals.length;
    }

    size_t getNumContigsExpected()
    {
        mixin(traceExecution);

        return mappedRegionsMask.intervals.length;
    }

    size_t getNumCorrectContigs()
    {
        return correctContigsPerIdentityLevel[0].length;
    }

    size_t getNumCorrectGaps()
    {
        return correctGapsPerIdentityLevel[0].length;
    }

    size_t[][identityLevels.length] getCorrectRegions(in ReferenceInterval[] intervals)
    {
        auto sortedResultAlignment = resultAlignment.assumeSorted!"a.contigA.id < b.contigA.id";

        size_t[][identityLevels.length] identicalLengthsPerLevel;
        foreach (identicalLengths; identicalLengthsPerLevel)
            identicalLengths.reserve(intervals.length);

        foreach (interval; parallel(intervals))
        {
            auto overlappingAlignments = sortedResultAlignment
                .equalRange(getDummyAC(interval));
            auto overlappedRegion = overlappingAlignments
                .map!(to!(ReferenceRegion, "contigA"))
                .fold!"a | b"(ReferenceRegion());
            size_t numDiffs = (ReferenceRegion(interval) - overlappedRegion).size;

            foreach (overlappingAlignment; overlappingAlignments)
            {
                auto alignmentBegin = overlappingAlignment.first.contigA.begin;
                auto alignmentEnd = overlappingAlignment.last.contigA.end;
                if (interval.end <= alignmentBegin || alignmentEnd <= interval.begin)
                    continue; // Ignore non-overlaps
                auto overlapBegin = interval.begin > alignmentBegin
                    ? interval.begin - alignmentBegin
                    : 0;
                auto overlapEnd = interval.end < alignmentEnd
                    ? interval.end - alignmentBegin
                    : alignmentEnd - alignmentBegin;

                if (overlappingAlignment.flags.complement)
                {
                    // Inversions are not allowed
                    numDiffs += overlapEnd - overlapBegin;
                }
                else
                {
                    auto exactAlignment = getOverlapAlignment(overlappingAlignment);
                    auto overlappingExactAlignment = exactAlignment.partial(overlapBegin, overlapEnd);

                    numDiffs += overlappingExactAlignment.score;
                }
            }

            foreach (i, identityLevel; identityLevels)
                if (numDiffs <= (1.0 - identityLevel) * interval.size)
                    synchronized identicalLengthsPerLevel[i] ~= interval.size;
        }

        return identicalLengthsPerLevel;
    }

    auto getOverlapAlignment(in AlignmentChain alignmentChain)
    {
        auto acFingerprint = fingerprint(&alignmentChain);

        return exactResultAlignments[acFingerprint];
    }

    static AlignmentChain getDummyAC(in ReferenceInterval refInterval)
    {
        return AlignmentChain(
            0,
            AlignmentChain.Contig(
                cast(id_t) refInterval.contigId,
                cast(coord_t) refInterval.end,
            ),
            AlignmentChain.Contig(0, 1),
            AlignmentChain.emptyFlags,
            [AlignmentChain.LocalAlignment(
                AlignmentChain.LocalAlignment.Locus(
                    cast(coord_t) refInterval.begin,
                    cast(coord_t) refInterval.end,
                ),
                AlignmentChain.LocalAlignment.Locus(
                    cast(coord_t) 0,
                    cast(coord_t) 1,
                ),
            )],
        );
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

    size_t getNumClosedGaps()
    {
        mixin(traceExecution);

        assert(referenceGaps.intervals.length == reconstructedGaps.length);
        logJsonDebug("closedGaps", zip(referenceGaps.intervals, reconstructedGaps)
            .filter!(pair => isGapClosed(pair))
            .map!"a[1].intervals"
            .array
            .toJson);

        return zip(referenceGaps.intervals, reconstructedGaps).count!(pair => isGapClosed(pair));
    }

    size_t getNumPartlyClosedGaps()
    {
        mixin(traceExecution);

        assert(referenceGaps.intervals.length == reconstructedGaps.length);
        logJsonDebug("partlyClosedGaps", zip(referenceGaps.intervals, reconstructedGaps)
            .filter!(pair => isGapPartlyClosed(pair))
            .map!"a[1].intervals"
            .array
            .toJson);

        return zip(referenceGaps.intervals, reconstructedGaps).count!(pair => isGapPartlyClosed(pair));
    }

    size_t getNumUnaddressedGaps()
    {
        mixin(traceExecution);

        assert(referenceGaps.intervals.length == reconstructedGaps.length);
        logJsonDebug("unaddressedGaps", zip(referenceGaps.intervals, reconstructedGaps)
            .filter!(pair => isGapUnaddressed(pair))
            .map!"a[1].intervals"
            .array
            .toJson);

        return zip(referenceGaps.intervals, reconstructedGaps).count!(pair => isGapUnaddressed(pair));
    }

    bool isGapClosed(Pair)(Pair gapAndInsertion)
    {
        auto gap = gapAndInsertion[0];
        auto insertion = gapAndInsertion[1];

        return insertion.size == gap.size;
    }

    bool isGapPartlyClosed(Pair)(Pair gapAndInsertion)
    {
        auto gap = gapAndInsertion[0];
        auto insertion = gapAndInsertion[1];

        return options.minInsertionLength <= insertion.size && insertion.size < gap.size;
    }

    bool isGapUnaddressed(Pair)(Pair gapAndInsertion)
    {
        auto gap = gapAndInsertion[0];
        auto insertion = gapAndInsertion[1];

        return insertion.size < min(options.minInsertionLength, gap.size);
    }

    size_t getNumBpsInClosedGaps()
    {
        mixin(traceExecution);

        return zip(referenceGaps.intervals, reconstructedGaps)
            .filter!(pair => isGapClosed(pair))
            .map!"a[1].size"
            .sum;
    }

    size_t getNumBpsInGaps()
    {
        mixin(traceExecution);

        return referenceGaps
            .intervals
            .filter!(gap => isInnerGap(gap))
            .map!(gap => gap.size)
            .sum;
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

        auto values = zip(referenceGaps.intervals, reconstructedGaps)
            .filter!(pair => isGapClosed(pair))
            .map!"a[0].size"
            .array;

        if (values.length == 0)
            return size_t.max;
        else
            return median(values);
    }

    ReferenceInterval getMaxClosedGap()
    {
        mixin(traceExecution);

        if (reconstructedGaps.length == 0)
            return ReferenceInterval(0, 0, 0);

        return zip(referenceGaps.intervals, reconstructedGaps)
            .filter!(pair => isGapClosed(pair))
            .map!"a[0]"
            .maxElement!"a.size";
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
            zip(referenceGaps.intervals, reconstructedGaps)
                .filter!(pair => isGapClosed(pair))
                .map!(pair => cast(coord_t) pair[0].size)
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
    size_t numBpsCorrect;
    size_t numTranslocatedGaps;
    size_t numCorrectGaps;
    size_t numContigsExpected;
    size_t numCorrectContigs;
    size_t numClosedGaps;
    size_t numPartlyClosedGaps;
    size_t numUnaddressedGaps;
    size_t numBpsInClosedGaps;
    size_t numBpsInGaps;
    size_t inputN50;
    size_t resultN50;
    size_t gapMedian;
    size_t closedGapMedian;
    ReferenceInterval maxClosedGap;
    Histogram!coord_t[identityLevels.length] correctGapLengthHistograms;
    Histogram!coord_t closedGapLengthHistogram;
    Histogram!coord_t gapLengthHistogram;

    string toJsonString() const
    {
        return [
            "numBpsExpected": numBpsExpected.toJson,
            "numBpsKnown": numBpsKnown.toJson,
            "numBpsResult": numBpsResult.toJson,
            "numBpsCorrect": numBpsCorrect.toJson,
            "numTranslocatedGaps": numTranslocatedGaps.toJson,
            "numCorrectGaps": numCorrectGaps.toJson,
            "numContigsExpected": numContigsExpected.toJson,
            "numCorrectContigs": numCorrectContigs.toJson,
            "numClosedGaps": numClosedGaps.toJson,
            "numPartlyClosedGaps": numPartlyClosedGaps.toJson,
            "numUnaddressedGaps": numUnaddressedGaps.toJson,
            "numBpsInClosedGaps": numBpsInClosedGaps.toJson,
            "numBpsInGaps": numBpsInGaps.toJson,
            "inputN50": inputN50.toJson,
            "resultN50": resultN50.toJson,
            "gapMedian": gapMedian.toJson,
            "closedGapMedian": closedGapMedian.toJson,
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
            format!"numCorrectContigs:     %*d"(columnWidth, numCorrectContigs),

            format!"numTranslocatedGaps:   %*d"(columnWidth, numTranslocatedGaps),
            format!"numClosedGaps:         %*d"(columnWidth, numClosedGaps),
            format!"numPartlyClosedGaps:   %*d"(columnWidth, numPartlyClosedGaps),
            format!"numUnaddressedGaps:    %*d"(columnWidth, numUnaddressedGaps),
            format!"numCorrectGaps:        %*d"(columnWidth, numCorrectGaps),

            format!"numBpsExpected:        %*d"(columnWidth, numBpsExpected),
            format!"numBpsKnown:           %*d"(columnWidth, numBpsKnown),
            format!"numBpsResult:          %*d"(columnWidth, numBpsResult),
            format!"numBpsCorrect:         %*d"(columnWidth, numBpsCorrect),
            format!"numBpsInGaps:          %*d"(columnWidth, numBpsInGaps),
            format!"numBpsInClosedGaps:    %*d"(columnWidth, numBpsInClosedGaps),
            format!"numBpsInNonClosedGaps: %*d"(columnWidth, numBpsInGaps - numBpsInClosedGaps),

            format!"inputN50:              %*d"(columnWidth, inputN50),
            format!"resultN50:             %*d"(columnWidth, resultN50),
            format!"gapMedian:             %*d"(columnWidth, gapMedian),
            format!"closedGapMedian:       %*d"(columnWidth, closedGapMedian),
            format!"maxClosedGap:          %*s"(columnWidth, toString(maxClosedGap)),
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

    const string toString(in ReferenceInterval interval)
    {
        return format!"[%d, %d) | len=%d"(interval.begin, interval.end, interval.size);
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
            numTranslocatedGaps,
            numBpsCorrect,
            numCorrectGaps,
            numContigsExpected,
            numCorrectContigs,
            numClosedGaps,
            numPartlyClosedGaps,
            numUnaddressedGaps,
            numBpsInClosedGaps,
            numBpsInGaps,
            inputN50,
            resultN50,
            gapMedian,
            closedGapMedian,
        ));
        auto strWidth = max(
            0, // dummy
            toString(maxClosedGap).length,
        );

        return max(numWidth, strWidth);
    }

    static protected size_t numWidth(Int)(Int i) nothrow
    {
        auto numWidth = cast(size_t) ceil(log10(i));

        return numWidth;
    }
}

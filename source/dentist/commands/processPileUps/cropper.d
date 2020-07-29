/**
    This package contains methods for cropping a pile up.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.processPileUps.cropper;

import dentist.common :
    ReadInterval,
    ReferenceInterval,
    ReferenceRegion,
    ReferencePoint,
    to;
import dentist.common.alignments :
    AlignmentChain,
    AlignmentLocationSeed,
    coord_t,
    getType,
    isExtension,
    PileUp,
    ReadAlignment,
    SeededAlignment,
    trace_point_t;
import dentist.dazzler : buildDbFile, getFastaSequences;
import dentist.util.algorithm : sliceBy;
import dentist.util.fasta : reverseComplement;
import dentist.util.log;
import dentist.util.math : ceil, ceildiv, RoundingMode;
import dentist.util.region : min, sup;
import std.algorithm :
    all,
    canFind,
    copy,
    countUntil,
    joiner,
    filter,
    fold,
    map,
    predSwitch,
    sort,
    sum,
    swap;
import std.array : array, minimallyInitializedArray;
import std.conv : to;
import std.format : format;
import std.path : buildPath;
import std.range :
    chain,
    chunks,
    iota,
    only,
    retro,
    takeExactly,
    zip;
import std.typecons : tuple;
import vibe.data.json : toJson = serializeToJson;


struct CropOptions
{
    string refDb;
    string readsDb;
    coord_t minAnchorLength;
    string tmpdir;
}

auto cropPileUp(PileUp pileUp, in ReferenceRegion mask, in CropOptions options)
{
    auto cropper = PileUpCropper(pileUp, mask, options);
    cropper.buildDb();

    return cropper.result;
}

private struct PileUpCropper
{
    PileUp pileUp;
    const(ReferenceRegion) repeatMask;
    const(CropOptions) options;
    private ReferencePoint[] croppingRefPositions;
    private AlignmentLocationSeed[] croppingSeeds;
    private string[] supportPatches;
    private string[] supportPatchesRevComp;
    private string croppedDb;

    @property auto result()
    {
        return tuple!(
            "db",
            "referencePositions",
            "seeds",
        )(
            croppedDb,
            croppingRefPositions,
            croppingSeeds,
        );
    }

    void buildDb()
    {
        enum extensionDbName = "pileup-%d.db";
        enum gapDbName = "pileup-%d-%d.db";

        fetchCroppingRefPositions();
        fetchSupportPatches();
        auto croppedSequences = pileUpWithSequence().map!(t => getCroppedSequence(t.expand));

        if (croppingRefPositions.length == 1)
            croppedDb = format!extensionDbName(croppingRefPositions[0].contigId);
        else
            croppedDb = format!gapDbName(
                croppingRefPositions[0].contigId,
                croppingRefPositions[1].contigId,
            );
        croppedDb = buildPath(options.tmpdir, croppedDb);

        buildDbFile(
            croppedDb,
            croppedSequences,
        );
    }

    private void fetchCroppingRefPositions()
    {
        auto result = getCroppingRefPositions();
        croppingRefPositions = result.refPositions;
        croppingSeeds = result.seeds;
        logJsonDebug("croppingRefPositions", croppingRefPositions.toJson);
        logJsonDebug("croppingSeeds", croppingSeeds.toJson);

        foreach (refPos; croppingRefPositions)
        {
            if (refPos.value == -1)
            {
                if (shouldLog(LogLevel.diagnostic))
                {
                    auto involvedContigs = croppingRefPositions.map!"a.contigId".array;
                    auto localRepeatMask = repeatMask
                        .intervals
                        .filter!(i => involvedContigs.canFind(i.contigId))
                        .array;
                    logJsonDiagnostic(
                        "info", "could not find a common trace point",
                        "croppingRefPositions", croppingRefPositions.toJson,
                        "croppingSeeds", croppingSeeds.toJson,
                        "tracePointDistance", pileUp.length > 0
                            ? pileUp[0][0].tracePointDistance
                            : 0,
                        "pileUp", [
                            "type": pileUp.getType.to!string.toJson,
                            "readAlignments": shouldLog(LogLevel.debug_)
                                ? pileUp.map!"a[]".array.toJson
                                : toJson(null),
                        ].toJson,
                        "repeatMask", localRepeatMask.toJson,
                    );
                }

                throw new Exception("could not find a common trace point");
            }
        }
    }

    /// Fetch pieces of the flanking contig(s) that can be appended to the
    /// cropped reads in case the sequence remaining after cropping is
    /// insufficient for later alignments, e.i. shorter than
    /// `options.minAnchorLength`.
    private void fetchSupportPatches()
    {
        auto contigSequences = getFastaSequences(
            options.refDb,
            croppingRefPositions.map!"a.contigId",
        ).map!idup.array;

        foreach (i, ref contigSequence; contigSequences)
        {
            auto seed = croppingSeeds[i];
            auto contigId = croppingRefPositions[i].contigId;
            auto pos = croppingRefPositions[i].value;
            auto patchInterval = ReferenceInterval(contigId);

            final switch(seed)
            {
                case AlignmentLocationSeed.front:
                    if (pos < options.minAnchorLength)
                    {
                        patchInterval.begin = pos;
                        patchInterval.end = options.minAnchorLength;
                    }
                    break;
                case AlignmentLocationSeed.back:
                    assert(pos <= contigSequence.length);
                    if (contigSequence.length - pos < options.minAnchorLength)
                    {
                        patchInterval.begin = contigSequence.length - options.minAnchorLength;
                        patchInterval.end = pos;
                    }
                    break;
            }

            contigSequence = contigSequence[patchInterval.begin .. patchInterval.end];
        }

        supportPatches = contigSequences;
        supportPatchesRevComp = contigSequences.map!reverseComplement.array;
    }

    private auto pileUpWithSequence()
    {
        auto readIds = pileUp.map!"a[0].contigB.id + 0".array;

        return zip(
            iota(pileUp.length),
            pileUp,
            getFastaSequences(options.readsDb, readIds),
        );
    }

    /**
        Get points on the reference where the pileUp should be cropped. Returns
        one common trace point for each involved contig.

        See_Also: `getCommonTracePoint`
    */
    private auto getCroppingRefPositions()
    {
        auto alignmentsByContig = pileUp.splitAlignmentsByContigA();
        auto commonTracePoints = alignmentsByContig
            .map!(acs => getCommonTracePoint(acs, repeatMask));
        auto contigAIds = alignmentsByContig.map!"a[0].contigA.id";
        auto seeds = alignmentsByContig.map!"a[0].seed".array;

        auto cropPos = zip(contigAIds, commonTracePoints)
            .map!(zipped => ReferencePoint(zipped.expand))
            .array;

        debug logJsonDebug(
            "alignmentsByContig", alignmentsByContig.toJson,
            "pileUp", [
                "type": pileUp.getType.to!string.toJson,
                "readAlignments": pileUp.map!"a[]".array.toJson,
            ],
            "cropPos", cropPos.toJson,
            "seeds", seeds.toJson,
        );

        return tuple!(
            "refPositions",
            "seeds",
        )(
            cropPos,
            seeds,
        );
    }

    private string getCroppedSequence(
        in size_t croppedDbIdx,
        in ReadAlignment readAlignment,
        in string readSequence,
    )
    {
        auto readCroppingSlice = getReadCroppingSlice(readAlignment);
        auto readPatches = getReadPatches(readAlignment);
        debug logJsonDebug(
            "croppedDbIdx", croppedDbIdx,
            "readCroppingSlice", readCroppingSlice.toJson,
        );

        return getCroppedReadAsFasta(
            readSequence,
            readCroppingSlice,
            readPatches.pre,
            readPatches.post,
        );
    }

    private auto getReadCroppingSlice(in ReadAlignment readAlignment)
    {
        auto readCroppingSlice = readAlignment[]
            .map!(alignment => alignment.getCroppingSlice(croppingRefPositions))
            .fold!"a & b";
        assert(!readCroppingSlice.empty, "invalid/empty read cropping slice");

        return readCroppingSlice;
    }

    private auto getReadPatches(in ReadAlignment readAlignment)
    {
        auto readPatches = readAlignment[]
            .map!(alignment => getSingleReadPatch(alignment))
            .array
            .sort
            .release;

        if (readPatches.length == 2)
            return tuple!("pre", "post")(readPatches[0].sequence, readPatches[1].sequence);
        else if (readPatches[0].readSeed == AlignmentLocationSeed.front)
            return tuple!("pre", "post")(readPatches[0].sequence, "");
        else
            return tuple!("pre", "post")("", readPatches[0].sequence);
    }

    auto getSingleReadPatch(in SeededAlignment alignment)
    {
        auto contigA = alignment.contigA.id;
        auto contigSeed = alignment.seed;
        auto complement = alignment.flags.complement;
        auto readSeed = (contigSeed == AlignmentLocationSeed.front) ^ complement
            ? AlignmentLocationSeed.back
            : AlignmentLocationSeed.front;
        auto cropIdx = croppingRefPositions.countUntil!(refPos => refPos.contigId == contigA);
        auto sequence = complement
            ? supportPatchesRevComp[cropIdx]
            : supportPatches[cropIdx];

        return tuple!("readSeed", "sequence")(readSeed, sequence);
    }

    private string getCroppedReadAsFasta(
        string readSequence,
        in ReadInterval readCroppingSlice,
        string prePatch,
        string postPatch,
    )
    {
        static enum fastaLineWidth = 100;
        auto fastaHeader = format!">read-%d/%d/%d_%d"(
            readCroppingSlice.readId,
            readCroppingSlice.readId,
            0,
            prePatch.length + readCroppingSlice.size + postPatch.length,
        );
        auto patchedSequence = chain(
            prePatch,
            readSequence[readCroppingSlice.begin .. readCroppingSlice.end],
            postPatch,
        );
        auto patchedSequenceLength = prePatch.length +
                                     readCroppingSlice.size +
                                     postPatch.length;
        auto fastaLength =
            fastaHeader.length + 1 +                      // header line + line break
            patchedSequenceLength +                       // sequence + ...
            ceildiv(patchedSequenceLength, fastaLineWidth);  // line breaks for sequence +
                                                          // new line at the end of file
        auto croppedFasta = minimallyInitializedArray!(ubyte[])(fastaLength);
        chain(
            fastaHeader[],
            "\n",
            patchedSequence
                .chunks(fastaLineWidth)
                .joiner("\n"),
            "\n",
        ).copy(croppedFasta);

        return cast(string) croppedFasta;
    }
}

/// Sort alignment chains of pileUp into groups with the same `contigA.id`.
private SeededAlignment[][] splitAlignmentsByContigA(PileUp pileUp)
{
    if (pileUp.length == 0)
    {
        return [];
    }

    auto orderedAlignments = pileUp
        .map!"a[]"
        .joiner
        .array
        .sort
        .release;
    auto alignmentsByContig = orderedAlignments.sliceBy!"a.contigA.id == b.contigA.id";

    return array(alignmentsByContig);
}

/// Returns a common trace points wrt. contigA that is not in mask.
private long getCommonTracePoint(
    in SeededAlignment[] alignments,
    in ReferenceRegion mask,
)
{
    static long _getCommonTracePoint(R)(R tracePointCandidates, ReferenceRegion tracePointRegion) pure
    {
        auto commonTracePoints = tracePointCandidates
            .filter!(c => c in tracePointRegion || c.value == tracePointRegion.intervals[$ - 1].end);

        return commonTracePoints.empty ? -1 : commonTracePoints.front.value.to!long;
    }

    auto contigA = alignments[0].contigA;
    auto locationSeed = alignments[0].seed;
    auto tracePointDistance = alignments[0].tracePointDistance;
    auto commonAlignmentRegion = alignments
        .map!(to!(ReferenceRegion, "contigA"))
        .fold!"a & b";
    auto unmaskedTracePointRegion = commonAlignmentRegion - mask;

    assert(alignments.all!(a => a.contigA == contigA && a.seed == locationSeed));
    debug logJsonDebug(
        "alignments", alignments.toJson,
        "commonAlignmentRegion", commonAlignmentRegion.toJson,
        "unmaskedTracePointRegion", unmaskedTracePointRegion.toJson,
    );

    foreach (ReferenceRegion tracePointRegion; [
        unmaskedTracePointRegion,
        commonAlignmentRegion,
    ])
    {
        if (tracePointRegion.empty)
            continue;

        auto tracePointMin = ceil(min(tracePointRegion), tracePointDistance);
        auto tracePointSup = ceil(sup(tracePointRegion), tracePointDistance);
        auto contigEndTracePoint = tracePointSup > contigA.length
            ? [contigA.length]
            : [];
        auto tracePointCandidates = chain(
            iota(tracePointMin, tracePointSup, tracePointDistance),
            contigEndTracePoint,
        ).map!(tracePointCandidate => ReferencePoint(contigA.id, tracePointCandidate));

        auto commonTracePoint = locationSeed == AlignmentLocationSeed.front
            ? _getCommonTracePoint(tracePointCandidates.retro, tracePointRegion)
            : _getCommonTracePoint(tracePointCandidates, tracePointRegion);

        if (commonTracePoint < 0)
            continue;

        return commonTracePoint;
    }

    return -1;
}

private ReadInterval getCroppingSlice(
    in SeededAlignment alignment,
    in ReferencePoint[] croppingRefPoints,
)
{
    auto locationSeed = alignment.seed;
    auto read = alignment.contigB;
    auto tracePointDistance = alignment.tracePointDistance;
    auto croppingRefPos = croppingRefPoints
        .filter!(p => p.contigId == alignment.contigA.id)
        .front
        .value;

    auto readCroppingPos = alignment
        .translateTracePoint(cast(coord_t) croppingRefPos, RoundingMode.floor)
        .contigB;
    size_t readBeginIdx;
    size_t readEndIdx;

    final switch (locationSeed)
    {
    case AlignmentLocationSeed.front:
        readBeginIdx = 0;
        readEndIdx = readCroppingPos;
        break;
    case AlignmentLocationSeed.back:
        readBeginIdx = readCroppingPos;
        readEndIdx = read.length;
        break;
    }

    if (alignment.complement)
    {
        swap(readBeginIdx, readEndIdx);
        readBeginIdx = read.length - readBeginIdx;
        readEndIdx = read.length - readEndIdx;
    }

    auto croppingSlice = ReadInterval(read.id, readBeginIdx, readEndIdx);

    debug logJsonDebug(
        "alignment", alignment.toJson,
        "croppingRefPoints", croppingRefPoints.toJson,
        "croppingSlice", croppingSlice.toJson,
    );

    return croppingSlice;
}

unittest
{
    alias Contig = AlignmentChain.Contig;
    alias Flags = AlignmentChain.Flags;
    alias emptyFlags = AlignmentChain.emptyFlags;
    alias complement = AlignmentChain.Flag.complement;
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias Locus = LocalAlignment.Locus;
    alias TracePoint = LocalAlignment.TracePoint;
    enum tracePointDistance = 100;

    {
        enum alignment = SeededAlignment(
            AlignmentChain(
                0,
                Contig(1, 2584),
                Contig(58024, 10570),
                Flags(complement),
                [LocalAlignment(
                    Locus(579, 2584),
                    Locus(0, 2158),
                    292,
                    [
                        TracePoint( 2,  23),
                        TracePoint(11, 109),
                        TracePoint(13, 109),
                        TracePoint(15, 107),
                        TracePoint(18, 107),
                        TracePoint(14, 103),
                        TracePoint(16, 106),
                        TracePoint(17, 106),
                        TracePoint( 9, 106),
                        TracePoint(14, 112),
                        TracePoint(16, 105),
                        TracePoint(16, 114),
                        TracePoint(10, 103),
                        TracePoint(14, 110),
                        TracePoint(15, 110),
                        TracePoint(15, 101),
                        TracePoint(17, 108),
                        TracePoint(17, 109),
                        TracePoint(15, 111),
                        TracePoint(17, 111),
                        TracePoint(11,  88),
                    ],
                )],
                tracePointDistance,
            ),
            AlignmentLocationSeed.back,
        );

        assert(
            getCroppingSlice(alignment, [ReferencePoint(1, 600)]) ==
            ReadInterval(58024, 0, 10547)
        );
        assert(
            getCroppingSlice(alignment, [ReferencePoint(1, 800)]) ==
            ReadInterval(58024, 0, 10329)
        );
    }
    {
        enum alignment = SeededAlignment(
            AlignmentChain(
                0,
                Contig(1, 2584),
                Contig(7194, 9366),
                emptyFlags,
                [LocalAlignment(
                    Locus(0, 724),
                    Locus(8612, 9366),
                    76,
                    [
                        TracePoint( 7, 107),
                        TracePoint(17, 106),
                        TracePoint(12,  96),
                        TracePoint( 9, 107),
                        TracePoint( 7, 102),
                        TracePoint(11, 105),
                        TracePoint(11, 105),
                        TracePoint( 2,  26),
                    ],
                )],
                tracePointDistance,
            ),
            AlignmentLocationSeed.front,
        );

        assert(
            getCroppingSlice(alignment, [ReferencePoint(1, 0)]) ==
            ReadInterval(7194, 0, 8612)
        );
        assert(
            getCroppingSlice(alignment, [ReferencePoint(1, 100)]) ==
            ReadInterval(7194, 0, 8719)
        );
        assert(
            getCroppingSlice(alignment, [ReferencePoint(1, 200)]) ==
            ReadInterval(7194, 0, 8825)
        );
    }
}

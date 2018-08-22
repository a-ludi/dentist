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
    ReferenceRegion,
    ReferencePoint,
    to;
import dentist.common.alignments :
    AlignmentChain,
    AlignmentLocationSeed,
    getType,
    PileUp,
    ReadAlignment,
    SeededAlignment,
    trace_point_t;
import dentist.dazzler : buildDamFile, getFastaSequences;
import dentist.util.algorithm : sliceBy;
import dentist.util.log;
import dentist.util.math : ceil, ceildiv;
import dentist.util.region : min, sup;
import std.algorithm :
    all,
    canFind,
    copy,
    joiner,
    filter,
    fold,
    map,
    sort,
    sum,
    swap;
import std.array : array, minimallyInitializedArray;
import std.conv : to;
import std.format : format;
import std.range :
    chain,
    chunks,
    iota,
    retro,
    zip;
import std.typecons : tuple;
import vibe.data.json : toJson = serializeToJson;


struct CropOptions
{
    string readsDb;
    trace_point_t tracePointDistance;
    string workdir;
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
    private string croppedDb;

    @property auto result()
    {
        return tuple!("db", "referencePositions")(croppedDb, croppingRefPositions);
    }

    void buildDb()
    {
        fetchCroppingRefPositions();
        auto croppedSequences = pileUpWithSequence().map!(t => getCroppedSequence(t.expand));
        croppedDb = buildDamFile(croppedSequences, options.workdir);
    }

    private void fetchCroppingRefPositions()
    {
        croppingRefPositions = getCroppingRefPositions();
        logJsonDebug("croppingRefPositions", croppingRefPositions.toJson);

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
                        "tracePointDistance", options.tracePointDistance,
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

    private auto pileUpWithSequence()
    {
        auto readIds = pileUp.map!"a[0].contigB.id + 0".array;

        return zip(
            iota(pileUp.length),
            pileUp,
            getFastaSequences(options.readsDb, readIds, options.workdir),
        );
    }

    /**
        Get points on the reference where the pileUp should be cropped. Returns
        one common trace point for each involved contig.

        See_Also: `getCommonTracePoint`
    */
    private ReferencePoint[] getCroppingRefPositions()
    {
        auto alignmentsByContig = pileUp.splitAlignmentsByContigA();
        auto commonTracePoints = alignmentsByContig
            .map!(acs => getCommonTracePoint(acs, repeatMask, options.tracePointDistance));
        auto contigAIds = alignmentsByContig.map!"a[0].contigA.id";

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
        );

        return cropPos;
    }

    private string getCroppedSequence(
        in size_t croppedDbIdx,
        in ReadAlignment readAlignment,
        in string readSequence,
    )
    {
        auto readCroppingSlice = getReadCroppingSlice(readAlignment);
        debug logJsonDebug(
            "croppedDbIdx", croppedDbIdx,
            "readCroppingSlice", readCroppingSlice.toJson,
        );

        return getCroppedReadAsFasta(
            readSequence,
            readCroppingSlice,
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

    private string getCroppedReadAsFasta(string readSequence, in ReadInterval readCroppingSlice)
    {
        static enum fastaLineWidth = 100;
        auto fastaHeader = format!">read-%d"(readCroppingSlice.readId);
        auto fastaLength =
            fastaHeader.length + 1 +                          // header line + line break
            readCroppingSlice.size +                          // sequence + ...
            ceildiv(readCroppingSlice.size, fastaLineWidth);  // line breaks for sequence +
                                                              // new line at the end of file
        auto croppedFasta = minimallyInitializedArray!(ubyte[])(fastaLength);
        chain(
            fastaHeader[],
            "\n",
            readSequence[readCroppingSlice.begin .. readCroppingSlice.end]
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

/// Returns a common trace points wrt. contigA that is not in mask and is on
/// the relevant half of the contig.
private long getCommonTracePoint(
    in SeededAlignment[] alignments,
    in ReferenceRegion mask,
    in size_t tracePointDistance
)
{
    static long getCommonTracePoint(R)(R tracePointCandidates, ReferenceRegion allowedTracePointRegion) pure
    {
        auto commonTracePoints = tracePointCandidates.filter!(c => c in allowedTracePointRegion);

        return commonTracePoints.empty ? -1 : commonTracePoints.front.value.to!long;
    }

    auto contigA = alignments[0].contigA;
    auto locationSeed = alignments[0].seed;
    auto commonAlignmentRegion = alignments
        .map!(to!(ReferenceRegion, "contigA"))
        .fold!"a & b";
    auto relevantContigHalf = locationSeed == AlignmentLocationSeed.front
        ? ReferenceRegion(contigA.id, 0, contigA.length / 2)
        : ReferenceRegion(contigA.id, contigA.length / 2, contigA.length);
    auto allowedTracePointRegion = (relevantContigHalf & commonAlignmentRegion) - mask;

    assert(alignments.all!(a => a.contigA == contigA && a.seed == locationSeed));
    debug logJsonDebug(
        "alignments", alignments.toJson,
        "commonAlignmentRegion", commonAlignmentRegion.toJson,
        "relevantContigHalf", relevantContigHalf.toJson,
        "allowedTracePointRegion", allowedTracePointRegion.toJson,
    );

    if (allowedTracePointRegion.empty)
    {
        return -1;
    }

    auto tracePointMin = ceil(min(allowedTracePointRegion), tracePointDistance);
    auto tracePointSup = ceil(sup(allowedTracePointRegion), tracePointDistance);
    auto tracePointCandidates = iota(tracePointMin, tracePointSup, tracePointDistance)
        .map!(tracePointCandidate => ReferencePoint(contigA.id, tracePointCandidate));

    return locationSeed == AlignmentLocationSeed.front
        ? getCommonTracePoint(tracePointCandidates.retro, allowedTracePointRegion)
        : getCommonTracePoint(tracePointCandidates, allowedTracePointRegion);
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

    size_t croppingTracePointIdx(in AlignmentChain.LocalAlignment localAlignment)
    {
        auto firstTracePointRefPos = ceil(localAlignment.contigA.begin, tracePointDistance);
        auto openInterval = locationSeed == AlignmentLocationSeed.front;
        assert(croppingRefPos >= firstTracePointRefPos);

        return (croppingRefPos - firstTracePointRefPos) / tracePointDistance + (openInterval ? 0 : 1);
    }

    bool coversCroppingRefPos(in AlignmentChain.LocalAlignment localAlignment)
    {
        return localAlignment.contigA.begin <= croppingRefPos
            && croppingRefPos < localAlignment.contigA.end;
    }

    auto coveringLocalAlignment = alignment
        .localAlignments
        .filter!coversCroppingRefPos
        .front;
    auto readCroppingPos =
        coveringLocalAlignment.contigB.begin +
        coveringLocalAlignment
            .tracePoints[0 .. croppingTracePointIdx(coveringLocalAlignment)]
            .map!"a.numBasePairs"
            .sum;
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

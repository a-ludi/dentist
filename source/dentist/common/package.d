/**
    This package holds common code for the `dentist` algorithm.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common;

import dentist.util.log;
import dentist.util.region : Region;
import std.algorithm : count, fold, map, max, sum;
import std.array : array;
import std.conv : to;
import std.math : floor;
import std.traits : TemplateOf;
import std.typecons : Flag;

public import dentist.common.alignments;
public import dentist.common.binio;
public import dentist.common.scaffold;


/// True iff building with testing command.
version (DentistTesting)
    enum isTesting = true;
else
    enum isTesting = false;

/// Evaluate to `value` if building with testing command;
/// otherwise to `typeof(value).init`.
template testingOnly(alias value)
{
    static if (isTesting)
        enum testingOnly = value;
    else
        enum testingOnly = typeof(value).init;
}

/// Thrown if some runtime error in the `dentist` algorithm occurs.
class DentistException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

/// A region of the reference aka. mask.
alias ReferenceRegion = Region!(size_t, size_t, "contigId");
/// An interval of a reference contig.
alias ReferenceInterval = ReferenceRegion.TaggedInterval;
/// A point on the reference.
alias ReferencePoint = ReferenceRegion.TaggedPoint;
/// A region of a read aka. mask.
alias ReadRegion = Region!(size_t, size_t, "readId");
/// An interval of a read contig.
alias ReadInterval = ReadRegion.TaggedInterval;
/// A point on a read.
alias ReadPoint = ReadRegion.TaggedPoint;

/// Returns the alignment region of alignmentChain.
R to(R, string contig = "contigA")(in AlignmentChain alignmentChain) pure
        if (__traits(isSame, TemplateOf!R, Region))
{
    auto contigId = mixin("alignmentChain." ~ contig ~ ".id");

    return R(alignmentChain
        .localAlignments
        .map!(la => R.TaggedInterval(
            contigId,
            mixin("la." ~ contig ~ ".begin"),
            mixin("la." ~ contig ~ ".end"),
        ))
        .array
    );
}

size_t insertionScore(Options)(
    in ReadAlignment readAlignment,
    in ReadAlignmentType pileUpType,
    in ReferenceRegion repeatMask,
    in Flag!"preferSpanning" preferSpanning,
    in Options options,
) pure
{
    immutable shortAnchorPenaltyMagnitude = AlignmentChain.maxScore / 512;
    immutable notSpanningPenaltyMagnitude = AlignmentChain.maxScore / 2;
    immutable improperAlignmentPenaltyMagnitude = AlignmentChain.maxScore / 8;

    long numAlignments = readAlignment.length;
    long expectedAlignmentCount = pileUpType == ReadAlignmentType.gap ? 2 : 1;
    auto alignmentAnchor = readAlignment[].map!(to!(ReferenceRegion, "contigA"))
        .fold!"a | b" - repeatMask;
    long avgAnchorSize = alignmentAnchor.size / numAlignments;
    debug long avgAlignmentLength = readAlignment[].map!"a.totalLength".sum / numAlignments;
    debug assert(avgAnchorSize <= avgAlignmentLength);
    long avgAlignmentScore = readAlignment[].map!"a.score".sum / numAlignments;
    long shortAnchorPenalty = floor(shortAnchorPenaltyMagnitude * (
            (options.goodAnchorLength + 1) / avgAnchorSize.to!float) ^^ 2).to!size_t;
    long notSpanningPenalty = preferSpanning
        ? (expectedAlignmentCount - numAlignments) * notSpanningPenaltyMagnitude
        : 0;
    long improperAlignmentPenalty = readAlignment[].count!"!a.isProper"
        * improperAlignmentPenaltyMagnitude / numAlignments;
    size_t score = max(0, (
          avgAlignmentScore
        - shortAnchorPenalty
        - notSpanningPenalty
        - improperAlignmentPenalty
    ));

    debug
    {
        size_t readId = readAlignment[0].contigB.id;
        auto contigIds = readAlignment[].map!"a.contigA.id".array;
        debug logJsonDebug(
            "readId", readId,
            "contigIds", contigIds.toJson,
            "expectedAlignmentCount", expectedAlignmentCount,
            "avgAnchorSize", avgAnchorSize,
            "avgAlignmentLength", avgAlignmentLength,
            "avgAlignmentScore", avgAlignmentScore,
            "shortAnchorPenalty", shortAnchorPenalty,
            "score", score,
        );
    }

    return score;
}

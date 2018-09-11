/**
    This is a collection of helpers for the alignment filtering.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.collectPileUps.filter;

import dentist.common : ReferenceRegion, to;
import dentist.common.alignments :
    AlignmentChain,
    haveEqualIds;
import dentist.util.math : NaturalNumberSet;
import std.algorithm :
    all,
    chunkBy,
    filter,
    joiner,
    map,
    sort,
    uniq;
import std.array : array;
import std.range :
    InputRange,
    inputRangeObject,
    walkLength;

interface AlignmentChainFilter
{
    AlignmentChain[] opCall(AlignmentChain[] alignmentChains);
}

abstract class ReadFilter : AlignmentChainFilter
{
    NaturalNumberSet* unusedReads;

    this(NaturalNumberSet* unusedReads)
    {
        this.unusedReads = unusedReads;
    }

    override AlignmentChain[] opCall(AlignmentChain[] alignmentChains)
    {
        NaturalNumberSet discardedReadIds;
        discardedReadIds.reserveFor(unusedReads.capacity);

        foreach (discardedAlignment; getDiscardedReadIds(alignmentChains))
        {
            auto discardedReadId = discardedAlignment.contigB.id;

            discardedReadIds.add(discardedReadId);
            unusedReads.remove(discardedReadId);
        }

        foreach (ref alignmentChain; alignmentChains)
        {
            auto wasReadDiscarded = discardedReadIds.has(alignmentChain.contigB.id);
            alignmentChain.disableIf(wasReadDiscarded);
        }

        return alignmentChains;
    }

    InputRange!(AlignmentChain) getDiscardedReadIds(AlignmentChain[] alignmentChains);
}

/// Discard improper alignments.
class ImproperAlignmentChainsFilter : AlignmentChainFilter
{
    override AlignmentChain[] opCall(AlignmentChain[] alignmentChains)
    {
        foreach (ref alignmentChain; alignmentChains)
        {
            alignmentChain.disableIf(!alignmentChain.isProper);
        }

        return alignmentChains;
    }
}

/// Discard read if it has an alignment that – extended with the
/// exceeding read sequence on either end – is fully contained in
/// a single contig.
class RedundantAlignmentChainsFilter : ReadFilter
{
    this(NaturalNumberSet* unusedReads)
    {
        super(unusedReads);
    }

    override InputRange!(AlignmentChain) getDiscardedReadIds(AlignmentChain[] alignmentChains)
    {
        assert(alignmentChains.map!"a.isProper || a.flags.disabled".all);

        return inputRangeObject(alignmentChains.filter!"!a.flags.disabled && a.isFullyContained");
    }
}

/// Discard read if it has an alignment that aligns in multiple
/// locations of one contig with approximately equal quality.
class AmbiguousAlignmentChainsFilter : ReadFilter
{
    /**
        Two scores are considered similar if the relative "error" of them is
        smaller than defaulMaxRelativeDiff.

        The relative error is defined as:

            relError(a, b) = abs(a - b)/max(a - b)

        **Implementation note:** all computations are done in integer
        arithmetic hence AlignmentChain.maxScore corresponds to 1 in the above
        equation.
    */
    static enum maxRelativeDiff = AlignmentChain.maxScore / 20; // 5% of larger value
    /**
        Two scores are considered similar if the absolute "error" of them is
        smaller than defaulMaxAbsoluteDiff.

        The relative error is defined as:

            absError(a, b) = abs(a - b)

        **Implementation note:** all computations are done in integer
        arithmetic hence all scores are lower than or equal to
        AlignmentChain.maxScore .
    */
    static enum maxAbsoluteDiff = AlignmentChain.maxScore / 100; // 1% wrt. score

    this(NaturalNumberSet* unusedReads)
    {
        super(unusedReads);
    }

    override InputRange!(AlignmentChain) getDiscardedReadIds(AlignmentChain[] alignmentChains)
    {
        assert(alignmentChains.map!"a.isProper || a.flags.disabled".all);

        return alignmentChains
            .filter!"!a.flags.disabled"
            .chunkBy!haveEqualIds
            .filter!isAmgiguouslyAlignedRead
            .joiner
            .inputRangeObject;
    }

    static bool isAmgiguouslyAlignedRead(Chunk)(Chunk readAlignments)
    {
        return readAlignments.save.walkLength(2) > 1 && haveSimilarScore(readAlignments);
    }

    static bool haveSimilarScore(Chunk)(Chunk readAlignments)
    {
        auto scores = readAlignments.map!"a.score".array.sort.release;

        return similarScore(scores[0], scores[1]);
    }

    protected static bool similarScore(size_t a, size_t b) pure
    {
        auto diff = a > b ? a - b : b - a;
        auto magnitude = a > b ? a : b;

        return diff < maxAbsoluteDiff || (diff * AlignmentChain.maxScore / magnitude) < maxRelativeDiff;
    }

    unittest
    {
        enum eps = 1;
        enum refValueSmall = 2 * maxAbsoluteDiff;
        enum refValueLarge = AlignmentChain.maxScore / 2;

        // test absolute part
        assert(similarScore(refValueSmall, refValueSmall));
        assert(similarScore(refValueSmall, refValueSmall + (maxAbsoluteDiff - 1)));
        assert(similarScore(refValueSmall, refValueSmall - (maxAbsoluteDiff - 1)));
        assert(!similarScore(refValueSmall, refValueSmall + (maxAbsoluteDiff + 1)));
        assert(!similarScore(refValueSmall, refValueSmall - (maxAbsoluteDiff + 1)));
        // test relative part
        assert(similarScore(refValueLarge, refValueLarge));
        assert(similarScore(refValueLarge,
                refValueLarge + refValueLarge * maxRelativeDiff / AlignmentChain.maxScore));
        assert(similarScore(refValueLarge,
                refValueLarge - refValueLarge * maxRelativeDiff / AlignmentChain.maxScore) + eps);
        assert(!similarScore(refValueLarge,
                refValueLarge + refValueLarge * 2 * maxRelativeDiff / AlignmentChain.maxScore));
        assert(!similarScore(refValueLarge,
                refValueLarge - refValueLarge * 2 * maxRelativeDiff / AlignmentChain.maxScore));
    }
}

class WeaklyAnchoredAlignmentChainsFilter : AlignmentChainFilter
{
    size_t minAnchorLength;
    const(ReferenceRegion) repetitiveRegions;

    this(const(ReferenceRegion) repetitiveRegions, size_t minAnchorLength)
    {
        this.repetitiveRegions = repetitiveRegions;
        this.minAnchorLength = minAnchorLength;
    }

    AlignmentChain[] opCall(AlignmentChain[] alignmentChains)
    {
        foreach (ref alignmentChain; alignmentChains)
        {
            alignmentChain.disableIf(isWeaklyAnchored(alignmentChain));
        }

        return alignmentChains;
    }

    bool isWeaklyAnchored(AlignmentChain alignment)
    {
        // Mark reads as weakly anchored if they are mostly anchored in a
        // repetitive region
        auto uniqueAlignmentRegion = alignment.to!(ReferenceRegion, "contigA") - repetitiveRegions;
        bool isWeaklyAnchored = uniqueAlignmentRegion.size <= minAnchorLength;

        return isWeaklyAnchored;
    }
}

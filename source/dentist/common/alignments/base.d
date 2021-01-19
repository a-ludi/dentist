/**
    Everything to handle local alignments and friends.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.alignments.base;

import dentist.common.scaffold :
    concatenatePayloads,
    ContigNode,
    ContigPart,
    Join;
import dentist.util.algorithm : cmpLexicographically, orderLexicographically;
import dentist.util.log;
import dentist.util.math : absdiff, ceildiv, floor, RoundingMode;
import core.exception : AssertError;
import std.algorithm :
    all,
    among,
    any,
    cache,
    canFind,
    chunkBy,
    copy,
    countUntil,
    cumulativeFold,
    equal,
    filter,
    find,
    isSorted,
    joiner,
    map,
    max,
    maxIndex,
    mean,
    min,
    minIndex,
    reverse,
    sort,
    sum,
    swap,
    SwapStrategy,
    uniq;
import std.array : appender, array, minimallyInitializedArray;
import std.conv : to;
import std.exception : assertNotThrown, assertThrown, enforce, ErrnoException;
import std.format : format;
import std.math : sgn;
import std.parallelism : defaultTaskPool = taskPool, TaskPool;
import std.range :
    assumeSorted,
    chain,
    chunks,
    ElementType,
    enumerate,
    InputRange,
    inputRangeObject,
    iota,
    isInputRange,
    only,
    radial,
    repeat,
    retro,
    slide,
    tee,
    takeNone,
    zip;
import std.range.primitives;
import std.string : capitalize, split;
import std.stdio : File, LockType;
import std.typecons : BitFlags, PhobosFlag = Flag, No, tuple, Tuple, Yes;
import std.traits : isArray, TemplateArgsOf, TemplateOf;
import vibe.data.json : Json, toJson = serializeToJson;

debug import std.stdio : writefln, writeln;

alias arithmetic_t = int;
alias coord_t = uint;
alias diff_t = uint;
alias id_t = uint;
alias trace_point_t = ushort;



struct Contig
{
    id_t id;
    coord_t length;
}


struct Locus
{
    coord_t begin;
    coord_t end;

    @property coord_t length() const pure nothrow
    {
        return end - begin;
    }
}


enum Flag : ubyte
{
    complement = 1 << 0,
    disabled = 1 << 1,
    alternateChain = 1 << 2,
}


alias Flags = BitFlags!Flag;


struct TracePoint
{
    trace_point_t numDiffs;
    trace_point_t numBasePairs;
}


struct TranslatedTracePoint
{
    coord_t contigA;
    coord_t contigB;
}


struct Trace
{
    Locus contigA;
    Locus contigB;
    trace_point_t tracePointDistance;
    const(TracePoint)[] tracePoints;


    TranslatedTracePoint translateTracePoint(string contig)(
        coord_t contigPos,
        RoundingMode roundingMode,
    ) const pure if (contig.among("contigA", "contigB"))
    {
        assert(mixin(contig ~ `.begin <= contigPos && contigPos <= ` ~ contig ~ `.end`));

        auto tracePointIndex = tracePointsUpTo!contig(contigPos, roundingMode);
        auto contigBPos = contigB.begin + tracePoints[0 .. tracePointIndex]
                .map!"a.numBasePairs"
                .sum;
        auto contigAPos = tracePointIndex == 0
            ? contigA.begin
            : tracePointIndex < tracePoints.length
                ? floor(contigA.begin, tracePointDistance) + cast(coord_t) (tracePointIndex * tracePointDistance)
                : contigA.end;

        return TranslatedTracePoint(contigAPos, contigBPos);
    }


    auto tracePointsUpTo(string contig)(
        coord_t contigAPos,
        RoundingMode roundingMode,
    ) const pure nothrow if (contig == "contigA")
    {
        assert(contigA.begin <= contigAPos && contigAPos <= contigA.end);

        auto firstTracePointRefPos = contigA.begin;
        auto secondTracePointRefPos = floor(firstTracePointRefPos, tracePointDistance) + tracePointDistance;
        auto secondFromLastTracePointRefPos = floor(contigA.end - 1, tracePointDistance);

        final switch (roundingMode)
        {
        case RoundingMode.floor:
            if (contigAPos < secondTracePointRefPos)
                return 0;
            if (contigAPos < contigA.end)
                return 1 + (contigAPos - secondTracePointRefPos) / tracePointDistance;
            else
                return tracePoints.length;
        case RoundingMode.round:
            assert(0, "unimplemented");
        case RoundingMode.ceil:
            if (firstTracePointRefPos == contigAPos)
                return 0;
            if (contigAPos <= secondTracePointRefPos)
                return 1;
            else if (contigAPos <= secondFromLastTracePointRefPos)
                return 1 + ceildiv(contigAPos - secondTracePointRefPos, tracePointDistance);
            else
                return tracePoints.length;
        }
    }

    unittest
    {
        auto trace = Trace(
            Locus(50, 2897),
            Locus(50, 2905),
            100,
            [TracePoint(1, 50), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100), TracePoint(0, 100),
             TracePoint(7, 105)],
        );
        coord_t contigPos = 79;
        auto roundingMode = RoundingMode.ceil;
        auto tpsUpTo = trace.tracePointsUpTo!"contigA"(contigPos, roundingMode);

        assert(0 <= tpsUpTo && tpsUpTo <= trace.tracePoints.length);
    }

    auto tracePointsUpTo(string contig)(
        coord_t contigBPos,
        RoundingMode roundingMode,
    ) const pure nothrow if (contig == "contigB")
    {
        assert(contigB.begin <= contigBPos && contigBPos <= contigB.end);

        if (contigBPos == contigB.begin)
            return 0;
        else if (contigBPos == contigB.end)
            return tracePoints.length;

        auto tracePointPositions = chain(only(contigB.begin), tracePoints.map!"a.numBasePairs")
            .cumulativeFold!"a + b"
            .enumerate;

        final switch (roundingMode)
        {
        case RoundingMode.floor:
            return tracePointPositions
                .find!(pair => contigBPos < pair.value)
                .front
                .index - 1;
        case RoundingMode.round:
            assert(0, "unimplemented");
        case RoundingMode.ceil:
            return tracePointPositions
                .find!(pair => contigBPos <= pair.value)
                .front
                .index;
        }
    }
}


/**
    Holds a chain of local alignments that form a compound alignment. An AlignmentChain should
    contain at least one element.
*/
struct AlignmentChain
{
    static struct LocalAlignment
    {
        alias Locus = .Locus;
        alias TracePoint = .TracePoint;

        Locus contigA;
        Locus contigB;
        diff_t numDiffs;
        TracePoint[] tracePoints;


        Trace getTrace(trace_point_t tracePointDistance) const pure nothrow @safe
        {
            return Trace(
                cast(Locus) contigA,
                cast(Locus) contigB,
                tracePointDistance,
                tracePoints,
            );
        }


        TranslatedTracePoint translateTracePoint(string contig = "contigA")(
            coord_t contigPos,
            trace_point_t tracePointDistance,
            RoundingMode roundingMode,
        ) const pure if (contig.among("contigA", "contigB"))
        {
            return getTrace(tracePointDistance).translateTracePoint!contig(contigPos, roundingMode);
        }

        /// Crops this local alignment from startingSeed to contigPos.
        void cropToTracePoint(string contig = "contigA")(
            in AlignmentLocationSeed startingSeed,
            in coord_t contigPos,
            in trace_point_t tracePointDistance,
            in RoundingMode roundingMode,
        ) pure if (contig.among("contigA", "contigB"))
        {
            auto index = tracePointsUpTo!contig(contigPos, tracePointDistance, roundingMode);
            auto newLocusBoundary = translateTracePoint!contig(
                contigPos,
                tracePointDistance,
                roundingMode,
            );

            final switch (startingSeed)
            {
            case AlignmentLocationSeed.front:
                contigA.end = newLocusBoundary.contigA;
                contigB.end = newLocusBoundary.contigB;
                tracePoints = tracePoints[0 .. index];
                break;
            case AlignmentLocationSeed.back:
                contigA.begin = newLocusBoundary.contigA;
                contigB.begin = newLocusBoundary.contigB;
                tracePoints = tracePoints[index .. $];
                break;
            }

            numDiffs = tracePoints.map!(tp => tp.numDiffs.to!diff_t).sum;
        }

        auto tracePointsUpTo(string contig)(
            coord_t contigPos,
            trace_point_t tracePointDistance,
            RoundingMode roundingMode,
        ) const pure nothrow if (contig.among("contigA", "contigB"))
        {
            return getTrace(tracePointDistance).tracePointsUpTo!contig(
                contigPos,
                roundingMode,
            );
        }
    }

    static enum maxScore = 2 ^^ 16;

    id_t id;
    Contig contigA;
    Contig contigB;
    Flags flags;
    LocalAlignment[] localAlignments;
    trace_point_t tracePointDistance;

    static @property AlignmentChain disabledInstance()
    {
        AlignmentChain ac;

        ac.flags |= Flag.disabled;

        return ac;
    }

    static foreach(flagName; __traits(allMembers, Flag))
    {
        mixin(format!(q"<
            static alias %1$s = PhobosFlag!"%2$s";

            @property PhobosFlag!"%2$s" %2$s() pure const nothrow @trusted
            {
                return cast(PhobosFlag!"%2$s") flags.%2$s;
            }

            @property void %2$s(PhobosFlag!"%2$s" %2$s) pure nothrow
            {
                flags.%2$s = %2$s;
            }
        >")(flagName.capitalize, flagName));
    }

    invariant
    {
        if (flags.disabled)
            return;

        assert(localAlignments.length >= 1, "empty chain is forbidden");
        foreach (la; localAlignments)
        {
            assert(0 <= la.contigA.begin && la.contigA.begin < la.contigA.end
                    && (contigA.length == 0 || la.contigA.end <= contigA.length), "non-sense alignment of contigA");
            assert(0 <= la.contigB.begin && la.contigB.begin < la.contigB.end
                    && (contigB.length == 0 || la.contigB.end <= contigB.length), "non-sense alignment of contigB");

            assert(tracePointDistance == 0 || la.tracePoints.length > 0, "missing trace points");
            if (tracePointDistance > 0)
            {
                coord_t traceLength = la.tracePoints.map!(tp => tp.numBasePairs.to!coord_t).sum;
                coord_t traceDiffs = la.tracePoints.map!(tp => tp.numDiffs.to!coord_t).sum;

                // TODO remove logging if fixed in LAdump (issue #23)
                if (la.numDiffs != traceDiffs)
                {
                    debug logJsonDebug(
                        "contigA", contigA.id,
                        "contigB", contigB.id,
                        "la.numDiffs", la.numDiffs,
                        "traceDiffs", traceDiffs,
                    );
                }
                // TODO make equality assertion if fixed in LAdump (issue #23)
                assert(la.numDiffs <= traceDiffs, "missing trace points");
                assert(la.contigB.end - la.contigB.begin == traceLength,
                        "trace distance does not match alignment");
            }
        }
    }

    unittest
    {
        with (LocalAlignment)
            {
                auto acZeroLength = AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), []);
                auto ac1 = AlignmentChain(1, Contig(1, 10), Contig(1, 10), Flags(),
                        [LocalAlignment(Locus(1, 1), Locus(1, 9), 0)]);
                auto ac2 = AlignmentChain(2, Contig(1, 10), Contig(1, 10), Flags(),
                        [LocalAlignment(Locus(1, 11), Locus(1, 9), 0)]);
                auto ac3 = AlignmentChain(3, Contig(1, 10), Contig(1, 10), Flags(),
                        [LocalAlignment(Locus(1, 9), Locus(1, 1), 0)]);
                auto ac4 = AlignmentChain(4, Contig(1, 10), Contig(1, 10), Flags(),
                        [LocalAlignment(Locus(1, 9), Locus(1, 11), 0)]);
                auto acFine = AlignmentChain(5, Contig(1, 10), Contig(1, 10),
                        Flags(), [LocalAlignment(Locus(1, 9), Locus(1, 9), 0)]);

                assertThrown!AssertError(acZeroLength.totalLength);
                assertThrown!AssertError(ac1.totalLength);
                assertThrown!AssertError(ac2.totalLength);
                assertThrown!AssertError(ac3.totalLength);
                assertThrown!AssertError(ac4.totalLength);
                assertNotThrown!AssertError(acFine.totalLength);
            }
    }

    @property ref const(LocalAlignment) first() const pure nothrow @safe
    {
        return localAlignments[0];
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto firstLA = LocalAlignment(Locus(1, 9), Locus(1, 9), 0);
                auto otherLA = LocalAlignment(Locus(9, 10), Locus(9, 10), 0);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [firstLA, otherLA]);

                assert(ac.first == firstLA);
            }
    }

    @property ref const(LocalAlignment) last() const pure nothrow @safe
    {
        return localAlignments[$ - 1];
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto lastLA = LocalAlignment(Locus(1, 9), Locus(1, 9), 0);
                auto otherLA = LocalAlignment(Locus(9, 10), Locus(9, 10), 0);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [otherLA, lastLA]);

                assert(ac.last == lastLA);
            }
    }

    PhobosFlag!"disabled" disableIf(lazy bool disable) pure
    {
        if (!flags.disabled)
        {
            flags.disabled = disable;
        }


        return disabled;
    }

    /// This alignment is called proper iff it starts and ends at a read boundary.
    @property bool isProper(coord_t allowance = 0) const pure nothrow @safe
    {
        return (beginsWith!"contigA" || beginsWith!"contigB") &&
               (endsWith!"contigA" || endsWith!"contigB");
    }

    @property bool beginsWith(string contig)(coord_t allowance = 0) const pure nothrow @safe
        if (contig.among("contigA", "contigB"))
    {
        return mixin("first."~contig~".begin <= allowance");
    }

    @property bool endsWith(string contig)(coord_t allowance = 0) const pure nothrow @safe
        if (contig.among("contigA", "contigB"))
    {
        return mixin("last."~contig~".end + allowance >= "~contig~".length");
    }

    /// Returns true iff this alignment covers `contig` completely.
    @property bool completelyCovers(string contig)(coord_t allowance = 0) const pure nothrow
        if (contig.among("contigA", "contigB"))
    {
        return beginsWith!contig(allowance) && endsWith!contig(allowance);
    }

    /**
        Returns true if the aligned read `contigB` (with extensions on either
        end) is fully contained in the reference `contigA`.

        According to the following 'definition' `contigA` is fully contained
        in `contigB` iff `x >= 0` and `y <= l_a`.
        ---
                0                  x   ua      va  y                        la
        contigA |------------------+---+-+---+-+---+-------------------------|
                                  / / /  | | |  \ \ \
                                 / / /   | | |   \ \ \
                                / / /    | | |    \ \ \
                               / / /     | | |     \ \ \
        contigB               |---+------+---+------+---|
                              0   ub                vb lb

        x = ua - (ub - 0) = ua - ub
        y = va + (lb - vb)
        ---
    */
    bool isFullyContained() const
    {
        if (first.contigB.begin > first.contigA.begin)
        {
            // x < 0; return early to avoid negative numbers in unsigned integers
            return false;
        }

        auto x = first.contigA.begin - first.contigB.begin;
        auto y = last.contigA.end + contigB.length - last.contigB.end;

        return 0 <= x && y < contigA.length;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(30, 35), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), Flags(), [la]);

                // read with extension align an contigA from 25 to 40
                assert(ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(10, 20), Locus(5, 10), 1);
                auto la2 = LocalAlignment(Locus(30, 40), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), Flags(), [la1, la2]);

                // read with extension align an contigA from 5 to 45
                assert(ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(0, 10), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), Flags(), [la]);

                // read with extension align an contigA from -5 to 15
                assert(!ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(40, 50), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), Flags(), [la]);

                // read with extension align an contigA from 35 to 55
                assert(!ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(0, 20), Locus(5, 10), 1);
                auto la2 = LocalAlignment(Locus(30, 50), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), Flags(), [la1, la2]);

                // read with extension align an contigA from -5 to 55
                assert(!ac.isFullyContained);
            }
    }

    @property size_t totalLength() const pure
    {
        return last.contigA.end - first.contigA.begin;
    }

    @property size_t coveredBases(string contig)() const pure
    {
        return localAlignments.map!("a." ~ contig ~ ".end - a." ~ contig ~ ".begin").sum;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [la1, la2]);

                assert(ac.totalLength == 9);
            }
    }

    @property size_t totalDiffs() const pure
    {
        return localAlignments.map!"a.numDiffs".sum;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [la1, la2]);

                assert(ac.totalDiffs == 3);
            }
    }

    @property size_t totalGapLength() const pure
    {
        return localAlignments
            .chunks(2)
            .map!(las => las.length < 2 ? 0 : absdiff(las[1].contigA.begin, las[0].contigA.end))
            .sum;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [la1, la2]);

                assert(ac.totalGapLength == 2);
            }
    }

    @property size_t numMatchingBps() const pure
    {
        return totalLength - (totalDiffs + totalGapLength);
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [la1, la2]);

                assert(ac.numMatchingBps == 9 - (3 + 2));
            }
    }

    @property size_t score() const pure
    {
        return numMatchingBps * maxScore / totalLength;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [la1, la2]);

                assert(ac.score == 4 * maxScore / 9);
            }
    }

    @property double averageErrorRate() const pure
    {
        return totalDiffs.to!double / coveredBases!"contigA".to!double;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [la1, la2]);

                assert(ac.averageErrorRate == 3.0 / 7.0);
            }
    }

    int compareIds(ref const AlignmentChain other) const pure nothrow
    {
        return cmpLexicographically!(
            typeof(this),
            ac => ac.contigA.id,
            ac => ac.contigB.id,
        )(this, other);
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(0, 1), Locus(0, 1), 1);
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [la]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 10), Flags(), [la]),
                    AlignmentChain(2, Contig(2, 10), Contig(1, 10), Flags(), [la]),
                    AlignmentChain(3, Contig(2, 10), Contig(2, 10), Flags(), [la]),
                ];

                foreach (i; 0 .. acs.length)
                    foreach (j; 0 .. acs.length)
                    {
                        auto compareValue = acs[i].compareIds(acs[j]);
                        alias errorMessage = (cmp) => format!"expected %s and %s to compare %s 0 but got %d"(
                                acs[i], acs[j], cmp, compareValue);

                        if (i < j)
                            assert(compareValue < 0, errorMessage("<"));
                        else if (i > j)
                            assert(compareValue > 0, errorMessage(">"));
                        else
                            assert(compareValue == 0, errorMessage("=="));
                    }
            }
    }

    int opCmp(ref const AlignmentChain other) const pure nothrow
    {
        return cmpLexicographically!(
            typeof(this),
            ac => ac.contigA.id,
            ac => ac.contigB.id,
            ac => ac.first.contigA.begin,
            ac => ac.first.contigB.begin,
            ac => ac.last.contigA.end,
            ac => ac.last.contigB.end,
        )(this, other);
    }

    unittest
    {
        // see compareIds
        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(0, 1), Locus(0, 1), 1);
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [la]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 10), Flags(), [la]),
                    AlignmentChain(2, Contig(2, 10), Contig(1, 10), Flags(), [la]),
                    AlignmentChain(3, Contig(2, 10), Contig(2, 10), Flags(), [la]),
                ];

                foreach (i; 0 .. acs.length)
                    foreach (j; 0 .. acs.length)
                    {
                        auto compareValue = acs[i].opCmp(acs[j]);
                        alias errorMessage = (cmp) => format!"expected %s and %s to compare %s 0 but got %d"(
                                acs[i], acs[j], cmp, compareValue);

                        if (i < j)
                            assert(compareValue < 0, errorMessage("<"));
                        else if (i > j)
                            assert(compareValue > 0, errorMessage(">"));
                        else
                            assert(compareValue == 0, errorMessage("=="));
                    }
            }

        // test non-id-related comparison
        with (Complement) with (LocalAlignment)
            {
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [
                        LocalAlignment(Locus(0, 2), Locus(0, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(1, Contig(1, 10), Contig(1, 10), Flags(), [
                        LocalAlignment(Locus(0, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(2, Contig(1, 10), Contig(1, 10), Flags(), [
                        LocalAlignment(Locus(1, 2), Locus(0, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(3, Contig(1, 10), Contig(1, 10), Flags(), [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(4, Contig(1, 10), Contig(1, 10), Flags(), [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 6), 1)
                    ]),
                    AlignmentChain(5, Contig(1, 10), Contig(1, 10), Flags(), [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 6), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(6, Contig(1, 10), Contig(1, 10), Flags(), [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 6), Locus(3, 6), 1)
                    ]),
                ];

                foreach (i; 0 .. acs.length)
                    foreach (j; 0 .. acs.length)
                    {
                        auto compareValue = acs[i].opCmp(acs[j]);
                        alias errorMessage = (cmp) => format!"expected %s and %s to compare %s 0 but got %d"(
                                acs[i], acs[j], cmp, compareValue);

                        if (i < j)
                            assert(compareValue < 0, errorMessage("<"));
                        else if (i > j)
                            assert(compareValue > 0, errorMessage(">"));
                        else
                            assert(compareValue == 0, errorMessage("=="));
                    }
            }
    }

    TranslatedTracePoint translateTracePoint(string contig = "contigA")(
        in coord_t contigPos,
        RoundingMode roundingMode,
    ) const pure if (contig.among("contigA", "contigB"))
    {
        auto index = coveringLocalAlignmentIndex!contig(contigPos, roundingMode);
        auto coveringLocalAlignment = localAlignments[index];

        return coveringLocalAlignment.translateTracePoint!contig(
            contigPos,
            tracePointDistance,
            roundingMode,
        );
    }

    unittest
    {
        alias Locus = LocalAlignment.Locus;
        alias TracePoint = LocalAlignment.TracePoint;
        enum tracePointDistance = 100;
        auto ac = AlignmentChain(
            0,
            Contig(1, 2584),
            Contig(58024, 10570),
            Flags(Flag.complement),
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
        );

        assert(ac.translateTracePoint(579, RoundingMode.floor) == TranslatedTracePoint(579, 0));
        assert(ac.translateTracePoint(599, RoundingMode.floor) == TranslatedTracePoint(579, 0));
        assert(ac.translateTracePoint(600, RoundingMode.floor) == TranslatedTracePoint(600, 23));
        assert(ac.translateTracePoint(699, RoundingMode.floor) == TranslatedTracePoint(600, 23));
        assert(
            ac.translateTracePoint(700, RoundingMode.ceil)
            ==
            ac.translateTracePoint(700, RoundingMode.floor)
        );
        assert(
            ac.translateTracePoint(699, RoundingMode.ceil)
            ==
            ac.translateTracePoint(701, RoundingMode.floor)
        );
        assert(ac.translateTracePoint(700, RoundingMode.floor) == TranslatedTracePoint(700, 23 + 109));
        assert(ac.translateTracePoint(799, RoundingMode.floor) == TranslatedTracePoint(700, 23 + 109));
        assert(ac.translateTracePoint(800, RoundingMode.floor) == TranslatedTracePoint(800, 23 + 109 + 109));
        assert(ac.translateTracePoint(899, RoundingMode.floor) == TranslatedTracePoint(800, 23 + 109 + 109));
        assert(ac.translateTracePoint(2583, RoundingMode.floor) == TranslatedTracePoint(2500, 2070));
        assert(ac.translateTracePoint(2584, RoundingMode.floor) == TranslatedTracePoint(2584, 2158));
        assertThrown!Exception(ac.translateTracePoint(578, RoundingMode.floor));
        assertThrown!Exception(ac.translateTracePoint(2585, RoundingMode.floor));
    }

    /// Crops this alignment chain from startingSeed to contigPos.
    void cropToTracePoint(string contig = "contigA")(
        in AlignmentLocationSeed startingSeed,
        in coord_t contigPos,
        in RoundingMode roundingMode,
    ) pure if (contig.among("contigA", "contigB"))
    {
        auto index = coveringLocalAlignmentIndex!contig(contigPos, roundingMode);

        localAlignments[index].cropToTracePoint!contig(
            startingSeed,
            contigPos,
            tracePointDistance,
            roundingMode,
        );
        auto isEmptyCroppedLA = localAlignments[index].contigA.length == 0 ||
                                localAlignments[index].contigB.length == 0;

        final switch (startingSeed)
        {
        case AlignmentLocationSeed.front:
            if (isEmptyCroppedLA)
                --index;
            localAlignments = localAlignments[0 .. index + 1];
            break;
        case AlignmentLocationSeed.back:
            if (isEmptyCroppedLA)
                ++index;
            localAlignments = localAlignments[index .. $];
            break;
        }

        if (localAlignments.length == 0)
            flags.disabled = true;
    }

    unittest
    {
        enum front = AlignmentLocationSeed.front;
        enum back = AlignmentLocationSeed.back;
        alias Locus = LocalAlignment.Locus;
        alias TracePoint = LocalAlignment.TracePoint;
        enum tracePointDistance = 100;
        alias getAC = () => AlignmentChain(
            0,
            Contig(1, 2584),
            Contig(58024, 10570),
            Flags(Flag.complement),
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
        );
        enum originalAC = getAC();

        {
            auto ac = getAC();
            ac.cropToTracePoint(front, 579, RoundingMode.floor);

            assert(ac.flags.disabled);
            assert(ac.localAlignments.length == 0);
        }
        {
            auto ac = getAC();
            ac.cropToTracePoint(front, 599, RoundingMode.floor);

            assert(ac.flags.disabled);
            assert(ac.localAlignments.length == 0);
        }
        {
            auto ac = getAC();
            ac.cropToTracePoint(front, 600, RoundingMode.floor);

            assert(!ac.flags.disabled);
            assert(ac.localAlignments.length == 1);
            assert(ac.first.contigA.begin == originalAC.first.contigA.begin);
            assert(ac.first.contigA.end == 600);
            assert(ac.first.contigB.begin == originalAC.first.contigB.begin);
            assert(ac.first.contigB.end == 23);
        }
        {
            auto ac = getAC();
            ac.cropToTracePoint(front, 699, RoundingMode.floor);

            assert(!ac.flags.disabled);
            assert(ac.localAlignments.length == 1);
            assert(ac.first.contigA.begin == originalAC.first.contigA.begin);
            assert(ac.first.contigA.end == 600);
            assert(ac.first.contigB.begin == originalAC.first.contigB.begin);
            assert(ac.first.contigB.end == 23);
        }
        {
            auto ac1 = getAC();
            ac1.cropToTracePoint(front, 699, RoundingMode.ceil);
            auto ac2 = getAC();
            ac2.cropToTracePoint(front, 701, RoundingMode.floor);

            assert(ac1 == ac2);
        }
        {
            auto ac = getAC();
            ac.cropToTracePoint(back, 2584, RoundingMode.ceil);

            assert(ac.flags.disabled);
            assert(ac.localAlignments.length == 0);
        }
        {
            auto ac = getAC();
            ac.cropToTracePoint(back, 2501, RoundingMode.ceil);

            assert(ac.flags.disabled);
            assert(ac.localAlignments.length == 0);
        }
        {
            auto ac = getAC();
            ac.cropToTracePoint(back, 2500, RoundingMode.ceil);
            import std.stdio;

            assert(!ac.flags.disabled);
            assert(ac.localAlignments.length == 1);
            assert(ac.last.contigA.end == originalAC.last.contigA.end);
            assert(ac.last.contigA.begin == 2500);
            assert(ac.last.contigB.end == originalAC.last.contigB.end);
            assert(ac.last.contigB.begin == 2070);
        }
        {
            auto ac = getAC();
            ac.cropToTracePoint(back, 2401, RoundingMode.ceil);

            assert(!ac.flags.disabled);
            assert(ac.localAlignments.length == 1);
            assert(ac.last.contigA.end == originalAC.last.contigA.end);
            assert(ac.last.contigA.begin == 2500);
            assert(ac.last.contigB.end == originalAC.last.contigB.end);
            assert(ac.last.contigB.begin == 2070);
        }
        {
            auto ac1 = getAC();
            ac1.cropToTracePoint(back, 2401, RoundingMode.ceil);
            auto ac2 = getAC();
            ac2.cropToTracePoint(back, 2500, RoundingMode.floor);

            assert(ac1 == ac2);
        }

        assertThrown!Exception(getAC().cropToTracePoint(front, 578, RoundingMode.floor));
        assertThrown!Exception(getAC().cropToTracePoint(front, 2585, RoundingMode.floor));
    }

    protected size_t coveringLocalAlignmentIndex(string contig = "contigA")(
        in coord_t contigPos,
        in RoundingMode roundingMode,
    ) const pure if (contig.among("contigA", "contigB"))
    {
        bool coversContigPos(in AlignmentChain.LocalAlignment localAlignment)
        {
            return mixin(`localAlignment.` ~ contig ~ `.begin <= contigPos
                && contigPos <= localAlignment.` ~ contig ~ `.end`);
        }

        auto index = localAlignments.countUntil!coversContigPos;
        enforce!Exception(
            index >= 0,
            "cannot translate coordinate due to lack of alignment coverage",
        );

        return index;
    }

    /**
        Generate a cartoon of this alignment relative to `contig`.

        Params:
            bpsPerChar =    Number of base pairs that one char represents.

        Returns: a cartoon of this alignment
    */
    static string cartoon(string contig)(in coord_t bpsPerChar, in AlignmentChain[] alignmentChains...)
    {
        if (alignmentChains.length == 0)
            return "";

        alias getContig = ac => mixin("ac." ~ contig);
        alias getCoords = la => mixin("la." ~ contig);

        auto referenceContig = getContig(alignmentChains[0]);

        enforce!Exception(
            alignmentChains.all!(ac => getContig(ac) == referenceContig),
            "all alignment chains must share the same reference contig",
        );

        InputRange!char cartoonLine(in AlignmentChain ac)
        {
            if (contig == "contigA" || !ac.flags.complement)
                return inputRangeObject(chain(
                    ' '.repeat(getCoords(ac.first).begin / bpsPerChar),
                    ac
                        .localAlignments
                        .slide!(No.withPartial)(2)
                        .map!(laPair => getCoords(laPair[1]).begin > getCoords(laPair[0]).end
                            ? chain(
                                '-'.repeat(ceildiv(getCoords(laPair[0]).length, bpsPerChar)),
                                '='.repeat((getCoords(laPair[1]).begin - getCoords(laPair[0]).end) / bpsPerChar),
                            )
                            : chain(
                                '-'.repeat(ceildiv(
                                    getCoords(laPair[0]).length - (getCoords(laPair[0]).end - getCoords(laPair[1]).begin),
                                    bpsPerChar,
                                )),
                                '='.repeat(0),
                            ),
                        )
                        .joiner,
                    '-'.repeat(ceildiv(getCoords(ac.last).length, bpsPerChar)),
                ));
            else
                return inputRangeObject(chain(
                    ' '.repeat((referenceContig.length - getCoords(ac.last).end) / bpsPerChar),
                    ac
                        .localAlignments
                        .retro
                        .slide!(No.withPartial)(2)
                        .map!(laPair => getCoords(laPair[0]).begin > getCoords(laPair[1]).end
                            ? chain(
                                '-'.repeat(ceildiv(getCoords(laPair[0]).length, bpsPerChar)),
                                '='.repeat((getCoords(laPair[0]).begin - getCoords(laPair[1]).end) / bpsPerChar),
                            )
                            : chain(
                                '-'.repeat(ceildiv(
                                    getCoords(laPair[0]).length - (getCoords(laPair[1]).end - getCoords(laPair[0]).begin),
                                    bpsPerChar,
                                )),
                                '='.repeat(0),
                            )
                        )
                        .joiner,
                    '-'.repeat(ceildiv(getCoords(ac.first).length, bpsPerChar)),
                ));
        }

        return chain(
            '-'.repeat(ceildiv(referenceContig.length, bpsPerChar)),
            only('\n'),
            alignmentChains.map!cartoonLine.joiner(only('\n')),
        ).to!string;
    }

    ///
    unittest
    {
        alias LocalAlignment = AlignmentChain.LocalAlignment;

        auto acs = [
            AlignmentChain(
                0,
                Contig(1, 10),
                Contig(1, 10),
                Flags(),
                [
                    LocalAlignment(Locus(0, 3), Locus(0, 3)),
                    LocalAlignment(Locus(4, 5), Locus(4, 5)),
                ]),
            AlignmentChain(
                1,
                Contig(1, 10),
                Contig(1, 10),
                Flags(Flag.complement),
                [
                    LocalAlignment(Locus(5, 8), Locus(0, 3)),
                    LocalAlignment(Locus(9, 10), Locus(4, 5)),
                ]),
        ];

        assert(cartoon!"contigA"(1, acs) == "----------\n" ~
                                            "---=-\n" ~
                                            "     ---=-");
        assert(cartoon!"contigB"(1, acs) == "----------\n" ~
                                            "---=-\n" ~
                                            "     -=---");
    }
}

bool idsPred(in AlignmentChain ac1, in AlignmentChain ac2) pure
{
    auto cmpValue = ac1.compareIds(ac2);

    return 0 != cmpValue && cmpValue < 0;
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;

    auto la = LocalAlignment(Locus(0, 1), Locus(0, 1), 1);
    auto acs = [
        AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [la]),
        AlignmentChain(1, Contig(1, 10), Contig(2, 10), Flags(), [la]),
        AlignmentChain(2, Contig(2, 10), Contig(1, 10), Flags(), [la]),
        AlignmentChain(3, Contig(2, 10), Contig(2, 10), Flags(), [la]),
    ];

    foreach (i; 0 .. acs.length)
        foreach (j; 0 .. acs.length)
        {
            auto compareValue = idsPred(acs[i], acs[j]);
            alias errorMessage = (expValue) => format!"expected idsPred(%s, %s) to be %s but got %s"(
                    acs[i], acs[j], expValue, compareValue);

            if (i < j)
                assert(compareValue, errorMessage(true));
            else
                assert(!compareValue, errorMessage(false));
        }
}

auto haveEqualIds(in AlignmentChain ac1, in AlignmentChain ac2) pure
{
    auto cmpValue = ac1.compareIds(ac2);

    return 0 == cmpValue;
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;

    auto sortedTestChains = [
        AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [LocalAlignment(Locus(0, 10), Locus(0, 1), 0)]),
        AlignmentChain(1, Contig(1, 10), Contig(2, 20), Flags(), [LocalAlignment(Locus(0, 10), Locus(0, 2), 0)]),
        AlignmentChain(2, Contig(1, 10), Contig(3, 30), Flags(), [LocalAlignment(Locus(0, 10), Locus(0, 3), 0)]),
        AlignmentChain(3, Contig(2, 20), Contig(1, 10), Flags(), [LocalAlignment(Locus(0, 20), Locus(0, 4), 0)]),
        AlignmentChain(4, Contig(2, 20), Contig(1, 10), Flags(), [LocalAlignment(Locus(0, 20), Locus(0, 5), 0)]),
        AlignmentChain(5, Contig(2, 20), Contig(3, 30), Flags(), [LocalAlignment(Locus(0, 20), Locus(0, 6), 0)]),
        AlignmentChain(6, Contig(3, 30), Contig(1, 10), Flags(), [LocalAlignment(Locus(0, 30), Locus(0, 7), 0)]),
        AlignmentChain(7, Contig(3, 30), Contig(2, 20), Flags(), [LocalAlignment(Locus(0, 30), Locus(0, 8), 0)]),
        AlignmentChain(8, Contig(3, 30), Contig(3, 30), Flags(), [LocalAlignment(Locus(0, 30), Locus(0, 9), 0)]),
    ];

    assert(sortedTestChains.chunkBy!haveEqualIds.equal!equal([
        sortedTestChains[0..1],
        sortedTestChains[1..2],
        sortedTestChains[2..3],
        sortedTestChains[3..5],
        sortedTestChains[5..6],
        sortedTestChains[6..7],
        sortedTestChains[7..8],
        sortedTestChains[8..9],
    ]));
}

auto equalIdsRange(in AlignmentChain[] acList, in id_t contigAID, in id_t contigBID) pure
{
    assert(isSorted!idsPred(acList));
    AlignmentChain needle;
    needle.contigA = Contig(contigAID, 1);
    needle.contigB = Contig(contigBID, 1);
    needle.localAlignments = [AlignmentChain.LocalAlignment(Locus(0, 1), Locus(0, 1))];

    return acList.assumeSorted!idsPred.equalRange(needle);
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;

    auto sortedTestChains = [
        AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [LocalAlignment(Locus(0, 1), Locus(0, 1), 0)]),
        AlignmentChain(1, Contig(1, 10), Contig(2, 20), Flags(), [LocalAlignment(Locus(0, 2), Locus(0, 2), 0)]),
        AlignmentChain(2, Contig(1, 10), Contig(3, 30), Flags(), [LocalAlignment(Locus(0, 3), Locus(0, 3), 0)]),
        AlignmentChain(3, Contig(2, 20), Contig(1, 10), Flags(), [LocalAlignment(Locus(0, 4), Locus(0, 4), 0)]),
        AlignmentChain(4, Contig(2, 20), Contig(1, 10), Flags(), [LocalAlignment(Locus(0, 5), Locus(0, 5), 0)]),
        AlignmentChain(5, Contig(2, 20), Contig(3, 30), Flags(), [LocalAlignment(Locus(0, 6), Locus(0, 6), 0)]),
        AlignmentChain(6, Contig(3, 30), Contig(1, 10), Flags(), [LocalAlignment(Locus(0, 7), Locus(0, 7), 0)]),
        AlignmentChain(7, Contig(3, 30), Contig(2, 20), Flags(), [LocalAlignment(Locus(0, 8), Locus(0, 8), 0)]),
        AlignmentChain(8, Contig(3, 30), Contig(3, 30), Flags(), [LocalAlignment(Locus(0, 9), Locus(0, 9), 0)]),
    ];

    assert(sortedTestChains.equalIdsRange(1, 1).equal(sortedTestChains[0 .. 1]));
    assert(sortedTestChains.equalIdsRange(2, 1).equal(sortedTestChains[3 .. 5]));
    assert(sortedTestChains.equalIdsRange(3, 1).equal(sortedTestChains[6 .. 7]));
    assert(sortedTestChains.equalIdsRange(42, 1337).equal(sortedTestChains[0 .. 0]));
}

/**
    Returns true iff ac1 begins before ac2 on the given contig (in).
*/
bool isBefore(string contig)(in AlignmentChain ac1, in AlignmentChain ac2) pure
        if (contig == "contigA" || contig == "contigB")
{
    assert(__traits(getMember, ac1, contig) == __traits(getMember, ac2, contig),
            "alignment chains do not belong to the same contig");

    static if (contig == "contigA")
    {
        return __traits(getMember, ac1.first, contig).begin < __traits(getMember,
                ac2.first, contig).begin;
    }
    else
    {

    }
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    enum complement = Flag.complement;

    auto acs = [
        AlignmentChain(0, Contig(1, 10), Contig(1, 10), Flags(), [LocalAlignment(Locus(1, 6), Locus(0, 1), 0)]),
        AlignmentChain(1, Contig(1, 10), Contig(2, 10), Flags(complement), [LocalAlignment(Locus(2, 6), Locus(0, 1), 0)]),
        AlignmentChain(2, Contig(1, 10), Contig(3, 10), Flags(), [LocalAlignment(Locus(3, 6), Locus(0, 1), 0)]),
        AlignmentChain(3, Contig(1, 10), Contig(4, 10), Flags(complement), [LocalAlignment(Locus(4, 6), Locus(0, 1), 0)]),
        AlignmentChain(4, Contig(1, 10), Contig(5, 10), Flags(), [LocalAlignment(Locus(5, 6), Locus(0, 1), 0)]),
    ];

    foreach (i; 0 .. acs.length)
        foreach (j; 0 .. acs.length)
        {
            auto compareValue = isBefore!"contigA"(acs[i], acs[j]);
            alias errorMessage = (expValue) => format!"expected isBefore!\"contigA\"(%s, %s) to be %s but got %s"(
                    acs[i], acs[j], expValue, compareValue);

            if (i < j)
                assert(compareValue, errorMessage(true));
            else
                assert(!compareValue, errorMessage(false));
        }
}

/// Calculates the coverage of the contigs by the given alignments. Only
/// contigs involved in the alignments are regarded.
double alignmentCoverage(in AlignmentChain[] alignments)
{
    static double coveredBases(T)(T alignmentsPerContig)
    {
        return alignmentsPerContig.map!(ac => ac.coveredBases!"contigA").sum;
    }

    auto alignmentsPerContig = alignments.chunkBy!"a.contigA.id == b.contigA.id";
    auto totalContigLength = alignmentsPerContig
        .save
        .map!"a.front.contigA.length"
        .sum;
    auto totalCoveredBases = alignmentsPerContig
        .save
        .map!coveredBases
        .sum;

    return totalCoveredBases.to!double / totalContigLength.to!double;
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    enum complement = Flag.complement;

    auto alignments = [
        AlignmentChain(
            0,
            Contig(1, 100),
            Contig(1, 50),
            Flags(),
            [
                LocalAlignment(
                    Locus(0, 10),
                    Locus(40, 50),
                    0
                ),
            ],
        ),
        AlignmentChain(
            1,
            Contig(1, 100),
            Contig(2, 30),
            Flags(complement),
            [
                LocalAlignment(
                    Locus(10, 20),
                    Locus(0, 10),
                    0
                ),
                LocalAlignment(
                    Locus(30, 40),
                    Locus(20, 30),
                    0
                ),
            ],
        ),
        AlignmentChain(
            2,
            Contig(1, 100),
            Contig(3, 20),
            Flags(),
            [
                LocalAlignment(
                    Locus(40, 60),
                    Locus(0, 20),
                    0
                ),
            ],
        ),
        AlignmentChain(
            3,
            Contig(1, 100),
            Contig(4, 50),
            Flags(complement),
            [
                LocalAlignment(
                    Locus(70, 100),
                    Locus(0, 30),
                    0
                ),
            ],
        ),
    ];

    assert(alignmentCoverage(alignments) == 80.0 / 100.0);
}


struct FlatLocalAlignment
{
    static struct FlatLocus
    {
        id_t id;
        coord_t length;
        coord_t begin;
        coord_t end;


        @property Contig contig() const pure nothrow @safe
        {
            return Contig(id, length);
        }


        @property void contig(Contig newContig) pure nothrow @safe
        {
            this.id = newContig.id;
            this.length = newContig.length;
        }


        @property Locus locus() const pure nothrow @safe
        {
            return Locus(begin, end);
        }


        @property void locus(Locus newLocus) pure nothrow @safe
        {
            this.begin = newLocus.begin;
            this.end = newLocus.end;
        }


        @property coord_t mappedLength() const pure nothrow @safe
        {
            return end - begin;
        }


        @property void boundedBegin(coord_t begin) pure nothrow @safe
        {
            this.begin = min(begin, length);
        }


        @property void boundedEnd(coord_t end) pure nothrow @safe
        {
            this.end = min(end, length);
        }


        bool beginsWithin(coord_t allowance) const pure nothrow @safe
        {
            return begin <= allowance;
        }


        bool endsWithin(coord_t allowance) const pure nothrow @safe
        {
            return end + allowance >= length;
        }


        bool isFullyContained(coord_t allowance) const pure nothrow @safe
        {
            return beginsWithin(allowance) && endsWithin(allowance);
        }
    }

    id_t id;
    FlatLocus contigA;
    FlatLocus contigB;
    Flags flags;
    trace_point_t tracePointDistance;
    TracePoint[] tracePoints;


    @property Trace trace() const pure nothrow @safe
    {
        return Trace(
            contigA.locus,
            contigB.locus,
            tracePointDistance,
            tracePoints,
        );
    }


    TranslatedTracePoint translateTracePoint(string contig = "contigA")(
        coord_t contigPos,
        RoundingMode roundingMode,
    ) const pure if (contig.among("contigA", "contigB"))
    {
        return trace.translateTracePoint!contig(contigPos, roundingMode);
    }


    /**
        Generate a cartoon of this alignment relative to `contig`.

        Params:
            bpsPerChar =    Number of base pairs that one char represents.

        Returns: a cartoon of this alignment
    */
    static string cartoon(string contig)(in coord_t bpsPerChar, in FlatLocalAlignment[] localAlignments...)
    {
        if (localAlignments.length == 0)
            return "";

        alias getContig = fla => mixin("fla." ~ contig ~ ".contig");
        alias getCoords = fla => mixin("fla." ~ contig ~ ".locus");

        auto referenceContig = getContig(localAlignments[0]);

        enforce!Exception(
            localAlignments.all!(ac => getContig(ac) == referenceContig),
            "all alignment chains must share the same reference contig",
        );

        auto cartoonLine(in FlatLocalAlignment fla)
        {
            auto skipBps = contig == "contigA" || !fla.flags.complement
                ? getCoords(fla).begin
                : referenceContig.length - getCoords(fla).end;

            return chain(
                ' '.repeat(skipBps / bpsPerChar),
                '-'.repeat(ceildiv(getCoords(fla).length, bpsPerChar)),
            );
        }

        return chain(
            '-'.repeat(ceildiv(referenceContig.length, bpsPerChar)),
            only('\n'),
            localAlignments.map!cartoonLine.joiner(only('\n')),
        ).to!string;
    }

    ///
    unittest
    {
        auto flas = [
            FlatLocalAlignment(
                0,
                FlatLocus(1, 10, 0, 3),
                FlatLocus(1, 10, 0, 3),
            ),
            FlatLocalAlignment(
                1,
                FlatLocus(1, 10, 4, 5),
                FlatLocus(1, 10, 4, 5),
            ),
            FlatLocalAlignment(
                2,
                FlatLocus(1, 10, 5, 8),
                FlatLocus(1, 10, 0, 3),
                Flags(Flag.complement),
            ),
            FlatLocalAlignment(
                3,
                FlatLocus(1, 10, 9, 10),
                FlatLocus(1, 10, 4, 5),
                Flags(Flag.complement),
            ),
        ];

        assert(cartoon!"contigA"(1, flas) == "----------\n" ~
                                             "---\n" ~
                                             "    -\n" ~
                                             "     ---\n" ~
                                             "         -");
        assert(cartoon!"contigB"(1, flas) == "----------\n" ~
                                             "---\n" ~
                                             "    -\n" ~
                                             "       ---\n" ~
                                             "     -");
    }
}


int cmpIdsAndComplement(ref const AlignmentChain lhs, ref const AlignmentChain rhs)
{
    return cmpLexicographically!(
        const(AlignmentChain),
        ac => ac.contigA.id,
        ac => ac.contigB.id,
        ac => ac.flags.complement,
    )(lhs, rhs);
}


/// Type of the read alignment.
enum AlignmentLocationSeed : ubyte
{
    front,
    back,
}

/**
    An alignment chain with a "seed", ie. hint for it's intended location.
*/
struct SeededAlignment
{
    AlignmentChain alignment;
    AlignmentLocationSeed seed;

    alias alignment this;

    invariant
    {
        assert(seed != AlignmentLocationSeed.front || isFrontExtension(alignment));
        assert(seed != AlignmentLocationSeed.back || isBackExtension(alignment));
    }

    int opCmp(ref const SeededAlignment other) const pure nothrow
    {
        return cmpLexicographically!(
            typeof(this),
            "a.alignment",
            "a.seed",
        )(this, other);
    }

    static InputRange!SeededAlignment from(AlignmentChain alignmentChain)
    {
        alias Seed = AlignmentLocationSeed;

        if (isFrontExtension(alignmentChain) && isBackExtension(alignmentChain))
        {
            return inputRangeObject(only(
                SeededAlignment(alignmentChain, Seed.front),
                SeededAlignment(alignmentChain, Seed.back),
            ));
        }
        else if (isFrontExtension(alignmentChain) && !isBackExtension(alignmentChain))
        {
            return inputRangeObject(only(SeededAlignment(alignmentChain, Seed.front)));
        }
        else if (isBackExtension(alignmentChain) && !isFrontExtension(alignmentChain))
        {
            return inputRangeObject(only(SeededAlignment(alignmentChain, Seed.back)));
        }
        else
        {
            return inputRangeObject(takeNone!(SeededAlignment[]));
        }
    }
}

/// Type of the read alignment.
static enum ReadAlignmentType
{
    front = 0,
    gap = 1,
    back = 2,
}

private bool isExtension(in AlignmentChain alignment) pure nothrow
{
    return isFrontExtension(alignment) ^ isBackExtension(alignment);
}

private bool isFrontExtension(in AlignmentChain alignment) pure nothrow
{
    auto readExtensionLength = alignment.first.contigB.begin;
    auto referenceExtensionLength = alignment.first.contigA.begin;

    return readExtensionLength > referenceExtensionLength;
}

private bool isBackExtension(in AlignmentChain alignment) pure nothrow
{
    auto readExtensionLength = alignment.contigB.length - alignment.last.contigB.end;
    auto referenceExtensionLength = alignment.contigA.length - alignment.last.contigA.end;

    return readExtensionLength > referenceExtensionLength;
}

/**
    Alignment of a read against the reference. This is either one or two
    alignment chains which belong to the same read and one or two reference
    contig(s).
*/
struct ReadAlignment
{
    SeededAlignment[2] _alignments;
    size_t _length;

    this(SeededAlignment[] alignments...)
    {
        if (1 <= alignments.length && alignments.length <= 2)
        {
            this._length = alignments.length;
            this._alignments[0 .. _length] = alignments[0 .. _length];
        }
        else
        {
            logJsonDebug(
                "info", format!"creating invalid read alignment with %d local alignments"(alignments.length),
                "alignments", alignments.toJson,
            );

            this._length = 0;
        }
    }

    @property size_t length() pure const nothrow
    {
        return _length;
    }

    @property size_t opDollar() pure const nothrow
    {
        return _length;
    }

    inout(SeededAlignment[]) opIndex() inout pure nothrow
    {
        return _alignments[0 .. _length];
    }

    inout(SeededAlignment) opIndex(T)(T idx) inout pure nothrow
    {
        return this[][idx];
    }

    unittest
    {
        auto ac1 = SeededAlignment(AlignmentChain(1), AlignmentLocationSeed.front);
        auto ac2 = SeededAlignment(AlignmentChain(2), AlignmentLocationSeed.back);

        auto ra1 = ReadAlignment(ac1);

        assert(ra1.length == 1);
        assert(ra1[] == [ac1]);
        assert(ra1[0] == ac1);
        assert(ra1[$ - 1] == ac1);

        auto ra2 = ReadAlignment(ac1, ac2);

        assert(ra2.length == 2);
        assert(ra2[] == [ac1, ac2]);
        assert(ra2[0] == ac1);
        assert(ra2[1] == ac2);
        assert(ra2[$ - 1] == ac2);
    }

    /**
        If readAlignment is a gap return true iff the first alignment's
        `contigA.id` is lower than the second alignment's; otherwise returns
        true.
    */
    @property bool isInOrder() const pure nothrow
    {
        return !isGap || _alignments[0].contigA.id < _alignments[1].contigA.id;
    }

    ReadAlignment getInOrder() pure nothrow
    {
        if (isGap && !isInOrder)
        {
            swap(_alignments[0], _alignments[1]);
        }

        return this;
    }

    /**
        Returns true iff the read alignment is valid, ie. it is either an
        extension or gap.
    */
    @property bool isValid() const pure nothrow
    {
        return isExtension ^ isGap;
    }

    /**
        Get the type of the read alignment.

        See_Also: `isFrontExtension`, `isBackExtension`, `isGap`
    */
    @property ReadAlignmentType type() const pure nothrow
    {
        assert(isValid, "invalid read alignment");

        if (isGap)
        {
            return ReadAlignmentType.gap;
        }
        else if (isFrontExtension)
        {
            return ReadAlignmentType.front;
        }
        else
        {
            return ReadAlignmentType.back;
        }
    }

    /**
        Returns true iff the read alignment is an extension, ie. it is a front or
        back extension.

        See_Also: `isFrontExtension`, `isBackExtension`
    */
    @property bool isExtension() const pure nothrow
    {
        // A single SeededAlignment is always a valid extension.
        return _length == 1;
    }

    /**
        Returns true iff the read alignment is an front extension, ie. it is an
        extension and reaches over the front of the reference contig.

        ---
        Case 1 (complement alignment):

                          0  rx
            ref           |--+->-+->-->-->--|
                             | | |
            read  |--<--<--<-+-<-+--|
                  0          ax

        Case 2 (non-complement alignment):

                          0  rx
            ref           |--+->-+->-->-->--|
                             | | |
            read  |-->-->-->-+->-+--|
                  0          ax
        ---
    */
    @property bool isFrontExtension() const pure nothrow
    {
        return _length == 1 && _alignments[0].seed == AlignmentLocationSeed.front;
    }

    /**
        Returns true iff the read alignment is an back extension, ie. it is an
        extension and reaches over the back of the reference contig.

        ---
        Case 1 (complement alignment):

                  0             ry lr
            ref   |-->-->-->-+->-+--|
                             | | |
            read          |--+-<-+-<--<--<--|
                          0     ay         la

        Case 2 (non-complement alignment):

                  0             ry lr
            ref   |-->-->-->-+->-+--|
                             | | |
            read          |--+->-+->-->-->--|
                          0     ay         la
        ---
    */
    @property bool isBackExtension() const pure nothrow
    {
        return _length == 1 && _alignments[0].seed == AlignmentLocationSeed.back;
    }

    /**
        Returns true iff the read alignment spans a gap, ie. two alignments of
        the same read on different reference contigs are involved.
    */
    @property bool isGap() const pure nothrow
    {
        return _length == 2 &&
            _alignments[0].contigA.id != _alignments[1].contigA.id &&
            _alignments[0].contigB.id == _alignments[1].contigB.id;
    }

    /**
        Returns true iff the read alignment spans a gap and the flanking
        contigs have in the same orientation (according to this alignment).
    */
    @property bool isParallel() const pure nothrow
    {
        return isGap &&
            _alignments[0].seed != _alignments[1].seed &&
            _alignments[0].complement == _alignments[1].complement;
    }

    /**
        Returns true iff the read alignment spans a gap and the flanking
        contigs have in different orientation (according to this alignment).
    */
    @property bool isAntiParallel() const pure nothrow
    {
        return isGap &&
            _alignments[0].seed == _alignments[1].seed &&
            _alignments[0].complement != _alignments[1].complement;
    }

    unittest
    {
        alias LocalAlignment = AlignmentChain.LocalAlignment;
        alias complement = Flag.complement;

        auto testData = [
            // TODO drop this case?
            //"innerAlignment": ReadAlignment(
            //    SeededAlignment.from(AlignmentChain(
            //        1,
            //        Contig(1, 100),
            //        Contig(1, 10),
            //        no,
            //        [
            //            LocalAlignment(
            //                Locus(10, 11),
            //                Locus(0, 1),
            //                0,
            //            ),
            //            LocalAlignment(
            //                Locus(19, 20),
            //                Locus(9, 10),
            //                0,
            //            ),
            //        ],
            //    )).front,
            //),
            // TODO drop this case?
            //"innerAlignmentComplement": ReadAlignment(
            //    SeededAlignment.from(AlignmentChain(
            //        2,
            //        Contig(1, 100),
            //        Contig(1, 10),
            //        yes,
            //        [
            //            LocalAlignment(
            //                Locus(10, 11),
            //                Locus(0, 1),
            //                0,
            //            ),
            //            LocalAlignment(
            //                Locus(19, 20),
            //                Locus(9, 10),
            //                0,
            //            ),
            //        ],
            //    )).front,
            //),
            "frontExtension": ReadAlignment(
                SeededAlignment.from(AlignmentChain(
                    3,
                    Contig(1, 100),
                    Contig(1, 10),
                    Flags(),
                    [
                        LocalAlignment(
                            Locus(2, 3),
                            Locus(5, 6),
                            0,
                        ),
                        LocalAlignment(
                            Locus(5, 6),
                            Locus(9, 10),
                            0,
                        ),
                    ],
                )).front,
            ),
            "frontExtensionComplement": ReadAlignment(
                SeededAlignment.from(AlignmentChain(
                    4,
                    Contig(1, 100),
                    Contig(1, 10),
                    Flags(complement),
                    [
                        LocalAlignment(
                            Locus(2, 3),
                            Locus(5, 6),
                            0,
                        ),
                        LocalAlignment(
                            Locus(5, 6),
                            Locus(9, 10),
                            0,
                        ),
                    ],
                )).front,
            ),
            "backExtension": ReadAlignment(
                SeededAlignment.from(AlignmentChain(
                    5,
                    Contig(1, 100),
                    Contig(1, 10),
                    Flags(),
                    [
                        LocalAlignment(
                            Locus(94, 95),
                            Locus(0, 1),
                            0,
                        ),
                        LocalAlignment(
                            Locus(97, 98),
                            Locus(4, 5),
                            0,
                        ),
                    ],
                )).front,
            ),
            "backExtensionComplement": ReadAlignment(
                SeededAlignment.from(AlignmentChain(
                    6,
                    Contig(1, 100),
                    Contig(1, 10),
                    Flags(complement),
                    [
                        LocalAlignment(
                            Locus(94, 95),
                            Locus(0, 1),
                            0,
                        ),
                        LocalAlignment(
                            Locus(97, 98),
                            Locus(4, 5),
                            0,
                        ),
                    ],
                )).front,
            ),
            "gapEnd2Front": ReadAlignment(
                SeededAlignment.from(AlignmentChain(
                    7,
                    Contig(1, 100),
                    Contig(1, 10),
                    Flags(),
                    [
                        LocalAlignment(
                            Locus(94, 95),
                            Locus(0, 1),
                            0,
                        ),
                        LocalAlignment(
                            Locus(97, 98),
                            Locus(4, 5),
                            0,
                        ),
                    ],
                )).front,
                SeededAlignment.from(AlignmentChain(
                    8,
                    Contig(2, 100),
                    Contig(1, 10),
                    Flags(),
                    [
                        LocalAlignment(
                            Locus(2, 3),
                            Locus(5, 6),
                            0,
                        ),
                        LocalAlignment(
                            Locus(5, 6),
                            Locus(9, 10),
                            0,
                        ),
                    ],
                )).front,
            ),
            "gapFront2EndComplement": ReadAlignment(
                SeededAlignment.from(AlignmentChain(
                    9,
                    Contig(2, 100),
                    Contig(1, 10),
                    Flags(complement),
                    [
                        LocalAlignment(
                            Locus(2, 3),
                            Locus(5, 6),
                            0,
                        ),
                        LocalAlignment(
                            Locus(5, 6),
                            Locus(9, 10),
                            0,
                        ),
                    ],
                )).front,
                SeededAlignment.from(AlignmentChain(
                    10,
                    Contig(1, 100),
                    Contig(1, 10),
                    Flags(complement),
                    [
                        LocalAlignment(
                            Locus(94, 95),
                            Locus(0, 1),
                            0,
                        ),
                        LocalAlignment(
                            Locus(97, 98),
                            Locus(4, 5),
                            0,
                        ),
                    ],
                )).front,
            ),
            "gapEnd2End": ReadAlignment(
                SeededAlignment.from(AlignmentChain(
                    11,
                    Contig(1, 100),
                    Contig(1, 10),
                    Flags(complement),
                    [
                        LocalAlignment(
                            Locus(94, 95),
                            Locus(0, 1),
                            0,
                        ),
                        LocalAlignment(
                            Locus(97, 98),
                            Locus(4, 5),
                            0,
                        ),
                    ],
                )).front,
                SeededAlignment.from(AlignmentChain(
                    12,
                    Contig(2, 100),
                    Contig(1, 10),
                    Flags(),
                    [
                        LocalAlignment(
                            Locus(94, 95),
                            Locus(0, 1),
                            0,
                        ),
                        LocalAlignment(
                            Locus(97, 98),
                            Locus(4, 5),
                            0,
                        ),
                    ],
                )).front,
            ),
            "gapBegin2Begin": ReadAlignment(
                SeededAlignment.from(AlignmentChain(
                    13,
                    Contig(2, 100),
                    Contig(1, 10),
                    Flags(),
                    [
                        LocalAlignment(
                            Locus(2, 3),
                            Locus(5, 6),
                            0,
                        ),
                        LocalAlignment(
                            Locus(5, 6),
                            Locus(9, 10),
                            0,
                        ),
                    ],
                )).front,
                SeededAlignment.from(AlignmentChain(
                    14,
                    Contig(1, 100),
                    Contig(1, 10),
                    Flags(complement),
                    [
                        LocalAlignment(
                            Locus(2, 3),
                            Locus(5, 6),
                            0,
                        ),
                        LocalAlignment(
                            Locus(5, 6),
                            Locus(9, 10),
                            0,
                        ),
                    ],
                )).front,
            ),
        ];
        //                               +--------- isInOrder
        //                               |+-------- isValid
        //                               ||+------- type
        //                               |||+------ isExtension
        //                               ||||+----- isFrontExt
        //                               |||||+---- isBackExt
        //                               ||||||+--- isGap
        //                               |||||||+-- isParallel
        //                               ||||||||+- isAntiParallel
        //                               |||||||||
        //                               |||||||||
        auto testCases = [
            //"innerAlignment":           "+.F......",
            //"innerAlignmentComplement": "+.F......",
            "frontExtension":           "++F++....",
            "frontExtensionComplement": "++F++....",
            "backExtension":            "++B+.+...",
            "backExtensionComplement":  "++B+.+...",
            "gapEnd2Front":             "++G...++.",
            "gapFront2EndComplement":   ".+G...++.",
            "gapEnd2End":               "++G...+.+",
            "gapBegin2Begin":           ".+G...+.+",
        ];

        alias getFailureMessage = (testCase, testFunction, expectedValue) => format!"expected %s.%s to be %s"(
                testCase, testFunction, expectedValue);
        alias toBool = (c) => c == '+' ? true : false;
        alias toRAT = (c) => c == 'F'
            ? ReadAlignmentType.front
            : c == 'B'
                ? ReadAlignmentType.back
                : ReadAlignmentType.gap;

        foreach (testCase, expectations; testCases)
        {
            auto readAlignment = testData[testCase];

            auto expIsInOrder = toBool(expectations[0]);
            auto expIsValid = toBool(expectations[1]);
            auto expType = toRAT(expectations[2]);
            auto expIsExtension = toBool(expectations[3]);
            auto expIsFrontExt = toBool(expectations[4]);
            auto expIsBackExt = toBool(expectations[5]);
            auto expIsGap = toBool(expectations[6]);
            auto expIsParallel = toBool(expectations[7]);
            auto expIsAntiParallel = toBool(expectations[8]);

            assert(expIsInOrder == readAlignment.isInOrder,
                    getFailureMessage(testCase, "isInOrder", expIsInOrder));
            assert(expIsValid == readAlignment.isValid,
                    getFailureMessage(testCase, "isValid", expIsValid));
            if (readAlignment.isValid)
                assert(expType == readAlignment.type,
                        getFailureMessage(testCase, "type", expType));
            else
                assertThrown!AssertError(readAlignment.type,
                        format!"expected type(%s) to throw"(testCase));
            assert(expIsExtension == readAlignment.isExtension,
                    getFailureMessage(testCase, "isExtension", expIsExtension));
            assert(expIsFrontExt == readAlignment.isFrontExtension,
                    getFailureMessage(testCase, "isFrontExtension", expIsFrontExt));
            assert(expIsBackExt == readAlignment.isBackExtension,
                    getFailureMessage(testCase, "isBackExtension", expIsBackExt));
            assert(expIsGap == readAlignment.isGap,
                    getFailureMessage(testCase, "isGap", expIsGap));
            assert(expIsParallel == readAlignment.isParallel,
                    getFailureMessage(testCase, "isParallel", expIsParallel));
            assert(expIsAntiParallel == readAlignment.isAntiParallel,
                    getFailureMessage(testCase, "isAntiParallel", expIsAntiParallel));
        }
    }

    double meanScore() const pure
    {
        return this[].map!"a.score".mean;
    }
}

/// Generate basic join from read alignment.
J makeJoin(J)(ReadAlignment readAlignment)
{
    final switch (readAlignment.type)
    {
    case ReadAlignmentType.front:
        return J(
            ContigNode(
                readAlignment[0].contigA.id,
                ContigPart.pre,
            ),
            ContigNode(
                readAlignment[0].contigA.id,
                ContigPart.begin,
            ),
        );
    case ReadAlignmentType.gap:
        alias contigPart = (locationSeed) => locationSeed == AlignmentLocationSeed.front
                    ? ContigPart.begin
                    : ContigPart.end;

        return J(
            ContigNode(
                readAlignment[0].contigA.id,
                contigPart(readAlignment[0].seed),
            ),
            ContigNode(
                readAlignment[1].contigA.id,
                contigPart(readAlignment[1].seed),
            ),
        );
    case ReadAlignmentType.back:
        return J(
            ContigNode(
                readAlignment[0].contigA.id,
                ContigPart.end,
            ),
            ContigNode(
                readAlignment[0].contigA.id,
                ContigPart.post,
            ),
        );
    }
}

/**
    A pile of read alignments belonging to the same gap/contig end.

    See_Also: Hit
*/
alias PileUp = ReadAlignment[];

/**
    Get the type of the read alignment.

    See_Also: `isFrontExtension`, `isBackExtension`, `isGap`
*/
ReadAlignmentType getType(in PileUp pileUp) pure nothrow
{
    assert(pileUp.isValid, "invalid read alignment");

    if (pileUp.isGap)
    {
        return ReadAlignmentType.gap;
    }
    else
    {
        assert(pileUp.isExtension);
        return pileUp[0].isFrontExtension
            ? ReadAlignmentType.front
            : ReadAlignmentType.back;
    }
}

bool isValid(in PileUp pileUp) pure nothrow
{
    return pileUp.isExtension ^ pileUp.isGap;
}

bool isExtension(in PileUp pileUp) pure nothrow
{
    if (pileUp.length > 0 && pileUp[0].isFrontExtension)
    {
        return pileUp.all!(readAlignment => readAlignment.isFrontExtension);
    }
    else if (pileUp.length > 0 && pileUp[0].isBackExtension)
    {
        return pileUp.all!(readAlignment => readAlignment.isBackExtension);
    }
    else
    {
        return false;
    }
}

bool isGap(in PileUp pileUp) pure nothrow
{
    return pileUp.any!(readAlignment => readAlignment.isGap);
}

auto isParallel(in PileUp pileUp) pure nothrow
{
    return pileUp.isGap && pileUp
        .filter!(readAlignment => readAlignment.isGap)
        .front
        .isParallel;
}

auto isAntiParallel(in PileUp pileUp) pure nothrow
{
    return pileUp.isGap && pileUp
        .filter!(readAlignment => readAlignment.isGap)
        .front
        .isAntiParallel;
}

Contig[] contigs(in PileUp pileUp) nothrow
{
    Contig[] contigs;
    contigs.reserve(2);

    foreach (readAlignment; pileUp)
    {
        if (contigs.length == 2)
            break;
        else if (contigs.length == 0)
            foreach (ac; readAlignment[])
                contigs ~= ac.contigA;
        else
            foreach (ac; readAlignment[])
                if (contigs[0] != ac.contigA)
                    contigs ~= ac.contigA;
    }

    return contigs;
}

/// Returns a list of pointers to all involved alignment chains.
AlignmentChain*[] getAlignmentRefs(PileUp pileUp) pure nothrow
{
    auto alignmentChainsAcc = appender!(AlignmentChain*[]);
    alignmentChainsAcc.reserve(pileUp.map!"a.length".sum);

    foreach (ref readAlignment; pileUp)
    {
        foreach (ref seededAlignment; readAlignment)
        {
            alignmentChainsAcc ~= &(seededAlignment.alignment);
        }
    }

    return alignmentChainsAcc.data;
}

///
unittest
{
    auto pileUp = [
        ReadAlignment(SeededAlignment(), SeededAlignment()), ReadAlignment(SeededAlignment())
    ];
    auto allAlignmentChains = pileUp.getAlignmentRefs();

    assert(allAlignmentChains.length == 3);
    assert(pileUp[0][0].id == 0);
    assert(pileUp[0][1].id == 0);
    assert(pileUp[1][0].id == 0);

    allAlignmentChains[0].id = 1;
    allAlignmentChains[1].id = 2;
    allAlignmentChains[2].id = 3;

    assert(pileUp[0][0].id == 1);
    assert(pileUp[0][1].id == 2);
    assert(pileUp[1][0].id == 3);
}

/// Converts the pileup into a simple JSON object for diagnostic purposes.
Json pileUpToSimpleJson(in PileUp pileUp)
{
    return [
        "type": pileUp.getType.to!string.toJson,
        "length": pileUp.length.toJson,
        "contigIds": pileUp
            .map!(ra => ra[].map!"0 + a.contigA.id".array)
            .joiner
            .array
            .sort
            .uniq
            .array
            .toJson,
    ].toJson;
}

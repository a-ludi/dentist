/**
    Everything to handle local alignments and friends.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.alignments;

import dentist.common.scaffold : buildScaffold, concatenatePayloads, ContigNode,
    ContigPart, discardAmbiguousJoins, Join, mergeExtensionsWithGaps;
import dentist.util.algorithm : cmpLexicographically, orderLexicographically;
import dentist.util.log;
import core.exception : AssertError;
import std.algorithm : all, any, canFind, chunkBy, equal, filter, isSorted,
    joiner, map, mean, min, sort, sum, swap, SwapStrategy;
import std.array : appender, array;
import std.conv : to;
import std.exception : assertNotThrown, assertThrown;
import std.format : format;
import std.math : sgn;
import std.range : assumeSorted, chain, chunks, InputRange, inputRangeObject,
    iota, only, slide, takeNone, zip;
import std.string : capitalize;
import std.typecons : BitFlags, PhobosFlag = Flag, No, Yes;
import vibe.data.json : toJson = serializeToJson;

debug import std.stdio : writeln;

alias arithmetic_t = int;
alias coord_t = uint;
alias diff_t = uint;
alias id_t = uint;
alias trace_point_t = ushort;

/**
    Holds a chain of local alignments that form a compound alignment. An AlignmentChain should
    contain at least one element.
*/
struct AlignmentChain
{
    static enum Flag : ubyte
    {
        complement = 1 << 0,
        disabled = 1 << 1,
    }

    static alias Flags = BitFlags!Flag;

    static struct LocalAlignment
    {
        static struct Locus
        {
            coord_t begin;
            coord_t end;
        }

        static struct TracePoint
        {
            trace_point_t numDiffs;
            trace_point_t numBasePairs;
        }

        Locus contigA;
        Locus contigB;
        diff_t numDiffs;
        TracePoint[] tracePoints;
    }

    static struct Contig
    {
        id_t id;
        coord_t length;
    }

    static immutable maxScore = 2 ^^ 16;
    static immutable Flags emptyFlags;

    id_t id;
    Contig contigA;
    Contig contigB;
    Flags flags;
    LocalAlignment[] localAlignments;
    trace_point_t tracePointDistance;

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
                    // dfmt off
                    debug logJsonDebug(
                        "contigA", contigA.id,
                        "contigB", contigB.id,
                        "la.numDiffs", la.numDiffs,
                        "traceDiffs", traceDiffs,
                    );
                    // dfmt on
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
                auto acZeroLength = AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, []);
                auto ac1 = AlignmentChain(1, Contig(1, 10), Contig(1, 10), emptyFlags,
                        [LocalAlignment(Locus(1, 1), Locus(1, 9), 0)]);
                auto ac2 = AlignmentChain(2, Contig(1, 10), Contig(1, 10), emptyFlags,
                        [LocalAlignment(Locus(1, 11), Locus(1, 9), 0)]);
                auto ac3 = AlignmentChain(3, Contig(1, 10), Contig(1, 10), emptyFlags,
                        [LocalAlignment(Locus(1, 9), Locus(1, 1), 0)]);
                auto ac4 = AlignmentChain(4, Contig(1, 10), Contig(1, 10), emptyFlags,
                        [LocalAlignment(Locus(1, 9), Locus(1, 11), 0)]);
                auto acFine = AlignmentChain(5, Contig(1, 10), Contig(1, 10),
                        emptyFlags, [LocalAlignment(Locus(1, 9), Locus(1, 9), 0)]);

                assertThrown!AssertError(acZeroLength.totalLength);
                assertThrown!AssertError(ac1.totalLength);
                assertThrown!AssertError(ac2.totalLength);
                assertThrown!AssertError(ac3.totalLength);
                assertThrown!AssertError(ac4.totalLength);
                assertNotThrown!AssertError(acFine.totalLength);
            }
    }

    @property ref const(LocalAlignment) first() const pure nothrow
    {
        return localAlignments[0];
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto firstLA = LocalAlignment(Locus(1, 9), Locus(1, 9), 0);
                auto otherLA = LocalAlignment(Locus(9, 10), Locus(9, 10), 0);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [firstLA, otherLA]);

                assert(ac.first == firstLA);
            }
    }

    @property ref const(LocalAlignment) last() const pure nothrow
    {
        return localAlignments[$ - 1];
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto lastLA = LocalAlignment(Locus(1, 9), Locus(1, 9), 0);
                auto otherLA = LocalAlignment(Locus(9, 10), Locus(9, 10), 0);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [otherLA, lastLA]);

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
    @property bool isProper() const pure nothrow
    {
        // dfmt off
        return (
            first.contigA.begin == 0 ||
            first.contigB.begin == 0
        )
        &&
        (
            last.contigA.end == contigA.length ||
            last.contigB.end == contigB.length
        );
        // dfmt on
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
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), emptyFlags, [la]);

                // read with extension align an contigA from 25 to 40
                assert(ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(10, 20), Locus(5, 10), 1);
                auto la2 = LocalAlignment(Locus(30, 40), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), emptyFlags, [la1, la2]);

                // read with extension align an contigA from 5 to 45
                assert(ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(0, 10), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), emptyFlags, [la]);

                // read with extension align an contigA from -5 to 15
                assert(!ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(40, 50), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), emptyFlags, [la]);

                // read with extension align an contigA from 35 to 55
                assert(!ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(0, 20), Locus(5, 10), 1);
                auto la2 = LocalAlignment(Locus(30, 50), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), emptyFlags, [la1, la2]);

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
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [la1, la2]);

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
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [la1, la2]);

                assert(ac.totalDiffs == 3);
            }
    }

    @property size_t totalGapLength() const pure
    {
        // dfmt off
        return localAlignments
            .chunks(2)
            .map!(las => las.length < 2 ? 0 : las[1].contigA.begin - las[0].contigA.end)
            .sum;
        // dfmt on
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [la1, la2]);

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
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [la1, la2]);

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
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [la1, la2]);

                assert(ac.score == 4 * maxScore / 9);
            }
    }

    int compareIds(ref const AlignmentChain other) const pure nothrow
    {
        // dfmt off
        return cmpLexicographically!(
            typeof(this),
            ac => ac.contigA.id,
            ac => ac.contigB.id,
        )(this, other);
        // dfmt on
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(0, 1), Locus(0, 1), 1);
                // dfmt off
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [la]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 10), emptyFlags, [la]),
                    AlignmentChain(2, Contig(2, 10), Contig(1, 10), emptyFlags, [la]),
                    AlignmentChain(3, Contig(2, 10), Contig(2, 10), emptyFlags, [la]),
                ];
                // dfmt on

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
        // dfmt off
        return cmpLexicographically!(
            typeof(this),
            ac => ac.contigA.id,
            ac => ac.contigB.id,
            ac => ac.first.contigA.begin,
            ac => ac.first.contigB.begin,
            ac => ac.last.contigA.end,
            ac => ac.last.contigB.end,
        )(this, other);
        // dfmt on
    }

    unittest
    {
        // see compareIds
        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(0, 1), Locus(0, 1), 1);
                // dfmt off
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [la]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 10), emptyFlags, [la]),
                    AlignmentChain(2, Contig(2, 10), Contig(1, 10), emptyFlags, [la]),
                    AlignmentChain(3, Contig(2, 10), Contig(2, 10), emptyFlags, [la]),
                ];
                // dfmt on

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
                // dfmt off
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [
                        LocalAlignment(Locus(0, 2), Locus(0, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(1, Contig(1, 10), Contig(1, 10), emptyFlags, [
                        LocalAlignment(Locus(0, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(2, Contig(1, 10), Contig(1, 10), emptyFlags, [
                        LocalAlignment(Locus(1, 2), Locus(0, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(3, Contig(1, 10), Contig(1, 10), emptyFlags, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(4, Contig(1, 10), Contig(1, 10), emptyFlags, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 6), 1)
                    ]),
                    AlignmentChain(5, Contig(1, 10), Contig(1, 10), emptyFlags, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 6), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(6, Contig(1, 10), Contig(1, 10), emptyFlags, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 6), Locus(3, 6), 1)
                    ]),
                ];
                // dfmt on

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
}

bool idsPred(in AlignmentChain ac1, in AlignmentChain ac2) pure
{
    auto cmpValue = ac1.compareIds(ac2);

    return 0 != cmpValue && cmpValue < 0;
}

unittest
{
    with (AlignmentChain) with (LocalAlignment) with (Complement)
            {
                auto la = LocalAlignment(Locus(0, 1), Locus(0, 1), 1);
                // dfmt off
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [la]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 10), emptyFlags, [la]),
                    AlignmentChain(2, Contig(2, 10), Contig(1, 10), emptyFlags, [la]),
                    AlignmentChain(3, Contig(2, 10), Contig(2, 10), emptyFlags, [la]),
                ];
                // dfmt on

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
}

auto haveEqualIds(in AlignmentChain ac1, in AlignmentChain ac2) pure
{
    auto cmpValue = ac1.compareIds(ac2);

    return 0 == cmpValue;
}

unittest
{
    with (AlignmentChain) with (Complement) with (LocalAlignment)
            {
                // dfmt off
                auto sortedTestChains = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [LocalAlignment(Locus(0, 10), Locus(0, 1), 0)]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 20), emptyFlags, [LocalAlignment(Locus(0, 10), Locus(0, 2), 0)]),
                    AlignmentChain(2, Contig(1, 10), Contig(3, 30), emptyFlags, [LocalAlignment(Locus(0, 10), Locus(0, 3), 0)]),
                    AlignmentChain(3, Contig(2, 20), Contig(1, 10), emptyFlags, [LocalAlignment(Locus(0, 20), Locus(0, 4), 0)]),
                    AlignmentChain(4, Contig(2, 20), Contig(1, 10), emptyFlags, [LocalAlignment(Locus(0, 20), Locus(0, 5), 0)]),
                    AlignmentChain(5, Contig(2, 20), Contig(3, 30), emptyFlags, [LocalAlignment(Locus(0, 20), Locus(0, 6), 0)]),
                    AlignmentChain(6, Contig(3, 30), Contig(1, 10), emptyFlags, [LocalAlignment(Locus(0, 30), Locus(0, 7), 0)]),
                    AlignmentChain(7, Contig(3, 30), Contig(2, 20), emptyFlags, [LocalAlignment(Locus(0, 30), Locus(0, 8), 0)]),
                    AlignmentChain(8, Contig(3, 30), Contig(3, 30), emptyFlags, [LocalAlignment(Locus(0, 30), Locus(0, 9), 0)]),
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
                // dfmt on
            }
}

auto equalIdsRange(in AlignmentChain[] acList, in id_t contigAID, in id_t contigBID) pure
{
    assert(isSorted!idsPred(acList));
    // dfmt off
    AlignmentChain needle = {
        contigA: AlignmentChain.Contig(contigAID, 1),
        contigB: AlignmentChain.Contig(contigBID, 1),
        localAlignments: [AlignmentChain.LocalAlignment(AlignmentChain.LocalAlignment.Locus(0, 1), AlignmentChain.LocalAlignment.Locus(0, 1))],
    };
    // dfmt on

    return acList.assumeSorted!idsPred.equalRange(needle);
}

unittest
{
    with (AlignmentChain) with (Complement) with (LocalAlignment)
            {
                // dfmt off
                auto sortedTestChains = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [LocalAlignment(Locus(0, 1), Locus(0, 1), 0)]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 20), emptyFlags, [LocalAlignment(Locus(0, 2), Locus(0, 2), 0)]),
                    AlignmentChain(2, Contig(1, 10), Contig(3, 30), emptyFlags, [LocalAlignment(Locus(0, 3), Locus(0, 3), 0)]),
                    AlignmentChain(3, Contig(2, 20), Contig(1, 10), emptyFlags, [LocalAlignment(Locus(0, 4), Locus(0, 4), 0)]),
                    AlignmentChain(4, Contig(2, 20), Contig(1, 10), emptyFlags, [LocalAlignment(Locus(0, 5), Locus(0, 5), 0)]),
                    AlignmentChain(5, Contig(2, 20), Contig(3, 30), emptyFlags, [LocalAlignment(Locus(0, 6), Locus(0, 6), 0)]),
                    AlignmentChain(6, Contig(3, 30), Contig(1, 10), emptyFlags, [LocalAlignment(Locus(0, 7), Locus(0, 7), 0)]),
                    AlignmentChain(7, Contig(3, 30), Contig(2, 20), emptyFlags, [LocalAlignment(Locus(0, 8), Locus(0, 8), 0)]),
                    AlignmentChain(8, Contig(3, 30), Contig(3, 30), emptyFlags, [LocalAlignment(Locus(0, 9), Locus(0, 9), 0)]),
                ];
                // dfmt on

                assert(sortedTestChains.equalIdsRange(1, 1).equal(sortedTestChains[0 .. 1]));
                assert(sortedTestChains.equalIdsRange(2, 1).equal(sortedTestChains[3 .. 5]));
                assert(sortedTestChains.equalIdsRange(3, 1).equal(sortedTestChains[6 .. 7]));
                assert(sortedTestChains.equalIdsRange(42, 1337).equal(sortedTestChains[0 .. 0]));
            }
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
    with (AlignmentChain) with (Flag) with (LocalAlignment)
            {
                // dfmt off
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), emptyFlags, [LocalAlignment(Locus(1, 6), Locus(0, 1), 0)]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 10), Flags(complement), [LocalAlignment(Locus(2, 6), Locus(0, 1), 0)]),
                    AlignmentChain(2, Contig(1, 10), Contig(3, 10), emptyFlags, [LocalAlignment(Locus(3, 6), Locus(0, 1), 0)]),
                    AlignmentChain(3, Contig(1, 10), Contig(4, 10), Flags(complement), [LocalAlignment(Locus(4, 6), Locus(0, 1), 0)]),
                    AlignmentChain(4, Contig(1, 10), Contig(5, 10), emptyFlags, [LocalAlignment(Locus(5, 6), Locus(0, 1), 0)]),
                ];
                // dfmt on

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
    // dfmt off
    auto totalContigLength = alignmentsPerContig
        .save
        .map!"a.front.contigA.length"
        .sum;
    auto totalCoveredBases = alignmentsPerContig
        .save
        .map!coveredBases
        .sum;
    // dfmt on

    return totalCoveredBases.to!double / totalContigLength.to!double;
}

unittest
{
    // dfmt off
    auto alignments = [
        AlignmentChain(
            0,
            AlignmentChain.Contig(1, 100),
            AlignmentChain.Contig(1, 50),
            AlignmentChain.emptyFlags,
            [
                AlignmentChain.LocalAlignment(
                    AlignmentChain.LocalAlignment.Locus(0, 10),
                    AlignmentChain.LocalAlignment.Locus(40, 50),
                    0
                ),
            ],
        ),
        AlignmentChain(
            1,
            AlignmentChain.Contig(1, 100),
            AlignmentChain.Contig(2, 30),
            AlignmentChain.Flags(AlignmentChain.Flag.complement),
            [
                AlignmentChain.LocalAlignment(
                    AlignmentChain.LocalAlignment.Locus(10, 20),
                    AlignmentChain.LocalAlignment.Locus(0, 10),
                    0
                ),
                AlignmentChain.LocalAlignment(
                    AlignmentChain.LocalAlignment.Locus(30, 40),
                    AlignmentChain.LocalAlignment.Locus(20, 30),
                    0
                ),
            ],
        ),
        AlignmentChain(
            2,
            AlignmentChain.Contig(1, 100),
            AlignmentChain.Contig(3, 20),
            AlignmentChain.emptyFlags,
            [
                AlignmentChain.LocalAlignment(
                    AlignmentChain.LocalAlignment.Locus(40, 60),
                    AlignmentChain.LocalAlignment.Locus(0, 20),
                    0
                ),
            ],
        ),
        AlignmentChain(
            3,
            AlignmentChain.Contig(1, 100),
            AlignmentChain.Contig(4, 50),
            AlignmentChain.Flags(AlignmentChain.Flag.complement),
            [
                AlignmentChain.LocalAlignment(
                    AlignmentChain.LocalAlignment.Locus(70, 100),
                    AlignmentChain.LocalAlignment.Locus(0, 30),
                    0
                ),
            ],
        ),
    ];
    // dfmt on

    assert(alignmentCoverage(alignments) == 80.0 / 100.0);
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
        // dfmt off
        return cmpLexicographically!(
            typeof(this),
            "a.alignment",
            "a.seed",
        )(this, other);
        // dfmt on
    }

    static InputRange!SeededAlignment from(AlignmentChain alignmentChain)
    {
        alias Seed = AlignmentLocationSeed;

        if (isFrontExtension(alignmentChain) && isBackExtension(alignmentChain))
        {
            // dfmt off
            return inputRangeObject(only(
                SeededAlignment(alignmentChain, Seed.front),
                SeededAlignment(alignmentChain, Seed.back),
            ));
            // dfmt on
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
            // dfmt off
            logJsonDebug(
                "info", format!"creating invalid read alignment with %d local alignments"(alignments.length),
                "alignments", alignments.toJson,
            );
            // dfmt on

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
        Returns true iff the read alignment spans a gap, ie. two alignments on
        different reference contigs with different individual extension type are
        involved.

        ---
        Case 1 (complement alignment):

                  0             ry lr   0             ry lr
            ref   |--<--<--<-+-<-+--|   |--+-<-+-<-+-<-+--|
                             | | |         | | |
            read          |--+->-+->-->-->-+->-+--|
                          0     ay         la

        Case 1 (complement alignment):

                  0             ry lr1  0 rx             lr2
            ref   |-->-->-->-+->-+--|   |--+->-+->-+->-+--|
                             | | |         | | |
            read          |--+->-+->-->-->-+->-+--|
                          0     ay         bx     la
        ---
    */
    @property bool isGap() const pure nothrow
    {
        return isParallel ^ isAntiParallel;
    }

    private @property bool joinsTwoContigs() const pure nothrow
    {
        // dfmt off
        return _length == 2 &&
            _alignments[0].contigA.id != _alignments[1].contigA.id &&
            _alignments[0].contigB.id == _alignments[1].contigB.id;
        // dfmt on
    }

    @property bool isParallel() const pure nothrow
    {
        // dfmt off
        return joinsTwoContigs &&
            _alignments[0].seed != _alignments[1].seed &&
            _alignments[0].complement == _alignments[1].complement;
        // dfmt on
    }

    @property bool isAntiParallel() const pure nothrow
    {
        // dfmt off
        return joinsTwoContigs &&
            _alignments[0].seed == _alignments[1].seed &&
            _alignments[0].complement != _alignments[1].complement;
        // dfmt on
    }

    unittest
    {
        with (AlignmentChain) with (LocalAlignment) with (Flag)
                {
                    // dfmt off
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
                                emptyFlags,
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
                                emptyFlags,
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
                                emptyFlags,
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
                                emptyFlags,
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
                                emptyFlags,
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
                                emptyFlags,
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
                    // dfmt on

                    alias getFailureMessage = (testCase, testFunction, expectedValue) => format!"expected %s.%s to be %s"(
                            testCase, testFunction, expectedValue);
                    alias toBool = (c) => c == '+' ? true : false;
                    // dfmt off
                    alias toRAT = (c) => c == 'F'
                        ? ReadAlignmentType.front
                        : c == 'B'
                            ? ReadAlignmentType.back
                            : ReadAlignmentType.gap;
                    // dfmt on

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
        // dfmt off
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
        // dfmt on
    case ReadAlignmentType.gap:
        // dfmt off
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
        // dfmt on
    case ReadAlignmentType.back:
        // dfmt off
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
        // dfmt on
    }
}

/// Generate join from read alignment.
J to(J : Join!(ReadAlignment[]))(ReadAlignment readAlignment)
{
    auto join = makeJoin!J(readAlignment);
    join.payload = [readAlignment];

    return join;
}

///
unittest
{
    with (AlignmentChain) with (LocalAlignment) with (Flag)
            {
                // dfmt off
                auto frontExtension = ReadAlignment(
                    SeededAlignment.from(AlignmentChain(
                        3,
                        Contig(1, 100),
                        Contig(1, 10),
                        emptyFlags,
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
                );
                auto backExtension = ReadAlignment(
                    SeededAlignment.from(AlignmentChain(
                        5,
                        Contig(1, 100),
                        Contig(1, 10),
                        emptyFlags,
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
                );
                auto gap = ReadAlignment(
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
                        emptyFlags,
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
                );
                // dfmt on

                auto join1 = frontExtension.to!(Join!(ReadAlignment[]));
                auto join2 = backExtension.to!(Join!(ReadAlignment[]));
                auto join3 = gap.to!(Join!(ReadAlignment[]));

                assert(join1.start == ContigNode(1, ContigPart.pre));
                assert(join1.end == ContigNode(1, ContigPart.begin));
                assert(join1.payload == [frontExtension]);

                assert(join2.start == ContigNode(1, ContigPart.end));
                assert(join2.end == ContigNode(1, ContigPart.post));
                assert(join2.payload == [backExtension]);

                assert(join3.start == ContigNode(1, ContigPart.end));
                assert(join3.end == ContigNode(2, ContigPart.end));
                assert(join3.payload == [gap]);
            }
}

/**
    A pile of read alignments belonging to the same gap/contig end.

    See_Also: Hit
*/
alias PileUp = ReadAlignment[];

PileUp[] buildPileUps(in size_t numReferenceContigs, AlignmentChain[] candidates)
{
    alias isSameRead = (a, b) => a.contigB.id == b.contigB.id;
    alias Payload = ReadAlignment[];
    alias ReadAlignmentJoin = Join!Payload;

    candidates.sort!("a.contigB.id < b.contigB.id", SwapStrategy.stable);

    // dfmt off
    auto readAlignmentJoins = candidates
        .filter!"!a.flags.disabled"
        .chunkBy!isSameRead
        .map!collectReadAlignments
        .joiner
        .filter!"a.isValid"
        .map!"a.getInOrder()"
        .map!(to!ReadAlignmentJoin);
    // dfmt on

    // dfmt off
    auto alignmentsScaffold = buildScaffold!(concatenatePayloads!Payload, Payload)(numReferenceContigs + 0, readAlignmentJoins)
        .discardAmbiguousJoins!Payload
        .mergeExtensionsWithGaps!("a ~ b", Payload);
    // dfmt on

    // dfmt off
    auto pileUps = alignmentsScaffold
        .edges
        .filter!"a.payload.length > 0"
        .map!"a.payload"
        .filter!(pileUp => pileUp.isValid)
        .array;
    // dfmt on

    return pileUps;
}

unittest
{
    // FIXME add test case with contig spanning read
    //               c1                       c2                       c3
    //        0    5   10   15         0    5   10   15         0    5   10   15
    // ref:   |---->---->----|.........|---->---->----|.........|---->---->----|
    // reads: .    .    .    :    .    :    .    .    :    .    :    .    .    .
    //        . #1 |---->---->--| .    :    .    .    :    .    :    .    .    .
    //        . #2 |----<----<--| .    :    .    .    :    .    :    .    .    .
    //        . #3 |---->---->----|    :    .    .    :    .    :    .    .    .
    //        .    . #4 |---->----|    :    .    .    :    .    :    .    .    .
    //        .    . #5 |---->---->---->----|    .    :    .    :    .    .    .
    //        .    . #6 |----<----<----<----|    .    :    .    :    .    .    .
    //        .    .    #7 |->---->---->---->--| .    :    .    :    .    .    .
    //        .    .    .    : #8 |---->---->--| .    :    .    :    .    .    .
    //        .    .    .    : #9 |----<----<--| .    :    .    :    .    .    .
    //        .    .    .    :#10 |---->---->----|    :    .    :    .    .    .
    //        .    .    .    :    #11 |----->----|    :    .    :    .    .    .
    //        .    .    .    :    .    :    .    .    :    .    :    .    .    .
    //        .    .    .    :    .    :#12 |---->---->--| .    :    .    .    .
    //        .    .    .    :    .    :#13 |----<----<--| .    :    .    .    .
    //        .    .    .    :    .    :#14 |---->---->----|    :    .    .    .
    //        .    .    .    :    .    :    .#15 |---->----|    :    .    .    .
    //        .    .    .    :    .    :    .#16 |---->---->---->----|    .    .
    //        .    .    .    :    .    :    .#17 |----<----<----<----|    .    .
    //        .    .    .    :    .    :    .   #18 |->---->---->---->--| .    .
    //        .    .    .    :    .    :    .    .    :#19 |---->---->--| .    .
    //        .    .    .    :    .    :    .    .    :#20 |----<----<--| .    .
    //        .    .    .    :    .    :    .    .    :#21 |---->---->----|    .
    //        .    .    .    :    .    :    .    .    :    #22 |----->----|    .
    import std.algorithm : clamp;

    with (AlignmentChain) with (LocalAlignment)
        {
            id_t alignmentChainId = 0;
            id_t contReadId = 0;
            ReadAlignment getDummyRead(id_t beginContigId, arithmetic_t beginIdx,
                    id_t endContigId, arithmetic_t endIdx, Complement complement)
            {
                static immutable contigLength = 16;
                static immutable gapLength = 9;
                static immutable numDiffs = 0;

                alignmentChainId += 2;
                auto readId = ++contReadId;
                // dfmt off
                auto readLength = beginContigId == endContigId
                    ? endIdx - beginIdx
                    : contigLength - beginIdx + gapLength + endIdx;
                // dfmt on
                auto flags = complement ? Flags(Flag.complement) : emptyFlags;
                coord_t firstReadBeginIdx;
                coord_t firstReadEndIdx;

                if (beginIdx < 0)
                {
                    firstReadBeginIdx = readLength - endIdx;
                    firstReadEndIdx = readLength;
                }
                else
                {
                    firstReadBeginIdx = 0;
                    firstReadEndIdx = contigLength - beginIdx;
                }

                beginIdx = clamp(beginIdx, 0, contigLength);
                endIdx = clamp(endIdx, 0, contigLength);

                if (beginContigId == endContigId)
                {
                    // dfmt off
                    return ReadAlignment(SeededAlignment.from(AlignmentChain(
                        alignmentChainId - 2,
                        Contig(beginContigId, contigLength),
                        Contig(readId, readLength),
                        flags,
                        [
                            LocalAlignment(
                                Locus(beginIdx, beginIdx + 1),
                                Locus(firstReadBeginIdx, firstReadBeginIdx + 1),
                                numDiffs,
                            ),
                            LocalAlignment(
                                Locus(endIdx - 1, endIdx),
                                Locus(firstReadEndIdx - 1, firstReadEndIdx),
                                numDiffs,
                            ),
                        ],
                    )).front);
                    // dfmt on
                }
                else
                {
                    auto secondReadBeginIdx = readLength - endIdx;
                    auto secondReadEndIdx = readLength;

                    // dfmt off
                    return ReadAlignment(
                        SeededAlignment.from(AlignmentChain(
                            alignmentChainId - 2,
                            Contig(beginContigId, contigLength),
                            Contig(readId, readLength),
                            flags,
                            [
                                LocalAlignment(
                                    Locus(beginIdx, beginIdx + 1),
                                    Locus(firstReadBeginIdx, firstReadBeginIdx + 1),
                                    numDiffs,
                                ),
                                LocalAlignment(
                                    Locus(contigLength - 1, contigLength),
                                    Locus(firstReadEndIdx - 1, firstReadEndIdx),
                                    numDiffs,
                                ),
                            ],
                        )).front,
                        SeededAlignment.from(AlignmentChain(
                            alignmentChainId - 1,
                            Contig(endContigId, contigLength),
                            Contig(readId, readLength),
                            flags,
                            [
                                LocalAlignment(
                                    Locus(0, 1),
                                    Locus(secondReadBeginIdx, secondReadBeginIdx + 1),
                                    numDiffs,
                                ),
                                LocalAlignment(
                                    Locus(endIdx - 1, endIdx),
                                    Locus(secondReadEndIdx - 1, secondReadEndIdx),
                                    numDiffs,
                                ),
                            ],
                        )).front,
                    );
                    // dfmt on
                }
            }

            id_t c1 = 1;
            id_t c2 = 2;
            id_t c3 = 3;
            // dfmt off
            auto pileUps = [
                [
                    getDummyRead(c1,  5, c1, 18, Complement.no),  //  #1
                    getDummyRead(c1,  5, c1, 18, Complement.yes), //  #2
                    getDummyRead(c1,  5, c1, 20, Complement.no),  //  #3
                    getDummyRead(c1, 10, c1, 20, Complement.no),  //  #4
                    getDummyRead(c1, 10, c2,  5, Complement.no),  //  #5
                    getDummyRead(c1, 10, c2,  5, Complement.yes), //  #6
                    getDummyRead(c1, 13, c2,  8, Complement.no),  //  #7
                    getDummyRead(c2, -5, c2,  8, Complement.no),  //  #8
                    getDummyRead(c2, -5, c2,  8, Complement.yes), //  #9
                    getDummyRead(c2, -5, c2, 10, Complement.no),  // #10
                    getDummyRead(c2, -1, c2, 10, Complement.no),  // #11
                ],
                [
                    getDummyRead(c2,  5, c2, 18, Complement.no),  // #12
                    getDummyRead(c2,  5, c2, 18, Complement.yes), // #13
                    getDummyRead(c2,  5, c2, 20, Complement.no),  // #14
                    getDummyRead(c2, 10, c2, 20, Complement.no),  // #15
                    getDummyRead(c2, 10, c3,  5, Complement.no),  // #16
                    getDummyRead(c2, 10, c3,  5, Complement.yes), // #17
                    getDummyRead(c2, 13, c3,  8, Complement.no),  // #18
                    getDummyRead(c3, -5, c3,  8, Complement.no),  // #19
                    getDummyRead(c3, -5, c3,  8, Complement.yes), // #20
                    getDummyRead(c3, -5, c3, 10, Complement.no),  // #21
                    getDummyRead(c3, -1, c3, 10, Complement.no),  // #22
                ],
            ];
            auto alignmentChains = pileUps
                .joiner
                .map!"a[]"
                .joiner
                .map!"a.alignment"
                .array;
            // dfmt on
            auto computedPileUps = buildPileUps(3, alignmentChains);

            foreach (pileUp, computedPileUp; zip(pileUps, computedPileUps))
            {
                // pileUp is subset of computedPileUp
                foreach (readAlignment; pileUp)
                {
                    assert(computedPileUp.canFind(readAlignment));
                }
                // computedPileUp is subset of pileUp
                foreach (readAlignment; computedPileUp)
                {
                    assert(pileUp.canFind(readAlignment));
                }
            }
        }
}

// Not meant for public usage.
ReadAlignment[] collectReadAlignments(Chunk)(Chunk sameReadAlignments)
{
    // dfmt off
    alias beginRelToContigB = (alignment) => alignment.complement
        ? alignment.contigB.length - alignment.last.contigB.end
        : alignment.first.contigB.begin;
    alias endRelToContigB = (alignment) => alignment.complement
        ? alignment.contigB.length - alignment.first.contigB.begin
        : alignment.last.contigB.end;
    alias seedRelToContigB = (alignment) => alignment.complement
        ? 0 - alignment.seed
        : 0 + alignment.seed;
    alias orderByLocusAndSeed = orderLexicographically!(
        SeededAlignment,
        beginRelToContigB,
        endRelToContigB,
        seedRelToContigB,
    );
    // dfmt on
    alias seededPartsOfOneAlignment = (a, b) => a.alignment == b.alignment && a.seed != b.seed;
    alias shareReadSequence = (a, b) => endRelToContigB(a) > beginRelToContigB(b);
    alias emptyRange = () => cast(ReadAlignment[])[];

    auto seededAlignments = sameReadAlignments.map!(SeededAlignment.from).joiner.array;
    seededAlignments.sort!orderByLocusAndSeed;

    if (seededAlignments.length == 0)
    {
        return emptyRange();
    }

    // Validate seeded alignments and discard accordingly.
    foreach (saPair; seededAlignments.slide!(No.withPartial)(2))
    {
        if (shareReadSequence(saPair[0], saPair[1])
                && !seededPartsOfOneAlignment(saPair[0], saPair[1]))
        {
            return emptyRange();
        }
    }

    // Collect read alignments
    bool startWithExtension = beginRelToContigB(seededAlignments[0]) > 0;
    size_t sliceStart = startWithExtension ? 1 : 0;

    // dfmt off
    auto remainingReadAlignments = iota(sliceStart, seededAlignments.length, 2)
        .map!(i => ReadAlignment(seededAlignments[i .. min(i + 2, $)]));
    auto readAlignments = startWithExtension
        ? chain(
            only(ReadAlignment(seededAlignments[0 .. 1])),
            remainingReadAlignments,
        ).array
        : remainingReadAlignments.array;
    // dfmt on

    if (readAlignments.any!"!a.isValid")
    {
        // If one read alignment is invalid we should not touch this read at all.
        return emptyRange();
    }

    return readAlignments;
}

unittest
{
    alias Contig = AlignmentChain.Contig;
    alias Flags = AlignmentChain.Flags;
    immutable emptyFlags = AlignmentChain.emptyFlags;
    immutable complement = AlignmentChain.Flag.complement;
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias Locus = LocalAlignment.Locus;
    alias Seed = AlignmentLocationSeed;

    {
        // Case 1:
        //
        // |-->-->-->--|  |-->-->-->--| |-->-->-->--|
        //
        //       |----->  <-----------> <-----|
        // dfmt off
        auto alignmentChains = [
            AlignmentChain(
                0,
                Contig(1, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(10, 20),
                    Locus(0, 10),
                )]
            ),
            AlignmentChain(
                1,
                Contig(2, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(0, 20),
                    Locus(20, 40),
                )]
            ),
            AlignmentChain(
                2,
                Contig(3, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(0, 10),
                    Locus(50, 60),
                )]
            ),
        ];
        assert(collectReadAlignments(alignmentChains).equal([
            ReadAlignment([
                SeededAlignment(alignmentChains[0], Seed.back),
                SeededAlignment(alignmentChains[1], Seed.front),
            ]),
            ReadAlignment([
                SeededAlignment(alignmentChains[1], Seed.back),
                SeededAlignment(alignmentChains[2], Seed.front),
            ]),
        ]));
        // dfmt on
    }
    {
        // Case 2:
        //
        // |-->-->-->--|  |--<--<--<--| |-->-->-->--|
        //
        //       |----->  <-----------> <-----|
        // dfmt off
        auto alignmentChains = [
            AlignmentChain(
                0,
                Contig(1, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(10, 20),
                    Locus(0, 10),
                )]
            ),
            AlignmentChain(
                1,
                Contig(2, 20),
                Contig(1, 60),
                Flags(complement),
                [LocalAlignment(
                    Locus(0, 20),
                    Locus(20, 40),
                )]
            ),
            AlignmentChain(
                2,
                Contig(3, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(0, 10),
                    Locus(50, 60),
                )]
            ),
        ];
        assert(collectReadAlignments(alignmentChains).equal([
            ReadAlignment([
                SeededAlignment(alignmentChains[0], Seed.back),
                SeededAlignment(alignmentChains[1], Seed.back),
            ]),
            ReadAlignment([
                SeededAlignment(alignmentChains[1], Seed.front),
                SeededAlignment(alignmentChains[2], Seed.front),
            ]),
        ]));
        // dfmt on
    }
    {
        // Case 3:
        //
        //                |--<--<--<--| |-->-->-->--|
        //
        //                <-----------> <-----|
        // dfmt off
        auto alignmentChains = [
            AlignmentChain(
                1,
                Contig(2, 20),
                Contig(1, 60),
                Flags(complement),
                [LocalAlignment(
                    Locus(0, 20),
                    Locus(20, 40),
                )]
            ),
            AlignmentChain(
                2,
                Contig(3, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(0, 10),
                    Locus(50, 60),
                )]
            ),
        ];
        assert(collectReadAlignments(alignmentChains).equal([
            ReadAlignment([
                SeededAlignment(alignmentChains[0], Seed.back),
            ]),
            ReadAlignment([
                SeededAlignment(alignmentChains[0], Seed.front),
                SeededAlignment(alignmentChains[1], Seed.front),
            ]),
        ]));
        // dfmt on
    }
    {
        // Case 4:
        //
        // |-->-->-->--|  |--<--<--<--|
        //
        //       |----->  <----------->
        // dfmt off
        auto alignmentChains = [
            AlignmentChain(
                0,
                Contig(1, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(10, 20),
                    Locus(0, 10),
                )]
            ),
            AlignmentChain(
                1,
                Contig(2, 20),
                Contig(1, 60),
                Flags(complement),
                [LocalAlignment(
                    Locus(0, 20),
                    Locus(20, 40),
                )]
            ),
        ];
        assert(collectReadAlignments(alignmentChains).equal([
            ReadAlignment([
                SeededAlignment(alignmentChains[0], Seed.back),
                SeededAlignment(alignmentChains[1], Seed.back),
            ]),
            ReadAlignment([
                SeededAlignment(alignmentChains[1], Seed.front),
            ]),
        ]));
        // dfmt on
    }
    {
        // Case 5:
        //
        //                |--<--<--<--|
        //
        //                <----------->
        // dfmt off
        auto alignmentChains = [
            AlignmentChain(
                1,
                Contig(2, 20),
                Contig(1, 60),
                Flags(complement),
                [LocalAlignment(
                    Locus(0, 20),
                    Locus(20, 40),
                )]
            ),
        ];
        assert(collectReadAlignments(alignmentChains).equal([
            ReadAlignment([
                SeededAlignment(alignmentChains[0], Seed.back),
            ]),
            ReadAlignment([
                SeededAlignment(alignmentChains[0], Seed.front),
            ]),
        ]));
        // dfmt on
    }
    {
        // Case 6:
        //
        // |-->-->-->--|
        // |-->-->-->--|  |--<--<--<--| |-->-->-->--|
        //
        //       |----->  <-----------> <-----|
        // dfmt off
        auto alignmentChains = [
            AlignmentChain(
                0,
                Contig(1, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(10, 20),
                    Locus(0, 10),
                )]
            ),
            AlignmentChain(
                1,
                Contig(2, 20),
                Contig(1, 60),
                Flags(complement),
                [LocalAlignment(
                    Locus(0, 20),
                    Locus(20, 40),
                )]
            ),
            AlignmentChain(
                2,
                Contig(3, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(0, 10),
                    Locus(50, 60),
                )]
            ),
            AlignmentChain(
                3,
                Contig(4, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(10, 20),
                    Locus(0, 10),
                )]
            ),
        ];
        // dfmt on
        assert(collectReadAlignments(alignmentChains).length == 0);
    }
    {
        // Case 7:
        //
        //                |-->-->-->--|
        // |-->-->-->--|  |--<--<--<--| |-->-->-->--|
        //
        //       |----->  <-----------> <-----|
        // dfmt off
        auto alignmentChains = [
            AlignmentChain(
                0,
                Contig(1, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(10, 20),
                    Locus(0, 10),
                )]
            ),
            AlignmentChain(
                1,
                Contig(2, 20),
                Contig(1, 60),
                Flags(complement),
                [LocalAlignment(
                    Locus(0, 20),
                    Locus(20, 40),
                )]
            ),
            AlignmentChain(
                2,
                Contig(3, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(0, 10),
                    Locus(50, 60),
                )]
            ),
            AlignmentChain(
                3,
                Contig(4, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(0, 20),
                    Locus(20, 40),
                )]
            ),
        ];
        // dfmt on
        assert(collectReadAlignments(alignmentChains).length == 0);
    }
    {
        // Case 8:
        //
        //                        |-->-->-->--|
        // |-->-->-->--|  |--<--<--<--| |-->-->-->--|
        //
        //       |----->  <-----------> <-----|
        // dfmt off
        auto alignmentChains = [
            AlignmentChain(
                0,
                Contig(1, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(10, 20),
                    Locus(0, 10),
                )]
            ),
            AlignmentChain(
                1,
                Contig(2, 20),
                Contig(1, 60),
                Flags(complement),
                [LocalAlignment(
                    Locus(0, 20),
                    Locus(20, 40),
                )]
            ),
            AlignmentChain(
                2,
                Contig(3, 20),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(0, 10),
                    Locus(50, 60),
                )]
            ),
            AlignmentChain(
                3,
                Contig(4, 30),
                Contig(1, 60),
                emptyFlags,
                [LocalAlignment(
                    Locus(0, 30),
                    Locus(30, 60),
                )]
            ),
        ];
        // dfmt on
        assert(collectReadAlignments(alignmentChains).length == 0);
    }
}

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
        // dfmt off
        return pileUp[0].isFrontExtension
            ? ReadAlignmentType.front
            : ReadAlignmentType.back;
        // dfmt on
    }
}

bool isValid(in PileUp pileUp) pure nothrow
{
    return pileUp.isExtension ^ pileUp.isGap;
}

bool isExtension(in PileUp pileUp) pure nothrow
{
    if (pileUp[0].isFrontExtension)
    {
        return pileUp.all!(readAlignment => readAlignment.isFrontExtension);
    }
    else if (pileUp[0].isBackExtension)
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
    // dfmt off
    return pileUp.isGap && pileUp
        .filter!(readAlignment => readAlignment.isGap)
        .front
        .isParallel;
    // dfmt on
}

auto isAntiParallel(in PileUp pileUp) pure nothrow
{
    // dfmt off
    return pileUp.isGap && pileUp
        .filter!(readAlignment => readAlignment.isGap)
        .front
        .isAntiParallel;
    // dfmt on
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

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
import std.array : appender, array, minimallyInitializedArray;
import std.conv : to;
import std.exception : assertNotThrown, assertThrown, enforce, ErrnoException;
import std.format : format;
import std.math : sgn;
import std.range : assumeSorted, chain, chunks, InputRange, inputRangeObject,
    iota, only, slide, takeNone, zip;
import std.string : capitalize, split;
import std.stdio : File, LockType;
import std.typecons : BitFlags, PhobosFlag = Flag, No, tuple, Tuple, Yes;
import std.traits : isArray, TemplateArgsOf, TemplateOf;
import vibe.data.json : toJson = serializeToJson;

debug import std.stdio : writefln, writeln;

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

class PileUpDbException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

struct PileUpDb
{
    private alias LocalAlignment = AlignmentChain.LocalAlignment;
    private alias TracePoint = LocalAlignment.TracePoint;
    private alias DbSlices = Tuple!(
        ArrayStorage!(StorageType!PileUp), "pileUps",
        ArrayStorage!(StorageType!ReadAlignment), "readAlignments",
        ArrayStorage!(StorageType!SeededAlignment), "seededAlignments",
        ArrayStorage!(StorageType!LocalAlignment), "localAlignments",
        ArrayStorage!(StorageType!TracePoint), "tracePoints",
    );

    private File pileUpDb;
    private PileUpDbIndex dbIndex;
    DbSlices dbSlices;

    static PileUpDb parse(in string dbFile)
    {
        auto pileUpDb = File(dbFile, "rb");
        pileUpDb.lock(LockType.read);

        return PileUpDb(pileUpDb);
    }

    PileUp[] opIndex()
    {
        ensureDbIndex();

        return readSlice(0, length);
    }

    PileUp opIndex(size_t i)
    {
        ensureDbIndex();
        enforce!PileUpDbException(
            i < length,
            format!"cannot read block %d in `%s`: out of bounds [0, %d)"(
                    i, pileUpDb.name, length)
        );

        return readSlice(i, i + 1)[0];
    }

    PileUp[] opIndex(size_t[2] slice)
    {
        auto from = slice[0];
        auto to = slice[1];
        ensureDbIndex();
        enforce!PileUpDbException(
            to <= length,
            format!"cannot read blocks %d-%d in `%s`: out of bounds [0, %d]"(
                    from, to, pileUpDb.name, length)
        );

        return readSlice(from, to);
    }

    size_t[2] opSlice(size_t dim)(size_t from, size_t to)
            if (dim == 0)
    {
        assert(from < to, "invalid slice");

        return [from, to];
    }

    @property size_t length()
    {
        ensureDbIndex();

        return (dbIndex.endPtr!PileUp - dbIndex.beginPtr!PileUp) / StorageType!PileUp.sizeof;
    }

    alias opDollar = length;

    private void ensureDbIndex()
    {
        if (dbIndex != dbIndex.init)
            return;

        dbIndex = readRecord!PileUpDbIndex();
    }

    private PileUp[] readSlice(size_t from, size_t to)
    {
        assert(from <= to && to <= length);

        // Step 1: determine memory requirements and DB slices
        dbSlices = getDbSlices(from, to);

        // Step 2: allocate minimally initialized memory for all blocks
        auto pileUps = minimallyInitializedArray!(PileUp[])(dbSlices.pileUps.length);
        auto readAlignments = minimallyInitializedArray!(ReadAlignment[])(dbSlices.readAlignments.length);
        auto seededAlignments = minimallyInitializedArray!(SeededAlignment[])(dbSlices.seededAlignments.length);
        auto localAlignments = minimallyInitializedArray!(LocalAlignment[])(dbSlices.localAlignments.length);
        auto tracePoints = minimallyInitializedArray!(TracePoint[])(dbSlices.tracePoints.length);

        // Step 3: parse each record for each block assigning already
        //         allocated array slices to the array fields
        parse(pileUps, readAlignments);
        parse(seededAlignments, localAlignments);
        parse(readAlignments, seededAlignments);
        parse(localAlignments, tracePoints);
        parse(tracePoints);

        return pileUps;
    }

    private DbSlices getDbSlices(size_t from, size_t to)
    {
        auto pileUps = dbIndex.pileUps[from .. to];
        auto firstPileUp = readRecordAt!(StorageType!PileUp)(pileUps[0]);
        auto lastPileUp = readRecordAt!(StorageType!PileUp)(pileUps[$ - 1]);

        auto readAlignments = ArrayStorage!(StorageType!ReadAlignment).fromPtrs(
            firstPileUp[0],
            lastPileUp[$],
        );
        auto firstReadAlignment = readRecordAt!(StorageType!ReadAlignment)(readAlignments[0]);
        auto lastReadAlignment = readRecordAt!(StorageType!ReadAlignment)(readAlignments[$ - 1]);

        auto seededAlignments = ArrayStorage!(StorageType!SeededAlignment).fromPtrs(
            firstReadAlignment[0],
            lastReadAlignment[$],
        );
        auto firstSeededAlignment = readRecordAt!(StorageType!SeededAlignment)(seededAlignments[0]);
        auto lastSeededAlignment = readRecordAt!(StorageType!SeededAlignment)(seededAlignments[$ - 1]);

        auto localAlignments = ArrayStorage!(StorageType!LocalAlignment).fromPtrs(
            firstSeededAlignment.localAlignments[0],
            lastSeededAlignment.localAlignments[$],
        );
        auto firstLocalAlignment = readRecordAt!(StorageType!LocalAlignment)(localAlignments[0]);
        auto lastLocalAlignment = readRecordAt!(StorageType!LocalAlignment)(localAlignments[$ - 1]);

        auto tracePoints = ArrayStorage!(StorageType!TracePoint).fromPtrs(
            firstLocalAlignment.tracePoints[0],
            lastLocalAlignment.tracePoints[$],
        );

        return DbSlices(
            pileUps,
            readAlignments,
            seededAlignments,
            localAlignments,
            tracePoints,
        );
    }

    private void parse(ref PileUp[] pileUps, ReadAlignment[] readAlignments)
    {
        static assert(PileUp.sizeof == StorageType!PileUp.sizeof);
        pileUpDb.seek(dbSlices.pileUps.ptr);
        pileUps = readRecords(pileUps);

        size_t[2] pileUpSlice;
        foreach (ref pileUp; pileUps)
        {
            auto pileUpStorage = *cast(StorageType!PileUp*) &pileUp;
            pileUpSlice[0] = pileUpSlice[1];
            pileUpSlice[1] += pileUpStorage.length;

            pileUp = readAlignments[pileUpSlice[0] .. pileUpSlice[1]];
        }
    }

    private void parse(
        ref SeededAlignment[] seededAlignments,
        AlignmentChain.LocalAlignment[] localAlignments
    )
    {
        // Parse `SeededAlignment`s
        alias Contig = AlignmentChain.Contig;
        pileUpDb.seek(dbSlices.seededAlignments.ptr);

        size_t[2] seededAlignmentsSlice;
        foreach (ref seededAlignment; seededAlignments)
        {
            auto seededAlignmentStorage = readRecord!(StorageType!SeededAlignment);

            seededAlignmentsSlice[0] = seededAlignmentsSlice[1];
            seededAlignmentsSlice[1] += seededAlignmentStorage.localAlignments.length;

            seededAlignment = SeededAlignment(
                AlignmentChain(
                    seededAlignmentStorage.id,
                    Contig(
                        seededAlignmentStorage.contigAId,
                        seededAlignmentStorage.contigALength,
                    ),
                    Contig(
                        seededAlignmentStorage.contigBId,
                        seededAlignmentStorage.contigBLength,
                    ),
                    seededAlignmentStorage.flags,
                    localAlignments[seededAlignmentsSlice[0] .. seededAlignmentsSlice[1]],
                    seededAlignmentStorage.tracePointDistance,
                ),
                seededAlignmentStorage.seed,
            );
        }
    }

    private void parse(ref ReadAlignment[] readAlignments, SeededAlignment[] seededAlignments)
    {
        pileUpDb.seek(dbSlices.readAlignments.ptr);

        size_t[2] readAlignmentsSlice;
        foreach (ref readAlignment; readAlignments)
        {
            auto readAlignmentStorage = readRecord!(StorageType!ReadAlignment);

            readAlignmentsSlice[0] = readAlignmentsSlice[1];
            readAlignmentsSlice[1] += readAlignmentStorage.length;

            readAlignment = ReadAlignment(
                seededAlignments[readAlignmentsSlice[0] .. readAlignmentsSlice[1]]
            );
        }
    }

    private void parse(
        ref AlignmentChain.LocalAlignment[] localAlignments,
        AlignmentChain.LocalAlignment.TracePoint[] tracePoints
    )
    {
        alias LocalAlignment = AlignmentChain.LocalAlignment;
        alias Locus = LocalAlignment.Locus;
        pileUpDb.seek(dbSlices.localAlignments.ptr);

        size_t[2] localAlignmentsSlice;
        foreach (ref localAlignment; localAlignments)
        {
            auto localAlignmentStorage = readRecord!(StorageType!LocalAlignment);

            localAlignmentsSlice[0] = localAlignmentsSlice[1];
            localAlignmentsSlice[1] += localAlignmentStorage.tracePoints.length;

            localAlignment = LocalAlignment(
                Locus(
                    localAlignmentStorage.contigABegin,
                    localAlignmentStorage.contigAEnd,
                ),
                Locus(
                    localAlignmentStorage.contigBBegin,
                    localAlignmentStorage.contigBEnd,
                ),
                localAlignmentStorage.numDiffs,
                tracePoints[localAlignmentsSlice[0] .. localAlignmentsSlice[1]],
            );
        }
    }

    private void parse(ref AlignmentChain.LocalAlignment.TracePoint[] tracePoints)
    {
        alias LocalAlignment = AlignmentChain.LocalAlignment;
        alias TracePoint = LocalAlignment.TracePoint;

        static assert(TracePoint.sizeof == StorageType!TracePoint.sizeof);
        pileUpDb.seek(dbSlices.tracePoints.ptr);
        tracePoints = readRecords(tracePoints);
    }

    private T readRecordAt(T)(size_t ptr)
    {
        pileUpDb.seek(ptr.to!long);

        return readRecord!T();
    }

    private T readRecord(T)()
    {
        return readRecords(new T[1])[0];
    }

    private T readRecords(T)(T records) if (isArray!T)
    {
        auto expectedLength = records.length;
        try
        {
            records = pileUpDb.rawRead(records);

            enforce!PileUpDbException(
                records.length == expectedLength,
                format!"malformed pile up DB `%s`: premature end of file"(
                        pileUpDb.name)
            );
        }
        catch (ErrnoException e)
        {
            throw new PileUpDbException(
                format!"malformed pile up DB `%s`: cannot read record of type %s: %s"(
                        pileUpDb.name, T.stringof, e.msg),
            );
        }

        return records;
    }
}

void writePileUpsDb(in PileUp[] pileUps, in string dbFile)
{
    auto pileUpDb = File(dbFile, "wb");
    pileUpDb.lock();
    scope (exit)
        pileUpDb.unlock();

    writePileUpsDb(pileUps, pileUpDb);
}

void writePileUpsDb(in PileUp[] pileUps, File pileUpDb)
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

    PileUpDbIndex dbIndex = buildPileUpDbIndex(pileUps);

    pileUpDb.rawWrite([dbIndex]);

    writePileUpsDbBlock!PileUp(pileUpDb, pileUps, dbIndex);
    writePileUpsDbBlock!ReadAlignment(pileUpDb, pileUps, dbIndex);
    writePileUpsDbBlock!SeededAlignment(pileUpDb, pileUps, dbIndex);
    writePileUpsDbBlock!LocalAlignment(pileUpDb, pileUps, dbIndex);
    writePileUpsDbBlock!TracePoint(pileUpDb, pileUps, dbIndex);
}

unittest
{
    import dentist.util.tempfile : mkstemp;
    import std.file : remove;

    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

    auto pileUps = getPileUpsTestData();

    enum numPileUps = 2;
    enum numReadAlignments = 5;
    enum numSeededAlignments = 7;
    enum numLocalAlignments = 8;
    enum numTracePoints = 393;
    enum totalDbSize =
        PileUpDbIndex.sizeof +
        StorageType!PileUp.sizeof * numPileUps +
        StorageType!ReadAlignment.sizeof * numReadAlignments +
        StorageType!SeededAlignment.sizeof * numSeededAlignments +
        StorageType!LocalAlignment.sizeof * numLocalAlignments +
        StorageType!TracePoint.sizeof * numTracePoints;

    auto tmpDb = mkstemp("./.unittest-XXXXXX");
    scope (exit)
    {
        tmpDb.file.close();
        remove(tmpDb.name);
    }

    writePileUpsDb(pileUps, tmpDb.file);
    tmpDb.file.sync();

    assert(tmpDb.file.size == totalDbSize);

    tmpDb.file.rewind();
    auto pileUpDb = PileUpDb(tmpDb.file);

    assert(pileUpDb[] == pileUps);
}

private PileUpDbIndex buildPileUpDbIndex(in PileUp[] pileUps) nothrow pure
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = AlignmentChain.LocalAlignment.TracePoint;

    PileUpDbIndex dbIndex;

    dbIndex.beginPtr!PileUp = PileUpDbIndex.sizeof;
    dbIndex.endPtr!PileUp = StorageType!PileUp.sizeof * pileUps.length;
    foreach (ref pileUp; pileUps)
    {
        dbIndex.endPtr!ReadAlignment += StorageType!ReadAlignment.sizeof * pileUp.length;
        foreach (ref readAlignment; pileUp)
        {
            dbIndex.endPtr!SeededAlignment += StorageType!SeededAlignment.sizeof * readAlignment.length;
            foreach (ref seededAlignment; readAlignment[])
            {
                dbIndex.endPtr!LocalAlignment += StorageType!LocalAlignment.sizeof * seededAlignment.localAlignments.length;
                foreach (ref localAlignment; seededAlignment.localAlignments)
                {
                    dbIndex.endPtr!TracePoint += StorageType!TracePoint.sizeof * localAlignment.tracePoints.length;
                }
            }
        }
    }

    dbIndex.readAlignmentsPtr += dbIndex.pileUpsPtr;
    dbIndex.seededAlignmentsPtr += dbIndex.readAlignmentsPtr;
    dbIndex.localAlignmentsPtr += dbIndex.seededAlignmentsPtr;
    dbIndex.tracePointsPtr += dbIndex.localAlignmentsPtr;
    dbIndex.eofPtr += dbIndex.tracePointsPtr;

    return dbIndex;
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = AlignmentChain.LocalAlignment.TracePoint;

    auto pileUps = getPileUpsTestData();
    auto dbIndex = buildPileUpDbIndex(pileUps);

    enum numPileUps = 2;
    enum numReadAlignments = 5;
    enum numSeededAlignments = 7;
    enum numLocalAlignments = 8;
    enum numTracePoints = 393;

    assert(dbIndex.pileUpsPtr == PileUpDbIndex.sizeof);
    assert(dbIndex.readAlignmentsPtr ==
            dbIndex.pileUpsPtr +
            StorageType!PileUp.sizeof * numPileUps);
    assert(dbIndex.seededAlignmentsPtr ==
            dbIndex.readAlignmentsPtr +
            StorageType!ReadAlignment.sizeof * numReadAlignments);
    assert(dbIndex.localAlignmentsPtr ==
            dbIndex.seededAlignmentsPtr +
            StorageType!SeededAlignment.sizeof * numSeededAlignments);
    assert(dbIndex.tracePointsPtr ==
            dbIndex.localAlignmentsPtr +
            StorageType!LocalAlignment.sizeof * numLocalAlignments);
    assert(dbIndex.eofPtr ==
            dbIndex.tracePointsPtr +
            StorageType!TracePoint.sizeof * numTracePoints);
}

void writePileUpsDbBlock(T)(ref File pileUpDb, in PileUp[] pileUps, in PileUpDbIndex dbIndex)
    if (is(T == PileUp))
{
    auto readAlignments = ArrayStorage!(StorageType!ReadAlignment)(dbIndex.beginPtr!ReadAlignment);

    version (assert)
    {
        auto storedPileUps = ArrayStorage!(StorageType!PileUp)(dbIndex.beginPtr!PileUp);
        assert(storedPileUps.ptr == pileUpDb.tell());
    }
    foreach (ref pileUp; pileUps)
    {
        readAlignments.length = pileUp.length;

        pileUpDb.rawWrite([readAlignments]);

        readAlignments.ptr = readAlignments[$];

        version (assert)
        {
            ++storedPileUps.length;
            assert(storedPileUps[$] == pileUpDb.tell());
        }
    }
}

void writePileUpsDbBlock(T)(ref File pileUpDb, in PileUp[] pileUps, in PileUpDbIndex dbIndex)
    if (is(T == ReadAlignment))
{
    auto seededAlignments = ArrayStorage!(StorageType!SeededAlignment)(dbIndex.beginPtr!SeededAlignment);

    version (assert)
    {
        auto readAlignments = ArrayStorage!(StorageType!ReadAlignment)(dbIndex.beginPtr!ReadAlignment);
        assert(readAlignments.ptr == pileUpDb.tell());
    }
    foreach (ref pileUp; pileUps)
    {
        foreach (ref readAlignment; pileUp)
        {
            seededAlignments.length = readAlignment.length;

            pileUpDb.rawWrite([seededAlignments]);

            seededAlignments.ptr = seededAlignments[$];

            version (assert)
            {
                ++readAlignments.length;
                assert(readAlignments[$] == pileUpDb.tell());
            }
        }
    }
}

void writePileUpsDbBlock(T)(ref File pileUpDb, in PileUp[] pileUps, in PileUpDbIndex dbIndex)
    if (is(T == SeededAlignment))
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;

    auto localAlignments = ArrayStorage!(StorageType!LocalAlignment)(dbIndex.beginPtr!LocalAlignment);

    version (assert)
    {
        auto seededAlignments = ArrayStorage!(StorageType!SeededAlignment)(dbIndex.beginPtr!SeededAlignment);
        assert(seededAlignments.ptr == pileUpDb.tell());
    }
    foreach (ref pileUp; pileUps)
    {
        foreach (ref readAlignment; pileUp)
        {
            foreach (ref seededAlignment; readAlignment[])
            {
                localAlignments.length = seededAlignment.localAlignments.length;

                pileUpDb.rawWrite([StorageType!SeededAlignment(
                    seededAlignment.id,
                    seededAlignment.contigA.id,
                    seededAlignment.contigA.length,
                    seededAlignment.contigB.id,
                    seededAlignment.contigB.length,
                    seededAlignment.flags,
                    localAlignments,
                    seededAlignment.tracePointDistance,
                    seededAlignment.seed,
                )]);

                localAlignments.ptr = localAlignments[$];

                version (assert)
                {
                    ++seededAlignments.length;
                  assert(seededAlignments[$] == pileUpDb.tell());
              }
            }
        }
    }
}

void writePileUpsDbBlock(T)(ref File pileUpDb, in PileUp[] pileUps, in PileUpDbIndex dbIndex)
    if (is(T == AlignmentChain.LocalAlignment))
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

    auto tracePoints = ArrayStorage!(StorageType!TracePoint)(dbIndex.beginPtr!TracePoint);

    version (assert)
    {
        auto localAlignments = ArrayStorage!(StorageType!LocalAlignment)(dbIndex.beginPtr!LocalAlignment);
        assert(localAlignments.ptr == pileUpDb.tell());
    }
    foreach (ref pileUp; pileUps)
    {
        foreach (ref readAlignment; pileUp)
        {
            foreach (ref seededAlignment; readAlignment[])
            {
                foreach (ref localAlignment; seededAlignment.localAlignments)
                {
                    tracePoints.length = localAlignment.tracePoints.length;

                    pileUpDb.rawWrite([StorageType!LocalAlignment(
                        localAlignment.contigA.begin,
                        localAlignment.contigA.end,
                        localAlignment.contigB.begin,
                        localAlignment.contigB.end,
                        localAlignment.numDiffs,
                        tracePoints,
                    )]);

                    tracePoints.ptr = tracePoints[$];

                    version (assert)
                    {
                        ++localAlignments.length;
                        assert(localAlignments[$] == pileUpDb.tell());
                    }
                }
            }
        }
    }
}

void writePileUpsDbBlock(T)(ref File pileUpDb, in PileUp[] pileUps, in PileUpDbIndex dbIndex)
    if (is(T == AlignmentChain.LocalAlignment.TracePoint))
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

    version (assert)
    {
        auto tracePoints = ArrayStorage!(StorageType!TracePoint)(dbIndex.beginPtr!TracePoint);
        assert(tracePoints.ptr == pileUpDb.tell());
    }
    foreach (ref pileUp; pileUps)
    {
        foreach (ref readAlignment; pileUp)
        {
            foreach (ref seededAlignment; readAlignment[])
            {
                foreach (ref localAlignment; seededAlignment.localAlignments)
                {
                    pileUpDb.rawWrite(cast(const(StorageType!TracePoint[])) localAlignment.tracePoints);

                    version (assert)
                    {
                        tracePoints.length += localAlignment.tracePoints.length;
                        assert(tracePoints[$] == pileUpDb.tell());
                    }
                }
            }
        }
    }
}

private struct PileUpDbIndex
{
    static struct EOF;

    private static template NextType(T)
    {
        static if (is(T == PileUp))
            alias NextType = ReadAlignment;
        else static if (is(T == ReadAlignment))
            alias NextType = SeededAlignment;
        else static if (is(T == SeededAlignment))
            alias NextType = AlignmentChain.LocalAlignment;
        else static if (is(T == AlignmentChain.LocalAlignment))
            alias NextType = AlignmentChain.LocalAlignment.TracePoint;
        else static if (is(T == AlignmentChain.LocalAlignment.TracePoint))
            alias NextType = EOF;
    }

    private static template fieldPtr(T)
    {
        static if (is(T == PileUp))
            alias fieldPtr = pileUpsPtr;
        else static if (is(T == ReadAlignment))
            alias fieldPtr = readAlignmentsPtr;
        else static if (is(T == SeededAlignment))
            alias fieldPtr = seededAlignmentsPtr;
        else static if (is(T == AlignmentChain.LocalAlignment))
            alias fieldPtr = localAlignmentsPtr;
        else static if (is(T == AlignmentChain.LocalAlignment.TracePoint))
            alias fieldPtr = tracePointsPtr;
        else static if (is(T == EOF))
            alias fieldPtr = eofPtr;
    }

    size_t pileUpsPtr;
    size_t readAlignmentsPtr;
    size_t seededAlignmentsPtr;
    size_t localAlignmentsPtr;
    size_t tracePointsPtr;
    size_t eofPtr;

    @property auto ref size_t beginPtr(T)() pure nothrow
    {
        return fieldPtr!T;
    }

    @property size_t beginPtr(T)() const pure nothrow
    {
        return fieldPtr!T;
    }

    @property auto ref size_t endPtr(T)() pure nothrow
    {
        return fieldPtr!(NextType!T);
    }

    @property size_t endPtr(T)() const pure nothrow
    {
        return fieldPtr!(NextType!T);
    }

    @property ArrayStorage!(StorageType!T) arrayStorage(T)() pure nothrow
    {
        return typeof(return).fromPtrs(beginPtr!T, endPtr!T);
    }

    @property alias pileUps = arrayStorage!PileUp;
    @property alias readAlignments = arrayStorage!PileUp;
    @property alias seededAlignments = arrayStorage!PileUp;
    @property alias localAlignments = arrayStorage!PileUp;
    @property alias tracePoints = arrayStorage!PileUp;
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = AlignmentChain.LocalAlignment.TracePoint;
    enum begin = 1;
    enum end = 2;
    enum modified = 3;


    {
        PileUpDbIndex dbIndex;

        dbIndex.pileUpsPtr = begin;
        dbIndex.readAlignmentsPtr = end;

        assert(dbIndex.beginPtr!PileUp == begin);
        assert(dbIndex.endPtr!PileUp == end);

        dbIndex.beginPtr!PileUp = modified;
        dbIndex.endPtr!PileUp = modified;

        assert(dbIndex.pileUpsPtr == modified);
        assert(dbIndex.readAlignmentsPtr == modified);
    }
    {
        PileUpDbIndex dbIndex;

        dbIndex.readAlignmentsPtr = begin;
        dbIndex.seededAlignmentsPtr = end;

        assert(dbIndex.beginPtr!ReadAlignment == begin);
        assert(dbIndex.endPtr!ReadAlignment == end);

        dbIndex.beginPtr!ReadAlignment = modified;
        dbIndex.endPtr!ReadAlignment = modified;

        assert(dbIndex.readAlignmentsPtr == modified);
        assert(dbIndex.seededAlignmentsPtr == modified);
    }
    {
        PileUpDbIndex dbIndex;

        dbIndex.seededAlignmentsPtr = begin;
        dbIndex.localAlignmentsPtr = end;

        assert(dbIndex.beginPtr!SeededAlignment == begin);
        assert(dbIndex.endPtr!SeededAlignment == end);

        dbIndex.beginPtr!SeededAlignment = modified;
        dbIndex.endPtr!SeededAlignment = modified;

        assert(dbIndex.seededAlignmentsPtr == modified);
        assert(dbIndex.localAlignmentsPtr == modified);
    }
    {
        PileUpDbIndex dbIndex;

        dbIndex.localAlignmentsPtr = begin;
        dbIndex.tracePointsPtr = end;

        assert(dbIndex.beginPtr!LocalAlignment == begin);
        assert(dbIndex.endPtr!LocalAlignment == end);

        dbIndex.beginPtr!LocalAlignment = modified;
        dbIndex.endPtr!LocalAlignment = modified;

        assert(dbIndex.localAlignmentsPtr == modified);
        assert(dbIndex.tracePointsPtr == modified);
    }
    {
        PileUpDbIndex dbIndex;

        dbIndex.tracePointsPtr = begin;
        dbIndex.eofPtr = end;

        assert(dbIndex.beginPtr!TracePoint == begin);
        assert(dbIndex.endPtr!TracePoint == end);

        dbIndex.beginPtr!TracePoint = modified;
        dbIndex.endPtr!TracePoint = modified;

        assert(dbIndex.tracePointsPtr == modified);
        assert(dbIndex.eofPtr == modified);
    }
}

private template StorageType(T)
{
    static if (is(T == PileUp))
        alias StorageType = ArrayStorage!(StorageType!ReadAlignment);
    else static if (is(T == ReadAlignment))
        alias StorageType = ArrayStorage!(StorageType!SeededAlignment);
    else static if (is(T == SeededAlignment))
        alias StorageType = SeededAlignmentStorage;
    else static if (is(T == AlignmentChain.LocalAlignment[]))
        alias StorageType = ArrayStorage!(StorageType!(AlignmentChain.LocalAlignment));
    else static if (is(T == AlignmentChain.LocalAlignment))
        alias StorageType = LocalAlignmentStorage;
    else static if (is(T == AlignmentChain.LocalAlignment.TracePoint[]))
        alias StorageType = ArrayStorage!(StorageType!(AlignmentChain.LocalAlignment.TracePoint));
    else static if (is(T == AlignmentChain.LocalAlignment.TracePoint))
        alias StorageType = TracePointStorage;
}

private struct ArrayStorage(T)
{
    enum elementSize = T.sizeof;
    size_t ptr;
    size_t length;

    ArrayStorage!T opIndex() const pure nothrow
    {
        return this;
    }

    size_t opIndex(size_t i) const pure nothrow
    {
        assert(i <= length, "out of bounds");

        return ptr + i * elementSize;
    }

    ArrayStorage!T opIndex(size_t[2] slice) const pure nothrow
    {
        auto from = slice[0];
        auto to = slice[1];
        assert(from < to && to <= length, "out of bounds");

        return typeof(return)(this[from], to - from);
    }

    size_t[2] opSlice(size_t dim)(size_t from, size_t to) const pure nothrow
            if (dim == 0)
    {
        assert(from < to, "slice `[from .. to]` must have `from < to`");

        return [from, to];
    }

    size_t opDollar() const pure nothrow
    {
        return length;
    }

    long indexOf(size_t ptr) const pure
    {
        assert((ptr.to!long - this.ptr.to!long) % elementSize.to!long == 0, "bad pointer alignment");

        return (ptr.to!long - this.ptr.to!long) / elementSize.to!long;
    }

    static ArrayStorage!T fromPtrs(size_t fromPtr, size_t toPtr)
    {
        assert((toPtr - fromPtr) % elementSize == 0, "bad pointer alignment");

        return typeof(return)(fromPtr, (toPtr - fromPtr) / elementSize);
    }
}

unittest
{
    alias T = int;
    static assert(ArrayStorage!T.elementSize == T.sizeof);
    enum basePtr = 1337;
    enum length = 42;

    ArrayStorage!T storage;
    storage.ptr = basePtr;
    storage.length = length;

    assert(storage[0] == basePtr);
    assert(storage[13] == basePtr + 13 * T.sizeof);
    assert(storage[$] == storage[length]);
    assertThrown!AssertError(storage[length + 1]);

    auto storageCopy = storage[];
    assert(storageCopy == storage);
    --storageCopy.length;
    assert(storageCopy.length != storage.length);

    storageCopy = storage[0 .. $];
    assert(storageCopy == storage);
    --storageCopy.length;
    assert(storageCopy.length != storage.length);

    auto storageSlice = storage[1 .. 13];
    assert(storageSlice.ptr == storage[1]);
    assert(storageSlice.length == 12);

    assert(ArrayStorage!T.fromPtrs(storage[1], storage[13]) == storageSlice);

    assert(storage.indexOf(basePtr + 13 * T.sizeof) == 13);
    assert(storage.indexOf(basePtr - 13 * T.sizeof) == -13);
    assertThrown!AssertError(storage.indexOf(basePtr + 1));
}

private struct SeededAlignmentStorage
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;

    id_t id;
    id_t contigAId;
    coord_t contigALength;
    id_t contigBId;
    coord_t contigBLength;
    AlignmentChain.Flags flags;
    StorageType!(LocalAlignment[]) localAlignments;
    trace_point_t tracePointDistance;
    AlignmentLocationSeed seed;
}

private struct LocalAlignmentStorage
{
    alias TracePoint = AlignmentChain.LocalAlignment.TracePoint;

    coord_t contigABegin;
    coord_t contigAEnd;
    coord_t contigBBegin;
    coord_t contigBEnd;
    diff_t numDiffs;
    StorageType!(TracePoint[]) tracePoints;
}

private struct TracePointStorage
{
    trace_point_t numDiffs;
    trace_point_t numBasePairs;
}

version (unittest)
{
    enum ushort defaultTracePointDistance = 100;
    PileUp[] getPileUpsTestData()
    {
        with (AlignmentChain) with (LocalAlignment) with (Flag)
            return [
                [
                    ReadAlignment(
                        SeededAlignment(
                            AlignmentChain(
                                9,
                                Contig(1, 8300),
                                Contig(1539, 6414),
                                Flags(complement),
                                [
                                    LocalAlignment(
                                        Locus(0, 6227),
                                        Locus(7, 6414),
                                        309,
                                        [
                                            TracePoint(10, 98),
                                            TracePoint(7, 107),
                                            TracePoint(7, 105),
                                            TracePoint(8, 108),
                                            TracePoint(4, 99),
                                            TracePoint(4, 102),
                                            TracePoint(5, 103),
                                            TracePoint(5, 105),
                                            TracePoint(1, 101),
                                            TracePoint(8, 107),
                                            TracePoint(6, 102),
                                            TracePoint(8, 103),
                                            TracePoint(7, 105),
                                            TracePoint(6, 105),
                                            TracePoint(2, 102),
                                            TracePoint(4, 104),
                                            TracePoint(5, 101),
                                            TracePoint(5, 99),
                                            TracePoint(0, 100),
                                            TracePoint(0, 100),
                                            TracePoint(10, 107),
                                            TracePoint(7, 103),
                                            TracePoint(5, 102),
                                            TracePoint(4, 101),
                                            TracePoint(2, 100),
                                            TracePoint(7, 106),
                                            TracePoint(1, 101),
                                            TracePoint(5, 101),
                                            TracePoint(7, 107),
                                            TracePoint(6, 106),
                                            TracePoint(3, 103),
                                            TracePoint(3, 103),
                                            TracePoint(7, 105),
                                            TracePoint(4, 103),
                                            TracePoint(5, 102),
                                            TracePoint(7, 103),
                                            TracePoint(7, 104),
                                            TracePoint(5, 105),
                                            TracePoint(5, 101),
                                            TracePoint(8, 106),
                                            TracePoint(4, 102),
                                            TracePoint(9, 104),
                                            TracePoint(2, 100),
                                            TracePoint(6, 102),
                                            TracePoint(8, 104),
                                            TracePoint(8, 106),
                                            TracePoint(6, 102),
                                            TracePoint(4, 102),
                                            TracePoint(2, 101),
                                            TracePoint(4, 102),
                                            TracePoint(4, 103),
                                            TracePoint(3, 103),
                                            TracePoint(4, 101),
                                            TracePoint(8, 108),
                                            TracePoint(2, 102),
                                            TracePoint(2, 101),
                                            TracePoint(3, 103),
                                            TracePoint(2, 100),
                                            TracePoint(3, 103),
                                            TracePoint(5, 102),
                                            TracePoint(6, 105),
                                            TracePoint(3, 98),
                                            TracePoint(1, 28),
                                        ],
                                    ),
                                ],
                                defaultTracePointDistance,
                            ),
                            AlignmentLocationSeed.front,
                        ),
                    ),
                    ReadAlignment(
                        SeededAlignment(
                            AlignmentChain(
                                17,
                                Contig(1, 8300),
                                Contig(3197, 8564),
                                Flags(complement),
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
                                        Locus(0, 8530),
                                        407,
                                        [
                                            TracePoint(6, 105),
                                            TracePoint(9, 108),
                                            TracePoint(7, 105),
                                            TracePoint(3, 103),
                                            TracePoint(2, 100),
                                            TracePoint(5, 100),
                                            TracePoint(6, 104),
                                            TracePoint(7, 104),
                                            TracePoint(5, 99),
                                            TracePoint(4, 104),
                                            TracePoint(2, 101),
                                            TracePoint(3, 103),
                                            TracePoint(6, 102),
                                            TracePoint(7, 102),
                                            TracePoint(8, 102),
                                            TracePoint(4, 104),
                                            TracePoint(4, 104),
                                            TracePoint(6, 99),
                                            TracePoint(7, 104),
                                            TracePoint(6, 105),
                                            TracePoint(3, 99),
                                            TracePoint(5, 102),
                                            TracePoint(3, 103),
                                            TracePoint(5, 103),
                                            TracePoint(4, 104),
                                            TracePoint(7, 104),
                                            TracePoint(6, 106),
                                            TracePoint(6, 105),
                                            TracePoint(7, 105),
                                            TracePoint(2, 100),
                                            TracePoint(2, 101),
                                            TracePoint(9, 101),
                                            TracePoint(6, 104),
                                            TracePoint(4, 102),
                                            TracePoint(6, 101),
                                            TracePoint(7, 104),
                                            TracePoint(8, 103),
                                            TracePoint(7, 107),
                                            TracePoint(3, 103),
                                            TracePoint(5, 103),
                                            TracePoint(3, 103),
                                            TracePoint(3, 101),
                                            TracePoint(7, 104),
                                            TracePoint(1, 101),
                                            TracePoint(8, 103),
                                            TracePoint(3, 103),
                                            TracePoint(1, 101),
                                            TracePoint(4, 103),
                                            TracePoint(4, 103),
                                            TracePoint(5, 104),
                                            TracePoint(4, 99),
                                            TracePoint(3, 103),
                                            TracePoint(5, 103),
                                            TracePoint(4, 102),
                                            TracePoint(3, 99),
                                            TracePoint(4, 98),
                                            TracePoint(8, 108),
                                            TracePoint(4, 104),
                                            TracePoint(6, 106),
                                            TracePoint(7, 106),
                                            TracePoint(7, 105),
                                            TracePoint(3, 103),
                                            TracePoint(4, 101),
                                            TracePoint(1, 101),
                                            TracePoint(6, 102),
                                            TracePoint(7, 104),
                                            TracePoint(1, 101),
                                            TracePoint(5, 98),
                                            TracePoint(8, 104),
                                            TracePoint(5, 104),
                                            TracePoint(4, 104),
                                            TracePoint(5, 102),
                                            TracePoint(4, 103),
                                            TracePoint(4, 100),
                                            TracePoint(3, 103),
                                            TracePoint(6, 104),
                                            TracePoint(7, 103),
                                            TracePoint(4, 104),
                                            TracePoint(5, 103),
                                            TracePoint(5, 102),
                                            TracePoint(2, 100),
                                            TracePoint(7, 102),
                                            TracePoint(5, 105),
                                        ],
                                    ),
                                ],
                                defaultTracePointDistance,
                            ),
                            AlignmentLocationSeed.front,
                        ),
                    ),
                ],
                [
                    ReadAlignment(
                        SeededAlignment(
                            AlignmentChain(
                                42,
                                Contig(1, 8300),
                                Contig(2325, 16069),
                                emptyFlags,
                                [
                                    LocalAlignment(
                                        Locus(3562, 8300),
                                        Locus(0, 4860),
                                        229,
                                        [
                                            TracePoint(3, 41),
                                            TracePoint(7, 102),
                                            TracePoint(2, 102),
                                            TracePoint(3, 103),
                                            TracePoint(4, 102),
                                            TracePoint(8, 108),
                                            TracePoint(2, 100),
                                            TracePoint(6, 102),
                                            TracePoint(6, 104),
                                            TracePoint(2, 102),
                                            TracePoint(4, 102),
                                            TracePoint(7, 102),
                                            TracePoint(3, 99),
                                            TracePoint(2, 102),
                                            TracePoint(3, 100),
                                            TracePoint(6, 100),
                                            TracePoint(3, 102),
                                            TracePoint(5, 104),
                                            TracePoint(5, 104),
                                            TracePoint(4, 101),
                                            TracePoint(5, 103),
                                            TracePoint(2, 102),
                                            TracePoint(7, 106),
                                            TracePoint(7, 103),
                                            TracePoint(3, 101),
                                            TracePoint(8, 106),
                                            TracePoint(6, 103),
                                            TracePoint(6, 103),
                                            TracePoint(4, 103),
                                            TracePoint(3, 102),
                                            TracePoint(2, 101),
                                            TracePoint(7, 106),
                                            TracePoint(4, 104),
                                            TracePoint(3, 101),
                                            TracePoint(8, 103),
                                            TracePoint(10, 102),
                                            TracePoint(8, 104),
                                            TracePoint(7, 102),
                                            TracePoint(3, 103),
                                            TracePoint(5, 105),
                                            TracePoint(6, 100),
                                            TracePoint(7, 105),
                                            TracePoint(3, 99),
                                            TracePoint(3, 100),
                                            TracePoint(5, 104),
                                            TracePoint(4, 104),
                                            TracePoint(3, 100),
                                            TracePoint(5, 103),
                                        ],
                                    ),
                                ],
                                defaultTracePointDistance,
                            ),
                            AlignmentLocationSeed.back,
                        ),
                        SeededAlignment(
                            AlignmentChain(
                                76,
                                Contig(2, 8350),
                                Contig(2325, 16069),
                                emptyFlags,
                                [
                                    LocalAlignment(
                                        Locus(0, 6798),
                                        Locus(9086, 16069),
                                        308,
                                        [
                                            TracePoint(5, 103),
                                            TracePoint(4, 101),
                                            TracePoint(8, 106),
                                            TracePoint(9, 104),
                                            TracePoint(9, 104),
                                            TracePoint(6, 106),
                                            TracePoint(3, 101),
                                            TracePoint(8, 104),
                                            TracePoint(4, 104),
                                            TracePoint(7, 101),
                                            TracePoint(4, 102),
                                            TracePoint(3, 100),
                                            TracePoint(7, 105),
                                            TracePoint(6, 99),
                                            TracePoint(3, 103),
                                            TracePoint(2, 102),
                                            TracePoint(6, 102),
                                            TracePoint(3, 103),
                                            TracePoint(4, 102),
                                            TracePoint(12, 110),
                                            TracePoint(7, 103),
                                            TracePoint(4, 104),
                                            TracePoint(4, 102),
                                            TracePoint(4, 100),
                                            TracePoint(5, 105),
                                            TracePoint(5, 100),
                                            TracePoint(6, 102),
                                            TracePoint(5, 105),
                                            TracePoint(3, 102),
                                            TracePoint(2, 102),
                                            TracePoint(6, 103),
                                            TracePoint(2, 100),
                                            TracePoint(4, 104),
                                            TracePoint(5, 103),
                                            TracePoint(6, 104),
                                            TracePoint(4, 104),
                                            TracePoint(3, 101),
                                            TracePoint(4, 104),
                                            TracePoint(6, 104),
                                            TracePoint(4, 104),
                                            TracePoint(1, 101),
                                            TracePoint(3, 103),
                                            TracePoint(8, 105),
                                            TracePoint(2, 102),
                                            TracePoint(4, 104),
                                            TracePoint(1, 99),
                                            TracePoint(6, 102),
                                            TracePoint(3, 101),
                                            TracePoint(3, 101),
                                            TracePoint(5, 105),
                                            TracePoint(1, 101),
                                            TracePoint(1, 100),
                                            TracePoint(4, 102),
                                            TracePoint(4, 101),
                                            TracePoint(3, 103),
                                            TracePoint(2, 101),
                                            TracePoint(5, 101),
                                            TracePoint(6, 103),
                                            TracePoint(6, 105),
                                            TracePoint(7, 106),
                                            TracePoint(8, 107),
                                            TracePoint(5, 103),
                                            TracePoint(5, 105),
                                            TracePoint(2, 102),
                                            TracePoint(2, 100),
                                            TracePoint(1, 100),
                                            TracePoint(5, 103),
                                            TracePoint(3, 99),
                                        ],
                                    ),
                                ],
                                defaultTracePointDistance,
                            ),
                            AlignmentLocationSeed.front,
                        ),
                   ),
                    ReadAlignment(
                        SeededAlignment(
                            AlignmentChain(
                                47,
                                Contig(1, 8300),
                                Contig(3332, 12946),
                                Flags(complement),
                                [
                                    LocalAlignment(
                                        Locus(4899, 8300),
                                        Locus(0, 3507),
                                        154,
                                        [
                                            TracePoint(0, 1),
                                            TracePoint(7, 105),
                                            TracePoint(6, 105),
                                            TracePoint(6, 106),
                                            TracePoint(7, 107),
                                            TracePoint(3, 102),
                                            TracePoint(7, 105),
                                            TracePoint(3, 103),
                                            TracePoint(2, 102),
                                            TracePoint(4, 100),
                                            TracePoint(3, 101),
                                            TracePoint(5, 105),
                                            TracePoint(3, 103),
                                            TracePoint(7, 101),
                                            TracePoint(3, 101),
                                            TracePoint(5, 102),
                                            TracePoint(7, 105),
                                            TracePoint(4, 104),
                                            TracePoint(3, 103),
                                            TracePoint(4, 104),
                                            TracePoint(3, 103),
                                            TracePoint(6, 102),
                                            TracePoint(4, 101),
                                            TracePoint(6, 104),
                                            TracePoint(2, 102),
                                            TracePoint(8, 104),
                                            TracePoint(4, 104),
                                            TracePoint(4, 103),
                                            TracePoint(1, 99),
                                            TracePoint(0, 100),
                                            TracePoint(9, 108),
                                            TracePoint(5, 104),
                                            TracePoint(5, 102),
                                            TracePoint(5, 103),
                                            TracePoint(3, 103),
                                        ],
                                    ),
                                ],
                                defaultTracePointDistance,
                            ),
                            AlignmentLocationSeed.back,
                        ),
                        SeededAlignment(
                            AlignmentChain(
                                89,
                                Contig(2, 8350),
                                Contig(3332, 12946),
                                Flags(complement),
                                [
                                    LocalAlignment(
                                        Locus(0, 5082),
                                        Locus(7726, 12946),
                                        244,
                                        [
                                            TracePoint(5, 104),
                                            TracePoint(2, 99),
                                            TracePoint(3, 103),
                                            TracePoint(7, 103),
                                            TracePoint(2, 98),
                                            TracePoint(2, 102),
                                            TracePoint(11, 107),
                                            TracePoint(4, 104),
                                            TracePoint(4, 102),
                                            TracePoint(8, 106),
                                            TracePoint(5, 105),
                                            TracePoint(11, 109),
                                            TracePoint(8, 101),
                                            TracePoint(5, 103),
                                            TracePoint(4, 104),
                                            TracePoint(2, 102),
                                            TracePoint(7, 105),
                                            TracePoint(3, 103),
                                            TracePoint(5, 101),
                                            TracePoint(3, 101),
                                            TracePoint(5, 102),
                                            TracePoint(4, 102),
                                            TracePoint(6, 103),
                                            TracePoint(5, 105),
                                            TracePoint(3, 99),
                                            TracePoint(5, 103),
                                            TracePoint(2, 100),
                                            TracePoint(7, 105),
                                            TracePoint(3, 101),
                                            TracePoint(5, 105),
                                            TracePoint(5, 103),
                                            TracePoint(5, 105),
                                            TracePoint(5, 105),
                                            TracePoint(7, 105),
                                            TracePoint(5, 102),
                                            TracePoint(3, 100),
                                            TracePoint(3, 102),
                                            TracePoint(4, 102),
                                            TracePoint(5, 103),
                                            TracePoint(5, 105),
                                            TracePoint(8, 100),
                                            TracePoint(6, 103),
                                            TracePoint(5, 102),
                                            TracePoint(6, 105),
                                            TracePoint(2, 98),
                                            TracePoint(2, 100),
                                            TracePoint(7, 104),
                                            TracePoint(5, 104),
                                            TracePoint(3, 103),
                                            TracePoint(2, 100),
                                            TracePoint(5, 82),
                                        ],
                                    ),
                                ],
                                defaultTracePointDistance,
                            ),
                            AlignmentLocationSeed.front,
                        )
                    ),
                    ReadAlignment(
                        SeededAlignment(
                            AlignmentChain(
                                43,
                                Contig(1, 8300),
                                Contig(3593, 6783),
                                emptyFlags,
                                [
                                    LocalAlignment(
                                        Locus(3943, 8300),
                                        Locus(0, 4471),
                                        213,
                                        [
                                            TracePoint(3, 56),
                                            TracePoint(7, 102),
                                            TracePoint(2, 102),
                                            TracePoint(7, 101),
                                            TracePoint(7, 103),
                                            TracePoint(9, 107),
                                            TracePoint(4, 104),
                                            TracePoint(4, 104),
                                            TracePoint(5, 100),
                                            TracePoint(4, 104),
                                            TracePoint(3, 101),
                                            TracePoint(5, 102),
                                            TracePoint(5, 101),
                                            TracePoint(2, 102),
                                            TracePoint(2, 98),
                                            TracePoint(6, 102),
                                            TracePoint(10, 109),
                                            TracePoint(1, 99),
                                            TracePoint(6, 102),
                                            TracePoint(1, 101),
                                            TracePoint(3, 101),
                                            TracePoint(5, 104),
                                            TracePoint(7, 101),
                                            TracePoint(7, 105),
                                            TracePoint(5, 102),
                                            TracePoint(0, 100),
                                            TracePoint(4, 104),
                                            TracePoint(10, 105),
                                            TracePoint(7, 103),
                                            TracePoint(11, 111),
                                            TracePoint(8, 105),
                                            TracePoint(4, 101),
                                            TracePoint(4, 103),
                                            TracePoint(3, 100),
                                            TracePoint(5, 104),
                                            TracePoint(4, 102),
                                            TracePoint(7, 105),
                                            TracePoint(5, 105),
                                            TracePoint(6, 104),
                                            TracePoint(5, 102),
                                            TracePoint(3, 103),
                                            TracePoint(3, 101),
                                            TracePoint(2, 100),
                                            TracePoint(2, 100),
                                        ],
                                    ),
                                ],
                                defaultTracePointDistance,
                            ),
                            AlignmentLocationSeed.back,
                        )
                    ),
                ],
            ];
    }
}

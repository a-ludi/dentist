/**
    This is he algorithm for building pile ups.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.collectPileUps.pileups;

import dentist.common.alignments :
    AlignmentChain,
    AlignmentLocationSeed,
    arithmetic_t,
    id_t,
    coord_t,
    getType,
    isValid,
    makeJoin,
    PileUp,
    ReadAlignment,
    SeededAlignment;
import dentist.common.binio : ArrayStorage;
import dentist.common.scaffold :
    buildScaffold,
    ContigNode,
    ContigPart,
    getUnkownJoin,
    isGap,
    isTranscendent,
    isReal,
    Join,
    mergeExtensionsWithGaps,
    removeNoneJoins,
    Scaffold;
import dentist.dazzler : GapSegment;
import dentist.util.algorithm : orderLexicographically;
import dentist.util.math :
    add,
    bulkAdd,
    filterEdges,
    mapEdges;
import dentist.util.log;
import std.algorithm :
    any,
    canFind,
    chunkBy,
    equal,
    filter,
    fold,
    joiner,
    map,
    min,
    minElement,
    sort,
    SwapStrategy;
import std.algorithm : equal;
import std.array : appender, array;
import std.bitmanip : bitsSet;
import std.conv : to;
import std.range :
    chain,
    enumerate,
    iota,
    only,
    slide,
    walkLength,
    zip;
import std.range.primitives;
import std.traits : EnumMembers;
import std.typecons :
    BitFlags,
    Flag,
    No,
    Yes;
import vibe.data.json : Json, toJson = serializeToJson;


private auto collectPileUps(Scaffold!ScaffoldPayload scaffold)
{
    return scaffold
        .edges
        .map!"a.payload"
        .filter!(payload => payload.types.pileUp)
        .map!(payload => payload.readAlignments)
        .filter!(pileUp => pileUp.length > 0)
        .filter!(pileUp => pileUp.isValid);
}

private void debugLogPileUps(string state, Scaffold!ScaffoldPayload scaffold)
{
    logJsonDebug(
        "state", state,
        "joins", scaffold
            .edges
            .map!joinToJson
            .array
            .toJson,
        "pileUps", collectPileUps(scaffold)
            .map!(pileUp => [
                "type": pileUp.getType.to!string.toJson,
                "readAlignments": pileUp.map!"a[]".array.toJson,
            ].toJson)
            .array
            .toJson,
    );
}

struct ScaffoldPayload
{
    static enum Type : ubyte
    {
        pileUp = 1 << 0,
        inputGap = 1 << 1,
    }

    BitFlags!Type types;
    ReadAlignment[] readAlignments;

    private this(Type type, ReadAlignment[] readAlignments = []) pure nothrow
    {
        this.types |= type;
        if (types.pileUp)
            this.readAlignments = readAlignments;
    }

    private this(BitFlags!Type types, ReadAlignment[] readAlignments) pure nothrow
    {
        this.types = types;
        this.readAlignments = readAlignments;
    }

    static ScaffoldPayload pileUp(ReadAlignment[] readAlignments) pure nothrow
    {
        return ScaffoldPayload(Type.pileUp, readAlignments);
    }

    void remove(Type type)() pure nothrow if (type == Type.pileUp)
    {
        types.pileUp = false;
        readAlignments = [];
    }

    static ScaffoldPayload inputGap() pure nothrow
    {
        return ScaffoldPayload(Type.inputGap);
    }

    void remove(Type type)() pure nothrow if (type == Type.inputGap)
    {
        types.inputGap = false;
    }

    static ScaffoldPayload merge(R)(R payloads)
        if (isInputRange!R && is(ElementType!R == ScaffoldPayload))
    {
        return ScaffoldPayload(
            payloads.save.map!"a.types".fold!"a | b",
            payloads.map!"a.readAlignments".joiner.array,
        );
    }

    static ScaffoldPayload merge(ScaffoldPayload[] payloads...)
    {
        return merge!(ScaffoldPayload[])(payloads);
    }

    /// Return the number of set types in this payload.
    size_t numTypes() const pure nothrow
    {
        return bitsSet(cast(ubyte) types).walkLength;
    }

    Json toJson()
    {
        return [
            "types": types.to!string.toJson,
            "pileUpType": types.pileUp && readAlignments.length > 0
                ? readAlignments.getType.to!string.toJson
                : Json(null),
            "readAlignments": shouldLog(LogLevel.debug_)
                ? readAlignments.map!"a[]".array.toJson
                : readAlignments.length.toJson,
        ].toJson;
    }
}

Join!ScaffoldPayload mergeJoins(Join!ScaffoldPayload[] joins...)
{
    assert(joins.length > 0);

    auto mergedJoin = joins[0];

    mergedJoin.payload = ScaffoldPayload.merge(joins.map!"a.payload");

    return mergedJoin;
}

static Join!ScaffoldPayload selectMeanest(Join!ScaffoldPayload[] joins...)
{
    return joins.minElement!"a.payload.numTypes";
}

PileUp[] build(Options)(
    in size_t numReferenceContigs,
    AlignmentChain[] candidates,
    GapSegment[] inputGaps,
    in Options options,
)
{
    alias isSameRead = (a, b) => a.contigB.id == b.contigB.id;
    alias ReadAlignmentJoin = Join!ScaffoldPayload;

    candidates.sort!("a.contigB.id < b.contigB.id", SwapStrategy.stable);

    auto readAlignmentJoins = candidates
        .filter!"!a.flags.disabled"
        .chunkBy!isSameRead
        .map!collectReadAlignments
        .joiner
        .filter!"a.isValid"
        .map!"a.getInOrder()"
        .map!makeScaffoldJoin;
    auto inputGapJoins = inputGaps
        .map!makeScaffoldJoin;

    auto alignmentsScaffold = buildScaffold!(mergeJoins, ScaffoldPayload)(
        numReferenceContigs + 0,
        chain(readAlignmentJoins, inputGapJoins),
    );
    debugLogPileUps("raw", alignmentsScaffold);
    alignmentsScaffold = alignmentsScaffold.discardAmbiguousJoins(
        options.bestPileUpMargin,
        options.existingGapBonus,
    );
    debugLogPileUps("unambiguous", alignmentsScaffold);
    alignmentsScaffold = alignmentsScaffold.enforceMinSpanningReads(options.minSpanningReads);
    debugLogPileUps("minSpanningEnforced", alignmentsScaffold);
    alignmentsScaffold = alignmentsScaffold.removeInputGaps();
    debugLogPileUps("inputGapsRemoved", alignmentsScaffold);
    alignmentsScaffold = alignmentsScaffold.mergeExtensionsWithGaps!(ScaffoldPayload.merge, ScaffoldPayload);
    debugLogPileUps("extensionsMerged", alignmentsScaffold);
    auto pileUps = collectPileUps(alignmentsScaffold).array;

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
            static struct Options
            {
                size_t minSpanningReads = 1;
                double bestPileUpMargin = 1.0;
                double existingGapBonus = 1.0;
            }

            id_t alignmentChainId = 0;
            id_t contReadId = 0;
            ReadAlignment getDummyRead(id_t beginContigId, arithmetic_t beginIdx,
                    id_t endContigId, arithmetic_t endIdx, Complement complement)
            {
                static enum contigLength = 16;
                static enum gapLength = 9;
                static enum numDiffs = 0;

                alignmentChainId += 2;
                auto readId = ++contReadId;
                auto readLength = beginContigId == endContigId
                    ? endIdx - beginIdx
                    : contigLength - beginIdx + gapLength + endIdx;
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
                }
                else
                {
                    auto secondReadBeginIdx = readLength - endIdx;
                    auto secondReadEndIdx = readLength;

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
                }
            }

            id_t c1 = 1;
            id_t c2 = 2;
            id_t c3 = 3;
            auto inputGaps = [
                GapSegment(
                    c1, // beginGlobalContigId
                    c2, // endGlobalContigId
                    0,  // scaffoldId
                    0,  // beginContigId
                    1,  // endContigId
                    15, // begin
                    25, // end
                ),
                GapSegment(
                    c2, // beginGlobalContigId
                    c3, // endGlobalContigId
                    0,  // scaffoldId
                    1,  // beginContigId
                    2,  // endContigId
                    40, // begin
                    50, // end
                ),
            ];
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
            auto computedPileUps = build(
                3,
                alignmentChains,
                inputGaps,
                Options(),
            );

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

/// Generate join from read alignment.
Join!ScaffoldPayload makeScaffoldJoin(ReadAlignment readAlignment)
{
    auto join = makeJoin!(typeof(return))(readAlignment);
    join.payload = ScaffoldPayload.pileUp([readAlignment]);

    return join;
}

///
unittest
{
    with (AlignmentChain) with (LocalAlignment) with (Flag)
            {
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
                auto inputGap = GapSegment(
                    1, // beginGlobalContigId
                    2, // endGlobalContigId
                    0,  // scaffoldId
                    0,  // beginContigId
                    1,  // endContigId
                    15, // begin
                    25, // end
                );

                auto join1 = makeScaffoldJoin(frontExtension);
                auto join2 = makeScaffoldJoin(backExtension);
                auto join3 = makeScaffoldJoin(gap);
                auto join4 = makeScaffoldJoin(inputGap);

                assert(join1.start == ContigNode(1, ContigPart.pre));
                assert(join1.end == ContigNode(1, ContigPart.begin));
                assert(join1.payload == ScaffoldPayload.pileUp([frontExtension]));

                assert(join2.start == ContigNode(1, ContigPart.end));
                assert(join2.end == ContigNode(1, ContigPart.post));
                assert(join2.payload == ScaffoldPayload.pileUp([backExtension]));

                assert(join3.start == ContigNode(1, ContigPart.end));
                assert(join3.end == ContigNode(2, ContigPart.end));
                assert(join3.payload == ScaffoldPayload.pileUp([gap]));

                assert(join4.start == ContigNode(1, ContigPart.end));
                assert(join4.end == ContigNode(2, ContigPart.begin));
                assert(join4.payload == ScaffoldPayload.inputGap());
            }
}

/// Generate join from inputGap.
Join!ScaffoldPayload makeScaffoldJoin(GapSegment inputGap)
{
    return typeof(return)(
        ContigNode(
            inputGap.beginGlobalContigId,
            ContigPart.end,
        ),
        ContigNode(
            inputGap.endGlobalContigId,
            ContigPart.begin,
        ),
        ScaffoldPayload.inputGap,
    );
}

// Not meant for public usage.
ReadAlignment[] collectReadAlignments(Chunk)(Chunk sameReadAlignments)
{
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
    alias seededPartsOfOneAlignment = (a, b) => a.alignment == b.alignment && a.seed != b.seed;
    alias shareReadSequence = (a, b) => endRelToContigB(a) > beginRelToContigB(b);
    alias emptyRange = () => cast(ReadAlignment[])[];

    auto seededAlignments = sameReadAlignments.map!(SeededAlignment.from).joiner.array;
    seededAlignments.sort!orderByLocusAndSeed;

    if (seededAlignments.length == 0)
    {
        return emptyRange();
    }

    // Validate seeded alignments
    version (assert)
        foreach (saPair; seededAlignments.slide!(No.withPartial)(2))
            assert(
                (
                    !shareReadSequence(saPair[0], saPair[1]) ||
                    seededPartsOfOneAlignment(saPair[0], saPair[1])
                ),
                "illegal overlap between seeded alignments",
            );

    // Collect read alignments
    bool startWithExtension = beginRelToContigB(seededAlignments[0]) > 0;
    size_t sliceStart = startWithExtension ? 1 : 0;

    auto remainingReadAlignments = iota(sliceStart, seededAlignments.length, 2)
        .map!(i => ReadAlignment(seededAlignments[i .. min(i + 2, $)]));
    auto readAlignments = startWithExtension
        ? chain(
            only(ReadAlignment(seededAlignments[0 .. 1])),
            remainingReadAlignments,
        ).array
        : remainingReadAlignments.array;

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
    enum emptyFlags = AlignmentChain.emptyFlags;
    enum complement = AlignmentChain.Flag.complement;
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias Locus = LocalAlignment.Locus;
    alias Seed = AlignmentLocationSeed;

    {
        // Case 1:
        //
        // |-->-->-->--|  |-->-->-->--| |-->-->-->--|
        //
        //       |----->  <-----------> <-----|
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
    }
    {
        // Case 2:
        //
        // |-->-->-->--|  |--<--<--<--| |-->-->-->--|
        //
        //       |----->  <-----------> <-----|
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
    }
    {
        // Case 3:
        //
        //                |--<--<--<--| |-->-->-->--|
        //
        //                <-----------> <-----|
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
    }
    {
        // Case 4:
        //
        // |-->-->-->--|  |--<--<--<--|
        //
        //       |----->  <----------->
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
    }
    {
        // Case 5:
        //
        //                |--<--<--<--|
        //
        //                <----------->
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
    }
}


/// This removes ambiguous gap insertions.
Scaffold!ScaffoldPayload discardAmbiguousJoins(
    Scaffold!ScaffoldPayload scaffold,
    in double bestPileUpMargin,
    in double existingGapBonus,
)
{
    auto incidentEdgesCache = scaffold.allIncidentEdges();
    auto removePileUpsAcc = appender!(Join!ScaffoldPayload[]);

    foreach (contigNode; scaffold.nodes)
    {
        assert(!contigNode.contigPart.isTranscendent || incidentEdgesCache[contigNode].length <= 1);

        if (contigNode.contigPart.isReal && incidentEdgesCache[contigNode].length > 2)
        {
            auto incidentGapJoins = incidentEdgesCache[contigNode]
                .filter!isGap
                .filter!"a.payload.types.pileUp"
                .array;

            if (incidentGapJoins.length > 1)
            {
                auto correctGapJoinIdx = incidentGapJoins.findCorrectGapJoin(
                    bestPileUpMargin,
                    existingGapBonus,
                );

                if (correctGapJoinIdx < incidentGapJoins.length)
                {
                    // Keep correct gap join for diagnostic output
                    auto correctGapJoin = incidentGapJoins[correctGapJoinIdx];

                    // Remove correct gap join from the list
                    incidentGapJoins = incidentGapJoins[0 .. correctGapJoinIdx] ~
                                       incidentGapJoins[correctGapJoinIdx + 1 .. $];

                    logJsonDiagnostic(
                        "info", "removing bad gap pile ups",
                        "sourceContigNode", contigNode.toJson,
                        "correctGapJoin", joinToJson(correctGapJoin),
                        "removedGapJoins", incidentGapJoins.map!joinToJson.array.toJson,
                    );
                }
                else
                {
                    logJsonDiagnostic(
                        "info", "skipping ambiguous gap pile ups",
                        "sourceContigNode", contigNode.toJson,
                        "removedGapJoins", incidentGapJoins.map!joinToJson.array.toJson,
                    );
                }

                // Mark bad/ambiguous pile ups for removal
                removePileUpsAcc ~= incidentGapJoins.filter!"a.payload.types.pileUp";
            }
        }
    }

    foreach (ref gapJoin; removePileUpsAcc.data)
        gapJoin.payload.remove!(ScaffoldPayload.Type.pileUp);

    scaffold.bulkAdd!selectMeanest(removePileUpsAcc.data);

    return removeNoneJoins!ScaffoldPayload(scaffold);
}

///
unittest
{
    ScaffoldPayload getDummyPayload(in size_t length, Flag!"hasInputGap" hasInputGap = No.hasInputGap)
    {
        auto payload = ScaffoldPayload.pileUp(new ReadAlignment[length]);

        if (hasInputGap)
            payload.types.inputGap = true;

        return payload;
    }

    //             contig 1      contig 2
    //
    //            o        o     o        o
    //                    / e1 e2 \      / e4
    //              o -> o ------- o -> o
    //               \        e3         \
    //                \                   \ e5 (strong evidence)
    // (input gap) e10 \   ____________   /
    //                  \ /   e6       \ /
    //    o <- o         o <- o         o <- o
    //                         \ e7 e8 /      \ e9
    //  o        o     o        o     o        o
    //
    //   contig 5       contig 4      contig 3
    //
    alias J = Join!ScaffoldPayload;
    alias S = Scaffold!ScaffoldPayload;
    alias CN = ContigNode;
    alias CP = ContigPart;
    auto scaffold = buildScaffold!(mergeJoins, ScaffoldPayload)(5, [
        J(CN(1, CP.end), CN(1, CP.post ), getDummyPayload(1)), // e1
        J(CN(1, CP.end), CN(1, CP.post ), getDummyPayload(1)), // e1
        J(CN(2, CP.pre), CN(2, CP.begin), getDummyPayload(1)), // e2
        J(CN(1, CP.end), CN(2, CP.begin), getDummyPayload(1)), // e3
        J(CN(2, CP.end), CN(2, CP.post ), getDummyPayload(1)), // e4
        J(CN(2, CP.end), CN(3, CP.end  ), getDummyPayload(2)), // e5
        J(CN(4, CP.end), CN(3, CP.end  ), getDummyPayload(1)), // e6
        J(CN(4, CP.end), CN(4, CP.post ), getDummyPayload(1)), // e7
        J(CN(3, CP.pre), CN(3, CP.begin), getDummyPayload(1)), // e8
        J(CN(3, CP.end), CN(3, CP.post ), getDummyPayload(1)), // e9
        J(CN(4, CP.end), CN(1, CP.begin), getDummyPayload(1, Yes.hasInputGap)), // e10
    ]).discardAmbiguousJoins(1.5, 2.0);
    //
    //   contig 1      contig 2
    //
    //            o        o     o        o
    //                    / e1 e2 \      / e4
    //              o -- o ------- o -- o
    //               \        e3         \ e5
    //            e10 \                  /
    //    o -- o       o -- o           o -- o
    //                       \ e7   e8 /      \ e9
    //  o        o   o        o       o        o
    //
    //   contig 5       contig 4      contig 3

    assert(J(CN(1, CP.end), CN(1, CP.post)) in scaffold); // e1
    assert(J(CN(2, CP.pre), CN(2, CP.begin)) in scaffold); // e2
    assert(J(CN(1, CP.end), CN(2, CP.begin)) in scaffold); // e3
    assert(J(CN(2, CP.end), CN(2, CP.post)) in scaffold); // e4
    assert(J(CN(2, CP.end), CN(3, CP.end)) in scaffold); // e5
    assert(J(CN(4, CP.end), CN(3, CP.end)) !in scaffold); // e6
    assert(J(CN(4, CP.end), CN(4, CP.post)) in scaffold); // e7
    assert(J(CN(3, CP.pre), CN(3, CP.begin)) in scaffold); // e8
    assert(J(CN(3, CP.end), CN(3, CP.post)) in scaffold); // e9
    assert(J(CN(4, CP.end), CN(1, CP.begin)) in scaffold); // e10

    assert(scaffold.get(J(CN(1, CP.end), CN(1, CP.post))).payload == getDummyPayload(2)); // e1
}

// Not meant for public usage.
auto joinToJson(Join!ScaffoldPayload join)
{
    return [
        "start": join.start.toJson,
        "end": join.end.toJson,
        "payload": join.payload.toJson,
    ].toJson;
}

size_t findCorrectGapJoin(
    Join!ScaffoldPayload[] incidentGapJoins,
    in double bestPileUpMargin,
    in double existingGapBonus,
)
{
    auto pileUpsLengths = incidentGapJoins
        .map!(gapJoin => gapJoin.payload.readAlignments.length * (gapJoin.payload.types.inputGap
            ? existingGapBonus
            : 1.0
        ))
        .enumerate
        .array;
    pileUpsLengths.sort!"a.value > b.value";
    auto largestPileUp = pileUpsLengths[0];
    auto sndLargestPileUp = pileUpsLengths[1];

    if (sndLargestPileUp.value * bestPileUpMargin < largestPileUp.value)
        return largestPileUp.index;
    else
        return size_t.max;
}

Scaffold!ScaffoldPayload enforceMinSpanningReads(Scaffold!ScaffoldPayload scaffold, size_t minSpanningReads)
{
    auto enforceMinSpanningReadsOnJoin(ref Join!ScaffoldPayload join)
    {
        if (
            join.payload.types.pileUp &&
            join.isGap &&
            join.payload.readAlignments.length < minSpanningReads
        )
            join.payload.remove!(ScaffoldPayload.Type.pileUp)();

        return join;
    }

    scaffold.mapEdges!enforceMinSpanningReadsOnJoin;

    return removeNoneJoins!ScaffoldPayload(scaffold);
}

Scaffold!ScaffoldPayload removeInputGaps(Scaffold!ScaffoldPayload scaffold)
{
    auto removeInputGapsOnJoin(ref Join!ScaffoldPayload join)
    {
        join.payload.remove!(ScaffoldPayload.Type.inputGap)();

        return join;
    }

    scaffold.mapEdges!removeInputGapsOnJoin;

    return removeNoneJoins!ScaffoldPayload(scaffold);
}

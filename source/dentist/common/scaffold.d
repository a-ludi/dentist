/**
    Work with scaffold graphs.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.scaffold;

import dentist.common.alignments : getType;
import dentist.util.log;
import dentist.util.math :
    add,
    bulkAdd,
    filterEdges,
    Graph,
    MissingNodeException,
    NaturalNumberSet;
import std.algorithm :
    count,
    equal,
    filter,
    fold,
    joiner,
    map,
    minElement,
    setDifference,
    sort,
    sum;
import std.array : appender, array;
import std.conv : to;
import std.functional : binaryFun;
import std.range :
    enumerate,
    iota,
    isForwardRange,
    only,
    refRange,
    retro,
    save,
    walkLength;
import std.typecons : Flag, No, Tuple, Yes;
import vibe.data.json : toJson = serializeToJson;

debug
{
    import std.conv : to;
    import std.stdio : writeln;
}

/// Each contig has four designated parts where joins can start or end.
static enum ContigPart : ubyte
{
    /// Designates a transcendent point *before* the contig where
    /// front extensions end.
    pre,
    /// Designates the begin of the contig.
    begin,
    /// Designates the end of the contig.
    end,
    /// Designates a transcendent point *after* the contig where
    /// back extensions end.
    post,
}

bool isReal(in ContigPart contigPart) pure nothrow
{
    return contigPart == ContigPart.begin || contigPart == ContigPart.end;
}

bool isTranscendent(in ContigPart contigPart) pure nothrow
{
    return contigPart == ContigPart.pre || contigPart == ContigPart.post;
}

/**
    A contig is represented by four `ContigNodes` in the scaffold graph.

    See_Also: ContigPart
*/
alias ContigNode = Tuple!(size_t, "contigId", ContigPart, "contigPart");
/// Represents a set of joins which conclude possibly several scaffolding
/// variants.
alias Scaffold(T) = Graph!(ContigNode, void, No.isDirected, T);
alias Join(T) = Scaffold!T.Edge;

Join!T sumPayloads(T)(Join!T[] joins...) nothrow
{
    assert(joins.length > 0);

    auto mergedJoin = joins[0];

    mergedJoin.payload = joins
        .map!"a.payload"
        .sum;

    return mergedJoin;
}

Join!T concatenatePayloads(T)(Join!T[] joins...) nothrow
{
    assert(joins.length > 0);

    auto mergedJoin = joins[0];

    mergedJoin.payload = joins
        .map!"a.payload"
        .joiner
        .array;

    return mergedJoin;
}

J dontJoin(J)(J j1, J j2) pure nothrow
{
    j1.payload = T.init;

    return j1;
}

/// Returns true iff join is a default edge of the scaffold graph.
bool isDefault(J)(in J join) pure nothrow
{
    return join.start.contigPart == ContigPart.begin &&
           join.end.contigPart == ContigPart.end &&
           join.start.contigId == join.end.contigId;
}

/// Returns true iff join is a unknown edge, ie. an edge for unknown sequence
/// (`n`s) of the scaffold graph.
bool isUnkown(J)(in J join) pure nothrow
{
    return join.start.contigId != join.end.contigId &&
        (join.start.contigPart != join.end.contigPart) &&
        join.start.contigPart.isTranscendent &&
        join.end.contigPart.isTranscendent;
}

/// Returns true iff join is a gap edge of the scaffold graph.
bool isGap(J)(in J join) pure nothrow
{
    return join.start.contigId != join.end.contigId &&
        (join.start.contigPart.isReal) &&
        (join.end.contigPart.isReal);
}

/// Returns true iff join is a gap edge and anti-parallel.
bool isAntiParallel(J)(in J join) pure nothrow
{
    return join.isGap && join.start.contigPart == join.end.contigPart;
}

/// Returns true iff join is a gap edge and parallel.
bool isParallel(J)(in J join) pure nothrow
{
    return join.isGap && join.start.contigPart != join.end.contigPart;
}

/// Returns true iff join is an extension edge of the scaffold graph.
bool isExtension(J)(in J join) pure nothrow
{
    return isFrontExtension(join) ^ isBackExtension(join);
}

/// Returns true iff join is a front extension edge of the scaffold graph.
bool isFrontExtension(J)(in J join) pure nothrow
{
    return join.start.contigId == join.end.contigId &&
        join.start.contigPart == ContigPart.pre &&
        join.end.contigPart == ContigPart.begin;
}

/// Returns true iff join is a back extension edge of the scaffold graph.
bool isBackExtension(J)(in J join) pure nothrow
{
    return join.start.contigId == join.end.contigId &&
        join.start.contigPart == ContigPart.end &&
        join.end.contigPart == ContigPart.post;
}

/// Returns true iff join is a valid edge of the scaffold graph.
bool isValid(J)(in J join) pure nothrow
{
    return isDefault(join) ^ isGap(join) ^ isExtension(join) ^ isUnkown(join);
}

/// Build a scaffold graph using `rawJoins`. This creates default edges and
/// inserts the rawJoins.
Scaffold!T buildScaffold(alias mergeMultiEdges, T, R)(in size_t numReferenceContigs, R rawJoins)
{
    auto scaffold = initScaffold!T(numReferenceContigs)
        .addJoins!(mergeMultiEdges, T)(rawJoins);

    return scaffold;
}

/// Creates a scaffold with all the default edges. Optionally specify a
/// function that produces the payloads.
///
/// See_Also: `getDefaultJoin`
Scaffold!T initScaffold(alias getPayload, T)(in size_t numReferenceContigs)
{
    auto contigIds = iota(1, numReferenceContigs + 1);
    auto contigNodes = contigIds
        .map!(contigId => only(
            ContigNode(contigId, ContigPart.pre),
            ContigNode(contigId, ContigPart.begin),
            ContigNode(contigId, ContigPart.end),
            ContigNode(contigId, ContigPart.post),
        ))
        .joiner
        .array;
    auto initialScaffold = Scaffold!T(contigNodes);

    static if (__traits(compiles, getPayload is null) && (getPayload is null))
    {
        alias createDefaultJoin = getDefaultJoin!T;
    }
    else
    {
        alias createDefaultJoin = getDefaultJoin!(getPayload, T);
    }

    initialScaffold.bulkAddForce(contigIds.map!createDefaultJoin.array);

    return initialScaffold;
}

/// ditto
Scaffold!T initScaffold(T)(in size_t numReferenceContigs)
{
    return initScaffold!(null, T)(numReferenceContigs);
}

/**
    Get the default join for contigId. Initialize payload with
    `getPayload(contigId)` if given.
*/
Join!T getDefaultJoin(T)(size_t contigId) pure nothrow
{
    return Join!T(
        ContigNode(contigId, ContigPart.begin),
        ContigNode(contigId, ContigPart.end),
    );
}

/// ditto
Join!T getDefaultJoin(alias getPayload, T)(size_t contigId) pure nothrow
{
    return Join!T(
        ContigNode(contigId, ContigPart.begin),
        ContigNode(contigId, ContigPart.end),
        getPayload(contigId),
    );
}

private Scaffold!T addJoins(alias mergeMultiEdges, T, R)(Scaffold!T scaffold, R rawJoins) if (isForwardRange!R)
{
    version (assert)
    {
        foreach (ref join; rawJoins.save)
        {
            assert(join.isValid && !join.isDefault);
        }
    }

    scaffold.bulkAdd!mergeMultiEdges(rawJoins);

    return scaffold;
}

/// Get join for a stretch of unknown sequence (`n`s).
Join!T getUnkownJoin(T)(size_t preContigId, size_t postContigId, T payload) pure nothrow
{
    assert(preContigId != postContigId);
    return Join!T(
        ContigNode(preContigId, ContigPart.post),
        ContigNode(postContigId, ContigPart.pre),
        payload,
    );
}

/// Normalizes unknown joins such that they join contigs or are removed as
/// applicable.
Scaffold!T normalizeUnkownJoins(T)(Scaffold!T scaffold)
{
    mixin(traceExecution);

    auto numUnkownJoins = scaffold.edges.count!isUnkown;
    auto newJoins = appender!(Join!T[]);
    newJoins.reserve(numUnkownJoins);
    auto removalAcc = appender!(Join!T[]);
    removalAcc.reserve(numUnkownJoins);
    auto degreesCache = scaffold.allDegrees();

    foreach (unkownJoin; scaffold.edges.filter!isUnkown)
    {
        auto preContigId = unkownJoin.start.contigId;
        auto preContigEnd = ContigNode(preContigId, ContigPart.end);
        auto postContigId = unkownJoin.end.contigId;
        auto postContigBegin = ContigNode(postContigId, ContigPart.begin);

        bool isPreContigUnconnected = degreesCache[preContigEnd] == 1;
        bool hasPreContigExtension = scaffold.has(Join!T(
            preContigEnd,
            unkownJoin.start,
        ));
        bool hasPreContigGap = !isPreContigUnconnected && !hasPreContigExtension;
        bool isPostContigUnconnected = degreesCache[postContigBegin] == 1;
        bool hasPostContigExtension = scaffold.has(Join!T(
            unkownJoin.end,
            postContigBegin,
        ));
        bool hasPostContigGap = !isPostContigUnconnected && !hasPostContigExtension;

        if (isPreContigUnconnected && isPostContigUnconnected)
        {
            newJoins ~= Join!T(
                preContigEnd,
                postContigBegin,
                unkownJoin.payload,
            );

            unkownJoin.payload = T.init;
            removalAcc ~= unkownJoin;
        }
        else if (isPreContigUnconnected && hasPostContigExtension)
        {
            newJoins ~= Join!T(
                preContigEnd,
                unkownJoin.end,
                unkownJoin.payload,
            );

            unkownJoin.payload = T.init;
            removalAcc ~= unkownJoin;
        }
        else if (hasPreContigExtension && isPostContigUnconnected)
        {
            newJoins ~= Join!T(
                unkownJoin.start,
                postContigBegin,
                unkownJoin.payload,
            );

            unkownJoin.payload = T.init;
            removalAcc ~= unkownJoin;
        }
        else if (hasPreContigGap || hasPostContigGap)
        {
            unkownJoin.payload = T.init;
            removalAcc ~= unkownJoin;
        }
    }

    scaffold.bulkAddForce(newJoins.data);
    foreach (noneJoin; removalAcc.data)
    {
        scaffold.add!(scaffold.ConflictStrategy.replace)(noneJoin);
    }

    return removeNoneJoins!T(scaffold);
}

///
unittest
{
    alias J = Join!int;
    alias S = Scaffold!int;
    alias CN = ContigNode;
    alias CP = ContigPart;

    //  Case 1:
    //
    //      o        oxxxxo        o   =>   o        o    o        o
    //                                 =>
    //        o -- o        o -- o     =>     o -- oxxxxxxxxo -- o
    auto scaffold1 = buildScaffold!(sumPayloads!int, int)(2, [
        getUnkownJoin!int(1, 2, 1),
    ]).normalizeUnkownJoins!int;
    assert(getDefaultJoin!int(1) in scaffold1);
    assert(getDefaultJoin!int(2) in scaffold1);
    assert(J(CN(1, CP.end), CN(2, CP.begin)) in scaffold1);
    assert(scaffold1.edges.walkLength == 3);

    //  Case 2:
    //
    //      o        oxxxxo        o  =>  o        o   xxxo        o
    //                     \          =>              /    \
    //        o -- o        o -- o    =>    o -- oxxxx      o -- o
    auto scaffold2 = buildScaffold!(sumPayloads!int, int)(2, [
        getUnkownJoin!int(1, 2, 1),
        J(CN(2, CP.pre), CN(2, CP.begin), 1),
    ]).normalizeUnkownJoins!int;
    assert(getDefaultJoin!int(1) in scaffold2);
    assert(getDefaultJoin!int(2) in scaffold2);
    assert(J(CN(1, CP.end), CN(2, CP.pre)) in scaffold2);
    assert(J(CN(2, CP.pre), CN(2, CP.begin)) in scaffold2);
    assert(scaffold2.edges.walkLength == 4);

    //  Case 3:
    //
    //      o        oxxxxo        o  =>  o        o    o        o
    //                                =>
    //        o -- o        o -- o    =>    o -- o        o -- o
    //                      |         =>                  |
    //                      o -- o    =>                  o -- o
    //                                =>
    //                    o        o  =>                o        o
    auto scaffold3 = buildScaffold!(sumPayloads!int, int)(3, [
        getUnkownJoin!int(1, 2, 1),
        J(CN(2, CP.begin), CN(3, CP.begin), 1),
    ]).normalizeUnkownJoins!int;
    assert(getDefaultJoin!int(1) in scaffold3);
    assert(getDefaultJoin!int(2) in scaffold3);
    assert(getDefaultJoin!int(3) in scaffold3);
    assert(J(CN(2, CP.begin), CN(3, CP.begin)) in scaffold3);
    assert(scaffold3.edges.walkLength == 4);

    //  Case 4:
    //
    //      o        oxxxxo        o  =>  o        oxxx   o        o
    //              /                 =>          /    \
    //        o -- o        o -- o    =>    o -- o      xxxxo -- o
    auto scaffold4 = buildScaffold!(sumPayloads!int, int)(2, [
        getUnkownJoin!int(1, 2, 1),
        J(CN(1, CP.end), CN(1, CP.post), 1),
    ]).normalizeUnkownJoins!int;
    assert(getDefaultJoin!int(1) in scaffold4);
    assert(getDefaultJoin!int(2) in scaffold4);
    assert(J(CN(1, CP.post), CN(2, CP.begin)) in scaffold4);
    assert(J(CN(1, CP.end), CN(1, CP.post)) in scaffold4);
    assert(scaffold4.edges.walkLength == 4);

    //  Case 5:
    //
    //      o        oxxxxo        o
    //              /      \
    //        o -- o        o -- o
    auto scaffold5 = buildScaffold!(sumPayloads!int, int)(2, [
        getUnkownJoin!int(1, 2, 1),
        J(CN(1, CP.end), CN(1, CP.post), 1),
        J(CN(2, CP.pre), CN(2, CP.begin), 1),
    ]).normalizeUnkownJoins!int;
    assert(getDefaultJoin!int(1) in scaffold5);
    assert(getDefaultJoin!int(2) in scaffold5);
    assert(getUnkownJoin!int(1, 2, 1) in scaffold5);
    assert(J(CN(1, CP.end), CN(1, CP.post)) in scaffold5);
    assert(J(CN(2, CP.pre), CN(2, CP.begin)) in scaffold5);
    assert(scaffold5.edges.walkLength == 5);

    //  Case 6:
    //
    //      o        oxxxxo        o  =>  o        o    o        o
    //              /                 =>          /
    //        o -- o        o -- o    =>    o -- o        o -- o
    //                      |         =>                  |
    //                      o -- o    =>                  o -- o
    //                                =>
    //                    o        o  =>                o        o
    auto scaffold6 = buildScaffold!(sumPayloads!int, int)(3, [
        getUnkownJoin!int(1, 2, 1),
        J(CN(1, CP.end), CN(1, CP.post), 1),
        J(CN(2, CP.begin), CN(3, CP.begin), 1),
    ]).normalizeUnkownJoins!int;
    assert(getDefaultJoin!int(1) in scaffold6);
    assert(getDefaultJoin!int(2) in scaffold6);
    assert(getDefaultJoin!int(3) in scaffold6);
    assert(J(CN(1, CP.end), CN(1, CP.post)) in scaffold6);
    assert(J(CN(2, CP.begin), CN(3, CP.begin)) in scaffold6);
    assert(scaffold6.edges.walkLength == 5);

    //  Case 7:
    //
    //      o        oxxxxo        o  =>  o        o    o        o
    //                                =>
    //        o -- o        o -- o    =>    o -- o        o -- o
    //             |                  =>         |
    //        o -- o                  =>    o -- o
    //                                =>
    //      o        o                =>  o        o
    auto scaffold7 = buildScaffold!(sumPayloads!int, int)(3, [
        getUnkownJoin!int(1, 2, 1),
        J(CN(1, CP.end), CN(3, CP.end), 1),
    ]).normalizeUnkownJoins!int;
    assert(getDefaultJoin!int(1) in scaffold7);
    assert(getDefaultJoin!int(2) in scaffold7);
    assert(getDefaultJoin!int(3) in scaffold7);
    assert(J(CN(1, CP.end), CN(3, CP.end)) in scaffold7);
    assert(scaffold7.edges.walkLength == 4);

    //  Case 8:
    //
    //      o        oxxxxo        o  =>  o        o    o        o
    //                     \          =>                 \
    //        o -- o        o -- o    =>    o -- o        o -- o
    //             |                  =>         |
    //        o -- o                  =>    o -- o
    //                                =>
    //      o        o                =>  o        o
    auto scaffold8 = buildScaffold!(sumPayloads!int, int)(3, [
        getUnkownJoin!int(1, 2, 1),
        J(CN(1, CP.end), CN(3, CP.end), 1),
        J(CN(2, CP.pre), CN(2, CP.begin), 1),
    ]).normalizeUnkownJoins!int;
    assert(getDefaultJoin!int(1) in scaffold8);
    assert(getDefaultJoin!int(2) in scaffold8);
    assert(getDefaultJoin!int(3) in scaffold8);
    assert(J(CN(1, CP.end), CN(3, CP.end)) in scaffold8);
    assert(J(CN(2, CP.pre), CN(2, CP.begin)) in scaffold8);
    assert(scaffold8.edges.walkLength == 5);

    //  Case 9:
    //
    //      o        oxxxxo        o  =>  o        o    o        o
    //                                =>
    //        o -- o        o -- o    =>    o -- o        o -- o
    //             |        |         =>         |        |
    //             o ------ o         =>         o ------ o
    //                                =>
    //           o            o       =>       o            o
    auto scaffold9 = buildScaffold!(sumPayloads!int, int)(3, [
        getUnkownJoin!int(1, 2, 1),
        J(CN(1, CP.end), CN(3, CP.begin), 1),
        J(CN(2, CP.begin), CN(3, CP.end), 1),
    ]).normalizeUnkownJoins!int;
    assert(getDefaultJoin!int(1) in scaffold9);
    assert(getDefaultJoin!int(2) in scaffold9);
    assert(getDefaultJoin!int(3) in scaffold9);
    assert(J(CN(1, CP.end), CN(3, CP.begin)) in scaffold9);
    assert(J(CN(2, CP.begin), CN(3, CP.end)) in scaffold9);
    assert(scaffold9.edges.walkLength == 5);
}

/// Determine which kinds of joins are allowed.
enum JoinPolicy
{
    /// Only join gaps inside of scaffolds (marked by `n`s in FASTA).
    scaffoldGaps,
    /// Join gaps inside of scaffolds (marked by `n`s in FASTA) and try to
    /// join scaffolds.
    scaffolds,
    /// Break input into contigs and re-scaffold everything; maintains scaffold gaps where new
    /// scaffolds are consistent.
    contigs,
}

/// Enforce joinPolicy in scaffold.
Scaffold!T enforceJoinPolicy(T)(Scaffold!T scaffold, in JoinPolicy joinPolicy)
{
    mixin(traceExecution);

    if (joinPolicy == JoinPolicy.contigs)
    {
        // Just do nothing; this is the default mode of operation.
        return scaffold;
    }

    alias orderByNodes = Scaffold!T.orderByNodes;

    assert(joinPolicy == JoinPolicy.scaffoldGaps || joinPolicy == JoinPolicy.scaffolds);

    auto allowedJoins = scaffold
        .edges
        .filter!isUnkown
        .map!(join => only(
            Join!T(
                ContigNode(join.start.contigId, ContigPart.end),
                ContigNode(join.start.contigId, ContigPart.post),
            ),
            Join!T(
                ContigNode(join.start.contigId, ContigPart.end),
                ContigNode(join.end.contigId, ContigPart.begin),
            ),
            Join!T(
                ContigNode(join.end.contigId, ContigPart.pre),
                ContigNode(join.end.contigId, ContigPart.begin),
            ),
        ))
        .joiner;
    auto gapJoins = scaffold.edges.filter!isGap;
    auto forbiddenJoins = setDifference!orderByNodes(gapJoins, allowedJoins).array;

    // NOTE: joins between scaffolds will be re-included into the graph later
    //       if joinPolicy == scaffolds.
    foreach (forbiddenJoin; forbiddenJoins)
    {
        forbiddenJoin.payload = T.init;
        scaffold.add!(scaffold.ConflictStrategy.replace)(forbiddenJoin);
    }

    scaffold = removeNoneJoins!T(scaffold);

    if (joinPolicy == JoinPolicy.scaffolds)
    {
        scaffold = normalizeUnkownJoins!T(scaffold);

        alias validJoin = (candidate) => scaffold.degree(candidate.start) == 1
            && scaffold.degree(candidate.end) == 1;
        foreach (scaffoldJoin; forbiddenJoins.filter!validJoin)
        {
            scaffold.add(scaffoldJoin);
        }

        if (shouldLog(LogLevel.info))
        {
            forbiddenJoins = forbiddenJoins.filter!(j => !validJoin(j)).array;
        }
    }

    logJsonInfo("forbiddenJoins", forbiddenJoins
        .filter!(gapJoin => scaffold.degree(gapJoin.start) == 1 && scaffold.degree(gapJoin.end) == 1)
        .map!(join => [
            "start": join.start.toJson,
            "end": join.end.toJson,
        ])
        .array
        .toJson);

    return scaffold;
}

/// Enforce joinPolicy in scaffold.
Scaffold!T removeExtensions(T)(Scaffold!T scaffold)
{
    auto extensionJoins = scaffold.edges.filter!isExtension;

    foreach (extensionJoin; extensionJoins)
    {
        extensionJoin.payload = T.init;
        scaffold.add!(scaffold.ConflictStrategy.replace)(extensionJoin);
    }

    return removeNoneJoins!T(scaffold);
}

/// Enforce joinPolicy in scaffold.
Scaffold!T removeSpanning(T)(Scaffold!T scaffold)
{
    auto spanningJoins = scaffold.edges.filter!isGap;

    foreach (spanningJoin; spanningJoins)
    {
        spanningJoin.payload = T.init;
        scaffold.add!(scaffold.ConflictStrategy.replace)(spanningJoin);
    }

    return removeNoneJoins!T(scaffold);
}

/// Remove marked edges from the graph. This always keeps the default edges.
Scaffold!T removeNoneJoins(T)(Scaffold!T scaffold)
{
    scaffold.filterEdges!(noneJoinFilter!T);

    return scaffold;
}

bool noneJoinFilter(T)(Join!T join)
{
    return isDefault(join) || join.payload != T.init;
}

/// Remove extension edges were they coincide with a gap edge combining their
/// payload. This is intended to build pile ups with all reads that contribute
/// to each gap.
Scaffold!T mergeExtensionsWithGaps(alias mergePayloads, T)(Scaffold!T scaffold)
{
    foreach (contigNode; scaffold.nodes)
    {
        assert(!contigNode.contigPart.isTranscendent || scaffold.degree(contigNode) <= 1);
        debug auto incidentEdges = scaffold.incidentEdges(contigNode).array;
        assert(scaffold.degree(contigNode) <= 3, "node degree must be <= 3");

        if (contigNode.contigPart.isReal && scaffold.degree(contigNode) == 3)
        {
            auto incidentJoins = scaffold.incidentEdges(contigNode)
                .filter!(j => !isDefault(j)).array;
            assert(incidentJoins.length == 2);
            // The gap join has real `contigPart`s on both ends.
            int gapJoinIdx = incidentJoins[0].target(contigNode).contigPart.isReal ? 0 : 1;
            auto gapJoin = incidentJoins[gapJoinIdx];
            auto extensionJoin = incidentJoins[$ - gapJoinIdx - 1];

            gapJoin.payload = binaryFun!mergePayloads(gapJoin.payload, extensionJoin.payload);
            extensionJoin.payload = T.init;

            scaffold.add!(scaffold.ConflictStrategy.replace)(gapJoin);
            scaffold.add!(scaffold.ConflictStrategy.replace)(extensionJoin);
        }
    }

    return removeNoneJoins!T(scaffold);
}

///
unittest
{
    alias J = Join!int;
    alias S = Scaffold!int;
    alias CN = ContigNode;
    alias CP = ContigPart;
    //   contig 1      contig 2
    //
    //  o        o     o        o
    //          / e1 e2 \      / e4
    //    o -- o ------- o -- o
    //              e3
    //
    //    o -- o         o -- o         o -- o
    //                         \ e7 e8 /      \ e9
    //  o        o     o        o     o        o
    //
    //   contig 5       contig 4      contig 3
    //
    auto scaffold = buildScaffold!(sumPayloads!int, int)(5, [
        J(CN(1, CP.end), CN(1, CP.post ), 1), // e1
        J(CN(1, CP.end), CN(1, CP.post ), 1), // e1
        J(CN(2, CP.pre), CN(2, CP.begin), 1), // e2
        J(CN(1, CP.end), CN(2, CP.begin), 1), // e3
        J(CN(2, CP.end), CN(2, CP.post ), 1), // e4
        J(CN(4, CP.end), CN(4, CP.post ), 1), // e7
        J(CN(3, CP.pre), CN(3, CP.begin), 1), // e8
        J(CN(3, CP.end), CN(3, CP.post ), 1), // e9
    ]).mergeExtensionsWithGaps!("a + b", int);

    assert(J(CN(1, CP.end), CN(1, CP.post)) !in scaffold); // e1
    assert(J(CN(2, CP.pre), CN(2, CP.begin)) !in scaffold); // e2
    assert(J(CN(1, CP.end), CN(2, CP.begin)) in scaffold); // e3
    assert(J(CN(2, CP.end), CN(2, CP.post)) in scaffold); // e4
    assert(J(CN(4, CP.end), CN(4, CP.post)) in scaffold); // e7
    assert(J(CN(3, CP.pre), CN(3, CP.begin)) in scaffold); // e8
    assert(J(CN(3, CP.end), CN(3, CP.post)) in scaffold); // e9

    assert(scaffold.get(J(CN(1, CP.end), CN(2, CP.begin))).payload == 4); // merged 2 * e1 + e2 + e3
}

/**
    Performs a linear walk through a scaffold graph starting in startNode.
    A linear walk is a sequence of adjacent joins where no node is visited
    twice unless the graph is cyclic in which case the first node will appear
    twice. The implementation requires the graph to have linear components,
    ie. for every node the degree must be at most two. If the component of
    `startNode` is cyclic then the walk will in `startNode` and the `isCyclic`
    flag will be set.

    The direction of the walk can be influenced by giving `firstJoin`.

    **Note:** if one wants to read the `isCyclic` flag it is required to use
    `std.range.refRange` in most cases.

    Returns: range of joins in the scaffold graph.
    Throws: MissingNodeException if any node is encountered that is not part
            of the graph.
*/
LinearWalk!T linearWalk(T)(
    Scaffold!T scaffold,
    ContigNode startNode,
    Scaffold!T.IncidentEdgesCache incidentEdgesCache = Scaffold!T.IncidentEdgesCache.init,
)
{
    return LinearWalk!T(scaffold, startNode, incidentEdgesCache);
}

/// ditto
LinearWalk!T linearWalk(T)(
    Scaffold!T scaffold,
    ContigNode startNode,
    Join!T firstJoin,
    Scaffold!T.IncidentEdgesCache incidentEdgesCache = Scaffold!T.IncidentEdgesCache.init,
)
{
    return LinearWalk!T(scaffold, startNode, firstJoin, incidentEdgesCache);
}

///
unittest
{
    alias Payload = int;
    alias J = Join!Payload;
    alias S = Scaffold!Payload;
    alias CN = ContigNode;
    alias CP = ContigPart;
    //   contig 1      contig 2
    //
    //  o        o     o        o
    //                         / e4
    //    o -- o ------- o -- o
    //              e3
    //
    //    o -- o         o -- o         o -- o
    //                         \ e7 e8 /      \ e9
    //  o        o     o        o     o        o
    //
    //   contig 5       contig 4      contig 3
    //
    auto joins1 = [
        J(CN(1, CP.end), CN(2, CP.begin), 1), // e3
        J(CN(2, CP.end), CN(2, CP.post ), 1), // e4
        J(CN(4, CP.end), CN(4, CP.post ), 1), // e7
        J(CN(3, CP.pre), CN(3, CP.begin), 1), // e8
        J(CN(3, CP.end), CN(3, CP.post ), 1), // e9
    ];
    auto scaffold1 = buildScaffold!(sumPayloads!Payload, Payload)(5, joins1);
    auto scaffold1Cache = scaffold1.allIncidentEdges();
    auto walks1 = [
        [
            getDefaultJoin!Payload(1),
            joins1[0],
            getDefaultJoin!Payload(2),
            joins1[1],
        ],
        [
            getDefaultJoin!Payload(4),
            joins1[2],
        ],
        [
            joins1[3],
            getDefaultJoin!Payload(3),
            joins1[4],
        ],
    ];

    alias getWalkStart = (walk) => walk[0].source(walk[0].getConnectingNode(walk[1]));

    foreach (walk; walks1)
    {
        auto reverseWalk = walk.retro.array;
        auto computedWalk = linearWalk!Payload(scaffold1, getWalkStart(walk));
        auto computedReverseWalk = linearWalk!Payload(scaffold1, getWalkStart(reverseWalk));

        assert(equal(walk[], refRange(&computedWalk)));
        assert(!computedWalk.isCyclic);
        assert(equal(reverseWalk[], refRange(&computedReverseWalk)));
        assert(!computedReverseWalk.isCyclic);
    }

    foreach (walk; walks1)
    {
        auto reverseWalk = walk.retro.array;
        auto computedWalk = linearWalk!Payload(scaffold1, getWalkStart(walk), scaffold1Cache);
        auto computedReverseWalk = linearWalk!Payload(scaffold1, getWalkStart(reverseWalk), scaffold1Cache);

        assert(equal(walk[], refRange(&computedWalk)));
        assert(!computedWalk.isCyclic);
        assert(equal(reverseWalk[], refRange(&computedReverseWalk)));
        assert(!computedReverseWalk.isCyclic);
    }

    //   contig 1      contig 2
    //
    //  o        o     o        o
    //
    //              e1
    //    o -- o ------- o -- o
    //     \_________________/
    //              e2
    //
    auto joins2 = [
        J(CN(1, CP.end), CN(2, CP.begin), 1), // e1
        J(CN(2, CP.end), CN(1, CP.begin ), 1), // e2
    ];
    auto scaffold2 = buildScaffold!(sumPayloads!Payload, Payload)(2, joins2);
    auto scaffold2Cache = scaffold2.allIncidentEdges();
    auto walk2 = [
        getDefaultJoin!Payload(1),
        joins2[0],
        getDefaultJoin!Payload(2),
        joins2[1],
    ];

    {
        auto computedWalk = linearWalk!Payload(scaffold2, getWalkStart(walk2), walk2[0]);
        auto reverseWalk2 = walk2.retro.array;
        auto computedReverseWalk2 = linearWalk!Payload(scaffold2,
                getWalkStart(reverseWalk2), reverseWalk2[0]);

        assert(equal(walk2[], refRange(&computedWalk)));
        assert(computedWalk.isCyclic);
        assert(equal(reverseWalk2[], refRange(&computedReverseWalk2)));
        assert(computedWalk.isCyclic);
    }

    {
        auto computedWalk = linearWalk!Payload(scaffold2, getWalkStart(walk2), walk2[0], scaffold2Cache);
        auto reverseWalk2 = walk2.retro.array;
        auto computedReverseWalk2 = linearWalk!Payload(scaffold2,
                getWalkStart(reverseWalk2), reverseWalk2[0], scaffold2Cache);

        assert(equal(walk2[], refRange(&computedWalk)));
        assert(computedWalk.isCyclic);
        assert(equal(reverseWalk2[], refRange(&computedReverseWalk2)));
        assert(computedWalk.isCyclic);
    }
}

struct LinearWalk(T)
{
    alias IncidentEdgesCache = Scaffold!T.IncidentEdgesCache;
    static enum emptyIncidentEdgesCache = IncidentEdgesCache.init;

    private Scaffold!T scaffold;
    private IncidentEdgesCache incidentEdgesCache;
    private size_t currentNodeIdx;
    private Join!T currentJoin;
    private bool isEmpty = false;
    private Flag!"isCyclic" isCyclic = No.isCyclic;
    private NaturalNumberSet visitedNodes;

    /// Start linear walk through a scaffold graph in startNode.
    this(
        Scaffold!T scaffold,
        ContigNode startNode,
        IncidentEdgesCache incidentEdgesCache = emptyIncidentEdgesCache,
    )
    {
        this.scaffold = scaffold;
        this.incidentEdgesCache = incidentEdgesCache == emptyIncidentEdgesCache
            ? scaffold.allIncidentEdges()
            : incidentEdgesCache;
        this.currentNode = startNode;
        this.visitedNodes.reserveFor(this.scaffold.nodes.length - 1);
        this.markVisited(this.currentNodeIdx);
        this.popFront();
    }

    /// Start linear walk through a scaffold graph in startNode.
    this(
        Scaffold!T scaffold,
        ContigNode startNode,
        Join!T firstJoin,
        IncidentEdgesCache incidentEdgesCache = emptyIncidentEdgesCache,
    )
    {
        this.scaffold = scaffold;
        this.incidentEdgesCache = incidentEdgesCache == emptyIncidentEdgesCache
            ? scaffold.allIncidentEdges()
            : incidentEdgesCache;
        this.currentNode = startNode;
        this.visitedNodes.reserveFor(this.scaffold.nodes.length - 1);
        this.markVisited(this.currentNodeIdx);
        this.currentJoin = firstJoin;
        currentNode = currentJoin.target(currentNode);
        this.markVisited(this.currentNodeIdx);
    }

    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty LinearWalk");
        assert(scaffold.degree(currentNode) <= 2, "fork in linear walk");

        // If isCyclic is set then the last edge of the cycle is being popped.
        // Thus we stop.
        if (isCyclic)
        {
            return endOfWalk();
        }

        auto candidateEdges = incidentEdgesCache[currentNode]
            .filter!(join => !visitedNodes.has(scaffold.indexOf(join.target(currentNode))));
        bool noSuccessorNodes = candidateEdges.empty;

        if (noSuccessorNodes)
        {
            // Iff the current node has more than one neighbor the graph has
            // a cycle.
            if (scaffold.degree(currentNode) > 1)
            {
                assert(scaffold.degree(currentNode) == 2);

                return lastEdgeOfCycle();
            }
            else
            {
                return endOfWalk();
            }
        }

        currentJoin = candidateEdges.front;
        currentNode = currentJoin.target(currentNode);
        markVisited(currentNodeIdx);
    }

    @property Join!T front()
    {
        assert(!empty, "Attempting to fetch the front of an empty LinearWalk");
        return currentJoin;
    }

    @property bool empty()
    {
        return isEmpty;
    }

    private void lastEdgeOfCycle()
    {
        isCyclic = Yes.isCyclic;
        // Find the missing edge.
        currentJoin = scaffold
            .incidentEdges(currentNode)
            .filter!(join => join != currentJoin)
            .front;
    }

    private void endOfWalk()
    {
        currentNodeIdx = scaffold.nodes.length;
        isEmpty = true;
    }

    private @property ContigNode currentNode()
    {
        return scaffold.nodes[currentNodeIdx];
    }

    private @property void currentNode(ContigNode node)
    {
        this.currentNodeIdx = this.scaffold.indexOf(node);
    }

    private void markVisited(size_t nodeIdx)
    {
        this.visitedNodes.add(nodeIdx);
    }
}

/// Get a range of `ContigNode`s where full contig walks should start.
auto scaffoldStarts(T)(
    Scaffold!T scaffold,
    Scaffold!T.IncidentEdgesCache incidentEdgesCache = Scaffold!T.IncidentEdgesCache.init,
)
{
    static struct ContigStarts
    {
        alias IncidentEdgesCache = Scaffold!T.IncidentEdgesCache;
        static enum emptyIncidentEdgesCache = IncidentEdgesCache.init;

        Scaffold!T scaffold;
        IncidentEdgesCache incidentEdgesCache;
        bool _empty = false;
        NaturalNumberSet unvisitedNodes;
        ContigNode currentContigStart;

        this(Scaffold!T scaffold, IncidentEdgesCache incidentEdgesCache = emptyIncidentEdgesCache)
        {
            this.scaffold = scaffold;
            this.incidentEdgesCache = incidentEdgesCache == emptyIncidentEdgesCache
                ? scaffold.allIncidentEdges()
                : incidentEdgesCache;
            unvisitedNodes.reserveFor(scaffold.nodes.length);

            foreach (nodeIdx; iota(scaffold.nodes.length))
            {
                unvisitedNodes.add(nodeIdx);
            }

            popFront();
        }

        void popFront()
        {
            ContigNode walkToEndNode(ContigNode lastNode, Join!T walkEdge)
            {
                auto nextNode = walkEdge.target(lastNode);
                unvisitedNodes.remove(scaffold.indexOf(nextNode));

                return nextNode;
            }

            assert(!empty, "Attempting to popFront an empty ContigStarts");

            if (unvisitedNodes.empty)
            {
                _empty = true;

                return;
            }

            auto unvisitedNodeIdx = unvisitedNodes.minElement();
            auto unvisitedNode = scaffold.nodes[unvisitedNodeIdx];
            auto unvisitedNodeDegree = incidentEdgesCache[unvisitedNode].length;
            unvisitedNodes.remove(unvisitedNodeIdx);

            // Ignore unconnected nodes.
            if (unvisitedNodeDegree > 0)
            {
                auto endNodes = incidentEdgesCache[unvisitedNode]
                    .map!(firstEdge => linearWalk!T(scaffold, unvisitedNode, firstEdge, incidentEdgesCache)
                        .fold!walkToEndNode(unvisitedNode));
                currentContigStart = unvisitedNodeDegree == 1
                    // If the start node has only one edge it is itself an end node.
                    ? minElement(endNodes, unvisitedNode)
                    : minElement(endNodes);
            }
            else
            {
                popFront();
            }
        }

        @property ContigNode front() pure nothrow
        {
            assert(!empty, "Attempting to fetch the front of an empty ContigStarts");
            return currentContigStart;
        }

        @property bool empty() const pure nothrow
        {
            return _empty;
        }
    }

    return ContigStarts(scaffold, incidentEdgesCache);
}

///
unittest
{
    alias Payload = int;
    alias J = Join!Payload;
    alias S = Scaffold!Payload;
    alias CN = ContigNode;
    alias CP = ContigPart;
    //   contig 1      contig 2
    //
    //  o        o     o        o
    //                         / e4
    //    o -- o ------- o -- o
    //              e3
    //
    //    o -- o         o -- o         o -- o
    //                         \ e7 e8 /      \ e9
    //  o        o     o        o     o        o
    //
    //   contig 5       contig 4      contig 3
    //
    auto scaffold1 = buildScaffold!(sumPayloads!Payload, Payload)(5, [
        J(CN(1, CP.end), CN(2, CP.begin), 1), // e3
        J(CN(2, CP.end), CN(2, CP.post ), 1), // e4
        J(CN(4, CP.end), CN(4, CP.post ), 1), // e7
        J(CN(3, CP.pre), CN(3, CP.begin), 1), // e8
        J(CN(3, CP.end), CN(3, CP.post ), 1), // e9
    ]);

    assert(equal(scaffoldStarts!Payload(scaffold1), [
        CN(1, CP.begin),
        CN(3, CP.pre),
        CN(4, CP.begin),
        CN(5, CP.begin),
    ]));

    //   contig 1      contig 2
    //
    //  o        o     o        o
    //
    //              e1
    //    o -- o ------- o -- o
    //     \_________________/
    //              e2
    //
    auto scaffold2 = buildScaffold!(sumPayloads!Payload, Payload)(2, [
        J(CN(1, CP.end), CN(2, CP.begin), 1), // e1
        J(CN(2, CP.end), CN(1, CP.begin ), 1), // e2
    ]);

    assert(equal(scaffoldStarts!Payload(scaffold2), [
        CN(1, CP.begin),
    ]));
}

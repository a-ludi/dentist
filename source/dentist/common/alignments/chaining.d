/**
    Implementation of local alignment chaining. The main function is
    `chainLocalAlignments`.

    The chaining algorithm works by solving a shortest path problem
    on a directed, node and edge weighted graph where each node is a local
    alignment and a directed edge between to nodes exists if they
    `areChainable`. The nodes give a bonus score relative to the amount of
    sequence covered by the local alignment whereas the edges give a penalty
    proportional to the "gap" between the two involved local alignments.

    This graph problem is reduced to a classical edge-weight single source
    shortest paths problem with one additional node `s` that is the source and
    is connected to all other nodes in the graph. The edges from `s` to `x`
    have a weight of `-alignmentScore(x)` and all other edges `(x, y)` have a
    weight of `chainScore(x, y)` which include the `-alignmentScore(y)` term
    that accounts for `y`'s node weight.

    `ChainingOptions` provides fine-grained control over the graph's structure
    by means of `ChainingOptions.maxIndelBps`, `ChainingOptions.maxChainGapBps`
    and `ChainingOptions.maxRelativeOverlap`.

    Selection of the final chains is controlled by `ChainingOptions.minScore`
    and `ChainingOptions.minRelativeScore`. See command line options for
    default values of these.

    See_also: `chainLocalAlignments`, `ChainingOptions`, `areChainable`,
        `alignmentScore`, `chainScore`

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.alignments.chaining;

import dentist.common.alignments.base;
import dentist.util.algorithm : cmpLexicographically;
import dentist.util.graphalgo :
    dagSingleSourceShortestPaths,
    connectedComponents;
import dentist.util.log;
import dentist.util.math :
    absdiff,
    NaturalNumberSet;
import std.algorithm :
    all,
    cache,
    chunkBy,
    copy,
    filter,
    joiner,
    map,
    max,
    maxIndex,
    min,
    minIndex,
    reverse,
    sum,
    sort;
import std.array : array;
import std.conv : to;
import std.exception : enforce;
import std.math : abs;
import std.range :
    enumerate,
    iota,
    only,
    tee;
import std.range.primitives;
import std.typecons : tuple, Tuple;
import vibe.data.json : Json, toJson = serializeToJson;

debug import std.stdio : writefln, writeln;


/// Options for that control the chaining algorithm.
struct ChainingOptions
{
    /// Maximum absolute distance between neighboring ends of local
    /// alignments.
    ///
    /// See_also: `indel`
    coord_t maxIndelBps;

    /// Maximum postively truncated distance between neighboring ends of local
    /// alignments.
    ///
    /// See_also: `gap`
    coord_t maxChainGapBps;

    /// Maximum absolute distance between neighboring ends of local
    /// alignments.
    ///
    /// See_also: `indel`
    double maxRelativeOverlap;

    /// Minimum chain score as fraction of the best score.
    ///
    /// See_also: `effectiveMinScore`, `chainScore`
    double minRelativeScore;

    /// Minimum chain score.
    ///
    /// See_also: `effectiveMinScore`, `chainScore`
    arithmetic_t minScore;


    /// Returns the effective minimum score taking into account
    /// `minRelativeScore` relative to `bestScore` and `minScore`.
    arithmetic_t effectiveMinScore(arithmetic_t bestScore) const pure @safe
    {
        return max(
            minScore,
            minRelativeScore * bestScore,
        ).to!arithmetic_t;
    }
}


/// Chain local alignments contained in `inputAlignments`.
auto chainLocalAlignments(R)(R inputAlignments, const ChainingOptions options)
    if (isInputRange!R && is(ElementType!R == FlatLocalAlignment))
{
    FlatLocalAlignment lastLA;
    lastLA.id = id_t.max;

    return inputAlignments
        // disregard disabled alignments
        .filter!"!a.flags.disabled"
        // make sure local alignments are ordered...
        .map!((la) {
            enforce(
                lastLA.id == id_t.max || cmpIds(lastLA, la) <= 0,
                "local alignments are not ordered properly",
            );

            lastLA = la;

            return la;
        })
        // process in blocks of alignments that share contig A and B IDs
        //
        // Note: this could be run in parallel as each call is independent.
        .chunkBy!sameIds
        .map!(chunk => buildAlignmentChains(chunk.array, options))
        .cache
        .joiner;
}


private auto buildAlignmentChains(
    FlatLocalAlignment[] flatLocalAlignments,
    const ChainingOptions options,
)
in (flatLocalAlignments.length > 0 && flatLocalAlignments.all!(la => sameIds(la, flatLocalAlignments[0])))
{
    // These functions define the directed and weighted graph with unsigned
    // integer nodes `x` and `y`. The graph defined by these functions does
    // NOT have the additional source node. This is added later when the
    // actual shortest path problem is solved
    alias _areChainable = (x, y) => areChainable(
        flatLocalAlignments[x],
        flatLocalAlignments[y],
        options,
    );
    alias _areChainableUndirected = (x, y) => _areChainable(x, y) || _areChainable(y, x);
    alias _alignmentScore = (x) => alignmentScore(
        flatLocalAlignments[x],
        options,
    );
    alias _chainScore = (x, y) => chainScore(
        flatLocalAlignments[x],
        flatLocalAlignments[y],
        options,
    );

    // split into independent, ie. non-chainable, components to:
    // 1. potentially produce all non-overlapping chains
    // 2. reduce runtime of `dagSingleSourceShortestPaths` which has O(n^^3)
    //    where n is the number of local alignments
    auto components = connectedComponents!_areChainableUndirected(flatLocalAlignments.length);

    debug (chaining)
    {
        import std.algorithm :
            find,
            canFind;

        enum x = 6491;
        enum y = 13019;

        if (x < flatLocalAlignments.length && y < flatLocalAlignments.length)
        {
            assert(_areChainable(x, y));
            assert(_areChainableUndirected(x, y));

            auto xComponent = components.find!(c => c.canFind(x));
            auto yComponent = components.find!(c => c.canFind(y));

            assert(xComponent == yComponent);
        }
    }

    logJsonDebug(
        "info", "buildAlignmentChains",
        "contigA", flatLocalAlignments[0].contigA.contig.toJson,
        "contigB", flatLocalAlignments[0].contigB.contig.toJson,
        "numComponents", components.length,
        "componentSizes", components
            .map!"a.length"
            .array
            .toJson,
    );

    auto selectedChains = components
        .enumerate
        // find best chain for each component
        .map!((enumComponent) {
            auto componentIdx = enumComponent.index;
            auto component = enumComponent.value;
            // Find best chain by means of a shortest paths problem. The
            // lambda functions introduce node 0 as the source node shifting
            // all other indices by +1 and implement the weighting scheme.
            // See `chainLocalAlignments` for a description of the full graph.
            auto ratedChains = dagSingleSourceShortestPaths!(
                (x, y) => (x == 0 && y > 0) ||
                          (x > 0 && y > 0 && _areChainable(component[x - 1], component[y - 1])),
                (x, y) => x == 0 && y > 0
                    ? -_alignmentScore(component[y - 1])
                    : _chainScore(component[x - 1], component[y - 1]),
            )(0, component.length + 1);

            // Extract and sort all computed distances. The list is sorted
            // from best to worst chain.
            auto sortedDistances = ratedChains
                .distances
                .enumerate
                .map!(x => cast(Tuple!(size_t, "endNode", arithmetic_t, "distance")) x)
                .array
                .sort!((x, y) => x.distance < y.distance);
            // Computed score threshold in terms of distance
            auto maxDistance = -options.effectiveMinScore(-sortedDistances[0].distance);

            size_t[] selectedEndNodes;
            selectedEndNodes.reserve(component.length);
            auto forbiddenNodes = NaturalNumberSet(component.length + 1);
            forbiddenNodes.add(0);
            bool[size_t] alternateEndNodes;

            // Find best (alternate) chain paths. An alternate chain path is
            // one that shares a prefix with another (better) chain path but
            // ends in a separate node.
            foreach (endNode, distance; sortedDistances)
                if (endNode !in forbiddenNodes && distance <= maxDistance)
                {
                    foreach (pathNode; ratedChains.reverseShortestPath(endNode))
                    {
                        if (pathNode > 0)
                        {
                            if (pathNode in forbiddenNodes)
                                alternateEndNodes[endNode] = true;
                            forbiddenNodes.add(pathNode);
                        }
                    }

                    selectedEndNodes ~= endNode;
                }

            // Construct chains alongside their scores from the end nodes.
            auto selectedChains = selectedEndNodes
                .map!((endNode) {
                    auto chainPath = ratedChains
                        .reverseShortestPath(endNode)
                        .array
                        .reverse[1 .. $];
                    auto additionalFlags = endNode in alternateEndNodes
                        ? Flags(Flag.alternateChain)
                        : Flags();

                    return tuple!("alignmentChain", "score")(
                        chainPath
                            .map!(node => flatLocalAlignments[component[node - 1]])
                            .composeAlignmentChain(additionalFlags),
                        -ratedChains.distance(endNode),
                    );
                })
                .array;

            logJsonDebug(
                "info", "selectedChains",
                "contigA", flatLocalAlignments[0].contigA.contig.toJson,
                "contigB", flatLocalAlignments[0].contigB.contig.toJson,
                "componentIdx", componentIdx,
                "numSelectedChains", selectedEndNodes.length,
                "bestScore", -sortedDistances[0].distance,
                "minScore", -maxDistance,
                "allScores", sortedDistances.toJson,
            );

            return selectedChains;
        })
        .cache
        .joiner
        .array;

    // determine minimum chain score
    auto bestChainIdx = selectedChains.maxIndex!"a.score < b.score";
    auto bestChain = selectedChains[bestChainIdx];
    auto minScore = options.effectiveMinScore(bestChain.score);

    // filter chains by minimum chain score
    auto acceptedChains = new AlignmentChain[selectedChains.length];
    auto bufferRest = selectedChains
        .filter!(acWithScore => minScore <= acWithScore.score)
        .map!(acWithScore => acWithScore.alignmentChain)
        .copy(acceptedChains);
    acceptedChains.length -= bufferRest.length;
    // make sure chains are sorted
    acceptedChains.sort();

    logJsonDebug(
        "info", "buildAlignmentChains",
        "contigA", flatLocalAlignments[0].contigA.contig.toJson,
        "contigB", flatLocalAlignments[0].contigB.contig.toJson,
        "numComponents", components.length,
        "numAcceptedChains", acceptedChains.length,
        "bestChainIdx", bestChainIdx,
        "bestChainScore", bestChain.score,
        "minScore", minScore,
    );

    return acceptedChains;
}


private AlignmentChain composeAlignmentChain(R)(R chainedLocalAlignments, Flags additionalFlags = Flags())
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;

    // compose local alignments into an alignment chain
    auto firstLA = chainedLocalAlignments.front;
    trace_point_t tracePointDistance = firstLA.tracePoints.length > 0
        ? firstLA.tracePointDistance
        : 0;
    auto alignmentChain = AlignmentChain(
        firstLA.id,
        firstLA.contigA.contig,
        firstLA.contigB.contig,
        firstLA.flags | additionalFlags,
        chainedLocalAlignments
            .map!(fla => LocalAlignment(
                fla.contigA.locus,
                fla.contigB.locus,
                fla.tracePoints.map!"a.numDiffs".sum,
                fla.tracePoints,
            ))
            .array,
        tracePointDistance,
    );

    return alignmentChain;
}


/// Returns the size of the gap between `x` and `y` on contig `seq`. This is
/// negative if they overlap. `x` and `y` are expected to be `areChainable`.
arithmetic_t gap(char seq)(const FlatLocalAlignment x, const FlatLocalAlignment y) pure nothrow @safe
{
    return mixin("cast(arithmetic_t) y.contig" ~ seq ~ ".begin - cast(arithmetic_t) x.contig" ~ seq ~ ".end");
}


/// Returns the maximum of the absolute gap sizes on contig A and B. `x` and
/// `y` are expected to be `areChainable`.
arithmetic_t maxAbsGap(const FlatLocalAlignment x, const FlatLocalAlignment y) pure nothrow @safe
{
    return max(abs(gap!'A'(x, y)), abs(gap!'B'(x, y)));
}


/// Returns the size of the overlap between `x` and `y` on contig `seq`,
/// i.e. the number of overlapping bases if they overlap and 0 otherwise.
arithmetic_t overlap(char seq)(const FlatLocalAlignment x, const FlatLocalAlignment y) pure nothrow @safe
{
    return max(0, -gap!seq(x, y));
}


/// Returns the number of bases covered by `la` on contig `seq`.
arithmetic_t length(char seq)(const FlatLocalAlignment la) pure nothrow @safe
{
    return mixin("cast(arithmetic_t) la.contig" ~ seq ~ ".end - la.contig" ~ seq ~ ".begin");
}

/// Returns the minimum of `length!'A'` and `length!'B'`.
arithmetic_t minLength(char seq)(const FlatLocalAlignment x, const FlatLocalAlignment y) pure nothrow @safe
{
    return min(length!seq(x), length!seq(y));
}


/// Returns the absolute difference between gap sizes on contig A and B.
arithmetic_t indel(const FlatLocalAlignment x, const FlatLocalAlignment y) pure nothrow @safe
{
    return absdiff(gap!'A'(x, y), gap!'B'(x, y));
}


/// Returns true if `x` and `y` have the same contig A and B IDs.
bool sameIds(
    const FlatLocalAlignment x,
    const FlatLocalAlignment y,
) pure nothrow @safe
{
    return x.contigA.id == y.contigA.id &&
           x.contigB.id == y.contigB.id;
}


/// Compares `x` and `y` by contig A and B ID in this order.
int cmpIds(const FlatLocalAlignment lhs, const FlatLocalAlignment rhs) pure nothrow @safe
{
    return cmpLexicographically!(
        const(FlatLocalAlignment),
        la => la.contigA.id,
        la => la.contigB.id,
    )(lhs, rhs);
}


// Returns true iff `x` may precede `y` in a chain.
// This relation is asymmetric and irreflexive but not transitive.
bool areChainable(
    const FlatLocalAlignment x,
    const FlatLocalAlignment y,
    const ChainingOptions options,
) pure nothrow @safe
    in (x.contigA.id == y.contigA.id && x.contigB.id == y.contigB.id, "local alignments must be grouped by ids")
    out (chainable; !chainable || !areChainable(y, x, options), "should be asymmetric")
    out (chainable; x.id != y.id || !chainable, "should be irreflexive")
{
    if (x.flags.complement != y.flags.complement)
        return false;

    return x.contigA.begin < y.contigA.begin && x.contigB.begin < y.contigB.begin && // chain only if both seqs advance
           indel(x, y) <= options.maxIndelBps &&
           maxAbsGap(x, y) <= options.maxChainGapBps &&
           (overlap!'A'(x, y) <= options.maxRelativeOverlap * minLength!'A'(x, y)) &&
           (overlap!'B'(x, y) <= options.maxRelativeOverlap * minLength!'B'(x, y));
}


/// Return the mean `length` of `x` on contig A and B.
arithmetic_t alignmentScore(
    ref const FlatLocalAlignment x,
    const ChainingOptions options,
) pure nothrow @safe
{
    return (length!'A'(x) + length!'B'(x))/2;
}


/// Returns `indel(x, y) + maxAbsGap(x, y)/10 - alignmentScore(y)`.
///
/// See_also: `indel`, `maxAbsGap`, `alignmentScore`
arithmetic_t chainScore(
    ref const FlatLocalAlignment x,
    ref const FlatLocalAlignment y,
    const ChainingOptions options,
) pure nothrow @safe
    in (areChainable(x, y, options))
{
    return indel(x, y) + maxAbsGap(x, y)/10 - alignmentScore(y, options);
}

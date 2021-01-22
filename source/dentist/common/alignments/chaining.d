/**
    Implementation of local alignment chaining.

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

struct ChainingOptions
{
    coord_t maxIndelBps;
    coord_t maxChainGapBps;
    double maxRelativeOverlap;
    double minRelativeScore;
    arithmetic_t minScore;


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
        .filter!"!a.flags.disabled"
        .map!((la) {
            // make sure local alignments are ordered...
            enforce(
                lastLA.id == id_t.max || cmpIds(lastLA, la) <= 0,
                "local alignments are not ordered properly",
            );

            lastLA = la;

            return la;
        })
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

    // find best chain for each component
    auto selectedChains = components
        .enumerate
        .map!((enumComponent) {
            auto componentIdx = enumComponent.index;
            auto component = enumComponent.value;
            // find best chain by means of a shortest paths problem
            auto ratedChains = dagSingleSourceShortestPaths!(
                (x, y) => (x == 0 && y > 0) ||
                          (x > 0 && y > 0 && _areChainable(component[x - 1], component[y - 1])),
                (x, y) => x == 0 && y > 0
                    ? -_alignmentScore(component[y - 1])
                    : _chainScore(component[x - 1], component[y - 1]),
            )(0, component.length + 1);

            auto sortedDistances = ratedChains
                .distances
                .enumerate
                .map!(x => cast(Tuple!(size_t, "endNode", arithmetic_t, "distance")) x)
                .array
                .sort!((x, y) => x.distance < y.distance);
            auto maxDistance = -options.effectiveMinScore(-sortedDistances[0].distance);

            size_t[] selectedEndNodes;
            selectedEndNodes.reserve(component.length);
            auto forbiddenNodes = NaturalNumberSet(component.length + 1);
            forbiddenNodes.add(0);
            bool[size_t] alternateEndNodes;

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


AlignmentChain composeAlignmentChain(R)(R chainedLocalAlignments, Flags additionalFlags = Flags())
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;

    // compose local alignments into an alignment chain
    auto firstLA = chainedLocalAlignments.front;
    auto tracePointDistance = firstLA.tracePointDistance;
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


arithmetic_t gap(char seq)(const FlatLocalAlignment x, const FlatLocalAlignment y) pure nothrow @safe
{
    return mixin("cast(arithmetic_t) y.contig" ~ seq ~ ".begin - cast(arithmetic_t) x.contig" ~ seq ~ ".end");
}


arithmetic_t maxAbsGap(const FlatLocalAlignment x, const FlatLocalAlignment y) pure nothrow @safe
{
    return max(abs(gap!'A'(x, y)), abs(gap!'B'(x, y)));
}


arithmetic_t overlap(char seq)(const FlatLocalAlignment x, const FlatLocalAlignment y) pure nothrow @safe
{
    return max(0, -gap!seq(x, y));
}


arithmetic_t length(char seq)(const FlatLocalAlignment la) pure nothrow @safe
{
    return mixin("cast(arithmetic_t) la.contig" ~ seq ~ ".end - la.contig" ~ seq ~ ".begin");
}


arithmetic_t minLength(char seq)(const FlatLocalAlignment x, const FlatLocalAlignment y) pure nothrow @safe
{
    return min(length!seq(x), length!seq(y));
}


arithmetic_t indel(const FlatLocalAlignment x, const FlatLocalAlignment y) pure nothrow @safe
{
    return absdiff(gap!'A'(x, y), gap!'B'(x, y));
}


bool sameIds(
    const FlatLocalAlignment x,
    const FlatLocalAlignment y,
) pure nothrow @safe
{
    return x.contigA.id == y.contigA.id &&
           x.contigB.id == y.contigB.id;
}


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


arithmetic_t alignmentScore(
    ref const FlatLocalAlignment x,
    const ChainingOptions options,
) pure nothrow @safe
{
    return (length!'A'(x) + length!'B'(x))/2;
}


arithmetic_t chainScore(
    ref const FlatLocalAlignment x,
    ref const FlatLocalAlignment y,
    const ChainingOptions options,
) pure nothrow @safe
    in (areChainable(x, y, options))
{
    return indel(x, y) + maxAbsGap(x, y)/10 - alignmentScore(y, options);
}

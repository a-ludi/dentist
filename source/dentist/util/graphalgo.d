/**
    Implementation of the Floyd–Warshall algorithm.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.graphalgo;

import dentist.util.math :
    absdiff,
    NaturalNumberSet;
import dentist.util.saturationmath;
import std.algorithm :
    any,
    copy,
    countUntil,
    map;
import std.array :
    appender,
    array,
    uninitializedArray;
import std.functional : binaryFun;
import std.range :
    enumerate,
    iota;
import std.typecons :
    Tuple,
    Yes;


/**
    Calculate connected components of the graph defined by `hasEdge`.

    Params:
        hasEdge = Binary predicate taking two nodes of type `size_t` which is
                  true iff the first node is adjacent to the second node.
        n =       Number of nodes in the graph. `hasEdge` must be defined for
                  every pair of integer in `0 .. n`.
    Returns: Array of components represented as arrays of node indices.
*/
size_t[][] connectedComponents(alias hasEdge)(size_t n)
{
    alias _hasEdge = binaryFun!hasEdge;

    auto unvisitedNodes = NaturalNumberSet(n, Yes.addAll);
    auto nodesBuffer = new size_t[n];
    auto components = appender!(size_t[][]);

    while (!unvisitedNodes.empty)
    {
        // discover startNode's component by depth-first search
        auto component = discoverComponent!_hasEdge(unvisitedNodes);
        // copy component indices to buffer
        auto restNodesBuffer = component
            .elements
            .copy(nodesBuffer);
        // append component to result
        components ~= nodesBuffer[0 .. $ - restNodesBuffer.length];
        // reduce node buffer
        nodesBuffer = restNodesBuffer;

    }

    return components.data;
}

///
unittest
{
    import std.algorithm :
        equal,
        min;

    //    _____________
    //   /             \
    // (0) --- (1) --- (2)     (3) --- (4)
    enum n = 5;
    alias connect = (u, v, x, y) => (u == x && v == y) || (u == y && v == x);
    alias hasEdge = (u, v) => connect(u, v, 0, 1) ||
                              connect(u, v, 1, 2) ||
                              connect(u, v, 2, 0) ||
                              connect(u, v, 3, 4);

    auto components = connectedComponents!hasEdge(n);

    assert(equal(components, [
        [0, 1, 2],
        [3, 4],
    ]));
}

unittest
{
    import std.algorithm :
        equal,
        min;

    //    _____________
    //   /             \
    // (0) --- (1) --- (2)     (3) --- (4)
    //   \_____________________/
    enum n = 5;
    alias connect = (u, v, x, y) => (u == x && v == y) || (u == y && v == x);
    alias hasEdge = (u, v) => connect(u, v, 0, 1) ||
                              connect(u, v, 1, 2) ||
                              connect(u, v, 2, 0) ||
                              connect(u, v, 0, 3) ||
                              connect(u, v, 3, 4);

    auto components = connectedComponents!hasEdge(n);

    import std.stdio;
    assert(equal(components, [
        [0, 1, 2, 3, 4],
    ]));
}


private NaturalNumberSet discoverComponent(alias hasEdge)(ref NaturalNumberSet nodes)
{
    assert(!nodes.empty, "cannot discoverComponent of an empty graph");

    // prepare component
    auto component = NaturalNumberSet(nodes.maxElement);
    // select start node
    auto currentNode = nodes.minElement;

    discoverComponent!hasEdge(nodes, currentNode, component);

    return component;
}


private void discoverComponent(alias hasEdge)(ref NaturalNumberSet nodes, size_t currentNode, ref NaturalNumberSet component)
{
    // move currentNode from available nodes to the component
    component.add(currentNode);
    nodes.remove(currentNode);

    // try to find successor of current node
    foreach (nextNode; nodes.elements)
    {
        if (hasEdge(currentNode, nextNode))
        {
            assert(
                hasEdge(nextNode, currentNode),
                "connectedComponents may be called only on an undirected graph",
            );
            // found successor -> recurse

            discoverComponent!hasEdge(nodes, nextNode, component);
        }
    }
}



struct FloydWarshallMatrix(weight_t)
{
    enum unconnectedWeight = saturatedInfinity!weight_t;
    enum noNext = size_t.max;

    private size_t _numNodes;
    private weight_t[] _dist;
    private size_t[] _next;


    @property size_t numNodes() const pure nothrow @safe
    {
        return _numNodes;
    }


    bool hasNegativeCycles() const pure nothrow @safe
    {
        return iota(_numNodes).any!(n => dist(n, n) < 0);
    }


    private size_t idx(size_t u, size_t v) const pure nothrow @safe
    in (u < numNodes && v < numNodes, "index out of bounds")
    {
        return u * numNodes + v;
    }


    @property ref weight_t dist(size_t u, size_t v) pure nothrow @safe
    {
        return _dist[idx(u, v)];
    }


    @property weight_t dist(size_t u, size_t v) const pure nothrow @safe
    {
        return _dist[idx(u, v)];
    }


    @property bool isConnected(size_t u, size_t v) const pure nothrow @safe
    {
        return dist(u, v) < unconnectedWeight;
    }


    @property ref size_t next(size_t u, size_t v) pure nothrow @safe
    {
        return _next[idx(u, v)];
    }


    @property size_t next(size_t u, size_t v) const pure nothrow @safe
    {
        return _next[idx(u, v)];
    }


    @property bool hasNext(size_t u, size_t v) const pure nothrow @safe
    {
        return next(u, v) != noNext;
    }


    static struct ShortestPath
    {
        private const(FloydWarshallMatrix!weight_t)* _matrix;
        private size_t _from;
        private size_t _to;
        private size_t _current;


        private this(const(FloydWarshallMatrix!weight_t)* matrix, size_t from, size_t to)
        {
            this._matrix = matrix;
            this._from = from;
            this._to = to;
            this._current = matrix !is null && matrix.isConnected(from, to)
                ? from
                : noNext;
        }


        @property const(FloydWarshallMatrix!weight_t) matrix() pure nothrow @safe
        {
            return *_matrix;
        }


        @property size_t from() const pure nothrow @safe
        {
            return _from;
        }


        @property size_t to() const pure nothrow @safe
        {
            return _to;
        }


        @property bool empty() const pure nothrow @safe
        {
            return _current == noNext;
        }


        @property size_t front() const pure nothrow @safe
        {
            assert(
                !empty,
                "Attempting to fetch the front of an empty FloydWarshallMatrix.ShortestPath",
            );

            return _current;
        }


        void popFront() pure nothrow @safe
        {
            assert(!empty, "Attempting to popFront an empty FloydWarshallMatrix.ShortestPath");

            this._current = _matrix !is null
                ? matrix.next(_current, _to)
                : noNext;
        }


        @property ShortestPath save() const pure nothrow @safe
        {
            return cast(typeof(return)) this;
        }
    }


    ShortestPath shortestPath(size_t from, size_t to) const pure nothrow
    {
        return ShortestPath(&this, from, to);
    }
}


private void updateBestConnections(weight_t)(
    size_t u,
    size_t v,
    weight_t d,
    ref size_t[2][] bestConnections,
    ref weight_t[] bestDists,
) nothrow @safe
{
    foreach (i, ref bestDist; bestDists)
    {
        if (d < bestDist)
        {
            // shift other best weights & connections one place further down
            foreach_reverse (j; i + 1 .. bestDists.length)
            {
                bestDists[j] = bestDists[j - 1];
                bestConnections[j] = bestConnections[j - 1];
            }
            // update current best weight & connection
            bestDist = d;
            bestConnections[i] = [u, v];

            return;
        }
    }
}


private auto floydWarshallMatrix(
    alias hasEdge,
    alias weight,
    weight_t = typeof(weight(size_t.init, size_t.init)),
)(size_t n, ref size_t[2][] bestConnections, ref weight_t[] bestDists)
    in (bestConnections.length <= n^^2 && bestConnections.length == bestDists.length)
{
    FloydWarshallMatrix!weight_t matrix;

    matrix._numNodes = n;
    matrix._dist = uninitializedArray!(weight_t[])(n * n);
    matrix._next = uninitializedArray!(size_t[])(n * n);
    foreach (ref c; bestConnections)
        c = [0, 0];
    bestDists[] = cast(weight_t) 0;

    foreach (u; 0 .. matrix.numNodes)
        foreach (v; 0 .. matrix.numNodes)
        {
            if (u == v)
            {
                matrix.dist(u, v) = 0;
                matrix.next(u, v) = matrix.noNext;
            }
            else if (hasEdge(u, v))
            {
                auto w = weight(u, v);
                matrix.dist(u, v) = w;
                matrix.next(u, v) = v;
                updateBestConnections(u, v, w, bestConnections, bestDists);
            }
            else
            {
                matrix.dist(u, v) = matrix.unconnectedWeight;
                matrix.next(u, v) = matrix.noNext;
            }
        }

    return matrix;
}


enum GraphType : ubyte
{
    /// Any weighted graph (directed or undirected). Negative cycles are
    /// allowed and may be detected after the algorithm finished.
    general,
    /// Directed acyclic graph. The algorithm can be faster by a constant
    /// factor by first sorting in topoloigcal order and than skipping
    /// irrelevant edges.
    DAG,
}


/**
    Calculate all shortest paths between all pairs of nodes. The functions
    `hasEdge` and `weight` define the graphs structure and weights,
    respectively. Nodes are represented as `size_t` integers.

    Params:
        hasEdge = Binary predicate taking two nodes of type `size_t` which is
                  true iff the first node is adjacent to the second node.
        weight =  Binary function taking two nodes of type `size_t` which
                  returns the weight of the edge between the first and the
                  second node. The function may be undefined if `hasEdge`
                  returns false for the given arguments.
        n =       Number of nodes in the graph. `hasEdge` must be defined for
                  every pair of integer in `0 .. n`.
        bestConnections =
                  If given, the array will be populated with the
                  `bestConnections.length` best connections, that
                  is the pairs of nodes with optimal distances.
    Returns: FloydWarshallMatrix
*/
auto shortestPathsFloydWarshall(
    alias hasEdge,
    alias weight,
    GraphType graphType = GraphType.general,
)(size_t n)
{
    size_t[2][] bestConnections;

    return shortestPathsFloydWarshall!(hasEdge, weight, graphType)(n, bestConnections);
}

/// ditto
auto shortestPathsFloydWarshall(
    alias hasEdge,
    alias weight,
    GraphType graphType = GraphType.general,
)(
    size_t n,
    ref size_t[2][] bestConnections,
)
{
    alias _hasEdge = binaryFun!hasEdge;
    alias _weight = binaryFun!weight;
    alias weight_t = typeof(_weight(size_t.init, size_t.init));

    static if (graphType == GraphType.DAG)
    {
        auto nodes = topologicalSort!_hasEdge(n);

        alias nodesSlice = (from, to) => nodes[from .. to].enumerate(from);
        alias kRange = () => nodesSlice(1, n - 1);
        alias uRange = (kIdx) => nodesSlice(0, kIdx);
        alias vRange = (kIdx, uIdx) => nodesSlice(kIdx + 1, n);
    }
    else
    {
        auto nodes = enumerate(iota(n));

        alias kRange = () => nodes;
        alias uRange = (kIdx) => nodes;
        alias vRange = (kIdx, uIdx) => nodes;
    }

    auto bestDists = uninitializedArray!(weight_t[])(bestConnections.length);
    auto matrix = floydWarshallMatrix!(_hasEdge, _weight)(
        n,
        bestConnections,
        bestDists,
    );

    if (n < 2)
        return matrix;

    foreach (kIdx, k; kRange())
        foreach (uIdx, u; uRange(kIdx))
            foreach (vIdx, v; vRange(kIdx, uIdx))
            {
                auto ukDist = matrix.dist(u, k);
                auto kvDist = matrix.dist(k, v);
                auto d = saturatedAdd(ukDist, kvDist);

                if (matrix.dist(u, v) > d)
                {
                    matrix.dist(u, v) = d;
                    matrix.next(u, v) = matrix.next(u, k);
                    updateBestConnections(u, v, d, bestConnections, bestDists);
                }
            }

    return matrix;
}

/// ditto
alias allPairsShortestPaths = shortestPathsFloydWarshall;

///
unittest
{
    import std.algorithm : equal;

    //    _____________   _____________
    //   /             v /             v
    // (0) --> (1) --> (2)     (3) --> (4)
    enum n = 5;
    alias hasEdge = (u, v) => (u + 1 == v && u != 2) ||
                              (u + 2 == v && u % 2 == 0);
    alias weight = (u, v) => 1;

    auto shortestPaths = shortestPathsFloydWarshall!(hasEdge, weight)(n);

    debug (2) shortestPaths.printConnections([[0, 4], [1, 4], [1, 3], [3, 4]]);

    assert(equal(shortestPaths.shortestPath(0, 4), [0, 2, 4]));
    assert(shortestPaths.dist(0, 4) == 2);
    assert(equal(shortestPaths.shortestPath(1, 4), [1, 2, 4]));
    assert(shortestPaths.dist(1, 4) == 2);
    assert(equal(shortestPaths.shortestPath(1, 3), size_t[].init));
    assert(!shortestPaths.isConnected(1, 3));
    assert(equal(shortestPaths.shortestPath(3, 4), [3, 4]));
    assert(shortestPaths.dist(3, 4) == 1);
}

///
unittest
{
    import std.algorithm : equal;

    //    _____-4______   _____-4______
    //   /             v /             v
    // (0) --> (1) --> (2)     (3) --> (4)
    //     -1      -1              -1
    enum n = 5;
    alias hasEdge = (u, v) => (u + 1 == v && u != 2) ||
                              (u + 2 == v && u % 2 == 0);
    alias weight = (u, v) => -(cast(long) u - cast(long) v)^^2;

    auto shortestPaths = shortestPathsFloydWarshall!(hasEdge, weight)(n);

    debug (2) shortestPaths.printConnections([[0, 4], [1, 4], [2, 4], [1, 3], [3, 4]]);

    assert(equal(shortestPaths.shortestPath(0, 4), [0, 2, 4]));
    assert(shortestPaths.dist(0, 4) == -8);
    assert(equal(shortestPaths.shortestPath(1, 4), [1, 2, 4]));
    assert(shortestPaths.dist(1, 4) == -5);
    assert(equal(shortestPaths.shortestPath(1, 3), size_t[].init));
    assert(!shortestPaths.isConnected(1, 3));
    assert(equal(shortestPaths.shortestPath(3, 4), [3, 4]));
    assert(shortestPaths.dist(3, 4) == -1);
}

///
unittest
{
    import std.algorithm : equal, map;

    //    _____-4______   _____-4______
    //   /             v /             v
    // (0) --> (1) --> (2)     (3) --> (4)
    //     -1      -1              -1
    enum n = 5;
    alias hasEdge = (u, v) => (u + 1 == v && u != 2) ||
                              (u + 2 == v && u % 2 == 0);
    alias weight = (u, v) => -(cast(long) u - cast(long) v)^^2;

    enum bestN = 7;
    auto bestConnections = new size_t[2][bestN];
    auto shortestPaths = shortestPathsFloydWarshall!(hasEdge, weight)(n, bestConnections);

    assert(equal(bestConnections, [
        [0, 4],
        [1, 4],
        [0, 2],
        [2, 4],
        [0, 1],
        [1, 2],
        [3, 4],
    ]));
    assert(equal(bestConnections.map!(c => shortestPaths.dist(c[0], c[1])), [
        -8,
        -5,
        -4,
        -4,
        -1,
        -1,
        -1,
    ]));
}

/// Optimize performance by choosing appropriate `graphType`
unittest
{
    import std.datetime.stopwatch;

    enum n = 5;
    alias hasEdge = (u, v) => (u + 1 == v && u != 2) ||
                              (u + 2 == v && u % 2 == 0);
    alias weight = (u, v) => -(cast(long) u - cast(long) v)^^2;

    alias shortestPathsGeneral = () => shortestPathsFloydWarshall!(
        hasEdge,
        weight,
        GraphType.general,
    )(n);
    alias shortestPathsDAG = () => shortestPathsFloydWarshall!(
        hasEdge,
        weight,
        GraphType.DAG,
    )(n);

    enum numRounds = 10_000;
    auto result = benchmark!(
        shortestPathsGeneral,  // => 157.029400ms
        shortestPathsDAG,      // =>  35.348500ms
    )(numRounds);

    debug (2)
    {
        import std.stdio : writefln;

        writefln!"Computed %d rounds:"(numRounds);
        writefln!"shortestPathsGeneral:  %fms"(result[0].total!"nsecs"/1e9*1e3);
        writefln!"shortestPathsDAG:      %fms"(result[1].total!"nsecs"/1e9*1e3);
    }
}

// functional test of different graph types
unittest
{
    enum n = 5;
    alias hasEdge = (u, v) => (u + 1 == v && u != 2) ||
                              (u + 2 == v && u % 2 == 0);
    alias weight = (u, v) => -(cast(long) u - cast(long) v)^^2;

    auto shortestPathsGeneral = shortestPathsFloydWarshall!(
        hasEdge,
        weight,
        GraphType.general,
    )(n);
    auto shortestPathsDAG = shortestPathsFloydWarshall!(
        hasEdge,
        weight,
        GraphType.DAG,
    )(n);

    debug (2)
    {
        import std.stdio : writefln;

        size_t[2][] connections = [[0, 4], [1, 4], [2, 4], [1, 3], [3, 4]];

        shortestPathsGeneral.printConnections(connections);
        writefln!"topologicalOrder: %(%d > %)"(topologicalSort!hasEdge(n));
        shortestPathsDAG.printConnections(connections);
    }

    assert(shortestPathsGeneral == shortestPathsDAG);
}

unittest
{
    enum n = 0;
    alias hasEdge = (u, v) => true;
    alias weight = (u, v) => 1;

    cast(void) shortestPathsFloydWarshall!(hasEdge, weight)(n);
}

unittest
{
    import std.algorithm : equal;

    enum n = 1;
    alias hasEdge = (u, v) => false;
    alias weight = (u, v) => 0;

    auto bestConnections = new size_t[2][1];
    auto shortestPaths = shortestPathsFloydWarshall!(hasEdge, weight)(n, bestConnections);

    assert(bestConnections == [[0, 0]]);
    assert(shortestPaths.dist(0, 0) == 0);
    assert(equal(shortestPaths.shortestPath(0, 0), [0]));
}


private version (unittest)
{
    void printConnections(weight_t)(
        const ref FloydWarshallMatrix!weight_t shortestPaths,
        size_t[2][] connections,
        string file = __FILE__,
        size_t line = __LINE__,
    )
    {
        import std.conv;
        import std.stdio;

        writefln!"--- BEGIN connections (%s:%d)"(file, line);
        foreach (uv; connections)
            writefln!"%d -> %d (%s): %-(%d -> %)"(
                uv[0],
                uv[1],
                shortestPaths.isConnected(uv[0], uv[1])
                    ? shortestPaths.dist(uv[0], uv[1]).to!string
                    : "-",
                shortestPaths.shortestPath(uv[0], uv[1]),
            );
        writefln!"--- END connections"();
    }
}


struct SingleSourceShortestPathsSolution(weight_t)
{
    enum unconnectedWeight = saturatedInfinity!weight_t;
    enum noPredecessor = size_t.max;

    size_t startNode;
    size_t[] topologicalOrder;
    private weight_t[] _distance;
    private size_t[] _predecessor;


    @property size_t numNodes() const pure nothrow @safe
    {
        return topologicalOrder.length;
    }


    private size_t originalNode(size_t u) const pure nothrow @safe
    {
        return topologicalOrder[u];
    }


    @property const(weight_t)[] distances() const pure nothrow @safe
    {
        return _distance[];
    }


    @property ref weight_t distance(size_t u) pure nothrow @safe
    {
        return _distance[u];
    }


    @property weight_t distance(size_t u) const pure nothrow @safe
    {
        return _distance[u];
    }


    @property bool isConnected(size_t u) const pure nothrow @safe
    {
        return distance(u) < unconnectedWeight;
    }


    @property ref size_t predecessor(size_t u) pure nothrow @safe
    {
        return _predecessor[u];
    }


    @property size_t predecessor(size_t u) const pure nothrow @safe
    {
        return _predecessor[u];
    }


    @property bool hasPredecessor(size_t u) const pure nothrow @safe
    {
        return predecessor(u) != noPredecessor;
    }


    static struct ReverseShortestPath
    {
        private const(SingleSourceShortestPathsSolution!weight_t)* _solution;
        private size_t _to;
        private size_t _current;


        private this(const(SingleSourceShortestPathsSolution!weight_t)* solution, size_t to)
        {
            this._solution = solution;
            this._to = to;
            this._current = solution !is null && solution.isConnected(to)
                ? to
                : noPredecessor;
        }


        @property const(SingleSourceShortestPathsSolution!weight_t) solution() pure nothrow @safe
        {
            return *_solution;
        }


        @property size_t from() const pure nothrow @safe
        {
            return _solution.startNode;
        }


        @property size_t to() const pure nothrow @safe
        {
            return _to;
        }


        @property bool empty() const pure nothrow @safe
        {
            return _current == noPredecessor;
        }


        @property size_t front() const pure nothrow @safe
        {
            assert(
                !empty,
                "Attempting to fetch the front of an empty SingleSourceShortestPathsSolution.ReverseShortestPath",
            );

            return _current;
        }


        void popFront() pure nothrow @safe
        {
            assert(!empty, "Attempting to popFront an empty SingleSourceShortestPathsSolution.ReverseShortestPath");

            this._current = _solution !is null
                ? solution.predecessor(_current)
                : noPredecessor;
        }
    }


    ReverseShortestPath reverseShortestPath(size_t to) const pure nothrow
    {
        return ReverseShortestPath(&this, to);
    }
}


/**
    Calculate all shortest paths in DAG starting at `start`. The
    functions `hasEdge` and `weight` define the graphs structure and
    weights, respectively. Nodes are represented as `size_t` integers.
    The graph must be directed and acyclic (DAG).

    Params:
        hasEdge = Binary predicate taking two nodes of type `size_t` which is
                  true iff the first node is adjacent to the second node.
        weight =  Binary function taking two nodes of type `size_t` which
                  returns the weight of the edge between the first and the
                  second node. The function may be undefined if `hasEdge`
                  returns false for the given arguments.
        n =       Number of nodes in the graph. `hasEdge` must be defined for
                  every pair of integer in `0 .. n`.
    Returns: SingleSourceShortestPathsSolution
*/
auto dagSingleSourceShortestPaths(alias hasEdge, alias weight)(size_t start, size_t n)
{
    alias _hasEdge = binaryFun!hasEdge;
    alias _weight = binaryFun!weight;
    alias weight_t = typeof(_weight(size_t.init, size_t.init));

    SingleSourceShortestPathsSolution!weight_t result;

    with (result)
    {
        // sort topological
        topologicalOrder = topologicalSort!_hasEdge(n);
        alias N = (u) => originalNode(u);

        _distance = uninitializedArray!(weight_t[])(n);
        _distance[] = saturatedInfinity!weight_t;
        _distance[start] = 0;
        _predecessor = uninitializedArray!(size_t[])(n);
        _predecessor[] = size_t.max;

        foreach (u; topologicalOrder.countUntil(start) .. n)
            foreach (v; u + 1 .. n)
                if (_hasEdge(N(u), N(v)))
                {
                    auto vDistance = distance(N(v));
                    auto uDistance = saturatedAdd(distance(N(u)), _weight(N(u), N(v)));

                    if (vDistance > uDistance)
                    {
                        distance(N(v)) = uDistance;
                        predecessor(N(v)) = N(u);
                    }
                }
    }

    return result;
}

///
unittest
{
    import std.algorithm : equal;

    //    _____________   _____________
    //   /             v /             v
    // (0) --> (1) --> (2)     (3) --> (4)
    enum n = 5;
    alias hasEdge = (u, v) => (u + 1 == v && u != 2) ||
                              (u + 2 == v && u % 2 == 0);
    alias weight = (u, v) => 1;

    auto shortestPaths = dagSingleSourceShortestPaths!(hasEdge, weight)(0, n);

    assert(equal(shortestPaths.reverseShortestPath(4), [4, 2, 0]));
    assert(shortestPaths.distance(4) == 2);
    assert(equal(shortestPaths.reverseShortestPath(2), [2, 0]));
    assert(shortestPaths.distance(2) == 1);
    assert(equal(shortestPaths.reverseShortestPath(1), [1, 0]));
    assert(shortestPaths.distance(1) == 1);
    assert(equal(shortestPaths.reverseShortestPath(3), size_t[].init));
    assert(!shortestPaths.isConnected(3));
}


auto topologicalSort(alias hasEdge)(size_t n)
{
    alias _hasEdge = binaryFun!hasEdge;

    // list that will contain the sorted nodes
    auto sortedNodes = new size_t[n];

    auto sortedNodesHead = sortedNodes[];
    void enqueueNode(size_t node)
    {
        sortedNodesHead[$ - 1] = node;
        --sortedNodesHead.length;
    }

    // keep track which nodes have been visited
    auto unvisitedNodes = NaturalNumberSet(n, Yes.addAll);
    auto temporaryVisitedNodes = NaturalNumberSet(n);

    void visit(size_t node)
    {
        if (node !in unvisitedNodes)
            // already visited
            return;

        if (node in temporaryVisitedNodes)
            // cycle detected
            throw new NoDAG();

        temporaryVisitedNodes.add(node);

        foreach (nextNode; unvisitedNodes.elements)
            if (_hasEdge(node, nextNode))
                visit(nextNode);

        temporaryVisitedNodes.remove(node);
        unvisitedNodes.remove(node);
        enqueueNode(node);
    }

    foreach (node; unvisitedNodes.elements)
        visit(node);

    return sortedNodes;
}

///
unittest
{
    import std.algorithm : equal;

    //    _____________   _____________
    //   /             v /             v
    // (0) --> (1) --> (2)     (3) --> (4)
    enum n = 5;
    alias hasEdge = (u, v) => (u + 1 == v && u != 2) ||
                              (u + 2 == v && u % 2 == 0);

    auto topologicalOrder = topologicalSort!hasEdge(n);

    assert(equal(topologicalOrder, [3, 0, 1, 2, 4]));
}

///
unittest
{
    import std.exception : assertThrown;

    //    _____________   _____________
    //   /             v /             v
    // (0) --> (1) --> (2)     (3) --> (4)
    //   ^_____________________________/
    enum n = 5;
    alias hasEdge = (u, v) => (u + 1 == v && u != 2) ||
                              (u + 2 == v && u % 2 == 0) ||
                              u == 4 && v == 0;

    assertThrown!NoDAG(topologicalSort!hasEdge(n));
}


/// Thrown if a cycle was detected.
class NoDAG : Exception
{
    this(string file = __FILE__, size_t line = __LINE__, Throwable next = null)
    {
        super("not a DAG: graph has cycles", file, line, next);
    }
}

/**
    Some additional mathematical functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.math;

import dentist.util.algorithm : cmpLexicographically, sliceBy;
import std.algorithm :
    all,
    among,
    copy,
    countUntil,
    cumulativeFold,
    filter,
    map,
    max,
    min,
    maxElement,
    sort,
    sum,
    swap,
    uniq;
import std.array : Appender, array;
import std.conv : to;
import std.exception : assertThrown;
import std.format : format;
import std.functional : binaryFun, unaryFun;
import std.range :
    assumeSorted,
    chain,
    ElementType,
    enumerate,
    isForwardRange,
    isRandomAccessRange,
    retro,
    save,
    walkLength;
import std.traits : isCallable, isIntegral, isNumeric;
import std.typecons : Flag, No, Yes;

debug import std.stdio : writeln;

/// Calculate the mean of range.
ElementType!Range mean(Range)(Range values) if (isForwardRange!Range)
{
    auto sum = values.save.sum;
    auto length = values.walkLength.to!(ElementType!Range);

    return sum / length;
}

unittest
{
    {
        auto values = [2, 4, 8];
        assert(values.mean == 4);
    }
    {
        auto values = [1.0, 2.0, 3.0, 4.0];
        assert(values.mean == 2.5);
    }
}

/// Calculate the median of range.
ElementType!Range median(Range)(Range values) if (__traits(compiles, sort(values)))
{
    assert(values.length > 0, "median is undefined for empty set");
    auto middleIdx = values.length / 2;
    auto useAverage = values.length > 1 && values.length % 2 == 0;
    auto sortedValues = values.sort;

    if (useAverage)
        return (sortedValues[middleIdx - 1] + sortedValues[middleIdx]) / 2;
    else
        return values[middleIdx];
}

unittest
{
    {
        auto values = [4, 2, 8];
        assert(values.median == 4);
    }
    {
        auto values = [4, 3, 2, 8];
        assert(values.median == 3);
    }
    {
        auto values = [4, 6, 2, 8];
        assert(values.median == 5);
    }
    {
        auto values = [2, 1, 3, 0, 4, 9, 8, 5, 6, 3, 9];
        assert(values.median == 4);
    }
    {
        auto values = [2.0, 1.0, 4.0, 3.0];
        assert(values.median == 2.5);
    }
    {
        auto values = [2.0, 1.0, 4.0, 3.0, 5.0];
        assert(values.median == 3.0);
    }
}

/// Calculate the Nxx (e.g. N50) of values.
ElementType!Range N(real xx, Range, Num)(Range values, Num totalSize) if (__traits(compiles, sort(values)))
{
    static assert(0 < xx && xx < 100, "N" ~ xx.to!string ~ " is undefined for empty set");
    assert(values.length > 0, "N" ~ xx.to!string ~ " is undefined for empty set");
    auto xxPercentile = xx/100.0 * totalSize;
    auto sortedValues = values.sort;
    auto targetIndex = sortedValues
        .retro
        .cumulativeFold!"a + b"(cast(ElementType!Range) 0)
        .countUntil!"a >= b"(xxPercentile);

    if (targetIndex == values.length)
        return 0;
    else
        return sortedValues[$ - targetIndex - 1];
}

unittest
{
    {
        auto totalSize = 54;
        auto values = [2, 3, 4, 5, 6, 7, 8, 9, 10];
        enum N50 = 8;
        enum N10 = 10;

        assert(N!50(values, totalSize) == N50);
        assert(N!10(values, totalSize) == N10);
    }
    {
        auto totalSize = 32;
        auto values = [2, 2, 2, 3, 3, 4, 8, 8];
        enum N50 = 8;

        assert(N!50(values, totalSize) == N50);
    }
}

enum RoundingMode : byte
{
    floor,
    round,
    ceil,
}

/**
    Round x upward according to base, ie. returns the next integer larger or
    equal to x which is divisible by base.

    Returns: x rounded upward according to base.
*/
Integer ceil(Integer)(in Integer x, in Integer base) pure nothrow
        if (isIntegral!Integer)
{
    return x % base == 0
        ? x
        : (x / base + 1) * base;
}

///
unittest
{
    assert(ceil(8, 10) == 10);
    assert(ceil(32, 16) == 32);
    assert(ceil(101, 100) == 200);
}

/**
    Round x downward according to base, ie. returns the next integer smaller or
    equal to x which is divisible by base.

    Returns: x rounded downward according to base.
*/
Integer floor(Integer)(in Integer x, in Integer base) pure nothrow
        if (isIntegral!Integer)
{
    return (x / base) * base;
}

///
unittest
{
    assert(floor(8, 10) == 0);
    assert(floor(32, 16) == 32);
    assert(floor(101, 100) == 100);
}

/// Returns the absolute difference between two numbers.
Num absdiff(Num)(in Num a, in Num b) pure nothrow if (isNumeric!Num)
{
    return a > b
        ? a - b
        : b - a;
}

///
unittest
{
    assert(absdiff(2UL, 3UL) == 1UL);
    assert(absdiff(-42, 13) == 55);
    assert(absdiff(2.5, 5) == 2.5);
}

/// Returns the result of `ceil(a / b)` but uses integer arithmetic only.
Integer ceildiv(Integer)(in Integer a, in Integer b) pure nothrow if (isIntegral!Integer)
{
    Integer resultSign = (a < 0) ^ (b < 0) ? -1 : 1;

    return resultSign < 0 || a % b == 0
        ? a / b
        : a / b + resultSign;
}

///
unittest
{
    assert(ceildiv(0, 3) == 0);
    assert(ceildiv(1UL, 3UL) == 1UL);
    assert(ceildiv(2L, 3L) == 1L);
    assert(ceildiv(3U, 3U) == 1U);
    assert(ceildiv(4, 3) == 2);
    assert(ceildiv(-4, 4) == -1);
    assert(ceildiv(-4, 3) == -1);
}

class EdgeExistsException : Exception
{
    pure nothrow @nogc @safe this(
        string file = __FILE__,
        size_t line = __LINE__,
        Throwable nextInChain = null,
    )
    {
        super("edge cannot be inserted: edge already exists", file, line, nextInChain);
    }
}

class MissingEdgeException : Exception
{
    pure nothrow @nogc @safe this(
        string file = __FILE__,
        size_t line = __LINE__,
        Throwable nextInChain = null,
    )
    {
        super("edge not found", file, line, nextInChain);
    }
}

class MissingNodeException : Exception
{
    pure nothrow @nogc @safe this(
        string file = __FILE__,
        size_t line = __LINE__,
        Throwable nextInChain = null,
    )
    {
        super("node not found", file, line, nextInChain);
    }
}

/// This structure represents a graph with optional edge
/// payloads. The graph is represented as a list of edges which is
/// particularly suited for sparse graphs. While the set of nodes is fixed the
/// set of edges is mutable.
struct Graph(Node, Weight = void, Flag!"isDirected" isDirected = No.isDirected, EdgePayload = void)
{
    static enum isWeighted = !is(Weight == void);
    static enum hasEdgePayload = !is(EdgePayload == void);

    static struct Edge
    {
        protected Node _start;
        protected Node _end;

        static if (isWeighted)
            Weight weight;

        static if (hasEdgePayload)
            EdgePayload payload;

        /// Construct an edge.
        this(Node start, Node end)
        {
            this._start = start;
            this._end = end;

            static if (!isDirected)
            {
                if (end < start)
                {
                    swap(this._start, this._end);
                }
            }
        }

        static if (isWeighted)
        {
            /// ditto
            this(Node start, Node end, Weight weight)
            {
                this(start, end);
                this.weight = weight;
            }
        }

        static if (hasEdgePayload && !is(EdgePayload : Weight))
        {
            /// ditto
            this(Node start, Node end, EdgePayload payload)
            {
                this(start, end);
                this.payload = payload;
            }
        }

        static if (isWeighted && hasEdgePayload)
        {
            /// ditto
            this(Node start, Node end, Weight weight, EdgePayload payload)
            {
                this(start, end);
                this.weight = weight;
                this.payload = payload;
            }
        }

        /// Get the start of this edge. For undirected graphs this is the
        /// smaller of both incident nodes.
        @property Node start() const pure nothrow
        {
            return _start;
        }

        /// Get the end of this edge. For undirected graphs this is the
        /// larger of both incident nodes.
        @property Node end() const pure nothrow
        {
            return _end;
        }

        /**
            Get target of this edge beginning at node `from`. For undirected
            graphs returns the other node of this edge.

            Throws: MissingNodeException if this edge does not start in node `from`.
        */
        Node target(Node from) const
        {
            static if (isDirected)
            {
                if (start == from)
                {
                    return end;
                }
                else
                {
                    throw new MissingNodeException();
                }
            }
            else
            {
                if (start == from)
                {
                    return end;
                }
                else if (end == from)
                {
                    return start;
                }
                else
                {
                    throw new MissingNodeException();
                }
            }
        }

        /**
            Get source of this edge beginning at node `from`. For undirected
            graphs returns the other node of this edge.

            Throws: MissingNodeException if this edge does not end in node `from`.
        */
        static if (isDirected)
        {
            Node source(Node from) const
            {
                if (end == from)
                {
                    return start;
                }
                else
                {
                    throw new MissingNodeException();
                }
            }
        }
        else
        {
            alias source = target;
        }

        /// Two edges are equal iff their incident nodes (and weight) are the
        /// same.
        bool opEquals(in Edge other) const pure nothrow
        {
            static if (isWeighted)
            {
                return this.start == other.start && this.end == other.end
                    && this.weight == other.weight;
            }
            else
            {
                return this.start == other.start && this.end == other.end;
            }
        }

        /// Orders edge lexicographically by `start`, `end`(, `weight`).
        int opCmp(in Edge other) const pure nothrow
        {
            static if (isWeighted)
            {
                return cmpLexicographically!(
                    typeof(this),
                    "a.start",
                    "a.end",
                    "a.weight",
                )(this, other);
            }
            else
            {
                return cmpLexicographically!(
                    typeof(this),
                    "a.start",
                    "a.end",
                )(this, other);
            }
        }

        private int compareNodes(in Edge other) const pure nothrow
        {
            return cmpLexicographically!(
                typeof(this),
                "a.start",
                "a.end",
            )(this, other);
        }

        /**
            Returns the node that is connects this edge with other edge. In
            case of undirected graphs this is just the common node of both
            edges; in directed case this is the end node of this edge if it
            matches the start node of other edge.

            Throws: MissingNodeException if the connecting node is undefined.
        */
        Node getConnectingNode(in Edge other) const
        {
            static if (isDirected)
            {
                if (this.end == other.start)
                {
                    return this.end;
                }
            }
            else
            {
                if (this.end == other.start || this.end == other.end)
                {
                    return this.end;
                }
                else if (this.start == other.start || this.start == other.end)
                {
                    return this.start;
                }
            }

            throw new MissingNodeException();
        }
    }

    static bool orderByNodes(in Edge a, in Edge b) nothrow pure
    {
        return a.compareNodes(b) < 0;
    }

    static bool groupByNodes(in Edge a, in Edge b) nothrow pure
    {
        return a.compareNodes(b) == 0;
    }

    /// Construct an edge for this graph.
    static Edge edge(T...)(T args)
    {
        return Edge(args);
    }

    protected Node[] _nodes;
    protected Appender!(Edge[]) _edges;

    /// The set (ordered list) of nodes.
    @property const(Node[]) nodes() const nothrow pure
    {
        return _nodes;
    }

    private @property void nodes(Node[] nodes)
    {
        nodes.sort();

        this._nodes = nodes.uniq.array;
    }

    /// Get the set (ordered list) of edges in this graph.
    @property auto edges() nothrow pure
    {
        return chain(_edges.data);
    }

    /// ditto
    @property auto edges() const nothrow pure
    {
        return chain(_edges.data);
    }

    /**
        Construct a graph from a set of nodes (and edges). Modifies `nodes`
        while sorting but releases it after construction.

        Throws: MissingNodeException if an edge has a node that is not present
                in this graph .
        Throws: EdgeExistsException if an edge already exists when trying
                inserting it.
    */
    this(Node[] nodes)
    {
        this.nodes = nodes;
    }

    /// ditto
    this(Node[] nodes, Edge[] edges)
    {
        this(nodes);

        _edges.reserve(edges.length);
        foreach (edge; edges)
        {
            add(this, edge);
        }
    }

    this(this)
    {
        _nodes = _nodes.dup;
    }

    /// Add a set of edges to this graph without any checks.
    void bulkAddForce(R)(R edges) if (isForwardRange!R && is(ElementType!R == Edge))
    {
        this._edges ~= edges;
        _edges.data.sort;
    }

    /// Add an edge to this graph.
    /// See_Also: `Edge add(Graph, Edge)`
    void opOpAssign(string op)(Edge edge) if (op == "~")
    {
        add(this, edge);
    }

    /// Some pre-defined conflict handlers for `add`.
    static struct ConflictStrategy
    {
        static if (isWeighted)
        {
            /// Return an edge with sum of both weights. If given payload will be
            /// kept from existingEdge .
            static Edge sumWeights(Edge existingEdge, Edge newEdge)
            {
                existingEdge.weight += newEdge.weight;

                return existingEdge;
            }

            ///
            unittest
            {
                auto g1 = Graph!(int, int)([1, 2]);
                alias CS = g1.ConflictStrategy;

                g1 ~= g1.edge(1, 2, 1);

                auto addedEdge = g1.add!(CS.sumWeights)(g1.edge(1, 2, 1));

                assert(addedEdge.weight == 2);
            }
        }

        /// Throw `EdgeExistsException`.
        static inout(Edge) error(inout(Edge) existingEdge, inout(Edge) newEdge)
        {
            throw new EdgeExistsException();
        }

        ///
        unittest
        {
            auto g1 = Graph!int([1, 2]);
            alias CS = g1.ConflictStrategy;

            g1 ~= g1.edge(1, 2);

            assertThrown!EdgeExistsException(g1.add!(CS.error)(g1.edge(1, 2)));
        }

        /// Replace the existingEdge by newEdge.
        static inout(Edge) replace(inout(Edge) existingEdge, inout(Edge) newEdge)
        {
            return newEdge;
        }

        ///
        unittest
        {
            auto g1 = Graph!(int, int)([1, 2]);
            alias CS = g1.ConflictStrategy;

            g1 ~= g1.edge(1, 2, 1);

            auto addedEdge = g1.add!(CS.replace)(g1.edge(1, 2, 2));

            assert(addedEdge.weight == 2);
        }

        /// Keep existingEdge – discard newEdge.
        static inout(Edge) keep(inout(Edge) existingEdge, inout(Edge) newEdge)
        {
            return existingEdge;
        }

        ///
        unittest
        {
            auto g1 = Graph!(int, int)([1, 2]);
            alias CS = g1.ConflictStrategy;

            g1 ~= g1.edge(1, 2, 1);

            auto addedEdge = g1.add!(CS.keep)(g1.edge(1, 2, 2));

            assert(addedEdge.weight == 1);
        }
    }

    /// Forcibly add an edge to this graph.
    protected Edge forceAdd(Edge edge)
    {
        _edges ~= edge;
        _edges.data.sort;

        return edge;
    }

    /// Replace an edge in this graph.
    protected Edge replaceEdge(in size_t edgeIdx, Edge newEdge)
    {
        auto shouldSort = _edges.data[edgeIdx] != newEdge;

        _edges.data[edgeIdx] = newEdge;

        if (shouldSort)
        {
            _edges.data.sort;
        }

        return newEdge;
    }

    /// Check if edge/node exists in this graph. Ignores the weight if weighted.
    bool opBinaryRight(string op)(in Node node) const pure nothrow if (op == "in")
    {
        auto sortedNodes = assumeSorted(nodes);

        return sortedNodes.contains(node);
    }

    /// ditto
    bool has(in Node node) const pure nothrow
    {
        return node in this;
    }

    /// Check if edge exists in this graph. Only the `start` and `end` node
    /// will be compared.
    bool opBinaryRight(string op)(in Edge edge) const pure nothrow if (op == "in")
    {
        auto sortedEdges = assumeSorted!orderByNodes(edges);

        return sortedEdges.contains(edge);
    }

    /// ditto
    bool has(in Edge edge) const pure nothrow
    {
        return edge in this;
    }

    /// Get the designated edge from this graph. Only the `start` and `end`
    /// node will be compared.
    auto get(in Edge edge)
    {
        auto sortedEdges = assumeSorted!orderByNodes(edges);
        auto existingEdges = sortedEdges.equalRange(edge);

        if (existingEdges.empty)
        {
            throw new MissingEdgeException();
        }
        else
        {
            return existingEdges.front;
        }
    }

    ///
    unittest
    {
        auto g1 = Graph!(int, int)([1, 2]);

        auto e1 = g1.edge(1, 2, 1);

        g1 ~= e1;

        assert(g1.get(g1.edge(1, 2)) == e1);
        assertThrown!MissingEdgeException(g1.get(g1.edge(1, 1)));
    }

    /// Returns the index of node `n` in the list of nodes.
    size_t indexOf(in Node n) const
    {
        auto sortedNodes = assumeSorted(nodes);
        auto tristectedNodes = sortedNodes.trisect(n);

        if (tristectedNodes[1].empty)
        {
            throw new MissingNodeException();
        }

        return tristectedNodes[0].length;
    }

    ///
    unittest
    {
        auto g1 = Graph!(int, int)([1, 2]);

        assert(g1.indexOf(1) == 0);
        assert(g1.indexOf(2) == 1);
        assertThrown!MissingNodeException(g1.indexOf(3));
    }

    /// Returns the index of node `n` in the list of nodes.
    size_t indexOf(in Edge edge) const
    {
        auto sortedEdges = assumeSorted!orderByNodes(edges);
        auto trisectedEdges = sortedEdges.trisect(edge);

        if (trisectedEdges[1].empty)
        {
            throw new MissingEdgeException();
        }

        return trisectedEdges[0].length;
    }

    ///
    unittest
    {
        auto g1 = Graph!(int, int)([1, 2]);

        auto e1 = g1.edge(1, 2, 1);

        g1 ~= e1;

        assert(g1.indexOf(g1.edge(1, 2)) == 0);
        assertThrown!MissingEdgeException(g1.indexOf(g1.edge(1, 1)));
    }

    static if (isDirected)
    {
        /// Returns a range of in/outgoing edges of node `n`.
        auto inEdges(Node n) nothrow pure
        {
            return _edges.data[].filter!(e => e.end == n);
        }

        /// ditto
        auto inEdges(Node n) const nothrow pure
        {
            return edges[].filter!(e => e.end == n);
        }

        /// ditto
        auto outEdges(Node n) nothrow pure
        {
            return _edges.data[].filter!(e => e.start == n);
        }

        /// ditto
        auto outEdges(Node n) const nothrow pure
        {
            return edges[].filter!(e => e.start == n);
        }

        ///
        unittest
        {
            import std.algorithm : equal;

            auto g1 = Graph!(int, void, Yes.isDirected)([1, 2, 3]);

            g1 ~= g1.edge(1, 1);
            g1 ~= g1.edge(1, 2);
            g1 ~= g1.edge(2, 2);
            g1 ~= g1.edge(2, 3);

            assert(g1.inEdges(1).equal([
                g1.edge(1, 1),
            ]));
            assert(g1.outEdges(1).equal([
                g1.edge(1, 1),
                g1.edge(1, 2),
            ]));
            assert(g1.inEdges(2).equal([
                g1.edge(1, 2),
                g1.edge(2, 2),
            ]));
            assert(g1.outEdges(2).equal([
                g1.edge(2, 2),
                g1.edge(2, 3),
            ]));
            assert(g1.inEdges(3).equal([
                g1.edge(2, 3),
            ]));
            assert(g1.outEdges(3).empty);
        }

        /// Get the in/out degree of node `n`.
        size_t inDegree(Node n) const nothrow pure
        {
            return inEdges(n).walkLength;
        }

        /// ditto
        size_t outDegree(Node n) const nothrow pure
        {
            return outEdges(n).walkLength;
        }

        ///
        unittest
        {
            auto g1 = Graph!(int, void, Yes.isDirected)([1, 2, 3]);

            g1 ~= g1.edge(1, 1);
            g1 ~= g1.edge(1, 2);
            g1 ~= g1.edge(2, 2);
            g1 ~= g1.edge(2, 3);

            assert(g1.inDegree(1) == 1);
            assert(g1.outDegree(1) == 2);
            assert(g1.inDegree(2) == 2);
            assert(g1.outDegree(2) == 2);
            assert(g1.inDegree(3) == 1);
            assert(g1.outDegree(3) == 0);
        }
    }
    else
    {
        /// Returns a range of all edges incident to node `n`.
        auto incidentEdges(Node n) nothrow pure
        {
            return _edges.data[].filter!(e => e.start == n || e.end == n);
        }

        /// ditto
        auto incidentEdges(Node n) const nothrow pure
        {
            return edges[].filter!(e => e.start == n || e.end == n);
        }

        /// ditto
        alias inEdges = incidentEdges;

        /// ditto
        alias outEdges = incidentEdges;

        ///
        unittest
        {
            import std.algorithm : equal;

            auto g1 = Graph!int([1, 2, 3]);

            g1 ~= g1.edge(1, 1);
            g1 ~= g1.edge(1, 2);
            g1 ~= g1.edge(2, 2);
            g1 ~= g1.edge(2, 3);

            assert(g1.incidentEdges(1).equal([
                g1.edge(1, 1),
                g1.edge(1, 2),
            ]));
            assert(g1.incidentEdges(2).equal([
                g1.edge(1, 2),
                g1.edge(2, 2),
                g1.edge(2, 3),
            ]));
            assert(g1.incidentEdges(3).equal([
                g1.edge(2, 3),
            ]));
        }

        IncidentEdgesCache allIncidentEdges()
        {
            return IncidentEdgesCache(this);
        }

        static struct IncidentEdgesCache
        {
            alias G = Graph!(Node, Weight, isDirected, EdgePayload);
            G graph;
            Edge[][] incidentEdges;

            this(G graph)
            {
                this.graph = graph;
                collectAllIncidentEdges();
            }

            private void collectAllIncidentEdges()
            {
                preallocateMemory();

                size_t startIdx;
                size_t endIdx;
                foreach (edge; graph._edges.data)
                {
                    if (graph._nodes[startIdx] < edge.start)
                        endIdx = startIdx;
                    while (graph._nodes[startIdx] < edge.start)
                        ++startIdx;
                    if (endIdx < startIdx)
                        endIdx = startIdx;
                    while (graph._nodes[endIdx] < edge.end)
                        ++endIdx;

                    incidentEdges[startIdx] ~= edge;
                    // Avoid double-counting of loops
                    if (startIdx != endIdx)
                        incidentEdges[endIdx] ~= edge;
                }
            }

            void preallocateMemory()
            {
                auto degreesCache = graph.allDegrees();
                Edge[] buffer;
                buffer.length = degreesCache.degrees.sum;
                incidentEdges.length = degreesCache.degrees.length;

                size_t sliceBegin;
                size_t startIdx;
                foreach (degree; degreesCache)
                {
                    incidentEdges[startIdx] = buffer[sliceBegin .. sliceBegin + degree];
                    incidentEdges[startIdx].length = 0;

                    sliceBegin += degree;
                    ++startIdx;
                }
            }

            ref inout(Edge[]) opIndex(in Node node) inout
            {
                return incidentEdges[graph.indexOf(node)];
            }

            ref inout(Edge[]) opIndex(in size_t nodeIdx) inout
            {
                return incidentEdges[nodeIdx];
            }

            int opApply(scope int delegate(Edge[]) yield)
            {
                int result = 0;

                foreach (currentIncidentEdges; incidentEdges)
                {
                    result = yield(currentIncidentEdges);
                    if (result)
                        break;
                }

                return result;
            }

            int opApply(scope int delegate(Node, Edge[]) yield)
            {
                int result = 0;

                foreach (i, currentIncidentEdges; incidentEdges)
                {
                    result = yield(graph._nodes[i], currentIncidentEdges);
                    if (result)
                        break;
                }

                return result;
            }
        }

        /// Get the `adjacencyList` of this graph where nodes are represented
        /// by their index in the nodes list.
        size_t[][] adjacencyList() const
        {
            size_t[][] _adjacencyList;
            _adjacencyList.length = nodes.length;
            size_t[] targetsBuffer;
            targetsBuffer.length = 2 * edges.length;

            foreach (i, node; _nodes)
            {
                auto bufferRest = edges
                    .filter!(e => e.start == node || e.end == node)
                    .map!(edge => indexOf(edge.target(node)))
                    .copy(targetsBuffer);
                _adjacencyList[i] = targetsBuffer[0 .. $ - bufferRest.length];
                _adjacencyList[i].sort;
                targetsBuffer = bufferRest;
            }

            return _adjacencyList;
        }

        ///
        unittest
        {
            auto g1 = Graph!int([1, 2, 3, 4]);

            g1 ~= g1.edge(1, 1);
            g1 ~= g1.edge(1, 2);
            g1 ~= g1.edge(2, 2);
            g1 ~= g1.edge(2, 3);
            g1 ~= g1.edge(2, 4);
            g1 ~= g1.edge(3, 4);

            assert(g1.adjacencyList() == [
                [0, 1],
                [0, 1, 2, 3],
                [1, 3],
                [1, 2],
            ]);
        }

        /// Get the degree of node `n`.
        size_t degree(Node n) const nothrow pure
        {
            return incidentEdges(n).walkLength;
        }

        /// ditto
        alias inDegree = degree;

        /// ditto
        alias outDegree = degree;

        DegreesCache allDegrees() const
        {
            return DegreesCache(this);
        }

        static struct DegreesCache
        {
            alias G = Graph!(Node, Weight, isDirected, EdgePayload);
            const(G) graph;
            size_t[] degrees;

            this(in G graph)
            {
                this.graph = graph;
                collectAllDegrees();
            }

            private void collectAllDegrees()
            {
                degrees.length = graph._nodes.length;

                size_t startIdx;
                size_t endIdx;
                foreach (edge; graph._edges.data)
                {
                    if (graph._nodes[startIdx] < edge.start)
                        endIdx = startIdx;
                    while (graph._nodes[startIdx] < edge.start)
                        ++startIdx;
                    if (endIdx < startIdx)
                        endIdx = startIdx;
                    while (graph._nodes[endIdx] < edge.end)
                        ++endIdx;

                    ++degrees[startIdx];
                    // Avoid double-counting of loops
                    if (startIdx != endIdx)
                        ++degrees[endIdx];
                }
            }

            size_t opIndex(in Node node) const
            {
                return degrees[graph.indexOf(node)];
            }

            size_t opIndex(in size_t nodeIdx) const
            {
                return degrees[nodeIdx];
            }

            int opApply(scope int delegate(size_t) yield) const
            {
                int result = 0;

                foreach (degree; degrees)
                {
                    result = yield(degree);
                    if (result)
                        break;
                }

                return result;
            }

            int opApply(scope int delegate(Node, size_t) yield) const
            {
                int result = 0;

                foreach (i, degree; degrees)
                {
                    result = yield(graph._nodes[i], degree);
                    if (result)
                        break;
                }

                return result;
            }
        }
    }
}

///
unittest
{
    //   +-+  +-+
    //   \ /  \ /
    //   (1)--(2)
    auto g1 = Graph!int([1, 2]);

    g1 ~= g1.edge(1, 1);
    g1 ~= g1.edge(1, 2);
    g1.add(g1.edge(2, 2));

    assert(g1.edge(1, 1) in g1);
    assert(g1.edge(1, 2) in g1);
    assert(g1.edge(2, 1) in g1);
    assert(g1.has(g1.edge(2, 2)));
    assert(g1.allDegrees().degrees == [2, 2]);
    assert(g1.allIncidentEdges().incidentEdges == [
        [g1.edge(1, 1), g1.edge(1, 2)],
        [g1.edge(1, 2), g1.edge(2, 2)],
    ]);

    //   0.5     0.5
    //   +-+     +-+
    //   \ /     \ /
    //   (1)-----(2)
    //       1.0
    auto g2 = Graph!(int, double)([1, 2]);

    g2 ~= g2.edge(1, 1, 0.5);
    g2 ~= g2.edge(1, 2, 1.0);
    g2.add(g2.edge(2, 2, 0.5));

    assert(g2.edge(1, 1) in g2);
    assert(g2.edge(1, 2) in g2);
    assert(g2.edge(2, 1) in g2);
    assert(g2.has(g2.edge(2, 2)));
    assert(g2.allDegrees().degrees == [2, 2]);
    assert(g2.allIncidentEdges().incidentEdges == [
        [g2.edge(1, 1, 0.5), g2.edge(1, 2, 1.0)],
        [g2.edge(1, 2, 1.0), g2.edge(2, 2, 0.5)],
    ]);

    //   0.5     0.5
    //   +-+     +-+
    //   \ v     v /
    //   (1)---->(2)
    //       1.0
    auto g3 = Graph!(int, double, Yes.isDirected)([1, 2]);

    g3 ~= g3.edge(1, 1, 0.5);
    g3 ~= g3.edge(1, 2, 1.0);
    g3.add(g3.edge(2, 2, 0.5));

    assert(g3.edge(1, 1) in g3);
    assert(g3.edge(1, 2) in g3);
    assert(!(g3.edge(2, 1) in g3));
    assert(g3.has(g3.edge(2, 2)));

    //   +-+   +-+
    //   \ v   v /
    //   (1)-->(2)
    auto g4 = Graph!(int, void, Yes.isDirected)([1, 2]);

    g4 ~= g4.edge(1, 1);
    g4 ~= g4.edge(1, 2);
    g4.add(g4.edge(2, 2));

    assert(g4.edge(1, 1) in g4);
    assert(g4.edge(1, 2) in g4);
    assert(!(g4.edge(2, 1) in g4));
    assert(g4.has(g4.edge(2, 2)));

    //   +-+  +-+
    //   \ /  \ /
    //   (1)--(2)
    //
    // payload(1, 1) = [1];
    // payload(1, 2) = [2];
    // payload(2, 2) = [3];
    auto g5 = Graph!(int, void, No.isDirected, int[])([1, 2]);

    g5 ~= g5.edge(1, 1, [1]);
    g5 ~= g5.edge(1, 2, [2]);
    g5.add(g5.edge(2, 2, [3]));

    assert(g5.edge(1, 1) in g5);
    assert(g5.get(g5.edge(1, 1)).payload == [1]);
    assert(g5.edge(1, 2) in g5);
    assert(g5.get(g5.edge(1, 2)).payload == [2]);
    assert(g5.edge(2, 1) in g5);
    assert(g5.get(g5.edge(2, 1)).payload == [2]);
    assert(g5.has(g5.edge(2, 2)));
    assert(g5.get(g5.edge(2, 2)).payload == [3]);
    assert(g5.allDegrees().degrees == [2, 2]);
    assert(g5.allIncidentEdges().incidentEdges == [
        [g5.edge(1, 1), g5.edge(1, 2)],
        [g5.edge(1, 2), g5.edge(2, 2)],
    ]);
}

///
unittest
{
    //     -1     1         1
    // (1)----(2)---(3) (4)---(5) (6)
    size_t[] contigs = [1, 2, 3, 4, 5, 6];
    auto contigGraph = Graph!(size_t, int)([1, 2, 3, 4, 5, 6]);

    contigGraph.add(contigGraph.edge(1, 2, -1));
    contigGraph.add(contigGraph.edge(2, 3, 1));
    contigGraph.add(contigGraph.edge(4, 5, 1));

    foreach (contig; contigs)
    {
        assert(contigGraph.degree(contig) <= 2);
    }
    assert(contigGraph.allDegrees().degrees == [1, 2, 1, 1, 1, 0]);
    assert(contigGraph.allIncidentEdges().incidentEdges == [
        [contigGraph.edge(1, 2, -1)],
        [contigGraph.edge(1, 2, -1), contigGraph.edge(2, 3, 1)],
        [contigGraph.edge(2, 3, 1)],
        [contigGraph.edge(4, 5, 1)],
        [contigGraph.edge(4, 5, 1)],
        [],
    ]);
}

/// Add a set of edges to this graph and merge mutli-edges using `merge`.
void bulkAdd(alias merge, G, R)(ref G graph, R edges)
        if (is(G : Graph!Params, Params...) && isForwardRange!R && is(ElementType!R == G.Edge))
{
    alias Edge = G.Edge;
    alias ReturnTypeMerge = typeof(merge(new Edge[0]));
    static assert(is(ReturnTypeMerge == Edge), "expected `Edge merge(Edge[] multiEdge)`");

    graph.bulkAddForce(edges);

    auto bufferRest = graph
        ._edges
        .data
        .sliceBy!(G.groupByNodes)
        .map!(unaryFun!merge)
        .copy(graph._edges.data);
    graph._edges.shrinkTo(graph._edges.data.length - bufferRest.length);
}

///
unittest
{
    auto g1 = Graph!(int, int)([1, 2]);

    static g1.Edge sumWeights(g1.Edge[] multiEdge)
    {
        auto sumOfWeights = multiEdge.map!"a.weight".sum;
        auto mergedEdge = multiEdge[0];
        mergedEdge.weight = sumOfWeights;

        return mergedEdge;
    }

    auto edges = [
        g1.edge(1, 2, 1),
        g1.edge(1, 2, 1),
        g1.edge(1, 2, 1),
        g1.edge(2, 3, 2),
        g1.edge(2, 3, 2),
        g1.edge(3, 4, 3),
    ];
    g1.bulkAdd!sumWeights(edges);
    assert(g1.edges == [
        g1.edge(1, 2, 3),
        g1.edge(2, 3, 4),
        g1.edge(3, 4, 3),
    ]);
}

/// Add an edge to this graph and handle existing edges with `handleConflict`.
/// The handler must have this signature `Edge handleConflict(Edge, Edge)`.
G.Edge add(alias handleConflict = 1337, G)(ref G graph, G.Edge edge)
        if (is(G : Graph!Params, Params...))
{
    static if (isCallable!handleConflict)
        alias handleConflict_ = binaryFun!handleConflict;
    else
        alias handleConflict_ = binaryFun!(G.ConflictStrategy.error);

    if (!graph.has(edge.start) || !graph.has(edge.end))
    {
        throw new MissingNodeException();
    }

    auto sortedEdges = assumeSorted!(G.orderByNodes)(graph._edges.data);
    auto trisectedEdges = sortedEdges.trisect(edge);
    auto existingEdges = trisectedEdges[1];
    auto existingEdgeIdx = trisectedEdges[0].length;

    if (existingEdges.empty)
    {
        return graph.forceAdd(edge);
    }
    else
    {
        auto newEdge = handleConflict_(existingEdges.front, edge);

        return graph.replaceEdge(existingEdgeIdx, newEdge);
    }
}

///
unittest
{
    auto g1 = Graph!(int, int)([1, 2]);

    auto e1 = g1.edge(1, 2, 1);
    auto e2 = g1.edge(1, 2, 2);

    g1 ~= e1;

    assertThrown!EdgeExistsException(g1.add(e2));

    with (g1.ConflictStrategy)
    {
        g1.add!replace(e2);

        assert(g1.get(g1.edge(1, 2)) == e2);

        g1.add!keep(e1);

        assert(g1.get(g1.edge(1, 2)) == e2);

        g1.add!sumWeights(e2);

        assert(g1.get(g1.edge(1, 2)).weight == 2 * e2.weight);
    }
}

void filterEdges(alias pred, G)(ref G graph) if (is(G : Graph!Params, Params...))
{
    auto bufferRest = graph
        ._edges
        .data
        .filter!pred
        .copy(graph._edges.data);
    graph._edges.shrinkTo(graph._edges.data.length - bufferRest.length);
}

void mapEdges(alias fun, G)(ref G graph) if (is(G : Graph!Params, Params...))
{
    foreach (ref edge; graph._edges.data)
        edge = unaryFun!fun(edge);

    graph._edges.data.sort();
}

class EmptySetException : Exception
{
    this(string msg)
    {
        super(msg);
    }
}

struct NaturalNumberSet
{
    private static enum partSize = 8 * size_t.sizeof;
    private static enum size_t firstBit = 1;
    private static enum size_t lastBit = firstBit << (partSize - 1);
    private static enum size_t emptyPart = 0;
    private static enum size_t fullPart = ~emptyPart;

    private size_t[] parts;
    private size_t nMax;

    this(size_t initialNumElements, Flag!"addAll" addAll = No.addAll)
    {
        reserveFor(initialNumElements);

        if (addAll)
        {
            foreach (i; 0 .. initialNumElements / partSize)
                parts[i] = fullPart;
            foreach (i; initialNumElements / partSize .. initialNumElements)
                add(i);
        }
    }

    static NaturalNumberSet create(size_t[] initialElements...)
    {
        auto set = NaturalNumberSet(initialElements.maxElement);

        foreach (i; initialElements)
            set.add(i);

        return set;
    }

    this(this)
    {
        parts = parts.dup;
    }

    private this(size_t[] parts)
    {
        this.parts = parts;
    }

    private bool inBounds(in size_t n) const pure nothrow
    {
        return n < nMax;
    }

    void reserveFor(in size_t n)
    {
        if (parts.length == 0)
        {
            parts.length = max(1, ceildiv(n, partSize));
            nMax = parts.length * partSize;
        }

        while (!inBounds(n))
        {
            parts.length *= 2;
            nMax = parts.length * partSize;
        }
    }

    @property size_t capacity() pure const nothrow
    {
        return nMax;
    }

    private size_t partIdx(in size_t n) const pure nothrow
    {
        return n / partSize;
    }

    private size_t idxInPart(in size_t n) const pure nothrow
    {
        return n % partSize;
    }

    private size_t itemMask(in size_t n) const pure nothrow
    {
        return firstBit << idxInPart(n);
    }

    static size_t inverse(in size_t n) pure nothrow
    {
        return n ^ fullPart;
    }

    void add(in size_t n)
    {
        reserveFor(n);

        parts[partIdx(n)] |= itemMask(n);
    }

    void remove(in size_t n)
    {
        if (!inBounds(n))
        {
            return;
        }

        parts[partIdx(n)] &= inverse(itemMask(n));
    }

    bool has(in size_t n) const pure nothrow
    {
        if (!inBounds(n))
        {
            return false;
        }

        return (parts[partIdx(n)] & itemMask(n)) != emptyPart;
    }

    bool empty() const pure nothrow
    {
        return parts.all!(part => part == emptyPart);
    }

    void clear() pure nothrow
    {
        foreach (ref part; parts)
            part = emptyPart;
    }

    bool opBinary(string op)(in NaturalNumberSet other) const pure nothrow if (op == "==")
    {
        auto numCommonParts = min(this.parts.length, other.parts.length);

        foreach (i; 0 .. numCommonParts)
        {
            if (this.parts[i] != other.parts[i])
                return false;
        }

        static bool hasEmptyTail(ref in NaturalNumberSet set, in size_t tailStart)
        {
            foreach (i; tailStart .. set.parts.length)
                if (set.parts[i] != emptyPart)
                    return false;

            return true;
        }

        if (this.parts.length > numCommonParts)
            return hasEmptyTail(this, numCommonParts);
        if (other.parts.length > numCommonParts)
            return hasEmptyTail(other, numCommonParts);

        return true;
    }

    bool opBinary(string op)(in NaturalNumberSet other) const pure nothrow if (op == "in")
    {
        auto numCommonParts = min(this.parts.length, other.parts.length);

        foreach (i; 0 .. numCommonParts)
            if ((this.parts[i] & other.parts[i]) != this.parts[i])
                return false;

        static bool hasEmptyTail(ref in NaturalNumberSet set, in size_t tailStart)
        {
            foreach (i; tailStart .. set.parts.length)
                if (set.parts[i] != emptyPart)
                    return false;

            return true;
        }

        if (this.parts.length > numCommonParts)
            return hasEmptyTail(this, numCommonParts);

        return true;
    }

    NaturalNumberSet opBinary(string op)(in NaturalNumberSet other) const pure nothrow if (op.among("|", "^", "&"))
    {
        NaturalNumberSet result;
        result.parts.length = max(this.parts.length, other.parts.length);
        result.nMax = max(this.nMax, other.nMax);

        auto numCommonParts = min(this.parts.length, other.parts.length);

        foreach (i; 0 .. numCommonParts)
            result.parts[i] = mixin("this.parts[i] " ~ op ~ " other.parts[i]");

        static if (op.among("|", "^"))
        {
            if (this.parts.length > numCommonParts)
                result.parts[numCommonParts .. $] = this.parts[numCommonParts .. $];
            if (other.parts.length > numCommonParts)
                result.parts[numCommonParts .. $] = other.parts[numCommonParts .. $];
        }

        return result;
    }

    bool intersects(in NaturalNumberSet other) const pure nothrow
    {
        auto numCommonParts = min(this.parts.length, other.parts.length);

        foreach (i; 0 .. numCommonParts)
        {
            if ((this.parts[i] & other.parts[i]) != emptyPart)
                return true;
        }

        return false;
    }

    @property size_t size() const pure nothrow
    {
        size_t numSetBits;

        foreach (i, part; parts)
        {
            size_t j = 0;

            while ((part >> j) != emptyPart)
            {
                while (((part >> j) & firstBit) != firstBit)
                    ++j;
                ++numSetBits;
                ++j;
            }
        }

        return numSetBits;
    }

    size_t minElement() const
    {
        foreach (i, part; parts)
        {
            if (part != emptyPart)
            {
                size_t j = 0;

                while (((part >> j) & firstBit) != firstBit)
                {
                    ++j;
                }

                return i * partSize + j;
            }
        }

        throw new EmptySetException("empty set has no minElement");
    }

    size_t maxElement() const
    {
        foreach (i, part; parts.retro.enumerate)
        {
            if (part != emptyPart)
            {
                size_t j = 0;

                while (((part << j) & lastBit) != lastBit)
                {
                    ++j;
                }

                return (parts.length - i - 1) * partSize + (partSize - j - 1);
            }
        }

        throw new EmptySetException("empty set has no maxElement");
    }

    unittest
    {
        foreach (i; 0 .. 2 * NaturalNumberSet.partSize)
        {
            NaturalNumberSet set;

            set.add(i + 5);
            set.add(i + 7);

            assert(set.minElement() == i + 5);
            assert(set.maxElement() == i + 7);
        }
    }

    /// Returns a range of the elements in this set. The elements are ordered
    /// ascending.
    @property auto elements() const pure nothrow
    {
        static struct ElementsRange
        {
            const NaturalNumberSet* set;
            bool _empty = false;
            size_t i = 0;
            size_t part;
            size_t j = 0;

            this(const NaturalNumberSet* set) pure nothrow
            {
                this.set = set;
                this._empty = set.empty;

                if (!this.empty)
                {
                    this.part = set.parts[i];
                    if (!set.has(front))
                    {
                        popFront();
                    }
                }
            }

            @property ElementsRange save() const pure nothrow
            {
                return this;
            }

            void popFront() pure nothrow
            {
                assert(!empty, "Attempting to popFront an empty elements range");
                ++j;

                while (shiftedPartEmpty)
                {
                    nextPart();

                    if (empty)
                    {
                        return;
                    }
                }

                while (((part >> j) & firstBit) != firstBit && !shiftedPartEmpty)
                {
                    ++j;
                }

                if (shiftedPartEmpty)
                {
                    popFront();
                }
            }

            @property size_t front() const pure nothrow
            {
                assert(!empty, "Attempting to fetch the front of an empty elements range");
                return i * partSize + j;
            }

            @property bool empty() const pure nothrow
            {
                return _empty;
            }

            private @property bool shiftedPartEmpty() const pure nothrow
            {
                return (part >> j) == emptyPart || j >= partSize;
            }

            private void nextPart() pure nothrow
            {
                // move to start of next part
                ++i;
                j = 0;

                if (i < set.parts.length)
                {
                    part = set.parts[i];
                }
                else
                {
                    _empty = true;
                }
            }
        }

        return ElementsRange(&this);
    }

    ///
    unittest
    {
        import std.algorithm : equal;
        import std.range : iota;

        NaturalNumberSet set;
        auto someNumbers = iota(set.partSize).filter!"a % 3 == 0";

        foreach (i; someNumbers)
        {
            set.add(i);
        }

        assert(equal(someNumbers, set.elements));
    }

    /// The set may be modified while iterating:
    unittest
    {
        import std.algorithm : equal;
        import std.range : iota;

        enum numElements = 64;
        auto set = NaturalNumberSet(numElements, Yes.addAll);

        foreach (i; set.elements)
        {
            if (i % 10 == 0)
                set.remove(i + 1);
        }

        auto expectedNumbers = iota(numElements).filter!"a == 0 || !((a - 1) % 10 == 0)";
        assert(equal(expectedNumbers, set.elements));
    }

    string toString() const pure
    {
        return format("[%(%d,%)]", this.elements);
    }
}

unittest
{
    NaturalNumberSet set;

    // add some numbers
    foreach (i; 0 .. set.partSize)
    {
        if (i % 2 == 0)
        {
            set.add(i);
        }
    }

    // force extension of set
    foreach (i; set.partSize .. 2 * set.partSize)
    {
        if (i % 3 == 0)
        {
            set.add(i);
        }
    }

    // validate presence
    foreach (i; 0 .. 2 * set.partSize)
    {
        if (i / set.partSize == 0 && i % 2 == 0)
        {
            assert(set.has(i));
        }
        else if (i / set.partSize == 1 && i % 3 == 0)
        {
            assert(set.has(i));
        }
        else
        {
            assert(!set.has(i));
        }
    }
}

/**
    Find all maximal connected components of a graph-like structure. The
    predicate `isConnected` will be evaluated `O(n^^2)` times in the
    worst-case and `Ω(n)` in the best case. In expectation it will be
    evaluated `θ(n*log(n))`.

    Params:
        isConnected =   binary predicate that evaluates to true iff two nodes,
                        represented as indices, are connected
        numNodes    =   total number of nodes in the graph

    Returns:    range of maxmimally connected components represented as
                `NaturalNumberSet`s
*/
auto findMaximallyConnectedComponents(alias isConnected)(in size_t numNodes)
{
    return MaximalConnectedComponents!(binaryFun!isConnected)(numNodes);
}

///
unittest
{
    import std.algorithm : equal;
    import std.range : only;

    alias modEqv(size_t m) = (a, b) => (a % m) == (b % m);
    alias clusterByThreshold(size_t t) = (a, b) => (a < t) == (b < t);

    assert(equal(
        findMaximallyConnectedComponents!(modEqv!5)(15),
        only(
            NaturalNumberSet.create(0, 5, 10),
            NaturalNumberSet.create(1, 6, 11),
            NaturalNumberSet.create(2, 7, 12),
            NaturalNumberSet.create(3, 8, 13),
            NaturalNumberSet.create(4, 9, 14),
        ),
    ));
    assert(equal(
        findMaximallyConnectedComponents!(modEqv!3)(15),
        only(
            NaturalNumberSet.create(0, 3, 6, 9, 12),
            NaturalNumberSet.create(1, 4, 7, 10, 13),
            NaturalNumberSet.create(2, 5, 8, 11, 14),
        ),
    ));
    assert(equal(
        findMaximallyConnectedComponents!(clusterByThreshold!10)(15),
        only(
            NaturalNumberSet.create(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
            NaturalNumberSet.create(10, 11, 12, 13, 14),
        ),
    ));
}

///
unittest
{
    import std.algorithm : equal;
    import std.range : only;

    auto connectivity = [
        [false, false, false, true ],
        [false, false, true , false],
        [false, true , false, false],
        [true , false, false, false],
    ];
    alias isConnected = (i, j) => connectivity[i][j];

    assert(equal(
        findMaximallyConnectedComponents!isConnected(4),
        only(
            NaturalNumberSet.create(0, 3),
            NaturalNumberSet.create(1, 2),
        ),
    ));
}

private struct MaximalConnectedComponents(alias isConnected)
{

    const(size_t) numNodes;
    NaturalNumberSet unvisited;
    NaturalNumberSet currentComponent;

    this(in size_t numNodes)
    {
        this.numNodes = numNodes;
        this.unvisited = NaturalNumberSet(numNodes, Yes.addAll);
        this.currentComponent = NaturalNumberSet(numNodes);

        if (!empty)
            popFront();
    }

    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty " ~ typeof(this).stringof);

        currentComponent.clear();

        if (unvisited.empty)
            return;

        auto seedNode = unvisited.minElement;

        maximizeConnectedComponent(seedNode);
    }

    private void maximizeConnectedComponent(size_t node)
    {
        currentComponent.add(node);
        unvisited.remove(node);

        foreach (nextNode; unvisited.elements)
            if (isConnected(node, nextNode))
                maximizeConnectedComponent(nextNode);
    }

    @property NaturalNumberSet front()
    {
        assert(!empty, "Attempting to fetch the front an empty " ~ typeof(this).stringof);

        return currentComponent;
    }

    @property bool empty() const pure nothrow
    {
        return unvisited.empty && currentComponent.empty;
    }
}

/**
    Find all maximal cliques in a graph represented by `adjacencyList`.
    The implementation is based on version 1 of the Bron-Kerbosch algorithm [1].

    [1]: Bron, C.; Kerbosch, J. (1973), "Algorithm 457: finding all cliques
         of an undirected graph", Communications of the ACM, 16 (9): 575–577,
         doi:10.1145/362342.362367.

    Returns: list of sets of nodes each representing a maximal clique
*/
auto findAllCliques(in size_t[][] adjacencyList)
{
    return BronKerboschVersion1(adjacencyList);
}

///
unittest
{
    auto g = Graph!int([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    g.add(g.edge(0, 1));
    g.add(g.edge(0, 2));
    g.add(g.edge(1, 2));
    g.add(g.edge(1, 7));
    g.add(g.edge(1, 8));
    g.add(g.edge(2, 3));
    g.add(g.edge(3, 4));
    g.add(g.edge(3, 5));
    g.add(g.edge(3, 6));
    g.add(g.edge(4, 5));
    g.add(g.edge(4, 6));
    g.add(g.edge(5, 6));
    g.add(g.edge(6, 7));
    g.add(g.edge(7, 8));

    auto cliques = array(findAllCliques(g.adjacencyList()));

    assert(cliques == [
        [0, 1, 2],
        [1, 7, 8],
        [2, 3],
        [3, 4, 5, 6],
        [6, 7],
        [9],
    ]);
}

private struct BronKerboschVersion1
{
    const size_t[][] adjacencyList;

    int opApply(scope int delegate(size_t[]) yield)
    {
        size_t[] clique;
        clique.reserve(adjacencyList.length);

        auto candidates = NaturalNumberSet(adjacencyList.length, Yes.addAll);
        auto not = NaturalNumberSet(adjacencyList.length);

        return extendClique(clique, candidates, not, yield);
    }

    private int extendClique(
        size_t[] clique,
        NaturalNumberSet candidates,
        NaturalNumberSet not,
        scope int delegate(size_t[]) yield,
    )
    {
        import std.stdio;

        if (not.empty && candidates.empty)
            return clique.length == 0 ? 0 : yield(clique);

        int result;

        foreach (candidate; candidates.elements)
        {
            clique ~= candidate;

            auto reducedCandidates = NaturalNumberSet(adjacencyList.length);
            auto reducedNot = NaturalNumberSet(adjacencyList.length);

            foreach (neighbourNode; adjacencyList[candidate])
            {
                if (candidates.has(neighbourNode))
                    reducedCandidates.add(neighbourNode);
                if (not.has(neighbourNode))
                    reducedNot.add(neighbourNode);
            }

            result = extendClique(clique, reducedCandidates, reducedNot, yield);

            if (result)
                return result;

            candidates.remove(candidate);
            not.add(candidate);
            --clique.length;
        }

        return result;
    }
}

/**
    Calculate a longest increasing subsequence of `sequence`. This subsequence
    is not necessarily contiguous, or unique. Given a `sequence` of `n`
    elements the algorithm uses `O(n log n)` evaluation of `pred`.

    See_Also: https://en.wikipedia.org/wiki/Longest_increasing_subsequence
*/
auto longestIncreasingSubsequence(alias pred = "a < b", Range)(Range sequence)
        if (isRandomAccessRange!Range)
{
    alias lessThan = binaryFun!pred;

    size_t[] subseqEnds;
    subseqEnds.length = sequence.length;
    size_t[] predecessors;
    predecessors.length = sequence.length;
    size_t subseqLength;

    foreach (i; 0 .. sequence.length)
    {
        // Binary search for the largest positive j < subseqLength
        // such that sequence[subseqEnds[j]] < sequence[i]
        long lo = 0;
        long hi = subseqLength - 1;
        auto pivot = sequence[i];
        assert(!lessThan(pivot, pivot), "`pred` is not anti-symmetric");

        while (lo <= hi)
        {
            auto mid = ceildiv(lo + hi, 2);

            if (lessThan(sequence[subseqEnds[mid]], pivot))
                lo = mid + 1;
            else
                hi = mid - 1;
        }

        // After searching, lo + 1 is the length of the longest prefix of
        // sequence[i]
        auto newSubseqLength = lo + 1;

        // The predecessor of sequence[i] is the last index of
        // the subsequence of length newSubseqLength - 1
        subseqEnds[lo] = i;
        if (lo > 0)
            predecessors[i] = subseqEnds[lo - 1];

        if (newSubseqLength > subseqLength)
            // If we found a subsequence longer than any we've
            // found yet, update subseqLength
            subseqLength = newSubseqLength;
    }

    auto subsequenceResult = subseqEnds[0 .. subseqLength];

    if (subseqLength > 0)
    {
        // Reconstruct the longest increasing subsequence
        // Note: reusing memory from now unused subseqEnds
        auto k = subseqEnds[subseqLength - 1];
        foreach_reverse (i; 0 .. subseqLength)
        {
            subsequenceResult[i] = k;
            k = predecessors[k];
        }
    }

    return subsequenceResult.map!(i => sequence[i]);
}

/// Example from Wikipedia
unittest
{
    import std.algorithm : equal;

    auto inputSequence = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];
    auto expectedOutput = [0, 2, 6, 9, 11, 15];

    assert(inputSequence.longestIncreasingSubsequence.equal(expectedOutput));
}

/// Example using a different `pred`
unittest
{
    import std.algorithm : equal;
    import std.range : retro;

    auto inputSequence = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];
    auto expectedOutput = [12, 10, 9, 5, 3];

    assert(inputSequence.longestIncreasingSubsequence!"a > b".equal(expectedOutput));
}

unittest
{
    import std.algorithm : equal;

    int[] inputSequence = [];
    int[] expectedOutput = [];

    assert(inputSequence.longestIncreasingSubsequence.equal(expectedOutput));
}

unittest
{
    import std.algorithm : equal;

    auto inputSequence = [1, 2, 3, 4, 5];
    auto expectedOutput = [1, 2, 3, 4, 5];

    assert(inputSequence.longestIncreasingSubsequence.equal(expectedOutput));
}

unittest
{
    import std.algorithm : equal;

    auto inputSequence = [2, 1, 3, 4, 5];
    auto expectedOutput = [1, 3, 4, 5];

    assert(inputSequence.longestIncreasingSubsequence.equal(expectedOutput));
}

unittest
{
    import std.algorithm : equal;

    auto inputSequence = [1, 2, 3, 5, 4];
    auto expectedOutput = [1, 2, 3, 4];

    assert(inputSequence.longestIncreasingSubsequence.equal(expectedOutput));
}

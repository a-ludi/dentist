/**
    Some additional range functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.range;

import std.functional : unaryFun;
import std.meta : AliasSeq;
import std.range : ElementType, isInputRange;
import std.typecons : tuple;

/**
    This range iterates over fixed-sized chunks of size chunkSize of a source
    range. Source must be an input range. chunkSize must be greater than zero.

    See Also: `std.range.chunks`
    Returns: Range of chunks, ie. `ElementType!Source[]`.
*/
auto arrayChunks(Source)(Source range, in size_t chunkSize) if (isInputRange!Source)
{
    alias Element = ElementType!Source;

    static struct ArrayChunks
    {
        private Source range;
        private const size_t chunkSize;
        private Element[] chunk;

        this(Source range, in size_t chunkSize)
        {
            this.range = range;
            this.chunkSize = chunkSize;
            this.popFront();
        }

        void popFront()
        {
            if (range.empty)
            {
                chunk = null;

                return;
            }

            chunk = new Element[chunkSize];

            foreach (i; 0 .. chunkSize)
            {
                chunk[i] = range.front;

                if (range.empty)
                {
                    chunk = chunk[0 .. i + 1];
                    break;
                }
                else
                {
                    range.popFront();
                }
            }
        }

        @property Element[] front()
        {
            return chunk;
        }

        @property bool empty()
        {
            return chunk is null;
        }
    }

    assert(chunkSize > 0, "chunkSize must be greater than zero");

    return ArrayChunks(range, chunkSize);
}

///
unittest
{
    import std.array : array;
    import std.range : iota;

    auto chunks = iota(10).arrayChunks(2);
    assert(chunks.array == [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9]]);
}

/// Generate a tuple of tuples of chunkSize.
template chunks(size_t chunkSize)
{
    auto chunks(T...)(T args) pure nothrow @safe if (args.length >= chunkSize)
    {
        return tuple(tuple(args[0 .. chunkSize]), chunks(args[chunkSize .. $]).expand);
    }

    auto chunks(T...)(T args) pure nothrow @safe
            if (0 < args.length && args.length < chunkSize)
    {
        return tuple(tuple(args[0 .. $]));
    }

    auto chunks(T...)(T args) pure nothrow @safe if (args.length == 0)
    {
        return tuple();
    }
}

///
unittest
{
    auto c1 = chunks!2(0, 1, 2, 3, 4, 5);

    assert(c1 == tuple(tuple(0, 1), tuple(2, 3), tuple(4, 5)));

    auto c2 = chunks!3(false, "1", 2.0, 3, '4', 5);

    assert(c2 == tuple(tuple(false, "1", 2.0), tuple(3, '4', 5)));

    immutable c4 = chunks!4(false, "1", 2.0, 3, '4', 5);

    static assert(c4 == tuple(tuple(false, "1", 2.0, 3), tuple('4', 5)));
}

/// Split a list of aliases into chunks.
template Chunks(size_t chunkSize, T...)
{
    static if (T.length >= chunkSize)
    {
        alias Chunks = AliasSeq!(Chunk!(T[0 .. chunkSize]), Chunks!(chunkSize, T[chunkSize .. $]));
    }
    else static if (0 < T.length && T.length < chunkSize)
    {
        alias Chunks = AliasSeq!(Chunk!(T[0 .. $]));
    }
    else static if (T.length == 0)
    {
        alias Chunks = AliasSeq!();
    }
    else
    {
        static assert(0);
    }
}

template Chunk(T...)
{
    struct Chunk
    {
        alias chunks = T;
    }
}

///
unittest
{
    alias c1 = Chunks!(2, AliasSeq!(int, int, int, int, int, int));

    static assert(is(c1 == AliasSeq!(Chunk!(int, int), Chunk!(int, int), Chunk!(int, int))));
    static foreach (pair; c1)
    {
        static foreach (type; pair.chunks)
        {
            static assert(is(type == int));
        }
    }
}

/*
    Build a comparator according to `pred`.
*/
template Comparator(pred...) if (pred.length == 1)
{
    /// Return comparison value akin to `opCmp`.
    int compare(S, T = S)(in S a, in T b)
    {
        alias _pred = unaryFun!pred;
        auto _a = _pred(a);
        auto _b = _pred(b);

        if (_a < _b)
            return -1;
        else if (_a == _b)
            return 0;
        else
            return 1;
    }

    /// Return `true` iff `a < b`.
    bool lt(S, T = S)(in S a, in T b)
    {
        return compare!(S, T)(a, b) < 0;
    }

    /// Return `true` iff `a <= b`.
    bool le(S, T = S)(in S a, in T b)
    {
        return compare!(S, T)(a, b) <= 0;
    }

    /// Return `true` iff `a == b`.
    bool eq(S, T = S)(in S a, in T b)
    {
        return compare!(S, T)(a, b) == 0;
    }

    /// Return `true` iff `a >= b`.
    bool ge(S, T = S)(in S a, in T b)
    {
        return compare!(S, T)(a, b) >= 0;
    }

    /// Return `true` iff `a > b`.
    bool gt(S, T = S)(in S a, in T b)
    {
        return compare!(S, T)(a, b) > 0;
    }
}

///
unittest
{
    alias compareSquares = Comparator!"a ^^ 2".compare;

    assert(compareSquares(1, 2) < 0);
    assert(compareSquares(1, -2) < 0);
    assert(compareSquares(-1, 1) == 0);
    assert(compareSquares(-2.0, 1) > 0);

    alias compareByLength = Comparator!"a.length".compare;

    assert(compareByLength([], [1]) < 0);
    assert(compareByLength([1, 2], [1]) > 0);
    assert(compareByLength([1, 2], ["1", "2"]) == 0);

    alias compareAbsInts = Comparator!("a > 0 ? a : -a").compare!(int);

    assert(compareSquares(1, 2) < 0);
    assert(compareSquares(1, -2) < 0);
    assert(compareSquares(-1, 1) == 0);
    assert(compareSquares(-2, 1) > 0);
    static assert(!__traits(compiles, compareAbsInts(-2.0, 1.0)));

    alias ltSquared = Comparator!("a ^^ 2").lt;

    assert(ltSquared(1, 2));
    assert(ltSquared(1, -2));
    assert(!ltSquared(-2, -1));

    alias eqSquared = Comparator!("a ^^ 2").eq;

    assert(eqSquared(1, 1));
    assert(eqSquared(1, -1));
    assert(!eqSquared(1, 2));
}

/// Take exactly `n` element from range. Throws an exception if range has not
/// enough  elements.
///
/// Throws: Exception if range has less than `n` elements.
ElementType!R[n] takeExactly(size_t n, R)(R range) if (isInputRange!R)
{
    import std.exception : enforce;

    ElementType!R[n] result;
    size_t i = 0;

    while (!range.empty && i < n)
    {
        result[i++] = range.front;
        range.popFront();
    }

    enforce!Exception(i == n, "not enough elements");

    return result;
}

///
unittest
{
    import std.exception : assertThrown;
    import std.range : iota;

    static assert(is(typeof(iota(10).takeExactly!5) == int[5]));
    assert(iota(10).takeExactly!5 == [0, 1, 2, 3, 4]);

    assertThrown!Exception(iota(2).takeExactly!5);
}

class WrapLinesImpl(R)
{
    R output;
    size_t lineWidth;
    size_t column;

    this(R output, size_t lineWidth)
    {
        this.output = output;
        this.lineWidth = lineWidth;
    }

    void put(inout(char) c)
    {
        import std.range.primitives;

        assert(c == '\n' || column < lineWidth);

        std.range.primitives.put(output, c);
        ++column;

        if (c == '\n')
        {
            column = 0;
        }

        if (column >= lineWidth)
        {
            put('\n');
        }
    }

    void put(string chunk)
    {
        foreach (c; chunk)
        {
            put(c);
        }
    }
}

auto wrapLines(R)(R output, size_t lineWidth)
{
    return new WrapLinesImpl!R(output, lineWidth);
}

unittest
{
    import std.range.primitives : put;

    auto outputBuffer = new dchar[12];
    auto output = wrapLines(outputBuffer, 10);

    put(output, "hello world");

    assert(outputBuffer == "hello worl\nd");
}

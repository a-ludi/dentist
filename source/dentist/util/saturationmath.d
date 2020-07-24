/**
    Efficient implementation of saturation math.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.saturationmath;

import std.traits;


/// Defines +/-infinity as interpreted by this module.
template saturatedInfinity(T) if (isNumeric!T)
{
    static if (isFloatingPoint!T)
        enum saturatedInfinity = T.infinity;
    else
        enum saturatedInfinity = T.max;
}

/// ditto
template saturatedMinusInfinity(T) if (isNumeric!T)
{
    static if (isFloatingPoint!T)
        enum saturatedInfinty = -T.infinity;
    else
        enum saturatedInfinty = T.min;
}


/// Computes the result of `a + b` limited to the natural min/max value of T.
/// Additionally the min/max value of T are treated as -/+infinity,
/// respectively.
T saturatedAdd(T)(T a, T b) if (isFloatingPoint!T)
{
    return a + b;
}

/// ditto
T saturatedAdd(T)(T a, T b) if (isIntegral!T && isUnsigned!T)
{
    // Implementation according to https://web.archive.org/web/20190213215419/https://locklessinc.com/articles/sat_arithmetic/
    auto sum = a + b;

    static assert(T.max == -cast(T) (true));
    static assert(0 == -cast(T) (false));
    // We got an overflow if sum is smaller than either argument; we set the
    // result to T.max in that case
    sum |= -cast(T) (sum < a);

    return sum;
}

unittest
{
    assert(saturatedAdd(42UL, 1337UL) == 42UL + 1337UL);
    assert(saturatedAdd(size_t.max, 0UL) == size_t.max);
    assert(saturatedAdd(size_t.max, 1UL) == size_t.max);
    assert(saturatedAdd(size_t.max, size_t.max) == size_t.max);
}

/// ditto
T saturatedAdd(T)(T a, T b) if (isIntegral!T && isSigned!T)
{
    // Implementation according to https://web.archive.org/web/20190213215419/https://locklessinc.com/articles/sat_arithmetic/
    alias U = Unsigned!T;

    auto ua = cast(U) a;
    auto ub = cast(U) b;
    auto res = ua + ub;

    // Calculate overflowed result. (Don't change the sign bit of ux)
    ua = (ua >> (8*T.sizeof - 1)) + T.max;

    // Force compiler to use cmovns instruction
    if (cast(T) ((ua ^ ub) | ~(ub ^ res)) >= 0)
        res = ua;

    // Treat T.max/T.min as +/-infinity
    //
    // T.max = 0111...111
    // T.min = 1000...000
    static assert(T.min == ~T.max);

    if (
        // We need to handle infinity only if different signs
        ((a ^ b) >> (8 * T.sizeof - 1)) &&
        (a == T.max || a == T.min || b == T.max || b == T.min)
    )
    {
        // Opposite infinity values add up to zero
        if ((a ^ b) == ~(cast(T) 0))
            res = 0;
        // Either a or b is +infinity (while the other is finite)
        else if ((a | b) == ~(cast(T) 0))
            res = T.max;
        // Either a or b is -infinity (while the other is finite)
        else if ((a & b) == (cast(T) 0))
            res = T.min;
    }

    return res;
}

/// ditto
unittest
{
    assert(saturatedAdd(42L, 1337L) == 42L + 1337L);
    assert(saturatedAdd(0, -4) == -4);
    assert(saturatedAdd(long.max, 0L) == long.max);
    assert(saturatedAdd(long.max, 1L) == long.max);
    assert(saturatedAdd(long.max, long.max) == long.max);
    assert(saturatedAdd(long.min, -0L) == long.min);
    assert(saturatedAdd(long.min, -1L) == long.min);
    assert(saturatedAdd(long.min, long.min) == long.min);

    assert(saturatedAdd(long.max, -1L) == long.max);
    assert(saturatedAdd(long.min, 1L) == long.min);
    assert(saturatedAdd(-1L, long.max) == long.max);
    assert(saturatedAdd(1L, long.min) == long.min);
    assert(saturatedAdd(long.max, long.min) == 0);
    assert(saturatedAdd(long.min, long.max) == 0);
}

// timed tests on pseudo-random numbers
unittest
{
    import std.algorithm;
    import std.datetime.stopwatch;
    import std.format;
    import std.range;
    import std.random;
    import std.stdio;

    auto rnd = Random(465864470);
    // Exclude min/max edge cases; they are tested separately
    alias randomInt = () => uniform!"()"(int.min, int.max, rnd);
    enum numRounds = 100_000;
    // Pre-compute random numbers
    auto randomInts = iota(2 * numRounds)
        .map!(i => randomInt())
        .array;
    // Pre-compute expected results
    auto sums = randomInts
        .chunks(2)
        .map!(pair => clamp(cast(long) pair[0] + cast(long) pair[1], int.min, int.max))
        .array;
    size_t i;
    int x, y;
    long sum;

    alias assertCorrectResult = () =>
        assert(saturatedAdd(x, y) == sum, format!"saturatedAdd(%d, %d) != %d"(x, y, sum));

    auto result = benchmark!(
        /* both values random */ {
            x = randomInts[i++ % randomInts.length];
            y = randomInts[i++ % randomInts.length];
            sum = sums[i / 2 - 1];

            assertCorrectResult();
        },
        /* first value -infinity */ {
            x = int.min;
            y = 42;//randomInts[i++ % randomInts.length];
            sum = int.min;

            assertCorrectResult();
        },
        /* first value +infinity */ {
            x = int.max;
            y = 42;//randomInts[i++ % randomInts.length];
            sum = int.max;

            assertCorrectResult();
        },
        /* second value -infinity */ {
            x = 42;//randomInts[i++ % randomInts.length];
            y = int.min;
            sum = int.min;

            assertCorrectResult();
        },
        /* second value +infinity */ {
            x = 42;//randomInts[i++ % randomInts.length];
            y = int.max;
            sum = int.max;

            assertCorrectResult();
        },
    )(numRounds);

    debug (2)
    {
        writefln!"Computed %d rounds:"(numRounds);
        writefln!"both values random:     %fms"(result[0].total!"nsecs"/1e9*1e3);
        writefln!"first value -infinity:  %fms"(result[1].total!"nsecs"/1e9*1e3);
        writefln!"first value +infinity:  %fms"(result[2].total!"nsecs"/1e9*1e3);
        writefln!"second value -infinity: %fms"(result[3].total!"nsecs"/1e9*1e3);
        writefln!"second value +infinity: %fms"(result[4].total!"nsecs"/1e9*1e3);
    }
}

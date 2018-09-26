/**
    Some additional string functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.string;

import std.algorithm :
    joiner,
    min,
    map;
import std.array :
    array,
    minimallyInitializedArray;
import std.conv : to;
import std.exception : basicExceptionCtors, enforce;
import std.functional : binaryFun;
import std.range : chain, cycle, take;
import std.range.primitives : hasSlicing;
import std.string : lineSplitter;
import std.traits : isSomeString;
import std.typecons :
    Flag,
    No,
    tuple,
    Yes;

/**
    Adds one level of indentation for a multi-line string. Adds indentSize
    spaces to each non-empty line.

    Returns: indented string
*/
S indent(S)(S str, in size_t indentSize = 4) if (isSomeString!S)
{
    enum lineSep = "\n";
    alias indentLine = (line) => chain(" ".cycle.take(line.length == 0 ? 0 : indentSize), line);

    return str[]
        .lineSplitter
        .map!indentLine
        .joiner(lineSep)
        .chain(str[$ - lineSep.length .. $] == lineSep ? lineSep : "")
        .array
        .to!S;
}

///
unittest
{
    assert("a\nb".indent == "    a\n    b");
    assert("a\nb\n\n".indent == "    a\n    b\n\n");
    assert("a\nb\n".indent(2) == "  a\n  b\n");
}

class AlignmentException : Exception
{
    ///
    mixin basicExceptionCtors;
}

/// One edit operation of the Needleman-Wunsch algorithm.
enum EditOp: byte
{
    substitution,
    deletetion,
    insertion,
}

alias score_t = uint;

/**
    Compute an alignment of `query` against `reference` using the
    Needleman-Wunsch algorithm with non-negative scores.

    Params:
        scoreFun =     calculate score for a 'substitution' at `i, j` using
                       `scoreFun(reference[i], reference[j])`
        reference =    Sequence to compare `query` against
        query =        Sequence to compare against `reference`
        indelPenalty = Penalize each indel with this value
        memoryLimit =  throw an error if the calculation would require more
                       than `memoryLimit` bytes.
    Throws: AlignmentException if the calculation would require more than
            `memoryLimit` bytes.

    See_Also: http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm
*/
auto findAlignment(
    alias scoreFun,
    S,
)(
    in S reference,
    in S query,
    in score_t indelPenalty,
    size_t memoryLimit = 2^^20,
)
{
    alias score = binaryFun!scoreFun;
    auto memoryRequired = score_t.sizeof * (
        // Space for the scoring matrix F
        ((reference.length + 1) * (query.length + 1)) +
        // Space for the edit path
        (reference.length + query.length));
    enforce!AlignmentException(
        memoryRequired <= memoryLimit,
        "memory limit exceeded; use short reference and/or query",
    );

    auto F = DPMatrix!score_t(reference.length + 1, query.length + 1);

    // Initialize scoring matrix
    foreach (i; 0 .. reference.length + 1)
        F[i, 0] = cast(score_t) i * indelPenalty;
    foreach (j; 0 .. query.length + 1)
        F[0, j] = cast(score_t) j * indelPenalty;

    // Compute scoring matrix by rows in a cache friendly manner
    foreach (i; 1 .. reference.length + 1)
        foreach (j; 1 .. query.length + 1)
            F[i, j] = min(
                F[i - 1, j - 1] + score(reference[i - 1], query[j - 1]),
                F[i - 1, j    ] + indelPenalty,
                F[i    , j - 1] + indelPenalty,
            );

    // Find the best alignment, ie. its edit path and score
    return tuple!(
        "score",
        "editPath",
    )(
        F[$ - 1, $ - 1],
        tracebackScoringMatrix(F),
    );
}

///
unittest
{
    alias score = (x, y) => x == y ? 0 : 1;
    auto alignment = findAlignment!score("ACGTC", "AGTC", 1);

        //   - A G T C
        // - 0 1 2 3 4
        // A 1 0 1 2 3
        // C 2 1 1 2 2
        // G 3 2 1 2 3
        // T 4 3 2 1 2
        // C 5 4 3 2 1
        //
        // ACGTC
        // | |||
        // A-GTC

    assert(alignment.score == 1);
    assert(alignment.editPath == [
        EditOp.substitution,
        EditOp.deletetion,
        EditOp.substitution,
        EditOp.substitution,
        EditOp.substitution,
    ]);
}

///
unittest
{
    alias score = (x, y) => x == y ? 0 : 1;
    auto alignment = findAlignment!score("GCATGCT", "GATTACA", 1);

    //   - G A T T A C A
    // - 0 1 2 3 4 5 6 7
    // G 1 0 1 2 3 4 5 6
    // C 2 1 2 3 4 5 4 5
    // A 3 2 1 2 3 4 5 4
    // T 4 3 2 1 2 3 6 5
    // G 5 4 3 2 3 4 5 6
    // C 6 5 4 3 4 5 4 5
    // T 7 6 5 4 3 6 5 6
    //
    // GCAT-GCT
    // | || *|*
    // G-ATTACA
    assert(alignment.score == 4);
    assert(alignment.editPath == [
        EditOp.substitution,
        EditOp.deletetion,
        EditOp.substitution,
        EditOp.substitution,
        EditOp.insertion,
        EditOp.substitution,
        EditOp.substitution,
        EditOp.substitution,
    ]);
}

// Find edit path of the best alignment
private EditOp[] tracebackScoringMatrix(in DPMatrix!score_t F)
{
    auto editPath = minimallyInitializedArray!(EditOp[])(F.size[0] + F.size[1]);
    size_t i = F.size[0] - 1;
    size_t j = F.size[1] - 1;
    size_t k = editPath.length;

    while (0 < i && 0 < j)
    {
        auto matchScore     = F[i - 1, j - 1];
        auto insertionScore = F[i    , j - 1];
        auto deletionScore  = F[i - 1, j    ];
        auto nextScore = min(
            matchScore,
            deletionScore,
            insertionScore,
        );

        assert(k > 0);
        switch (nextScore)
        {
        case matchScore:
            editPath[--k] = EditOp.substitution;
            --i;
            --j;
            break;
        case insertionScore:
            editPath[--k] = EditOp.insertion;
            --j;
            break;
        case deletionScore:
            editPath[--k] = EditOp.deletetion;
            --i;
            break;
        default:
            assert(0, "corrupted scoring matrix");
        }
    }

    return editPath[k .. $];
}

unittest
{
    auto F = DPMatrix!score_t(6, 5);
    F.elements = [
        // -  A  G  T  C
           0, 1, 2, 3, 4, // -
           1, 0, 1, 2, 3, // A
           2, 1, 1, 2, 2, // C
           3, 2, 1, 2, 3, // G
           4, 3, 2, 1, 2, // T
           5, 4, 3, 2, 1, // C
    ];

    assert(F.tracebackScoringMatrix == [
        EditOp.substitution,
        EditOp.deletetion,
        EditOp.substitution,
        EditOp.substitution,
        EditOp.substitution,
    ]);
}

private struct DPMatrix(T)
{
    size_t[2] size;
    T[] elements;

    this(in size_t n, in size_t m)
    {
        this.size[0] = n;
        this.size[1] = m;
        this.elements = minimallyInitializedArray!(T[])(n * m);
    }

    unittest
    {
        auto M = DPMatrix!int(2, 3);

        assert(M.size[0] == 2);
        assert(M.size[1] == 3);

        foreach (i; 0 .. M.size[0])
            foreach (j; 0 .. M.size[1])
                M[i, j] = cast(int) (i + j);

        assert(M.elements == [
            0, 1, 2,
            1, 2, 3,
        ]);
    }

    size_t opDollar(size_t dim)() const pure nothrow if (dim < 2)
    {
        return size[dim];
    }

    auto ref T opIndex(size_t i, size_t j) pure nothrow
    {
        assert(i < size[0] && j < size[1], "index out of bounds");

        return elements[i * size[1] + j];
    }

    const(T) opIndex(size_t i, size_t j) const pure nothrow
    {
        assert(i < size[0] && j < size[1], "index out of bounds");

        return elements[i * size[1] + j];
    }

    unittest
    {
        auto M = DPMatrix!int(2, 3);
        M.elements = [
            0, 1, 2,
            1, 2, 3,
        ];

        assert(M[0, 0] == 0 + 0);
        assert(M[0, 1] == 0 + 1);
        assert(M[0, 2] == 0 + 2);
        assert(M[1, 0] == 1 + 0);
        assert(M[1, 1] == 1 + 1);
        assert(M[1, 2] == 1 + 2);
    }

    int opApply(scope int delegate(size_t i, size_t j, ref T) yield)
    {
        mixin(opApplyImpl!(No.reverse));
    }

    int opApplyReverse(scope int delegate(size_t i, size_t j, ref T) yield)
    {
        mixin(opApplyImpl!(Yes.reverse));
    }

    unittest
    {
        auto M = DPMatrix!int(2, 3);

        foreach (i, j, ref elem; M)
            elem = cast(int) (i + j);

        assert(M.elements == [
            0, 1, 2,
            1, 2, 3,
        ]);
    }

    int opApply(scope int delegate(size_t i, size_t j, in T) yield) const
    {
        mixin(opApplyImpl!(No.reverse));
    }

    int opApplyReverse(scope int delegate(size_t i, size_t j, in T) yield) const
    {
        mixin(opApplyImpl!(Yes.reverse));
    }

    unittest
    {
        auto M = DPMatrix!int(2, 3);
        M.elements = [
            0, 1, 2,
            1, 2, 3,
        ];

        foreach (i, j, elem; cast(const(typeof(M))) M)
        {
            assert(M[i, j] == elem);
            assert(elem == i + j);
        }
    }

    private static enum opApplyImpl(Flag!"reverse" reverse) = `
        int result = 0;

        foreach` ~ (reverse ? `_reverse` : ``) ~ ` (i, ref element; elements)
        {
            result = yield(i / size[1], i % size[1], element);

            if (result)
                break;
        }

        return result;
    `;
}

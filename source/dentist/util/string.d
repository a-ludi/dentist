/**
    Some additional string functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.string;

import std.algorithm :
    countUntil,
    joiner,
    min,
    map;
import std.array :
    appender,
    array,
    minimallyInitializedArray;
import std.conv : to;
import std.exception : basicExceptionCtors, enforce;
import std.functional : binaryFun;
import std.range :
    chain,
    chunks,
    cycle,
    only,
    take,
    zip;
import std.range.primitives : hasSlicing;
import std.string : join, lineSplitter;
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

/// Represents an alignment of two sequences.
struct SequenceAlignment(S, alias scoreFun = "a == b ? 0 : 1")
{
    alias getScore = binaryFun!scoreFun;

    score_t score;
    EditOp[] editPath;
    S reference;
    S query;
    score_t indelPenalty;
    Flag!"freeShift" freeShift;

    bool isValid() const pure nothrow
    {
        score_t computedScore;
        size_t i, j;
        foreach (k, editOp; editPath)
        {
            if (i > reference.length || j > query.length)
                return false;

            final switch (editOp)
            {
            case EditOp.substitution:
                if (i == reference.length || j == query.length)
                    return false;
                computedScore += getScore(reference[i], query[j]);
                ++i;
                ++j;
                break;
            case EditOp.deletetion:
                if (i == reference.length)
                    return false;
                computedScore += indelPenalty;
                ++i;
                break;
            case EditOp.insertion:
                if (j == query.length)
                    return false;
                computedScore += indelPenalty;
                ++j;
                break;
            }
        }
        if (freeShift)
            computedScore -= editPath.countUntil(EditOp.substitution);

        return
            i == reference.length &&
            j == query.length &&
            computedScore == score;
    }

    /**
        Get a partial alignment with respect to `reference`.
    */
    auto partial(in size_t begin, in size_t end) inout pure nothrow
    {
        assert(0 <= begin && begin <= end && end <= reference.length, "index out of bounds");

        if (end == begin)
            return typeof(this)(0, [], reference, query, indelPenalty);
        else if (0 == begin && end == reference.length)
            return this;

        score_t newScore;
        size_t editBegin, editEnd;
        size_t queryBegin, queryEnd;
        size_t i, j;
        foreach (k, editOp; editPath)
        {
            if (i == begin)
            {
                newScore = 0;
                editBegin = k;
                queryBegin = j;
            }

            final switch (editOp)
            {
            case EditOp.substitution:
                newScore += getScore(reference[i], query[j]);
                ++i;
                ++j;
                break;
            case EditOp.deletetion:
                newScore += indelPenalty;
                ++i;
                break;
            case EditOp.insertion:
                newScore += indelPenalty;
                ++j;
                break;
            }

            if (i >= end)
            {
                editEnd = k + 1;
                queryEnd = j;
                break;
            }
        }

        auto partialAlignment = typeof(this)(
            newScore,
            editPath[editBegin .. editEnd],
            reference[begin .. end],
            query[queryBegin .. queryEnd],
            indelPenalty,
            freeShift,
        );
        assert(partialAlignment.isValid());

        return partialAlignment;
    }

    auto opDollar() const pure nothrow
    {
        return reference.length;
    }

    /// ditto
    auto opIndex(in size_t[2] slice) inout pure nothrow
    {
        return partial(slice[0], slice[1]);
    }

    ///
    unittest
    {
        enum indelPenalty = 1;
        auto alignment = findAlignment("GCATGCT", "GATTACA", indelPenalty);

        assert(alignment.score == 4);
        assert(alignment.toString ==
            "GCAT-GCT\n" ~
            "| || *|*\n" ~
            "G-ATTACA");

        auto partialAlignment = alignment[1 .. 5];

        assert(partialAlignment.score == 3);
        assert(partialAlignment.toString ==
            "CAT-G\n" ~
            " || *\n" ~
            "-ATTA");
    }

    size_t[2] opSlice(size_t dim)(in size_t begin, in size_t end) const pure nothrow
    {
        return [begin, end];
    }

    /// Get a string representation of this alignment. Visual alignment breaks
    /// unless elements of the sequences convert to single chars via `to!string`.
    string toString(
        alias matchSymbol = '|',
        alias substitutionSymbol = '*',
        alias indelSymbol = ' ',
        alias gapSymbol = '-',
    )(in size_t width = 0) const pure
    {
        auto referenceLine = appender!string;
        referenceLine.reserve(editPath.length);
        auto compareLine = appender!string;
        compareLine.reserve(editPath.length);
        auto queryLine = appender!string;
        queryLine.reserve(editPath.length);

        size_t i, j;
        foreach (editOp; editPath)
        {
            final switch (editOp)
            {
            case EditOp.substitution:
                referenceLine ~= reference[i].to!string;
                compareLine ~= reference[i] == query[j]
                    ? matchSymbol
                    : substitutionSymbol;
                queryLine ~= query[j].to!string;
                ++i;
                ++j;
                break;
            case EditOp.deletetion:
                referenceLine ~= reference[i].to!string;
                compareLine ~= indelSymbol;
                queryLine ~= gapSymbol;
                ++i;
                break;
            case EditOp.insertion:
                referenceLine ~= gapSymbol;
                compareLine ~= indelSymbol;
                queryLine ~= query[j].to!string;
                ++j;
                break;
            }
        }

        if (width == 0)
            return only(
                referenceLine.data.to!string,
                compareLine.data.to!string,
                queryLine.data.to!string,
            ).join('\n');
        else
            return zip(
                referenceLine.data.to!string.chunks(width),
                compareLine.data.to!string.chunks(width),
                queryLine.data.to!string.chunks(width),
            )
                .map!(lineTriple => only(lineTriple.expand)
                    .joiner("\n"))
                .joiner("\n\n")
                .to!string;
    }

    ///
    unittest
    {
        auto alignment = findAlignment("ACGTC", "AGTC", 1);
        assert(alignment.toString ==
            "ACGTC\n" ~
            "| |||\n" ~
            "A-GTC");
    }

    ///
    unittest
    {
        auto alignment = findAlignment("GCATGCT", "GATTACA", 1);

        assert(alignment.toString ==
            "GCAT-GCT\n" ~
            "| || *|*\n" ~
            "G-ATTACA");
    }
}

/**
    Compute an alignment of `query` against `reference` using the
    Needleman-Wunsch algorithm with non-negative scores and constant
    indel penalty. Optionally, the `freeShift` mode may be activated
    as to allow large indels at the beginning and end of the alignment.

    **Implementation Notes:** The current implementation needs
    `O(reference.length * query.length)` in time and memory. As the
    memory requirement easily exceeds available memory it can be
    limited for now. This may change in future and an implementation
    using `O(max(reference.length, query.length))` memory will be
    silently selected for large inputs.

    Params:
        scoreFun =     calculate score for a 'substitution' at `i, j` using
                       `scoreFun(reference[i], reference[j])`
        reference =    Sequence to compare `query` against
        query =        Sequence to compare against `reference`
        indelPenalty = Penalize each indel with this value
        freeShift =    Allow indels at the beginning and end of the alignment
        memoryLimit =  throw an error if the calculation would require more
                       than `memoryLimit` bytes.
    Throws: AlignmentException if the calculation would require more than
            `memoryLimit` bytes.

    See_Also: http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm
*/
SequenceAlignment!(const(S), scoreFun) findAlignment(
    alias scoreFun = "a == b ? 0 : 1",
    S,
)(
    in S reference,
    in S query,
    in score_t indelPenalty,
    Flag!"freeShift" freeShift = No.freeShift,
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
        F[i, 0] = freeShift ? 0 : cast(score_t) i * indelPenalty;
    foreach (j; 0 .. query.length + 1)
        F[0, j] = freeShift ? 0 : cast(score_t) j * indelPenalty;

    // Compute scoring matrix by rows in a cache friendly manner
    foreach (i; 1 .. reference.length + 1)
        foreach (j; 1 .. query.length + 1)
            F[i, j] = min(
                F[i - 1, j - 1] + score(reference[i - 1], query[j - 1]),
                F[i - 1, j    ] + indelPenalty,
                F[i    , j - 1] + indelPenalty,
            );

    return typeof(return)(
        F[$ - 1, $ - 1],
        tracebackScoringMatrix(F),
        reference,
        query,
        indelPenalty,
        freeShift,
    );
}

///
unittest
{
    auto alignment = findAlignment("ACGTC", "AGTC", 1);

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
    auto alignment = findAlignment("GCATGCT", "GATTACA", 1);

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

///
unittest
{
    auto alignment = findAlignment!"a == b ? 0 : 1"(
        "tgaggacagaagggtcataggtttaattctggtcacaggcacattcctgg" ~
        "gttgcaggtttgatctccacctggtcggggcacatgca",
        "tgaggacagaagggtcatggtttaattctggtcacaggcacattcctggg" ~
        "ttgtaggctcaattcccacccggtcggggccacatgca",
        1,
    );

    assert(alignment.toString(50) ==
        "tgaggacagaagggtcataggtttaattctggtcacaggcacattcctgg\n" ~
        "|||||||||||||||||| |||||||||||||||||||||||||||||||\n" ~
        "tgaggacagaagggtcat-ggtttaattctggtcacaggcacattcctgg\n" ~
        "\n" ~
        "gttgcaggtttgatctcc-acctggtcggggc-acatgca\n" ~
        "||||*|||*|**|| ||| |||*||||||||| |||||||\n" ~
        "gttgtaggctcaat-tcccacccggtcggggccacatgca");
}

///
unittest
{
    auto alignment = findAlignment(
        "tgaggacagaagggtcataggtttaattctggtcacaggcacattcctgg" ~
        "gttgcaggtttgatctccacctggtcggggcacatgca",
        "aatgaacagctgaggacagaagggtcatggtttaattctggtcacaggca" ~
        "cattcctgggttgtaggctcaattcccacccggtcggggccacatgca",
        1,
        Yes.freeShift,
    );

    assert(alignment.toString(50) ==
        "----------tgaggacagaagggtcataggtttaattctggtcacaggc\n" ~
        "          |||||||||||||||||| |||||||||||||||||||||\n" ~
        "aatgaacagctgaggacagaagggtcat-ggtttaattctggtcacaggc\n" ~
        "\n" ~
        "acattcctgggttgcaggtttgatctcc-acctggtcggggc-acatgca\n" ~
        "||||||||||||||*|||*|**|| ||| |||*||||||||| |||||||\n" ~
        "acattcctgggttgtaggctcaat-tcccacccggtcggggccacatgca");
}

///
unittest
{
    auto alignment = findAlignment(
        "tatcctcaggtgaggactaacaacaaatatatatatatttatatctaaca" ~
        "acatatgatttctaaaatttcaaaaatcttaaggctgaattaat",
        "tatcctcaggtgaggcttaacaacaaatatatatatactgtaatatctaa" ~
        "caacatatgattctaaaatttcaaaatgcttaaaggtctga",
        1,
        Yes.freeShift,
    );

    assert(alignment.toString(50) ==
        "tatcctcaggtgaggact-aacaacaaatatatatata-ttta-tatcta\n" ~
        "||||||||||||||| || ||||||||||||||||||| |*|| ||||||\n" ~
        "tatcctcaggtgagg-cttaacaacaaatatatatatactgtaatatcta\n" ~
        "\n" ~
        "acaacatatgatttctaaaatttcaaaaatcttaa-gg-ctgaattaat\n" ~
        "||||||||||||| ||||||||||||||**||||| || ||||      \n" ~
        "acaacatatgatt-ctaaaatttcaaaatgcttaaaggtctga------");
}

unittest
{
    auto alignment = findAlignment(
        "tgaggacagaagggtcataggtttaattctggtcacaggcacattcctgg" ~
        "gttgcaggtttgatctccacctggtcggggcacatgcaggaggcaaccaa" ~
        "tcaatgtgtctctttctcagtgatgtttcttctctctgtctctctccccc" ~
        "ctctcttctactctctctgaaaaatagatggaaaaaatatcctcaggtga" ~
        "ggactaacaacaaatatatatatatttatatctaacaacatatgatttct" ~
        "aaaatttcaaaaatcttaaggctgaattaat",
        "aatgaacagctgaggacagaagggtcatggtttaattctggtcacaggca" ~
        "cattcctgggttgtaggctcaattcccacccggtcggggccacatgcaga" ~
        "aggcaaccatcaatgtgtctctttcaagtgatgtttcttctctctgtcta" ~
        "ctccccctccctatctactctctgggaaacttatggaaaaaatatcctca" ~
        "ggtgaggcttaacaacaaatatatatatactgtaatatctaacaacatat" ~
        "gattctaaaatttcaaaatgcttaaaggtctga",
        1,
        Yes.freeShift,
    );

    assert(alignment.toString(50) ==
        "----------tgaggacagaagggtcataggtttaattctggtcacaggc\n" ~
        "          |||||||||||||||||| |||||||||||||||||||||\n" ~
        "aatgaacagctgaggacagaagggtcat-ggtttaattctggtcacaggc\n" ~
        "\n" ~
        "acattcctgggttgcaggtttgatctcc-acctggtcggggc-acatgca\n" ~
        "||||||||||||||*|||*|**|| ||| |||*||||||||| |||||||\n" ~
        "acattcctgggttgtaggctcaat-tcccacccggtcggggccacatgca\n" ~
        "\n" ~
        "ggaggcaaccaatcaatgtgtctctttctcagtgatgtttcttctctctg\n" ~
        "|*||||||||| |||||||||||||||| *||||||||||||||||||||\n" ~
        "gaaggcaacca-tcaatgtgtctctttc-aagtgatgtttcttctctctg\n" ~
        "\n" ~
        "tctctctcccccctctct-tctactctctctgaaaaatagatggaaaaaa\n" ~
        "||| *||||||| ||*|| ||||||||||**||||**||  |||||||||\n" ~
        "tct-actccccc-tccctatctactctctgggaaactta--tggaaaaaa\n" ~
        "\n" ~
        "tatcctcaggtgaggact-aacaacaaatatatatata-ttta-tatcta\n" ~
        "||||||||||||||| || ||||||||||||||||||| |*|| ||||||\n" ~
        "tatcctcaggtgagg-cttaacaacaaatatatatatactgtaatatcta\n" ~
        "\n" ~
        "acaacatatgatttctaaaatttcaaaaatcttaa-gg-ctgaattaat\n" ~
        "||||||||||||| ||||||||||||||**||||| || ||||      \n" ~
        "acaacatatgatt-ctaaaatttcaaaatgcttaaaggtctga------");
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

    while (0 < i)
    {
        assert(k > 0);
        editPath[--k] = EditOp.deletetion;
        --i;
    }

    while (0 < j)
    {
        assert(k > 0);
        editPath[--k] = EditOp.insertion;
        --j;
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

    this(in size_t n, in size_t m, ref T[] buffer)
    {
        this.size[0] = n;
        this.size[1] = m;
        auto numElements = n * m;
        assert(buffer.length <= numElements);
        this.elements = buffer;
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

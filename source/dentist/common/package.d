/**
    This package holds common code for the `dentist` algorithm.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common;

import dentist.util.region : Region;
import std.algorithm : map;
import std.array : array;
import std.traits : TemplateOf;

public import dentist.common.alignments;
public import dentist.common.scaffold;

/// Thrown if some runtime error in the `dentist` algorithm occurs.
class DentistException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

/// A region of the reference aka. mask.
alias ReferenceRegion = Region!(size_t, size_t, "contigId");
/// An interval of a reference contig.
alias ReferenceInterval = ReferenceRegion.TaggedInterval;
/// A point on the reference.
alias ReferencePoint = ReferenceRegion.TaggedPoint;
/// A region of a read aka. mask.
alias ReadRegion = Region!(size_t, size_t, "readId");
/// An interval of a read contig.
alias ReadInterval = ReadRegion.TaggedInterval;
/// A point on a read.
alias ReadPoint = ReadRegion.TaggedPoint;

/// Returns the alignment region of alignmentChain.
R to(R, string contig = "contigA")(in AlignmentChain alignmentChain) pure
        if (__traits(isSame, TemplateOf!R, Region))
{
    auto contigId = mixin("alignmentChain." ~ contig ~ ".id");

    return R(alignmentChain
        .localAlignments
        .map!(la => R.TaggedInterval(
            contigId,
            mixin("la." ~ contig ~ ".begin"),
            mixin("la." ~ contig ~ ".end"),
        ))
        .array
    );
}

/**
    This package holds common code for the `dentist` algorithm.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common;

import dentist.util.log;
import dentist.util.region : Region;
import std.algorithm : count, fold, map, max, sum;
import std.array : array;
import std.conv : to;
import std.format : format;
import std.math : floor;
import std.traits : TemplateOf;
import std.typecons : Flag;

public import dentist.common.alignments;
public import dentist.common.binio;
public import dentist.common.scaffold;


/// True iff building with testing commands.
version (DentistTesting)
    enum isTesting = true;
else
    enum isTesting = false;

/// Evaluate to `value` if building with testing command;
/// otherwise to `typeof(value).init`.
template testingOnly(alias value)
{
    static if (isTesting)
        enum testingOnly = value;
    else
        enum testingOnly = typeof(value).init;
}

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

/// A point on the output assembly.
struct OutputCoordinate
{
    static enum OriginType : ubyte
    {
        global,
        contig,
        scaffold,
        scaffoldContig,
    }

    id_t scaffoldId;
    id_t contigId;
    coord_t coord;

    @property coord_t idx() const pure nothrow
    {
        return coord - 1;
    }

    @property OriginType originType() const pure nothrow
    {
        if (scaffoldId == 0 && contigId == 0)
            return OriginType.global;
        else if (scaffoldId == 0 && contigId > 0)
            return OriginType.contig;
        else if (scaffoldId > 0 && contigId == 0)
            return OriginType.scaffold;
        else
            return OriginType.scaffoldContig;
    }

    string toString() const
    {
        final switch(originType)
        {
        case OriginType.global:
            return format!`%d`(coord);
        case OriginType.contig:
            return format!`contig/%d/%d`(contigId, coord);
        case OriginType.scaffold:
            return format!`scaffold/%d/%d`(scaffoldId, coord);
        case OriginType.scaffoldContig:
            return format!`scaffold/%d/contig/%d/%d`(scaffoldId, contigId, coord);
        }
    }
}

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

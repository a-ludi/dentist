/**
    This package holds common code for the DENTIST algorithm.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common;

import dentist.util.log;
import dentist.util.region : Region;
import std.algorithm :
    among,
    count,
    fold,
    map,
    max,
    sum;
import std.array : array;
import std.conv : to;
import std.exception : enforce;
import std.format : format;
import std.math : floor;
import std.traits : TemplateOf;
import std.typecons : Flag;
import vibe.data.json : Json;

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


/// Thrown if some runtime error in the DENTIST algorithm occurs.
class DentistException : Exception
{
    /// Auxiliary data describing the circumstances under which the exception
    /// occurred.
    Json payload;

    /**
        Construct a new `DentistException`.

        Params:
            msg      = The message for the exception.
            payload  = Additional information for the exception.
            file     = The file where the exception occurred.
            line     = The line number where the exception occurred.
            next     = The previous exception in the chain of exceptions, if any.
    */
    this(string msg, string file = __FILE__, size_t line = __LINE__,
         Throwable next = null) @nogc @safe pure nothrow
    {
        super(msg, file, line, next);
    }

    /// ditto
    this(string msg, Throwable next, string file = __FILE__,
         size_t line = __LINE__) @nogc @safe pure nothrow
    {
        super(msg, file, line, next);
    }

    /// ditto
    this(string msg, Json payload, string file = __FILE__, size_t line = __LINE__,
         Throwable next = null) @nogc @safe pure nothrow
    {
        super(msg, file, line, next);
        this.payload = payload;
    }

    /// ditto
    this(string msg, Json payload, Throwable next, string file = __FILE__,
         size_t line = __LINE__) @nogc @safe pure nothrow
    {
        super(msg, file, line, next);
        this.payload = payload;
    }
}


/**
    Enforces that the given value is true. If the given value is false, a
    `DentistException` is thrown.

    Params:
        value    = Test value.
        msg      = The message for the exception.
        payload  = Additional information for the exception.
        file     = The file where the exception occurred.
        line     = The line number where the exception occurred.
    See_also: `std.exception.enforce`
    Returns:  `value`, if `cast(bool) value` is true. Otherwise,
              `new DentistException(msg, payload)` is thrown.
*/
T dentistEnforce(T)(
    T value,
    lazy string msg,
    lazy Json payload = Json(),
    string file = __FILE__,
    size_t line = __LINE__,
)
{
    return enforce(value, new DentistException(
        msg,
        payload,
        file,
        line,
    ));
}

/// A region, interval or point of the reference/read.
///
/// A point is a zero-based coordinate (property `value`) on the contig/read
/// specified by `contigId`/`readId`.
///
/// An interval is a pair of zero-based coordinates describing the right-open
/// interval `[begin, end)` on the contig/read specified by
/// `contigId`/`readId`.
///
/// A region is a collection of disjunctive intervals. When creating or
/// expanding a region the intervals are normalized, i.e. merged if they
/// overlap. This is the internal representation of masks.
///
/// See_also: `dentist.util.region.Region`,
///     `dentist.util.region.Region.TaggedInterval`,
///     `dentist.util.region.Region.TaggedPoint`
alias ReferenceRegion = Region!(size_t, size_t, "contigId");
/// ditto
alias ReferenceInterval = ReferenceRegion.TaggedInterval;
/// ditto
alias ReferencePoint = ReferenceRegion.TaggedPoint;
/// ditto
alias ReadRegion = Region!(size_t, size_t, "readId");
/// ditto
alias ReadInterval = ReadRegion.TaggedInterval;
/// ditto
alias ReadPoint = ReadRegion.TaggedPoint;


/// A point on the output assembly. This is used in the `translate-coords`
/// command.
///
/// See_also: `dentist.commands.translateCoords`
struct OutputCoordinate
{
    /// Type of the coordinate origin or reference system.
    static enum OriginType : ubyte
    {
        /// Not implemented. Global base pair position counting all ACGTN's.
        global,
        /// Not implemented. Position on a specific contig (only ACGT).
        contig,
        /// Position on specific contig counting all ACGTN's.
        scaffold,
        /// Not implemented. Position on a specific contig of a specific
        /// scaffold (only ACGT).
        scaffoldContig,
    }

    /// One-based scaffold ID; a value of zero signifies absence.
    id_t scaffoldId;
    /// One-based contig ID; a value of zero signifies absence.
    id_t contigId;
    /// One-based coordinate.
    coord_t coord;


    /// Zero-based coordinate.
    @property coord_t idx() const pure nothrow
    {
        return coord - 1;
    }


    /// Return the origin type or reference system of this coordinate.
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


    /// Encode this coordinate in a string with format
    /// `[scaffold/<scaff>/][contig/<contig>/]<coord>`.
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


/// Returns the alignment region of `alignmentChain`.
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


/**
    Get the interval that the alignment covers. This method does returns a
    single interval from the first to the last matching base pair. It takes
    complementary alignments into account when an interval on `contigB` is
    requested

    Params:
        Interval  = interval type with (at least) three fields `contigId`,
                   `begin` and `end`
        contig    = either `"contigA"` or `"contigB"`
        alignment = alignment chain

    Returns: interval from the first to the last matching base pair.
*/
Interval toInterval(Interval, string contig)(in AlignmentChain alignment)
{
    static assert(contig.among("contigA", "contigB"), "invalid contig name");

    static if (contig == "contigA")
    {
        return Interval(
            alignment.contigA.id + 0,
            alignment.first.contigA.begin + 0,
            alignment.last.contigA.end + 0,
        );
    }
    else
    {
        static assert(contig == "contigB");

        if (alignment.flags.complement)
            return Interval(
                alignment.contigB.id + 0,
                alignment.contigB.length - alignment.last.contigB.end,
                alignment.contigB.length - alignment.first.contigB.begin,
            );
        else
            return Interval(
                alignment.contigB.id + 0,
                alignment.first.contigB.begin,
                alignment.last.contigB.end,
            );
    }
}

/**
    Defines a Region and common operation with these. A Region is a set of
    tagged intervals where differently tagged intervals are distinct.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.region;

import std.algorithm : all, cmp, filter, map, max, min, sort, sum;
import std.array : appender, array, join;
import std.exception : assertThrown;
import std.format : format;
import std.functional : unaryFun;
import std.range : assumeSorted, dropExactly, ElementType, isInputRange, only, retro;
import std.traits : isNumeric, Unqual;

/// Returns the type of the property `tag` of `T`.
template TagType(T)
{
    private T instance;

    alias TagType = typeof(instance.tag);
}

/// Checks if T has a property `tag` implicitly convertible to `Tag` – if given.
template isTaggable(T)
{
    enum isTaggable = is(TagType!T);
}

unittest
{
    struct Taggable
    {
        int tag;
    }

    static assert(isTaggable!Taggable);
    static assert(!isTaggable!int);
}

unittest
{
    struct Taggable
    {
        int tag;
    }

    static assert(isTaggable!Taggable);
    static assert(!isTaggable!int);
}

/// Thrown if two operands require the same tag but different were provided.
static class MismatchingTagsException(Tag) : Exception
{
    const(Tag[2]) tags;

    this(in Tag tagA, in Tag tagB)
    {
        this.tags = [tagA, tagB];

        super(format!"mismatching interval tags: %s"(this.tags));
    }
}

/*
    Throws an exception if tags do not match.

    Throws: MismatchingTagsException if tags do not match.
*/
void enforceMatchingTags(Taggable)(in Taggable taggableA, in Taggable taggableB) pure
        if (isTaggable!Taggable)
{
    alias Tag = TagType!Taggable;

    if (taggableA.tag != taggableB.tag)
    {
        throw new MismatchingTagsException!Tag(taggableA.tag, taggableB.tag);
    }
}

/// Thrown if regions is unexpectedly empty.
static class EmptyRegionException : Exception
{
    this()
    {
        super("empty region");
    }
}

/*
    Throws an exception if region is empty.

    Throws: EmptyRegionException if region is empty.
*/
void enforceNonEmpty(R)(in R region) if (is(R : Region!Args, Args...))
{
    if (region.empty)
    {
        throw new EmptyRegionException();
    }
}

/**
    A Region is a set of tagged intervals where differently tagged intervals are distinct.
*/
struct Region(Number, Tag, string tagAlias = null, Tag emptyTag = Tag.init)
{
    static assert(isNumeric!Number, "interval limits must be numeric");

    static if (__traits(compiles, Number.infinity))
    {
        static enum numberSup = Number.infinity;
    }
    else
    {
        static enum numberSup = Number.max;
    }

    /**
        This represents a single tagged point.
    */
    static struct TaggedPoint
    {
        Tag tag = emptyTag;
        Number value;
        static if (!(tagAlias is null) && tagAlias != "tag")
        {
            mixin("alias " ~ tagAlias ~ " = tag;");
        }

        bool opBinary(string op)(auto ref const TaggedPoint other) const pure nothrow
                if (op == "==")
        {
            return this.tag == other.tag &&
                   this.value == other.value;
        }

        int opCmp(in TaggedPoint other) const pure nothrow
        {
            return cmp(
                only(this.tag, this.value),
                only(other.tag, other.value),
            );
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TP = R.TaggedPoint;

            assert(TP(0, 10) > TP(0, 0));
            assert(TP(0, 10) > TP(0, 5));
            assert(TP(0, 10) < TP(0, 12));
            assert(TP(0, 10) < TP(0, 15));
            assert(TP(0, 10) == TP(0, 10));
            assert(TP(0, 10) < TP(0, 25));
            assert(TP(0, 10) < TP(1, 25));
        }

        /// Returns true iff `this` is in `interval`.
        bool opBinary(string op)(in TaggedInterval interval) const pure nothrow if (op == "in")
        {
            return this.tag == interval.tag &&
                   interval.begin <= this.value &&
                   this.value < interval.end;
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TP = R.TaggedPoint;
            alias TI = R.TaggedInterval;

            enum interval = TI(0, 10, 20);

            assert(TP(0, 10)  in interval);
            assert(TP(1, 10) !in interval);
            assert(TP(0,  0) !in interval);
            assert(TP(0, 20) !in interval);
        }
    }

    /**
        This is a right-open interval `[begin, end)` tagged with `tag`.
        If `tagAlias` is given then the tag may be access as a property of
        that name.
    */
    static struct TaggedInterval
    {
        Tag tag = emptyTag;
        Number begin;
        Number end;
        static if (!(tagAlias is null) && tagAlias != "tag")
        {
            mixin("alias " ~ tagAlias ~ " = tag;");
        }

        invariant
        {
            assert(begin <= end, "begin must be less than or equal to end");
        }

        /// Returns the size of this interval.
        @property Number size() pure const nothrow
        {
            return this.end - this.begin;
        }

        ///
        unittest
        {
            assert(Region!(int, int).TaggedInterval().size == 0);
            assert(TaggedInterval(1, 10, 20).size == 10);
            assert(TaggedInterval(2, 20, 40).size == 20);
            assert(TaggedInterval(3, 30, 60).size == 30);
        }

        /// Returns true iff the interval is empty. An interval is empty iff
        /// `begin == end`.
        @property bool empty() pure const nothrow
        {
            return begin >= end;
        }

        ///
        unittest
        {
            assert(TaggedInterval().empty);
            assert(!TaggedInterval(1, 10, 20).empty);
            assert(!TaggedInterval(2, 20, 40).empty);
            assert(TaggedInterval(3, 60, 60).empty);
        }

        bool opBinary(string op)(auto ref const TaggedInterval other) const pure nothrow
                if (op == "==")
        {
            return (this.empty && other.empty) || (
                this.tag == other.tag &&
                this.begin == other.begin &&
                this.end == other.end
            );
        }

        /**
            Returns the convex hull of the intervals.

            Throws: MismatchingTagsException if `tag`s differ.
        */
        static TaggedInterval convexHull(in TaggedInterval[] intervals...) pure
        {
            if (intervals.length == 0)
            {
                return TaggedInterval();
            }

            TaggedInterval convexHullInterval = intervals[0];

            foreach (interval; intervals[1 .. $])
            {
                enforceMatchingTags(convexHullInterval, interval);

                if (convexHullInterval.empty)
                {
                    convexHullInterval = interval;
                }
                else if (interval.empty)
                {
                    continue;
                }
                else
                {
                    convexHullInterval.begin = min(convexHullInterval.begin, interval.begin);
                    convexHullInterval.end = max(convexHullInterval.end, interval.end);
                }
            }

            return convexHullInterval.empty
                ? TaggedInterval()
                : convexHullInterval;
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(TI.convexHull(TI(0, 10, 20), TI(0, 0, 5)) == TI(0, 0, 20));
            assert(TI.convexHull(TI(0, 10, 20), TI(0, 5, 15)) == TI(0, 5, 20));
            assert(TI.convexHull(TI(0, 10, 20), TI(0, 12, 18)) == TI(0, 10, 20));
            assert(TI.convexHull(TI(0, 10, 20), TI(0, 10, 20)) == TI(0, 10, 20));
            assert(TI.convexHull(TI(0, 10, 20), TI(0, 15, 25)) == TI(0, 10, 25));
            assert(TI.convexHull(TI(0, 10, 20), TI(0, 25, 30)) == TI(0, 10, 30));
            assertThrown!(MismatchingTagsException!int)(TI.convexHull(TI(0, 10, 20), TI(1, 25, 30)));
        }

        /// Returns the intersection of both intervals; empty if `tag`s differ.
        TaggedInterval opBinary(string op)(in TaggedInterval other) const pure nothrow
                if (op == "&")
        {
            if (this.tag != other.tag)
            {
                return TaggedInterval();
            }

            auto newBegin = max(this.begin, other.begin);
            auto newEnd = min(this.end, other.end);

            if (newBegin > newEnd)
            {
                return TaggedInterval();
            }

            return TaggedInterval(
                tag,
                newBegin,
                newEnd,
            );
        }

        /// ditto
        TaggedInterval opOpAssign(string op)(in TaggedInterval other) if (op == "&")
        {
            {
                auto tmp = this & other;

                this.tag = tmp.tag;
                this.begin = tmp.begin;
                this.end = tmp.end;

                return this;
            }
        }
        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert((TI(0, 10, 20) & TI(0, 0, 5)).empty);
            assert((TI(0, 10, 20) & TI(0, 5, 15)) == TI(0, 10, 15));
            assert((TI(0, 10, 20) & TI(0, 12, 18)) == TI(0, 12, 18));
            assert((TI(0, 10, 20) & TI(0, 10, 20)) == TI(0, 10, 20));
            assert((TI(0, 10, 20) & TI(0, 15, 25)) == TI(0, 15, 20));
            assert((TI(0, 10, 20) & TI(0, 25, 30)).empty);
            assert((TI(0, 10, 20) & TI(1, 25, 30)).empty);
        }

        /// Returns the difference of both intervals.
        Region opBinary(string op)(in TaggedInterval other) const if (op == "-")
        {
            auto intersection = this & other;

            if (intersection.empty)
            {
                TaggedInterval thisCopy = this;

                return Region(this.empty ? [] : [thisCopy]);
            }

            return Region(only(
                TaggedInterval(
                    tag,
                    this.begin,
                    intersection.begin,
                ),
                TaggedInterval(
                    tag,
                    intersection.end,
                    this.end,
                ),
            ).filter!"!a.empty".array);
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(TI(0, 10, 20) - TI(0, 0, 5) == R([TI(0, 10, 20)]));
            assert(TI(0, 10, 20) - TI(0, 5, 15) == R([TI(0, 15, 20)]));
            assert(TI(0, 10, 20) - TI(0, 12, 18) == R([TI(0, 10, 12), TI(0, 18, 20)]));
            assert(TI(0, 10, 20) - TI(0, 10, 20) == R([]));
            assert(TI(0, 10, 20) - TI(0, 15, 25) == R([TI(0, 10, 15)]));
            assert(TI(0, 10, 20) - TI(0, 25, 30) == R([TI(0, 10, 20)]));
            assert(TI(0, 10, 20) - TI(1, 25, 30) == R([TI(0, 10, 20)]));
        }

        int opCmp(in TaggedInterval other) const pure nothrow
        {
            return cmp(
                only(this.tag, this.begin, this.end),
                only(other.tag, other.begin, other.end),
            );
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(TI(0, 10, 20) > TI(0, 0, 5));
            assert(TI(0, 10, 20) > TI(0, 5, 15));
            assert(TI(0, 10, 20) < TI(0, 12, 18));
            assert(TI(0, 10, 20) < TI(0, 15, 25));
            assert(TI(0, 10, 20) == TI(0, 10, 20));
            assert(TI(0, 10, 20) < TI(0, 25, 30));
            assert(TI(0, 10, 20) < TI(1, 25, 30));
        }

        /// Returns true iff the tagged intervals intersect.
        bool intersects(in TaggedInterval other) const pure nothrow
        {
            return !(this & other).empty;
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(!TI(0, 10, 20).intersects(TI(0, 0, 5)));
            assert(TI(0, 10, 20).intersects(TI(0, 5, 15)));
            assert(TI(0, 10, 20).intersects(TI(0, 12, 18)));
            assert(TI(0, 10, 20).intersects(TI(0, 15, 25)));
            assert(TI(0, 10, 20).intersects(TI(0, 10, 20)));
            assert(!TI(0, 10, 20).intersects(TI(0, 25, 30)));
            assert(!TI(0, 10, 20).intersects(TI(1, 25, 30)));
        }

        /// Returns true iff this interval contains other.
        bool contains(in TaggedInterval other) const pure nothrow
        {
            return (this & other) == other;
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(!TI(0, 10, 20).contains(TI(0, 0, 5)));
            assert(!TI(0, 10, 20).contains(TI(0, 5, 15)));
            assert(TI(0, 10, 20).contains(TI(0, 12, 18)));
            assert(!TI(0, 10, 20).contains(TI(0, 15, 25)));
            assert(TI(0, 10, 20).contains(TI(0, 10, 20)));
            assert(!TI(0, 10, 20).contains(TI(0, 25, 30)));
            assert(!TI(0, 10, 20).contains(TI(1, 25, 30)));
        }

        /// Returns true iff `this` is a subset of `other`, ie. fully included _in_.
        bool opBinary(string op)(in TaggedInterval other) const pure nothrow if (op == "in")
        {
            return (this & other) == this;
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(TI(0, 0, 5) !in TI(0, 10, 20));
            assert(TI(0, 5, 15) !in TI(0, 10, 20));
            assert(TI(0, 12, 18) in TI(0, 10, 20));
            assert(TI(0, 15, 25) !in TI(0, 10, 20));
            assert(TI(0, 10, 20) in TI(0, 10, 20));
            assert(TI(0, 25, 30) !in TI(0, 10, 20));
            assert(TI(1, 12, 18) !in TI(0, 10, 20));
        }

        /// Returns true iff the tagged intervals do not intersect and `this < other`.
        bool isStrictlyBefore(in TaggedInterval other) const pure nothrow
        {
            return !this.intersects(other) && this < other;
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(!TI(0, 10, 20).isStrictlyBefore(TI(0, 0, 5)));
            assert(!TI(0, 10, 20).isStrictlyBefore(TI(0, 5, 15)));
            assert(!TI(0, 10, 20).isStrictlyBefore(TI(0, 12, 18)));
            assert(!TI(0, 10, 20).isStrictlyBefore(TI(0, 15, 25)));
            assert(!TI(0, 10, 20).isStrictlyBefore(TI(0, 10, 20)));
            assert(TI(0, 10, 20).isStrictlyBefore(TI(0, 25, 30)));
            assert(TI(0, 10, 20).isStrictlyBefore(TI(1, 25, 30)));
        }

        /// Returns true iff the tagged intervals do not intersect and `this > other`.
        bool isStrictlyAfter(in TaggedInterval other) const pure nothrow
        {
            return !this.intersects(other) && this > other;
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(TI(0, 10, 20).isStrictlyAfter(TI(0, 0, 5)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(0, 5, 15)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(0, 12, 18)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(0, 15, 25)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(0, 10, 20)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(0, 25, 30)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(1, 25, 30)));
        }
    }

    ///
    unittest
    {
        static enum emptyTag = 42;
        alias R = Region!(int, int, "bucketId", emptyTag);
        alias TI = R.TaggedInterval;

        TI emptyInterval;

        // Default constructor produces empty interval.
        assert((emptyInterval).empty);
        assert(emptyInterval.tag == emptyTag);

        auto ti1 = TI(1, 0, 10);

        // The tag can be aliased:
        assert(ti1.tag == ti1.bucketId);

        auto ti2 = TI(1, 5, 15);

        // Tagged intervals with the same tag behave like regular intervals:
        assert((ti1 & ti2) == TI(1, 5, 10));

        auto ti3 = TI(2, 0, 10);

        // Tagged intervals with different tags are distinct:
        assert((ti1 & ti3).empty);
    }

    TaggedInterval[] _intervals;

    this(TaggedInterval[] intervals)
    {
        this._intervals = intervals;
        this.normalize();
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        auto region = R([TI(0, 15, 20), TI(0, 5, 10), TI(0, 0, 10)]);

        // Intervals get implicitly normalized.
        assert(region.intervals == [TI(0, 0, 10), TI(0, 15, 20)]);
    }

    this(TaggedInterval interval)
    {
        this._intervals = [interval];
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        R region = TI(0, 15, 20);

        assert(region.intervals == [TI(0, 15, 20)]);
    }

    this(Tag tag, Number begin, Number end)
    {
        this._intervals = [TaggedInterval(tag, begin, end)];
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        auto region = R(0, 15, 20);

        assert(region.intervals == [TI(0, 15, 20)]);
    }

    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        auto tis = [TI(0, 15, 20)];
        auto region = R(tis);
        auto regionDup = region;

        assert(region == regionDup);

        // Changing the original does not affect duplicate
        region |= R(TI(0, 25, 30));

        assert(region != regionDup);
    }

    /// Return a list of the tagged intervals in this region.
    @property const(TaggedInterval)[] intervals() const pure nothrow
    {
        return _intervals;
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        R emptyRegion;
        auto region1 = R(0, 5, 10);
        auto region2 = R([TI(0, 5, 10), TI(0, 15, 20)]);

        assert(emptyRegion.intervals == []);
        assert(region1.intervals == [TI(0, 5, 10)]);
        assert(region2.intervals == [TI(0, 5, 10), TI(0, 15, 20)]);
    }

    /// Return a list of the tagged intervals in this region.
    TaggedInterval[] releaseIntervals() pure nothrow
    {
        auto intervals = _intervals;

        this._intervals = [];

        return intervals;
    }

    /// Returns the size of this region.
    Number size() pure const nothrow
    {
        return _intervals.map!"a.size".sum;
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        R emptyRegion;
        auto region1 = R(0, 0, 10);
        auto region2 = R([TI(0, 0, 10), TI(0, 20, 30)]);
        auto region3 = R([TI(0, 0, 20), TI(0, 10, 30)]);

        assert(emptyRegion.size == 0);
        assert(region1.size == 10);
        assert(region2.size == 20);
        assert(region3.size == 30);
    }

    /// Returns true iff the region is empty.
    bool empty() pure const nothrow
    {
        return _intervals.length == 0;
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        R emptyRegion1;
        auto emptyRegion2 = R([TI(0, 0, 0), TI(0, 10, 10)]);
        auto emptyRegion3 = R([TI(0, 0, 0), TI(1, 0, 0)]);
        auto region1 = R(0, 0, 10);

        assert(emptyRegion1.empty);
        assert(emptyRegion2.empty);
        assert(emptyRegion3.empty);
        assert(!region1.empty);
    }

    protected void normalize()
    {
        if (_intervals.length == 0)
        {
            return;
        }

        _intervals.sort();
        TaggedInterval accInterval = _intervals[0];
        size_t insertIdx = 0;

        alias intervalsTouch = (lhs, rhs) => lhs.tag == rhs.tag &&
                                             (lhs.end == rhs.begin || lhs.begin == rhs.end);

        foreach (i, intervalB; _intervals[1 .. $])
        {
            if (intervalB.empty)
            {
                continue;
            }
            else if (accInterval.intersects(intervalB) || intervalsTouch(accInterval, intervalB))
            {
                // If two intervals intersect or touch their union is the same as the convex hull of both.
                accInterval = TaggedInterval.convexHull(accInterval, intervalB);
            }
            else
            {
                if (!accInterval.empty)
                {
                    _intervals[insertIdx++] = accInterval;
                }
                accInterval = intervalB;
            }
        }

        if (!accInterval.empty)
        {
            _intervals[insertIdx++] = accInterval;
        }
        _intervals.length = insertIdx;
    }

    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        auto region = R([
            TI(3, 8349, 8600),
            TI(3, 8349, 8349),
            TI(3, 8349, 8850),
            TI(3, 8349, 9100),
            TI(3, 8349, 8349),
            TI(3, 8349, 9350),
            TI(3, 8349, 9600),
            TI(3, 8349, 8349),
            TI(3, 8349, 9900),
            TI(3, 8349, 10150),
            TI(3, 8349, 8349),
            TI(3, 8349, 10400),
            TI(3, 8349, 10650),
            TI(3, 8349, 10800),
            TI(3, 8499, 10800),
            TI(3, 8749, 10800),
            TI(3, 8749, 8749),
            TI(3, 8999, 10800),
            TI(3, 9249, 10800),
            TI(3, 9549, 10800),
            TI(3, 9799, 10800),
            TI(3, 10049, 10800),
            TI(3, 10299, 10800),
            TI(3, 10549, 10800),
        ]);
        auto normalizedRegion = R([
            TI(3, 8349, 10800),
        ]);

        assert(region == normalizedRegion);
    }

    /// Computes the union of all tagged intervals.
    Region opBinary(string op)(in Region other) const if (op == "|")
    {
        return Region(this._intervals.dup ~ other._intervals.dup);
    }

    /// ditto
    Region opBinary(string op)(in TaggedInterval other) const if (op == "|")
    {
        return Region(this._intervals.dup ~ [cast(TaggedInterval) other]);
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        assert((R(0, 10, 20) | R(0, 0, 5)) == R([TI(0, 10, 20), TI(0, 0, 5)]));
        assert((R(0, 10, 20) | R(0, 5, 15)) == R(0, 5, 20));
        assert((R(0, 10, 20) | R(0, 12, 18)) == R(0, 10, 20));
        assert((R(0, 10, 20) | R(0, 10, 20)) == R(0, 10, 20));
        assert((R(0, 10, 20) | R(0, 15, 25)) == R(0, 10, 25));
        assert((R(0, 10, 20) | R(0, 25, 30)) == R([TI(0, 10, 20), TI(0, 25, 30)]));
        assert((R(0, 10, 20) | R(1, 25, 30)) == R([TI(0, 10, 20), TI(1, 25, 30)]));
    }

    /// Computes the intersection of the two regions.
    Region opBinary(string op)(in Region other) const if (op == "&")
    {
        if (this.empty || other.empty)
        {
            return Region();
        }

        auto intersectionAcc = appender!(TaggedInterval[]);
        intersectionAcc.reserve(this._intervals.length + other._intervals.length);

        size_t lhsIdx = 0;
        size_t lhsLength = this._intervals.length;
        TaggedInterval lhsInterval;
        size_t rhsIdx = 0;
        size_t rhsLength = other._intervals.length;
        TaggedInterval rhsInterval;

        while (lhsIdx < lhsLength && rhsIdx < rhsLength)
        {
            lhsInterval = this._intervals[lhsIdx];
            rhsInterval = other._intervals[rhsIdx];

            if (lhsInterval.isStrictlyBefore(rhsInterval))
            {
                ++lhsIdx;
            }
            else if (rhsInterval.isStrictlyBefore(lhsInterval))
            {
                ++rhsIdx;
            }
            else
            {
                assert(lhsInterval.intersects(rhsInterval));

                intersectionAcc ~= lhsInterval & rhsInterval;

                if (lhsIdx + 1 < lhsLength && rhsIdx + 1 < rhsLength)
                {
                    if (this._intervals[lhsIdx + 1].begin < other._intervals[rhsIdx + 1].begin)
                    {
                        ++lhsIdx;
                    }
                    else
                    {
                        ++rhsIdx;
                    }
                }
                else if (lhsIdx + 1 < lhsLength)
                {
                    ++lhsIdx;
                }
                else
                {
                    ++rhsIdx;
                }
            }
        }

        return Region(intersectionAcc.data);
    }

    /// ditto
    Region opBinary(string op)(in TaggedInterval other) const if (op == "&")
    {
        return this & Region([other]);
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        assert((R(0, 10, 20) & R(0, 0, 5)) == R([]));
        assert((R(0, 10, 20) & R(0, 5, 15)) == R(0, 10, 15));
        assert((R(0, 10, 20) & R(0, 12, 18)) == R(0, 12, 18));
        assert((R(0, 10, 20) & R(0, 10, 20)) == R(0, 10, 20));
        assert((R(0, 10, 20) & R(0, 15, 25)) == R(0, 15, 20));
        assert((R(0, 10, 20) & R(0, 25, 30)) == R([]));
        assert((R(0, 10, 20) & R(1, 25, 30)) == R([]));
        // R1:       [-------)   [-------)   [-------)
        // R2:             [-------)   [-------)   [-------)
        // R1 & R2:        [-)   [-)   [-)   [-)   [-)
        assert((R([
            TI(0, 0, 30),
            TI(0, 40, 70),
            TI(0, 80, 110),
        ]) & R([
            TI(0, 20, 50),
            TI(0, 60, 90),
            TI(0, 100, 130),
        ])) == R([
            TI(0, 20, 30),
            TI(0, 40, 50),
            TI(0, 60, 70),
            TI(0, 80, 90),
            TI(0, 100, 110),
        ]));
    }

    Region opBinary(string op)(in TaggedInterval interval) const if (op == "-")
    {
        if (interval.empty)
        {
            return Region(_intervals.dup);
        }

        auto differenceAcc = appender!(TaggedInterval[]);
        differenceAcc.reserve(_intervals.length + 1);

        auto trisection = _intervals
            .assumeSorted!"a.isStrictlyBefore(b)"
            .trisect(interval);
        auto copyHeadEnd = trisection[0].length;
        auto copyTailBegin = trisection[0].length + trisection[1].length;

        differenceAcc ~= _intervals[0 .. copyHeadEnd];
        foreach (lhsInterval; trisection[1])
        {
            assert(lhsInterval.intersects(interval));
            auto tmpDifference = lhsInterval - interval;

            differenceAcc ~= tmpDifference._intervals;
        }
        differenceAcc ~= _intervals[copyTailBegin .. $];

        return Region(differenceAcc.data);
    }

    Region opBinary(string op)(in Region other) const if (op == "-")
    {
        if (other.empty)
        {
            return Region(this._intervals.dup);
        }

        auto otherDifferenceCandidates = getDifferenceCandidates(other).assumeSorted!"a.isStrictlyBefore(b)";
        auto differenceAcc = appender!(TaggedInterval[]);
        differenceAcc.reserve(this._intervals.length + otherDifferenceCandidates.length);

        foreach (lhsInterval; this._intervals)
        {
            auto intersectingRhs = otherDifferenceCandidates.equalRange(lhsInterval);
            auto tmpDifference = Region(lhsInterval);

            foreach (rhsInterval; intersectingRhs)
            {
                if (tmpDifference.empty
                        || tmpDifference._intervals[$ - 1].isStrictlyBefore(rhsInterval))
                {
                    // Remaining rhsItervals will not intersect anymore.
                    break;
                }

                tmpDifference -= rhsInterval;
            }

            if (!tmpDifference.empty)
            {
                differenceAcc ~= tmpDifference._intervals;
            }
        }

        return Region(differenceAcc.data);
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        assert((R(0, 10, 20) - R(0, 0, 5)) == R(0, 10, 20));
        assert((R(0, 10, 20) - R(0, 5, 15)) == R(0, 15, 20));
        assert((R(0, 10, 20) - R(0, 12, 18)) == R([TI(0, 10, 12), TI(0, 18, 20)]));
        assert((R(0, 10, 20) - R(0, 10, 20)).empty);
        assert((R(0, 10, 20) - R(0, 15, 25)) == R(0, 10, 15));
        assert((R(0, 10, 20) - R(0, 25, 30)) == R(0, 10, 20));
        assert((R(0, 10, 20) - R(1, 25, 30)) == R(0, 10, 20));
    }

    private auto getDifferenceCandidates(in Region other) const pure nothrow
    {
        auto otherIntervals = other._intervals.assumeSorted!"a.isStrictlyBefore(b)";
        auto sliceBegin = otherIntervals.lowerBound(this._intervals[0]).length;
        auto sliceEnd = other._intervals.length - otherIntervals.upperBound(this._intervals[$ - 1]).length;

        if (sliceBegin < sliceEnd)
        {
            return other._intervals[sliceBegin .. sliceEnd];
        }
        else
        {
            return cast(typeof(other._intervals)) [];
        }
    }

    Region opOpAssign(string op, T)(in T other)
            if (is(T : Region) || is(T : TaggedInterval))
    {
        static if (op == "|" || op == "&" || op == "-")
        {
            mixin("auto tmp = this " ~ op ~ " other;");

            this._intervals = tmp._intervals;

            return this;
        }
        else
        {
            static assert(0, "unsupported operator: " ~ op);
        }
    }

    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        R accRegion;
        auto inputRegion1 = R([
            TI(3, 8349, 8600),
            TI(3, 8349, 8850),
            TI(3, 8349, 9100),
            TI(3, 8349, 9350),
            TI(3, 8349, 9600),
            TI(3, 8349, 9900),
            TI(3, 8349, 10150),
            TI(3, 8349, 10400),
            TI(3, 8349, 10650),
            TI(3, 8349, 10800),
            TI(3, 8499, 10800),
            TI(3, 8749, 10800),
            TI(3, 8999, 10800),
            TI(3, 9249, 10800),
            TI(3, 9549, 10800),
            TI(3, 9799, 10800),
            TI(3, 10049, 10800),
            TI(3, 10299, 10800),
            TI(3, 10549, 10800),
        ]);
        auto inputRegion2 = R([
            TI(3, 2297, 11371),
        ]);
        auto expectedResult = R([
            TI(3, 2297, 11371),
        ]);

        accRegion |= inputRegion1;

        assert(accRegion == inputRegion1);

        accRegion |= inputRegion2;

        assert(accRegion == expectedResult);
        assert((inputRegion1 | inputRegion2) == expectedResult);
    }

    /// Returns true iff point is in this region.
    bool opBinaryRight(string op)(in TaggedPoint point) const pure nothrow
            if (op == "in")
    {
        if (this.empty)
        {
            return false;
        }

        auto needle = TaggedInterval(point.tag, point.value, numberSup);
        auto sortedIntervals = intervals.assumeSorted;
        auto candidateIntervals = sortedIntervals.lowerBound(needle).retro;

        if (candidateIntervals.empty)
        {
            return false;
        }

        auto candidateInterval = candidateIntervals.front;

        return candidateInterval.begin <= point.value && point.value < candidateInterval.end;
    }

    /// ditto
    alias includes = opBinaryRight!"in";

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;
        alias TP = R.TaggedPoint;

        R emptyRegion;
        auto region = R([TI(0, 0, 10), TI(1, 0, 10)]);

        assert(TP(0, 0) !in emptyRegion);
        assert(TP(0, 5) !in emptyRegion);
        assert(TP(0, 10) !in emptyRegion);
        assert(TP(0, 20) !in emptyRegion);
        assert(TP(0, 0) in region);
        assert(TP(0, 5) in region);
        assert(TP(0, 10) !in region);
        assert(TP(0, 20) !in region);
        assert(TP(1, 0) in region);
        assert(TP(1, 5) in region);
        assert(TP(1, 10) !in region);
        assert(TP(1, 20) !in region);
    }
}

unittest
{
    Region!(int, int) r;
}

/**
    Returns true iff `thing` is empty

    See_Also: Region.empty, Region.TaggedInterval.empty
*/
bool empty(T)(in T thing) pure nothrow
        if (is(T : Region!Args, Args...) || is(T : Region!Args.TaggedInterval, Args...))
{
    return thing.empty;
}

///
unittest
{
    alias R = Region!(int, int);
    alias TI = R.TaggedInterval;

    R emptyRegion;
    TI emptyTI;
    auto region = R(0, 0, 10);
    auto ti = TI(0, 0, 10);

    assert(empty(emptyRegion));
    assert(empty(emptyTI));
    assert(!empty(region));
    assert(!empty(ti));
}

/**
    Returns the union of all elements.

    See_Also: Region.opBinary!"|", Region.TaggedInterval.opBinary!"|"
*/
auto union_(Range)(Range regions)
        if (isInputRange!Range && is(ElementType!Range : Region!Args, Args...))
{
    alias Region = Unqual!(ElementType!Range);

    return Region(regions
        .map!"a._intervals.dup"
        .join);
}

///
unittest
{
    alias R = Region!(int, int);
    alias TI = R.TaggedInterval;

    R emptyRegion;
    TI emptyTI;
    auto region = R(0, 0, 10);
    auto ti = TI(0, 0, 10);

    assert(empty(emptyRegion));
    assert(empty(emptyTI));
    assert(!empty(region));
    assert(!empty(ti));
}

/**
    Returns the minimum/supremum point or convex hull of the intervals. Both
    minimum and supremum are undefined for empty regions but the convex hull
    is not.

    Throws: MismatchingTagsException if `tag`s differ.
    Throws: EmptyRegionException if region is empty.
*/
auto min(R)(R region) if (is(R : Region!Args, Args...))
{
    alias TI = R.TaggedInterval;

    enforceNonEmpty(region);
    auto convexHull = TI.convexHull(region.intervals[0], region.intervals[$ - 1]);

    return convexHull.begin;
}

/// ditto
auto sup(R)(R region) if (is(R : Region!Args, Args...))
{
    alias TI = R.TaggedInterval;

    enforceNonEmpty(region);
    auto convexHull = TI.convexHull(region.intervals[0], region.intervals[$ - 1]);

    return convexHull.end;
}

/// ditto
auto convexHull(R)(R region) if (is(R : Region!Args, Args...))
{
    alias TI = R.TaggedInterval;

    if (region.empty)
        return TI();

    auto convexHull = TI.convexHull(region.intervals[0], region.intervals[$ - 1]);

    return convexHull;
}

///
unittest
{
    alias R = Region!(int, int);
    alias TI = R.TaggedInterval;

    R emptyRegion;
    auto region1 = R([TI(0, 0, 10), TI(0, 20, 30)]);
    auto region2 = R([TI(0, 0, 10), TI(1, 0, 10)]);

    assert(min(region1) == 0);
    assert(sup(region1) == 30);
    assertThrown!EmptyRegionException(min(emptyRegion));
    assertThrown!EmptyRegionException(sup(emptyRegion));
    assertThrown!(MismatchingTagsException!int)(min(region2));
    assertThrown!(MismatchingTagsException!int)(sup(region2));
}

struct Tiling(R, N, T)
{
    R region;
    N totalOverlap;
    T[] elements;
}

auto findTilings(alias toInterval, T, N)(T[] elements, in N maxLocalOverlap, in N maxGlobalOverlap = N.max)
{
    alias interval = unaryFun!toInterval;
    alias Region = typeof(interval(elements[0]) - interval(elements[0]));
    alias Num = typeof(interval(elements[0]).size);
    alias overlapEachOther = (lhs, rhs) => interval(lhs).intersects(interval(rhs));

    auto tilingsAcc = appender!(Tiling!(Region, Num, T)[]);

    void findAllowedTilings(T[] candidates, Region region, T[] included, in Num globalOverlap)
    {
        if (globalOverlap > maxGlobalOverlap)
            return;

        tilingsAcc ~= Tiling!(Region, Num, T)(
            region,
            globalOverlap,
            included,
        );

        size_t numExpansions;
        foreach (i, candidate; candidates)
        {
            auto candidateInterval = interval(candidate);
            auto newRegion = region | candidateInterval;
            auto localOverlap = region.size + candidateInterval.size - newRegion.size;

            if (localOverlap <= maxLocalOverlap)
            {
                findAllowedTilings(
                    candidates.dropExactly(i + 1),
                    newRegion,
                    included ~ [candidate],
                    globalOverlap + localOverlap,
                );
                ++numExpansions;
            }
        }
    }

    findAllowedTilings(
        elements,
        Region(),
        [],
        Num.init,
    );

    assert(tilingsAcc.data.length > 0);

    return tilingsAcc.data;
}

unittest
{
    alias R = Region!(int, int);
    alias TI = R.TaggedInterval;
    auto elements = [
        [1, 5],
        [3, 9],
        [7, 10],
        [1, 10],
    ];
    auto tilings = findTilings!(
        interval => TI(0, interval[0], interval[1])
    )(
        elements,
        2,
        4,
    );

    import std.stdio;

    assert(tilings == [
        Tiling!(R, int, int[])(
            R([]),
            0,
            [],
        ),
        Tiling!(R, int, int[])(
            R([TI(0, 1, 5)]),
            0,
            [[1, 5]],
        ),
        Tiling!(R, int, int[])(
            R([TI(0, 1, 9)]),
            2,
            [[1, 5], [3, 9]],
        ),
        Tiling!(R, int, int[])(
            R([TI(0, 1, 10)]),
            4,
            [[1, 5], [3, 9], [7, 10]],
        ),
        Tiling!(R, int, int[])(
            R([TI(0, 1, 5), TI(0, 7, 10)]),
            0,
            [[1, 5], [7, 10]],
        ),
        Tiling!(R, int, int[])(
            R([TI(0, 3, 9)]),
            0,
            [[3, 9]],
        ),
        Tiling!(R, int, int[])(
            R([TI(0, 3, 10)]),
            2,
            [[3, 9], [7, 10]],
        ),
        Tiling!(R, int, int[])(
            R([TI(0, 7, 10)]),
            0,
            [[7, 10]],
        ),
        Tiling!(R, int, int[])(
            R([TI(0, 1, 10)]),
            0,
            [[1, 10]],
        ),
]);
}

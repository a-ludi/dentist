/**
    Collection of common structure and functions for working with binary data.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.binio.common;

import dentist.util.log;
import dentist.util.math : ceil, ceildiv;
import std.algorithm :
    all,
    among,
    joiner,
    map,
    permutations;
import std.array : minimallyInitializedArray;
import std.bitmanip : bitfields;
import std.conv : to;
import std.exception : enforce, ErrnoException;
import std.format : format;
import std.process : environment;
import std.range :
    chunks,
    dropExactly,
    empty,
    enumerate,
    front,
    popFront,
    retro,
    takeExactly;
import std.stdio : File;
import std.string : toLower;
import std.traits :
    isArray,
    isSomeChar,
    isSomeString,
    Unqual;
import std.typecons : BitFlags, Flag, No, Yes;
import std.utf : UTFException;


/// Place a lock on `file` if the underlying file system supports it. This
/// ignores and logs errors encountered while trying to acquire the lock.
/// File locking may disabled entirely by setting the environment variable
/// `SKIP_FILE_LOCKING=1`.
///
/// Environment:
/// $(UL
///     $(LI `SKIP_FILE_LOCKING`: disable file locking by setting to `1`)
/// )
void lockIfPossible(ref File file) nothrow
{
    enum skipFileLocking = "SKIP_FILE_LOCKING";

    try
    {
        if (environment.get(skipFileLocking, null) == "1")
        {
            logJsonInfo(
                "info", "skipping file locking on user's request",
                "file", file.name,
            );

            return;
        }
    }
    catch (UTFException e)
    {
        logJsonError(
            "info", "ignoring environment variable " ~ skipFileLocking ~ ": " ~
                    "contains invalid UTF characters",
            "file", file.name,
            "error", e.message().to!string,
        );
    }
    catch (Exception e)
    {
        logJsonError(
            "info", "ignoring environment variable " ~ skipFileLocking ~ ": " ~
                    e.message().to!string,
            "file", file.name,
            "error", e.message().to!string,
        );
    }

    try
    {
        file.lock();
    }
    catch (Exception e)
    {
        logJsonDebug(
            "info", "could not lock file `%s`",
            "file", file.name,
            "error", e.message().to!string,
        );
    }
}


/// Thrown if an error occurs during binary I/O operations.
class BinaryIOException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}


/// Provides array access to the fields defined by `fieldPtr!T`. The latter
/// method must be implemented by the struct where this is mixed in.
mixin template DbIndex(T...)
{
    static struct EOF;


    /// Return the pointer to the begin of the `T`s array.
    @property auto ref size_t beginPtr(T)() pure nothrow
    {
        return fieldPtr!T;
    }

    /// ditto
    @property size_t beginPtr(T)() const pure nothrow
    {
        return fieldPtr!T;
    }


    /// Return the pointer to the end of the `T`s array.
    @property auto ref size_t endPtr(T)() pure nothrow
    {
        return fieldPtr!(NextType!T);
    }

    /// ditto
    @property size_t endPtr(T)() const pure nothrow
    {
        return fieldPtr!(NextType!T);
    }


    /// Return an `ArrayStorage` constructed from `beginPtr!T` and `endPtr!T`.
    @property ArrayStorage!(StorageType!T) arrayStorage(T)() const pure nothrow
    {
        return typeof(return).fromPtrs(beginPtr!T, endPtr!T);
    }
}


/// Read a single record of type `T` at offset `ptr` in `dbFile`.
T readRecordAt(T)(File dbFile, size_t ptr)
{
    dbFile.seek(ptr.to!long);

    return readRecord!T(dbFile);
}


/// Read a single record of type `T` at the current offset in `dbFile`.
T readRecord(T)(File dbFile)
{
    return readRecords(dbFile, new T[1])[0];
}


/// Read an array of records at the current offset in `dbFile`.
///
/// Throws:
///     `BinaryIOException` if an `ErrnoException` occurs or the file ends
///         before all records are read.
T readRecords(T)(File dbFile, T records) if (isArray!T)
{
    if (records.length == 0)
        return records;

    auto expectedLength = records.length;
    try
    {
        records = dbFile.rawRead(records);

        enforce!BinaryIOException(
            records.length == expectedLength,
            format!"malformed DB `%s`: premature end of file"(
                    dbFile.name)
        );
    }
    catch (ErrnoException e)
    {
        throw new BinaryIOException(
            format!"malformed DB `%s`: cannot read record of type %s: %s"(
                    dbFile.name, T.stringof, e.msg),
        );
    }

    return records;
}


/// An array-like structure of pointers that is used to access the data in a
/// binary file. In this context, pointers are expected to be offsets in the
/// binary file.
struct ArrayStorage(T)
{
    /// Number of bytes consumed by a single `T`.
    enum elementSize = T.sizeof;
    /// Base pointer to the first array element.
    size_t ptr;
    /// Number of `T` items in this array.
    size_t length;


    /// Return a shallow copy if `this`.
    ArrayStorage!T opIndex() const pure nothrow
    {
        return this;
    }


    /// Return the pointer to element `0 <= i <= length`. Contrary to regular
    /// arrays, the element at `i == length` can be accessed. This is useful
    /// the construct slices (see `opIndex(size_t[2])`).
    size_t opIndex(size_t i) const pure nothrow
    {
        assert(i <= length, "out of bounds");

        return ptr + i * elementSize;
    }


    /// Return a slice of this `ArrayStorage`. An empty slice is not allowed.
    ArrayStorage!T opIndex(size_t[2] slice) const pure nothrow
    {
        auto from = slice[0];
        auto to = slice[1];
        assert(from < to && to <= length, "out of bounds");

        return typeof(return)(this[from], to - from);
    }

    size_t[2] opSlice(size_t dim)(size_t from, size_t to) const pure nothrow
        if (dim == 0)
    {
        assert(from < to, "slice `[from .. to]` must have `from < to`");

        return [from, to];
    }


    /// Returns `length`.
    size_t opDollar() const pure nothrow
    {
        return length;
    }


    /// Computes the index corresponding to `ptr` relative to this array.
    /// The result may be out of bounds. `ptr` must have a valid alignment
    /// relative to `this.ptr`.
    long indexOf(size_t ptr) const pure
    {
        assert((ptr.to!long - this.ptr.to!long) % elementSize.to!long == 0, "bad pointer alignment");

        return (ptr.to!long - this.ptr.to!long) / elementSize.to!long;
    }


    /// Construct an `ArrayStorage` from given pointers. The pointers must
    /// be properly aligned.
    static ArrayStorage!T fromPtrs(size_t fromPtr, size_t toPtr)
    {
        assert((toPtr - fromPtr) % elementSize == 0, "bad pointer alignment");

        return typeof(return)(fromPtr, (toPtr - fromPtr) / elementSize);
    }
}

unittest
{
    import core.exception : AssertError;
    import std.exception : assertThrown;

    alias T = int;
    static assert(ArrayStorage!T.elementSize == T.sizeof);
    enum basePtr = 1337;
    enum length = 42;

    ArrayStorage!T storage;
    storage.ptr = basePtr;
    storage.length = length;

    assert(storage[0] == basePtr);
    assert(storage[13] == basePtr + 13 * T.sizeof);
    assert(storage[$] == storage[length]);
    assertThrown!AssertError(storage[length + 1]);

    auto storageCopy = storage[];
    assert(storageCopy == storage);
    --storageCopy.length;
    assert(storageCopy.length != storage.length);

    storageCopy = storage[0 .. $];
    assert(storageCopy == storage);
    --storageCopy.length;
    assert(storageCopy.length != storage.length);

    auto storageSlice = storage[1 .. 13];
    assert(storageSlice.ptr == storage[1]);
    assert(storageSlice.length == 12);

    assert(ArrayStorage!T.fromPtrs(storage[1], storage[13]) == storageSlice);

    assert(storage.indexOf(basePtr + 13 * T.sizeof) == 13);
    assert(storage.indexOf(basePtr - 13 * T.sizeof) == -13);
    assertThrown!AssertError(storage.indexOf(basePtr + 1));
}


/// 2-bit encoding of bases `a`, `c`, `t`, `g`.
enum CompressedBase : ubyte
{
    a = 0b00,
    c = 0b01,
    t = 0b10,
    g = 0b11,
}

/// Four `CompressedBase`s packed into one byte. Provides array-like access
/// to these four bases.
struct CompressedBaseQuad
{
    private
    {
        mixin(bitfields!(
            CompressedBase, "_0", 2,
            CompressedBase, "_1", 2,
            CompressedBase, "_2", 2,
            CompressedBase, "_3", 2,
        ));
    }

    /// Number of bases – always four.
    enum length = 4;

    /// ditto
    alias opDollar = length;


    /// Construct a dynamic array of the bases. This requires memory
    /// allocation and takes up the size of the array struct and additional
    /// four bytes – one per base.
    CompressedBase[] opIndex() const pure nothrow
    {
        return [_0, _1, _2, _3];
    }


    /// Access base at `i`.
    CompressedBase opIndex(size_t i) const pure nothrow
    {
        switch (i)
        {
        case 0:
            return _0;
        case 1:
            return _1;
        case 2:
            return _2;
        case 3:
            return _3;
        default:
            assert(0, "index out of bounds: " ~ i.to!string);
        }
    }

    /// Set base at `i` to `base`.
    void opIndexAssign(CompressedBase base, size_t i) pure nothrow
    {
        switch (i)
        {
        case 0:
            _0 = base;
            return;
        case 1:
            _1 = base;
            return;
        case 2:
            _2 = base;
            return;
        case 3:
            _3 = base;
            return;
        default:
            assert(0, "index out of bounds: " ~ i.to!string);
        }
    }
}


/// A string of bases encoded in `CompressedBaseQuad`s. The sequence may have
/// any length – not just multiples of four.
struct CompressedSequence
{
    // The sequence is stored in quads of bases but may have a length that is
    // not a multiple of four. The actual sequence is a slice of the underlying
    // quad buffer defined by an `_baseOffset` (in bases) and the number of
    // bases `numBases`.
    //
    // The `_baseOffset` is actually required to enable slicing without copy.

    private CompressedBaseQuad[] _data;
    private ubyte _baseOffset;
    private size_t numBases;


    /// Compress ASCII-encoded string of `acgt`s.
    ///
    /// Throws: `Exception` if `fastaString` contains characters other than
    ///     `acgt`.
    static CompressedSequence from(S)(S fastaString) if (isSomeString!S)
    {
        auto compressedLength = ceil(fastaString.length, 4) / 4;
        auto data = minimallyInitializedArray!(CompressedBaseQuad[])(compressedLength);

        foreach (i, baseQuad; fastaString.chunks(4).enumerate)
        {
            static foreach (j; 0 .. 4)
            {{
                if (!baseQuad.empty)
                {
                    data[i][j] = convert(baseQuad.front);
                    baseQuad.popFront();
                }

            }}
        }

        auto compressedSequence = CompressedSequence(data);
        compressedSequence.numBases = fastaString.length;

        return compressedSequence;
    }

    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        static foreach (i, base; testSequence)
        {
            mixin("assert(cs[" ~ i.stringof ~ "] == CompressedBase." ~ base.toLower.to!char ~ ");");
        }
    }


    /// Access underlying buffer directly. This is intended for binary storage.
    @property const(CompressedBaseQuad[]) data() const pure nothrow
    {
        return _data;
    }


    /// Return a lazy range of the stored bases.
    ///
    /// Bugs:
    /// $(UL
    ///     $(LI this procedure requires many allocations because the
    ///         `CompressedBaseQuad`s are turned into ranges with `a[]`.)
    /// )
    @property auto bases(T : CompressedBase, Flag!"reverse" reverse = No.reverse)() const pure nothrow
    {
        static if (reverse)
            return _data
                .retro
                .map!"a[]"
                .map!retro
                .joiner
                .dropExactly(4 * _data.length - length - _baseOffset)
                .takeExactly(length);
        else
            return _data
                .map!"a[]"
                .joiner
                .dropExactly(_baseOffset)
                .takeExactly(length);
    }

    unittest
    {
        import std.algorithm : equal;

        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        with (CompressedBase)
        {
            assert(cs.bases!CompressedBase.equal([
                a, t, g, c, c, a, a, c, t, a, c, t, t, t, g, a, a, c, g, c, g, c,
                c, g, c, a, a, g, g, c, a, c, a, g, g, t, g, c, g, c, c, t,
            ]));
            assert(cs.bases!(CompressedBase, Yes.reverse).equal([
                t, c, c, g, c, g, t, g, g, a, c, a, c, g, g, a, a, c, g, c, c, g, c, g, c, a, a,
                g, t, t, t, c, a, t, c, a, a, c, c, g, t, a,
            ]));
        }
    }


    /// ditto
    @property auto bases(C, Flag!"reverse" reverse = No.reverse)() const pure nothrow if (isSomeChar!C)
    {
        return bases!(CompressedBase, reverse).map!(cb => convert!C(cb));
    }

    unittest
    {
        import std.algorithm : equal;

        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        assert(cs.bases!dchar.equal(testSequence.toLower));
        assert(cs.bases!(dchar, Yes.reverse).equal(testSequence.toLower.retro));
    }


    /// Decompress sequence into an ASCII-encoded string.
    S to(S)() const pure nothrow if (isSomeString!S)
    {
        cast(void) is(S == Char[], Char);
        alias C = Unqual!Char;
        auto fastaSequence = minimallyInitializedArray!(C[])(length);

        foreach (fastaIdx; 0 .. length)
        {
            fastaSequence[fastaIdx] = convert!C(this[fastaIdx]);
        }

        return cast(S) fastaSequence;
    }

    /// ditto
    alias toString = to!string;

    ///
    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        assert(cs.to!string == testSequence.toLower);
    }


    /// Internal property that is exposed for binary storage.
    @property ubyte baseOffset() const pure nothrow
    {
        return _baseOffset;
    }


    /// Returns the number of bases in sequence.
    @property size_t length() const pure nothrow
    {
        return numBases;
    }

    /// ditto
    alias opDollar = length;

    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        assert(cs.length == testSequence.length);
        assert(cs[1 .. $ - 2].length == testSequence.length - 3);
    }


    /// Access base at index `fastaIdx`.
    CompressedBase opIndex(size_t fastaIdx) const pure nothrow
    {
        assert(fastaIdx < length, "index out o bounds" ~ fastaIdx.to!string);

        fastaIdx += _baseOffset;
        auto i = fastaIdx / 4;
        auto j = fastaIdx % 4;

        return data[i][j];
    }

    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        assert(cs[0] == CompressedBase.a);
        assert(cs[1] == CompressedBase.t);
        assert(cs[2] == CompressedBase.g);
        assert(cs[3] == CompressedBase.c);
        assert(cs[$ - 4] == CompressedBase.g);
        assert(cs[$ - 3] == CompressedBase.c);
        assert(cs[$ - 2] == CompressedBase.c);
        assert(cs[$ - 1] == CompressedBase.t);

        assert(cs[1 .. $ - 2][0] == CompressedBase.t);
        assert(cs[1 .. $ - 2][1] == CompressedBase.g);
        assert(cs[1 .. $ - 2][2] == CompressedBase.c);
        assert(cs[1 .. $ - 2][3] == CompressedBase.c);
        assert(cs[1 .. $ - 2][$ - 4] == CompressedBase.g);
        assert(cs[1 .. $ - 2][$ - 3] == CompressedBase.c);
        assert(cs[1 .. $ - 2][$ - 2] == CompressedBase.g);
        assert(cs[1 .. $ - 2][$ - 1] == CompressedBase.c);
    }


    /// Return a slice of this `CompressedSequence` in single base coordinates.
    inout(CompressedSequence) opIndex(size_t[2] slice) inout pure nothrow
    {
        assert(0 <= slice[0] && slice[0] <= slice[1] && slice[1] <= length, "index out of bounds");


        if (slice[1] == slice[0])
            return typeof(return)();

        size_t dataStart = slice[0] / 4;
        size_t dataEnd = ceildiv(slice[1], 4);

        return typeof(return)(
            this._data[dataStart .. dataEnd],
            cast(ubyte) slice[0] % 4,
            slice[1] - slice[0],
        );
    }

    ///
    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        assert(cs[1 .. $ - 2].to!string == testSequence[1 .. $ - 2].toLower);
        assert(cs[0 .. $].to!string == testSequence.toLower);
        assert(cs[0 .. $] == cs);
    }


    size_t[2] opSlice(size_t dim)(in size_t start, in size_t end) const pure nothrow if (dim == 0)
    {
        return [start, end];
    }


    /// Number of bytes required to store the sequence. This does not include
    /// meta data.
    @property size_t compressedLength() const pure nothrow
    {
        return data.length;
    }

    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        assert(cs.compressedLength == 11);
    }


    /// ditto
    static size_t compressedLength(S)(in S fastaSequence) if (isSomeString!S)
    {
        assert(canConvert(fastaSequence));

        return ceildiv(fastaSequence.length, 4);
    }

    ///
    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";

        assert(compressedLength(testSequence) == 11);
    }


    /// Returns true if `fastaSequence` contains only characters `acgt`
    /// ignoring casing.
    static size_t canConvert(S)(in S fastaSequence) if (isSomeString!S)
    {
        return fastaSequence.all!(c => c.toLower.among!('a', 'c', 'g', 't'));
    }

    ///
    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        enum falseSequence = "atgccaactactttgaNNNNNnnnnnAGGCACAGGTGCGCCT";

        assert(canConvert(testSequence));
        assert(!canConvert(falseSequence));
    }


    /// Convert a single ASCII-encoded `fastaBase` to `CompressedBase` and
    /// vice versa.
    ///
    /// Throws: `Exception` if not `canConvert(fastaBase)`
    /// See_also: `canConvert`
    static CompressedBase convert(C)(in C fastaBase) if (isSomeChar!C)
    {
        switch (fastaBase.toLower)
        {
            case 'a':
                return CompressedBase.a;
            case 'c':
                return CompressedBase.c;
            case 't':
                return CompressedBase.t;
            case 'g':
                return CompressedBase.g;
            default:
                throw new Exception("invalid base: " ~ fastaBase.to!string);
        }
    }

    ///
    unittest
    {
        enum bases = "actgACTG";

        static foreach (i, base; bases)
        {
            mixin("assert(convert(base) == CompressedBase." ~ base.toLower.to!char ~ ");");
        }
    }

    /// ditto
    static C convert(C)(in CompressedBase compressedBase)
    {
        final switch (compressedBase)
        {
            case CompressedBase.a:
                return 'a';
            case CompressedBase.c:
                return 'c';
            case CompressedBase.t:
                return 't';
            case CompressedBase.g:
                return 'g';
        }
    }

    ///
    unittest
    {
        enum bases = "actgACTG";

        static foreach (i, base; bases)
        {
            mixin("assert(convert!dchar(CompressedBase." ~ base.toLower.to!char ~ ") == base.toLower);");
        }
    }
};


/**
    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.binio._base;

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


void lockIfPossible(ref File file) nothrow
{
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

class BinaryIOException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

mixin template DbIndex(T...)
{
    static struct EOF;

    @property auto ref size_t beginPtr(T)() pure nothrow
    {
        return fieldPtr!T;
    }

    @property size_t beginPtr(T)() const pure nothrow
    {
        return fieldPtr!T;
    }

    @property auto ref size_t endPtr(T)() pure nothrow
    {
        return fieldPtr!(NextType!T);
    }

    @property size_t endPtr(T)() const pure nothrow
    {
        return fieldPtr!(NextType!T);
    }

    @property ArrayStorage!(StorageType!T) arrayStorage(T)() const pure nothrow
    {
        return typeof(return).fromPtrs(beginPtr!T, endPtr!T);
    }
}

T readRecordAt(T)(File dbFile, size_t ptr)
{
    dbFile.seek(ptr.to!long);

    return readRecord!T(dbFile);
}

T readRecord(T)(File dbFile)
{
    return readRecords(dbFile, new T[1])[0];
}

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

struct ArrayStorage(T)
{
    enum elementSize = T.sizeof;
    size_t ptr;
    size_t length;

    ArrayStorage!T opIndex() const pure nothrow
    {
        return this;
    }

    size_t opIndex(size_t i) const pure nothrow
    {
        assert(i <= length, "out of bounds");

        return ptr + i * elementSize;
    }

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

    size_t opDollar() const pure nothrow
    {
        return length;
    }

    long indexOf(size_t ptr) const pure
    {
        assert((ptr.to!long - this.ptr.to!long) % elementSize.to!long == 0, "bad pointer alignment");

        return (ptr.to!long - this.ptr.to!long) / elementSize.to!long;
    }

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

enum CompressedBase : ubyte
{
    a = 0b00,
    c = 0b01,
    t = 0b10,
    g = 0b11,
}

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

    enum length = 4;
    alias opDollar = length;

    CompressedBase[] opIndex() const pure nothrow
    {
        return [_0, _1, _2, _3];
    }

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

struct CompressedSequence
{
    private CompressedBaseQuad[] _data;
    private ubyte _baseOffset;
    private size_t numBases;

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

    @property const(CompressedBaseQuad[]) data() const pure nothrow
    {
        return _data;
    }

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

    alias toString = to!string;

    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        assert(cs.to!string == testSequence.toLower);
    }

    @property ubyte baseOffset() const pure nothrow
    {
        return _baseOffset;
    }

    @property size_t length() const pure nothrow
    {
        return numBases;
    }

    alias opDollar = length;

    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        assert(cs.length == testSequence.length);
        assert(cs[1 .. $ - 2].length == testSequence.length - 3);
    }

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

    size_t[2] opSlice(size_t dim)(in size_t start, in size_t end) const pure nothrow if (dim == 0)
    {
        return [start, end];
    }

    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        auto cs = CompressedSequence.from(testSequence);

        assert(cs[1 .. $ - 2].to!string == testSequence[1 .. $ - 2].toLower);
        assert(cs[0 .. $].to!string == testSequence.toLower);
        assert(cs[0 .. $] == cs);
    }

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

    static size_t compressedLength(S)(in S fastaSequence) if (isSomeString!S)
    {
        assert(canConvert(fastaSequence));

        return ceildiv(fastaSequence.length, 4);
    }

    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";

        assert(compressedLength(testSequence) == 11);
    }

    static size_t canConvert(S)(in S fastaSequence) if (isSomeString!S)
    {
        return fastaSequence.all!(c => c.toLower.among!('a', 'c', 'g', 't'));
    }

    unittest
    {
        enum testSequence = "atgccaactactttgaacgcgCCGCAAGGCACAGGTGCGCCT";
        enum falseSequence = "atgccaactactttgaNNNNNnnnnnAGGCACAGGTGCGCCT";

        assert(canConvert(testSequence));
        assert(!canConvert(falseSequence));
    }

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

    unittest
    {
        enum bases = "actgACTG";

        static foreach (i, base; bases)
        {
            mixin("assert(convert(base) == CompressedBase." ~ base.toLower.to!char ~ ");");
        }
    }

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

    unittest
    {
        enum bases = "actgACTG";

        static foreach (i, base; bases)
        {
            mixin("assert(convert!dchar(CompressedBase." ~ base.toLower.to!char ~ ") == base.toLower);");
        }
    }
};


/**
    This package contains methods to handle the proprietary binary data
    container for `Insertion`s.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.binio.insertiondb;

import core.exception : AssertError;
import dentist.common : ReferencePoint;
import dentist.common.alignments : AlignmentChain;
import dentist.common.binio._base :
    ArrayStorage,
    CompressedBaseQuad,
    CompressedSequence,
    DbIndex,
    lockIfPossible,
    readRecord,
    readRecordAt,
    readRecords;
import dentist.common.insertions :
    Insertion,
    InsertionInfo,
    SpliceSite;
import dentist.common.scaffold :
    ContigNode,
    ContigPart;
import std.array : minimallyInitializedArray;
import std.conv : to;
import std.exception : assertThrown, enforce, ErrnoException;
import std.format : format;
import std.range :
    ElementType,
    empty,
    front,
    hasLength,
    isForwardRange,
    isInputRange,
    popFront,
    save;
import std.stdio : File;
import std.traits : isArray;
import std.typecons : tuple, Tuple;

version (unittest) import dentist.common.binio._testdata :
    getInsertionsTestData,
    numCompressedBaseQuads,
    numInsertions,
    numSpliceSites;


class InsertionDbException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

struct InsertionDb
{
    private alias DbSlices = Tuple!(
        ArrayStorage!(StorageType!Insertion), "insertions",
        ArrayStorage!(StorageType!CompressedBaseQuad), "compressedBaseQuads",
        ArrayStorage!(StorageType!SpliceSite), "spliceSites",
    );

    private File file;
    private InsertionDbIndex index;
    private DbSlices slices;

    @property auto insertions() const pure nothrow
    {
        return index.insertions;
    }

    @property auto compressedBaseQuads() const pure nothrow
    {
        return index.compressedBaseQuads;
    }

    @property auto spliceSites() const pure nothrow
    {
        return index.spliceSites;
    }

    static InsertionDb parse(in string dbFile)
    {
        auto file = File(dbFile, "rb");
        lockIfPossible(file);
        auto db = InsertionDb(file);
        db.ensureDbIndex();

        return db;
    }

    void releaseDb()
    {
        file.close();
    }

    Insertion[] opIndex()
    {
        ensureDbIndex();

        return readSlice(0, length);
    }

    Insertion opIndex(size_t i)
    {
        ensureDbIndex();
        enforce!InsertionDbException(
            i < length,
            format!"cannot read block %d in `%s`: out of bounds [0, %d)"(
                    i, file.name, length)
        );

        return readSlice(i, i + 1)[0];
    }

    Insertion[] opIndex(size_t[2] slice)
    {
        auto from = slice[0];
        auto to = slice[1];
        ensureDbIndex();
        enforce!InsertionDbException(
            to <= length,
            format!"cannot read blocks %d-%d in `%s`: out of bounds [0, %d]"(
                    from, to, file.name, length)
        );

        return readSlice(from, to);
    }

    size_t[2] opSlice(size_t dim)(size_t from, size_t to)
            if (dim == 0)
    {
        assert(from < to, "invalid slice");

        return [from, to];
    }

    @property size_t length()
    {
        ensureDbIndex();

        return index.insertions.length;
    }

    alias opDollar = length;

    private void ensureDbIndex()
    {
        if (index != index.init)
            return;

        index = file.readRecord!InsertionDbIndex();
    }

    private Insertion[] readSlice(size_t from, size_t to)
    {
        assert(from <= to && to <= length);

        // Step 1: determine memory requirements and DB slices
        slices = getSlices(from, to);

        // Step 2: allocate minimally initialized memory for all blocks
        auto insertions = minimallyInitializedArray!(Insertion[])(slices.insertions.length);
        auto compressedBaseQuads = minimallyInitializedArray!(CompressedBaseQuad[])(slices.compressedBaseQuads.length);
        auto spliceSites = minimallyInitializedArray!(SpliceSite[])(slices.spliceSites.length);

        // Step 3: parse each record for each block assigning already
        //         allocated array slices to the array fields
        parse(insertions, compressedBaseQuads, spliceSites);
        parse(compressedBaseQuads);
        parse(spliceSites);

        return insertions;
    }

    private DbSlices getSlices(size_t from, size_t to)
    {
        auto insertions = index.insertions[from .. to];
        auto firstInsertion = file.readRecordAt!(StorageType!Insertion)(insertions[0]);
        auto lastInsertion = file.readRecordAt!(StorageType!Insertion)(insertions[$ - 1]);

        auto compressedBaseQuads = ArrayStorage!(StorageType!CompressedBaseQuad).fromPtrs(
            firstInsertion.sequence[0],
            lastInsertion.sequence[$],
        );

        auto spliceSites = ArrayStorage!(StorageType!SpliceSite).fromPtrs(
            firstInsertion.spliceSites[0],
            lastInsertion.spliceSites[$],
        );

        return DbSlices(
            insertions,
            compressedBaseQuads,
            spliceSites,
        );
    }

    private void parse(
        ref Insertion[] insertions,
        CompressedBaseQuad[] compressedBaseQuads,
        SpliceSite[] spliceSites,
    )
    {
        file.seek(slices.insertions.ptr);

        size_t[2] compressedBaseQuadsSlice;
        size_t[2] spliceSitesSlice;

        foreach (ref insertion; insertions)
        {
            auto insertionStorage = file.readRecord!(StorageType!Insertion);

            compressedBaseQuadsSlice[0] = compressedBaseQuadsSlice[1];
            compressedBaseQuadsSlice[1] += insertionStorage.sequence.length;

            spliceSitesSlice[0] = spliceSitesSlice[1];
            spliceSitesSlice[1] += insertionStorage.spliceSites.length;

            insertion = Insertion(
                insertionStorage.start,
                insertionStorage.end,
                InsertionInfo(
                    CompressedSequence(
                        compressedBaseQuads[
                            compressedBaseQuadsSlice[0] .. compressedBaseQuadsSlice[1]
                        ],
                        insertionStorage.baseOffset,
                        insertionStorage.sequenceLength,
                    ),
                    insertionStorage.contigLength,
                    spliceSites[spliceSitesSlice[0] .. spliceSitesSlice[1]],
                ),
            );
        }
    }

    private void parse(
        ref SpliceSite[] spliceSites,
    )
    {
        static assert(SpliceSite.sizeof == StorageType!SpliceSite.sizeof);
        file.seek(slices.spliceSites.ptr);
        spliceSites = file.readRecords(spliceSites);
    }

    private void parse(
        ref CompressedBaseQuad[] compressedBaseQuads,
    )
    {
        static assert(CompressedBaseQuad.sizeof == StorageType!CompressedBaseQuad.sizeof);
        file.seek(slices.compressedBaseQuads.ptr);
        compressedBaseQuads = file.readRecords(compressedBaseQuads);
    }

    static void write(R)(in string dbFile, R insertions)
            if (isForwardRange!R && hasLength!R && is(ElementType!R : const(Insertion)))
    {
        auto writer = InsertionDbFileWriter!R(File(dbFile, "wb"), insertions);

        lockIfPossible(writer.file);
        writer.writeToFile();
    }
}

unittest
{
    import dentist.util.tempfile : mkstemp;
    import std.file : remove;

    auto insertions = getInsertionsTestData();

    enum totalDbSize =
        InsertionDbIndex.sizeof +
        StorageType!Insertion.sizeof * numInsertions +
        StorageType!CompressedBaseQuad.sizeof * numCompressedBaseQuads +
        StorageType!SpliceSite.sizeof * numSpliceSites;

    auto tmpDb = mkstemp("./.unittest-XXXXXX");
    scope (exit)
    {
        tmpDb.file.close();
        remove(tmpDb.name);
    }

    InsertionDbFileWriter!(Insertion[])(tmpDb.file, insertions).writeToFile();
    tmpDb.file.sync();

    assert(tmpDb.file.size == totalDbSize);

    tmpDb.file.rewind();
    auto insertionDb = InsertionDb(tmpDb.file);

    assert(insertionDb[] == insertions);
}

private struct InsertionDbFileWriter(R)
        if (isForwardRange!R && hasLength!R && is(ElementType!R : const(Insertion)))
{
    File file;
    R insertions;
    InsertionDbIndex index;

    void writeToFile()
    {
        index = InsertionDbIndex.from(insertions.save);

        file.rawWrite([index]);
        writeBlock!Insertion();
        writeBlock!CompressedBaseQuad();
        writeBlock!SpliceSite();
    }

    void writeBlock(T : Insertion)()
    {
        auto compressedBaseQuads = index.compressedBaseQuads;
        compressedBaseQuads.length = 0;
        auto spliceSites = index.spliceSites;
        spliceSites.length = 0;

        version (assert)
        {
            auto insertions = index.insertions;
            insertions.length = 0;
            assert(insertions.ptr == file.tell());
        }
        foreach (insertion; this.insertions.save)
        {
            compressedBaseQuads.length = insertion.payload.sequence.compressedLength;
            spliceSites.length = insertion.payload.spliceSites.length;
            auto insertionStorage = InsertionStorage(
                insertion.start,
                insertion.end,
                insertion.payload.sequence.baseOffset,
                insertion.payload.sequence.length,
                compressedBaseQuads,
                insertion.payload.contigLength,
                spliceSites,
            );

            file.rawWrite([insertionStorage]);

            compressedBaseQuads.ptr = compressedBaseQuads[$];
            spliceSites.ptr = spliceSites[$];
            version (assert)
            {
                ++insertions.length;
                assert(insertions[$] == file.tell());
            }
        }
    }

    void writeBlock(T : CompressedBaseQuad)()
    {
        version (assert)
        {
            auto compressedBaseQuads = index.compressedBaseQuads;
            compressedBaseQuads.length = 0;
            assert(compressedBaseQuads.ptr == file.tell());
        }
        foreach (insertion; this.insertions.save)
        {
            static assert(CompressedBaseQuad.sizeof == StorageType!CompressedBaseQuad.sizeof);
            file.rawWrite(insertion.payload.sequence.data);

            version (assert)
            {
                compressedBaseQuads.length += insertion.payload.sequence.compressedLength;
                assert(compressedBaseQuads[$] == file.tell());
            }
        }
    }

    void writeBlock(T : SpliceSite)()
    {
        version (assert)
        {
            auto spliceSites = index.spliceSites;
            spliceSites.length = 0;
            assert(spliceSites.ptr == file.tell());
        }
        foreach (insertion; this.insertions.save)
        {
            static assert(SpliceSite.sizeof == StorageType!SpliceSite.sizeof);
            file.rawWrite(insertion.payload.spliceSites);

            version (assert)
            {
                spliceSites.length += insertion.payload.spliceSites.length;
                assert(spliceSites[$] == file.tell());
            }
        }
    }
}

private struct InsertionDbIndex
{
    mixin DbIndex;

    private static template NextType(T)
    {
        static if (is(T == Insertion))
            alias NextType = CompressedBaseQuad;
        else static if (is(T == CompressedBaseQuad))
            alias NextType = SpliceSite;
        else static if (is(T == SpliceSite))
            alias NextType = EOF;
    }

    private static template fieldPtr(T)
    {
        static if (is(T == Insertion))
            alias fieldPtr = insertionsPtr;
        else static if (is(T == CompressedBaseQuad))
            alias fieldPtr = compressedBaseQuadsPtr;
        else static if (is(T == SpliceSite))
            alias fieldPtr = spliceSitesPtr;
        else static if (is(T == EOF))
            alias fieldPtr = eofPtr;
    }

    size_t insertionsPtr;
    size_t compressedBaseQuadsPtr;
    size_t spliceSitesPtr;
    size_t eofPtr;

    @property alias insertions = arrayStorage!Insertion;
    @property alias compressedBaseQuads = arrayStorage!CompressedBaseQuad;
    @property alias spliceSites = arrayStorage!SpliceSite;

    static InsertionDbIndex from(R)(R insertions) nothrow pure
            if (isInputRange!R && hasLength!R && is(ElementType!R : const(Insertion)))
    {
        InsertionDbIndex index;

        index.beginPtr!Insertion = InsertionDbIndex.sizeof;
        index.endPtr!Insertion = StorageType!Insertion.sizeof * insertions.length;
        foreach (insertion; insertions)
        {
            index.endPtr!CompressedBaseQuad += StorageType!CompressedBaseQuad.sizeof *
                    insertion.payload.sequence.compressedLength;
            index.endPtr!SpliceSite += StorageType!SpliceSite.sizeof *
                    insertion.payload.spliceSites.length;
        }

        index.compressedBaseQuadsPtr += index.insertionsPtr;
        index.spliceSitesPtr += index.compressedBaseQuadsPtr;
        index.eofPtr += index.spliceSitesPtr;

        return index;
    }

    unittest
    {
        auto dbIndex = InsertionDbIndex.from(getInsertionsTestData());

        assert(dbIndex.insertionsPtr == InsertionDbIndex.sizeof);
        assert(dbIndex.compressedBaseQuadsPtr ==
                dbIndex.insertionsPtr +
                StorageType!Insertion.sizeof * numInsertions);
        assert(dbIndex.spliceSitesPtr ==
                dbIndex.compressedBaseQuadsPtr +
                StorageType!CompressedBaseQuad.sizeof * numCompressedBaseQuads);
        assert(dbIndex.eofPtr ==
                dbIndex.spliceSitesPtr +
                StorageType!SpliceSite.sizeof * numSpliceSites);
    }
}

unittest
{
    enum begin = 1;
    enum end = 2;
    enum modified = 3;

    {
        InsertionDbIndex dbIndex;

        dbIndex.insertionsPtr = begin;
        dbIndex.compressedBaseQuadsPtr = end;

        assert(dbIndex.beginPtr!Insertion == begin);
        assert(dbIndex.endPtr!Insertion == end);

        dbIndex.beginPtr!Insertion = modified;
        dbIndex.endPtr!Insertion = modified;

        assert(dbIndex.insertionsPtr == modified);
        assert(dbIndex.compressedBaseQuadsPtr == modified);
    }
    {
        InsertionDbIndex dbIndex;

        dbIndex.compressedBaseQuadsPtr = begin;
        dbIndex.spliceSitesPtr = end;

        assert(dbIndex.beginPtr!CompressedBaseQuad == begin);
        assert(dbIndex.endPtr!CompressedBaseQuad == end);

        dbIndex.beginPtr!CompressedBaseQuad = modified;
        dbIndex.endPtr!CompressedBaseQuad = modified;

        assert(dbIndex.compressedBaseQuadsPtr == modified);
        assert(dbIndex.spliceSitesPtr == modified);
    }
    {
        InsertionDbIndex dbIndex;

        dbIndex.spliceSitesPtr = begin;
        dbIndex.eofPtr = end;

        assert(dbIndex.beginPtr!SpliceSite == begin);
        assert(dbIndex.endPtr!SpliceSite == end);

        dbIndex.beginPtr!SpliceSite = modified;
        dbIndex.endPtr!SpliceSite = modified;

        assert(dbIndex.spliceSitesPtr == modified);
        assert(dbIndex.eofPtr == modified);
    }
}

private template StorageType(T)
{
    static if (is(T == Insertion))
        alias StorageType = InsertionStorage;
    else static if (is(T == CompressedBaseQuad[]))
        alias StorageType = ArrayStorage!(StorageType!CompressedBaseQuad);
    else static if (is(T == CompressedBaseQuad))
        alias StorageType = CompressedBaseQuad;
    else static if (is(T == SpliceSite))
        alias StorageType = SpliceSite;
    else static if (is(T == SpliceSite[]))
        alias StorageType = ArrayStorage!(StorageType!SpliceSite);
}

private struct InsertionStorage
{
    ContigNode start;
    ContigNode end;
    ubyte baseOffset;
    size_t sequenceLength;
    StorageType!(CompressedBaseQuad[]) sequence;
    size_t contigLength;
    StorageType!(SpliceSite[]) spliceSites;
}

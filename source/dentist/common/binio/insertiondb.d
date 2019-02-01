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
import dentist.common.alignments :
    AlignmentChain,
    AlignmentLocationSeed,
    coord_t,
    diff_t,
    id_t,
    SeededAlignment,
    trace_point_t;
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
    InsertionInfo;
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

version (unittest) import dentist.common.binio._testdata.insertiondb :
    getInsertionsTestData,
    numInsertions,
    numCompressedBaseQuads,
    numOverlaps,
    numLocalAlignments,
    numTracePoints;


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
    private alias LocalAlignment = AlignmentChain.LocalAlignment;
    private alias TracePoint = LocalAlignment.TracePoint;
    private alias DbSlices = Tuple!(
        ArrayStorage!(StorageType!Insertion), "insertions",
        ArrayStorage!(StorageType!CompressedBaseQuad), "compressedBaseQuads",
        ArrayStorage!(StorageType!SeededAlignment), "overlaps",
        ArrayStorage!(StorageType!LocalAlignment), "localAlignments",
        ArrayStorage!(StorageType!TracePoint), "tracePoints",
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

    @property auto overlaps() const pure nothrow
    {
        return index.overlaps;
    }

    @property auto localAlignments() const pure nothrow
    {
        return index.localAlignments;
    }

    @property auto tracePoints() const pure nothrow
    {
        return index.tracePoints;
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

        if (from == to)
            return [];

        // Step 1: determine memory requirements and DB slices
        slices = getSlices(from, to);

        // Step 2: allocate minimally initialized memory for all blocks
        auto insertions = minimallyInitializedArray!(Insertion[])(slices.insertions.length);
        auto compressedBaseQuads = minimallyInitializedArray!(CompressedBaseQuad[])(slices.compressedBaseQuads.length);
        auto overlaps = minimallyInitializedArray!(SeededAlignment[])(slices.overlaps.length);
        auto localAlignments = minimallyInitializedArray!(LocalAlignment[])(slices.localAlignments.length);
        auto tracePoints = minimallyInitializedArray!(TracePoint[])(slices.tracePoints.length);

        // Step 3: parse each record for each block assigning already
        //         allocated array slices to the array fields
        parse(insertions, compressedBaseQuads, overlaps);
        parse(compressedBaseQuads);
        parse(overlaps, localAlignments);
        parse(localAlignments, tracePoints);
        parse(tracePoints);

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

        auto overlaps = ArrayStorage!(StorageType!SeededAlignment).fromPtrs(
            firstInsertion.overlaps[0],
            lastInsertion.overlaps[$],
        );
        auto firstSeededAlignment = file.readRecordAt!(StorageType!SeededAlignment)(overlaps[0]);
        auto lastSeededAlignment = file.readRecordAt!(StorageType!SeededAlignment)(overlaps[$ - 1]);

        auto localAlignments = ArrayStorage!(StorageType!LocalAlignment).fromPtrs(
            firstSeededAlignment.localAlignments[0],
            lastSeededAlignment.localAlignments[$],
        );
        auto firstLocalAlignment = file.readRecordAt!(StorageType!LocalAlignment)(localAlignments[0]);
        auto lastLocalAlignment = file.readRecordAt!(StorageType!LocalAlignment)(localAlignments[$ - 1]);

        auto tracePoints = ArrayStorage!(StorageType!TracePoint).fromPtrs(
            firstLocalAlignment.tracePoints[0],
            lastLocalAlignment.tracePoints[$],
        );

        return DbSlices(
            insertions,
            compressedBaseQuads,
            overlaps,
            localAlignments,
            tracePoints,
        );
    }

    private void parse(
        ref Insertion[] insertions,
        CompressedBaseQuad[] compressedBaseQuads,
        SeededAlignment[] overlaps,
    )
    {
        file.seek(slices.insertions.ptr);

        size_t[2] compressedBaseQuadsSlice;
        size_t[2] overlapsSlice;

        foreach (ref insertion; insertions)
        {
            auto insertionStorage = file.readRecord!(StorageType!Insertion);

            compressedBaseQuadsSlice[0] = compressedBaseQuadsSlice[1];
            compressedBaseQuadsSlice[1] += insertionStorage.sequence.length;

            overlapsSlice[0] = overlapsSlice[1];
            overlapsSlice[1] += insertionStorage.overlaps.length;

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
                    overlaps[overlapsSlice[0] .. overlapsSlice[1]],
                ),
            );
        }
    }

    private void parse(
        ref CompressedBaseQuad[] compressedBaseQuads,
    )
    {
        static assert(CompressedBaseQuad.sizeof == StorageType!CompressedBaseQuad.sizeof);
        file.seek(slices.compressedBaseQuads.ptr);
        compressedBaseQuads = file.readRecords(compressedBaseQuads);
    }

    private void parse(
        ref SeededAlignment[] overlaps,
        AlignmentChain.LocalAlignment[] localAlignments
    )
    {
        // Parse `SeededAlignment`s
        alias Contig = AlignmentChain.Contig;
        file.seek(slices.overlaps.ptr);

        size_t[2] localAlignmentsSlice;
        foreach (ref overlap; overlaps)
        {
            auto overlapStorage = file.readRecord!(StorageType!SeededAlignment);

            localAlignmentsSlice[0] = localAlignmentsSlice[1];
            localAlignmentsSlice[1] += overlapStorage.localAlignments.length;

            overlap = SeededAlignment(
                AlignmentChain(
                    overlapStorage.id,
                    Contig(
                        overlapStorage.contigAId,
                        overlapStorage.contigALength,
                    ),
                    Contig(
                        overlapStorage.contigBId,
                        overlapStorage.contigBLength,
                    ),
                    overlapStorage.flags,
                    localAlignments[localAlignmentsSlice[0] .. localAlignmentsSlice[1]],
                    overlapStorage.tracePointDistance,
                ),
                overlapStorage.seed,
            );
        }
    }

    private void parse(
        ref AlignmentChain.LocalAlignment[] localAlignments,
        AlignmentChain.LocalAlignment.TracePoint[] tracePoints
    )
    {
        alias LocalAlignment = AlignmentChain.LocalAlignment;
        alias Locus = LocalAlignment.Locus;
        file.seek(slices.localAlignments.ptr);

        size_t[2] localAlignmentsSlice;
        foreach (ref localAlignment; localAlignments)
        {
            auto localAlignmentStorage = file.readRecord!(StorageType!LocalAlignment);

            localAlignmentsSlice[0] = localAlignmentsSlice[1];
            localAlignmentsSlice[1] += localAlignmentStorage.tracePoints.length;

            localAlignment = LocalAlignment(
                Locus(
                    localAlignmentStorage.contigABegin,
                    localAlignmentStorage.contigAEnd,
                ),
                Locus(
                    localAlignmentStorage.contigBBegin,
                    localAlignmentStorage.contigBEnd,
                ),
                localAlignmentStorage.numDiffs,
                tracePoints[localAlignmentsSlice[0] .. localAlignmentsSlice[1]],
            );
        }
    }

    private void parse(ref AlignmentChain.LocalAlignment.TracePoint[] tracePoints)
    {
        alias LocalAlignment = AlignmentChain.LocalAlignment;
        alias TracePoint = LocalAlignment.TracePoint;

        static assert(TracePoint.sizeof == StorageType!TracePoint.sizeof);
        file.seek(slices.tracePoints.ptr);
        tracePoints = file.readRecords(tracePoints);
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

    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

    auto insertions = getInsertionsTestData();

    enum totalDbSize =
        InsertionDbIndex.sizeof +
        StorageType!Insertion.sizeof * numInsertions +
        StorageType!CompressedBaseQuad.sizeof * numCompressedBaseQuads +
        StorageType!SeededAlignment.sizeof * numOverlaps +
        StorageType!LocalAlignment.sizeof * numLocalAlignments +
        StorageType!TracePoint.sizeof * numTracePoints;

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
    private alias LocalAlignment = AlignmentChain.LocalAlignment;
    private alias TracePoint = LocalAlignment.TracePoint;

    File file;
    R insertions;
    InsertionDbIndex index;

    void writeToFile()
    {
        index = InsertionDbIndex.from(insertions.save);

        file.rawWrite([index]);
        writeBlock!Insertion();
        writeBlock!CompressedBaseQuad();
        writeBlock!SeededAlignment();
        writeBlock!LocalAlignment();
        writeBlock!TracePoint();
    }

    void writeBlock(T : Insertion)()
    {
        auto compressedBaseQuads = index.compressedBaseQuads;
        compressedBaseQuads.length = 0;
        auto overlaps = index.overlaps;
        overlaps.length = 0;

        version (assert)
        {
            auto insertions = index.insertions;
            insertions.length = 0;
            assert(insertions.ptr == file.tell());
        }
        foreach (insertion; this.insertions.save)
        {
            compressedBaseQuads.length = insertion.payload.sequence.compressedLength;
            overlaps.length = insertion.payload.overlaps.length;
            auto insertionStorage = InsertionStorage(
                insertion.start,
                insertion.end,
                insertion.payload.sequence.baseOffset,
                insertion.payload.sequence.length,
                compressedBaseQuads,
                insertion.payload.contigLength,
                overlaps,
            );

            file.rawWrite([insertionStorage]);

            compressedBaseQuads.ptr = compressedBaseQuads[$];
            overlaps.ptr = overlaps[$];
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

    void writeBlock(T : SeededAlignment)()
    {
        auto localAlignments = index.localAlignments;
        localAlignments.length = 0;

        version (assert)
        {
            auto overlaps = index.overlaps;
            overlaps.length = 0;
            assert(overlaps.ptr == file.tell());
        }

        foreach (insertion; this.insertions.save)
        {
            foreach (overlap; insertion.payload.overlaps)
            {
                localAlignments.length = overlap.localAlignments.length;

                auto overlapStorage = SeededAlignmentStorage(
                    overlap.id,
                    overlap.contigA.id,
                    overlap.contigA.length,
                    overlap.contigB.id,
                    overlap.contigB.length,
                    overlap.flags,
                    localAlignments,
                    overlap.tracePointDistance,
                    overlap.seed,
                );

                file.rawWrite([overlapStorage]);

                localAlignments.ptr = localAlignments[$];
                version (assert)
                {
                    ++overlaps.length;
                    assert(overlaps[$] == file.tell());
                }
            }
        }
    }

    void writeBlock(T : LocalAlignment)()
    {
        auto tracePoints = index.tracePoints;
        tracePoints.length = 0;

        version (assert)
        {
            auto localAlignments = index.localAlignments;
            localAlignments.length = 0;
            assert(localAlignments.ptr == file.tell());
        }

        foreach (insertion; this.insertions.save)
        {
            foreach (overlap; insertion.payload.overlaps)
            {
                foreach (localAlignment; overlap.localAlignments)
                {
                    tracePoints.length = localAlignment.tracePoints.length;

                    auto localAlignmentStorage = LocalAlignmentStorage(
                        localAlignment.contigA.begin,
                        localAlignment.contigA.end,
                        localAlignment.contigB.begin,
                        localAlignment.contigB.end,
                        localAlignment.numDiffs,
                        tracePoints,
                    );

                    file.rawWrite([localAlignmentStorage]);

                    tracePoints.ptr = tracePoints[$];
                    version (assert)
                    {
                        ++localAlignments.length;
                        assert(localAlignments[$] == file.tell());
                    }
                }
            }
        }
    }

    void writeBlock(T : TracePoint)()
    {
        version (assert)
        {
            auto tracePoints = index.tracePoints;
            tracePoints.length = 0;
            assert(tracePoints.ptr == file.tell());
        }

        foreach (insertion; this.insertions.save)
        {
            foreach (overlap; insertion.payload.overlaps)
            {
                foreach (localAlignment; overlap.localAlignments)
                {
                    static assert(TracePoint.sizeof == StorageType!TracePoint.sizeof);
                    file.rawWrite(localAlignment.tracePoints);

                    version (assert)
                    {
                        tracePoints.length += localAlignment.tracePoints.length;
                        assert(tracePoints[$] == file.tell());
                    }
                }
            }
        }
    }
}

private struct InsertionDbIndex
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

    mixin DbIndex;

    private static template NextType(T)
    {
        static if (is(T == Insertion))
            alias NextType = CompressedBaseQuad;
        else static if (is(T == CompressedBaseQuad))
            alias NextType = SeededAlignment;
        else static if (is(T == SeededAlignment))
            alias NextType = LocalAlignment;
        else static if (is(T == LocalAlignment))
            alias NextType = TracePoint;
        else static if (is(T == TracePoint))
            alias NextType = EOF;
    }

    private static template fieldPtr(T)
    {
        static if (is(T == Insertion))
            alias fieldPtr = insertionsPtr;
        else static if (is(T == CompressedBaseQuad))
            alias fieldPtr = compressedBaseQuadsPtr;
        else static if (is(T == SeededAlignment))
            alias fieldPtr = overlapsPtr;
        else static if (is(T == LocalAlignment))
            alias fieldPtr = localAlignmentsPtr;
        else static if (is(T == TracePoint))
            alias fieldPtr = tracePointsPtr;
        else static if (is(T == EOF))
            alias fieldPtr = eofPtr;
    }

    size_t insertionsPtr;
    size_t compressedBaseQuadsPtr;
    size_t overlapsPtr;
    size_t localAlignmentsPtr;
    size_t tracePointsPtr;
    size_t eofPtr;

    @property alias insertions = arrayStorage!Insertion;
    @property alias compressedBaseQuads = arrayStorage!CompressedBaseQuad;
    @property alias overlaps = arrayStorage!SeededAlignment;
    @property alias localAlignments = arrayStorage!LocalAlignment;
    @property alias tracePoints = arrayStorage!TracePoint;

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
            index.endPtr!SeededAlignment += StorageType!SeededAlignment.sizeof *
                    insertion.payload.overlaps.length;

            foreach (overlap; insertion.payload.overlaps)
            {
                index.endPtr!LocalAlignment += StorageType!LocalAlignment.sizeof *
                        overlap.localAlignments.length;

                foreach (localAlignment; overlap.localAlignments)
                    index.endPtr!TracePoint += StorageType!TracePoint.sizeof *
                            localAlignment.tracePoints.length;
            }
        }

        index.compressedBaseQuadsPtr += index.insertionsPtr;
        index.overlapsPtr += index.compressedBaseQuadsPtr;
        index.localAlignmentsPtr += index.overlapsPtr;
        index.tracePointsPtr += index.localAlignmentsPtr;
        index.eofPtr += index.tracePointsPtr;

        return index;
    }

    unittest
    {
        auto dbIndex = InsertionDbIndex.from(getInsertionsTestData());

        assert(dbIndex.insertionsPtr == InsertionDbIndex.sizeof);
        assert(dbIndex.compressedBaseQuadsPtr ==
                dbIndex.insertionsPtr +
                StorageType!Insertion.sizeof * numInsertions);
        assert(dbIndex.overlapsPtr ==
                dbIndex.compressedBaseQuadsPtr +
                StorageType!CompressedBaseQuad.sizeof * numCompressedBaseQuads);
        assert(dbIndex.localAlignmentsPtr ==
                dbIndex.overlapsPtr +
                StorageType!SeededAlignment.sizeof * numOverlaps);
        assert(dbIndex.tracePointsPtr ==
                dbIndex.localAlignmentsPtr +
                StorageType!LocalAlignment.sizeof * numLocalAlignments);
        assert(dbIndex.eofPtr ==
                dbIndex.tracePointsPtr +
                StorageType!TracePoint.sizeof * numTracePoints);
    }
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

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
        dbIndex.overlapsPtr = end;

        assert(dbIndex.beginPtr!CompressedBaseQuad == begin);
        assert(dbIndex.endPtr!CompressedBaseQuad == end);

        dbIndex.beginPtr!CompressedBaseQuad = modified;
        dbIndex.endPtr!CompressedBaseQuad = modified;

        assert(dbIndex.compressedBaseQuadsPtr == modified);
        assert(dbIndex.overlapsPtr == modified);
    }
    {
        InsertionDbIndex dbIndex;

        dbIndex.overlapsPtr = begin;
        dbIndex.localAlignmentsPtr = end;

        assert(dbIndex.beginPtr!SeededAlignment == begin);
        assert(dbIndex.endPtr!SeededAlignment == end);

        dbIndex.beginPtr!SeededAlignment = modified;
        dbIndex.endPtr!SeededAlignment = modified;

        assert(dbIndex.overlapsPtr == modified);
        assert(dbIndex.localAlignmentsPtr == modified);
    }
    {
        InsertionDbIndex dbIndex;

        dbIndex.localAlignmentsPtr = begin;
        dbIndex.tracePointsPtr = end;

        assert(dbIndex.beginPtr!LocalAlignment == begin);
        assert(dbIndex.endPtr!LocalAlignment == end);

        dbIndex.beginPtr!LocalAlignment = modified;
        dbIndex.endPtr!LocalAlignment = modified;

        assert(dbIndex.localAlignmentsPtr == modified);
        assert(dbIndex.tracePointsPtr == modified);
    }
    {
        InsertionDbIndex dbIndex;

        dbIndex.tracePointsPtr = begin;
        dbIndex.eofPtr = end;

        assert(dbIndex.beginPtr!TracePoint == begin);
        assert(dbIndex.endPtr!TracePoint == end);

        dbIndex.beginPtr!TracePoint = modified;
        dbIndex.endPtr!TracePoint = modified;

        assert(dbIndex.tracePointsPtr == modified);
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
    else static if (is(T == SeededAlignment[]))
        alias StorageType = ArrayStorage!(StorageType!SeededAlignment);
    else static if (is(T == SeededAlignment))
        alias StorageType = SeededAlignmentStorage;
    else static if (is(T == AlignmentChain.LocalAlignment[]))
        alias StorageType = ArrayStorage!(StorageType!(AlignmentChain.LocalAlignment));
    else static if (is(T == AlignmentChain.LocalAlignment))
        alias StorageType = LocalAlignmentStorage;
    else static if (is(T == AlignmentChain.LocalAlignment.TracePoint[]))
        alias StorageType = ArrayStorage!(StorageType!(AlignmentChain.LocalAlignment.TracePoint));
    else static if (is(T == AlignmentChain.LocalAlignment.TracePoint))
        alias StorageType = TracePointStorage;
}

private struct InsertionStorage
{
    ContigNode start;
    ContigNode end;
    ubyte baseOffset;
    size_t sequenceLength;
    StorageType!(CompressedBaseQuad[]) sequence;
    size_t contigLength;
    StorageType!(SeededAlignment[]) overlaps;
}

private struct SeededAlignmentStorage
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;

    id_t id;
    id_t contigAId;
    coord_t contigALength;
    id_t contigBId;
    coord_t contigBLength;
    AlignmentChain.Flags flags;
    StorageType!(LocalAlignment[]) localAlignments;
    trace_point_t tracePointDistance;
    AlignmentLocationSeed seed;
}

private struct LocalAlignmentStorage
{
    alias TracePoint = AlignmentChain.LocalAlignment.TracePoint;

    coord_t contigABegin;
    coord_t contigAEnd;
    coord_t contigBBegin;
    coord_t contigBEnd;
    diff_t numDiffs;
    StorageType!(TracePoint[]) tracePoints;
}

private struct TracePointStorage
{
    trace_point_t numDiffs;
    trace_point_t numBasePairs;
}

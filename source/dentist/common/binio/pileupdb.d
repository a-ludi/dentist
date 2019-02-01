/**
    This package contains methods to handle the proprietary binary data
    container for `PileUp`s.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.binio.pileupdb;

import core.exception : AssertError;
import dentist.common.alignments :
    AlignmentChain,
    AlignmentLocationSeed,
    coord_t,
    diff_t,
    id_t,
    PileUp,
    ReadAlignment,
    SeededAlignment,
    trace_point_t;
import dentist.common.binio._base :
    ArrayStorage,
    DbIndex,
    lockIfPossible,
    readRecord,
    readRecordAt,
    readRecords;
import std.array : minimallyInitializedArray;
import std.conv : to;
import std.exception : assertThrown, enforce, ErrnoException;
import std.format : format;
import std.typecons : tuple, Tuple;
import std.stdio : File;

version (unittest) import dentist.common.binio._testdata.pileupdb :
    getPileUpsTestData,
    numLocalAlignments,
    numPileUps,
    numReadAlignments,
    numSeededAlignments,
    numTracePoints;

class PileUpDbException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

struct PileUpDb
{
    private alias LocalAlignment = AlignmentChain.LocalAlignment;
    private alias TracePoint = LocalAlignment.TracePoint;
    private alias DbSlices = Tuple!(
        ArrayStorage!(StorageType!PileUp), "pileUps",
        ArrayStorage!(StorageType!ReadAlignment), "readAlignments",
        ArrayStorage!(StorageType!SeededAlignment), "seededAlignments",
        ArrayStorage!(StorageType!LocalAlignment), "localAlignments",
        ArrayStorage!(StorageType!TracePoint), "tracePoints",
    );

    private File pileUpDb;
    private PileUpDbIndex dbIndex;
    private DbSlices dbSlices;

    @property auto pileUps() const pure nothrow
    {
        return dbIndex.pileUps;
    }
    @property auto readAlignments() const pure nothrow
    {
        return dbIndex.readAlignments;
    }
    @property auto seededAlignments() const pure nothrow
    {
        return dbIndex.seededAlignments;
    }
    @property auto localAlignments() const pure nothrow
    {
        return dbIndex.localAlignments;
    }
    @property auto tracePoints() const pure nothrow
    {
        return dbIndex.tracePoints;
    }

    static PileUpDb parse(in string dbFile)
    {
        auto pileUpDb = File(dbFile, "rb");
        lockIfPossible(pileUpDb);

        auto db = PileUpDb(pileUpDb);
        db.ensureDbIndex();

        return db;
    }

    void releaseDb()
    {
        pileUpDb.close();
    }

    PileUp[] opIndex()
    {
        ensureDbIndex();

        return readSlice(0, length);
    }

    PileUp opIndex(size_t i)
    {
        ensureDbIndex();
        enforce!PileUpDbException(
            i < length,
            format!"cannot read block %d in `%s`: out of bounds [0, %d)"(
                    i, pileUpDb.name, length)
        );

        return readSlice(i, i + 1)[0];
    }

    PileUp[] opIndex(size_t[2] slice)
    {
        auto from = slice[0];
        auto to = slice[1];
        ensureDbIndex();
        enforce!PileUpDbException(
            from <= to && to <= length,
            format!"cannot read blocks %d-%d in `%s`: out of bounds [0, %d]"(
                    from, to, pileUpDb.name, length)
        );

        return readSlice(from, to);
    }

    size_t[2] opSlice(size_t dim)(size_t from, size_t to)
            if (dim == 0)
    {
        assert(from <= to, "invalid slice");

        return [from, to];
    }

    @property size_t length()
    {
        ensureDbIndex();

        return dbIndex.pileUps.length;
    }

    alias opDollar = length;

    private void ensureDbIndex()
    {
        if (dbIndex != dbIndex.init)
            return;

        dbIndex = pileUpDb.readRecord!PileUpDbIndex();
    }

    private PileUp[] readSlice(size_t from, size_t to)
    {
        assert(from <= to && to <= length);

        if (from == to)
            return [];

        // Step 1: determine memory requirements and DB slices
        dbSlices = getDbSlices(from, to);

        // Step 2: allocate minimally initialized memory for all blocks
        auto pileUps = minimallyInitializedArray!(PileUp[])(dbSlices.pileUps.length);
        auto readAlignments = minimallyInitializedArray!(ReadAlignment[])(dbSlices.readAlignments.length);
        auto seededAlignments = minimallyInitializedArray!(SeededAlignment[])(dbSlices.seededAlignments.length);
        auto localAlignments = minimallyInitializedArray!(LocalAlignment[])(dbSlices.localAlignments.length);
        auto tracePoints = minimallyInitializedArray!(TracePoint[])(dbSlices.tracePoints.length);

        // Step 3: parse each record for each block assigning already
        //         allocated array slices to the array fields
        parse(pileUps, readAlignments);
        parse(seededAlignments, localAlignments);
        parse(readAlignments, seededAlignments);
        parse(localAlignments, tracePoints);
        parse(tracePoints);

        return pileUps;
    }

    private DbSlices getDbSlices(size_t from, size_t to)
    {
        auto pileUps = dbIndex.pileUps[from .. to];
        auto firstPileUp = pileUpDb.readRecordAt!(StorageType!PileUp)(pileUps[0]);
        auto lastPileUp = pileUpDb.readRecordAt!(StorageType!PileUp)(pileUps[$ - 1]);

        auto readAlignments = ArrayStorage!(StorageType!ReadAlignment).fromPtrs(
            firstPileUp[0],
            lastPileUp[$],
        );
        auto firstReadAlignment = pileUpDb.readRecordAt!(StorageType!ReadAlignment)(readAlignments[0]);
        auto lastReadAlignment = pileUpDb.readRecordAt!(StorageType!ReadAlignment)(readAlignments[$ - 1]);

        auto seededAlignments = ArrayStorage!(StorageType!SeededAlignment).fromPtrs(
            firstReadAlignment[0],
            lastReadAlignment[$],
        );
        auto firstSeededAlignment = pileUpDb.readRecordAt!(StorageType!SeededAlignment)(seededAlignments[0]);
        auto lastSeededAlignment = pileUpDb.readRecordAt!(StorageType!SeededAlignment)(seededAlignments[$ - 1]);

        auto localAlignments = ArrayStorage!(StorageType!LocalAlignment).fromPtrs(
            firstSeededAlignment.localAlignments[0],
            lastSeededAlignment.localAlignments[$],
        );
        auto firstLocalAlignment = pileUpDb.readRecordAt!(StorageType!LocalAlignment)(localAlignments[0]);
        auto lastLocalAlignment = pileUpDb.readRecordAt!(StorageType!LocalAlignment)(localAlignments[$ - 1]);

        auto tracePoints = ArrayStorage!(StorageType!TracePoint).fromPtrs(
            firstLocalAlignment.tracePoints[0],
            lastLocalAlignment.tracePoints[$],
        );

        return DbSlices(
            pileUps,
            readAlignments,
            seededAlignments,
            localAlignments,
            tracePoints,
        );
    }

    private void parse(ref PileUp[] pileUps, ReadAlignment[] readAlignments)
    {
        static assert(PileUp.sizeof == StorageType!PileUp.sizeof);
        pileUpDb.seek(dbSlices.pileUps.ptr);
        pileUps = pileUpDb.readRecords(pileUps);

        size_t[2] pileUpSlice;
        foreach (ref pileUp; pileUps)
        {
            auto pileUpStorage = *cast(StorageType!PileUp*) &pileUp;
            pileUpSlice[0] = pileUpSlice[1];
            pileUpSlice[1] += pileUpStorage.length;

            pileUp = readAlignments[pileUpSlice[0] .. pileUpSlice[1]];
        }
    }

    private void parse(
        ref SeededAlignment[] seededAlignments,
        AlignmentChain.LocalAlignment[] localAlignments
    )
    {
        // Parse `SeededAlignment`s
        alias Contig = AlignmentChain.Contig;
        pileUpDb.seek(dbSlices.seededAlignments.ptr);

        size_t[2] seededAlignmentsSlice;
        foreach (ref seededAlignment; seededAlignments)
        {
            auto seededAlignmentStorage = pileUpDb.readRecord!(StorageType!SeededAlignment);

            seededAlignmentsSlice[0] = seededAlignmentsSlice[1];
            seededAlignmentsSlice[1] += seededAlignmentStorage.localAlignments.length;

            seededAlignment = SeededAlignment(
                AlignmentChain(
                    seededAlignmentStorage.id,
                    Contig(
                        seededAlignmentStorage.contigAId,
                        seededAlignmentStorage.contigALength,
                    ),
                    Contig(
                        seededAlignmentStorage.contigBId,
                        seededAlignmentStorage.contigBLength,
                    ),
                    seededAlignmentStorage.flags,
                    localAlignments[seededAlignmentsSlice[0] .. seededAlignmentsSlice[1]],
                    seededAlignmentStorage.tracePointDistance,
                ),
                seededAlignmentStorage.seed,
            );
        }
    }

    private void parse(ref ReadAlignment[] readAlignments, SeededAlignment[] seededAlignments)
    {
        pileUpDb.seek(dbSlices.readAlignments.ptr);

        size_t[2] readAlignmentsSlice;
        foreach (ref readAlignment; readAlignments)
        {
            auto readAlignmentStorage = pileUpDb.readRecord!(StorageType!ReadAlignment);

            readAlignmentsSlice[0] = readAlignmentsSlice[1];
            readAlignmentsSlice[1] += readAlignmentStorage.length;

            readAlignment = ReadAlignment(
                seededAlignments[readAlignmentsSlice[0] .. readAlignmentsSlice[1]]
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
        pileUpDb.seek(dbSlices.localAlignments.ptr);

        size_t[2] localAlignmentsSlice;
        foreach (ref localAlignment; localAlignments)
        {
            auto localAlignmentStorage = pileUpDb.readRecord!(StorageType!LocalAlignment);

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
        pileUpDb.seek(dbSlices.tracePoints.ptr);
        tracePoints = pileUpDb.readRecords(tracePoints);
    }
}

void writePileUpsDb(in PileUp[] pileUps, in string dbFile)
{
    auto pileUpDb = File(dbFile, "wb");
    lockIfPossible(pileUpDb);

    writePileUpsDb(pileUps, pileUpDb);
}

void writePileUpsDb(in PileUp[] pileUps, File pileUpDb)
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

    PileUpDbIndex dbIndex = buildPileUpDbIndex(pileUps);

    pileUpDb.rawWrite([dbIndex]);

    writePileUpsDbBlock!PileUp(pileUpDb, pileUps, dbIndex);
    writePileUpsDbBlock!ReadAlignment(pileUpDb, pileUps, dbIndex);
    writePileUpsDbBlock!SeededAlignment(pileUpDb, pileUps, dbIndex);
    writePileUpsDbBlock!LocalAlignment(pileUpDb, pileUps, dbIndex);
    writePileUpsDbBlock!TracePoint(pileUpDb, pileUps, dbIndex);
}

unittest
{
    import dentist.util.tempfile : mkstemp;
    import std.file : remove;

    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

    auto pileUps = getPileUpsTestData();

    enum totalDbSize =
        PileUpDbIndex.sizeof +
        StorageType!PileUp.sizeof * numPileUps +
        StorageType!ReadAlignment.sizeof * numReadAlignments +
        StorageType!SeededAlignment.sizeof * numSeededAlignments +
        StorageType!LocalAlignment.sizeof * numLocalAlignments +
        StorageType!TracePoint.sizeof * numTracePoints;

    auto tmpDb = mkstemp("./.unittest-XXXXXX");
    scope (exit)
    {
        tmpDb.file.close();
        remove(tmpDb.name);
    }

    writePileUpsDb(pileUps, tmpDb.file);
    tmpDb.file.sync();

    assert(tmpDb.file.size == totalDbSize);

    tmpDb.file.rewind();
    auto pileUpDb = PileUpDb(tmpDb.file);

    assert(pileUpDb[] == pileUps);
}

private PileUpDbIndex buildPileUpDbIndex(in PileUp[] pileUps) nothrow pure
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = AlignmentChain.LocalAlignment.TracePoint;

    PileUpDbIndex dbIndex;

    dbIndex.beginPtr!PileUp = PileUpDbIndex.sizeof;
    dbIndex.endPtr!PileUp = StorageType!PileUp.sizeof * pileUps.length;
    foreach (ref pileUp; pileUps)
    {
        dbIndex.endPtr!ReadAlignment += StorageType!ReadAlignment.sizeof * pileUp.length;
        foreach (ref readAlignment; pileUp)
        {
            dbIndex.endPtr!SeededAlignment += StorageType!SeededAlignment.sizeof * readAlignment.length;
            foreach (ref seededAlignment; readAlignment[])
            {
                dbIndex.endPtr!LocalAlignment += StorageType!LocalAlignment.sizeof * seededAlignment.localAlignments.length;
                foreach (ref localAlignment; seededAlignment.localAlignments)
                {
                    dbIndex.endPtr!TracePoint += StorageType!TracePoint.sizeof * localAlignment.tracePoints.length;
                }
            }
        }
    }

    dbIndex.readAlignmentsPtr += dbIndex.pileUpsPtr;
    dbIndex.seededAlignmentsPtr += dbIndex.readAlignmentsPtr;
    dbIndex.localAlignmentsPtr += dbIndex.seededAlignmentsPtr;
    dbIndex.tracePointsPtr += dbIndex.localAlignmentsPtr;
    dbIndex.eofPtr += dbIndex.tracePointsPtr;

    return dbIndex;
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = AlignmentChain.LocalAlignment.TracePoint;

    auto pileUps = getPileUpsTestData();
    auto dbIndex = buildPileUpDbIndex(pileUps);

    assert(dbIndex.pileUpsPtr == PileUpDbIndex.sizeof);
    assert(dbIndex.readAlignmentsPtr ==
            dbIndex.pileUpsPtr +
            StorageType!PileUp.sizeof * numPileUps);
    assert(dbIndex.seededAlignmentsPtr ==
            dbIndex.readAlignmentsPtr +
            StorageType!ReadAlignment.sizeof * numReadAlignments);
    assert(dbIndex.localAlignmentsPtr ==
            dbIndex.seededAlignmentsPtr +
            StorageType!SeededAlignment.sizeof * numSeededAlignments);
    assert(dbIndex.tracePointsPtr ==
            dbIndex.localAlignmentsPtr +
            StorageType!LocalAlignment.sizeof * numLocalAlignments);
    assert(dbIndex.eofPtr ==
            dbIndex.tracePointsPtr +
            StorageType!TracePoint.sizeof * numTracePoints);
}

void writePileUpsDbBlock(T)(ref File pileUpDb, in PileUp[] pileUps, in PileUpDbIndex dbIndex)
    if (is(T == PileUp))
{
    auto readAlignments = ArrayStorage!(StorageType!ReadAlignment)(dbIndex.beginPtr!ReadAlignment);

    version (assert)
    {
        auto storedPileUps = ArrayStorage!(StorageType!PileUp)(dbIndex.beginPtr!PileUp);
        assert(storedPileUps.ptr == pileUpDb.tell());
    }
    foreach (ref pileUp; pileUps)
    {
        readAlignments.length = pileUp.length;

        pileUpDb.rawWrite([readAlignments]);

        readAlignments.ptr = readAlignments[$];

        version (assert)
        {
            ++storedPileUps.length;
            assert(storedPileUps[$] == pileUpDb.tell());
        }
    }
}

void writePileUpsDbBlock(T)(ref File pileUpDb, in PileUp[] pileUps, in PileUpDbIndex dbIndex)
    if (is(T == ReadAlignment))
{
    auto seededAlignments = ArrayStorage!(StorageType!SeededAlignment)(dbIndex.beginPtr!SeededAlignment);

    version (assert)
    {
        auto readAlignments = ArrayStorage!(StorageType!ReadAlignment)(dbIndex.beginPtr!ReadAlignment);
        assert(readAlignments.ptr == pileUpDb.tell());
    }
    foreach (ref pileUp; pileUps)
    {
        foreach (ref readAlignment; pileUp)
        {
            seededAlignments.length = readAlignment.length;

            pileUpDb.rawWrite([seededAlignments]);

            seededAlignments.ptr = seededAlignments[$];

            version (assert)
            {
                ++readAlignments.length;
                assert(readAlignments[$] == pileUpDb.tell());
            }
        }
    }
}

void writePileUpsDbBlock(T)(ref File pileUpDb, in PileUp[] pileUps, in PileUpDbIndex dbIndex)
    if (is(T == SeededAlignment))
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;

    auto localAlignments = ArrayStorage!(StorageType!LocalAlignment)(dbIndex.beginPtr!LocalAlignment);

    version (assert)
    {
        auto seededAlignments = ArrayStorage!(StorageType!SeededAlignment)(dbIndex.beginPtr!SeededAlignment);
        assert(seededAlignments.ptr == pileUpDb.tell());
    }
    foreach (ref pileUp; pileUps)
    {
        foreach (ref readAlignment; pileUp)
        {
            foreach (ref seededAlignment; readAlignment[])
            {
                localAlignments.length = seededAlignment.localAlignments.length;

                pileUpDb.rawWrite([StorageType!SeededAlignment(
                    seededAlignment.id,
                    seededAlignment.contigA.id,
                    seededAlignment.contigA.length,
                    seededAlignment.contigB.id,
                    seededAlignment.contigB.length,
                    seededAlignment.flags,
                    localAlignments,
                    seededAlignment.tracePointDistance,
                    seededAlignment.seed,
                )]);

                localAlignments.ptr = localAlignments[$];

                version (assert)
                {
                    ++seededAlignments.length;
                  assert(seededAlignments[$] == pileUpDb.tell());
              }
            }
        }
    }
}

void writePileUpsDbBlock(T)(ref File pileUpDb, in PileUp[] pileUps, in PileUpDbIndex dbIndex)
    if (is(T == AlignmentChain.LocalAlignment))
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

    auto tracePoints = ArrayStorage!(StorageType!TracePoint)(dbIndex.beginPtr!TracePoint);

    version (assert)
    {
        auto localAlignments = ArrayStorage!(StorageType!LocalAlignment)(dbIndex.beginPtr!LocalAlignment);
        assert(localAlignments.ptr == pileUpDb.tell());
    }
    foreach (ref pileUp; pileUps)
    {
        foreach (ref readAlignment; pileUp)
        {
            foreach (ref seededAlignment; readAlignment[])
            {
                foreach (ref localAlignment; seededAlignment.localAlignments)
                {
                    tracePoints.length = localAlignment.tracePoints.length;

                    pileUpDb.rawWrite([StorageType!LocalAlignment(
                        localAlignment.contigA.begin,
                        localAlignment.contigA.end,
                        localAlignment.contigB.begin,
                        localAlignment.contigB.end,
                        localAlignment.numDiffs,
                        tracePoints,
                    )]);

                    tracePoints.ptr = tracePoints[$];

                    version (assert)
                    {
                        ++localAlignments.length;
                        assert(localAlignments[$] == pileUpDb.tell());
                    }
                }
            }
        }
    }
}

void writePileUpsDbBlock(T)(ref File pileUpDb, in PileUp[] pileUps, in PileUpDbIndex dbIndex)
    if (is(T == AlignmentChain.LocalAlignment.TracePoint))
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = LocalAlignment.TracePoint;

    version (assert)
    {
        auto tracePoints = ArrayStorage!(StorageType!TracePoint)(dbIndex.beginPtr!TracePoint);
        assert(tracePoints.ptr == pileUpDb.tell());
    }
    foreach (ref pileUp; pileUps)
    {
        foreach (ref readAlignment; pileUp)
        {
            foreach (ref seededAlignment; readAlignment[])
            {
                foreach (ref localAlignment; seededAlignment.localAlignments)
                {
                    pileUpDb.rawWrite(cast(const(StorageType!TracePoint[])) localAlignment.tracePoints);

                    version (assert)
                    {
                        tracePoints.length += localAlignment.tracePoints.length;
                        assert(tracePoints[$] == pileUpDb.tell());
                    }
                }
            }
        }
    }
}

private struct PileUpDbIndex
{
    mixin DbIndex;

    private static template NextType(T)
    {
        static if (is(T == PileUp))
            alias NextType = ReadAlignment;
        else static if (is(T == ReadAlignment))
            alias NextType = SeededAlignment;
        else static if (is(T == SeededAlignment))
            alias NextType = AlignmentChain.LocalAlignment;
        else static if (is(T == AlignmentChain.LocalAlignment))
            alias NextType = AlignmentChain.LocalAlignment.TracePoint;
        else static if (is(T == AlignmentChain.LocalAlignment.TracePoint))
            alias NextType = EOF;
    }

    private static template fieldPtr(T)
    {
        static if (is(T == PileUp))
            alias fieldPtr = pileUpsPtr;
        else static if (is(T == ReadAlignment))
            alias fieldPtr = readAlignmentsPtr;
        else static if (is(T == SeededAlignment))
            alias fieldPtr = seededAlignmentsPtr;
        else static if (is(T == AlignmentChain.LocalAlignment))
            alias fieldPtr = localAlignmentsPtr;
        else static if (is(T == AlignmentChain.LocalAlignment.TracePoint))
            alias fieldPtr = tracePointsPtr;
        else static if (is(T == EOF))
            alias fieldPtr = eofPtr;
    }

    size_t pileUpsPtr;
    size_t readAlignmentsPtr;
    size_t seededAlignmentsPtr;
    size_t localAlignmentsPtr;
    size_t tracePointsPtr;
    size_t eofPtr;

    @property alias pileUps = arrayStorage!PileUp;
    @property alias readAlignments = arrayStorage!ReadAlignment;
    @property alias seededAlignments = arrayStorage!SeededAlignment;
    @property alias localAlignments = arrayStorage!(AlignmentChain.LocalAlignment);
    @property alias tracePoints = arrayStorage!(AlignmentChain.LocalAlignment.TracePoint);
}

unittest
{
    alias LocalAlignment = AlignmentChain.LocalAlignment;
    alias TracePoint = AlignmentChain.LocalAlignment.TracePoint;
    enum begin = 1;
    enum end = 2;
    enum modified = 3;


    {
        PileUpDbIndex dbIndex;

        dbIndex.pileUpsPtr = begin;
        dbIndex.readAlignmentsPtr = end;

        assert(dbIndex.beginPtr!PileUp == begin);
        assert(dbIndex.endPtr!PileUp == end);

        dbIndex.beginPtr!PileUp = modified;
        dbIndex.endPtr!PileUp = modified;

        assert(dbIndex.pileUpsPtr == modified);
        assert(dbIndex.readAlignmentsPtr == modified);
    }
    {
        PileUpDbIndex dbIndex;

        dbIndex.readAlignmentsPtr = begin;
        dbIndex.seededAlignmentsPtr = end;

        assert(dbIndex.beginPtr!ReadAlignment == begin);
        assert(dbIndex.endPtr!ReadAlignment == end);

        dbIndex.beginPtr!ReadAlignment = modified;
        dbIndex.endPtr!ReadAlignment = modified;

        assert(dbIndex.readAlignmentsPtr == modified);
        assert(dbIndex.seededAlignmentsPtr == modified);
    }
    {
        PileUpDbIndex dbIndex;

        dbIndex.seededAlignmentsPtr = begin;
        dbIndex.localAlignmentsPtr = end;

        assert(dbIndex.beginPtr!SeededAlignment == begin);
        assert(dbIndex.endPtr!SeededAlignment == end);

        dbIndex.beginPtr!SeededAlignment = modified;
        dbIndex.endPtr!SeededAlignment = modified;

        assert(dbIndex.seededAlignmentsPtr == modified);
        assert(dbIndex.localAlignmentsPtr == modified);
    }
    {
        PileUpDbIndex dbIndex;

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
        PileUpDbIndex dbIndex;

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
    static if (is(T == PileUp))
        alias StorageType = ArrayStorage!(StorageType!ReadAlignment);
    else static if (is(T == ReadAlignment))
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

/**
    This is the `processPileUps` command of DENTIST.

    Command_Summary:

    ---
    Filter and process pile ups into insertions for gap closing.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.processPileUps;

package(dentist) enum summary = "
    Filter and process pile ups into insertions for gap closing.
";

import dentist.commandline : OptionsFor;
import dentist.commands.collectPileUps.filter : filterContainedAlignmentChains;
import dentist.commands.processPileUps.cropper : CropOptions, cropPileUp;
import dentist.common :
    dentistEnforce,
    DentistException,
    id_t,
    ReferenceInterval,
    ReferencePoint,
    ReferenceRegion;
import dentist.common.alignments :
    AlignmentChain,
    AlignmentLocationSeed,
    ChainingOptions,
    Contig,
    contigs,
    coord_t,
    getAlignmentRefs,
    getType,
    isExtension,
    isGap,
    makeJoin,
    PileUp,
    pileUpToSimpleJson,
    ReadAlignment,
    SeededAlignment,
    toChar;
import dentist.common.binio :
    CompressedSequence,
    InsertionDb,
    PileUpDb;
import dentist.common.commands : DentistCommand;
import dentist.common.insertions :
    Insertion,
    InsertionInfo;
import dentist.common.scaffold :
    ContigNode,
    getDefaultJoin,
    isParallel;
import dentist.util.containers : hashSet, HashSet;
import dentist.util.log;
import dentist.util.math : absdiff;
import dentist.dazzler :
    computeQVs,
    dbdust,
    dbEmpty,
    dbSubset,
    chainLocalAlignments,
    DBdumpOptions,
    DbRecord,
    filterPileUpAlignments,
    filterLocalAlignments,
    getAlignments,
    getDalignment,
    getDbRecords,
    getConsensus,
    getFastaSequence,
    lasEmpty,
    minQVCoverage,
    readMask,
    writeMask;
import std.algorithm :
    canFind,
    count,
    countUntil,
    cumulativeFold,
    equal,
    filter,
    find,
    joiner,
    map,
    max,
    maxElement,
    mean,
    min,
    merge,
    sort,
    swap,
    uniq;
import std.array : array, minimallyInitializedArray;
import std.conv : to;
import std.file : exists;
import std.format : format;
import std.parallelism : parallel, taskPool;
import std.path : buildPath;
import std.range :
    chain,
    drop,
    enumerate,
    evenChunks,
    iota,
    only,
    retro,
    zip;
import std.range.primitives : empty, front, popFront;
import std.typecons : tuple, Tuple, Yes;
import vibe.data.json : toJson = serializeToJson;


/// Options for the `processPileUps` command.
alias Options = OptionsFor!(DentistCommand.processPileUps);

/// Execute the `processPileUps` command with `options`.
void execute(in Options options)
{
    auto processor = new PileUpsProcessor(options);

    processor.run();
}

/// This class comprises the `processPileUps` step of the DENTIST algorithm
class PileUpsProcessor
{
    protected const Options options;
    protected PileUp[] pileUps;
    ReferenceRegion repeatMask;
    Insertion[] insertions;

    this(in ref Options options)
    {
        this.options = options;
        this.pileUps.reserve(options.numPileUps);
        this.insertions = minimallyInitializedArray!(Insertion[])(options.numPileUps);
    }

    void run()
    {
        mixin(traceExecution);

        readPileUps();
        readRepeatMask();

        foreach (i, pileUp; parallel(pileUps))
            processPileUp(i, pileUp);

        insertions.sort();
        dropEmptyInsertions();
        writeInsertions();
    }

    protected void processPileUp(size_t i, PileUp pileUp)
    {
        auto processor = new PileUpProcessor(options, repeatMask);

        processor.run(i, pileUp, &insertions[i]);
    }

    protected void readPileUps()
    {
        mixin(traceExecution);

        auto pileUpDb = PileUpDb.parse(options.pileUpsFile);

        foreach (pileUpBatch; options.pileUpBatches)
            pileUps ~= pileUpDb[pileUpBatch[0] .. pileUpBatch[1]];
    }

    protected void readRepeatMask()
    {
        mixin(traceExecution);

        foreach (mask; options.repeatMasks)
            repeatMask |= ReferenceRegion(readMask!ReferenceInterval(
                options.refDb,
                mask,
            ));
    }

    protected void dropEmptyInsertions()
    {
        insertions = insertions.find!(ins => ins.start.contigId != 0);
    }

    protected void writeInsertions()
    {
        mixin(traceExecution);

        // Use this output to quickly generate test data
        debug if (shouldLog(LogLevel.debug_))
            printInsertions(insertions);

        InsertionDb.write(options.insertionsFile, insertions);
    }
}


alias Enumerated(T) = Tuple!(
    size_t, "index",
    T, "value",
);


/// This class processes a single pileup.
protected class PileUpProcessor
{
    const(Options) options;
    const(ReferenceRegion) originalRepeatMask;
    ReferenceRegion repeatMask;

    protected const id_t[] pileUpIdMapping;
    protected id_t pileUpId;
    protected PileUp pileUp;
    protected Insertion* resultInsertion;
    protected string croppedDb;
    protected Contig[] pileUpContigs;
    protected ReferencePoint[] croppingPositions;
    protected AlignmentLocationSeed[] croppingSeeds;
    protected string pileUpAlignment;
    protected HashSet!id_t allowedReferenceReadIds;
    protected id_t[] referenceReadCandidateIndices;
    protected size_t referenceReadIdx;
    protected id_t referenceReadTry;
    protected string consensusDb;
    protected AlignmentChain[] postConsensusAlignment;
    protected ReadAlignment insertionAlignment;
    protected CompressedSequence insertionSequence;
    protected Insertion insertion;

    this(in Options options, in ReferenceRegion repeatMask)
    {
        this.options = options;
        this.originalRepeatMask = repeatMask;
        this.pileUpIdMapping = options
            .pileUpBatches
            .map!(batch => iota(batch[0], batch[1]))
            .joiner
            .array;
    }

    void run(size_t pileUpIdx, PileUp pileUp, Insertion* resultInsertion)
    {
        mixin(traceExecution);

        this.pileUpId = pileUpIdMapping[pileUpIdx];
        this.pileUp = pileUp;
        this.resultInsertion = resultInsertion;
        this.pileUpContigs = pileUp.contigs();

        if (
            (!options.onlyFlags.extending && pileUp.isExtension) ||
            (!options.onlyFlags.spanning && pileUp.isGap)
        )
        {
            logJsonWarn(
                "event", "pileUpSkipped",
                "info", "skipping pile up due to --only",
                "pileUpId", this.pileUpId,
                "pileUp", pileUp.pileUpToSimpleJson(),
            );

            return;
        }

        logJsonDiagnostic(
            "info", "processing pile up",
            "pileUpId", this.pileUpId,
            "pileUp", pileUp.pileUpToSimpleJson(),
        );

        processPileUp();
    }

    protected void processPileUp()
    {
        mixin(traceExecution);

        try
        {
            if (shouldSkipSmallPileUp())
                return;

            if (shouldProcessSingularPileUp())
            {
                postConsensusAlignment = pileUp[0][].map!"a.alignment".array;

                logJsonInfo(
                    "event", "singularPileUp",
                    "info", "using single read instead of pile up",
                    "pileUpId", pileUpId,
                    "pileUp", pileUp.pileUpToSimpleJson,
                    "readId", pileUp[0][0].contigB.id,
                );
            }
            else
            {
                reduceRepeatMaskToFlankingContigs();
                crop();
                adjustRepeatMaskToMakeMappingPossible();
                selectAllowedReferenceReadIds();
                computeQVs();
                findReferenceReadCandidates();
                while (selectReferenceRead(referenceReadTry++))
                {
                    try
                    {
                        computeConsensus();
                        break;
                    }
                    catch (Exception e)
                    {
                        logJsonDiagnostic(
                            "event", "consensusFailed",
                            "info", "computing reference-based consensus failed",
                            "reason", "error",
                            "error", e.message.to!string,
                            "pileUpId", pileUpId,
                            "pileUp", pileUp.pileUpToSimpleJson,
                            "referenceReadIdx", referenceReadIdx,
                        );
                    }
                }

                dentistEnforce(
                    referenceReadIdx < size_t.max,
                    "no valid reference read found",
                    [
                        "referenceReadCandidateIds": referenceReadCandidateIndices
                            .map!(idx => pileUp[idx][0].contigB.id)
                            .array,
                    ].toJson
                );

                alignConsensusToFlankingContigs();
            }

            getInsertionAlignment();
            getInsertionSequence();

            *resultInsertion = makeInsertion();
        }
        catch(DentistException e)
        {
            logJsonWarn(
                "event", "pileUpSkipped",
                "info", "skipping pile up due to errors",
                "reason", "error",
                "error", e.message.to!string,
                "errorPayload", e.payload,
                "pileUpId", pileUpId,
                "pileUp", pileUp.pileUpToSimpleJson,
            );
        }
        catch(Exception e)
        {
            logJsonWarn(
                "event", "pileUpSkipped",
                "info", "skipping pile up due to errors",
                "reason", "error",
                "error", e.message.to!string,
                "pileUpId", pileUpId,
                "pileUp", pileUp.pileUpToSimpleJson,
            );
        }
    }

    protected bool shouldProcessSingularPileUp() const nothrow
    {
        return options.allowSingleReads && pileUp.length == 1;
    }

    protected bool shouldSkipSmallPileUp() const nothrow
    {
        if (pileUp.length < options.minReadsPerPileUp && !shouldProcessSingularPileUp())
        {
            logJsonWarn(
                "event", "pileUpSkipped",
                "info", "skipping pile up due to `minReadsPerPileUp`",
                "reason", "minReadsPerPileUp",
                "pileUpId", pileUpId,
                "pileUp", pileUp.pileUpToSimpleJson,
            );

            return true;
        }

        return false;
    }

    protected void reduceRepeatMaskToFlankingContigs()
    {
        auto pileUpContigsRegion = ReferenceRegion(
            pileUpContigs
                .map!(contig => ReferenceInterval(contig.id, 0, contig.length))
                .array
        );
        repeatMask = originalRepeatMask & pileUpContigsRegion;
    }

    protected void crop()
    {
        auto croppingResult = cropPileUp(pileUp, repeatMask, CropOptions(
            options.refDb,
            options.readsDb,
            options.minAnchorLength,
            options.tmpdir,
        ));

        croppedDb = croppingResult.db;
        croppingPositions = croppingResult.referencePositions;
        croppingSeeds = croppingResult.seeds;

        if (pileUpContigs[0].id != croppingPositions[0].contigId)
            swap(pileUpContigs[0], pileUpContigs[1]);
    }

    protected void adjustRepeatMaskToMakeMappingPossible()
    {
        auto croppings = zip(pileUpContigs, croppingPositions, croppingSeeds)
            .map!(pair => tuple!(
                "contig",
                "position",
                "seed",
            )(
                pair[0],
                pair[2] == AlignmentLocationSeed.front
                    ? max(pair[1].value, options.minAnchorLength)
                    : min(pair[1].value, pair[0].length - options.minAnchorLength),
                pair[2],
            ));


        foreach (cropping; croppings)
        {
            auto contigInterval = ReferenceInterval(cropping.contig.id, 0, cropping.contig.length);
            auto croppedInterval = contigInterval;

            if (cropping.seed == AlignmentLocationSeed.front)
                croppedInterval.end = cropping.position;
            else
                croppedInterval.begin = cropping.position;

            auto localInvertedMask = ReferenceRegion(croppedInterval) - repeatMask;
            auto numUnmaskedBps = localInvertedMask.size;

            // Remove mask from contig if insufficient number of anchor bps
            if (numUnmaskedBps < options.minAnchorLength)
                repeatMask -= contigInterval;
        }
    }

    protected void selectAllowedReferenceReadIds()
    {
        // NOTE pileUp is not modified but the read alignments need to be assignable.
        allowedReferenceReadIds = hashSet(
            (cast(PileUp) pileUp)
                .enumerate
                .filter!(enumRA =>
                    enumRA.value.length == croppingPositions.length &&
                    enumRA.value[].map!"a.contigA.id".equal(croppingPositions.map!"a.contigId"))
                .map!(enumRA => cast(id_t) (enumRA.index + 1))
        );
    }

    protected void computeQVs()
    {
        dbdust(croppedDb, options.pileUpDustOptions);

        auto rawPileUpAlignment = getDalignment(
            croppedDb,
            options.pileUpAlignmentOptions,
            options.tmpdir,
        );
        auto hqPileUpAlignment = filterLocalAlignments!(
            a => a.averageErrorRate <= options.maxAlignmentError
        )(croppedDb, rawPileUpAlignment);

        dentistEnforce(
            !lasEmpty(hqPileUpAlignment),
            "empty pileup alignment",
        );

        auto chainedPileUpAlignment = chainLocalAlignments(
            croppedDb,
            hqPileUpAlignment,
            options.chainingOptions,
        );

        auto coverage = cast(id_t) allowedReferenceReadIds.size;

        if (coverage < minQVCoverage && pileUp.length >= minQVCoverage)
            coverage = minQVCoverage;

        .computeQVs(croppedDb, chainedPileUpAlignment, coverage);

        pileUpAlignment = filterPileUpAlignments(
            croppedDb,
            chainedPileUpAlignment,
            options.properAlignmentAllowance,
            Yes.forceFlat,
        );

        dentistEnforce(
            !lasEmpty(pileUpAlignment),
            "empty pileup alignment after filtering",
        );
    }

    protected void findReferenceReadCandidates()
    {
        auto croppedReads = getDbRecords(croppedDb, [
            DBdumpOptions.readNumber,
            DBdumpOptions.intrinsicQualityVector,
        ]).filter!(read => read.readNumber in allowedReferenceReadIds).array;

        auto hist = new size_t[DbRecord.maxQV];
        size_t histTotal;
        foreach (read; croppedReads)
            foreach (qv; read.intrinsicQVs)
                if (qv < DbRecord.maxQV)
                {
                    ++hist[qv];
                    ++histTotal;
                }

        auto badThres = cast(size_t) (options.badFraction * histTotal);
        auto badQV = DbRecord.maxQV - 1 - hist
            .retro
            .cumulativeFold!"a + b"(0UL)
            .countUntil!"a >= b"(badThres);

        logJsonInfo(
            "badQV", badQV,
            "pileUpId", pileUpId,
            "pileUp", pileUp.pileUpToSimpleJson,
        );

        auto readsWithScores = croppedReads.map!(read => tuple!(
            "numBadQVs",
            "meanQV",
            "readId",
        )(
            read.intrinsicQVs.count!"a >= b"(badQV),
            read.intrinsicQVs.mean,
            read.readNumber,
        )).array;

        sort(readsWithScores);

        logJsonDebug(
            "readsWithScores", readsWithScores.toJson,
            "pileUpId", pileUpId,
            "pileUp", pileUp.pileUpToSimpleJson,
        );

        referenceReadCandidateIndices = readsWithScores
            .map!(rws => cast(id_t) (rws.readId - 1))
            .array;
    }

    protected bool selectReferenceRead(in id_t referenceReadTry)
    {
        referenceReadIdx = bestReferenceReadIndex(referenceReadTry);

        if (referenceReadIdx == size_t.max)
            return false;

        logJsonDiagnostic(
            "referenceReadIdx", referenceReadIdx,
            "pileUpId", pileUpId,
            "pileUp", pileUp.pileUpToSimpleJson,
        );
        assert(referenceRead.length == croppingPositions.length);

        return true;
    }

    protected @property inout(ReadAlignment) referenceRead() inout
    {
        return pileUp[referenceReadIdx];
    }

    protected size_t bestReferenceReadIndex(in id_t skip) const pure
    {
        if (skip >= referenceReadCandidateIndices.length)
            return size_t.max;

        return referenceReadCandidateIndices[skip];
    }

    protected void computeConsensus()
    {
        mixin(traceExecution);

        consensusDb = getConsensus(
            croppedDb,
            pileUpAlignment,
            referenceReadIdx + 1,
            options.consensusOptions,
        );

        dentistEnforce(
            !dbEmpty(consensusDb),
            "consensus could not be computed",
            [
                "consensusDb": consensusDb.toJson,
                "referenceRead": (referenceReadIdx + 1).toJson,
            ].toJson,
        );
    }

    protected void alignConsensusToFlankingContigs()
    {
        mixin(traceExecution);

        auto flankingContigIds = croppingPositions.map!"a.contigId".array;
        auto flankingContigIdsWithSeeds = zip(flankingContigIds, croppingSeeds)
            .map!(pair => tuple(
                pair[0],
                pair[1].toChar
            ));
        auto flankingContigsDbName = buildPath(
            options.tmpdir,
            format!"contigs-%-(%(%d%c%)-%).dam"(flankingContigIdsWithSeeds),
        );
        auto flankingContigsDb = dbSubset(
            flankingContigsDbName,
            options.refDb,
            flankingContigIds,
            options.consensusOptions,
        );
        auto flankingContigsRepeatMask = repeatMask
            .intervals
            .filter!(interval => flankingContigIds.canFind(interval.contigId))
            .map!(interval => ReferenceInterval(
                1 + flankingContigIds.countUntil(interval.contigId),
                interval.begin,
                interval.end,
            ))
            .array;
        writeMask(
            flankingContigsDb,
            options.flankingContigsRepeatMaskName,
            flankingContigsRepeatMask,
        );
        dbdust(flankingContigsDb, options.consensusOptions.dbdustOptions);

        postConsensusAlignment = getAlignments(
            flankingContigsDb,
            consensusDb,
            getDalignment(
                flankingContigsDb,
                consensusDb,
                options.postConsensusAlignmentOptions,
                options.tmpdir,
            ),
            Yes.includeTracePoints,
        );

        postConsensusAlignment = filterContainedAlignmentChains(postConsensusAlignment);

        foreach (ref ac; postConsensusAlignment)
        {
            alias leftComplementMismatch = () =>
                ac.contigA.id == referenceRead[0].contigA.id &&
                ac.flags.complement != referenceRead[0].flags.complement;
            alias rightComplementMismatch = () =>
                croppingPositions.length > 1 &&
                ac.contigA.id == referenceRead[1].contigA.id &&
                ac.flags.complement != referenceRead[1].flags.complement;

            // Insert correct `contigId`s
            ac.contigA.id = croppingPositions[ac.contigA.id - 1]
                .contigId
                .to!id_t;
            ac.disableIf(
                !ac.isProper(options.properAlignmentAllowance) ||
                leftComplementMismatch() ||
                rightComplementMismatch()
            );
        }

        dentistEnforce(
            postConsensusAlignment.canFind!"!a.flags.disabled",
            "consensus does not align properly to flanking contig(s)",
            ["consensusDb": consensusDb].toJson,
        );
    }

    protected void getInsertionAlignment()
    {
        mixin(traceExecution);

        SeededAlignment[2] insertionAlignmentBuffer;

        foreach (ref croppingPos; croppingPositions)
        {
            alias contigsMatch = ac => ac.contigA.id == croppingPos.contigId;

            auto refReadFlankAlignmentIdx = referenceRead[].countUntil!contigsMatch;
            assert(refReadFlankAlignmentIdx >= 0);
            auto alignmentSeed = referenceRead[refReadFlankAlignmentIdx].seed;

            alias isProperInsertionOverlap = ac =>
                alignmentSeed == AlignmentLocationSeed.front
                    ? ac.first.contigA.begin <= options.properAlignmentAllowance &&
                      absdiff(ac.last.contigB.end, ac.contigB.length) <= options.properAlignmentAllowance
                    : absdiff(ac.last.contigA.end, ac.contigA.length) <= options.properAlignmentAllowance &&
                      ac.first.contigB.begin <= options.properAlignmentAllowance;

            auto flankAlignments = postConsensusAlignment
                .filter!(ac => !ac.flags.disabled)
                .filter!contigsMatch
                .filter!isProperInsertionOverlap;

            dentistEnforce(
                !flankAlignments.empty,
                format!"consensus does not align to flanking contig %d"(croppingPos.contigId),
                ["consensusDb": consensusDb].toJson,
            );

            auto flankAlignment = SeededAlignment(flankAlignments.front, alignmentSeed);
            insertionAlignmentBuffer[refReadFlankAlignmentIdx] = flankAlignment;

            flankAlignments.popFront();
            dentistEnforce(
                flankAlignments.empty,
                format!"consensus ambiguously aligns to flanking contig %d"(croppingPos.contigId),
                [
                    "consensusDb": consensusDb.toJson,
                    "alignments": chain([flankAlignment.alignment], flankAlignments)
                        .array
                        .toJson,
                ].toJson,
            );
        }

        insertionAlignment = ReadAlignment(insertionAlignmentBuffer[0 .. referenceRead.length]);

        dentistEnforce(
            insertionAlignment.isValid,
            "consensus alignment is invalid",
            ["consensusDb": consensusDb].toJson,
        );
        dentistEnforce(
            (
                insertionAlignment.type == referenceRead.type &&
                insertionAlignment.isParallel == referenceRead.isParallel
            ),
            format!"consensus alignment has an unexpected type: %s%s"(
                insertionAlignment.type,
                insertionAlignment.isGap
                    ? insertionAlignment.isParallel
                        ? " (parallel)"
                        : " (anti-parallel)"
                    : "",
            ),
            ["consensusDb": consensusDb].toJson,
        );
    }

    protected void getInsertionSequence()
    {
        if (shouldProcessSingularPileUp())
        {
            auto fastaSequence = getFastaSequence(
                options.readsDb,
                pileUp[0][0].contigB.id,
            );
            insertionSequence = CompressedSequence.from(fastaSequence);
        }
        else
        {
            auto fastaSequence = getFastaSequence(consensusDb, 1);
            insertionSequence = CompressedSequence.from(fastaSequence);
            assert(fastaSequence.length == insertionSequence.length);
        }
    }

    protected Insertion makeInsertion()
    {
        auto insertion = makeJoin!Insertion(referenceRead);
        insertion.payload = InsertionInfo(
            insertionSequence,
            0,
            insertionAlignment[],
            pileUp.map!"a[0].contigB.id".array.sort.release,
        );

        assert(insertion.isParallel == insertionAlignment.isParallel);
        if (insertionAlignment.length == 2)
            assert(insertionAlignment[0].contigB == insertionAlignment[1].contigB);
        assert(insertionSequence.length == insertionAlignment[0].contigB.length);

        return insertion;
    }
}

debug private void printInsertions(in Insertion[] insertions)
{
    import std.stdio : writefln;

    foreach (insertion; insertions)
    {
        writefln!"Insertion(";
        writefln!"    ContigNode(%d, %s),"(insertion.start.contigId, insertion.start.contigPart.to!string);
        writefln!"    ContigNode(%d, %s),"(insertion.end.contigId, insertion.end.contigPart.to!string);
        writefln!"    InsertionInfo(";
        writefln!`        CompressedSequence.from("%s"),`(insertion.payload.sequence.to!string);
        writefln!"        %d,"(insertion.payload.contigLength);
        writefln!"        [";

        foreach (overlap; insertion.payload.overlaps)
        {

        writefln!"            SeededAlignment(";
        writefln!"                AlignmentChain(";
        writefln!"                    %d,"(overlap.id);
        writefln!"                    Contig(%d, %d),"(overlap.contigA.id, overlap.contigA.length);
        writefln!"                    Contig(%d, %d),"(overlap.contigB.id, overlap.contigB.length);
        writefln!"                    %s,"(overlap.flags.complement ? "Flags(complement)" : "emptyFlags");
        writefln!"                    [";

        foreach (localAlignment; overlap.localAlignments)
        {
        writefln!"                        LocalAlignment(";
        writefln!"                            Locus(%d, %d),"(localAlignment.contigA.begin, localAlignment.contigA.end);
        writefln!"                            Locus(%d, %d),"(localAlignment.contigB.begin, localAlignment.contigB.end);
        writefln!"                            %d,"(localAlignment.numDiffs);
        writefln!"                            [";

        foreach (tracePoint; localAlignment.tracePoints)
        {
        writefln!"                                TracePoint(%d, %d),"(tracePoint.numDiffs, tracePoint.numBasePairs);
        }

        writefln!"                            ],";
        writefln!"                        ),";
        }

        writefln!"                    ],";
        writefln!"                    %d,"(overlap.tracePointDistance);
        writefln!"                ),";
        writefln!"                AlignmentLocationSeed.%s,"(overlap.seed.to!string);
        writefln!"            ),";

        }

        writefln!"        ],";
        writefln!"        [";

        foreach (readId; insertion.payload.readIds)
        {
        writefln!"            %d"(readId);
        }

        writefln!"        ],";
        writefln!"    ),";
        writefln!"),";
    }
}

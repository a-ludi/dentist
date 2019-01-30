/**
    This is the `processPileUps` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.processPileUps;

import dentist.commandline : DentistCommand, OptionsFor;
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
    getAlignmentRefs,
    getType,
    isExtension,
    makeJoin,
    PileUp,
    pileUpToSimpleJson,
    ReadAlignment,
    SeededAlignment;
import dentist.common.binio :
    CompressedSequence,
    InsertionDb,
    PileUpDb;
import dentist.common.insertions :
    Insertion,
    InsertionInfo,
    SpliceSite;
import dentist.common.scaffold : ContigNode, getDefaultJoin;
import dentist.util.log;
import dentist.dazzler :
    attachTracePoints,
    dbEmpty,
    dbSubset,
    getAlignments,
    getDalignment,
    getConsensus,
    getFastaSequence,
    readMask;
import std.algorithm :
    canFind,
    countUntil,
    equal,
    filter,
    find,
    joiner,
    map,
    maxElement,
    merge,
    sort,
    uniq;
import std.array : array;
import std.conv : to;
import std.format : format;
import std.parallelism : parallel, taskPool;
import std.range : enumerate, evenChunks, only, zip;
import std.range.primitives : empty, front, popFront;
import std.typecons : Yes;
import vibe.data.json : toJson = serializeToJson;


/// Options for the `processPileUps` command.
alias Options = OptionsFor!(DentistCommand.processPileUps);

/// Execute the `processPileUps` command with `options`.
void execute(Options)(in Options options)
{
    auto processor = new PileUpsProcessor(options);

    processor.run();
}

/// This class comprises the `processPileUps` step of the `dentist` algorithm
class PileUpsProcessor
{
    protected const Options options;
    protected PileUp[] pileUps;
    ReferenceRegion repeatMask;
    Insertion[] insertions;

    this(in ref Options options)
    {
        this.options = options;
        this.insertions.length = options.pileUpBatchSize;
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
        pileUps = pileUpDb[options.pileUpBatch[0] .. options.pileUpBatch[1]];
    }

    protected void readRepeatMask()
    {
        mixin(traceExecution);

        repeatMask = ReferenceRegion(readMask!ReferenceInterval(
            options.refDb,
            options.repeatMask,
            options.workdir,
        ));

        foreach (mask; options.additionalMasks)
            repeatMask |= ReferenceRegion(readMask!ReferenceInterval(
                options.refDb,
                mask,
                options.workdir,
            ));
    }

    protected void dropEmptyInsertions()
    {
        insertions = insertions.find!(ins => ins.start.contigId != 0);
    }

    protected void writeInsertions()
    {
        mixin(traceExecution);

        InsertionDb.write(options.insertionsFile, insertions);
    }
}

/// This class processes a single pileup.
protected class PileUpProcessor
{
    const(Options) options;
    const(ReferenceRegion) repeatMask;

    protected size_t pileUpId;
    protected PileUp pileUp;
    protected Insertion* resultInsertion;
    protected string croppedDb;
    protected ReferencePoint[] croppingPositions;
    protected size_t referenceReadIdx;
    protected string consensusDb;
    protected AlignmentChain[] postConsensusAlignment;
    protected ReadAlignment insertionAlignment;
    protected CompressedSequence insertionSequence;
    protected Insertion insertion;

    this(in Options options, in ReferenceRegion repeatMask)
    {
        this.options = options;
        this.repeatMask = repeatMask;
        this.pileUpId = 0 + pileUpId;
        this.pileUp = pileUp;
        this.resultInsertion = resultInsertion;
    }

    void run(size_t pileUpId, PileUp pileUp, Insertion* resultInsertion)
    {
        this.pileUpId = 0 + pileUpId;
        this.pileUp = fetchTracePoints(pileUp);
        this.resultInsertion = resultInsertion;

        processPileUp();
    }

    protected ref PileUp fetchTracePoints(ref PileUp pileUp)
    {
        mixin(traceExecution);

        auto allAlignmentChains = pileUp.getAlignmentRefs();
        allAlignmentChains.sort!"*a < *b";
        allAlignmentChains.attachTracePoints(
            options.refDb,
            options.readsDb,
            options.readsAlignmentFile,
            options.tracePointDistance,
            options.workdir
        );

        return pileUp;
    }

    protected void processPileUp()
    {
        mixin(traceExecution);

        try
        {
            if (shouldSkipSmallPileUp())
                return;

            crop();
            selectReferenceRead();
            computeConsensus();
            alignConsensusToFlankingContigs();
            getInsertionAlignment();
            getInsertionSequence();

            if (pileUp.isExtension && shouldSkipShortExtension())
                return;

            *resultInsertion = makeInsertion();
        }
        catch(DentistException e)
        {
            logJsonWarn(
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
                "info", "skipping pile up due to errors",
                "reason", "error",
                "error", e.message.to!string,
                "pileUpId", pileUpId,
                "pileUp", pileUp.pileUpToSimpleJson,
            );
        }
    }

    protected bool shouldSkipSmallPileUp() const nothrow
    {
        if (pileUp.length < options.minReadsPerPileUp)
        {
            logJsonInfo(
                "info", "skipping pile up due to `minReadsPerPileUp`",
                "reason", "minReadsPerPileUp",
                "pileUpId", pileUpId,
                "pileUp", pileUp.pileUpToSimpleJson,
            );

            return true;
        }

        return false;
    }

    protected void crop()
    {
        auto croppingResult = cropPileUp(pileUp, repeatMask, CropOptions(
            options.readsDb,
            options.tracePointDistance,
            options.workdir,
        ));

        croppedDb = croppingResult.db;
        croppingPositions = croppingResult.referencePositions;
    }

    protected void selectReferenceRead()
    {
        referenceReadIdx = bestReadAlignmentIndex(pileUp, croppingPositions);

        assert(referenceRead.length == croppingPositions.length);
    }

    protected @property inout(ReadAlignment) referenceRead() inout
    {
        return pileUp[referenceReadIdx];
    }

    protected size_t bestReadAlignmentIndex(
        in PileUp pileUp,
        in ReferencePoint[] referencePositions,
    ) const pure
    {
        // NOTE pileUp is not modified but the read alignments need to be assignable.
        return (cast(PileUp) pileUp)
            .enumerate
            .filter!(enumReadAlignment =>
                enumReadAlignment.value.length == referencePositions.length &&
                enumReadAlignment.value[].map!"a.contigA.id".equal(referencePositions.map!"a.contigId"))
            .maxElement!(enumReadAlignment => enumReadAlignment.value.meanScore)
            .index;
    }

    protected void computeConsensus()
    {
        consensusDb = getConsensus(
            croppedDb,
            referenceReadIdx + 1,
            options.consensusOptions,
        );

        dentistEnforce(
            !dbEmpty(consensusDb, options.workdir),
            "consensus could not be computed",
            [
                "consensusDb": consensusDb.toJson,
                "referenceRead": (referenceReadIdx + 1).toJson,
            ].toJson,
        );
    }

    protected void alignConsensusToFlankingContigs()
    {
        // FIXME crop flanks to expected matching region + a small margin
        auto flankingContigsDb = dbSubset(
            options.refDb,
            croppingPositions.map!"a.contigId",
            options.consensusOptions,
        );
        postConsensusAlignment = getAlignments(
            flankingContigsDb,
            consensusDb,
            getDalignment(
                flankingContigsDb,
                consensusDb,
                options.postConsensusAlignmentOptions,
                options.workdir,
            ),
            options.workdir,
            options.tracePointDistance,
        );

        foreach (ref ac; postConsensusAlignment)
        {
            // Insert correct `contigId`s
            ac.contigA.id = croppingPositions[ac.contigA.id - 1]
                .contigId
                .to!id_t;
            ac.disableIf(!ac.isProper);
        }

        dentistEnforce(
            postConsensusAlignment.canFind!"!a.flags.disabled",
            "consensus does not align to flanking contig(s)",
            ["consensusDb": consensusDb].toJson,
        );
    }

    protected void getInsertionAlignment()
    {
        SeededAlignment[2] insertionAlignmentBuffer;

        foreach (ref croppingPos; croppingPositions)
        {
            alias contigsMatch = ac => ac.contigA.id == croppingPos.contigId;

            auto refReadFlankAlignmentIdx = referenceRead[].countUntil!contigsMatch;
            assert(refReadFlankAlignmentIdx >= 0);
            auto alignmentSeed = referenceRead[refReadFlankAlignmentIdx].seed;

            alias isProperInsertionOverlap = ac =>
                alignmentSeed == AlignmentLocationSeed.front
                    ? ac.first.contigA.begin == 0 && ac.last.contigB.end == ac.contigB.length
                    : ac.last.contigA.end == ac.contigA.length && ac.first.contigB.begin == 0;

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
                ["consensusDb": consensusDb].toJson,
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
        auto fastaSequence = getFastaSequence(consensusDb, 1, options.workdir);
        insertionSequence = CompressedSequence.from(fastaSequence);
    }

    protected bool shouldSkipShortExtension() const
    {
        assert(croppingPositions.length == 1);

        auto refPos = croppingPositions[0].value;
        auto refLength = referenceRead[0].contigA.length;
        ulong extensionLength = referenceRead.isFrontExtension
            ? insertionSequence.length.to!ulong - refPos.to!ulong
            : insertionSequence.length.to!ulong - (refLength - refPos).to!ulong;

        if (extensionLength < options.minExtensionLength)
        {
            logJsonInfo(
                "info", "skipping pile up due to `minExtensionLength`",
                "reason", "minExtensionLength",
                "pileUpId", pileUpId,
                "pileUp", pileUp.pileUpToSimpleJson,
            );

            return true;
        }

        return false;
    }

    protected Insertion makeInsertion()
    {
        auto insertion = makeJoin!Insertion(referenceRead);
        // FIXME reduce payload to `insertionSequence`
        insertion.payload = InsertionInfo(
            insertionSequence,
            0,
            zip(croppingPositions, referenceRead[].map!"a.seed", referenceRead[].map!"a.flags")
                .map!(spliceSite => SpliceSite(spliceSite.expand))
                .array,
        );

        return insertion;
    }
}

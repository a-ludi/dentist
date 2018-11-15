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
    ReferenceInterval,
    ReferencePoint,
    ReferenceRegion;
import dentist.common.alignments :
    AlignmentChain,
    getAlignmentRefs,
    getType,
    isExtension,
    makeJoin,
    PileUp,
    ReadAlignment;
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
    getConsensus,
    getFastaSequences,
    readMask;
import std.algorithm :
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
import std.exception : enforce;
import std.parallelism : parallel, taskPool;
import std.range : enumerate, evenChunks, only, zip;
import std.typecons : Yes;
import vibe.data.json : toJson = serializeToJson;


/// Options for the `processPileUps` command.
alias Options = OptionsFor!(DentistCommand.processPileUps);

/// Execute the `processPileUps` command with `options`.
void execute(Options)(in Options options)
{
    auto processor = new PileUpProcessor(options);

    processor.run();
}

/// This class comprises the `processPileUps` step of the `dentist` algorithm
class PileUpProcessor
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
        {
            fetchTracePoints(pileUp);

            processPileUp(pileUp, &insertions[i]);
        }

        insertions.sort();
        dropEmptyInsertions();
        debug logJsonDebug("insertions", insertions
            .map!(ins => [
                "start": [
                    "contigId": ins.start.contigId.toJson,
                    "contigPart": ins.start.contigPart.to!string.toJson,
                ],
                "end": [
                    "contigId": ins.end.contigId.toJson,
                    "contigPart": ins.end.contigPart.to!string.toJson,
                ],
                "payload": [
                    "sequence": ins.payload.sequence.to!string.toJson,
                    "contigLength": ins.payload.contigLength.toJson,
                    "spliceSites": ins.payload.spliceSites.toJson,
                ],
            ])
            .array
            .toJson);

        writeInsertions();
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

    protected void processPileUp(PileUp pileUp, Insertion* resultInsertion) const
    {
        mixin(traceExecution);

        try
        {
            if (pileUp.length < options.minReadsPerPileUp)
            {
                return;
            }

            auto croppingResult = cropPileUp(pileUp, repeatMask, CropOptions(
                options.readsDb,
                options.tracePointDistance,
                options.workdir,
            ));
            auto referenceReadIdx = bestReadAlignmentIndex(pileUp, croppingResult.referencePositions);
            auto referenceRead = pileUp[referenceReadIdx];
            assert(referenceRead.length == croppingResult.referencePositions.length);

            auto consensusDb = getConsensus(
                croppingResult.db,
                referenceReadIdx + 1,
                options.consensusOptions
            );
            auto insertSequences = getFastaSequences(consensusDb, only(1), options.workdir);
            enforce!Exception(!insertSequences.empty, "consensus could not be computed");
            auto insertSequence = insertSequences.front;

            if (
                pileUp.isExtension &&
                shouldSkipShortExtension(croppingResult, insertSequence, referenceRead)
            )
            {
                return;
            }

            auto compressedSequence = CompressedSequence.from(insertSequence);
            *resultInsertion = makeInsertions(
                referenceRead,
                compressedSequence,
                croppingResult.referencePositions,
            );
        }
        catch(Exception e)
        {
            logJsonWarn(
                "info", "skipping pile up due to errors",
                "error", e.message().to!string,
                "pileUp", [
                    "type": pileUp.getType.to!string.toJson,
                    "contigIds": pileUp
                        .map!(ra => ra[].map!"a.contigA.id".array)
                        .joiner
                        .array
                        .sort
                        .uniq
                        .array
                        .toJson,
                ],
            );
        }
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

    protected bool shouldSkipShortExtension(T)(
        T croppingResult,
        in string insertSequence,
        in ReadAlignment referenceRead,
    ) const
    {
        assert(croppingResult.referencePositions.length == 1);

        auto refPos = croppingResult.referencePositions[0].value;
        auto refLength = referenceRead[0].contigA.length;
        ulong extensionLength = referenceRead.isFrontExtension
            ? insertSequence.length.to!ulong - refPos.to!ulong
            : insertSequence.length.to!ulong - (refLength - refPos).to!ulong;

        return extensionLength < options.minExtensionLength;
    }

    protected static Insertion makeInsertions(
        ReadAlignment referenceRead,
        CompressedSequence compressedSequence,
        ReferencePoint[] referencePositions,
    )
    {
        auto insertion = makeJoin!Insertion(referenceRead);
        insertion.payload = InsertionInfo(
            compressedSequence,
            0,
            zip(referencePositions, referenceRead[].map!"a.seed", referenceRead[].map!"a.flags")
                .map!(spliceSite => SpliceSite(spliceSite.expand))
                .array,
        );

        return insertion;
    }
}

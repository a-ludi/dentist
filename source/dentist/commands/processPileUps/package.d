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
    insertionScore,
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
    find,
    joiner,
    map,
    maxIndex,
    merge,
    sort,
    uniq;
import std.array : array;
import std.conv : to;
import std.exception : enforce;
import std.parallelism : parallel, taskPool;
import std.range : evenChunks, only, zip;
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
    static enum maxInsertionsPerPileup = 3;

    protected const Options options;
    protected PileUp[] pileUps;
    ReferenceRegion repeatMask;
    Insertion[] insertions;

    this(in ref Options options)
    {
        this.options = options;
        this.insertions.length = maxInsertionsPerPileup * options.pileUpBatchSize;
    }

    void run()
    {
        readPileUps();
        readRepeatMask();

        foreach (i, pileUp; parallel(pileUps))
        {
            fetchTracePoints(pileUp);

            auto baseIdx = maxInsertionsPerPileup * i;
            processPileUp(pileUp, insertions[baseIdx .. baseIdx + maxInsertionsPerPileup]);
        }

        import std.datetime.stopwatch;

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

        InsertionDb.write(options.insertionsFile, insertions);
    }

    protected void readPileUps()
    {
        auto pileUpDb = PileUpDb.parse(options.pileUpsFile);
        pileUps = pileUpDb[options.pileUpBatch[0] .. options.pileUpBatch[1]];
    }

    protected void readRepeatMask()
    {
        repeatMask = ReferenceRegion(readMask!ReferenceInterval(
            options.refDb,
            options.repeatMask,
            options.workdir,
        ));
    }

    protected ref PileUp fetchTracePoints(ref PileUp pileUp)
    {
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

    protected void processPileUp(PileUp pileUp, Insertion[] insertionsBuffer)
    {
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
            auto referenceReadIdx = bestReadAlignmentIndex(pileUp);
            auto referenceRead = pileUp[referenceReadIdx];
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
            insertionsBuffer[0] = makeInsertions(
                referenceRead,
                compressedSequence,
                croppingResult.referencePositions,
            );
            insertionsBuffer[1 .. 1 + croppingResult.referencePositions.length] =
                makeFlankingContigSlices(croppingResult.referencePositions)[];
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

    protected size_t bestReadAlignmentIndex(in PileUp pileUp) pure
    {
        auto pileUpType = pileUp.getType;

        return pileUp
            .map!(readAlignment => insertionScore(
                readAlignment,
                pileUpType,
                repeatMask,
                Yes.preferSpanning,
                options,
            ))
            .maxIndex;
    }

    protected bool shouldSkipShortExtension(T)(
        T croppingResult,
        in string insertSequence,
        in ReadAlignment referenceRead,
    )
    {
        assert(croppingResult.referencePositions.length == 1);

        auto refPos = croppingResult.referencePositions[0].value;
        auto refLength = referenceRead[0].contigA.length;
        ulong extensionLength = referenceRead.isFrontExtension
            ? insertSequence.length.to!ulong - refPos.to!ulong
            : insertSequence.length.to!ulong - (refLength - refPos).to!ulong;

        return extensionLength < options.minExtensionLength;
    }

    protected Insertion makeInsertions(
        ReadAlignment referenceRead,
        CompressedSequence compressedSequence,
        ReferencePoint[] referencePositions,
    )
    {
        auto insertion = makeJoin!Insertion(referenceRead);
        insertion.payload = InsertionInfo(
            compressedSequence,
            0,
            zip(referencePositions, referenceRead[].map!"a.flags")
                .map!(spliceSite => cast(SpliceSite) spliceSite)
                .array,
        );

        return insertion;
    }

    protected Insertion[] makeFlankingContigSlices(ref ReferencePoint[] referencePositions)
    {
        return referencePositions
            .map!(referencePosition => {
                auto contigEdge = getDefaultJoin!InsertionInfo(referencePosition.contigId);
                contigEdge.payload.spliceSites = [
                    SpliceSite(referencePosition, AlignmentChain.emptyFlags)
                ];

                return contigEdge;
            })
            .map!"a()"
            .array;
    }
}

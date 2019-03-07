/**
    This is the `output` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.output;

import dentist.commandline : OptionsFor;
import dentist.common :
    ReferenceInterval,
    ReferenceRegion,
    ReferencePoint;
import dentist.common.alignments :
    AlignmentChain,
    AlignmentLocationSeed,
    coord_t,
    id_t,
    SeededAlignment,
    trace_point_t;
import dentist.common.commands : DentistCommand;
import dentist.common.binio :
    CompressedSequence,
    InsertionDb;
import dentist.common.insertions :
    getCroppingPosition,
    getInfoForExistingContig,
    getInfoForGap,
    getInfoForNewSequenceInsertion,
    Insertion,
    InsertionInfo,
    isOutputGap,
    isValidInsertion,
    OutputScaffold;
import dentist.common.scaffold :
    ContigNode,
    ContigPart,
    scaffoldStarts,
    enforceJoinPolicy,
    getUnkownJoin,
    initScaffold,
    isAntiParallel,
    isDefault,
    isExtension,
    isFrontExtension,
    isGap,
    linearWalk,
    normalizeUnkownJoins,
    removeExtensions,
    removeSpanning;
import dentist.dazzler :
    ContigSegment,
    GapSegment,
    getNumContigs,
    getFastaSequence,
    getScaffoldStructure,
    ScaffoldSegment;
import dentist.util.algorithm : replaceInPlace;
import dentist.util.fasta : complement, reverseComplementer;
import dentist.util.log;
import dentist.util.math :
    absdiff,
    add,
    bulkAdd,
    ceil,
    filterEdges,
    floor,
    mean,
    RoundingMode;
import dentist.util.range : wrapLines;
import std.algorithm :
    all,
    canFind,
    copy,
    count,
    countUntil,
    filter,
    find,
    fold,
    joiner,
    map,
    max,
    maxElement,
    min,
    predSwitch,
    swap,
    swapAt;
import std.array : array;
import std.ascii : toUpper;
import std.conv : to;
import std.format : format;
import std.range :
    enumerate,
    dropExactly,
    only,
    repeat,
    takeExactly,
    walkLength;
import std.range.primitives : empty, front, popFront, save;
import std.stdio : File, stderr, stdout;
import std.typecons : Flag, No, tuple, Tuple, Yes;
import vibe.data.json : toJson = serializeToJson;


/// Options for the `output` command.
alias Options = OptionsFor!(DentistCommand.output);


/// Execute the `output` command with `options`.
void execute(Options)(in Options options)
{
    auto writer = new AssemblyWriter(options);

    writer.run();
}


class AssemblyWriter
{
    alias FastaWriter = typeof(wrapLines(stdout.lockingTextWriter, 0));

    protected const Options options;
    protected const(ScaffoldSegment)[] scaffoldStructure;
    size_t numReferenceContigs;
    size_t[] contigLengths;
    OutputScaffold assemblyGraph;
    OutputScaffold.IncidentEdgesCache incidentEdgesCache;
    ContigNode[] scaffoldStartNodes;
    File assemblyFile;
    FastaWriter writer;

    this(in ref Options options)
    {
        this.options = options;
        this.assemblyFile = options.assemblyFile is null
            ? stdout
            : File(options.assemblyFile, "w");
        this.writer = wrapLines(assemblyFile.lockingTextWriter, options.fastaLineWidth);
    }

    void run()
    {
        mixin(traceExecution);

        init();
        buildAssemblyGraph();
        scaffoldStartNodes = scaffoldStarts!InsertionInfo(assemblyGraph, incidentEdgesCache).array;
        logStatistics();

        foreach (startNode; scaffoldStartNodes)
            writeNewScaffold(startNode);
    }

    protected void init()
    {
        mixin(traceExecution);

        numReferenceContigs = getNumContigs(options.refDb, options.workdir);
        scaffoldStructure = getScaffoldStructure(options.refDb).array;
        contigLengths = scaffoldStructure
            .filter!(part => part.peek!ContigSegment !is null)
            .map!(contigPart => contigPart.get!ContigSegment)
            .map!(contigPart => contigPart.end - contigPart.begin + 0)
            .array;
    }


    protected void buildAssemblyGraph()
    {
        mixin(traceExecution);

        auto insertionDb = InsertionDb.parse(options.insertionsFile);
        auto insertions = insertionDb[];

        assemblyGraph = initScaffold!(
            contigId => InsertionInfo(CompressedSequence(), contigLengths[contigId - 1], []),
            InsertionInfo,
        )(numReferenceContigs);
        assemblyGraph.bulkAdd!(joins => mergeInsertions(joins))(insertions);
        appendUnkownJoins();

        if (options.onlyFlags.extending)
        {
            assemblyGraph.filterEdges!(insertion => skipShortExtension(insertion));
        }
        else
        {
            assemblyGraph = assemblyGraph.removeExtensions!InsertionInfo();
        }

        if (!options.onlyFlags.spanning)
        {
            assemblyGraph = assemblyGraph.removeSpanning!InsertionInfo();
        }

        assemblyGraph = assemblyGraph
            .enforceJoinPolicy!InsertionInfo(options.joinPolicy)
            .normalizeUnkownJoins!InsertionInfo()
            .fixCropping(options.tracePointDistance);
        incidentEdgesCache = assemblyGraph.allIncidentEdges();

        if (options.assemblyGraphFile !is null)
            InsertionDb.write(options.assemblyGraphFile, assemblyGraph.edges);
    }

    protected void appendUnkownJoins()
    {
        auto unkownJoins = scaffoldStructure[]
            .filter!(part => part.peek!GapSegment !is null)
            .map!(gapPart => gapPart.get!GapSegment)
            .map!(gapPart => getUnkownJoin(
                gapPart.beginGlobalContigId,
                gapPart.endGlobalContigId,
                InsertionInfo(CompressedSequence(), gapPart.length, []),
            ));
        assemblyGraph.bulkAddForce(unkownJoins);
    }

    protected Flag!"keepInsertion" skipShortExtension(Insertion insertion) const
    {
        assert(insertion.payload.overlaps.length == 1);

        auto extensionOverlap = insertion.payload.overlaps[0];
        auto extensionLength = extensionOverlap.seed == AlignmentLocationSeed.front
            ? extensionOverlap.first.contigB.begin
            : extensionOverlap.contigB.length - extensionOverlap.last.contigB.end;

        if (extensionLength < options.minExtensionLength)
        {
            logJsonInfo(
                "info", "skipping pile up due to `minExtensionLength`",
                "reason", "minExtensionLength",
                "flankingContigId", extensionOverlap.contigA.id,
            );

            return No.keepInsertion;
        }

        return Yes.keepInsertion;
    }

    void logStatistics()
    {
        mixin(traceExecution);

        // Due to several transformations of the assembly graph the predicates
        // for "join types" must be combined to yield the expected result.
        alias _isSpannedGap = j => j.isGap && !(j.isOutputGap || j.isDefault);
        alias _isExtension = j => !j.isGap && !(j.isOutputGap || j.isDefault);
        alias _isRemaingGap = isOutputGap;
        alias _isExistingContig = isDefault;

        logJsonInfo(
            "numGraphEdges", assemblyGraph.edges.length,
            "numScaffolds", scaffoldStartNodes.length,
            "numSpannedGaps", assemblyGraph.edges.count!_isSpannedGap,
            "numExtensions", assemblyGraph.edges.count!_isExtension,
            "numRemainingGaps", assemblyGraph.edges.count!_isRemaingGap,
            "numExistingContigs", assemblyGraph.edges.count!_isExistingContig,
        );

        if (shouldLog(LogLevel.diagnostic))
            logSparseInsertionWalks();

        debug logJsonDebug(
            "insertionWalks", scaffoldStartNodes
                .map!(startNode => linearWalk!InsertionInfo(assemblyGraph, startNode, incidentEdgesCache)
                    .map!(join => [
                        "start": join.start.toJson,
                        "end": join.end.toJson,
                        "payload": [
                            "sequence": join.payload.sequence.to!string.toJson,
                            "contigLength": join.payload.contigLength.toJson,
                            "overlaps": join.payload.overlaps.toJson,
                        ].toJson,
                    ])
                    .array
                )
                .array
                .toJson
        );
    }

    void logSparseInsertionWalks()
    {
        stderr.write(`{"scaffolds":[`);

        foreach (i, startNode; scaffoldStartNodes)
        {
            if (i == 0)
                stderr.write(`[`);
            else
                stderr.write(`,[`);

            foreach (enumJoin; enumerate(linearWalk!InsertionInfo(assemblyGraph, startNode, incidentEdgesCache)))
            {
                auto j = enumJoin.index;
                auto join = enumJoin.value;

                if (j > 0)
                    stderr.write(`,`);

                if (join.isDefault)
                    stderr.writef!`{"action":"copy","contigId":%d}`(join.start.contigId);
                else if (join.isOutputGap)
                    stderr.writef!`{"action":"gap","length":%d}`(join.payload.contigLength);
                else if (join.isGap)
                    stderr.writef!`{"action":"insert","contigIds":[%d,%d],"length":%d}`(
                        join.start.contigId,
                        join.end.contigId,
                        join.payload.sequence.length
                    );
                else if (join.isExtension)
                    stderr.writef!`{"action":"insert","contigIds":[%d],"contigPart":"%s","length":%d}`(
                        join.start.contigId,
                        join.isFrontExtension ? "front" : "back",
                        join.payload.sequence.length
                    );
                else
                    assert(0, "unexpected join type");
            }

            stderr.write(`]`);
        }

        stderr.writeln(`]}`);
    }

    void writeNewScaffold(ContigNode startNode)
    {
        mixin(traceExecution);

        auto globalComplement = false;
        auto insertionBegin = startNode;

        logJsonDebug(
            "info", "writing scaffold",
            "scaffoldId", startNode.contigId,
        );

        writeHeader(startNode);
        foreach (currentInsertion; linearWalk!InsertionInfo(assemblyGraph, startNode, incidentEdgesCache))
        {
            writeInsertion(
                insertionBegin,
                currentInsertion,
                globalComplement,
            );

            insertionBegin = currentInsertion.target(insertionBegin);
            if (currentInsertion.isAntiParallel)
                globalComplement = !globalComplement;
        }
        "\n".copy(writer);
    }

    Insertion mergeInsertions(Insertion[] insertionsChunk)
    {
        assert(insertionsChunk.length > 0);

        if (insertionsChunk[0].isDefault)
        {
            auto merged = insertionsChunk[0];

            merged.payload.contigLength = insertionsChunk
                .map!"a.payload.contigLength"
                .maxElement;
            assert(insertionsChunk
                .all!(a => a.payload.contigLength == 0 ||
                           a.payload.contigLength == merged.payload.contigLength));
            merged.payload.overlaps = insertionsChunk
                .map!"a.payload.overlaps"
                .joiner
                .array;

            return merged;
        }
        else
        {
            assert(insertionsChunk.length == 1);

            return insertionsChunk[0];
        }
    }

    protected void writeHeader(in ContigNode begin)
    {
        format!">scaffold-%d\n"(begin.contigId).copy(writer);
    }

    protected void writeInsertion(
        in ContigNode begin,
        in Insertion insertion,
        in bool globalComplement,
    )
    {
        assert(insertion.isValidInsertion, "invalid insertion");

        if (insertion.isDefault)
            writeExistingContig(begin, insertion, globalComplement);
        else if (insertion.isOutputGap)
            writeGap(insertion);
        else
            writeNewSequenceInsertion(begin, insertion, globalComplement);
    }

    protected void writeExistingContig(
        in ContigNode begin,
        in Insertion insertion,
        in bool globalComplement,
    )
    {
        auto insertionInfo = getInfoForExistingContig(begin, insertion, globalComplement);
        auto contigSequence = getFastaSequence(options.refDb, insertionInfo.contigId, options.workdir);
        auto croppedContigSequence = contigSequence[insertionInfo.cropping.begin .. insertionInfo.cropping.end];

        logJsonDebug(
            "info", "writing contig insertion",
            "contigId", insertionInfo.contigId,
            "contigLength", insertion.payload.contigLength,
            "spliceStart", insertionInfo.cropping.begin,
            "spliceEnd", insertionInfo.cropping.end,
            "overlaps", insertion.payload.overlaps.toJson,
            "start", insertion.start.toJson,
            "end", insertion.end.toJson,
        );

        if (insertionInfo.complement)
            croppedContigSequence.reverseComplementer.copy(writer);
        else
            croppedContigSequence.copy(writer);
    }

    protected void writeGap(
        in Insertion insertion,
    )
    {
        auto insertionInfo = getInfoForGap(insertion);

        logJsonDebug(
            "info", "writing gap",
            "gapLength", insertionInfo.length,
            "start", insertion.start.toJson,
            "end", insertion.end.toJson,
        );

        enum char unkownBase = 'n';
        unkownBase
            .repeat
            .takeExactly(insertionInfo.length)
            .copy(writer);
    }

    protected void writeNewSequenceInsertion(
        in ContigNode begin,
        in Insertion insertion,
        in bool globalComplement,
    )
    {
        auto insertionInfo = getInfoForNewSequenceInsertion(begin, insertion, globalComplement);

        logJsonDebug(
            "info", "writing new sequence insertion",
            "type", insertion.isGap ? "gap" : "extension",
            "spliceStart", insertionInfo.cropping.begin,
            "spliceEnd", insertionInfo.cropping.end,
            "insertionLength", insertionInfo.length,
            "isAntiParallel", insertion.isAntiParallel,
            "localComplement", insertionInfo.complement,
            "start", insertion.start.toJson,
            "end", insertion.end.toJson,
        );

        alias highlightIfRequested = (char c) => options.noHighlightInsertions ? c : toUpper(c);

        if (insertionInfo.complement)
            insertionInfo
                .sequence
                .bases!(char, Yes.reverse)
                .map!(complement!char)
                .dropExactly(insertionInfo.cropping.begin)
                .takeExactly(insertionInfo.cropping.size)
                .map!highlightIfRequested
                .copy(writer);
        else
            insertionInfo
                .sequence
                .bases!char
                .dropExactly(insertionInfo.cropping.begin)
                .takeExactly(insertionInfo.cropping.size)
                .map!highlightIfRequested
                .copy(writer);

    }
}


/// Remove contig cropping where no new sequence is to be inserted and adjust
/// cropping where cropped regions overlap.
OutputScaffold fixCropping(
    OutputScaffold scaffold,
    trace_point_t tracePointDistance,
)
{
    mixin(traceExecution);

    alias replace = OutputScaffold.ConflictStrategy.replace;
    auto contigJoins = scaffold.edges.filter!isDefault;
    auto incidentEdgesCache = scaffold.allIncidentEdges();

    foreach (contigJoin; contigJoins)
    {
        auto insertionUpdated = transferCroppingFromIncidentJoins(contigJoin, incidentEdgesCache);

        version (assert)
        {
            auto overlaps = contigJoin.payload.overlaps;

            assert(overlaps.length <= 2);
            if (overlaps.length == 2)
            {
                auto cropPos0 = getCroppingPosition!"contigA"(overlaps[0]);
                auto cropPos1 = getCroppingPosition!"contigA"(overlaps[1]);

                assert(
                    cropPos0 == cropPos1 ||
                    (cropPos0 < cropPos1) == (overlaps[0].seed < overlaps[1].seed)
                );
            }
        }

        if (insertionUpdated)
            scaffold.add!replace(contigJoin);
    }


    return scaffold;
}

Flag!"insertionUpdated" transferCroppingFromIncidentJoins(
    ref Insertion contigJoin,
    OutputScaffold.IncidentEdgesCache incidentEdgesCache,
)
{
    SeededAlignment overlapsFromIncidentJoins(ContigNode contigNode)
    {
        alias isSomeInsertion = insertion =>
            !insertion.isOutputGap && (insertion.isGap || insertion.isExtension);

        auto incidentJoins = incidentEdgesCache[contigNode].filter!isSomeInsertion;

        if (incidentJoins.empty)
            return SeededAlignment(AlignmentChain.disabledInstance);

        assert(incidentJoins.save.walkLength(2) == 1, "non-linearized scaffold");
        auto overlapsForContigNode = incidentJoins
            .front
            .payload
            .overlaps
            .filter!(overlap => overlap.contigA.id == contigNode.contigId);
        assert(!overlapsForContigNode.empty, "missing splice site");
        assert(overlapsForContigNode.save.walkLength(2) <= 1, "too many splice sites");

        return overlapsForContigNode.front;
    }

    contigJoin.payload.overlaps = only(contigJoin.start, contigJoin.end)
        .map!overlapsFromIncidentJoins
        .filter!"!a.flags.disabled"
        .array;

    return cast(typeof(return)) (contigJoin.payload.overlaps.length > 0);
}

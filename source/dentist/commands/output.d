/**
    This is the `output` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.output;

import dentist.commandline : DentistCommand, OptionsFor;
import dentist.common : ReferenceInterval, ReferencePoint;
import dentist.common.alignments :
    AlignmentChain,
    AlignmentLocationSeed,
    id_t,
    trace_point_t;
import dentist.common.binio :
    CompressedSequence,
    InsertionDb;
import dentist.common.insertions :
    getInfoForExistingContig,
    getInfoForGap,
    getInfoForNewSequenceInsertion,
    Insertion,
    InsertionInfo,
    isOutputGap,
    isValidInsertion,
    OutputScaffold,
    SpliceSite;
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
    removeExtensions;
import dentist.dazzler :
    ContigSegment,
    GapSegment,
    getNumContigs,
    getFastaSequence,
    getScaffoldStructure,
    ScaffoldSegment;
import dentist.util.fasta : complement, reverseComplementer;
import dentist.util.log;
import dentist.util.math :
    absdiff,
    add,
    bulkAdd,
    ceil,
    floor,
    mean;
import dentist.util.range : wrapLines;
import std.algorithm :
    all,
    canFind,
    copy,
    count,
    filter,
    find,
    fold,
    joiner,
    map,
    maxElement,
    swap;
import std.array : array;
import std.conv : to;
import std.format : format;
import std.range :
    enumerate,
    only,
    repeat,
    takeExactly,
    walkLength;
import std.range.primitives : empty, front, popFront, save;
import std.stdio : File, stderr, stdout;
import std.typecons : Flag, No, tuple, Yes;
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
        scaffoldStructure = getScaffoldStructure(options.refDb, options).array;
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

        if (!options.shouldExtendContigs)
        {
            assemblyGraph = assemblyGraph.removeExtensions!InsertionInfo();
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
                            "spliceSites": join.payload.spliceSites.toJson,
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
            merged.payload.spliceSites = insertionsChunk
                .map!"a.payload.spliceSites"
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
        contigSequence = contigSequence[insertionInfo.spliceStart .. insertionInfo.spliceEnd];

        logJsonDebug(
            "info", "writing contig insertion",
            "contigId", insertionInfo.contigId,
            "contigLength", insertion.payload.contigLength,
            "spliceStart", insertionInfo.spliceStart,
            "spliceEnd", insertionInfo.spliceEnd,
            "spliceSites", insertion.payload.spliceSites.toJson,
            "start", insertion.start.toJson,
            "end", insertion.end.toJson,
        );

        if (insertionInfo.complement)
            contigSequence.reverseComplementer.copy(writer);
        else
            contigSequence.copy(writer);
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
            "insertionLength", insertion.payload.sequence.length,
            "isAntiParallel", insertion.isAntiParallel,
            "localComplement", insertion.payload.spliceSites[0].flags.complement,
            "start", insertion.start.toJson,
            "end", insertion.end.toJson,
        );

        if (insertionInfo.complement)
            insertionInfo.sequence.bases!(char, Yes.reverse).map!(complement!char).copy(writer);
        else
            insertionInfo.sequence.bases!char.copy(writer);

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
        auto insertionUpdated1 = removeDegradingCropping(contigJoin, incidentEdgesCache);
        auto insertionUpdated2 = resolveCroppingConflicts(scaffold, contigJoin, incidentEdgesCache,
                                                         tracePointDistance);

        if (insertionUpdated1 || insertionUpdated2)
            scaffold.add!replace(contigJoin);
    }

    return scaffold;
}

Flag!"insertionUpdated" removeDegradingCropping(
    ref Insertion contigJoin,
    OutputScaffold.IncidentEdgesCache incidentEdgesCache,
)
{
    typeof(return) insertionUpdated;

    foreach (contigNode; [contigJoin.start, contigJoin.end])
    {
        auto shouldInsertNewSequence = incidentEdgesCache[contigNode]
            .canFind!(insertion => !insertion.isOutputGap && (insertion.isGap || insertion.isExtension));

        if (!shouldInsertNewSequence)
        {
            auto contigLength = contigJoin.payload.contigLength;
            auto newSpliceSites = contigJoin
                .payload
                .spliceSites
                .filter!(spliceSite =>
                    (contigNode.contigPart == ContigPart.begin) ^
                    (spliceSite.alignmentSeed == AlignmentLocationSeed.front))
                .array;
            if (newSpliceSites.length < contigJoin.payload.spliceSites.length)
            {
                contigJoin.payload.spliceSites = newSpliceSites;
                insertionUpdated = Yes.insertionUpdated;
            }
        }
    }

    return insertionUpdated;
}

Flag!"insertionUpdated" resolveCroppingConflicts(
    ref OutputScaffold scaffold,
    ref Insertion contigJoin,
    OutputScaffold.IncidentEdgesCache incidentEdgesCache,
    trace_point_t tracePointDistance,
)
{
    auto spliceSites = contigJoin.payload.spliceSites;

    if (spliceSites.length <= 1)
        return No.insertionUpdated;

    auto frontResult = resolveCroppingConflictsAt!(AlignmentLocationSeed.front)(
        scaffold,
        contigJoin,
        incidentEdgesCache,
        tracePointDistance,
    );
    auto backResult = resolveCroppingConflictsAt!(AlignmentLocationSeed.back)(
        scaffold,
        contigJoin,
        incidentEdgesCache,
        tracePointDistance,
    );

    auto newSpliceSites = only(frontResult, backResult)
        .map!"a[1]"
        .filter!"a.croppingRefPosition.contigId > 0"
        .array;
    auto insertionUpdated = only(frontResult, backResult)
        .map!"a[0]"
        .fold!"a | b";

    if (insertionUpdated)
        contigJoin.payload.spliceSites = newSpliceSites;

    return insertionUpdated;
}

auto resolveCroppingConflictsAt(AlignmentLocationSeed seed)(
    ref OutputScaffold scaffold,
    ref Insertion contigJoin,
    OutputScaffold.IncidentEdgesCache incidentEdgesCache,
    trace_point_t tracePointDistance,
)
{
    auto seedSpliceSites = contigJoin
        .payload
        .spliceSites
        .filter!(spliceSite => spliceSite.alignmentSeed == seed);
    auto numSeedSpliceSites = seedSpliceSites.walkLength(2);

    if (numSeedSpliceSites == 0)
        return tuple(No.insertionUpdated, SpliceSite());
    else if (numSeedSpliceSites == 1)
        return tuple(No.insertionUpdated, seedSpliceSites.front);

    auto commonCroppingPostition = seedSpliceSites
        .map!"a.croppingRefPosition.value"
        .mean;

    static if (seed == AlignmentLocationSeed.front)
        commonCroppingPostition = floor(commonCroppingPostition, tracePointDistance);
    else
        commonCroppingPostition = ceil(commonCroppingPostition, tracePointDistance);

    Flag!"insertionUpdated" insertionUpdated;

    foreach (ref spliceSite; seedSpliceSites)
    {
        auto croppingDiff = absdiff(
            commonCroppingPostition,
            spliceSite.croppingRefPosition.value,
        );

        if (croppingDiff > 0)
        {
            alias replace = OutputScaffold.ConflictStrategy.replace;

            auto contigNode = seed == AlignmentLocationSeed.front
                ? contigJoin.start
                : contigJoin.end;
            auto incidentInsertion = incidentEdgesCache[contigNode]
                .find!(insertion => !insertion.isOutputGap && (insertion.isGap || insertion.isExtension))
                .front;
            incidentInsertion.payload.sequence = seed == AlignmentLocationSeed.front
                ? incidentInsertion.payload.sequence[0 .. $ - croppingDiff]
                : incidentInsertion.payload.sequence[croppingDiff .. $];

            scaffold.add!replace(incidentInsertion);
            insertionUpdated |= Yes.insertionUpdated;
        }
    }

    return tuple(
        insertionUpdated,
        SpliceSite(
            ReferencePoint(
                seedSpliceSites.front.croppingRefPosition.contigId,
                commonCroppingPostition,
            ),
            seed,
            seedSpliceSites.front.flags,
        ),
    );
}

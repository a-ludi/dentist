/**
    This is the `output` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.output;

import dentist.commandline : DentistCommand, OptionsFor;
import dentist.common : ReferencePoint;
import dentist.common.binio :
    CompressedSequence,
    InsertionDb;
import dentist.common.alignments : AlignmentChain, id_t;
import dentist.common.insertions :
    Insertion,
    InsertionInfo,
    isOutputGap,
    OutputScaffold,
    SpliceSite;
import dentist.common.scaffold :
    ContigNode,
    ContigPart,
    contigStarts,
    enforceJoinPolicy,
    getUnkownJoin,
    initScaffold,
    isAntiParallel,
    isDefault,
    isExtension,
    isGap,
    isValid,
    linearWalk,
    normalizeUnkownJoins,
    removeExtensions;
import dentist.dazzler :
    ContigSegment,
    GapSegment,
    getNumContigs,
    getFastaSequences,
    getScaffoldStructure,
    ScaffoldSegment;
import dentist.util.fasta : complement, reverseComplementer;
import dentist.util.log;
import dentist.util.range : wrapLines;
import std.algorithm :
    canFind,
    copy,
    filter,
    joiner,
    map,
    maxElement,
    swap;
import std.array : array;
import std.conv : to;
import std.format : format;
import std.range :
    only,
    repeat,
    takeExactly;
import std.range.primitives : empty, front, popFront, save;
import std.stdio : File, stdout;
import std.typecons : Yes;
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

        foreach (startNode; contigStarts!InsertionInfo(assemblyGraph))
            writeNewContig(startNode);

        debug logJsonDebug(
            "insertionWalks", contigStarts!InsertionInfo(assemblyGraph)
                .map!(startNode => linearWalk!InsertionInfo(assemblyGraph, startNode)
                    .map!(join => [
                        "start": join.start.toJson,
                        "end": join.end.toJson,
                        "payload": join.payload.toJson,
                    ])
                    .array
                )
                .array
                .toJson
        );
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
        assemblyGraph.bulkAdd!mergeInsertions(insertions);
        appendUnkownJoins();

        if (!options.shouldExtendContigs)
        {
            assemblyGraph = assemblyGraph.removeExtensions!InsertionInfo();
        }

        assemblyGraph = assemblyGraph
            .enforceJoinPolicy!InsertionInfo(options.joinPolicy)
            .normalizeUnkownJoins!InsertionInfo()
            .fixContigCropping();
    }


    protected void appendUnkownJoins()
    {
        auto unkownJoins = scaffoldStructure[]
            .filter!(part => part.peek!GapSegment !is null)
            .map!(gapPart => gapPart.get!GapSegment)
            .map!(gapPart => getUnkownJoin(
                gapPart.beginGlobalContigId,
                gapPart.endGlobalContigId,
                InsertionInfo(
                    CompressedSequence(),
                    gapPart.end - gapPart.begin,
                    [
                        SpliceSite(ReferencePoint(
                            gapPart.beginGlobalContigId,
                            gapPart.begin,
                        ), AlignmentChain.emptyFlags),
                        SpliceSite(ReferencePoint(
                            gapPart.beginGlobalContigId,
                            gapPart.end,
                        ), AlignmentChain.emptyFlags),
                    ]
                ),
            ));
        assemblyGraph.bulkAddForce(unkownJoins);
    }

    void writeNewContig(ContigNode startNode)
    {
        mixin(traceExecution);

        auto globalComplement = false;
        auto insertionBegin = startNode;

        writeHeader(startNode);
        foreach (currentInsertion; linearWalk!InsertionInfo(assemblyGraph, startNode))
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

    static Insertion mergeInsertions(Insertion[] insertionsChunk)
    {
        assert(insertionsChunk.length > 0);

        if (insertionsChunk[0].isDefault)
        {
            auto merged = insertionsChunk[0];
            merged.payload.contigLength = insertionsChunk
                .map!"a.payload.contigLength"
                .maxElement;
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
        format!"> contig %d\n"(begin.contigId).copy(writer);
    }

    protected void writeInsertion(
        in ContigNode begin,
        in Insertion insertion,
        in bool globalComplement,
    )
    {
        assert(insertion.isValid, "invalid insertion");

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
        auto spliceSites = insertion.payload.spliceSites;
        auto contigId = begin.contigId;
        auto contigLength = insertion.payload.contigLength;
        size_t spliceStart;
        size_t spliceEnd;

        assert(contigLength > 0);

        switch (spliceSites.length)
        {
        case 0:
            spliceStart = 0;
            spliceEnd = contigLength;
            break;
        case 1:
            auto splicePosition = spliceSites[0].croppingRefPosition.value;

            if (splicePosition < contigLength / 2)
            {
                spliceStart = splicePosition;
                spliceEnd = contigLength;
            }
            else
            {
                spliceStart = 0;
                spliceEnd = splicePosition;
            }
            break;
        case 2:
            assert(spliceSites.length == 2);
            assert(spliceSites[0].croppingRefPosition.contigId
                    == spliceSites[1].croppingRefPosition.contigId);

            spliceStart = spliceSites[0].croppingRefPosition.value;
            spliceEnd = spliceSites[1].croppingRefPosition.value;

            if (spliceEnd < spliceStart)
            {
                swap(spliceStart, spliceEnd);
            }
            break;
        default:
            assert(0, "too many spliceSites");
        }

        assert(globalComplement != (begin < insertion.target(begin)));

        auto contigSequence = getFastaSequences(options.refDb, [contigId], options.workdir).front;
        contigSequence = contigSequence[spliceStart .. spliceEnd];

        logJsonDebug(
            "info", "writing contig insertion",
            "contigId", contigId,
            "contigLength", contigLength,
            "spliceStart", spliceStart,
            "spliceEnd", spliceEnd,
            "spliceSites", spliceSites.toJson,
        );

        if (globalComplement)
            contigSequence.reverseComplementer.copy(writer);
        else
            contigSequence.copy(writer);
    }

    protected void writeGap(
        in Insertion insertion,
    )
    {
        assert(insertion.payload.spliceSites.length == 2);

        logJsonDebug(
            "info", "writing gap",
            "gapLength", insertion.payload.contigLength,
        );

        enum char unkownBase = 'n';
        unkownBase
            .repeat
            .takeExactly(insertion.payload.contigLength)
            .copy(writer);
    }

    protected void writeNewSequenceInsertion(
        in ContigNode begin,
        in Insertion insertion,
        in bool globalComplement,
    )
    {
        auto spliceSites = insertion.payload.spliceSites;
        auto effectiveComplement = spliceSites[0].flags.complement ^ globalComplement;

        assert(
            (insertion.isExtension && spliceSites.length == 1) ^
            (insertion.isGap && spliceSites.length == 2)
        );

        auto newSequence = insertion.payload.sequence;

        if (effectiveComplement)
            newSequence.bases!(char, Yes.reverse).map!(complement!char).copy(writer);
        else
            newSequence.bases!char.copy(writer);

    }
}


/// Remove contig cropping where no new sequence is to be inserted.
OutputScaffold fixContigCropping(OutputScaffold scaffold)
{

    alias replace = OutputScaffold.ConflictStrategy.replace;
    auto contigJoins = scaffold.edges.filter!isDefault;

    foreach (contigJoin; contigJoins)
    {
        bool insertionUpdated;

        foreach (contigNode; [contigJoin.start, contigJoin.end])
        {
            auto shouldInsertNewSequence = scaffold
                .incidentEdges(contigNode)
                .canFind!(insertion => !insertion.isOutputGap && (insertion.isGap || insertion.isExtension));

            if (!shouldInsertNewSequence)
            {
                auto contigLength = contigJoin.payload.contigLength;
                auto newSpliceSites = contigJoin
                    .payload
                    .spliceSites
                    .filter!(spliceSite => contigNode.contigPart == ContigPart.begin
                        ? !(spliceSite.croppingRefPosition.value < contigLength / 2)
                        : !(spliceSite.croppingRefPosition.value >= contigLength / 2))
                    .array;
                if (newSpliceSites.length < contigJoin.payload.spliceSites.length)
                {
                    contigJoin.payload.spliceSites = newSpliceSites;
                    insertionUpdated = true;
                }
            }
        }

        if (insertionUpdated)
        {
            scaffold.add!replace(contigJoin);
        }
    }

    return scaffold;
}

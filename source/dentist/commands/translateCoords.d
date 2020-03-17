/**
    This is the `translateCoords` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.translateCoords;

import dentist.commandline : OptionsFor;
import dentist.common :
    DentistException,
    OutputCoordinate;
import dentist.common.alignments : arithmetic_t, coord_t, id_t;
import dentist.common.binio : InsertionDb;
import dentist.common.commands : DentistCommand;
import dentist.common.insertions :
    getInfoForExistingContig,
    getInfoForGap,
    getInfoForNewSequenceInsertion,
    Insertion,
    InsertionInfo,
    isValidInsertion,
    isOutputGap,
    OutputScaffold;
import dentist.common.scaffold :
    ContigNode,
    isAntiParallel,
    isDefault,
    linearWalk,
    scaffoldStarts;
import dentist.util.algorithm : uniqInPlace;
import dentist.util.log;
import std.algorithm :
    find,
    map,
    maxElement,
    min,
    sort,
    uniq;
import std.array : array;
import std.conv : to;
import std.exception : enforce;
import std.format : format;
import std.range :
    chain,
    dropExactly,
    only;
import std.range.primitives;
import std.stdio : writefln, writeln;
import vibe.data.json :
    toJson = serializeToJson,
    toJsonString = serializeToPrettyJson;


/// Options for the `translateCoords` command.
alias Options = OptionsFor!(DentistCommand.translateCoords);

/// Execute the `translateCoords` command with `options`.
void execute(in Options options)
{
    auto translator = new CoordinateTranslator(options);

    translator.run();
}

enum SequenceType
{
    existing,
    insertion,
    unkown,
}

struct ReferenceCoordinate
{
    id_t contigId;
    arithmetic_t contigCoord;
    coord_t contigLength;
    id_t[] contigIds;
    bool isReverseComplement;
    SequenceType sequenceType;
}

class CoordinateTranslator
{
    alias OriginType = OutputCoordinate.OriginType;

    const(Options) options;
    OutputScaffold assemblyGraph;
    OutputScaffold.IncidentEdgesCache incidentEdgesCache;
    ContigNode[] scaffoldStartNodes;

    this(in Options options)
    {
        this.options = options;
    }

    void run()
    {
        mixin(traceExecution);

        init();

        foreach (i, outputCoordinate; options.outputCoordinates)
        {
            try
            {
                auto walker = Walker(
                    outputCoordinate,
                    assemblyGraph,
                    incidentEdgesCache,
                    scaffoldStartNodes,
                );

                walker.walkToCoordinate();

                if (options.useJson)
                    writeJson(walker.refCoord);
                else
                {
                    if (i > 0)
                        writeln();
                    writeTabular(walker.refCoord);
                }
            }
            catch (DentistException e)
            {
                logJsonError(
                    "info", e.msg,
                    "coord", outputCoordinate.toString(),
                );
            }
        }
    }

    protected void init()
    {
        mixin(traceExecution);

        auto insertionDb = InsertionDb.parse(options.assemblyGraphFile);
        auto insertions = insertionDb[];
        auto nodes = chain(
            insertions.map!"a.start".uniq,
            insertions.map!"a.end".uniq,
        )
            .array
            .sort
            .release;
        uniqInPlace(nodes);

        assemblyGraph = OutputScaffold(nodes);
        assemblyGraph.bulkAddForce(insertions);
        incidentEdgesCache = assemblyGraph.allIncidentEdges();
        scaffoldStartNodes = scaffoldStarts!InsertionInfo(assemblyGraph, incidentEdgesCache).array;
    }

    protected static void writeTabular(ReferenceCoordinate refCoord)
    {
        auto stats = [
            format!`%d`(refCoord.contigId),
            format!`%d`(refCoord.contigCoord),
            format!`%d`(refCoord.contigLength),
            format!`%-(%d,%)`(refCoord.contigIds),
            refCoord.isReverseComplement ? "yes" : "no",
            refCoord.sequenceType.to!string,
        ];
        auto columnWidth = stats.map!"a.length".maxElement;

        writefln!"refContigId:        %*s"(columnWidth, stats[0]);
        writefln!"refContigCoord:     %*s"(columnWidth, stats[1]);
        writefln!"refContigLength:    %*s"(columnWidth, stats[2]);
        writefln!"flankingRefContigs: %*s"(columnWidth, stats[3]);
        writefln!"reverseComplement:  %*s"(columnWidth, stats[4]);
        writefln!"sequenceType:       %*s"(columnWidth, stats[5]);
    }

    protected static void writeJson(ReferenceCoordinate refCoord)
    {
        writeln(toJsonString([
            "refContig": [
                "id": refCoord.contigId,
                "coord": refCoord.contigCoord,
                "length": refCoord.contigLength,
            ].toJson,
            "flankingRefContigs": refCoord.contigIds.toJson,
            "isReverseComplement": refCoord.isReverseComplement.toJson,
            "sequenceType": refCoord.sequenceType.to!string.toJson,
        ]));
    }
}

private struct Walker
{
    alias OriginType = OutputCoordinate.OriginType;

    const(OutputCoordinate) outCoord;
    OutputScaffold assemblyGraph;
    OutputScaffold.IncidentEdgesCache incidentEdgesCache;
    ContigNode[] scaffoldStartNodes;

    coord_t numWalkedBasePairs;
    id_t lastContigId;
    arithmetic_t lastContigStart;
    coord_t lastContigLength;
    SequenceType lastWalkedSequenceType;
    ReferenceCoordinate refCoord;

    protected void walkToCoordinate()
    {
        final switch (outCoord.originType)
        {
        case OriginType.scaffold:
            auto scaffoldStart = getStartOfScaffold(outCoord.scaffoldId);

            walkScaffold(scaffoldStart);
            break;
        case OriginType.scaffoldContig:
            assert(0, "unimplemented");
        case OriginType.global:
            assert(0, "unimplemented");
        case OriginType.contig:
            assert(0, "unimplemented");
        }
    }

    protected ContigNode getStartOfScaffold(id_t scaffoldId)
    {
        alias byContigId = (node, needleId) => node.contigId == needleId;

        auto candidates = scaffoldStartNodes.find!byContigId(scaffoldId);

        enforce!DentistException(!candidates.empty, "scaffold does not exist");

        return candidates.front;
    }

    protected void walkScaffold(ContigNode startNode)
    {
        mixin(traceExecution);

        ContigNode insertionBegin = startNode;
        Insertion lastInsertion;
        auto globalComplement = false;

        foreach (currentInsertion; linearWalk!InsertionInfo(assemblyGraph, startNode, incidentEdgesCache))
        {
            lastInsertion = currentInsertion;

            walkInsertion(
                insertionBegin,
                currentInsertion,
                globalComplement,
            );

            if (numWalkedBasePairs >= outCoord.coord)
                break;

            insertionBegin = currentInsertion.target(insertionBegin);
            if (currentInsertion.isAntiParallel)
                globalComplement = !globalComplement;
        }

        enforce!DentistException(
            numWalkedBasePairs == outCoord.coord,
            "requested coordinate is out of bounds",
        );

        refCoord.contigId = lastContigId;
        refCoord.contigLength = lastContigLength;
        refCoord.contigCoord = numWalkedBasePairs - lastContigStart;
        if (globalComplement)
            refCoord.contigCoord = lastContigLength - refCoord.contigCoord;
        refCoord.contigIds = only(
            cast(id_t) lastInsertion.start.contigId,
            cast(id_t) lastInsertion.end.contigId,
        ).uniq.array;
        refCoord.sequenceType = lastWalkedSequenceType;
        refCoord.isReverseComplement = globalComplement;
    }

    protected void walkInsertion(
        in ContigNode begin,
        in Insertion insertion,
        in bool globalComplement,
    )
    {
        assert(insertion.isValidInsertion, "invalid insertion");
        size_t length;

        if (insertion.isDefault)
        {
            auto info = getInfoForExistingContig(begin, insertion, globalComplement);

            lastContigId = info.contigId;
            lastContigStart = numWalkedBasePairs;
            lastContigLength = cast(coord_t) info.contigLength;
            length = info.length;
            lastWalkedSequenceType = SequenceType.existing;
        }
        else if (insertion.isOutputGap)
        {
            length = getInfoForGap(insertion).length;
            lastWalkedSequenceType = SequenceType.unkown;
        }
        else
        {
            length = getInfoForNewSequenceInsertion(begin, insertion, globalComplement).length;
            lastWalkedSequenceType = SequenceType.insertion;
        }

        numWalkedBasePairs = min(
            numWalkedBasePairs + cast(coord_t) length,
            outCoord.coord,
        );
    }
}

/**
    This is the `translateCoords` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.translateCoords;

import dentist.common : isTesting;

static if (isTesting):

import dentist.commandline : OptionsFor, TestingCommand;
import dentist.common :
    DentistException,
    OutputCoordinate;
import dentist.common.alignments : arithmetic_t, coord_t, id_t;
import dentist.common.binio : InsertionDb;
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
    joiner,
    map,
    min,
    sort,
    uniq;
import std.array : array;
import std.exception : enforce;
import std.range :
    dropExactly,
    only;
import std.range.primitives;
import std.stdio : writeln;
import vibe.data.json :
    serializeToJsonString,
    toJson = serializeToJson,
    toJsonString = serializeToPrettyJson;


/// Options for the `translateCoords` command.
alias Options = OptionsFor!(TestingCommand.translateCoords);

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

    this(in Options options)
    {
        this.options = options;
        this.outCoord = options.outputCoordinate;
    }

    void run()
    {
        mixin(traceExecution);

        init();
        walkToCoordinate();

        writeln(refCoord);
    }

    protected void init()
    {
        mixin(traceExecution);

        auto insertionDb = InsertionDb.parse(options.assemblyGraphFile);
        auto insertions = insertionDb[];
        auto nodes = insertions
            .map!(ins => only(ins.start, ins.end))
            .joiner
            .array
            .sort
            .release;
        uniqInPlace(nodes);

        assemblyGraph = OutputScaffold(nodes);
        assemblyGraph.bulkAddForce(insertions);
        incidentEdgesCache = assemblyGraph.allIncidentEdges();
        scaffoldStartNodes = scaffoldStarts!InsertionInfo(assemblyGraph, incidentEdgesCache).array;
    }

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

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
    isCyclic,
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


enum AGPComponentType : string
{
    /// Active Finishing
    activeFinishing = "A",
    /// Draft HTG (often phase1 and phase2 are called Draft, whether or not they have the draft keyword).
    draftHTG = "D",
    /// Finished HTG (phase3)
    finishedHTG = "F",
    /// Whole Genome Finishing
    wholeGenomeFinishing = "G",
    /// Other sequence (typically means no HTG keyword)
    otherSequence = "O",
    /// Pre Draft
    preDraft = "P",
    /// WGS contig
    wgsContig = "W",
    /// gap with specified size
    gapWithSpecifiedSize = "N",
    /// gap of unknown size, defaulting to 100 bases.
    gapOfUnknownSize = "U",
}


enum AGPLinkageEvidence : string
{
    /// used when no linkage is being asserted (column 8b is ‘no’)
    na = "na",
    /// paired sequences from the two ends of a DNA fragment, mate-pairs and molecular-barcoding.
    pairedEnds = "paired-ends",
    /// alignment to a reference genome within the same genus.
    alignGenus = "align_genus",
    /// alignment to a reference genome within another genus.
    alignXgenus = "align_xgenus",
    /// alignment to a transcript from the same species.
    alignTrnscpt = "align_trnscpt",
    /// sequence on both sides of the gap is derived from the same clone, but the gap is not spanned by paired-ends. The adjacent sequence contigs have unknown order and orientation.
    withinClone = "within_clone",
    /// linkage is provided by a clone contig in the tiling path (TPF). For example, a gap where there is a known clone, but there is not yet sequence for that clone.
    cloneContig = "clone_contig",
    /// linkage asserted using a non-sequence based map such as RH, linkage, fingerprint or optical.
    map = "map",
    /// PCR using primers on both sides of the gap.
    pcr = "pcr",
    /// ligation of segments of DNA that were brought into proximity in chromatin (Hi-C and related technologies).
    proximityLigation = "proximity_ligation",
    /// strobe sequencing.
    strobe = "strobe",
    /// used only for gaps of type contamination and when converting old AGPs that lack a field for linkage evidence into the new format.
    unspecified = "unspecified",
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
    File resultFile;
    FastaWriter writer;
    File agpFile;
    string currentScaffold;
    id_t currentScaffoldPartId;
    coord_t currentScaffoldCoord;

    this(in ref Options options)
    {
        this.options = options;
        this.resultFile = options.resultFile is null
            ? stdout
            : File(options.resultFile, "w");
        this.writer = wrapLines(resultFile.lockingTextWriter, options.fastaLineWidth);
        if (options.agpFile !is null)
            this.agpFile = File(options.agpFile, "w");
    }

    void run()
    {
        mixin(traceExecution);

        init();
        buildAssemblyGraph();
        scaffoldStartNodes = scaffoldStarts!InsertionInfo(assemblyGraph, incidentEdgesCache).array;
        logStatistics();
        if (agpFile.isOpen)
            writeAGPHeader();

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

    void writeAGPHeader()
    {
        import dentist.swinfo;

        assert(agpFile.isOpen);
        agpFile.writefln!"##agp-version\t%s"(options.agpVersion);
        agpFile.writefln!"# TOOL: %s %s"(executableName, version_);
        agpFile.writefln!"# INPUT_ASSEMBLY: %s"(options.refDb);
        agpFile.writefln!"# object\tobject_beg\tobject_end\tpart_number\tcomponent_type\tcomponent_id/gap_length\tcomponent_beg/gap_type\tcomponent_end/linkage\torientation\tlinkage_evidence"();
    }

    void logSparseInsertionWalks()
    {
        auto oldLogLevel = getLogLevel();
        setLogLevel(LogLevel.fatal);
        scope (exit)
            setLogLevel(oldLogLevel);

        stderr.write(`{"scaffolds":[`);

        foreach (i, startNode; scaffoldStartNodes)
        {
            if (i == 0)
                stderr.write(`[`);
            else
                stderr.write(`,[`);

            auto insertionBegin = startNode;
            coord_t currentPosition;

            foreach (enumJoin; enumerate(linearWalk!InsertionInfo(assemblyGraph, startNode, incidentEdgesCache)))
            {
                auto j = enumJoin.index;
                auto join = enumJoin.value;

                if (j > 0)
                    stderr.write(`,`);

                if (join.isDefault)
                {
                    auto insertionInfo = getInfoForExistingContig(insertionBegin, join, false);

                    stderr.writef!`{"action":"copy","contigId":%d,"length":%d,"position":%d}`(
                        join.start.contigId,
                        insertionInfo.length,
                        currentPosition,
                    );
                    currentPosition += insertionInfo.length;
                }
                else if (join.isOutputGap)
                {
                    auto insertionInfo = getInfoForGap(join);

                    stderr.writef!`{"action":"gap","length":%d,"position":%d}`(
                        insertionInfo.length,
                        currentPosition,
                    );
                    currentPosition += insertionInfo.length;
                }
                else if (join.isGap)
                {
                    auto insertionInfo = getInfoForNewSequenceInsertion(insertionBegin, join, false);

                    stderr.writef!`{"action":"insert","contigs":[{"id":%d,"pos":"%s"},{"id":%d,"pos":"%s"}],"length":%d,"position":%d}`(
                        join.start.contigId,
                        join.start.contigPart.to!string,
                        join.end.contigId,
                        join.end.contigPart.to!string,
                        insertionInfo.length,
                        currentPosition,
                    );
                    currentPosition += insertionInfo.length;
                }
                else if (join.isExtension)
                {
                    auto insertionInfo = getInfoForNewSequenceInsertion(insertionBegin, join, false);

                    stderr.writef!`{"action":"insert","contigs":[{"id":%d,"pos":"%s"}],"length":%d,"position":%d}`(
                        join.start.contigId,
                        join.isFrontExtension ? "front" : "back",
                        join.payload.sequence.length,
                        currentPosition,
                    );
                    currentPosition += insertionInfo.length;
                }
                else
                {
                    assert(0, "unexpected join type");
                }

                insertionBegin = join.target(insertionBegin);
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

        currentScaffold = scaffoldHeader(startNode, isCyclic!InsertionInfo(assemblyGraph, startNode, incidentEdgesCache));
        currentScaffoldPartId = 1;
        currentScaffoldCoord = 1;
        writeHeader();
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
            ++currentScaffoldPartId;
            currentScaffoldCoord += currentInsertion.payload.contigLength - 1;
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

    protected string scaffoldHeader(in ContigNode begin, in Flag!"isCyclic" isCyclic)
    {
        if (isCyclic)
            return format!"circular-scaffold-%d"(begin.contigId);
        else
            return format!"scaffold-%d"(begin.contigId);
    }

    protected void writeHeader()
    {
        format!">%s\n"(currentScaffold).copy(writer);
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

        if (agpFile.isOpen)
            agpFile.writeln(only(
                currentScaffold,
                to!string(currentScaffoldCoord),
                to!string(currentScaffoldCoord + insertion.payload.contigLength - 1),
                to!string(currentScaffoldPartId),
                cast(string) AGPComponentType.wgsContig,
                to!string(insertionInfo.contigId),
                to!string(insertionInfo.cropping.begin),
                to!string(insertionInfo.cropping.end),
                insertionInfo.complement ? "+" : "-",
                cast(string) AGPLinkageEvidence.na,
            ).joiner("\t"));

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

        if (agpFile.isOpen)
            agpFile.writeln(only(
                currentScaffold,
                to!string(currentScaffoldCoord),
                to!string(currentScaffoldCoord + insertion.payload.contigLength - 1),
                to!string(currentScaffoldPartId),
                cast(string) AGPComponentType.gapWithSpecifiedSize,
                to!string(insertion.payload.contigLength),
                to!string("scaffold"),
                to!string("yes"),
                to!string("na"),
                cast(string) AGPLinkageEvidence.unspecified,
            ).joiner("\t"));

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

        if (agpFile.isOpen)
            agpFile.writeln(only(
                currentScaffold,
                to!string(currentScaffoldCoord),
                to!string(currentScaffoldCoord + insertion.payload.contigLength - 1),
                to!string(currentScaffoldPartId),
                cast(string) AGPComponentType.otherSequence,
                format!"reads-%(%d-%)"(insertion.payload.readIds),
                to!string(insertionInfo.cropping.begin),
                to!string(insertionInfo.cropping.end),
                to!string(insertionInfo.complement ? '+' : '-'),
                cast(string) AGPLinkageEvidence.cloneContig,
            ).joiner("\t"));

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

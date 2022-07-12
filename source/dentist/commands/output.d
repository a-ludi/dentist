/**
    This is the `output` command of DENTIST.

    Command_Summary:

    ---
    Generate the output assembly by closing gaps.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.output;

package(dentist) enum summary = "
    Generate the output assembly by closing gaps.
";

import dentist.commandline : OptionsFor;
import dentist.common :
    isTesting,
    ReferenceInterval,
    ReferenceRegion,
    ReferencePoint;
import dentist.common.alignments :
    AlignmentChain,
    AlignmentLocationSeed,
    coord_t,
    id_t,
    ReadAlignmentType,
    SeededAlignment,
    trace_point_t;
static if (isTesting)
    import dentist.commands.checkResults :
        Complement,
        ContigAlignmentsCache,
        ContigMapping,
        DuplicateQueryContig;
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
    isBackExtension,
    isGap,
    linearWalk,
    normalizeUnkownJoins,
    removeBlacklisted;
import dentist.dazzler :
    ContigSegment,
    GapSegment,
    getDbHeaders,
    getNumContigs,
    getFastaSequence,
    getScaffoldStructure,
    ScaffoldSegment;
import dentist.util.algorithm : replaceInPlace, sliceUntil;
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
    NaturalNumberSet,
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
    sort,
    swap,
    swapAt;
import std.array : array, join, split;
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
import std.ascii : isWhite;
import vibe.data.json : Json, toJson = serializeToJson;


/// Options for the `output` command.
alias Options = OptionsFor!(DentistCommand.output);


/// Execute the `output` command with `options`.
void execute(in Options options)
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
    const(ContigSegment)[] contigs;
    string[] readFastaIds;
    bool[size_t[2]] skipGaps;
    StringUniqifier!size_t getUniqScaffold;
    OutputScaffold assemblyGraph;
    OutputScaffold.IncidentEdgesCache incidentEdgesCache;
    ContigNode[] scaffoldStartNodes;
    File resultFile;
    FastaWriter writer;
    File agpFile;
    File closedGapsBedFile;
    string currentScaffold;
    id_t currentScaffoldPartId;
    coord_t currentScaffoldCoord;
    coord_t nextScaffoldCoord;
    id_t currentContigId;
    coord_t currentContigCoord;
    coord_t nextContigCoord;
    static if (isTesting)
    {
        ContigAlignmentsCache contigAlignmentsCache;
        ContigMapping[] contigAlignments;
    }

    this(in ref Options options)
    {
        this.options = options;
        this.resultFile = options.resultFile is null
            ? stdout
            : File(options.resultFile, "w");
        this.writer = wrapLines(resultFile.lockingTextWriter, options.fastaLineWidth);
        if (options.agpFile !is null)
            this.agpFile = File(options.agpFile, "w");
        if (options.closedGapsBedFile !is null)
            this.closedGapsBedFile = File(options.closedGapsBedFile, "w");
        static if (isTesting)
            if (options.contigAlignmentsCache !is null)
            {
                with (this.contigAlignmentsCache)
                {
                    contigAlignmentsCache = options.contigAlignmentsCache;
                    dbA = options.resultFile;
                    dbB = options.refDb;
                }
            }
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

        gotoNextContig();

        static if (isTesting)
            if (contigAlignmentsCache.contigAlignmentsCache !is null)
                contigAlignmentsCache.write(
                    contigAlignments,
                    NaturalNumberSet(),
                    Yes.forceOverwrite,
                );
    }

    protected void init()
    {
        mixin(traceExecution);

        numReferenceContigs = getNumContigs(options.refDb);
        scaffoldStructure = getScaffoldStructure(options.refDb).array;
        contigs = scaffoldStructure
            .filter!(part => part.peek!ContigSegment !is null)
            .map!(contigPart => contigPart.get!ContigSegment)
            .array;
        readFastaIds = !agpFile.isOpen || options.agpDazzler || options.agpSkipReadIds
            ? []
            : getDbHeaders(options.readsDb);

        // remove leading `>` and take only the prefix up to a whitespace
        foreach (ref header; readFastaIds)
            header = header[1 .. $].sliceUntil!isWhite;

        foreach (skipGap; options.skipGaps)
        {

            skipGaps[[
                cast(size_t) skipGap[0],
                cast(size_t) skipGap[1],
            ]] = true;
        }
        static if (isTesting)
            contigAlignments.reserve(numReferenceContigs);
    }


    protected void buildAssemblyGraph()
    {
        mixin(traceExecution);

        auto insertionDb = InsertionDb.parse(options.insertionsFile);
        auto allInsertions = insertionDb[].sort;
        auto insertions = allInsertions
            .enumerate
            .filter!(enumInsertion => !enumInsertion.value.isExtension || (
                options.onlyFlags.extending &&
                skipShortExtension(enumInsertion.expand)
            ))
            .filter!(enumInsertion => !enumInsertion.value.isGap || options.onlyFlags.spanning)
            .filter!(enumInsertion => ensureHighQualityConsensus(enumInsertion.expand))
            .map!"a.value";

        assemblyGraph = initScaffold!(
            contigId => InsertionInfo(CompressedSequence(), contigs[contigId - 1].length, []),
            InsertionInfo,
        )(numReferenceContigs);
        assemblyGraph.bulkAdd!(joins => mergeInsertions(joins))(insertions);
        appendUnkownJoins();

        Insertion[] blacklistedInsertions;
        Insertion[] forbiddenInsertions;
        assemblyGraph = assemblyGraph
            .removeBlacklisted!InsertionInfo(skipGaps, blacklistedInsertions)
            .enforceJoinPolicy!InsertionInfo(options.joinPolicy, forbiddenInsertions)
            .normalizeUnkownJoins!InsertionInfo()
            .fixCropping();
        incidentEdgesCache = assemblyGraph.allIncidentEdges();

        foreach (forbiddenInsertion; forbiddenInsertions)
            logJsonInfo(
                "info", "skipping pile up due to `joinPolicy`",
                "event", "insertionSkipped",
                "reason", "joinPolicy",
                "insertionId", allInsertions.lowerBound(forbiddenInsertion).length,
                "insertion", forbiddenInsertion.insertionToSimpleJson,
            );

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

    protected Flag!"keepInsertion" skipShortExtension(size_t insertionId, Insertion insertion) const
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
                "event", "insertionSkipped",
                "reason", "minExtensionLength",
                "insertionId", insertionId,
                "insertion", insertion.insertionToSimpleJson,
            );

            return No.keepInsertion;
        }

        return Yes.keepInsertion;
    }

    protected Flag!"keepInsertion" ensureHighQualityConsensus(size_t insertionId, Insertion insertion) const
    {
        foreach (consensusAlignment; insertion.payload.overlaps)
        {
            if (consensusAlignment.averageErrorRate > options.maxInsertionError)
            {
                logJsonWarn(
                    "info", "skipping insertion due to `maxInsertionError`",
                    "event", "insertionSkipped",
                    "reason", "maxInsertionError",
                    "insertionId", insertionId,
                    "insertion", insertion.insertionToSimpleJson,
                    "skippedContigId", consensusAlignment.contigA.id,
                    "averageErrorRate", consensusAlignment.averageErrorRate,
                    "maxErrorRate", options.maxInsertionError,
                );

                return No.keepInsertion;
            }
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

    void writeAGPContig(
        id_t contigId,
        size_t contigBegin,
        size_t contigEnd,
        bool complement,
        AGPLinkageEvidence linkageEvidence = AGPLinkageEvidence.na,
    )
    in (agpFile.isOpen)
    {
        const contig = options.agpDazzler
            ? ContigSegment(contigId)
            : contigs[contigId - 1];
        assert(contig.globalContigId == contigId);

        writeAGPComponent(
            AGPComponentType.wgsContig,
            options.agpDazzler
                ? contigId.to!string
                : contig.header[1 .. $].split("\t")[0],
            contig.begin + contigBegin,
            contig.begin + contigEnd,
            complement,
            linkageEvidence,
        );
    }

    void writeAGPInsertion(
        const id_t[] readIds,
        size_t insertionBegin,
        size_t insertionEnd,
        bool complement,
        AGPLinkageEvidence linkageEvidence = AGPLinkageEvidence.cloneContig,
    )
    in (agpFile.isOpen)
    {
        writeAGPComponent(
            AGPComponentType.otherSequence,
            options.agpSkipReadIds
                ? format!"%d reads"(readIds.length)
                : (options.agpDazzler
                    ? format!"reads-%(%d-%)"(readIds[])
                    : readIds[]
                        .map!(readId => readFastaIds[readId - 1])
                        .join(" ")
                ),
            insertionBegin,
            insertionEnd,
            complement,
            linkageEvidence,
        );
    }

    void writeAGPComponent(
        AGPComponentType componentType,
        string componentId,
        size_t componentBeg,
        size_t componentEnd,
        bool complement,
        AGPLinkageEvidence linkageEvidence,
    )
    in (agpFile.isOpen)
    {
        writeAGPObject();
        // component_type
        agpFile.writef!"%s\t"(cast(string) componentType);
        // component_id/gap_length
        agpFile.writef!"%s\t"(componentId);
        // component_beg/gap_type
        agpFile.writef!"%d\t"(componentBeg);
        // component_end/linkage
        agpFile.writef!"%d\t"(componentEnd);
        // orientation
        agpFile.writef!"%s\t"(complement ? "+" : "-");
        // linkage_evidence
        agpFile.writef!"%s\n"(cast(string) linkageEvidence);
    }

    void writeAGPGap(
        size_t gapLength,
        string gapType = "scaffold",
        AGPLinkageEvidence linkageEvidence = AGPLinkageEvidence.unspecified,
    )
    in (agpFile.isOpen)
    {
        writeAGPObject();
        // component_type
        agpFile.writef!"%s\t"(cast(string) AGPComponentType.gapWithSpecifiedSize);
        // component_id/gap_length
        agpFile.writef!"%d\t"(gapLength);
        // component_beg/gap_type
        agpFile.writef!"%s\t"(gapType);
        // component_end/linkage
        agpFile.writef!"%s\t"("yes");
        // orientation
        agpFile.writef!"%s\t"("na");
        // linkage_evidence
        agpFile.writef!"%s\n"(cast(string) linkageEvidence);
    }

    void writeAGPObject()
    in (agpFile.isOpen)
    {
        string object = currentScaffold.split("\t")[0];
        coord_t objectBeg = currentScaffoldCoord;
        coord_t objectEnd = nextScaffoldCoord - 1;
        id_t partNumber = currentScaffoldPartId;

        agpFile.writef!"%s\t%d\t%d\t%d\t"(object, objectBeg, objectEnd, partNumber);
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
        gotoNextContig();
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
            currentScaffoldCoord = nextScaffoldCoord;
            currentContigCoord = nextContigCoord;
        }
        "\n".copy(writer);
    }

    void gotoNextContig()
    {
        static if (isTesting)
            foreach_reverse (ref contigAlignment; contigAlignments)
            {
                if (contigAlignment.reference.contigId == currentContigId)
                    contigAlignment.referenceContigLength = currentContigCoord;
                else
                    break;
            }

        ++currentContigId;
        currentContigCoord = 0;
        nextContigCoord = 0;
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
        const contigId = begin.contigId;
        // get original scaffold ID and make sure there is no tab contained
        const origScaffId = contigs[contigId - 1].header[1 .. $].sliceUntil('\t');
        // get a uniquified version of `origScaffId`
        const uniqScaffId = getUniqScaffold(contigId, origScaffId);
        // append additional information
        const circularSuffix = isCyclic ? "\tisCyclic" : "";
        const header = format!"%s\tscaffold-%d%s"(
            uniqScaffId,
            contigId,
            circularSuffix,
        );

        return header;
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
        auto contigSequence = getFastaSequence(options.refDb, insertionInfo.contigId);
        auto croppedContigSequence = contigSequence[insertionInfo.cropping.begin .. insertionInfo.cropping.end];
        nextScaffoldCoord = currentScaffoldCoord + cast(coord_t) insertionInfo.length;
        nextContigCoord = currentContigCoord + cast(coord_t) insertionInfo.length;

        if (agpFile.isOpen)
            writeAGPContig(
                insertionInfo.contigId,
                insertionInfo.cropping.begin,
                insertionInfo.cropping.end,
                insertionInfo.complement,
            );

        static if (isTesting)
            contigAlignments ~= ContigMapping(
                ReferenceInterval(
                    currentContigId,
                    currentContigCoord,
                    nextContigCoord,
                ),
                0,
                insertionInfo.contigId,
                DuplicateQueryContig.no,
                cast(Complement) insertionInfo.complement,
            );

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
        nextScaffoldCoord = currentScaffoldCoord + cast(coord_t) insertionInfo.length;
        nextContigCoord = currentContigCoord + cast(coord_t) insertionInfo.length;

        if (agpFile.isOpen)
            writeAGPGap(insertion.payload.contigLength);

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

        gotoNextContig();
    }

    protected void writeNewSequenceInsertion(
        in ContigNode begin,
        in Insertion insertion,
        in bool globalComplement,
    )
    {
        auto insertionInfo = getInfoForNewSequenceInsertion(begin, insertion, globalComplement);
        nextScaffoldCoord = currentScaffoldCoord + cast(coord_t) insertionInfo.length;
        nextContigCoord = currentContigCoord + cast(coord_t) insertionInfo.length;
        auto leftContigId = begin.contigId;
        auto rightContigId = insertion.target(begin).contigId;

        if (agpFile.isOpen)
            writeAGPInsertion(
                insertion.payload.readIds,
                insertionInfo.cropping.begin,
                insertionInfo.cropping.end,
                insertionInfo.complement,
            );

        if (closedGapsBedFile.isOpen)
            closedGapsBedFile.writeln(only(
                currentScaffold.split("\t")[0],
                to!string(currentScaffoldCoord - 1),
                to!string(nextScaffoldCoord),
                format!("%s-%d-%d|%s-%(%d-%)")(
                    Options.contigsExtraName,
                    leftContigId,
                    rightContigId,
                    Options.readsExtraName,
                    insertion.payload.readIds,
                ),
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
                .dropExactly(insertionInfo.sequence.length - insertionInfo.cropping.end)
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

/// Converts the pileup into a simple JSON object for diagnostic purposes.
Json insertionToSimpleJson(in Insertion insertion)
{
    return [
        "type": insertion.getType.to!string.toJson,
        "numReads": insertion.payload.readIds.length.toJson,
        "contigIds": insertion
            .payload
            .overlaps
            .map!(sa => sa.contigA.id + 0)
            .array
            .sort
            .release
            .toJson,
    ].toJson;
}

private ReadAlignmentType getType(in Insertion insertion)
{
    if (insertion.isGap)
        return ReadAlignmentType.gap;
    else if (insertion.isFrontExtension)
        return ReadAlignmentType.front;
    else if (insertion.isBackExtension)
        return ReadAlignmentType.back;
    else
        assert(0, "illegal insertion type");
}


struct StringUniqifier(K, string fmt = "%s-%d")
{
    protected size_t[string] dupCounts;
    protected string[K] cache;

    string opCall(K key, string label)
    {
        // check cache first
        auto cacheHit = key in cache;
        if (cacheHit)
            // found ulabel in cache
            return *cacheHit;

        // get the number of duplciates of `label`
        auto dupCount = dupCounts.get(label, 0);
        auto ulabel = dupCount == 0
            // use plain `label` if there are no duplicates
            ? label
            // create a unique label from `fmt`
            : format!fmt(label, dupCount);

        while (ulabel in dupCounts)
            ulabel = format!fmt(label, ++dupCount);

        // add `ulabel` to cache
        cache[key] = ulabel;
        // update dupcliate counts
        dupCounts[label] = dupCount + 1;

        return ulabel;
    }
}

unittest
{
    StringUniqifier!int uniq;

    assert(uniq(1, "A") == "A");
    assert(uniq(1, "A") == "A");
    assert(uniq(2, "A") == "A-1");
    assert(uniq(2, "A") == "A-1");
    assert(uniq(3, "A") == "A-2");
    assert(uniq(4, "B") == "B");
}

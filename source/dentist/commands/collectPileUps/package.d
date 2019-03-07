/**
    This is the `collectPileUps` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.collectPileUps;

import dentist.commandline : OptionsFor;
import dentist.commands.collectPileUps.filter :
    AmbiguousAlignmentChainsFilter,
    ImproperAlignmentChainsFilter,
    RedundantAlignmentChainsFilter,
    WeaklyAnchoredAlignmentChainsFilter;
import dentist.common :
    DentistException,
    ReferenceInterval,
    ReferenceRegion;
import dentist.common.alignments :
    AlignmentChain,
    getType,
    PileUp;
import dentist.common.commands : DentistCommand;
import dentist.common.binio : writePileUpsDb;
import dentist.dazzler :
    GapSegment,
    getAlignments,
    getNumContigs,
    getScaffoldStructure,
    readMask;
import dentist.util.log;
import dentist.util.math : NaturalNumberSet;
import std.algorithm :
    canFind,
    count,
    filter,
    map,
    sum;
import std.array : array;
import std.conv : to;
import std.exception : enforce;
import std.typecons : tuple, Yes;
import vibe.data.json : toJson = serializeToJson;

/// Options for the `collectPileUps` command.
alias Options = OptionsFor!(DentistCommand.collectPileUps);

/// Execute the `collectPileUps` command with `options`.
void execute(in Options options)
{
    auto collector = new PileUpCollector(options);

    collector.run();
}

/// This class comprises the `collectPileUps` step of the `dentist` algorithm
class PileUpCollector
{
    protected const Options options;
    protected size_t numReferenceContigs;
    protected size_t numReads;
    protected AlignmentChain[] readsAlignment;
    protected ReferenceRegion repetitiveRegions;
    protected GapSegment[] inputGaps;
    protected NaturalNumberSet unusedReads;

    this(in ref Options options)
    {
        this.options = options;
    }

    void run()
    {
        mixin(traceExecution);

        readInputs();
        filterAlignments();
        auto pileUps = buildPileUps();
        writePileUps(pileUps);
    }

    protected void readInputs()
    {
        mixin(traceExecution);

        numReferenceContigs = getNumContigs(options.refDb, options.workdir);
        numReads = getNumContigs(options.readsDb, options.workdir);
        unusedReads = NaturalNumberSet(numReads, Yes.addAll);
        inputGaps = getScaffoldStructure(options.refDb)
            .filter!(part => part.peek!GapSegment !is null)
            .map!(gapPart => gapPart.get!GapSegment)
            .array;

        logJsonInfo(
            "numReferenceContigs", numReferenceContigs,
            "numReads", numReads,
        );
        readsAlignment = getAlignments(
            options.refDb,
            options.readsDb,
            options.readsAlignmentFile,
            options.workdir,
            options.tracePointDistance,
        );

        foreach (mask; options.repeatMasks)
            repetitiveRegions |= ReferenceRegion(readMask!ReferenceInterval(
                options.refDb,
                mask,
                options.workdir,
            ));

        enforce!DentistException(readsAlignment.length > 0, "empty ref vs. reads alignment");
    }

    protected void filterAlignments()
    {
        mixin(traceExecution);

        auto filters = tuple(
            new ImproperAlignmentChainsFilter(),
            new WeaklyAnchoredAlignmentChainsFilter(repetitiveRegions, options.minAnchorLength),
            new AmbiguousAlignmentChainsFilter(&unusedReads),
            new RedundantAlignmentChainsFilter(&unusedReads),
        );
        logJsonDiagnostic(
            "filterStage", "Input",
            "readsAlignment", shouldLog(LogLevel.debug_)
                ? readsAlignment.toJson
                : toJson(null),
            "numAlignmentChains", readsAlignment.count!"!a.flags.disabled",
        );

        foreach (i, filter; filters)
        {
            readsAlignment = filter(readsAlignment);

            logJsonDiagnostic(
                "filterStage", typeof(filter).stringof,
                "readsAlignment", shouldLog(LogLevel.debug_)
                    ? readsAlignment.toJson
                    : toJson(null),
                "numAlignmentChains", readsAlignment.count!"!a.flags.disabled",
            );
        }

        enforce!DentistException(
            readsAlignment.canFind!"!a.flags.disabled",
            "no alignment chains left after filtering",
        );
    }

    protected PileUp[] buildPileUps()
    {
        mixin(traceExecution);

        import dentist.commands.collectPileUps.pileups : build;
        auto pileUps = build(
            numReferenceContigs,
            readsAlignment,
            inputGaps,
            options,
        );

        logJsonDebug("pileUps", pileUps
            .map!(pileUp => toJson([
                "type": toJson(pileUp.getType.to!string),
                "readAlignments": pileUp.map!"a[]".array.toJson,
            ]))
            .array
            .toJson);
        logJsonInfo(
            "numPileUps", pileUps.length,
            "numAlignmentChains", pileUps.map!"a[].length".sum,
        );

        return pileUps;
    }

    protected void writePileUps(PileUp[] pileUps)
    {
        mixin(traceExecution);

        writePileUpsDb(pileUps, options.pileUpsFile);
    }
}

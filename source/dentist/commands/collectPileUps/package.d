/**
    This is the `collectPileUps` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.collectPileUps;

import dentist.commandline : DentistCommand, OptionsFor;
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
import dentist.common.binio : writePileUpsDb;
import dentist.dazzler :
    getAlignments,
    getNumContigs,
    readMask;
import dentist.util.log;
import dentist.util.math : NaturalNumberSet;
import std.algorithm : count, isSorted, map, sum;
import std.array : array;
import std.conv : to;
import std.exception : enforce;
import std.typecons : tuple;
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
    protected NaturalNumberSet unusedReads;

    this(in ref Options options)
    {
        this.options = options;
    }

    void run()
    {
        mixin(traceExecution);

        readInputs();
        initUnusedReads();
        filterAlignments();
        auto pileUps = buildPileUps();
        writePileUps(pileUps);
    }

    protected void readInputs()
    {
        mixin(traceExecution);

        numReferenceContigs = getNumContigs(options.refDb, options.workdir);
        numReads = getNumContigs(options.readsDb, options.workdir);
        logJsonInfo(
            "numReferenceContigs", numReferenceContigs,
            "numReads", numReads,
        );
        readsAlignment = getAlignments(
            options.refDb,
            options.readsDb,
            options.readsAlignmentFile,
            options.workdir,
        );
        repetitiveRegions = ReferenceRegion(readMask!ReferenceInterval(
            options.refDb,
            options.repeatMask,
            options.workdir,
        ));

        foreach (mask; options.additionalMasks)
            repetitiveRegions |= ReferenceRegion(readMask!ReferenceInterval(
                options.refDb,
                mask,
                options.workdir,
            ));

        enforce!DentistException(readsAlignment.length > 0, "empty ref vs. reads alignment");
    }

    protected void initUnusedReads()
    {
        mixin(traceExecution);

        unusedReads.reserveFor(numReads);
        foreach (readId; 1 .. numReads + 1)
        {
            unusedReads.add(readId);
        }
    }

    protected void filterAlignments()
    {
        mixin(traceExecution);

        auto filters = tuple(
            new WeaklyAnchoredAlignmentChainsFilter(repetitiveRegions, options.minAnchorLength),
            new ImproperAlignmentChainsFilter(),
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
            assert(isSorted(readsAlignment));

            logJsonDiagnostic(
                "filterStage", typeof(filter).stringof,
                "readsAlignment", shouldLog(LogLevel.debug_)
                    ? readsAlignment.toJson
                    : toJson(null),
                "numAlignmentChains", readsAlignment.count!"!a.flags.disabled",
            );
        }
    }

    protected PileUp[] buildPileUps()
    {
        mixin(traceExecution);

        import dentist.commands.collectPileUps.pileups : build;
        auto pileUps = build(
            numReferenceContigs,
            readsAlignment,
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

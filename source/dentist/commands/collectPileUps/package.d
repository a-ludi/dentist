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
import dentist.commands.collectPileUps.maskassessment :
    assessmentStage,
    BadAlignmentCoverageAssessor;
import dentist.common :
    DentistException,
    ReferenceRegion;
import dentist.common.alignments :
    AlignmentChain,
    alignmentCoverage,
    getAlignmentRefs,
    getType,
    PileUp;
import dentist.common.binio : writePileUpsDb;
import dentist.dazzler :
    attachTracePoints,
    getAlignments,
    getNumContigs,
    writeMask;
import dentist.util.log;
import dentist.util.math : ceildiv, NaturalNumberSet;
import dstats.distrib : invPoissonCDF;
import std.algorithm : count, isSorted, map, sort, sum;
import std.array : array;
import std.conv : to;
import std.exception : enforce;
import std.format : format;
import std.parallelism : parallel, taskPool;
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
    protected AlignmentChain[] selfAlignment;
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
        assessRepeatStructure();
        writeMask(
            options.refDb,
            options.repeatMask,
            this.repetitiveRegions.intervals,
            options.workdir,
        );
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
        selfAlignment = getAlignments(
            options.refDb,
            options.selfAlignmentFile,
            options.workdir,
        );
        readsAlignment = getAlignments(
            options.refDb,
            options.readsDb,
            options.readsAlignmentFile,
            options.workdir,
        );

        enforce!DentistException(selfAlignment.length > 0, "empty self-alignment");
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

    protected void assessRepeatStructure()
    {
        mixin(traceExecution);

        auto selfCoverageInfo = getCoverageInfo(selfAlignment);
        auto readsCoverageInfo = getCoverageInfo(readsAlignment);
        logJsonDebug(
            "selfCoverage", selfCoverageInfo.intrinsicCoverage,
            "selfCoverageConfidenceInterval", selfCoverageInfo.confidenceInterval.toJson,
            "readsCoverage", readsCoverageInfo.intrinsicCoverage,
            "readsCoverageConfidenceInterval", readsCoverageInfo.confidenceInterval.toJson,
        );

        auto assessmentStages = tuple(
            assessmentStage(
                "self-alignment",
                new BadAlignmentCoverageAssessor(selfCoverageInfo.confidenceInterval.expand),
                selfAlignment,
            ),
            assessmentStage(
                "reads-alignment",
                new BadAlignmentCoverageAssessor(readsCoverageInfo.confidenceInterval.expand),
                readsAlignment,
            ),
        );

        foreach (stage; assessmentStages)
        {
            auto repetitiveRegions = stage.assessor(stage.input);

            logJsonDiagnostic(
                "assessor", stage.name,
                "repetitiveRegions", shouldLog(LogLevel.debug_)
                    ? repetitiveRegions.intervals.toJson
                    : toJson(null),
                "numRepetitiveRegions", repetitiveRegions.intervals.length,
            );

            if (shouldLog(LogLevel.debug_))
            {
                auto maskName = format!"%s-%s"(options.repeatMask, stage.name);

                writeMask(
                    options.refDb,
                    maskName,
                    repetitiveRegions.intervals,
                    options.workdir,
                );
            }

            this.repetitiveRegions |= repetitiveRegions;
        }

        logJsonDiagnostic(
            "assessor", "finalResult",
            "repetitiveRegions", shouldLog(LogLevel.debug_)
                ? this.repetitiveRegions.intervals.toJson
                : toJson(null),
            "numRepetitiveRegions", this.repetitiveRegions.intervals.length,
        );
    }

    protected auto getCoverageInfo(in AlignmentChain[] alignments)
    {
        auto alphaHalf = (1 - options.confidence) / 2;
        auto intrinsicCoverage = alignmentCoverage(alignments);
        auto confidenceInterval = tuple(
            invPoissonCDF(alphaHalf, intrinsicCoverage),
            invPoissonCDF(1 - alphaHalf, intrinsicCoverage),
        );

        return tuple!(
            "intrinsicCoverage",
            "confidenceInterval",
        )(
            intrinsicCoverage,
            confidenceInterval,
        );
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
        auto pileUps = build(numReferenceContigs, readsAlignment);

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

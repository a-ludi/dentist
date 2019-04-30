/**
    This is the `translocateGaps` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.translocateGaps;

import dentist.common : isTesting;

static if (isTesting):

import dentist.commandline : OptionsFor;
import dentist.common :
    ReferenceInterval,
    ReferenceRegion,
    toInterval;
import dentist.common.alignments :
    AlignmentChain,
    id_t;
import dentist.common.commands : TestingCommand;
import dentist.dazzler :
    getAlignments,
    getNumContigs,
    getFastaSequence,
    writeMask;
import dentist.util.log;
import dentist.util.range : wrapLines;
import std.algorithm :
    cache,
    copy,
    filter,
    joiner,
    map;
import std.array : array;
import std.format : format;
import std.range :
    assumeSorted,
    chain,
    chunks,
    enumerate,
    iota,
    only,
    repeat,
    slide,
    takeExactly;
import std.stdio : File, stdout;
import std.traits : Unqual;
import std.typecons : No;
import vibe.data.json : toJson = serializeToJson;


/// Options for the `collectPileUps` command.
alias Options = OptionsFor!(TestingCommand.translocateGaps);

/// Execute the `translocateGaps` command with `options`.
void execute(in Options options)
{
    auto translocator = Translocator(options);

    return translocator.run();
}

private struct Translocator
{
    alias FastaWriter = typeof(wrapLines(stdout.lockingTextWriter, 0));

    const(Options) options;
    size_t numContigsAssembly2;
    AlignmentChain[] alignments;
    ReferenceRegion mappedRegions;
    File resultFile;
    FastaWriter writer;

    this(in Options options)
    {
        this.options = options;
        this.resultFile = options.resultFile is null
            ? stdout
            : File(options.resultFile, "w");
        this.writer = wrapLines(resultFile.lockingTextWriter, options.fastaLineWidth);
    }

    void run()
    {
        mixin(traceExecution);

        init();
        writeOutputAssembly();
    }

    protected void init()
    {
        mixin(traceExecution);

        alignments = getAlignments(
            options.trueAssemblyDb,
            options.shortReadAssemblyDb,
            options.shortReadAssemblyAlignmentFile,
            options.workdir,
        );
        mappedRegions = ReferenceRegion(alignments
            .filter!"a.isProper"
            .map!(ac => ac.toInterval!(ReferenceInterval, "contigA"))
            .filter!"a.size > 0"
            .array);
        if (mappedRegions.intervals.length > 0)
            mappedRegions = ReferenceRegion(mappedRegions
                .intervals
                // Remove small gaps
                .slide!(No.withPartial)(2)
                .map!(gapPair => gapPair[0].contigId == gapPair[1].contigId
                    ? (gapPair[1].begin - gapPair[0].end >= options.minGapSize
                        ? gapPair[0]
                        : ReferenceInterval.convexHull(gapPair[0], gapPair[1])
                    )
                    : gapPair[0])
                .chain(only(mappedRegions.intervals[$ - 1]))
                .map!(mappedInterval => cast(Unqual!(typeof(mappedInterval))) mappedInterval)
                .array
            );
        if (mappedRegions.intervals.length > 0)
            mappedRegions = ReferenceRegion(mappedRegions
                .intervals
                // Remove small contigs
                .filter!(mappedInterval => mappedInterval.size >= options.minContigSize)
                .map!(mappedInterval => cast(Unqual!(typeof(mappedInterval))) mappedInterval)
                .array
            );

        writeMask(
            options.trueAssemblyDb,
            options.mappedRegionsMask,
            mappedRegions.intervals,
            options.workdir,
        );
    }

    protected void writeOutputAssembly()
    {
        mixin(traceExecution);

        enum dchar unknownBase = 'n';
        auto numRefContigs = getNumContigs(options.trueAssemblyDb, options.workdir);
        auto mappedRegions = mappedRegions.intervals.assumeSorted!"a.contigId < b.contigId";
        ReferenceInterval needle;

        foreach (id_t contigId; 1 .. numRefContigs + 1)
        {
            needle.contigId = contigId;
            auto contigMappedRegions = mappedRegions.equalRange(needle);

            if (contigMappedRegions.length == 0)
                continue;

            auto contigSequence = getFastaSequence(
                options.trueAssemblyDb,
                contigId,
                options.workdir,
            );
            // Prepend needle to produce the first contig.
            auto keepRegionsPairs = chain(only(needle), contigMappedRegions)
                .slide!(No.withPartial)(2);

            getScaffoldHeader(contigId).copy(writer);
            foreach (keepRegions; keepRegionsPairs)
            {
                if (keepRegions[0].end > 0)
                    repeat(unknownBase)
                        .takeExactly(keepRegions[1].begin - keepRegions[0].end)
                        .copy(writer);

                contigSequence[keepRegions[1].begin .. keepRegions[1].end].copy(writer);
            }

            "\n".copy(writer);
        }
    }

    static protected string getScaffoldHeader(in size_t scaffoldId)
    {
        return format!">translocated_gaps_%d\n"(scaffoldId);
    }
}

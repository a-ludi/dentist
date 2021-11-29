/**
    This is the `filterMask` command of DENTIST.

    Command_Summary:

    ---
    Filter a Dazzler mask. See `dentist filter-mask --help` for a list of
    available filters.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.filterMask;

package(dentist) enum summary = "
    Filter a Dazzler mask. See `dentist filter-mask --help` for a list of
    available filters.
";

import dentist.commandline : OptionsFor;
import dentist.common :
    ReferenceInterval,
    ReferenceRegion,
    toInterval;
import dentist.common.alignments :
    AlignmentChain,
    id_t;
import dentist.common.commands : DentistCommand;
import dentist.dazzler :
    readMask,
    writeMask;
import dentist.util.log;
import dentist.util.range : wrapLines;
import dentist.util.algorithm : filterInPlace;
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


/// Options for the `filterMask` command.
alias Options = OptionsFor!(DentistCommand.filterMask);


/// Execute the `filterMask` command with `options`.
void execute(in Options options)
{
    auto builder = PartialAssemblyBuilder(options);

    return builder.run();
}

private struct PartialAssemblyBuilder
{
    alias FastaWriter = typeof(wrapLines(stdout.lockingTextWriter, 0));

    const(Options) options;
    ReferenceRegion mask;

    this(in Options options)
    {
        this.options = options;
    }

    void run()
    {
        mixin(traceExecution);

        readInputMask();
        runFilters();
        writeOutputMask();
    }

    protected void readInputMask()
    {
        mixin(traceExecution);

        mask = ReferenceRegion(readMask!ReferenceInterval(
            options.refDb,
            options.inMask,
        ));
    }


    protected void runFilters()
    {
        removeSmallGaps();
        removeSmallContigs();
    }


    protected void writeOutputMask()
    {
        mixin(traceExecution);

        writeMask(
            options.refDb,
            options.outMask,
            mask.intervals,
        );
    }

    protected void removeSmallGaps()
    {
        scope mappedIntervals = mask.releaseIntervals();

        size_t accIdx;
        foreach (i, currentInterval; mappedIntervals)
        {
            auto accInterval = &mappedIntervals[accIdx];

            if (
                accInterval.contigId == currentInterval.contigId &&
                accInterval.end + options.minGapSize >= currentInterval.begin
            )
                // extend accInterval
                accInterval.end = currentInterval.end;
            else if (accIdx + 1 < mappedIntervals.length)
                // start next accInterval
                mappedIntervals[++accIdx] = currentInterval;
        }

        mappedIntervals = mappedIntervals[0 .. accIdx + 1];

        mask = ReferenceRegion(mappedIntervals);
    }


    protected void removeSmallContigs()
    {
        scope mappedIntervals = mask.releaseIntervals();

        mappedIntervals = filterInPlace!(
            mappedInterval => mappedInterval.size >= options.minIntervalSize
        )(mappedIntervals);

        mask = ReferenceRegion(mappedIntervals);
    }
}

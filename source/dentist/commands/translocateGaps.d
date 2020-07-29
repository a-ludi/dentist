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
import dentist.common.alignments : AlignmentChain;
import dentist.common.commands : TestingCommand;
import dentist.dazzler :
    getAlignments,
    writeMask;
import dentist.util.log;
import std.algorithm :
    filter,
    map;
import std.array : array;


/// Options for the `translocateGaps` command.
alias Options = OptionsFor!(TestingCommand.translocateGaps);

/// Execute the `translocateGaps` command with `options`.
void execute(in Options options)
{
    auto translocator = GapTranslocator(options);

    return translocator.run();
}

private struct GapTranslocator
{
    const(Options) options;
    AlignmentChain[] alignments;
    ReferenceRegion mappedRegions;

    this(in Options options)
    {
        this.options = options;
    }

    void run()
    {
        mixin(traceExecution);

        readAlignments();
        buildMask();
        writeMask(
            options.trueAssemblyDb,
            options.mappedRegionsMask,
            mappedRegions.intervals,
        );
    }

    protected void readAlignments()
    {
        mixin(traceExecution);

        alignments = getAlignments(
            options.trueAssemblyDb,
            options.shortReadAssemblyDb,
            options.shortReadAssemblyAlignmentFile,
        );
    }

    protected void buildMask()
    {
        mixin(traceExecution);

        mappedRegions = ReferenceRegion(alignments
            .filter!(a => a.isProper(options.properAlignmentAllowance))
            .map!(ac => ac.toInterval!(ReferenceInterval, "contigA"))
            .filter!"a.size > 0"
            .array);
    }
}

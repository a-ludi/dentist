/**
    This is the `chainLocalAlignments` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.chainLocalAlignments;

import dentist.commandline : OptionsFor;
import dentist.common.alignments :
    AlignmentChain,
    ChainingOptions,
    chainLocalAlignments;
import dentist.common.commands : DentistCommand;
import dentist.dazzler :
    AlignmentHeader,
    BufferMode,
    getFlatLocalAlignments,
    writeAlignments;
import dentist.util.log;
import std.algorithm : sort;
import std.range : tee;
import std.typecons : Yes;


/// Options for the `chainLocalAlignments` command.
alias Options = OptionsFor!(DentistCommand.chainLocalAlignments);

/// Execute the `chainLocalAlignments` command with `options`.
void execute(in Options options)
{
    auto chainer = new CLIChainer(options);

    return chainer.run();
}


class CLIChainer
{
    alias FlatLocalAlignments = typeof(getFlatLocalAlignments("", ""));

    const Options options;
    const ChainingOptions chainingOptions;
    FlatLocalAlignments alignments;
    ProgressMeter progress;


    this(const Options options)
    {
        this.options = options;
        this.chainingOptions = options.chainingOptions;
        this.progress = options.createProgressMeter();
    }


    void run()
    {
        mixin(traceExecution);

        readAlignments();
        chainLocalAlignments();
    }


    void readAlignments()
    {
        mixin(traceExecution);

        alignments = getFlatLocalAlignments(
            options.refDb,
            options.readsDb,
            options.dbAlignmentFile,
            BufferMode.preallocated,
        );
    }


    void chainLocalAlignments()
    {
        mixin(traceExecution);

        progress.start();

        auto chainedAlignments = alignments
            .chainLocalAlignments(chainingOptions)
            .tee!(_ => progress.tick());

        options.chainedAlignments.writeAlignments(chainedAlignments);

        progress.stop();
    }
}

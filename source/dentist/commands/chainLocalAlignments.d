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
    chainLocalAlignments,
    cmpIdsAndComplement;
import dentist.common.commands : DentistCommand;
import dentist.dazzler :
    AlignmentHeader,
    getAlignments,
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
    const Options options;
    AlignmentChain[] alignments;
    AlignmentHeader headerData;
    ProgressMeter progress;


    this(const Options options)
    {
        this.options = options;
        this.progress = options.createProgressMeter();
    }


    void run()
    {
        mixin(traceExecution);

        readAlignments();
        sortAlignments();
        inferHeaderData();
        chainLocalAlignments();
    }


    void readAlignments()
    {
        mixin(traceExecution);

        alignments = getAlignments(
            options.refDb,
            options.readsDb,
            options.dbAlignmentFile,
            Yes.includeTracePoints
        );
    }


    void sortAlignments()
    {
        mixin(traceExecution);

        // Order required for faster chaining because local aligments with
        // different complements must not be chained.
        sort!((a, b) => cmpIdsAndComplement(a, b) < 0)(alignments);
    }


    void inferHeaderData()
    {
        mixin(traceExecution);

        headerData = AlignmentHeader.inferFrom(alignments);
    }


    void chainLocalAlignments()
    {
        mixin(traceExecution);

        progress.start();

        auto chainedAlignments = alignments
            .chainLocalAlignments()
            .tee!(_ => progress.tick());

        options.chainedAlignments.writeAlignments(chainedAlignments, headerData);

        progress.stop();
    }
}

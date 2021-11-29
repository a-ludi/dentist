/**
    This is the `buildPartialAssembly` command of DENTIST.

    Command_Summary:

    ---
    Build a partial assembly from a mask.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.buildPartialAssembly;

package(dentist) enum summary = "Build a partial assembly from a mask.";

import dentist.common : isTesting;

static if (isTesting):

import dentist.commandline : OptionsFor;
import dentist.common :
    ReferenceInterval,
    ReferenceRegion;
import dentist.common.alignments : id_t;
import dentist.common.commands : TestingCommand;
import dentist.dazzler :
    getNumContigs,
    getFastaSequence,
    readMask;
import dentist.util.log;
import dentist.util.range : wrapLines;
import std.algorithm : copy;
import std.format : format;
import std.range :
    assumeSorted,
    chain,
    only,
    repeat,
    slide,
    takeExactly;
import std.stdio : File, stdout;
import std.typecons : No;


/// Options for the `buildPartialAssembly` command.
alias Options = OptionsFor!(TestingCommand.buildPartialAssembly);


/// Execute the `buildPartialAssembly` command with `options`.
void execute(in Options options)
{
    auto translocator = Translocator(options);

    return translocator.run();
}

private struct Translocator
{
    alias FastaWriter = typeof(wrapLines(stdout.lockingTextWriter, 0));

    const(Options) options;
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

        mappedRegions = ReferenceRegion(readMask!ReferenceInterval(
            options.trueAssemblyDb,
            options.mappedRegionsMask,
        ));
        writeOutputAssembly();
    }

    protected void writeOutputAssembly()
    {
        mixin(traceExecution);

        enum dchar unknownBase = 'n';
        auto numRefContigs = getNumContigs(options.trueAssemblyDb);
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
            );
            // Prepend needle to produce the first contig.
            auto mappedRegionsPairs = chain(only(needle), contigMappedRegions)
                .slide!(No.withPartial)(2);

            getScaffoldHeader(contigId).copy(writer);
            foreach (mappedPair; mappedRegionsPairs)
            {
                if (mappedPair[0].end > 0)
                    repeat(unknownBase)
                        .takeExactly(mappedPair[1].begin - mappedPair[0].end)
                        .copy(writer);

                contigSequence[mappedPair[1].begin .. mappedPair[1].end].copy(writer);
            }

            "\n".copy(writer);
        }
    }

    static protected string getScaffoldHeader(in size_t scaffoldId)
    {
        return format!">translocated_gaps_%d\n"(scaffoldId);
    }
}

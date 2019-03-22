/**
    This is the `generateDazzlerOptions` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.generateDazzlerOptions;

import dentist.common : isTesting;
import std.stdio : writefln, writeln;

/// Execute the `generateDazzlerOptions` command with `options`.
void execute(Options)(in Options options)
{
    static if (isTesting)
    {
        writeln("# short vs true assembly alignment options (consider using `HPC.daligner`)");
        writefln!"daligner %-(%s %) <true-assembly> <short-read-assembly>"(options.shortVsTrueAssemblyAlignmentOptions);
    }

    writeln("# self alignment options (consider using `HPC.daligner`)");
    writefln!"daligner %-(%s %) <reference> <reference>"(options.selfAlignmentOptions);
    writeln("# ref vs reads alignment options (consider using `HPC.damapper`)");
    writefln!"damapper %-(%s %) <reference> <reads>"(options.refVsReadsAlignmentOptions);
}

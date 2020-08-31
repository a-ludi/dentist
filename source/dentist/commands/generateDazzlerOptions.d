/**
    This is the `generateDazzlerOptions` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.generateDazzlerOptions;

import dentist.common : isTesting;
import dentist.common.external : ExternalDependency;
import std.stdio : writefln, writeln;

/// Execute the `generateDazzlerOptions` command with `options`.
void execute(Options)(in Options options)
{
    writeln("# self alignment options (consider using `HPC.daligner`)");
    writefln!"daligner %-(%s %) <reference> <reference>"(options.selfAlignmentOptions);
    writeln("# tandem alignment options (consider using `HPC.daligner`)");
    writefln!"datander %-(%s %) <reference>"(options.tandemAlignmentOptions);
    writeln("# ref vs reads alignment options (consider using `HPC.daligner`)");
    writefln!"daligner %-(%s %) <reference> <reads>"(options.refVsReadsAlignmentOptions);
}

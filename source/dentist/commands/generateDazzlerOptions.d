/**
    This is the `generate-dazzler-options` command of DENTIST.

    Command_Summary:

    ---
    Outputs advice on how to call some of the Dazzler tools. The output is
    used by the Snakemake workflow to generate the correct commands.
    ---

    Example:
    ---
    # DBdust reference
    DBdust  <reference>
    # self alignment options (consider using `HPC.daligner`)
    daligner -I -T{threads} -l500 -e0.7 -mdust -P/tmp/ludwig <reference> <reference>
    # tandem alignment options (consider using `HPC.daligner`)
    datander -T{threads} -s126 -l500 -e0.7 -P/tmp/ludwig <reference>
    # ref vs reads alignment options (consider using `HPC.damapper`)
    damapper -C -T{threads} -e0.7 -P/tmp/ludwig <reference> <reads>
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.generateDazzlerOptions;

package(dentist) enum summary = "
    Outputs advice on how to call some of the Dazzler tools. The output is
    used by the Snakemake workflow to generate the correct commands.
";

import dentist.commandline : OptionsFor;
import dentist.common.commands : DentistCommand;
import dentist.common.external : ExternalDependency;
import std.stdio : writefln, writeln;


/// Options for the `generate-dazzler-options` command.
alias Options = OptionsFor!(DentistCommand.generateDazzlerOptions);


/// Execute the `generate-dazzler-options` command with `options`.
void execute(in Options options)
{
    writeln("# DBdust reference");
    writefln!"DBdust %-(%s %) <reference>"(options.additionalRefDustOptions);
    writeln("# self alignment options (consider using `HPC.daligner`)");
    writefln!"daligner %-(%s %) <reference> <reference>"(options.selfAlignmentOptions);
    writeln("# tandem alignment options (consider using `HPC.daligner`)");
    writefln!"datander %-(%s %) <reference>"(options.tandemAlignmentOptions);
    writeln("# ref vs reads alignment options (consider using `HPC.damapper`)");
    writefln!"damapper %-(%s %) <reference> <reads>"(options.refVsReadsAlignmentOptions);
}

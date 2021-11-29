/**
    This is the `mergeMasks` command of DENTIST.

    Command_Summary:

    ---
    Merge several masks into a single one with non-overlapping intervals. This
    can also be used to merge overlapping intervals in a single mask.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.mergeMasks;

package(dentist) enum summary = "
    Merge several masks into a single one with non-overlapping intervals. This
    can also be used to merge overlapping intervals in a single mask.
";

import dentist.commandline : OptionsFor;
import dentist.common :
    ReferenceInterval,
    ReferenceRegion;
import dentist.common.alignments :
    coord_t,
    id_t,
    FlatLocalAlignment;
import dentist.common.commands : DentistCommand;
import dentist.util.log;
import dentist.dazzler :
    readMask,
    writeMask;


/// Options for the `mergeMasks` command.
alias Options = OptionsFor!(DentistCommand.mergeMasks);


/// Execute the `mergeMasks` command with `options`.
void execute(in Options options)
{
    mixin(traceExecution);

    ReferenceRegion outMask;

    foreach (inMask; options.inMasks)
        outMask |= ReferenceRegion(readMask!ReferenceInterval(
            options.refDb,
            inMask,
        ));

    writeMask(options.refDb, options.outMask, outMask.intervals);
}

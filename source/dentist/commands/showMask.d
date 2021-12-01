/**
    This is the `show-mask` command of DENTIST.

    Command_Summary:

    ---
    Show a short summary of a mask. Can also be used to convert a mask to JSON
    by increasing verbosity two times.
    ---

    Example:

    Table output (default).

    ---
    name:                      dust
    numIntervals:           3043971
    numMaskedBases:       141388227

    name:                     tan-H
    numIntervals:             27765
    numMaskedBases:       124147972

    name:            dentist-self-H
    numIntervals:            212770
    numMaskedBases:       532462661

    name:           dentist-reads-H
    numIntervals:             31679
    numMaskedBases:        35656876

    name:                __merged__
    numIntervals:           2612751
    numMaskedBases:       759305412
    ---

    Example:

    JSON output (`--json`, `-j`).

    ---
    [
        {
            "name": "dust",
            "numIntervals": 3043971,
            "numMaskedBases": 141388227
        },
        {
            "name": "tan-H",
            "numIntervals": 27765,
            "numMaskedBases": 124147972
        },
        {
            "name": "dentist-self-H",
            "numIntervals": 212770,
            "numMaskedBases": 532462661
        },
        {
            "name": "dentist-reads-H",
            "numIntervals": 31679,
            "numMaskedBases": 35656876
        },
        {
            "name": "__merged__",
            "numIntervals": 2612751,
            "numMaskedBases": 759305412
        }
    ]
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.showMask;

package(dentist) enum summary = "
    Show a short summary of a mask. Can also be used to convert a mask to JSON
    by increasing verbosity two times.
";

import dentist.commandline : OptionsFor;
import dentist.common :
    ReferenceInterval,
    ReferenceRegion;
import dentist.common.commands : DentistCommand;
import dentist.util.log;
import dentist.dazzler : readMask;
import std.algorithm :
    map,
    max,
    maxElement,
    sum;
import std.math : log10, lrint, FloatingPointControl;
import std.stdio : writefln, writeln, stderr;
import vibe.data.json :
    serializeToJsonString,
    toJson = serializeToJson,
    toJsonString = serializeToPrettyJson;


/// Options for the `show-mask` command.
alias Options = OptionsFor!(DentistCommand.showMask);


/// Execute the `show-mask` command with `options`.
void execute(in Options options)
{
    auto masks = options.masks;
    ReferenceRegion mergedMask;
    Stats[] statsList;

    foreach (mask; masks)
    {
        auto maskRegion = ReferenceRegion(readMask!ReferenceInterval(
            options.refDb,
            mask,
        ));

        if (shouldLog(LogLevel.debug_))
            stderr.writeln(serializeToJsonString(maskRegion.intervals));
        statsList ~= statsFor(mask, maskRegion);

        if (masks.length > 1)
            mergedMask |= maskRegion;
    }

    if (masks.length > 1)
    {
        if (shouldLog(LogLevel.debug_))
            stderr.writeln(serializeToJsonString(mergedMask.intervals));
        statsList ~= statsFor("__merged__", mergedMask);
    }

    if (options.useJson)
        writeln(statsList.toJsonString);
    else
        writeTabular(statsList);
}

struct Stats
{
    string name;
    size_t numIntervals;
    size_t numMaskedBases;

    size_t columnWidth() const nothrow
    {
        FloatingPointControl fpCtrl;
        fpCtrl.rounding = FloatingPointControl.roundUp;
        auto numWidth = lrint(log10(max(
            numIntervals,
            numMaskedBases,
        )));
        numWidth = max(numWidth, name.length);

        return numWidth;
    }
}

Stats statsFor(string name, ReferenceRegion maskRegion)
{
    return Stats(
        name,
        maskRegion.intervals.length,
        maskRegion.intervals.map!"a.size".sum,
    );
}

void writeTabular(Stats[] statsList)
{
    auto columnWidth = statsList.map!"a.columnWidth".maxElement;

    foreach (i, stats; statsList)
    {
        writefln!"name:           %*s"(columnWidth, stats.name);
        writefln!"numIntervals:   %*d"(columnWidth, stats.numIntervals);
        writefln!"numMaskedBases: %*d"(columnWidth, stats.numMaskedBases);

        if (i < statsList.length - 1)
            writeln();
    }
}

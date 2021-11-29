/**
    This is the `showPileUps` command of DENTIST.

    Command_Summary:

    ---
    Show a short summary of the pile ups. Can also be used to convert pile υps
    to JSON by increasing verbosity two times.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.showPileUps;

package(dentist) enum summary = "
    Show a short summary of the pile ups. Can also be used to convert pile υps
    to JSON by increasing verbosity two times.
";

import dentist.commandline : OptionsFor;
import dentist.common.alignments : getType;
import dentist.common.binio : PileUpDb;
import dentist.common.commands : DentistCommand;
import dentist.util.log;
import std.algorithm : map, max;
import std.array : array;
import std.conv : to;
import std.file : getSize;
import std.math : log10, lrint, FloatingPointControl;
import std.stdio : stderr, writefln, writeln;
import std.typecons : tuple;
import vibe.data.json :
    toJson = serializeToJson,
    toJsonString = serializeToPrettyJson,
    serializeToJsonString;


/// Options for the `showPileUps` command.
alias Options = OptionsFor!(DentistCommand.showPileUps);


/// Execute the `showPileUps` command with `options`.
void execute(in Options options)
{
    size_t totalDbSize = options.pileUpsFile.getSize();
    auto pileUpDb = PileUpDb.parse(options.pileUpsFile);
    if (!shouldLog(LogLevel.debug_))
        pileUpDb.releaseDb();

    auto stats = Stats(
        totalDbSize,
        pileUpDb.pileUps.length,
        pileUpDb.readAlignments.length,
        pileUpDb.seededAlignments.length,
        pileUpDb.localAlignments.length,
        pileUpDb.tracePoints.length,
    );

    if (options.useJson)
        writeln(stats.toJsonString);
    else
        writeTabular(stats);

    if (shouldLog(LogLevel.debug_))
    {
        foreach (pileUp; pileUpDb[])
        {
            stderr.writeln(serializeToJsonString([
                "type": pileUp.getType.to!string.toJson,
                "readAlignments": pileUp.map!"a[]".array.toJson,
            ]));
        }
    }
}

struct Stats
{
    size_t totalDbSize;
    size_t numPileUps;
    size_t numReadAlignments;
    size_t numSeededAlignments;
    size_t numLocalAlignments;
    size_t numTracePoints;

    size_t columnWidth() const nothrow
    {
        FloatingPointControl fpCtrl;
        fpCtrl.rounding = FloatingPointControl.roundUp;
        auto numWidth = lrint(log10(max(
            totalDbSize,
            numPileUps,
            numReadAlignments,
            numSeededAlignments,
            numLocalAlignments,
            numTracePoints,
        )));

        return numWidth;
    }
}

void writeTabular(Stats stats)
{
    auto columnWidth = stats.columnWidth();

    writefln!"totalDbSize:         %*d bytes"(columnWidth, stats.totalDbSize);
    writefln!"numPileUps:          %*d"(columnWidth, stats.numPileUps);
    writefln!"numReadAlignments:   %*d"(columnWidth, stats.numReadAlignments);
    writefln!"numSeededAlignments: %*d"(columnWidth, stats.numSeededAlignments);
    writefln!"numLocalAlignments:  %*d"(columnWidth, stats.numLocalAlignments);
    writefln!"numTracePoints:      %*d"(columnWidth, stats.numTracePoints);
}

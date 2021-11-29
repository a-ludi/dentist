/**
    This is the `showInsertions` command of DENTIST.

    Command_Summary:

    ---
    Show a short summary of the insertions. Can also be used to convert
    insertions to JSON by increasing verbosity two times.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.showInsertions;

package(dentist) enum summary = "
    Show a short summary of the insertions. Can also be used to convert
    insertions to JSON by increasing verbosity two times.
";

import dentist.commandline : OptionsFor;
import dentist.common.binio : InsertionDb;
import dentist.common.commands : DentistCommand;
import dentist.util.log;
import std.algorithm : map, max;
import std.array : array;
import std.file : getSize;
import std.math : log10, lrint, FloatingPointControl;
import std.stdio : stderr, writefln, writeln;
import std.typecons : tuple;
import vibe.data.json :
    serializeToJsonString,
    toJson = serializeToJson,
    toJsonString = serializeToPrettyJson;


/// Options for the `showInsertions` command.
alias Options = OptionsFor!(DentistCommand.showInsertions);


/// Execute the `showInsertions` command with `options`.
void execute(in Options options)
{
    size_t totalDbSize = options.insertionsFile.getSize();
    auto insertionDb = InsertionDb.parse(options.insertionsFile);
    if (!shouldLog(LogLevel.debug_))
        insertionDb.releaseDb();

    auto stats = Stats(
        totalDbSize,
        insertionDb.insertions.length,
        insertionDb.compressedBaseQuads.length,
        insertionDb.overlaps.length,
        insertionDb.localAlignments.length,
        insertionDb.tracePoints.length,
        insertionDb.readIds.length,
    );

    if (options.useJson)
        writeln(stats.toJsonString);
    else
        writeTabular(stats);

    if (shouldLog(LogLevel.debug_))
    {
        foreach (join; insertionDb[])
        {
            stderr.writeln(serializeToJsonString([
                "start": join.start.toJson,
                "end": join.end.toJson,
                "payload": [
                    "sequence": join.payload.sequence.to!string.toJson,
                    "contigLength": join.payload.contigLength.toJson,
                    "overlaps": join.payload.overlaps.toJson,
                    "readIds": join.payload.readIds.toJson,
                ].toJson,
            ]));
        }
    }
}

struct Stats
{
    size_t totalDbSize;
    size_t numInsertions;
    size_t numCompressedBaseQuads;
    size_t numOverlaps;
    size_t numLocalAlignments;
    size_t numTracePoints;
    size_t numReadIds;

    size_t columnWidth() const nothrow
    {
        FloatingPointControl fpCtrl;
        fpCtrl.rounding = FloatingPointControl.roundUp;
        auto numWidth = lrint(log10(max(
            totalDbSize,
            numInsertions,
            numCompressedBaseQuads,
            numOverlaps,
            numLocalAlignments,
            numTracePoints,
            numReadIds,
        )));

        return numWidth;
    }
}

void writeTabular(Stats stats)
{
    auto columnWidth = stats.columnWidth();

    writefln!"totalDbSize:            %*d bytes"(columnWidth, stats.totalDbSize);
    writefln!"numInsertions:          %*d"(columnWidth, stats.numInsertions);
    writefln!"numCompressedBaseQuads: %*d"(columnWidth, stats.numCompressedBaseQuads);
    writefln!"numOverlaps:            %*d"(columnWidth, stats.numOverlaps);
    writefln!"numLocalAlignments:     %*d"(columnWidth, stats.numLocalAlignments);
    writefln!"numTracePoints:         %*d"(columnWidth, stats.numTracePoints);
    writefln!"numReadIds:             %*d"(columnWidth, stats.numReadIds);
}

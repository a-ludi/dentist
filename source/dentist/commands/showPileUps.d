/**
    This is the `showPileUps` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.showPileUps;

import dentist.common.binio : PileUpDb;
import dentist.common.alignments : getType;
import dentist.util.log;
import std.algorithm : map, max;
import std.array : array;
import std.conv : to;
import std.file : getSize;
import std.math : log10, lrint, FloatingPointControl;
import std.stdio : writefln, writeln;
import std.typecons : tuple;
import vibe.data.json : toJson = serializeToJson, toJsonString = serializeToPrettyJson;

/// Execute the `showPileUps` command with `options`.
void execute(Options)(in Options options)
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

    logJsonDebug("pileUps", pileUpDb[]
        .map!(pileUp => [
            "type": pileUp.getType.to!string.toJson,
            "readAlignments": pileUp.map!"a[]".array.toJson,
        ].toJson)
        .array
        .toJson);
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

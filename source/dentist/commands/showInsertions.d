/**
    This is the `showInsertions` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.showInsertions;

import dentist.common.binio : InsertionDb;
import dentist.util.log;
import std.algorithm : map, max;
import std.array : array;
import std.file : getSize;
import std.math : log10, lrint, FloatingPointControl;
import std.stdio : writefln, writeln;
import std.typecons : tuple;
import vibe.data.json : toJson = serializeToJson, toJsonString = serializeToPrettyJson;

/// Execute the `showInsertions` command with `options`.
void execute(Options)(in Options options)
{
    size_t totalDbSize = options.insertionsFile.getSize();
    auto insertionDb = InsertionDb.parse(options.insertionsFile);
    if (!shouldLog(LogLevel.debug_))
        insertionDb.releaseDb();

    auto stats = Stats(
        totalDbSize,
        insertionDb.insertions.length,
        insertionDb.compressedBaseQuads.length,
        insertionDb.spliceSites.length,
    );

    if (options.useJson)
        writeln(stats.toJsonString);
    else
        writeTabular(stats);

    logJsonDebug("insertions", insertionDb[]
        .map!(join => [
            "start": join.start.toJson,
            "end": join.end.toJson,
            "payload": [
                "sequence": join.payload.sequence.to!string.toJson,
                "contigLength": join.payload.contigLength.toJson,
                "spliceSites": join.payload.spliceSites.toJson,
            ].toJson,
        ])
        .array
        .toJson);
}

struct Stats
{
    size_t totalDbSize;
    size_t numInsertions;
    size_t numCompressedBaseQuads;
    size_t numSpliceSites;

    size_t columnWidth() const nothrow
    {
        FloatingPointControl fpCtrl;
        fpCtrl.rounding = FloatingPointControl.roundUp;
        auto numWidth = lrint(log10(max(
            totalDbSize,
            numInsertions,
            numCompressedBaseQuads,
            numSpliceSites,
        )));
        numWidth += numWidth / 3;

        return numWidth;
    }
}

void writeTabular(Stats stats)
{
    auto columnWidth = stats.columnWidth();

    writefln!"totalDbSize:            %*d bytes"(columnWidth, stats.totalDbSize);
    writefln!"numInsertions:          %*d"(columnWidth, stats.numInsertions);
    writefln!"numCompressedBaseQuads: %*d"(columnWidth, stats.numCompressedBaseQuads);
    writefln!"numSpliceSites:         %*d"(columnWidth, stats.numSpliceSites);
}

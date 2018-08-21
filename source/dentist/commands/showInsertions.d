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
import vibe.data.json : toJson = serializeToJson;

/// Execute the `showInsertions` command with `options`.
void execute(Options)(in Options options)
{
    size_t totalDbSize = options.insertionsFile.getSize();
    auto insertionDb = InsertionDb.parse(options.insertionsFile);
    if (!shouldLog(LogLevel.debug_))
        insertionDb.releaseDb();

    auto stats = tuple(
        totalDbSize,
        insertionDb.insertions.length,
        insertionDb.compressedBaseQuads.length,
        insertionDb.spliceSites.length,
    );

    FloatingPointControl fpCtrl;
    fpCtrl.rounding = FloatingPointControl.roundUp;
    auto numWidth = lrint(log10(max(stats.expand)));
    numWidth += numWidth / 3;

    writefln!"totalDbSize:            %*,d bytes"(numWidth, stats[0]);
    writefln!"numInsertions:          %*,d"(numWidth, stats[1]);
    writefln!"numCompressedBaseQuads: %*,d"(numWidth, stats[2]);
    writefln!"numSpliceSites:         %*,d"(numWidth, stats[3]);

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

/**
    This is the `showPileUps` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.showPileUps;

import dentist.common.binio : PileUpDb;
import std.algorithm : max;
import std.file : getSize;
import std.math : log10, lrint, FloatingPointControl;
import std.stdio : writefln, writeln;
import std.typecons : tuple;

/// Execute the `showPileUps` command with `options`.
void execute(Options)(in Options options)
{
    size_t totalDbSize = options.pileUpsFile.getSize();
    auto pileUpDb = PileUpDb.parse(options.pileUpsFile);
    pileUpDb.releaseDb();

    auto stats = tuple(
        totalDbSize,
        pileUpDb.pileUps.length,
        pileUpDb.readAlignments.length,
        pileUpDb.seededAlignments.length,
        pileUpDb.localAlignments.length,
        pileUpDb.tracePoints.length,
    );

    FloatingPointControl fpCtrl;
    fpCtrl.rounding = FloatingPointControl.roundUp;
    auto numWidth = lrint(log10(max(stats.expand)));

    writefln!"totalDbSize:         %*d"(numWidth, stats[0]);
    writefln!"numPileUps:          %*d"(numWidth, stats[1]);
    writefln!"numReadAlignments:   %*d"(numWidth, stats[2]);
    writefln!"numSeededAlignments: %*d"(numWidth, stats[3]);
    writefln!"numLocalAlignments:  %*d"(numWidth, stats[4]);
    writefln!"numTracePoints:      %*d"(numWidth, stats[5]);
}

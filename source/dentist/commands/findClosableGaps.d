/**
    This is the `findClosableGaps` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.findClosableGaps;

import dentist.commandline : OptionsFor, TestingCommand;
import dentist.common :
    isTesting,
    ReferenceInterval,
    ReferenceRegion;
import dentist.common.alignments :
    AlignmentChain,
    coord_t,
    id_t;
import dentist.dazzler :
    ContigSegment,
    GapSegment,
    getScaffoldStructure,
    readMask,
    ScaffoldSegment;
import dentist.util.algorithm : cmpLexicographically;
import dentist.util.log;
import std.algorithm :
    copy,
    filter,
    map,
    sort,
    swap;
import std.array : array;
import std.format : format, formattedRead;
import std.range : assumeSorted, enumerate;
import std.stdio : File, writeln;
import std.string : capitalize;
import std.typecons : PhobosFlag = Flag, tuple;
import vibe.data.json : toJson = serializeToJson, serializeToJsonString;


static if (isTesting):

alias Options = OptionsFor!(TestingCommand.findClosableGaps);

/// Execute the `findClosableGaps` command with `options`.
void execute(in Options options)
{
    auto finder = ClosableGapsFinder(options);

    finder.run();
}

struct ClosableGapsFinder
{
    const(Options) options;
    const(ScaffoldSegment)[] trueAssemblyScaffoldStructure;
    ReferenceRegion mappedRegionsMask;
    TrueAlignment[] trueAlignments;
    ClosableGap[] closableGaps;

    void run()
    {
        init();
        findClosableGaps();

        writeln(serializeToJsonString(closableGaps));
    }

    void init()
    {
        trueAssemblyScaffoldStructure = getScaffoldStructure(options.trueAssemblyDb).array;
        mappedRegionsMask = ReferenceRegion(readMask!ReferenceInterval(
            options.trueAssemblyDb,
            options.mappedRegionsMask,
            options.workdir,
        ));
        trueAlignments = readTrueAlignments(options.readsMap);

        logJsonDebug(
            "trueAssemblyScaffoldStructure", trueAssemblyScaffoldStructure
                .map!(scaffoldPart => scaffoldPart.peek!GapSegment !is null
                    ? scaffoldPart.get!GapSegment.toJson
                    : scaffoldPart.get!ContigSegment.toJson)
                .array
                .toJson,
            "mappedRegionsMask", mappedRegionsMask.intervals.toJson,
            "trueAlignments", trueAlignments.toJson,
        );
    }

    void findClosableGaps()
    {
        alias isContigPart = contigPart => contigPart.peek!ContigSegment !is null;
        alias mkContigPart = contigPart => contigPart.get!ContigSegment;

        foreach (contigPart; trueAssemblyScaffoldStructure.filter!isContigPart.map!mkContigPart)
        {
            auto scaffoldRegion = ReferenceRegion(ReferenceInterval(
                contigPart.globalContigId,
                0,
                contigPart.length,
            ));
            auto scaffoldGaps = scaffoldRegion - mappedRegionsMask;
            closableGaps.length += scaffoldGaps.intervals.length;

            foreach (i, gap; scaffoldGaps.intervals)
            {
                closableGaps[i].fromContig = cast(id_t) i;
                closableGaps[i].toContig = cast(id_t) (i + 1);
                closableGaps[i].gapSize = cast(coord_t) gap.size;

                foreach (j, read; trueAlignments)
                    if (
                        read.begin < gap.begin &&
                        (gap.begin - read.begin) >= options.minAnchorLength &&
                        gap.end < read.end &&
                        (read.end - gap.end) >= options.minAnchorLength
                    )
                        closableGaps[i].spanningReads ~= cast(id_t) (j + 1);
            }
        }

        auto bufferRest = closableGaps
            .filter!(closableGap => closableGap.spanningReads.length >= options.minSpanningReads)
            .copy(closableGaps);
        closableGaps = closableGaps[0 .. $ - bufferRest.length];
    }
}

struct ClosableGap
{
    id_t fromContig;
    id_t toContig;
    coord_t gapSize;
    id_t[] spanningReads;
}

struct TrueAlignment
{
    static alias Flag = AlignmentChain.Flag;
    static alias Flags = AlignmentChain.Flags;

    id_t scaffoldId;
    coord_t begin;
    coord_t end;
    Flags flags;

    static foreach(flagName; __traits(allMembers, Flag))
    {
        mixin(format!(q"<
            static alias %1$s = PhobosFlag!"%2$s";

            @property PhobosFlag!"%2$s" %2$s() pure const nothrow @trusted
            {
                return cast(PhobosFlag!"%2$s") flags.%2$s;
            }

            @property void %2$s(PhobosFlag!"%2$s" %2$s) pure nothrow
            {
                flags.%2$s = %2$s;
            }
        >")(flagName.capitalize, flagName));
    }

    int opCmp(in TrueAlignment other) const pure nothrow
    {
        return cmpLexicographically!(
            TrueAlignment,
            "a.scaffoldId",
            "a.begin",
            "a.end",
        )(this, other);
    }
}

TrueAlignment[] readTrueAlignments(in string readsMap)
{
    auto trueAlignments = File(readsMap)
        .byLine
        .map!parseTrueAlignment
        .array;
    trueAlignments.sort;

    return trueAlignments;
}

// Not meant for public usage.
TrueAlignment parseTrueAlignment(S)(S line)
{
    TrueAlignment trueAlignment;

    line.formattedRead!" %d %d %d"(
        trueAlignment.scaffoldId,
        trueAlignment.begin,
        trueAlignment.end,
    );

    if (trueAlignment.begin > trueAlignment.end)
    {
        swap(trueAlignment.begin, trueAlignment.end);
        trueAlignment.flags |= TrueAlignment.Flag.complement;
    }

    return trueAlignment;
}

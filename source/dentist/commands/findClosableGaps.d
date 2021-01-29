/**
    This is the `findClosableGaps` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.findClosableGaps;

import dentist.commandline : OptionsFor;
import dentist.common :
    isTesting,
    ReferenceInterval,
    ReferenceRegion;
import dentist.common.alignments :
    AlignmentChain,
    AlignmentFlag = Flag,
    AlignmentFlags = Flags,
    coord_t,
    id_t;
import dentist.common.commands : TestingCommand;
import dentist.dazzler :
    DBdumpOptions,
    DbRecord,
    getDbRecords,
    readMask;
import dentist.util.algorithm :
    cmpLexicographically,
    filterInPlace;
import dentist.util.log;
import std.algorithm :
    copy,
    count,
    filter,
    map,
    sort,
    swap;
import std.array : array;
import std.format : format, formattedRead;
import std.range :
    assumeSorted,
    enumerate,
    slide;
import std.stdio : File, writeln;
import std.string : capitalize;
import std.typecons :
    No,
    PhobosFlag = Flag,
    tuple;
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
    /// Contig locations of the base assembly (trueAssembly)
    DbRecord[] baseContigs;
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
        baseContigs = getDbRecords(options.trueAssemblyDb, [
            DBdumpOptions.readNumber,
            DBdumpOptions.originalHeader,
        ]).array;

        auto mappedIntervals = readMask!ReferenceInterval(
            options.trueAssemblyDb,
            options.mappedRegionsMask,
        );
        mappedIntervals.filterInPlace!(mapped => mapped.size >= options.contigCutoff);
        mappedRegionsMask = ReferenceRegion(mappedIntervals);
        trueAlignments = readTrueAlignments(options.readsMap);

        logJsonDebug(
            "baseContigs", baseContigs.toJson,
            "mappedRegionsMask", mappedRegionsMask.intervals.toJson,
            "trueAlignments", trueAlignments.toJson,
        );
    }

    void findClosableGaps()
    {
        closableGaps = new ClosableGap[mappedRegionsMask.intervals.length - 1];
        auto mappedRegionsPairs = mappedRegionsMask
            .intervals
            .slide!(No.withPartial)(2)
            .enumerate;

        foreach (i, mappedPair; mappedRegionsPairs)
        {
            if (mappedPair[0].contigId != mappedPair[1].contigId)
                // skip if pair is not on the same contig, i.e. not a gap
                continue;

            auto gap = ReferenceInterval(
                mappedPair[0].contigId,
                mappedPair[0].end,
                mappedPair[1].begin,
            );
            auto baseContig = baseContigs[gap.contigId - 1];
            assert(baseContig.contigId == gap.contigId);
            // contigIdx starts at 0 on each scaffold -> count zeros
            // but subtract one because scaffolds are zero-indexed
            auto scaffoldId = cast(id_t) baseContigs[0 .. baseContig.contigId]
                .count!"a.location.contigIdx == 0" - 1;

            auto trueAlignments = this.trueAlignments
                .assumeSorted!"a.scaffoldId < b.scaffoldId"
                // find all TrueAlignments for baseContig
                .equalRange(TrueAlignment(scaffoldId))
                .filter!(read => baseContig.location.begin <= read.begin && read.end <= baseContig.location.end)
                // adjust to contig-coordinates
                .map!(read => TrueAlignment(
                    read.scaffoldId,
                    cast(coord_t) (read.begin - baseContig.location.begin),
                    cast(coord_t) (read.end - baseContig.location.begin),
                    read.flags,
                    read.readId,
                ))
                .array;

            closableGaps[i].fromContig = cast(id_t) (i + 1);
            closableGaps[i].toContig = cast(id_t) (i + 2);
            closableGaps[i].gapSize = cast(coord_t) gap.size;
            closableGaps[i].mappedInterval = gap;

            foreach (read; trueAlignments)
                if (
                    read.begin + options.minAnchorLength <= gap.begin &&
                    gap.end + options.minAnchorLength <= read.end
                )
                    // read spans gap including a context of
                    // minAnchorLength on either side
                    closableGaps[i].spanningReads ~= read.readId;

            closableGaps[i].spanningReads.sort;
        }

        closableGaps.filterInPlace!(
            closableGap => closableGap.spanningReads.length >= options.minSpanningReads
        );
    }
}

struct ClosableGap
{
    id_t fromContig;
    id_t toContig;
    coord_t gapSize;
    ReferenceInterval mappedInterval;
    id_t[] spanningReads;
}

struct TrueAlignment
{
    id_t scaffoldId;
    coord_t begin;
    coord_t end;
    AlignmentFlags flags;
    id_t readId;

    static foreach(flagName; __traits(allMembers, AlignmentFlag))
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
            "a.readId",
        )(this, other);
    }
}

TrueAlignment[] readTrueAlignments(in string readsMap)
{
    auto trueAlignments = File(readsMap)
        .byLine
        .enumerate(1)
        .map!parseTrueAlignment
        .array;
    trueAlignments.sort;

    return trueAlignments;
}

// Not meant for public usage.
TrueAlignment parseTrueAlignment(EnumLine)(EnumLine enumLine)
{
    auto line = enumLine.value;
    TrueAlignment trueAlignment;
    trueAlignment.readId = enumLine.index;

    line.formattedRead!" %d %d %d"(
        trueAlignment.scaffoldId,
        trueAlignment.begin,
        trueAlignment.end,
    );

    if (trueAlignment.begin > trueAlignment.end)
    {
        swap(trueAlignment.begin, trueAlignment.end);
        trueAlignment.flags |= AlignmentFlag.complement;
    }

    return trueAlignment;
}

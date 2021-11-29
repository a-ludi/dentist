/**
    This is the `findClosableGaps` command of DENTIST.

    Command_Summary:

    ---
    Find which gaps are closable, i.e. the true alignment of the reads
    provides a sufficient number spanning reads including some amount of
    anchor sequence.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.findClosableGaps;

package(dentist) enum summary = "
    Find which gaps are closable, i.e. the true alignment of the reads
    provides a sufficient number spanning reads including some amount of
    anchor sequence.
";

import dentist.commandline : OptionsFor;
import dentist.common :
    isTesting,
    ReferenceInterval;
import dentist.common.alignments :
    AlignmentChain,
    coord_t,
    id_t;
import dentist.common.commands : DentistCommand;
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
    tuple;
import vibe.data.json : toJson = serializeToJson, serializeToJsonString;


static if (isTesting):

/// Options for the `findClosableGaps` command.
alias Options = OptionsFor!(DentistCommand.findClosableGaps);

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
    ReferenceInterval[] mappedContigs;
    ReadSample[] readSamples;
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

        mappedContigs = readMask!ReferenceInterval(
            options.trueAssemblyDb,
            options.mappedRegionsMask,
        );
        mappedContigs.filterInPlace!(mapped => mapped.size >= options.contigCutoff);
        readSamples = readReadSamples(options.readsMap);

        logJsonDebug(
            "baseContigs", baseContigs.toJson,
            "mappedContigs", mappedContigs.toJson,
            "readSamples", readSamples.toJson,
        );
    }

    void findClosableGaps()
    {
        closableGaps = new ClosableGap[mappedContigs.length - 1];
        auto mappedContigPairs = mappedContigs
            .slide!(No.withPartial)(2)
            .enumerate;

        foreach (i, contigPair; mappedContigPairs)
        {
            if (contigPair[0].contigId != contigPair[1].contigId)
                // skip if pair is not on the same contig, i.e. not a gap
                continue;

            auto gap = ReferenceInterval(
                contigPair[0].contigId,
                contigPair[0].end,
                contigPair[1].begin,
            );
            auto baseContig = baseContigs[gap.contigId - 1];
            assert(baseContig.contigId == gap.contigId);
            // contigIdx starts at 0 on each scaffold -> count zeros
            // but subtract one because scaffolds are zero-indexed
            auto scaffoldId = cast(id_t) baseContigs[0 .. baseContig.contigId]
                .count!"a.location.contigIdx == 0" - 1;

            auto readSamples = this.readSamples
                .assumeSorted!"a.scaffoldId < b.scaffoldId"
                // find all ReadSamples for baseContig
                .equalRange(ReadSample(scaffoldId))
                .filter!(read => baseContig.location.begin <= read.begin && read.end <= baseContig.location.end)
                // adjust to contig-coordinates
                .map!(read => ReadSample(
                    read.scaffoldId,
                    cast(coord_t) (read.begin - baseContig.location.begin),
                    cast(coord_t) (read.end - baseContig.location.begin),
                    read.complement,
                    read.readId,
                ))
                .array;

            closableGaps[i].fromContig = cast(id_t) (i + 1);
            closableGaps[i].toContig = cast(id_t) (i + 2);
            closableGaps[i].gapSize = cast(coord_t) gap.size;
            closableGaps[i].mappedInterval = gap;

            foreach (read; readSamples)
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

struct ReadSample
{
    id_t scaffoldId;
    coord_t begin;
    coord_t end;
    bool complement;
    id_t readId;

    int opCmp(in ReadSample other) const pure nothrow
    {
        return cmpLexicographically!(
            ReadSample,
            "a.scaffoldId",
            "a.begin",
            "a.end",
            "a.readId",
        )(this, other);
    }
}

ReadSample[] readReadSamples(in string readsMap)
{
    auto readSamples = File(readsMap)
        .byLine
        .enumerate(1)
        .map!parseReadSample
        .array;
    readSamples.sort;

    return readSamples;
}

// Not meant for public usage.
ReadSample parseReadSample(EnumLine)(EnumLine enumLine)
{
    auto line = enumLine.value;
    ReadSample readSample;
    readSample.readId = enumLine.index;

    line.formattedRead!" %d %d %d"(
        readSample.scaffoldId,
        readSample.begin,
        readSample.end,
    );

    if (readSample.begin > readSample.end)
    {
        swap(readSample.begin, readSample.end);
        readSample.complement = true;
    }

    return readSample;
}

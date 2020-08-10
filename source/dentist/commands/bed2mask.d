/**
    This is the `bed2mask` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.bed2mask;

import dentist.commandline : OptionsFor;
import dentist.common : ReferenceInterval;
import dentist.common.alignments :
    coord_t,
    id_t;
import dentist.common.commands : DentistCommand;
import dentist.dazzler :
    ContigSegment,
    getScaffoldStructure,
    writeMask;
import dentist.util.log;
import std.algorithm :
    chunkBy,
    filter,
    joiner,
    min,
    map,
    max;
import std.array :
    array,
    split;
import std.conv : to;
import std.range : tee;
import std.range.primitives;
import std.stdio : File;
import vibe.data.json : toJson = serializeToJson;


/// Options for the `bed2mask` command.
alias Options = OptionsFor!(DentistCommand.bed2mask);


/// Execute the `bed2mask` command with `options`.
void execute(Options)(in Options options)
{
    auto contigsByScaffold = getContigsByScaffold(options.refDb);
    auto bedFile = options.openBedFile;
    auto mask = readBedFile(bedFile, contigsByScaffold);

    writeDazzlerMask(
        options.refDb,
        options.outMask,
        mask,
    );
}


ContigSegment[][string] getContigsByScaffold(string refDb)
{
    ContigSegment[][string] contigsByScaffold;

    auto contigChunks = getScaffoldStructure(refDb)
            .filter!(part => part.peek!ContigSegment !is null)
            .map!(contigPart => contigPart.get!ContigSegment)
            .chunkBy!"a.scaffoldId == b.scaffoldId";

    foreach (contigChunk; contigChunks)
        contigsByScaffold[contigChunk.front.header[1 .. $]] = contigChunk.array;

    return contigsByScaffold;
}


ReferenceInterval[] readBedFile(File bedFile, ContigSegment[][string] contigsByScaffold)
{
    mixin(traceExecution);

    return bedFile
        .byLine
        .filter!(line => line.length > 0)
        .map!((line) {
            auto fields = line.split('\t');
            scope scaffoldName = cast(string) fields[0];
            auto begin = fields[1].to!coord_t;
            auto end = fields[2].to!coord_t;
            auto affectedContigs = contigsByScaffold.getOverlappingContigs(
                scaffoldName,
                begin,
                end,
            );

            return affectedContigs.map!(affectedContig => ReferenceInterval(
                affectedContig.globalContigId,
                max(begin, affectedContig.begin) - affectedContig.begin,
                min(end, affectedContig.end) - affectedContig.begin,
            )).tee!(interval => assert(&interval));
        })
        .joiner
        .array;
}


ContigSegment[] getOverlappingContigs(
    scope ContigSegment[][string] contigsByScaffold,
    scope string scaffoldName,
    coord_t begin,
    coord_t end,
)
{
    auto overlappingContigs = contigsByScaffold[scaffoldName];

    while (overlappingContigs.front.end <= begin)
        overlappingContigs.popFront();
    while (overlappingContigs.back.begin >= end)
        overlappingContigs.popBack();

    return overlappingContigs;
}


void writeDazzlerMask(string refDb, string maskName, ReferenceInterval[] mask)
{
    mixin(traceExecution);

    writeMask(refDb, maskName, mask);
}

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
    dazzExtra,
    getScaffoldStructure,
    writeDazzExtra,
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
import std.conv :
    ConvException,
    to;
import std.exception : enforce;
import std.format : format;
import std.range :
    chain,
    enumerate,
    only,
    tee;
import std.range.primitives;
import std.regex :
    ctRegex,
    matchFirst;
import std.stdio : File;
import std.typecons : Tuple;
import vibe.data.json : toJson = serializeToJson;


/// Options for the `bed2mask` command.
alias Options = OptionsFor!(DentistCommand.bed2mask);


/// Execute the `bed2mask` command with `options`.
void execute(Options)(in Options options)
{
    auto contigsByScaffold = getContigsByScaffold(options.refDb);
    auto bedFile = options.openBedFile;
    auto augmentedMask = readBedFile(
        bedFile,
        options.bedFileName,
        contigsByScaffold,
        cast(bool) options.dataComments,
    );

    writeDazzlerMask(
        options.refDb,
        options.outMask,
        augmentedMask,
        cast(bool) options.dataComments,
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


alias AugmentedReferenceInterval = Tuple!(
    ReferenceInterval, "interval",
    id_t[], "contigIds",
    id_t[], "readIds",
);


AugmentedReferenceInterval[] readBedFile(
    File bedFile,
    string bedFileName,
    ContigSegment[][string] contigsByScaffold,
    bool parseDataComments,
)
{
    mixin(traceExecution);

    return bedFile
        .byLine
        .enumerate(1)
        .filter!(enumLine => enumLine.value.length > 0)
        .map!((enumLine) {
            auto lineNumber = enumLine.index;
            auto line = enumLine.value;
            auto fields = line.split('\t');
            scope scaffoldName = cast(string) fields[0];
            auto begin = fields[1].to!coord_t;
            auto end = fields[2].to!coord_t;
            auto affectedContigs = contigsByScaffold.getOverlappingContigs(
                scaffoldName,
                begin,
                end,
            );

            id_t[] contigIds;
            id_t[] readIds;

            if (parseDataComments && fields.length >= 4)
                parseDataComment(fields[3], bedFileName, lineNumber, contigIds, readIds);

            return affectedContigs.map!(affectedContig => AugmentedReferenceInterval(
                ReferenceInterval(
                    affectedContig.globalContigId,
                    max(begin, affectedContig.begin) - affectedContig.begin,
                    min(end, affectedContig.end) - affectedContig.begin,
                ),
                contigIds,
                readIds,
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


void parseDataComment(
    const char[] comment,
    string filename,
    size_t lineNumber,
    ref id_t[] contigIds,
    ref id_t[] readIds,
)
{
    string[] errors;

    void parseFields(T)(ref T[] dest, const char[][] fields)
    {
        try
        {
            dest = fields[1 .. $].map!(to!T).array;
        }
        catch (ConvException e)
        {
            dest = [];
            errors ~= "could not parse " ~ T.stringof ~ " in " ~ fields[0].idup;
        }
    }

    auto parts = comment.split('|');

    foreach (part; parts)
    {
        auto partFields = part.split('-');
        auto partName = partFields[0];

        switch (partName)
        {
            case Options.contigsExtraName:
                if (partFields.length == 3)
                {
                    parseFields(contigIds, partFields);
                }
                else
                {
                    errors ~= Options.contigsExtraName ~ " must have exactly two entries";

                    break;
                }
                break;
            case Options.readsExtraName:
                parseFields(readIds, partFields);
                break;
            default:
                break;
        }
    }

    enforce(
        errors.length == 0,
        format!(
            "ill-formatted data comment in BED file %s:%d:\n%-(  - %s\n%)"
        )(filename, lineNumber, errors),
    );
}


void writeDazzlerMask(
    string refDb,
    string maskName,
    AugmentedReferenceInterval[] augmentedMask,
    bool hasDataComments,
)
{
    mixin(traceExecution);

    writeMask(refDb, maskName, augmentedMask.map!"a.interval");

    if (hasDataComments)
    {
        auto contigsExtra = dazzExtra(Options.contigsExtraName, augmentedMask
            .map!"a.contigIds"
            .joiner
            .map!"cast(long) a"
            .array);
        writeDazzExtra(refDb, maskName, contigsExtra);

        auto readsExtra = dazzExtra(Options.readsExtraName, augmentedMask
            .map!(a => chain(
                only(a.readIds.length.to!id_t),
                a.readIds,
            ))
            .joiner
            .map!"cast(long) a"
            .array);
        writeDazzExtra(refDb, maskName, readsExtra);
    }
}

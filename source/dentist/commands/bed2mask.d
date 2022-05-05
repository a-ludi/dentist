/**
    The `bed2mask` command creates a Dazzler mask from a BED file.

    Command_Summary:

    ---
    Convert a BED file to a Dazzler mask. Implements 'data comments' -- a
    special feature for DENTIST.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.bed2mask;

package(dentist) enum summary = "
    Convert a BED file to a Dazzler mask. Implements 'data comments' -- a
    special feature for DENTIST.
";

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
void execute(in Options options)
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


/// Returns an associative array containing an array of `ContigSegment`s for
/// each scaffold identified by its FASTA header (without leading `>`).
ContigSegment[][string] getContigsByScaffold(string refDb)
{
    ContigSegment[][string] contigsByScaffold;

    auto contigChunks = getScaffoldStructure(refDb)
            .filter!(part => part.peek!ContigSegment !is null)
            .map!(contigPart => contigPart.get!ContigSegment)
            .chunkBy!"a.scaffoldId == b.scaffoldId";

    foreach (contigChunk; contigChunks)
        contigsByScaffold[contigChunk.front.header.split("\t")[0][1 .. $]] = contigChunk.array;

    return contigsByScaffold;
}


/// `ReferenceInterval` with optional supplementary data.
///
/// See_also: `parseDataComment`
alias AugmentedReferenceInterval = Tuple!(
    ReferenceInterval, "interval",
    id_t[], "contigIds",
    id_t[], "readIds",
);


/// Main procedure that converts BED entries into an array of
/// `AugmentedReferenceInterval`s. The returned intervals have Dazzler
/// coordinates and contain data from the comment column if
/// `parseDataComments`.
///
/// Params:
///     bedFile = `File` that contains the BED data.
///     bedFileName = file name that is displayed in error messages
///     contigsByScaffold = index of Dazzler DB (see `getContigsByScaffold`)
///     parseDataComments = expected 4th column in BED to be a data comment if true.
/// See_also: `getContigsByScaffold`, `parseDataComment`
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
            // Convert BED file line-by-line ...
            auto lineNumber = enumLine.index;
            auto line = enumLine.value;
            // Split line into columns
            auto fields = line.split('\t');
            scope scaffoldName = cast(string) fields[0];
            auto begin = fields[1].to!coord_t;
            auto end = fields[2].to!coord_t;
            // Find contigs that overlap with the interval `[begin, end]` in
            // scaffold coordinates
            auto affectedContigs = contigsByScaffold.getOverlappingContigs(
                scaffoldName,
                begin,
                end,
            );

            id_t[] contigIds;
            id_t[] readIds;

            if (parseDataComments && fields.length >= 4)
                // Parse contigIds and readIds from comment column
                parseDataComment(fields[3], bedFileName, lineNumber, contigIds, readIds);

            // Generate list of `AugmentedReferenceInterval`s
            return affectedContigs.map!(affectedContig => AugmentedReferenceInterval(
                // Interval with scaffold-global contig ID and scaffold coordinates
                ReferenceInterval(
                    affectedContig.globalContigId,
                    // Crop contig interval to `[begin, end]`
                    max(begin, affectedContig.begin) - affectedContig.begin,
                    min(end, affectedContig.end) - affectedContig.begin,
                ),
                contigIds,
                readIds,
            )).tee!(interval => assert(&interval));
        })
        // Concatenate lists from each line
        .joiner
        .array;
}


/// Find contigs that overlap with the interval `[begin, end]` on `scaffold`.
///
/// See_also: `getContigsByScaffold`
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


/// Parse data from `comment` into `contigIds` and `readIds`, respectively.
/// If multiple parts of the same type are given later parts overwrite
/// previous parts.
///
/// Grammar:
/// ---
/// <data-comment>   ::== <part> | <part> "|" <data-comment>
/// <part>           ::== <contigs-part> | <reads-part>
/// <contigs-part>   ::== "contigs-" <id> "-" <id>
/// <reads-part>     ::== "reads" <id-list>
/// <id-list>        ::== "-" <id> | "-" <id> <id-list>
/// <id>             ::== <positive-digit> | <positive-digit> <digits>
/// <digits>         ::== <digit> | <digit> <digits>
/// <positive-digit> ::== "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9"
/// <digit>          ::== "0" | <positive-digit>
/// ---
///
/// Params:
///     comment = comment column (4th) from BED file
///     filename = filename to report in errors
///     lineNumber = line number to report in errors
///     contigIds = destination for contig IDs
///     readIds = destination for read IDs
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


/// Create Dazzler mask of intervals in `augmentedMask` adding extra tracks
/// if `hasDataComments` is true.
///
/// Overwrites existing mask without asking.
///
/// Params:
///     refDb = Dazzler DB for which the track is created
///     maskName = name of the mask
///     augmentedMask = mask data
///     hasDataComments = create mask track extras for contig and read IDs
///
/// See_also: `dentist.dazzler.writeMask`, `dentist.dazzler.dazzExtra`,
///     `dentist.dazzler.writeDazzExtra`
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

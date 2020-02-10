/**
    Everything to handle insertions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.insertions;

import dentist.common :
    dentistEnforce,
    ReadInterval,
    ReferenceInterval;
import dentist.common.alignments :
    AlignmentLocationSeed,
    coord_t,
    id_t,
    SeededAlignment;
import dentist.common.binio : CompressedSequence;
import dentist.common.scaffold :
    ContigNode,
    ContigPart,
    isDefault,
    isExtension,
    isGap,
    Scaffold;
import dentist.util.log;
import dentist.util.math : add;
import std.algorithm :
    among,
    canFind,
    filter;
import std.array : array;
import std.typecons : tuple;
import vibe.data.json : toJson = serializeToJson;


/// This characterizes an insertion.
struct InsertionInfo
{
    CompressedSequence sequence;
    size_t contigLength;
    SeededAlignment[] overlaps;
    id_t[] readIds;
}

/// This is used to collect all sub-sequences of the output.
alias OutputScaffold = Scaffold!InsertionInfo;

/// This characterizes an insertion.
alias Insertion = OutputScaffold.Edge;

/// Returns true iff a sequence of `n`s should be written.
bool isOutputGap(in Insertion insertion)
{
    return !insertion.hasSequence && !insertion.isDefault;
}

bool isValidInsertion(in Insertion insertion)
{
    return
        insertion.isDefault ^
        insertion.isOutputGap ^
        (insertion.isGap && insertion.hasSequence) ^
        (insertion.isExtension && insertion.hasSequence);
}

bool hasSequence(in Insertion insertion)
{
    return insertion.payload.sequence.length > 0;
}

coord_t getCroppingPosition(string contig)(in SeededAlignment overlap)
    if (contig.among("contigA", "contigB"))
{
    final switch (overlap.seed)
    {
        case AlignmentLocationSeed.front:
            return mixin("overlap.first." ~ contig ~ ".begin");
        case AlignmentLocationSeed.back:
            return mixin("overlap.last." ~ contig ~ ".end");
    }
}

auto getInfoForExistingContig(in ContigNode begin, in Insertion insertion, in bool globalComplement)
{
    auto overlaps = insertion.payload.overlaps;
    auto contigId = cast(id_t) begin.contigId;
    auto contigLength = insertion.payload.contigLength;
    auto slice = ReferenceInterval(contigId, 0, contigLength);

    assert(contigLength > 0);

    switch (overlaps.length)
    {
    case 0:
        // Keep whole contig
        break;
    case 1:
        auto overlap = overlaps[0];

        final switch (overlap.seed)
        {
        case AlignmentLocationSeed.front:
            slice.begin = getCroppingPosition!"contigA"(overlap);
            break;
        case AlignmentLocationSeed.back:
            slice.end = getCroppingPosition!"contigA"(overlap);
            break;
        }
        break;
    case 2:
        assert(overlaps[0].contigA.id == overlaps[1].contigA.id);
        assert(overlaps[0].seed != overlaps[1].seed);

        if (overlaps[0].seed < overlaps[1].seed)
        {
            slice.begin = getCroppingPosition!"contigA"(overlaps[0]);
            slice.end = getCroppingPosition!"contigA"(overlaps[1]);
        }
        else
        {
            slice.begin = getCroppingPosition!"contigA"(overlaps[1]);
            slice.end = getCroppingPosition!"contigA"(overlaps[0]);
        }
        break;
    default:
        dentistEnforce(0, "too many splice sites");
    }

    assert(
        // walk contig in correct orientation
        globalComplement != (begin < insertion.target(begin)) ||
        // the contig is circular
        begin.contigId == insertion.target(begin).contigId
    );

    return tuple!(
        "contigId",
        "contigLength",
        "cropping",
        "length",
        "complement",
    )(
        contigId,
        contigLength,
        slice,
        slice.size,
        globalComplement,
    );
}

auto getInfoForGap(in Insertion insertion)
{
    return tuple!("length")(insertion.payload.contigLength);
}

auto getInfoForNewSequenceInsertion(
    in ContigNode begin,
    in Insertion insertion,
    in bool globalComplement,
)
{
    auto overlaps = insertion.payload.overlaps;
    auto firstOverlap = overlaps[0].contigA.id == begin.contigId ? overlaps[0] : overlaps[1];
    auto effectiveComplement = firstOverlap.flags.complement ^ globalComplement;

    assert(
        (insertion.isExtension && overlaps.length == 1) ^
        (insertion.isGap && overlaps.length == 2)
    );

    auto slice = ReadInterval(
        firstOverlap.contigB.id,
        0,
        firstOverlap.contigB.length,
    );
    switch (overlaps.length)
    {
    case 1:
        final switch (firstOverlap.seed)
        {
        case AlignmentLocationSeed.front:
            slice.end = getCroppingPosition!"contigB"(firstOverlap);
            break;
        case AlignmentLocationSeed.back:
            slice.begin = getCroppingPosition!"contigB"(firstOverlap);
            break;
        }
        break;
    case 2:
        if (overlaps[0].seed > overlaps[1].seed)
        {
            slice.begin = getCroppingPosition!"contigB"(overlaps[0]);
            slice.end = getCroppingPosition!"contigB"(overlaps[1]);
        }
        else
        {
            slice.begin = getCroppingPosition!"contigB"(overlaps[1]);
            slice.end = getCroppingPosition!"contigB"(overlaps[0]);
        }
        if (slice.end < slice.begin)
        {
            logJsonWarn(
                "info", "adjusting invalid interval",
                "cropping", [
                    "begin": slice.begin,
                    "end": slice.end,
                ].toJson,
                "sequence", insertion.payload.sequence.to!string,
            );
            slice.end = slice.begin;
        }
        break;
    default:
        assert(0);
    }

    return tuple!(
        "sequence",
        "cropping",
        "length",
        "complement",
    )(
        insertion.payload.sequence,
        slice,
        slice.size,
        effectiveComplement,
    );
}

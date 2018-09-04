/**
    Everything to handle insertions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.insertions;

import dentist.common : ReferencePoint;
import dentist.common.alignments : AlignmentChain, id_t;
import dentist.common.binio : CompressedSequence;
import dentist.common.scaffold :
    ContigNode,
    ContigPart,
    isDefault,
    isExtension,
    isGap,
    Scaffold;
import std.algorithm :
    canFind,
    filter,
    swap;
import std.array : array;
import std.typecons : Flag, tuple, Tuple;

/// Information about the point where the two sequences should be spliced.
alias SpliceSite = Tuple!(
    ReferencePoint, "croppingRefPosition",
    AlignmentChain.Flags, "flags",
);

/// This characterizes an insertion.
alias InsertionInfo = Tuple!(
    CompressedSequence, "sequence",
    size_t, "contigLength",
    SpliceSite[], "spliceSites",
);

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

Insertion concatenateSpliceSites(Insertion existingJoin, Insertion newJoin)
{
    if (existingJoin.payload.sequence.length == 0)
    {
        existingJoin.payload.sequence = newJoin.payload.sequence;
    }
    existingJoin.payload.spliceSites ~= newJoin.payload.spliceSites;

    return existingJoin;
}

/// Remove contig cropping where no new sequence is to be inserted.
OutputScaffold fixContigCropping(OutputScaffold scaffold)
{
    alias replace = OutputScaffold.ConflictStrategy.replace;
    auto contigJoins = scaffold.edges.filter!isDefault;

    foreach (contigJoin; contigJoins)
    {
        bool insertionUpdated;

        foreach (contigNode; [contigJoin.start, contigJoin.end])
        {
            auto shouldInsertNewSequence = scaffold
                .incidentEdges(contigNode)
                .canFind!(insertion => !insertion.isOutputGap && (insertion.isGap || insertion.isExtension));

            if (!shouldInsertNewSequence)
            {
                auto contigLength = contigJoin.payload.contigLength;
                auto newSpliceSites = contigJoin
                    .payload
                    .spliceSites
                    .filter!(spliceSite => contigNode.contigPart == ContigPart.begin
                        ? !(spliceSite.croppingRefPosition.value < contigLength / 2)
                        : !(spliceSite.croppingRefPosition.value >= contigLength / 2))
                    .array;
                if (newSpliceSites.length < contigJoin.payload.spliceSites.length)
                {
                    contigJoin.payload.spliceSites = newSpliceSites;
                    insertionUpdated = true;
                }
            }
        }

        if (insertionUpdated)
        {
            scaffold.add!replace(contigJoin);
        }
    }

    return scaffold;
}

auto getInfoForExistingContig(in ContigNode begin, in Insertion insertion, in bool globalComplement)
{
    auto spliceSites = insertion.payload.spliceSites;
    auto contigId = cast(id_t) begin.contigId;
    auto contigLength = insertion.payload.contigLength;
    size_t spliceStart;
    size_t spliceEnd;

    assert(contigLength > 0);

    switch (spliceSites.length)
    {
    case 0:
        spliceStart = 0;
        spliceEnd = contigLength;
        break;
    case 1:
        auto splicePosition = spliceSites[0].croppingRefPosition.value;

        if (splicePosition < contigLength / 2)
        {
            spliceStart = splicePosition;
            spliceEnd = contigLength;
        }
        else
        {
            spliceStart = 0;
            spliceEnd = splicePosition;
        }
        break;
    case 2:
        assert(spliceSites.length == 2);
        assert(spliceSites[0].croppingRefPosition.contigId
                == spliceSites[1].croppingRefPosition.contigId);

        spliceStart = spliceSites[0].croppingRefPosition.value;
        spliceEnd = spliceSites[1].croppingRefPosition.value;

        if (spliceEnd < spliceStart)
        {
            swap(spliceStart, spliceEnd);
        }
        break;
    default:
        assert(0, "too many spliceSites");
    }

    assert(globalComplement != (begin < insertion.target(begin)));

    return tuple!(
        "contigId",
        "contigLength",
        "spliceStart",
        "spliceEnd",
        "length",
        "complement",
    )(
        contigId,
        contigLength,
        spliceStart,
        spliceEnd,
        spliceEnd - spliceStart,
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
    auto spliceSites = insertion.payload.spliceSites;
    auto effectiveComplement = spliceSites[0].flags.complement ^ globalComplement;

    assert(
        (insertion.isExtension && spliceSites.length == 1) ^
        (insertion.isGap && spliceSites.length == 2)
    );

    return tuple!(
        "sequence",
        "length",
        "complement",
    )(
        insertion.payload.sequence,
        insertion.payload.sequence.length,
        effectiveComplement,
    );
}

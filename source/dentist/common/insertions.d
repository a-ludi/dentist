/**
    Everything to handle insertions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.insertions;

import dentist.common : ReferencePoint;
import dentist.common.alignments : AlignmentChain;
import dentist.common.binio : CompressedSequence;
import dentist.common.scaffold :
    ContigPart,
    isDefault,
    isExtension,
    isGap,
    Scaffold;
import std.algorithm :
    canFind,
    filter;
import std.array : array;
import std.typecons : Flag, Tuple;

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

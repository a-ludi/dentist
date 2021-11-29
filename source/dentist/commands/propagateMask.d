/**
    This is the `propagateMask` command of DENTIST.

    Command_Summary:

    ---
    Propagate masked regions through the provided alignment. That means the
    mask is first transferred to the B-contigs/reads according to the given
    alignments.

    The default workflow is to first propagate from the reference assembly to
    the reads and then back again to the reference. Propagating, once again,
    to the reads will produce a complete repeat mask on the reads.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.propagateMask;

package(dentist) enum summary = "
    Propagate masked regions through the provided alignment. That means the
    mask is first transferred to the B-contigs/reads according to the given
    alignments.

    The default workflow is to first propagate from the reference assembly to
    the reads and then back again to the reference. Propagating, once again,
    to the reads will produce a complete repeat mask on the reads.
";

import dentist.commandline : OptionsFor;
import dentist.common :
    ReferenceInterval,
    ReferenceRegion;
import dentist.common.alignments :
    coord_t,
    id_t,
    FlatLocalAlignment;
import dentist.common.commands : DentistCommand;
import dentist.util.algorithm : sliceBy;
import dentist.util.log;
import dentist.util.math :
    ceildiv,
    RoundingMode;
import dentist.dazzler :
    AlignmentHeader,
    BufferMode,
    getFlatLocalAlignments,
    getNumContigs,
    readMask,
    writeMask;
import std.algorithm :
    chunkBy,
    copy,
    map,
    max,
    min,
    sort,
    swap;
import std.array : array, appender;
import std.format : format;
import std.range :
    assumeSorted,
    enumerate,
    tee;
import std.range.primitives;
import std.typecons :
    tuple,
    Yes;
import vibe.data.json : toJson = serializeToJson;


/// Options for the `propagateMask` command.
alias Options = OptionsFor!(DentistCommand.propagateMask);


/// Alias for readbility.
alias QueryInterval = ReferenceInterval;
/// ditto
alias QueryRegion = ReferenceRegion;


/// Execute the `propagateMask` command with `options`.
void execute(in Options options)
{
    return new MaskPropagator(options).run();
}


class MaskPropagator
{

    const Options options;
    string destinationDb;
    id_t numSourceContigs;
    ReferenceRegion inputMask;
    QueryInterval[][] outputMaskByContigs;
    QueryRegion outputMask;


    this(const Options options)
    {
        this.options = options;
        this.destinationDb = options.hasReadsDb
            ? options.readsDb
            : options.refDb;
    }


    void run()
    {
        mixin(traceExecution);

        init();

        auto bufferRest = getLocalAlignmentsByContig()
            .map!(localAlignments => propagateMaskPerContig(inputMask.intervals, localAlignments))
            .copy(outputMaskByContigs);
        outputMaskByContigs.length -= bufferRest.length;

        mergeMasks();
        writeMask(destinationDb, options.repeatMask, outputMask.intervals);
    }


    void init()
    {
        numSourceContigs = getNumContigs(options.refDb);
        outputMaskByContigs = new QueryInterval[][numSourceContigs];
        readMasks();
    }


    void readMasks()
    {
        mixin(traceExecution);

        foreach (mask; options.repeatMasks)
            this.inputMask |= ReferenceRegion(readMask!ReferenceInterval(options.refDb, mask));
    }


    auto getLocalAlignmentsByContig()
    {
        auto localAlignments = getFlatLocalAlignments(
            options.refDb,
            options.readsDb,
            options.dbAlignmentFile,
            BufferMode.preallocated,
        ).array;

        return localAlignments.sliceBy!((a, b) => a.contigA.id == b.contigA.id);
    }


    static QueryInterval[] propagateMaskPerContig(
        const(ReferenceInterval)[] inputMask,
        FlatLocalAlignment[] localAlignments,
    )
    {
        mixin(traceExecution);

        auto contigId = localAlignments[0].contigA.id;
        const(ReferenceInterval)[] maskIntervals = inputMask
            .assumeSorted!"a.contigId < b.contigId"
            .equalRange(ReferenceInterval(contigId))
            .release;
        if (maskIntervals.length > 0)
            sort!"a.contigA.begin < b.contigA.begin"(localAlignments);

        size_t maskIdx;
        size_t laIdx;
        auto propagatedMaskAcc = appender!(QueryInterval[]);
        propagatedMaskAcc.reserve(localAlignments.length);

        while (maskIdx < maskIntervals.length && laIdx < localAlignments.length)
        {
            auto currentMaskInterval = maskIntervals[maskIdx];
            auto currentLA = localAlignments[laIdx];

            if (currentMaskInterval.end <= currentLA.contigA.begin)
            {
                ++maskIdx; // mask is lagging behind -> catch up!
            }
            else if (currentLA.contigA.end <= currentMaskInterval.begin)
            {
                ++laIdx; // local alignment is lagging behind -> catch up!
            }
            else
            {
                auto intersectionIntervals = getIntersectionIntervals(
                    currentLA,
                    maskIntervals[maskIdx .. $],
                );
                auto propagatedIntervals = propagateIntervals(currentLA, intersectionIntervals);
                propagatedMaskAcc ~= propagatedIntervals;
                ++laIdx;
            }
        }

        logJsonInfo(
            "info", format!"propagated mask of contig %d"(contigId),
            "refContigId", contigId,
            "numInputMaskIntervals", maskIntervals.length,
            "numLocalAlignment", localAlignments.length,
            "numPropagatedMaskIntervals", propagatedMaskAcc.data.length,
        );

        return propagatedMaskAcc.data;
    }


    /// Identify maskIntervals intersecting currentLA and return a range with
    /// first and last interval cropped to the extent of currentLA.
    static auto getIntersectionIntervals(
        const FlatLocalAlignment currentLA,
        const(ReferenceInterval)[] maskIntervals,
    )
    {
        auto intersectingIntervals = getIntersectingIntervals(currentLA, maskIntervals);

        const ReferenceInterval cropToIntervalOfLocalAlignment(E)(E enumInterval) pure nothrow @safe
            out(i; currentLA.contigA.begin <= i.begin && i.end <= currentLA.contigA.end)
        {
            if (0 < enumInterval.index && enumInterval.index < intersectingIntervals.length - 1)
                return enumInterval.value;

            auto interval = ReferenceInterval(
                enumInterval.value.contigId,
                enumInterval.value.begin,
                enumInterval.value.end,
            );

            if (enumInterval.index == 0)
                interval.begin = max(interval.begin, currentLA.contigA.begin);
            if (enumInterval.index == intersectingIntervals.length - 1)
                interval.end = min(interval.end, currentLA.contigA.end);

            return interval;
        }

        return intersectingIntervals
            .enumerate
            .map!(e => cropToIntervalOfLocalAlignment(e));
    }


    /// Identify all mask intervals intersecting currentLA.
    static const(ReferenceInterval)[] getIntersectingIntervals(
        const FlatLocalAlignment currentLA,
        const(ReferenceInterval)[] maskIntervals,
    )
        in (
            currentLA.contigA.begin < maskIntervals[0].end &&
            maskIntervals[0].begin < currentLA.contigA.end,
            "first mask interval and local alignment should intersect",
        )
    {
        auto intersectionEnd = ReferenceInterval(0, currentLA.contigA.end);
        auto sortedMaskIntervals = maskIntervals.assumeSorted!"a.begin < b.begin";
        auto intersectingIntervals = sortedMaskIntervals.lowerBound(intersectionEnd);

        return intersectingIntervals.release;
    }


    /// Propagate intersectionIntervals by means of currentLA and return the
    /// resulting query intervals.
    static QueryInterval[] propagateIntervals(R)(
        const FlatLocalAlignment currentLA,
        R intersectionIntervals,
    )
    {
        auto propagatedIntervals = new QueryInterval[intersectionIntervals.length];

        foreach (i, intersectingInterval; enumerate(intersectionIntervals))
        {
            assert(currentLA.contigA.begin <= intersectingInterval.begin && intersectingInterval.end <= currentLA.contigA.end);
            auto propagatedInterval = QueryInterval(
                currentLA.contigB.id,
                currentLA
                    .translateTracePoint(cast(coord_t) intersectingInterval.begin, RoundingMode.floor)
                    .contigB,
                currentLA
                    .translateTracePoint(cast(coord_t) intersectingInterval.end, RoundingMode.ceil)
                    .contigB,
            );
            assert(propagatedInterval.begin <= propagatedInterval.end);

            if (currentLA.flags.complement)
            {
                propagatedInterval.begin = currentLA.contigB.length - propagatedInterval.begin;
                propagatedInterval.end = currentLA.contigB.length - propagatedInterval.end;
                swap(propagatedInterval.begin, propagatedInterval.end);
            }

            propagatedIntervals[i] = propagatedInterval;
        }

        return propagatedIntervals;
    }

    void mergeMasks()
    {
        mixin(traceExecution);

        foreach (outputMaskByContig; outputMaskByContigs)
            outputMask |= QueryRegion(outputMaskByContig);
    }
}

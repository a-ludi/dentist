/**
    This is the `mask-repetitive-regions` command of DENTIST.

    Command_Summary:

    ---
    Mask regions that have a alignment coverage that is out of bounds.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.maskRepetitiveRegions;

package(dentist) enum summary = "Mask regions that have a alignment coverage that is out of bounds.";

import dentist.commandline : OptionsFor;
import dentist.common :
    ReferenceInterval,
    ReferenceRegion;
import dentist.common.alignments :
    AlignmentChain,
    AlignmentFlags = Flags,
    Contig,
    coord_t,
    id_t,
    Locus;
import dentist.common.commands : DentistCommand;
import dentist.dazzler :
    alignmentChainPacker,
    AlignmentHeader,
    BufferMode,
    DBdumpOptions,
    getDbRecords,
    getFlatLocalAlignments,
    LocalAlignmentReader,
    writeMask;
import dentist.util.log;
import std.algorithm :
    filter,
    joiner,
    map,
    predSwitch,
    sort,
    uniq;
import std.array : appender, array, uninitializedArray;
import std.conv : to;
import std.format : format;
import std.range :
    chain,
    only;
import std.range.primitives :
    empty,
    ElementType,
    isInputRange;
import std.typecons : Flag, No, Tuple, Yes;
import vibe.data.json : toJson = serializeToJson;


/// Options for the `mask-repetitive-regions` command.
alias Options = OptionsFor!(DentistCommand.maskRepetitiveRegions);


/// Execute the `mask-repetitive-regions` command with `options`.
void execute(in Options options)
{
    auto assessor = new RepeatMaskAssessor(options);

    assessor.run();
}

enum AlignmentType : ubyte
{
    self,
    reads,
}

/// This class comprises the `maskRepetitiveRegions` step of the DENTIST algorithm
class RepeatMaskAssessor
{
    protected const Options options;
    protected AlignmentType alignmentType;
    protected AlignmentHeader alignmentHeader;
    protected LocalAlignmentReader alignment;
    protected ReferenceRegion repetitiveRegions;
    protected ReferenceRegion repetitiveRegionsImproper;

    this(in ref Options options)
    {
        this.options = options;
    }

    void run()
    {
        mixin(traceExecution);

        readInputs();
        assessRepeatStructure();
        writeRepeatMask();
    }

    protected void readInputs()
    {
        mixin(traceExecution);

        alignmentType = options.readsDb is null
            ? AlignmentType.self
            : AlignmentType.reads;
        alignmentHeader = AlignmentHeader.inferFrom(options.dbAlignmentFile);

        final switch (alignmentType)
        {
            case AlignmentType.self:
                alignment = getFlatLocalAlignments(
                    options.refDb,
                    options.dbAlignmentFile,
                );
                break;
            case AlignmentType.reads:
                alignment = getFlatLocalAlignments(
                    options.refDb,
                    options.readsDb,
                    options.dbAlignmentFile,
                );
                break;
        }

        if (alignment.empty)
            logJsonWarn("info", "empty " ~ alignmentType.to!string ~ "-alignment");
    }

    void assessRepeatStructure()
    {
        mixin(traceExecution);

        id_t[2] coverageBounds = predSwitch(alignmentType,
            AlignmentType.self, options.coverageBoundsSelf,
            AlignmentType.reads, options.coverageBoundsReads,
        );
        auto repeatAssessor = new BadAlignmentCoverageAssessor(
            coverageBounds[0],
            coverageBounds[1]
        );

        repetitiveRegions = repeatAssessor(
            alignmentIntervals(),
            contigIntervals(),
        );

        logJsonDiagnostic(
            "alignmentType", alignmentType.to!string,
            "repetitiveRegions", shouldLog(LogLevel.debug_)
                ? repetitiveRegions.intervals.toJson
                : toJson(null),
            "numRepetitiveRegions", repetitiveRegions.intervals.length,
        );

        if (alignmentType == AlignmentType.reads)
        {
            auto improperRepeatAssessor = new BadAlignmentCoverageAssessor(
                options.improperCoverageBoundsReads[0],
                options.improperCoverageBoundsReads[1]
            );

            alignment.reset();
            repetitiveRegionsImproper = improperRepeatAssessor(
                alignmentIntervals(Yes.improperOnly),
                contigIntervals(),
            );

            logJsonDiagnostic(
                "alignmentType", alignmentType.to!string,
                "repetitiveRegionsImproper", shouldLog(LogLevel.debug_)
                    ? repetitiveRegionsImproper.intervals.toJson
                    : toJson(null),
                "numRepetitiveRegionsImproper", repetitiveRegionsImproper.intervals.length,
            );
        }
    }

    protected auto alignmentIntervals(Flag!"improperOnly" improperOnly = No.improperOnly)
    {
        alignment.reset();
        static AlignmentChain.LocalAlignment[] localAlignmentBuffer;

        if (localAlignmentBuffer.length == 0)
            localAlignmentBuffer = uninitializedArray!(AlignmentChain.LocalAlignment[])(
                alignmentHeader.maxLocalAlignments,
            );

        return alignment
            .alignmentChainPacker(BufferMode.overwrite, localAlignmentBuffer)
            .filter!(ac => !improperOnly || !ac.isProper(options.properAlignmentAllowance))
            .map!(ac => ReferenceInterval(
                ac.contigA.id,
                ac.first.contigA.begin,
                ac.last.contigA.end,
            ));
    }

    protected auto contigIntervals()
    {
        return getDbRecords(options.refDb, [DBdumpOptions.readNumber, DBdumpOptions.originalHeader])
            .map!(contig => ReferenceInterval(
                contig.contigId,
                0,
                contig.location.length,
            ));
    }

    protected void writeRepeatMask()
    {
        mixin(traceExecution);

        auto mergedMask = repetitiveRegions | repetitiveRegionsImproper;

        writeMask(
            options.refDb,
            options.repeatMask,
            mergedMask.intervals,
        );

        if (alignmentType == AlignmentType.reads && options.debugRepeatMasks)
        {
            writeMask(
                options.refDb,
                options.repeatMask ~ "-all",
                repetitiveRegions.intervals,
            );
            writeMask(
                options.refDb,
                options.repeatMask ~ "-improper",
                repetitiveRegionsImproper.intervals,
            );
        }
    }
}


/**
    Mask reference regions where the alignment coverage is not within set
    limits. This helps to identify repetitive or bad quality regions.
*/
class BadAlignmentCoverageAssessor
{
    double lowerLimit;
    double upperLimit;

    /// Create an assessor with these limits.
    this(double lowerLimit, double upperLimit)
    {
        this.lowerLimit = lowerLimit;
        this.upperLimit = upperLimit;
    }

    version (unittest)
    {
        /**
            Generate a list of alignment chains for testing:

            ```
                                 c1                             c2                    c3
                   0    5   10   15   20   25   30       0    5   10   15      0    5   10   15
            ref:   [-----------------------------)       [--------------)      [--------------)
            reads: .    .    .    .    .    .    .       .    .    .    .      .    .    .    .
                   . #1 [------------) .    .    .   #12 [--) .    .    .  #23 .[--).    .    .
                   . #2 [------------) .    .    .   #13 [--) .    .    .  #24 . [--)    .    .
                   . #3 [--------------)    .    .   #14 [----)    .    .  #25 .  [--)   .    .
                   .    . #4 [---------)    .    .   #15 [----)    .    .  #26 .   [--)  .    .
                   .    . #5 [-------------------)   #16 [--------------)  #27 .    [--) .    .
                   .    . #6 [-------------------)   #17 [--------------)  #28 .    .[--).    .
                   .    .    #7 [----------------)   #18 [--------------)  #29 .    . [--)    .
                   .    .    .    . #8 [---------)       .#19 [---------)  #30 .    .  [--)   .
                   .    .    .    . #9 [---------)       .#20 [---------)  #31 .    .   [--)  .
                   .    .    .    .#10 [---------)       .#21 [---------)  #32 .    .    [--) .
                   .    .    .    .    #11 [-----)       .    #22 [-----)  #33 .    .    .[--).
                   .    .    .    .    .    .    .       .    .    .    .      .    .    .    .
            mask:  [====)    [=======) [=========)       [==) [=========)      [==) .    . [==)
                   .    .    .    .    .    .    .       .    .    .    .      .    .    .    .
            cov:   ^    .    .    .    .    .    .       .    .    .    .      .    .    .    .
                   |    .    .    .    .    .    .       .    .    .    .      .    .    .    .
                   |    .    .  +----+ .   +----+.       +-+  .   +----+.      .    .    .    .
                   |    .    +--+####| +---+####|.       |#|  +---+####|.      .    .    .    .
                 5 |.........|#######+-+########|.       |#+--+########|.      .    .    .    .
                   |    .    |##################|.       |#############|.      .    .    .    .
                   |    +----+##################|.       |#############|.      .  +-------+   .
                   |    |#######################|.       |#############|.      . ++#######++  .
                   |    |#######################|.       |#############|.      .++#########++ .
                 0 +----+-----------------------+--------+-------------+--------+-----------+---->
            ```
        */
        private static ReferenceInterval[] getTestAlignmentIntervals()
        {
            return [
                ReferenceInterval(1,  5, 18), //  #1
                ReferenceInterval(1,  5, 18), //  #2
                ReferenceInterval(1,  5, 20), //  #3
                ReferenceInterval(1, 10, 20), //  #4
                ReferenceInterval(1, 10, 30), //  #5
                ReferenceInterval(1, 10, 30), //  #6
                ReferenceInterval(1, 13, 30), //  #7
                ReferenceInterval(1, 20, 30), //  #8
                ReferenceInterval(1, 20, 30), //  #9
                ReferenceInterval(1, 20, 30), // #10
                ReferenceInterval(1, 24, 30), // #11
                ReferenceInterval(2,  0,  3), // #12
                ReferenceInterval(2,  0,  3), // #13
                ReferenceInterval(2,  0,  5), // #14
                ReferenceInterval(2,  0,  5), // #15
                ReferenceInterval(2,  0, 15), // #16
                ReferenceInterval(2,  0, 15), // #17
                ReferenceInterval(2,  0, 15), // #18
                ReferenceInterval(2,  5, 15), // #19
                ReferenceInterval(2,  5, 15), // #20
                ReferenceInterval(2,  5, 15), // #21
                ReferenceInterval(2,  9, 15), // #22
                ReferenceInterval(3,  1,  4), // #23
                ReferenceInterval(3,  2,  5), // #24
                ReferenceInterval(3,  3,  6), // #25
                ReferenceInterval(3,  4,  7), // #26
                ReferenceInterval(3,  5,  8), // #27
                ReferenceInterval(3,  6,  9), // #28
                ReferenceInterval(3,  7, 10), // #29
                ReferenceInterval(3,  8, 11), // #30
                ReferenceInterval(3,  9, 12), // #31
                ReferenceInterval(3, 10, 13), // #32
                ReferenceInterval(3, 11, 14), // #33
            ];
        }

        private static ReferenceInterval[] getTestContigIntervals()
        {
            return [
                ReferenceInterval(1, 0, 30),
                ReferenceInterval(2, 0, 15),
                ReferenceInterval(3, 0, 15),
            ];
        }
    }

    /// Apply the assessor to the given set of alignment.
    ReferenceRegion opCall(R1, R2)(R1 alignmentIntervals, R2 contigIntervals)
        if (
            isInputRange!R1 && is(ElementType!R1 == ReferenceInterval) &&
            isInputRange!R2 && is(ElementType!R2 == ReferenceInterval)
        )
    {
        if (alignmentIntervals.empty)
            return ReferenceRegion();

        static enum OK = CoverageZone.ok;
        auto maskAcc = appender!(ReferenceInterval[]);
        auto masker = Masker();
        auto changeEvents = coverageChanges(alignmentIntervals, contigIntervals);
        auto lastEvent = changeEvents.front;

        debug logJsonDebug("changeEvents", changeEvents.array.toJson);
        foreach (event; changeEvents)
        {
            auto currentZone = coverageZone(event.currentCoverage);
            auto newZone = coverageZone(event.newCoverage);

            if (masker.isMasking && event.contigId != lastEvent.contigId)
            {
                maskAcc ~= masker.finish(lastEvent.position);
            }

            if (
                !masker.isMasking &&
                (newZone != OK || (currentZone == OK && currentZone != newZone))
            )
            {
                masker.start(event.contigId, event.position);
            }
            else if (masker.isMasking && currentZone != OK && newZone == OK)
            {
                maskAcc ~= masker.finish(event.position);
            }

            lastEvent = event;
        }

        if (masker.isMasking)
        {
            maskAcc ~= masker.finish(lastEvent.position);
        }

        return ReferenceRegion(maskAcc.data);
    }

    ///
    unittest
    {
        auto alignmentIntervals = getTestAlignmentIntervals();
        auto contigIntervals = getTestContigIntervals();
        alias CoverageChange = CoverageChangeRange.CoverageChange;

        auto assessor = new BadAlignmentCoverageAssessor(3, 5);

        assert(assessor(alignmentIntervals, contigIntervals) == ReferenceRegion([
            ReferenceInterval(1,  0,  5),
            ReferenceInterval(1, 10, 18),
            ReferenceInterval(1, 20, 30),
            ReferenceInterval(2,  0,  3),
            ReferenceInterval(2,  5, 15),
            ReferenceInterval(3,  0,  3),
            ReferenceInterval(3, 12, 15),
        ]));
    }

    private CoverageZone coverageZone(in int coverage) pure nothrow
    {
        return coverage < lowerLimit
            ? CoverageZone.low
            : coverage > upperLimit
                ? CoverageZone.high
                : CoverageZone.ok;
    }
}

enum CoverageZone
{
    low,
    ok,
    high,
}

struct Masker
{
    private bool _isMasking = false;
    private size_t contigId;
    private size_t maskStart;

    void start(in size_t contigId, in size_t maskStart) pure nothrow
    {
        this._isMasking = true;
        this.contigId = contigId;
        this.maskStart = maskStart;
    }

    ReferenceInterval finish(in size_t maskEnd) pure nothrow
    {
        this._isMasking = false;

        return ReferenceInterval(
            this.contigId,
            this.maskStart,
            maskEnd,
        );
    }

    @property bool isMasking() const pure nothrow
    {
        return this._isMasking;
    }
}

/// Transforms a range of alignment chains into a range of coverage
/// change events.
struct CoverageChangeRange
{
    static alias AlignmentEvent = Tuple!
    (
        size_t, "contigId",
        size_t, "position",
        int, "diff",
    );

    static struct CoverageChange
    {
        size_t contigId;
        size_t position;
        int currentCoverage;
        int newCoverage;

        bool hasChanged() const pure nothrow
        {
            return currentCoverage != newCoverage;
        }
    }

    AlignmentEvent[] alignmentEvents;
    size_t currentEventIdx = 0;
    CoverageChange _front;

    private this(AlignmentEvent[] alignmentEvents)
    {
        this.alignmentEvents = alignmentEvents;
        if (!empty)
        {
            // Advance to first event.
            popFront();
        }
    }

    private static CoverageChangeRange create(R1, R2)(R1 alignmentIntervals, R2 contigIntervals)
        if (
            isInputRange!R1 && is(ElementType!R1 == ReferenceInterval) &&
            isInputRange!R2 && is(ElementType!R2 == ReferenceInterval)
        )
    {
        if (alignmentIntervals.empty)
            return CoverageChangeRange();

        auto alignmentEvents = alignmentIntervals
            .map!(alignmentInterval => only(
                AlignmentEvent(alignmentInterval.contigId, alignmentInterval.begin, 1),
                AlignmentEvent(alignmentInterval.contigId, alignmentInterval.end, -1),
            ))
            .joiner;
        auto contigBoundaryEvents = contigIntervals
            .map!(contigInterval => only(
                AlignmentEvent(contigInterval.contigId, 0, 0),
                AlignmentEvent(contigInterval.contigId, contigInterval.size, 0),
            ))
            .joiner;
        auto changeEvents = chain(alignmentEvents, contigBoundaryEvents).array;
        changeEvents.sort();

        return CoverageChangeRange(changeEvents);
    }

    void popFront()
    {
        _front.currentCoverage = _front.newCoverage;

        if (currentEventIdx == alignmentEvents.length)
        {
            // There is no more data; just pop the front and return;
            assert(_front.currentCoverage == 0, "coverage should drop to zero in the end");
            ++currentEventIdx;

            return;
        }

        auto currentEvent = alignmentEvents[currentEventIdx];
        auto eventAcc = currentEvent;
        eventAcc.diff = 0;

        // Collect all events for one position.
        for (; currentEventIdx < alignmentEvents.length; ++currentEventIdx)
        {
            currentEvent = alignmentEvents[currentEventIdx];

            if (
                eventAcc.contigId != currentEvent.contigId ||
                eventAcc.position != currentEvent.position
            )
            {
                break;
            }

            eventAcc.diff += currentEvent.diff;
        }

        assert(eventAcc.contigId == _front.contigId || _front.currentCoverage == 0,
                "coverage should drop to zero between contigs");
        _front.contigId = eventAcc.contigId;
        _front.position = eventAcc.position;
        _front.newCoverage = _front.currentCoverage + eventAcc.diff;
    }

    @property bool empty()
    {
        return currentEventIdx > alignmentEvents.length;
    }

    @property CoverageChange front()
    {
        return _front;
    }

    @property CoverageChangeRange save() pure nothrow @safe
    {
        // nothing to be done
        return this;
    }
}

CoverageChangeRange coverageChanges(R1, R2)(R1 alignmentIntervals, R2 contigIntervals)
    if (
        isInputRange!R1 && is(ElementType!R1 == ReferenceInterval) &&
        isInputRange!R2 && is(ElementType!R2 == ReferenceInterval)
    )
{
    return CoverageChangeRange.create(alignmentIntervals, contigIntervals);
}

unittest
{
    import std.algorithm : equal;

    auto alignmentIntervals = BadAlignmentCoverageAssessor.getTestAlignmentIntervals();
    auto contigIntervals = BadAlignmentCoverageAssessor.getTestContigIntervals();

    alias CoverageChange = CoverageChangeRange.CoverageChange;

    assert(coverageChanges(alignmentIntervals, contigIntervals).equal([
        CoverageChange(1,  0, 0, 0),
        CoverageChange(1,  5, 0, 3),
        CoverageChange(1, 10, 3, 6),
        CoverageChange(1, 13, 6, 7),
        CoverageChange(1, 18, 7, 5),
        CoverageChange(1, 20, 5, 6),
        CoverageChange(1, 24, 6, 7),
        CoverageChange(1, 30, 7, 0),
        CoverageChange(2,  0, 0, 7),
        CoverageChange(2,  3, 7, 5),
        CoverageChange(2,  5, 5, 6),
        CoverageChange(2,  9, 6, 7),
        CoverageChange(2, 15, 7, 0),
        CoverageChange(3,  0, 0, 0),
        CoverageChange(3,  1, 0, 1),
        CoverageChange(3,  2, 1, 2),
        CoverageChange(3,  3, 2, 3),
        CoverageChange(3,  4, 3, 3),
        CoverageChange(3,  5, 3, 3),
        CoverageChange(3,  6, 3, 3),
        CoverageChange(3,  7, 3, 3),
        CoverageChange(3,  8, 3, 3),
        CoverageChange(3,  9, 3, 3),
        CoverageChange(3, 10, 3, 3),
        CoverageChange(3, 11, 3, 3),
        CoverageChange(3, 12, 3, 2),
        CoverageChange(3, 13, 2, 1),
        CoverageChange(3, 14, 1, 0),
        CoverageChange(3, 15, 0, 0),
    ]));
}

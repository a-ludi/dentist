/**
    This is a collection of helpers for the assessment of the repeat mask.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.collectPileUps.maskassessment;

import dentist.common.alignments :
    id_t,
    coord_t,
    AlignmentChain;
import dentist.common : ReferenceInterval, ReferenceRegion;
import std.algorithm :
    equal,
    joiner,
    map,
    sort,
    uniq;
import std.array : appender, array;
import std.exception : enforce;
import std.typecons : Tuple;
import std.range :
    chain,
    empty,
    ElementType,
    front,
    isInputRange,
    only,
    popFront;

// see `assessmentStage` below.
alias MaskAssessmentStage(Assessor : RepeatAssessor) = Tuple!(
    string, "name",
    Assessor, "assessor",
    const(AlignmentChain[]), "input",
);

/// Construct a stage of the repeat mask assessment pipeline.
MaskAssessmentStage!Assessor assessmentStage(Assessor : RepeatAssessor)(
    string name,
    Assessor assessor,
    const(AlignmentChain[]) input,
)
{
    return typeof(return)(name, assessor, input);
}

/**
    A `RepeatAssessor` should generate a repeat mask derived from a set of
    alignment chains.
*/
interface RepeatAssessor
{
    ReferenceRegion opCall(in AlignmentChain[] input);
}

/**
    Mask reference regions where the alignment coverage is not within set
    limits. This helps to identify repetitive or bad quality regions.
*/
class BadAlignmentCoverageAssessor : RepeatAssessor
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
        private static AlignmentChain[] getTestAlignments()
        {
            with (AlignmentChain) with (LocalAlignment)
                {
                    id_t alignmentChainId = 0;
                    id_t contReadId = 0;
                    AlignmentChain getDummyAlignment(id_t contigId,
                            coord_t contigLength, coord_t beginIdx, coord_t endIdx)
                    {
                        return AlignmentChain(
                            ++alignmentChainId,
                            Contig(contigId, contigLength),
                            Contig(alignmentChainId, endIdx),
                            emptyFlags,
                            [
                                LocalAlignment(
                                    Locus(beginIdx, beginIdx + 1),
                                    Locus(0, 1),
                                    0,
                                ),
                                LocalAlignment(
                                    Locus(endIdx - 1, endIdx),
                                    Locus(2, 3),
                                    0,
                                ),
                            ],
                        );
                    }

                    return [
                        getDummyAlignment(1, 30,  5, 18), //  #1
                        getDummyAlignment(1, 30,  5, 18), //  #2
                        getDummyAlignment(1, 30,  5, 20), //  #3
                        getDummyAlignment(1, 30, 10, 20), //  #4
                        getDummyAlignment(1, 30, 10, 30), //  #5
                        getDummyAlignment(1, 30, 10, 30), //  #6
                        getDummyAlignment(1, 30, 13, 30), //  #7
                        getDummyAlignment(1, 30, 20, 30), //  #8
                        getDummyAlignment(1, 30, 20, 30), //  #9
                        getDummyAlignment(1, 30, 20, 30), // #10
                        getDummyAlignment(1, 30, 24, 30), // #11
                        getDummyAlignment(2, 15,  0,  3), // #12
                        getDummyAlignment(2, 15,  0,  3), // #13
                        getDummyAlignment(2, 15,  0,  5), // #14
                        getDummyAlignment(2, 15,  0,  5), // #15
                        getDummyAlignment(2, 15,  0, 15), // #16
                        getDummyAlignment(2, 15,  0, 15), // #17
                        getDummyAlignment(2, 15,  0, 15), // #18
                        getDummyAlignment(2, 15,  5, 15), // #19
                        getDummyAlignment(2, 15,  5, 15), // #20
                        getDummyAlignment(2, 15,  5, 15), // #21
                        getDummyAlignment(2, 15,  9, 15), // #22
                        getDummyAlignment(3, 15,  1,  4), // #23
                        getDummyAlignment(3, 15,  2,  5), // #24
                        getDummyAlignment(3, 15,  3,  6), // #25
                        getDummyAlignment(3, 15,  4,  7), // #26
                        getDummyAlignment(3, 15,  5,  8), // #27
                        getDummyAlignment(3, 15,  6,  9), // #28
                        getDummyAlignment(3, 15,  7, 10), // #29
                        getDummyAlignment(3, 15,  8, 11), // #30
                        getDummyAlignment(3, 15,  9, 12), // #31
                        getDummyAlignment(3, 15, 10, 13), // #32
                        getDummyAlignment(3, 15, 11, 14), // #33
                    ];
                }
        }
    }

    /// Apply the assessor to the given set of alignment.
    override ReferenceRegion opCall(const(AlignmentChain[]) alignments)
    {
        if (alignments.length == 0)
        {
            return ReferenceRegion();
        }

        static immutable OK = CoverageZone.ok;
        auto maskAcc = appender!(ReferenceInterval[]);
        auto masker = Masker();
        auto changeEvents = coverageChanges(alignments);
        auto lastEvent = changeEvents.front;

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
        auto alignments = getTestAlignments();
        alias CoverageChange = CoverageChangeRange.CoverageChange;

        auto assessor = new BadAlignmentCoverageAssessor(3, 5);

        assert(assessor(alignments) == ReferenceRegion([
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

private enum CoverageZone
{
    low,
    ok,
    high,
}

private struct Masker
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
private struct CoverageChangeRange
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

    private static CoverageChangeRange create(Range)(Range alignments)
            if (isInputRange!Range && is(ElementType!Range : const(AlignmentChain)))
    {
        if (alignments.length == 0)
        {
            return CoverageChangeRange();
        }

        auto alignmentEvents = alignments
            .map!(alignment => only(
                AlignmentEvent(alignment.contigA.id, alignment.first.contigA.begin, 1),
                AlignmentEvent(alignment.contigA.id, alignment.last.contigA.end, -1),
            ))
            .joiner;
        auto contigBoundaryEvents = alignments
            .map!"a.contigA"
            .uniq
            .map!(contig => only(
                AlignmentEvent(contig.id, 0, 0),
                AlignmentEvent(contig.id, contig.length, 0),
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
}

CoverageChangeRange coverageChanges(Range)(Range alignments)
{
    return CoverageChangeRange.create(alignments);
}

unittest
{
    auto alignments = BadAlignmentCoverageAssessor.getTestAlignments();

    alias CoverageChange = CoverageChangeRange.CoverageChange;

    assert(coverageChanges(alignments).equal([
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

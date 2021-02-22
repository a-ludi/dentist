/**
    This is the `checkScaffolding` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.checkScaffolding;

import dentist.common : isTesting;

static if (isTesting):

import dentist.commandline : OptionsFor;
import dentist.commands.checkResults :
    Complement,
    ContigAlignmentsCache,
    ContigMapping;
import dentist.common :
    dentistEnforce,
    ReferenceInterval,
    ReferenceRegion;
import dentist.common.alignments.base;
import dentist.common.binio.insertiondb : InsertionDb;
import dentist.common.commands : TestingCommand;
import dentist.common.insertions :
    getInfoForExistingContig,
    getInfoForNewSequenceInsertion,
    InsertionInfo,
    isOutputGap;
import dentist.common.scaffold :
    buildScaffold,
    isAntiParallel,
    isDefault,
    linearWalk,
    scaffoldStarts;
import dentist.dazzler :
    DBdumpOptions,
    DbRecord,
    getContigCutoff,
    getDbRecords,
    readDazzExtra,
    readMask;
import dentist.util.log;
import dentist.util.range : arrayChunks;
import std.algorithm :
    equal,
    chunkBy,
    copy,
    filter,
    find,
    isSorted,
    max,
    maxElement,
    min,
    minElement,
    map,
    reverse,
    sort,
    until;
import std.array : array;
import std.conv : to;
import std.format : format;
import std.range :
    assumeSorted,
    chain,
    enumerate,
    iota,
    only,
    slide;
import std.range.primitives;
import std.stdio : writeln;
import std.typecons : No;
import vibe.data.json :
    Json,
    toJson = serializeToJson,
    toJsonCompressed = serializeToJsonString,
    toJsonString = serializeToPrettyJson;


/// Options for the `collectPileUps` command.
alias Options = OptionsFor!(TestingCommand.checkScaffolding);

/// Execute the `checkScaffolding` command with `options`.
void execute(in Options options)
{
    auto analyzer = ResultAnalyzer(options);

    analyzer.run();
}


bool referenceOrder(const ContigMapping lhs, const ContigMapping rhs) pure nothrow @safe
{
    return lhs.reference < rhs.reference;
}


bool queryOrder(const ContigMapping lhs, const ContigMapping rhs) pure nothrow @safe
{
    return lhs.queryContigId < rhs.queryContigId;
}


private struct ResultAnalyzer
{
    static enum JoinState : ubyte
    {
        /// Ignore this join
        ignored,
        /// The involved contigs are joined as in the test assembly.
        correct,
        /// Two scaffold ends are joined.
        novel,
        /// The join contradicts the test assembly.
        broken,
    }

    static struct JoinSummary
    {
        JoinState state;
        ContigMapping lhsContigMapping;
        ContigMapping rhsContigMapping;
        ContigMapping[] skippedContigMappings;


        string toJson() const
        {
            return format!`{"state":"%s","lhsContigMapping":%s,"rhsContigMapping":%s,"skippedContigMappings":%s}`(
                state.to!string,
                lhsContigMapping.toJsonCompressed,
                rhsContigMapping.toJsonCompressed,
                skippedContigMappings.toJsonCompressed,
            );
        }
    }

    static struct JoinedContigs
    {
        id_t lhs;
        id_t rhs;
    }

    enum dbdumpOptions = [DBdumpOptions.readNumber, DBdumpOptions.originalHeader];

    const(Options) options;
    protected ReferenceInterval[] mappedContigs;
    protected DbRecord[] trueScaffoldStructure;
    protected DbRecord[] resultScaffoldStructure;
    protected ContigMapping[] knownContigMappings;
    protected size_t[] gapStartIndices;
    protected ContigMapping[] contigMappings;
    protected JoinSummary[] joinSummaries;

    void run()
    {
        mixin(traceExecution);

        init();
        analyzeJoins();
        writeScaffolingReportJson();
    }

    void init()
    {
        mixin(traceExecution);

        auto contigCutoff = getContigCutoff(options.refDb);
        mappedContigs = readMask!ReferenceInterval(
            options.trueAssemblyDb,
            options.mappedRegionsMask,
        ).filter!(interval => interval.size >= contigCutoff).array;
        enforceMappedAndTestAssemblyContigLengthsMatch();
        trueScaffoldStructure = getDbRecords(options.trueAssemblyDb, dbdumpOptions).array;
        resultScaffoldStructure = getDbRecords(options.resultDb, dbdumpOptions).array;

        auto contigAlignmentsCache = ContigAlignmentsCache(
            options.contigAlignmentsCache,
            options.resultDb,
            options.refDb,
        );

        dentistEnforce(contigAlignmentsCache.isValid, "invalid contig-alignments-cache");

        contigMappings = contigAlignmentsCache.read();
        contigMappings.sort!queryOrder;
        fakeContigMappingsFromScaffolding();

        logJsonDebug("contigMappings", contigMappings.toJson);
    }

    void enforceMappedAndTestAssemblyContigLengthsMatch()
    {
        auto mappedContigLengths = mappedContigs.map!"a.size";
        auto testAssemblyContigLengths = getDbRecords(options.refDb, dbdumpOptions)
            .map!"a.location.length";

        dentistEnforce(
            equal(mappedContigLengths, testAssemblyContigLengths),
            "lengths of mapped contigs differ from lengths of test assembly " ~
            "contigs. Check if the true assembly and mapped mask match the " ~
            "given test assembly; if not there mgith be a bug in " ~
            "check-scaffolding."
        );
    }

    void fakeContigMappingsFromScaffolding()
    {
        mixin(traceExecution);

        auto insertions = InsertionDb.parse(options.assemblyGraphFile)[];
        auto assemblyGraph = buildScaffold(insertions);
        auto incidentEdgesCache = assemblyGraph.allIncidentEdges();
        auto scaffoldStartNodes = scaffoldStarts!InsertionInfo(assemblyGraph, incidentEdgesCache).array;

        knownContigMappings.length = 0;
        knownContigMappings.reserve(mappedContigs.length);
        gapStartIndices.length = 0;
        gapStartIndices.reserve(mappedContigs.length);

        id_t currentContigId = 1;
        // for every scaffold in the graph
        foreach (startNode; scaffoldStartNodes)
        {
            auto globalComplement = false;
            auto insertionBegin = startNode;
            coord_t currentContigCoord;

            alias getFinalContigLength = () =>
                getDbRecord(resultScaffoldStructure, currentContigId).location.length;

            // walk the scaffold path
            foreach (insertion; linearWalk!InsertionInfo(assemblyGraph, startNode, incidentEdgesCache))
            {
                if (insertion.isDefault)
                {
                    // this is an input contig copied to the result
                    auto insertionInfo = getInfoForExistingContig(insertionBegin, insertion, globalComplement);
                    auto nextContigCoord = currentContigCoord + cast(coord_t) insertionInfo.length;

                    knownContigMappings ~= ContigMapping(
                        ReferenceInterval(
                            currentContigId,
                            currentContigCoord,
                            nextContigCoord,
                        ),
                        getFinalContigLength(),
                        insertionInfo.contigId,
                        No.duplicateQueryContig, // we don't really know – does not matter
                        cast(Complement) insertionInfo.complement,
                        0.0, // should be a perfect match :P
                    );

                    currentContigCoord = nextContigCoord;
                }
                else if (insertion.isOutputGap)
                {
                    assert(currentContigCoord == getFinalContigLength());
                    // goto next contig
                    ++currentContigId;
                    currentContigCoord = 0;
                }
                else
                {
                    // close the gap
                    assert(currentContigCoord > 0);
                    auto insertionInfo = getInfoForNewSequenceInsertion(insertionBegin, insertion, globalComplement);
                    currentContigCoord += cast(coord_t) insertionInfo.length;
                    assert(knownContigMappings.length > 0);
                    gapStartIndices ~= knownContigMappings.length - 1;
                }

                insertionBegin = insertion.target(insertionBegin);
                if (insertion.isAntiParallel)
                    globalComplement = !globalComplement;
            }

            assert(currentContigCoord == getFinalContigLength());
            ++currentContigId;
        }

        knownContigMappings.sort!referenceOrder;
    }

    void analyzeJoins()
    {
        mixin(traceExecution);

        joinSummaries.reserve(gapStartIndices.length);
        foreach (gapStartIndex; gapStartIndices)
        {
            auto lhsContigMapping = knownContigMappings[gapStartIndex];
            auto rhsContigMapping = knownContigMappings[gapStartIndex + 1];

            if (!onSameResultContig(lhsContigMapping, rhsContigMapping))
                // moving to the next contig (scaffold gaps are the same
                // as in the reference assembly and can be ignored in this
                // analysis)
                continue;

            JoinSummary joinSummary;
            joinSummary.lhsContigMapping = lhsContigMapping;
            joinSummary.rhsContigMapping = rhsContigMapping;

            if (adjacentInTrueAssembly(lhsContigMapping, rhsContigMapping))
            {
                joinSummary.state = JoinState.correct;
            }
            else if (maybeAdjacentInTrueAssembly(lhsContigMapping, rhsContigMapping))
            {
                if (skippedContigsArePresent(
                    lhsContigMapping,
                    rhsContigMapping,
                    joinSummary.skippedContigMappings
                ))
                    joinSummary.state = JoinState.correct;
                else
                    joinSummary.state = JoinState.broken;
            }
            else if (
                endOfTrueAssemblyScaffold(lhsContigMapping) &&
                endOfTrueAssemblyScaffold(rhsContigMapping)
            )
            {
                joinSummary.state = JoinState.novel;
            }
            else
            {
                joinSummary.state = JoinState.broken;
            }

            joinSummaries ~= joinSummary;
        }
    }

    bool onSameResultContig(const ContigMapping lhs, const ContigMapping rhs) const pure nothrow @safe
    {
        return lhs.reference.contigId == rhs.reference.contigId;
    }

    bool areOrderedInTheResult(const ContigMapping lhs, const ContigMapping rhs) const pure nothrow @safe
    {
        return areOrderedInTheResult(lhs.reference, rhs.reference);
    }

    bool areOrderedInTheResult(const ReferenceInterval lhs, const ReferenceInterval rhs) const pure nothrow @safe
    {
        return lhs < rhs && (lhs & rhs).size <= options.properAlignmentAllowance;
    }

    bool nonoverlappingInResult(const ContigMapping lhs, const ContigMapping rhs) const pure nothrow @safe
    {
        return (lhs.reference & rhs.reference).empty;
    }

    bool adjacentInTrueAssembly(const ContigMapping lhs, const ContigMapping rhs) const pure nothrow @safe
    {
        auto lhsContig = mappedContigs[lhs.queryContigId - 1];
        auto lhsDbRecord = getDbRecord(trueScaffoldStructure, lhsContig.contigId);
        auto rhsContig = mappedContigs[rhs.queryContigId - 1];
        auto rhsDbRecord = getDbRecord(trueScaffoldStructure, rhsContig.contigId);

        id_t lhsIncrement = lhs.complement ? 0 : 1;
        id_t rhsIncrement = 1 - lhsIncrement;

        return (
            // part of the same ground-truth scaffold
            lhsDbRecord.header == rhsDbRecord.header &&
            // same orientation
            lhs.complement == rhs.complement &&
            // correct contig sequence
            lhs.queryContigId + lhsIncrement == rhs.queryContigId + rhsIncrement
        );
    }

    bool maybeAdjacentInTrueAssembly(const ContigMapping lhs, const ContigMapping rhs) const pure nothrow @safe
    {
        auto lhsContig = mappedContigs[lhs.queryContigId - 1];
        auto lhsDbRecord = getDbRecord(trueScaffoldStructure, lhsContig.contigId);
        auto rhsContig = mappedContigs[rhs.queryContigId - 1];
        auto rhsDbRecord = getDbRecord(trueScaffoldStructure, rhsContig.contigId);

        return (
            // part of the same ground-truth scaffold
            lhsDbRecord.header == rhsDbRecord.header &&
            // same orientation
            lhs.complement == rhs.complement &&
            // correct contig order but possibly with skipped contigs
            lhs.complement
                ? lhs.queryContigId > rhs.queryContigId
                : lhs.queryContigId < rhs.queryContigId
        );
    }

    bool skippedContigsArePresent(
        const ContigMapping lhs,
        const ContigMapping rhs,
        out ContigMapping[] skippedContigMappings,
    ) const pure @safe
    {
        auto gap = ReferenceInterval(
            lhs.reference.contigId,
            lhs.reference.end,
            rhs.reference.begin,
        );
        // reduce contigMappings to the skipped query contigs
        auto candidateContigMappings = contigMappings
            .assumeSorted!queryOrder
            .upperBound(lhs)
            .lowerBound(rhs)
            .release();
        alias ContigIdsRange = typeof(iota(arithmetic_t.init, arithmetic_t.init, arithmetic_t.init));
        auto skippedContigIds = lhs.complement
            ? ContigIdsRange(lhs.queryContigId - 1, rhs.queryContigId, -1)
            : ContigIdsRange(lhs.queryContigId + 1, rhs.queryContigId, 1);
        ContigMapping prevContigMapping = knownContigMappings[lhs.queryContigId - 1];
        skippedContigMappings.reserve(skippedContigIds.length);
        foreach (skippedContigId; skippedContigIds)
        {
            auto validJoinsFinder = candidateContigMappings
                .find!(cm =>
                    // query contig must match
                    cm.queryContigId == skippedContigId &&
                    // the mapping must be located inside the gap (with some allowance)
                    (cm.reference - gap).size <= options.properAlignmentAllowance &&
                    // the scaffolding must match the ground-truth
                    adjacentInTrueAssembly(prevContigMapping, cm) &&
                    // the ordering in the result assembly must be as expected
                    areOrderedInTheResult(prevContigMapping, cm)
                );

            if (validJoinsFinder.empty)
            {
                // could not find a valid join -> gap broken :-(
                skippedContigMappings = candidateContigMappings.dup;

                return false;
            }

            // found a valid join -> take it and go to the next skipped contig
            skippedContigMappings ~= validJoinsFinder.front;
            prevContigMapping = validJoinsFinder.front;
        }

        // got chain of valid joins; let's see if it connects to the rhs
        return adjacentInTrueAssembly(prevContigMapping, rhs);
    }

    bool endOfTrueAssemblyScaffold(const ContigMapping contigMapping) const pure nothrow @safe
    {
        auto queryContigId = contigMapping.queryContigId;
        auto thisContigId = mappedContigs[queryContigId - 1].contigId;
        auto prevContigId = queryContigId > 1
            ? mappedContigs[queryContigId - 2].contigId
            : cast(id_t) 0;
        auto nextContigId = queryContigId < mappedContigs.length
            ? mappedContigs[queryContigId].contigId
            : thisContigId + 1;
        auto thisContigScaffold = getDbRecord(trueScaffoldStructure, thisContigId).header;
        auto prevContigScaffold = 0 < prevContigId
            ? getDbRecord(trueScaffoldStructure, prevContigId).header
            : null;
        auto nextContigScaffold = nextContigId < trueScaffoldStructure.length
            ? getDbRecord(trueScaffoldStructure, nextContigId).header
            : null;

        return (
            // contig is the first on its scaffold
            prevContigScaffold != thisContigScaffold ||
            // contig is the last on its scaffold
            thisContigScaffold != nextContigScaffold
        );
    }

    bool endOfResultAssemblyContig(const ContigMapping contigMapping) const pure nothrow @safe
    {
        return contigMapping.reference.begin <= options.properAlignmentAllowance ||
               contigMapping.reference.end + options.properAlignmentAllowance >= contigMapping.referenceContigLength;
    }

    const(DbRecord) getDbRecord(const DbRecord[] dbRecords, size_t contigId) const pure nothrow @safe
    {
        auto dbRecord = dbRecords[contigId - 1];
        assert(dbRecord.contigId == contigId);

        return dbRecord;
    }

    void writeScaffolingReportJson()
    {
        mixin(traceExecution);

        foreach (joinSummary; joinSummaries)
            if (joinSummary.state != JoinState.ignored)
                writeln(joinSummary.toJson);
    }
}

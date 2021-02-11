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
import dentist.common : dentistEnforce;
import dentist.common.alignments.base;
import dentist.common.commands : TestingCommand;
import dentist.dazzler :
    DBdumpOptions,
    DbRecord,
    getDbRecords;
import dentist.util.log;
import std.algorithm :
    sort;
import std.array : array;
import std.conv : to;
import std.format : format;
import std.range :
    enumerate,
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


        string toJson() const
        {
            return format!`{"state":"%s","lhsContigMapping":%s,"rhsContigMapping":%s}`(
                state.to!string,
                lhsContigMapping.toJsonCompressed,
                rhsContigMapping.toJsonCompressed,
            );
        }
    }

    const(Options) options;
    protected DbRecord[] expectedScaffoldStructure;
    protected DbRecord[] resultScaffoldStructure;
    protected ContigMapping[] contigAlignments;
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

        enum dbdumpOptions = [DBdumpOptions.readNumber, DBdumpOptions.originalHeader];
        expectedScaffoldStructure = getDbRecords(options.refDb, dbdumpOptions).array;
        resultScaffoldStructure = getDbRecords(options.resultDb, dbdumpOptions).array;

        auto contigAlignmentsCache = ContigAlignmentsCache(
            options.contigAlignmentsCache,
            options.resultDb,
            options.refDb,
        );

        dentistEnforce(contigAlignmentsCache.isValid, "invalid contig-alignments-cache");

        contigAlignments = contigAlignmentsCache.read();
    }

    void analyzeJoins()
    {
        mixin(traceExecution);

        contigAlignments.sort!((a, b) => a.reference < b.reference);
        joinSummaries.length = contigAlignments.length;
        foreach (i, contigMappingPair; contigAlignments.slide!(No.withPartial)(2).enumerate)
        {
            auto lhsContigMapping = contigMappingPair[0];
            auto rhsContigMapping = contigMappingPair[1];

            if (!onSameResultContig(lhsContigMapping, rhsContigMapping))
                continue; // joinState == ignored

            auto joinSummary = &joinSummaries[i];
            joinSummary.lhsContigMapping = lhsContigMapping;
            joinSummary.rhsContigMapping = rhsContigMapping;

            if (adjacentInTestAssembly(lhsContigMapping, rhsContigMapping))
            {
                joinSummary.state = JoinState.correct;
            }
            else if (
                endOfTestAssemblyScaffold(lhsContigMapping) &&
                endOfTestAssemblyScaffold(rhsContigMapping)
            )
            {
                joinSummary.state = JoinState.novel;
            }
            else
            {
                joinSummary.state = JoinState.broken;
            }
        }
    }

    bool onSameResultContig(const ContigMapping lhs, const ContigMapping rhs) const pure nothrow @safe
    {
        return lhs.reference.contigId == rhs.reference.contigId;
    }

    bool nonoverlappingInResult(const ContigMapping lhs, const ContigMapping rhs) const pure nothrow @safe
    {
        return (lhs.reference & rhs.reference).empty;
    }

    bool adjacentInTestAssembly(const ContigMapping lhs, const ContigMapping rhs) const pure nothrow @safe
    {
        auto lhsDbRecord = getDbRecord(expectedScaffoldStructure, lhs.queryContigId);
        auto rhsDbRecord = getDbRecord(expectedScaffoldStructure, rhs.queryContigId);

        id_t lhsIncrement = lhs.complement ? 0 : 1;
        id_t rhsIncrement = 1 - lhsIncrement;

        return (
            lhsDbRecord.header == rhsDbRecord.header &&
            lhs.complement == rhs.complement &&
            lhs.queryContigId + lhsIncrement == rhs.queryContigId + rhsIncrement
        );
    }

    bool endOfTestAssemblyScaffold(const ContigMapping contigMapping) const pure nothrow @safe
    {
        auto dbRecord = getDbRecord(expectedScaffoldStructure, contigMapping.queryContigId);
        auto nextDbRecord = contigMapping.queryContigId < expectedScaffoldStructure.length
            ? getDbRecord(expectedScaffoldStructure, contigMapping.queryContigId + 1)
            : DbRecord();

        return (
            // contig is the first on its scaffold
            dbRecord.location.contigIdx == 0 ||
            // next contig is the first on its scaffold or undefined
            nextDbRecord.location.contigIdx == 0
        );
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

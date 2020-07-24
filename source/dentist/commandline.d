/**
    Defines the behavior of the `dentist` command line client.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commandline;

import darg :
    ArgParseError,
    ArgParseHelp,
    Argument,
    ArgumentsParser,
    handleArg,
    Help,
    helpString,
    MetaVar,
    Multiplicity,
    Option,
    OptionFlag,
    parseArgs,
    usageString;
import dentist.common :
    isTesting,
    OutputCoordinate,
    testingOnly;
import dentist.common.alignments :
    coord_t,
    id_t,
    trace_point_t;
import dentist.common.commands :
    DentistCommand,
    dentistCommands,
    TestingCommand;
import dentist.common.configfile :
    retroInitFromConfig,
    fromBytes,
    SizeUnit,
    toBytes;
import dentist.common.binio : PileUpDb;
import dentist.common.external :
    ExternalDependency,
    externalDependencies;
import dentist.common.scaffold : JoinPolicy;
import dentist.dazzler :
    DaccordOptions,
    DalignerOptions,
    DamapperOptions,
    DatanderOptions,
    dbdustMaskName,
    getHiddenDbFiles,
    getMaskFiles,
    getNumContigs,
    getTracePointDistance,
    lasEmpty,
    LasFilterAlignmentsOptions;
import dentist.util.process : isExecutable;
import dentist.swinfo :
    copyright,
    description,
    executableName,
    gitCommit,
    license,
    version_;
import dentist.util.algorithm : staticPredSwitch;
import dentist.util.log;
import dentist.util.tempfile : mkdtemp;
import dentist.util.string : dashCaseCT, toString;
import std.algorithm :
    among,
    canFind,
    each,
    endsWith,
    filter,
    find,
    map,
    max,
    sort,
    startsWith,
    sum,
    swap;
import std.array :
    array,
    split;
import std.conv;
import std.exception :
    basicExceptionCtors,
    enforce,
    ErrnoException;
import std.file :
    exists,
    FileException,
    getcwd,
    isDir,
    mkdirRecurse,
    tempDir,
    remove,
    rmdirRecurse;
import std.format :
    format,
    formattedRead;
import std.math :
    ceil,
    exp,
    floor,
    log_e = log;
import std.meta :
    Alias,
    AliasSeq,
    allSatisfy,
    staticMap,
    staticSort;
import std.parallelism :
    defaultPoolThreads,
    totalCPUs;
import std.path :
    absolutePath,
    buildPath;
import std.process :
    Config,
    execute;
import std.range :
    ElementType,
    enumerate,
    only,
    takeOne;
import std.regex :
    ctRegex,
    matchFirst;
import std.stdio :
    File,
    stderr;
import std.string :
    join,
    lineSplitter,
    tr,
    wrap;
import std.traits :
    arity,
    EnumMembers,
    getSymbolsByUDA,
    getUDAs,
    isArray,
    isAssignable,
    isCallable,
    isDynamicArray,
    isFloatingPoint,
    isIntegral,
    isSomeString,
    isStaticArray,
    isUnsigned,
    Parameters,
    ReturnType;
import std.typecons :
    BitFlags,
    tuple,
    Yes;
import transforms : camelCase;
import vibe.data.json :
    serializeToJsonString,
    toJson = serializeToJson;


/// Possible returns codes of the command line execution.
enum ReturnCode
{
    ok,
    commandlineError,
    runtimeError,
}

/// Start `dentist` with the given set of arguments.
ReturnCode run(in string[] args)
{
    if (args.length == 1)
    {
        printBaseHelp();

        return ReturnCode.commandlineError;
    }

    switch (args[1])
    {
    case "--version":
        printVersion();

        return ReturnCode.ok;
    case "-h":
        goto case;
    case "--help":
        printBaseHelp();

        return ReturnCode.ok;
    case "-d":
        goto case;
    case "--dependencies":
        printExternalDependencies();

        return ReturnCode.ok;
    case "--usage":
        stderr.write(usageString!BaseOptions(executableName));

        return ReturnCode.ok;
    default:
        break;
    }

    string commandName;
    try
    {
        commandName = parseCommandName(args);
    }
    catch (Exception e)
    {
        stderr.writeln("Error: " ~ (shouldLog(LogLevel.diagnostic)
            ? e.to!string
            : e.msg));
        stderr.writeln();
        stderr.write(usageString!BaseOptions(executableName));

        return ReturnCode.commandlineError;
    }

    auto commandWithArgs = args[1 .. $];
    DentistCommand command = parseArgs!BaseOptions([commandName]).command;

    try
    {
        final switch (command)
        {
            static foreach (caseCommand; EnumMembers!DentistCommand)
            {
                case caseCommand:
                    return runCommand!caseCommand(commandWithArgs);
            }
        }
    }
    catch (Exception e)
    {
        stderr.writeln(shouldLog(LogLevel.diagnostic)
            ? e.to!string
            : e.msg);

        return ReturnCode.runtimeError;
    }
}

unittest
{
    import std.stdio : File;

    auto _stderr = stderr;
    stderr = File("/dev/null", "w");

    scope (exit)
    {
        stderr = _stderr;
    }

    assert(run([executableName, "--help"]) == ReturnCode.ok);
    assert(run([executableName, "--usage"]) == ReturnCode.ok);
    assert(run([executableName, "--version"]) == ReturnCode.ok);
    assert(run([executableName, "foobar"]) == ReturnCode.commandlineError);
    assert(run([executableName, "--foo"]) == ReturnCode.commandlineError);
    assert(run([executableName]) == ReturnCode.commandlineError);
}

void printExternalDependencies()
{
    import std.stdio : writefln;

    static assert(externalDependencies.length > 0);

    writefln!"%-(%s\n%)"(externalDependencies);
}

void assertExternalToolsAvailable()
{
    static assert(externalDependencies.length > 0);

    auto missingExternalDependencies = externalDependencies
        .filter!(extDep => !isExecutable(extDep.executable))
        .array;

    enforce!CLIException(
        missingExternalDependencies.length == 0,
        format!"missing external tools:\n%-(- %s\n%)\n\nCheck your PATH and/or install the required software."(
            missingExternalDependencies,
        ),
    );
}

string parseCommandName(in string[] args)
{
    enforce!CLIException(!args[1].startsWith("-"), format!"Missing <command> '%s'"(args[1]));

    auto candidates = only(dentistCommands).filter!(cmd => cmd.startsWith(args[1]));

    enforce!CLIException(!candidates.empty, format!"Unkown <command> '%s'"(args[1]));

    auto dashCaseCommand = candidates.front;

    candidates.popFront();
    enforce!CLIException(candidates.empty, format!"Ambiguous <command> '%s'"(args[1]));

    return dashCaseCommand.tr("-", "_").camelCase;
}

private void printBaseHelp()
{
    stderr.write(usageString!BaseOptions(executableName));
    stderr.writeln();
    stderr.writeln(description);
    stderr.writeln();
    stderr.write(helpString!BaseOptions);
}

private void printVersion()
{
    stderr.writeln(format!"%s %s (commit %s)"(executableName, version_, gitCommit));
    stderr.writeln();
    stderr.write(copyright);
    stderr.writeln();
    stderr.write(license);
}

/// The set of options common to all stages.
mixin template HelpOption()
{
    @Option("help", "h")
    @Help("Prints this help.")
    OptionFlag help;

    @Option("usage")
    @Help("Print a short command summary.")
    void requestUsage() pure
    {
        enforce!UsageRequested(false, "usage requested");
    }
}

class UsageRequested : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

/// Options for the different commands.
struct OptionsFor(DentistCommand _command)
{
    enum command = _command;
    enum commandName = dentistCommands[_command];

    static enum needWorkdir = command.among(
        TestingCommand.translocateGaps,
        TestingCommand.findClosableGaps,
        DentistCommand.maskRepetitiveRegions,
        DentistCommand.showMask,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
        DentistCommand.output,
        TestingCommand.checkResults,
    );

    @Option()
    string executableVersion = version_;

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
        DentistCommand.chainLocalAlignments,
    ))
    {
        @ArgumentsParser
        auto parseArguments(const(string)[] leftOver)
        {
            alias referenceSymbol = Alias!(__traits(getMember, this, "refDb"));
            enum referenceUDA = getUDAs!(referenceSymbol, Argument)[0];

            enforce!ArgParseError(leftOver.length >= numArguments.lowerBound, referenceUDA.multiplicityError(0));
            enforce!ArgParseError(leftOver.length <= numArguments.upperBound, "Missing positional arguments.");

            auto hasReadsDb = (leftOver.length == numArguments.upperBound);
            handleArg!"refDb"(this, leftOver[0]);
            leftOver = leftOver[1 .. $];

            if (hasReadsDb)
            {
                handleArg!"readsDb"(this, leftOver[0]);
                leftOver = leftOver[1 .. $];
            }

            handleArg!"dbAlignmentFile"(this, leftOver[0]);
            leftOver = leftOver[1 .. $];

            foreach (member; __traits(allMembers, typeof(this)))
            {
                alias symbol = Alias!(__traits(getMember, this, member));
                alias argUDAs = getUDAs!(symbol, Argument);

                static if (
                    argUDAs.length > 0 &&
                    !member.among("refDb", "readsDb", "dbAlignmentFile")
                )
                {
                    handleArg!member(this, leftOver[0]);
                    leftOver = leftOver[1 .. $];
                }
            }

            return this;
        }
    }

    static if (command.among(
        TestingCommand.translocateGaps,
        TestingCommand.buildPartialAssembly,
        TestingCommand.findClosableGaps,
        TestingCommand.checkResults,
    ))
    {
        @Argument("<in:true-assembly>")
        @Help("the 'true' assembly in .dam format")
        @Validate!(validateDB!".dam")
        string trueAssemblyDb;
    }

    static if (command.among(
        TestingCommand.translocateGaps,
    ))
    {
        @Argument("<in:short-read-assembly>")
        @Help("short-read assembly in .dam format")
        @Validate!(validateDB!".dam")
        string shortReadAssemblyDb;
    }

    static if (command.among(
        TestingCommand.filterMask,
        DentistCommand.maskRepetitiveRegions,
        DentistCommand.chainLocalAlignments,
        DentistCommand.showMask,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
        DentistCommand.output,
        TestingCommand.checkResults,
    ))
    {
        @Argument("<in:reference>")
        @Help("reference assembly in .dam format")
        @Validate!(validateDB)
        string refDb;


        @property id_t numReferenceContigs() inout
        {
            static id_t _numReferenceContigs;

            if (_numReferenceContigs == 0)
                _numReferenceContigs = getNumContigs(refDb, null);

            return _numReferenceContigs;
        }
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
        DentistCommand.chainLocalAlignments,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        static if (command.among(
            DentistCommand.maskRepetitiveRegions,
            DentistCommand.chainLocalAlignments,
        ))
            enum argReadsMultiplicity = Multiplicity.optional;
        else
            enum argReadsMultiplicity = 1;

        @Argument("<in:reads>", argReadsMultiplicity)
        @Help("set of PacBio reads in .db/.dam format")
        @Validate!(validateReadsDb)
        string readsDb;

        static void validateReadsDb(string readsDb)
        {
            if (argReadsMultiplicity == 1 || readsDb !is null)
                validateDB(readsDb);
        }

        @property bool hasReadsDb() const pure nothrow
        {
            return readsDb !is null;
        }
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Argument("<in:result>")
        @Help("result assembly in .dam format")
        @Validate!(validateDB!".dam")
        string resultDb;
    }

    static if (command.among(
        TestingCommand.translocateGaps,
    ))
    {
        @Argument("<in:short-vs-true-alignment>")
        @Help(q"{
            locals alignments of the short-read assembly against the 'true'
            assembly in form of a .las file as produced by `daligner`
        }")
        @Validate!((value, options) => validateLasFile(value, options.trueAssemblyDb, options.shortReadAssemblyDb))
        string shortReadAssemblyAlignmentFile;
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
        DentistCommand.chainLocalAlignments,
    ))
    {
        @Argument("<in:alignment>")
        @Help("self-alignment of the reference assembly or reads vs. reference alignment")
        @Validate!((value, options) => options.readsDb is null
            ? validateLasFile(value, options.refDb)
            : validateLasFile(value, options.refDb, options.readsDb))
        string dbAlignmentFile;
    }

    static if (command.among(
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @Argument("<in:ref-vs-reads-alignment>")
        @Help(q"{
            alignments chains of the reads against the reference in form of a .las
            file as produced by `damapper`
        }")
        @Validate!((value, options) => validateLasFile(value, options.refDb, options.readsDb))
        string readsAlignmentFile;
    }

    static if (command.among(
        DentistCommand.showPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @Argument("<in:pile-ups>")
        @Help("read pile ups from <pile-ups>")
        @Validate!validateFileExists
        string pileUpsFile;
    }

    static if (command.among(
        TestingCommand.filterMask,
    ))
    {
        @Argument("<in:input-mask>")
        @Help("filter Dazzler mask <input-mask>")
        @Validate!((value, options) => validateInputMask(options.refDb, value))
        string inMask;
    }

    static if (command.among(
        DentistCommand.showMask,
    ))
    {
        @Argument("<in:repeat-mask>", Multiplicity.oneOrMore)
        @Help("read Dazzler mask <repeat-mask>")
        @Validate!((values, options) => validateInputMasks(options.refDb, values))
        string[] masks;
    }

    static if (command.among(
        TestingCommand.findClosableGaps,
        TestingCommand.buildPartialAssembly,
        TestingCommand.checkResults,
    ))
    {
        @Argument("<in:mapped-regions-mask>")
        @Help("read regions that were kept aka. output contigs from the Dazzler mask.")
        @Validate!((value, options) => validateInputMask(options.trueAssemblyDb, value))
        string mappedRegionsMask;
    }


    static if (command.among(
        TestingCommand.findClosableGaps,
    ))
    {
        @Argument("<in:reads-map>")
        @Help(q"{
            true alignment of the reads aka. reads map; this is produced by the
            `-M` option of the `simulator` utility of the Dazzler tools.
        }")
        @Validate!validateFileExists
        string readsMap;
    }

    static if (command.among(
        DentistCommand.showInsertions,
        DentistCommand.output,
    ))
    {
        @Argument("<in:insertions>")
        @Help("read insertion information from <insertions> generated by the `merge-insertions` command")
        @Validate!validateFileExists
        string insertionsFile;
    }

    static if (command.among(
        DentistCommand.translateCoords,
    ))
    {
        @Argument("<in:scaffolding>")
        @Help(q"{
            read the assembly graph from <scaffolding> generate
            (see `--scaffolding` of the `output` command)
        }")
        @Validate!validateFileExists
        string assemblyGraphFile;
    }

    static if (command.among(
        DentistCommand.translateCoords,
    ))
    {
        @Argument("<coord-string>", Multiplicity.oneOrMore)
        @Help(q"{
            translate coordinate(s) given by <coord-string> of the result into
            coordinates on the reference. Coordinates are always 1-based.
            A <coord-string> the format `scaffold/<uint:scaffold-id>/<uint:coord>`
            which describes a coordinate on `>scaffold-<scaffold-id>` starting
            a the first base pair of the scaffold
        }")
        @Validate!validateCoordStrings
        string[] coordStrings;

        OutputCoordinate[] outputCoordinates;

        @PostValidate()
        void hookParseCoordStrings()
        {
            outputCoordinates.length = coordStrings.length;
            foreach (i, coordString; coordStrings)
                outputCoordinates[i] = parseCoordString(coordString);
        }
    }

    static if (command.among(
        TestingCommand.translocateGaps,
    ))
    {
        @Argument("<out:mapped-regions-mask>")
        @Help("write regions that were kept aka. output contigs into a Dazzler mask.")
        @Validate!((value, options) => validateOutputMask(options.trueAssemblyDb, value))
        string mappedRegionsMask;
    }

    static if (command.among(
        TestingCommand.filterMask,
    ))
    {
        @Argument("<out:filtered-mask>")
        @Help("write filtered Dazzler mask to <filtered-mask>")
        @Validate!((value, options) => validateOutputMask(options.refDb, value))
        string outMask;
    }

    static if (command.among(
        DentistCommand.collectPileUps,
    ))
    {
        @Argument("<out:pile-ups>")
        @Help("write inferred pile ups into <pile-ups>")
        @Validate!validateFileWritable
        string pileUpsFile;
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Argument("<out:repeat-mask>")
        @Help("write inferred repeat mask into a Dazzler mask.")
        @Validate!((value, options) => validateOutputMask(options.refDb, value))
        string repeatMask;
    }

    static if (command.among(
        DentistCommand.chainLocalAlignments,
    ))
    {
        @Argument("<out:chained-las>")
        @Help("write alignment chains to <chained-las>.")
        @Validate!(value => validateFileWritable(value))
        string chainedAlignments;
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Argument("<out:insertions>")
        @Help("write insertion information into <insertions>")
        @Validate!validateFileWritable
        string insertionsFile;
    }

    static if (command.among(
        DentistCommand.mergeInsertions,
    ))
    {
        @Argument("<out:insertions>")
        @Help("write merged insertion information to <insertions>")
        @Validate!validateFileWritable
        string mergedInsertionsFile;

        @Argument("<in:partitioned-insertions>", Multiplicity.oneOrMore)
        @Help("merge insertion information from <partitioned-insertions>... generated by the `processPileUps` command")
        @Validate!validateFilesExist
        @Validate!(value => value.length >= 2)
        string[] insertionsFiles;
    }

    static if (command.among(
        TestingCommand.buildPartialAssembly,
    ))
    {
        @Argument("<out:test-assembly>", Multiplicity.optional)
        @Help("write output assembly to <test-assembly> (default: stdout)")
        @Validate!(value => (value is null).execUnless!(() => validateFileWritable(value)))
        string resultFile;
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        @Argument("<out:result>", Multiplicity.optional)
        @Help("write gap-closed assembly to <result> (default: stdout)")
        @Validate!(value => (value is null).execUnless!(() => validateFileWritable(value)))
        string resultFile;
    }

    mixin HelpOption;

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("batch", "b")
        @MetaVar("<from>..<to>")
        @Help(q"{
            process only a subset of the gaps in the given range (excluding <to>).
            <from> and <to> are zero-based indices for the contigs of the
            reference assembly or <to> may be `$` to indicate the end of the
            reference.
        }")
        void parseReferenceContigBatch(string batchString) pure
        {
            try
            {
                if (batchString.endsWith("$"))
                {
                    batchString.formattedRead!"%d..$"(referenceContigBatch[0]);
                    referenceContigBatch[1] = id_t.max;
                }
                else
                {
                    batchString.formattedRead!"%d..%d"(referenceContigBatch[0], referenceContigBatch[1]);
                }
            }
            catch (Exception e)
            {
                throw new CLIException("ill-formatted batch range");
            }
        }

        @Option()
        @Validate!validateReferenceContigBatchRange
        id_t[2] referenceContigBatch;

        static void validateReferenceContigBatchRange(id_t[2] referenceContigBatch, OptionsFor!command options)
        {
            auto from = referenceContigBatch[0];
            auto to = referenceContigBatch[1];

            enforce!CLIException(
                referenceContigBatch == referenceContigBatch.init ||
                (0 <= from && from < to && (to == id_t.max || to <= options.numReferenceContigs)),
                format!"invalid batch range; check that 0 <= <from> < <to> <= %d"(
                        options.numReferenceContigs)
            );
        }

        @PostValidate()
        void hookEnsurePresenceOfBatchRange()
        {
            if (referenceContigBatch == referenceContigBatch.init || referenceContigBatch[1] == id_t.max)
            {
                referenceContigBatch[1] = numReferenceContigs;
            }
        }

        @property id_t referenceContigBatchSize() const pure nothrow
        {
            return referenceContigBatch[1] - referenceContigBatch[0];
        }
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        enum agpVersion = "2.1";

        @Option("agp")
        @Help(format!"write AGP v%s file that describes the output assembly"(agpVersion))
        string agpFile;
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("allow-single-reads")
        @Help("allow using single reads instead of consensus sequence for gap closing")
        OptionFlag allowSingleReads;
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("bad-fraction")
        @MetaVar("<frac>")
        @Help("
            Intrinsic QVs are categorized as \"bad\" if they are greater or equal to the best QV
            of the worst <frac> trace point intervals.
        ")
        @(Validate!(value => enforce!CLIException(
            0.0 <= value && value < 0.5,
            "--bad-fraction must be within [1, 0.5)")
        ))
        double badFraction = 0.08;
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("batch", "b")
        @MetaVar("<idx-spec>[,<idx-spec>...]")
        @Help(q"{
            process only a subset of the pile ups. <pile-up-ids> is a
            comma-separated list of <idx-spec>. Each
            <id-specifications> is either a single integer <idx> or a range
            <from>..<to>. <idx>, <from> and <to> are zero-based indices into
            the pile up DB. The range is right-open, i.e. index <to> is
            excluded. <to> may be a dollar-sign (`$`) to indicate the end of
            the pile up DB.
        }")
        void parsePileUpBatch(string batchString) pure
        {
            foreach (idxSpec; batchString.split(","))
                try
                    pileUpBatches ~= parsePileUpIdxSpec(idxSpec);
                catch (Exception e)
                    throw new CLIException("ill-formatted batch range");
        }

        static id_t[2] parsePileUpIdxSpec(string idxSpec) pure
        {
            id_t[2] batch;

            if (!idxSpec.canFind('.'))
            {
                idxSpec.formattedRead!"%d"(batch[0]);
                batch[1] = batch[0] + 1;
            }
            if (idxSpec.endsWith("$"))
            {
                idxSpec.formattedRead!"%d..$"(batch[0]);
                batch[1] = id_t.max;
            }
            else
            {
                idxSpec.formattedRead!"%d..%d"(batch[0], batch[1]);
            }

            return batch;
        }

        @property id_t pileUpLength() inout
        {
            static id_t numPileUps;

            if (numPileUps == 0)
            {
                numPileUps = PileUpDb.parse(pileUpsFile).length.to!id_t;
            }

            return numPileUps;
        }


        @Option()
        @Validate!validatePileUpBatches
        id_t[2][] pileUpBatches;

        static void validatePileUpBatches(id_t[2][] pileUpBatches, OptionsFor!command options)
        {
            foreach (pileUpBatch; pileUpBatches)
                validatePileUpBatchRange(pileUpBatch, options);

            validatePileUpBatchesDontIntersect(pileUpBatches, options);
        }

        static void validatePileUpBatchRange(id_t[2] pileUpBatch, OptionsFor!command options)
        {
            auto from = pileUpBatch[0];
            auto to = pileUpBatch[1];

            enforce!CLIException(
                pileUpBatch == pileUpBatch.init ||
                (0 <= from && from < to && (to == id_t.max || to <= options.pileUpLength)),
                format!"invalid batch range; check that 0 <= <from> < <to> <= %d"(
                        options.pileUpLength)
            );
        }

        static void validatePileUpBatchesDontIntersect(id_t[2][] pileUpBatches, OptionsFor!command options)
        {
            static bool intersect(id_t[2] batch1, id_t[2] batch2)
            {
                if (batch1[0] > batch2[0])
                    swap(batch1, batch2);

                return batch1[1] > batch2[0];
            }

            foreach (i, batch1; pileUpBatches)
                foreach (j, batch2; pileUpBatches[i + 1 .. $])
                    enforce!CLIException(
                        !intersect(batch1, batch2),
                        format!"invalid --batch: <idx-spec>'s at indices %d and %d intersect"(i, j),
                    );
        }

        @PostValidate()
        void hookEnsurePresenceOfBatchRanges()
        {
            foreach (ref pileUpBatch; pileUpBatches)
                if (pileUpBatch == pileUpBatch.init || pileUpBatch[1] == id_t.max)
                    pileUpBatch[1] = pileUpLength;
        }

        @PostValidate()
        void hookOptimizeBatchRanges()
        {
            pileUpBatches.sort!"a[0] < b[0] || (a[0] == b[0] && a[1] < b[1])";

            size_t accIdx;
            foreach (ref pileUpBatch; pileUpBatches[1 .. $])
            {
                if (pileUpBatch[0] <= pileUpBatches[accIdx][1])
                    pileUpBatches[accIdx][1] = max(pileUpBatch[1], pileUpBatches[accIdx][1]);
                else
                    pileUpBatches[++accIdx] = pileUpBatch;
            }

            pileUpBatches.length = accIdx + 1;
        }

        @property id_t numPileUps() const pure nothrow
        {
            return pileUpBatches
                .map!(pileUpBatch => pileUpBatch[1] - pileUpBatch[0])
                .sum;
        }
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        @Option("skip-gaps")
        @MetaVar("<gap-spec>[,<gap-spec>...]")
        @Help(q"{
            Do not close the specified gaps. Each <gap-spec> is a pair
            of contig IDs <contigA>-<contigA> meaning that the specified
            contigs should not be closed. They will still be joined by a
            prexisting gap.
        }")
        void parseSkipGaps(string skipGapsString) pure
        {
            foreach (gapSpec; skipGapsString.split(","))
                try
                    skipGaps ~= parseGapSpec(gapSpec);
                catch (Exception e)
                    throw new CLIException("ill-formatted <gap-spec>; should be <contigA>-<contigA>");
        }

        static id_t[2] parseGapSpec(string gapSpec) pure
        {
            id_t[2] gap;

            gapSpec.formattedRead!"%d-%d"(gap[0], gap[1]);

            return gap;
        }

        @Option()
        @(Validate!validateSkipGaps)
        id_t[2][] skipGaps;

        static void validateSkipGaps(id_t[2][] skipGaps, OptionsFor!command options)
        {
            foreach (skipGap; skipGaps)
                validatePileUpSkipGap(skipGap, options);
        }

        static void validatePileUpSkipGap(id_t[2] skipGap, OptionsFor!command options)
        {
            auto numReferenceContigs = options.numReferenceContigs;
            auto contigA = skipGap[0];
            auto contigB = skipGap[1];

            enforce!CLIException(
                0 < contigA && contigA <= numReferenceContigs,
                format!"invalid <gap-spec>: <contigA> == %d is out of bounds [1, %d]"(contigA, numReferenceContigs),
            );
            enforce!CLIException(
                0 < contigB && contigB <= numReferenceContigs,
                format!"invalid <gap-spec>: <contigB> == %d is out of bounds [1, %d]"(contigB, numReferenceContigs),
            );
            enforce!CLIException(
                contigA != contigB,
                "invalid <gap-spec>: <contigA> mut not be equal to <contigB>",
            );
        }

        @PostValidate()
        void hookSortSkipGaps()
        {
            foreach (ref skipGap; skipGaps)
                if (skipGap[0] > skipGap[1])
                    swap(skipGap[0], skipGap[1]);

            sort(skipGaps);
        }
    }

    static if (command.among(
        DentistCommand.collectPileUps,
    ))
    {
        @Option("best-pile-up-margin")
        @Help(format!q"{
            given a set of of conflicting gap closing candidates, if the
            largest has <double> times more reads than the second largest it
            is considered unique. If a candidates would close gap in the
            reference assembly marked by `n`s the number reads is multipled by
            --existing-gap-bonus. (default: %s)
        }"(defaultValue!bestPileUpMargin))
        @Validate!(value => enforce!CLIException(value > 1.0, "--best-pile-up-margin must be greater than 1.0"))
        double bestPileUpMargin = 3.0;
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("bucket-size", "B")
        @Help(format!q"{
            bucket size of the gap length histogram; use 0 to disable (default: %d)
        }"(defaultValue!bucketSize))
        coord_t bucketSize = 500;
    }

    static if (command.among(
        TestingCommand.checkResults,
        testingOnly!(DentistCommand.output),
    ))
    {
        @Option("cache-contig-alignments")
        @Help(format!(command == TestingCommand.checkResults
            ? "if given results of contig alignments will be cached as JSON
               (default: %s)"
            : "if given the contig location will be cached as JSON faking the
               effect of the same option in `check-results`. NOTE: the result
               has to amended manually to be fully valid. (default: %s)"
        )(defaultValue!contigAlignmentsCache))
        string contigAlignmentsCache;
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("cache-only")
        @Help("stop execution after writing the contig alignments cache")
        @Validate!((value, options) => enforce!CLIException(
            !value || options.contigAlignmentsCache !is null,
            "requires --cache-contig-alignments",
        ))
        OptionFlag cacheOnly;
    }

    static if (command == DentistCommand.validateConfig)
    {
        @Argument("<in:config>")
        @Help(configHelpString)
        @Validate!validateFileExists
        string configFile;
    }
    else
    {
        @Option("config")
        @MetaVar("<config-json>")
        @Help(configHelpString)
        @Validate!(value => (value is null).execUnless!(() => validateFileExists(value)))
        string configFile;

        @PreValidate(Priority.high)
        void retroMergeConfig()
        {
            if (configFile is null)
                return;

            this.retroInitFromConfig(configFile);
        }
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("crop-ambiguous")
        @MetaVar("<num>")
        @Help(format!q"{
            crop <num> bp from both ends of each reference contig when searching for exact copies
            in the reference itself in order to identify ambiguous contigs. If comparing different
            gap closing tools use the same value for all tools. (default: %d)
        }"(defaultValue!cropAmbiguous))
        // This is the amount PBJelly may modify
        coord_t cropAmbiguous = 100;
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("crop-alignment")
        @MetaVar("<num>")
        @Help(format!q"{
            crop <num> bp from both ends of each reference contig when searching for exact copies
            in the output of the gap closer. Keep this as low as possible, ie. if the contigs are
            not modified use zero. (default: %d)
        }"(defaultValue!cropAlignment))
        @Validate!((value, options) => enforce!CLIException(value <= options.cropAmbiguous,
                                                            "must be <= --crop-ambiguous"))
        coord_t cropAlignment = 0;
    }

    enum configHelpString = "
        provide configuration values in a JSON file. See README.md for
        usage and examples.
    ";

    static if (command.among(
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
        TestingCommand.checkResults,
    ))
    {
        @Option("auxiliary-threads", "aux-threads")
        @MetaVar("num-threads")
        @Help("
            use <num-threads> threads for auxiliary tools like `daligner`,
            `damapper` and `daccord`
            (defaults to floor(totalCpus / <threads>) )
        ")
        uint numAuxiliaryThreads;

        @PostValidate(Priority.low)
        void hookInitDaccordThreads()
        {
            if (numAuxiliaryThreads == 0)
            {
                numAuxiliaryThreads = totalCPUs / numThreads;
            }
        }
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("dust", "d")
        @MetaVar("<string>")
        @Help("
            Dazzler mask for low complexity regions. Uses " ~ defaultValue!dustMask ~ "
            by default but only if present.
        ")
        @Validate!((value, options) => value == defaultValue!dustMask || validateInputMask(options.trueAssemblyDb, value))
        string dustMask = "dust";

        @PostValidate(Priority.medium)
        void fixDefaultDustMask()
        {
            if (dustMask != defaultValue!dustMask)
                return;

            try
            {
                validateInputMask(trueAssemblyDb, dustMask);
            }
            catch(CLIException e)
            {
                dustMask = null;
            }
        }
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        @Option("scaffolding")
        @MetaVar("<insertions-db>")
        @Help("write the assembly scaffold to <insertions-db>; use `show-insertions` to inspect the result")
        @Validate!(value => (value is null).execUnless!(() => validateFileWritable(value)))
        string assemblyGraphFile;
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("gap-details")
        @MetaVar("<json>")
        @Help("write the summary for all gaps to a JSON file <json>")
        @Validate!(value => (value is null).execUnless!(() => validateFileWritable(value)))
        string gapDetailsJson;
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("gap-details-context")
        @MetaVar("<num>")
        @Help("add <num> base pairs of context to inserted sequence in gap details")
        coord_t gapDetailsContext;
    }

    static if (command.among(
        DentistCommand.collectPileUps,
    ))
    {
        @Option("debug-pile-ups")
        @MetaVar("<db-stem>")
        @Help("write pile ups of intermediate steps to `<db-stem>.<state>.db`")
        @Validate!(value => (value is null).execUnless!(() => validateFileWritable(value)))
        string intermediatePileUpsStem;
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Option("debug-repeat-masks")
        @Help("(only for reads-mask) write mask components into additional masks `<repeat-mask>-<component-type>`")
        OptionFlag debugRepeatMasks;
    }

    static if (command.among(
        TestingCommand.buildPartialAssembly,
        DentistCommand.output,
    ))
    {
        @Option("fasta-line-width", "w")
        @Help(format!"line width for ouput FASTA (default: %d)"(defaultValue!fastaLineWidth))
        @Validate!(value => enforce!CLIException(value > 0, "fasta line width must be greater than zero"))
        size_t fastaLineWidth = 50;
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        @Option("join-policy")
        @Help(format!q"{
            allow only joins (gap filling) in the given mode:
            `scaffoldGaps` (only join gaps inside of scaffolds –
            marked by `n`s in FASTA),
            `scaffolds` (join gaps inside of scaffolds and try to join scaffolds),
            `contigs` (break input into contigs and re-scaffold everything;
            maintains scaffold gaps where new scaffolds are consistent)
            (default: `%s`)
        }"(defaultValue!joinPolicy))
        JoinPolicy joinPolicy = JoinPolicy.scaffoldGaps;
    }

    static if (command.among(
        DentistCommand.showMask,
        DentistCommand.showPileUps,
        DentistCommand.showInsertions,
        DentistCommand.translateCoords,
        TestingCommand.checkResults,
    ))
    {
        @Option("json", "j")
        @Help("if given write the information in JSON format")
        OptionFlag useJson;
    }

    static if (needWorkdir)
    {
        @Option("keep-temp", "k")
        @Help("keep the temporary files; outputs the exact location")
        OptionFlag keepTemp;
    }

    static if (command.among(
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @Option("mask", "m")
        @MetaVar("<string>...")
        @Help("
            Dazzler masks for repetitive regions (at least one required;
            generate with `mask-repetitive-regions` command)
        ")
        void addMask(string mask) pure
        {
            repeatMasks ~= mask;
        }

        @Option()
        @Validate!((values, options) => validateInputMasks(options.refDb, values))
        @Validate!(value => value.length > 0)
        string[] repeatMasks;
    }

    static if (command.among(
        DentistCommand.collectPileUps,
    ))
    {
        enum numBubblesEscapeNodes = 2;
        enum nodePerContig = 2;
        enum numIntermediateContigs = 3;

        /// Consider cyclic subgraphs of up to the size when detecting
        /// _bubbles_ aka. _skipping pile ups.
        @Option()
        ubyte maxBubbleSize = numBubblesEscapeNodes + nodePerContig * numIntermediateContigs;

        /// Run the solver at most this number of times
        @Option()
        ubyte maxBubbleResolverIterations = 1 + numIntermediateContigs;
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Option("max-coverage-reads")
        @MetaVar("<uint>")
        @Help(q"{
            this is used to derive a repeat mask from the ref vs. reads alignment;
            if the alignment coverage is larger than <uint> it will be
            considered repetitive; a default value is derived from --read-coverage;
            both options are mutually exclusive
        }")
        id_t maxCoverageReads;

        @Option()
        id_t[2] coverageBoundsReads;

        @PostValidate(Priority.medium)
        void setCoverageBoundsReads()
        {
            if (!hasReadsDb)
                return;

            enforce!CLIException(
                maxCoverageReads != maxCoverageReads.init ||
                readCoverage != readCoverage.init,
                "must provide either --read-coverage or --max-coverage-reads",
            );
            enforce!CLIException(
                (maxCoverageReads != maxCoverageReads.init) ^
                (readCoverage != readCoverage.init),
                "must not provide both --read-coverage and --max-coverage-reads",
            );

            id_t upperBound(double x)
            {
                enum aReads = 1.65;
                enum bReads = 0.1650612;
                enum cReads = 5.9354533;

                return to!id_t(x / log_e(log_e(log_e(bReads * x + cReads)/log_e(aReads))));
            }

            if (readCoverage != readCoverage.init)
                maxCoverageReads = upperBound(readCoverage);

            coverageBoundsReads = [0, maxCoverageReads];
        }
    }


    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Option("max-improper-coverage-reads")
        @MetaVar("<uint>")
        @Help("
            this is used to derive a repeat mask from the ref vs. reads alignment;
            if the coverage of improper alignments is larger than <uint> it will be
            considered repetitive; a default value is derived from --read-coverage;
            both options are mutually exclusive
        ")
        id_t maxImproperCoverageReads;

        @Option()
        id_t[2] improperCoverageBoundsReads;

        @PostValidate(Priority.medium)
        void setImproperCoverageBoundsReads()
        {
            if (!hasReadsDb)
                return;

            enforce!CLIException(
                maxImproperCoverageReads != maxImproperCoverageReads.init ||
                readCoverage != readCoverage.init,
                "must provide either --read-coverage or --max-improper-coverage-reads",
            );
            enforce!CLIException(
                (maxImproperCoverageReads != maxImproperCoverageReads.init) ^
                (readCoverage != readCoverage.init),
                "must not provide both --read-coverage and --max-improper-coverage-reads",
            );

            id_t upperBound(double x)
            {
                enum a = 0.5;
                enum b = 0.1875;
                enum c = 8.0;

                // This is a smooth version of `max(4, a*x)`
                return to!id_t(a*x + exp(b*(c - x)));
            }

            if (readCoverage != readCoverage.init)
                maxImproperCoverageReads = upperBound(readCoverage);

            improperCoverageBoundsReads = [0, maxImproperCoverageReads];
        }
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Option("max-coverage-self")
        @MetaVar("<uint>")
        @Help(format!q"{
            this is used to derive a repeat mask from the self alignment;
            if the alignment coverage larger than <uint> it will be
            considered repetitive (default: %d)
        }"(defaultValue!maxCoverageSelf))
        @Validate!(validatePositive!("max-coverage-self", id_t))
        id_t maxCoverageSelf = 4;

        @Option()
        id_t[2] coverageBoundsSelf;

        @PostValidate(Priority.medium)
        void setCoverageBoundsSelf()
        {
            if (hasReadsDb)
                return;

            coverageBoundsSelf = [0, maxCoverageSelf];
        }
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        @Option("max-insertion-error")
        @Help(format!"
            insertion and existing contigs must match with less error than <double> (default: %s)
        "(defaultValue!maxInsertionError))
        @Validate!(value => enforce!CLIException(
            0.0 < value && value <= 0.3,
            "maximum insertion error rate must be in (0, 0.3]"
        ))
        double maxInsertionError = 1e-2;
    }

    static if (command.among(
        TestingCommand.findClosableGaps,
        DentistCommand.generateDazzlerOptions,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @Option("min-anchor-length")
        @Help(format!q"{
            alignment need to have at least this length of unique anchoring sequence (default: %d)
        }"(defaultValue!minAnchorLength))
        @Validate!(value => enforce!CLIException(value > 0, "minimum anchor length must be greater than zero"))
        @Validate!(
            (value, options) => enforce!CLIException(
                value > options.tracePointDistance,
                "minimum anchor length should be greater than trace point spacing of *.las file"
            ),
            is(typeof(OptionsFor!command().tracePointDistance)),
        )
        coord_t minAnchorLength = 500;
    }

    static if (command.among(
        TestingCommand.filterMask,
    ))
    {
        @Option("min-interval-size")
        @Help(format!"
            minimum size for mask intervals (default: %d)
        "(defaultValue!minIntervalSize))
        coord_t minIntervalSize;
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        @Option("min-extension-length")
        @Help(format!q"{
            extensions must have at least <ulong> bps of consensus to be inserted (default: %d)
        }"(defaultValue!minExtensionLength))
        @Validate!(value => enforce!CLIException(value > 0, "minimum extension length must be greater than zero"))
        size_t minExtensionLength = 100;
    }

    static if (command.among(
        TestingCommand.filterMask,
    ))
    {
        @Option("min-gap-size")
        @Help(format!"
            minimum size for gaps between mask intervals (default: %d)
        "(defaultValue!minGapSize))
        coord_t minGapSize;
    }

    static enum defaultMinSpanningReads = 3;

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("min-reads-per-pile-up")
        @Help(format!q"{
            pile ups must have at least <ulong> reads to be processed (default: %d)
        }"(defaultValue!minReadsPerPileUp))
        @Validate!(value => enforce!CLIException(value > 0, "min reads per pile up must be greater than zero"))
        size_t minReadsPerPileUp = defaultMinSpanningReads;
    }

    static if (command.among(
        TestingCommand.findClosableGaps,
        DentistCommand.collectPileUps,
    ))
    {
        @Option("min-spanning-reads", "s")
        @Help(format!q"{
            require at least <ulong> spanning reads to close a gap (default: %d)
        }"(defaultValue!minSpanningReads))
        size_t minSpanningReads = defaultMinSpanningReads;
    }

    static if (command.among(
        DentistCommand.collectPileUps,
    ))
    {
        @Option("no-merge-extension")
        @Help("Do not merge extension reads into spanning pile ups.")
        OptionFlag noMergeExtensions;

        @property bool mergeExtensions() const pure nothrow
        {
            return !noMergeExtensions;
        }
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        @Option("no-highlight-insertions", "H")
        @Help("
            turn off highlighting (upper case) of inserted sequences in the
            FASTA output
        ")
        OptionFlag noHighlightInsertions;
    }

    static if (command.among(
        DentistCommand.processPileUps,
        DentistCommand.output,
    ))
    {
        @Option("only")
        @Help(format!"
            only process/output insertions of the given type. Note, extending
            insertions are experimental and may produce invalid results.
            (default: spanning)
        ")
        OnlyFlag onlyFlag;

        enum OnlyFlag
        {
            spanning = 1 << 0,
            extending = 1 << 1,
            both = spanning | extending,
        }

        alias OnlyFlags = BitFlags!(OnlyFlag, Yes.unsafe);

        OnlyFlags onlyFlags;

        @PostValidate()
        void initOnlyFlags() pure nothrow
        {
            onlyFlags = OnlyFlags(onlyFlag);
        }
    }

    static if (command.among(
        DentistCommand.chainLocalAlignments,
    ))
    {
        @Option("progress")
        @Help("Print regular status reports on the progress.")
        OptionFlag printProgress;

        @Option("progress-every")
        @MetaVar("<msecs>")
        @Help(format!"
            Print status reports every <msecs>. (default: %d)
        "(defaultValue!printProgressEvery))
        uint printProgressEvery = 500;

        @Option("progress-format")
        @MetaVar("<format>")
        @Help(format!"
            Use <format> for status report lines where <format> is either
            `human` or `json`. The former prints a status line that updates
            regularly while the latter prints a full JSON record per line with
            every update (default: %s)
        "(defaultValue!progressFormat))
        ProgressMeter.Format progressFormat = ProgressMeter.Format.human;


        ProgressMeter createProgressMeter() const @safe
        {
            import std.typecons : Flag;

            ProgressMeter progress;

            progress.silent = cast(Flag!"silent") !printProgress;
            progress.unit = ProgressMeter.Unit.auto_;
            progress.printEveryMsecs = printProgressEvery;
            progress.format = progressFormat;

            return progress;
        }
    }

    static if (command.among(
        TestingCommand.translocateGaps,
        DentistCommand.maskRepetitiveRegions,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @Option("proper-alignment-allowance")
        @MetaVar("num")
        @Help("
            An alignment is called proper if it is end-to-end with at most <num> bp allowance.
            (default: --trace-point-spacing)
        ")
        coord_t properAlignmentAllowance;

        @PostValidate(Priority.low)
        void hookEnsurePresenceOfProperAlignmentAllowance()
        {
            if (properAlignmentAllowance > 0)
                return;

            static if (command == TestingCommand.translocateGaps)
                properAlignmentAllowance = 1_000;
            else
                properAlignmentAllowance = tracePointDistance;
        }
    }


    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Option("read-coverage", "C")
        @Help(q"{
            this is used to provide good default values for --max-coverage-reads;
            both options are mutually exclusive
        }")
        double readCoverage;
    }

    static if (command.among(
        DentistCommand.generateDazzlerOptions,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @Option("reads-error")
        @Help(format!"
            estimated error rate in reads (default: %s)
        "(defaultValue!readsErrorRate))
        @Validate!(value => enforce!CLIException(
            0.0 < value && value <= 0.3,
            "reads error rate must be in (0, 0.3]"
        ))
        double readsErrorRate = .15;
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("recover-imperfect-contigs", "R")
        @Help("try to recover imperfect contigs")
        OptionFlag recoverImperfectContigs;

        @Option()
        double maxImperfectContigError = 0.015;

        @property string[] recoverImperfectContigsAlignmentOptions() const
        {
            return [
                DamapperOptions.symmetric,
                DamapperOptions.oneDirection,
                DamapperOptions.numThreads ~ numAuxiliaryThreads.to!string,
                DamapperOptions.kMerSize ~ "32",
                format!(DamapperOptions.averageCorrelationRate ~ "%f")(
                    1.0 - maxImperfectContigError,
                ),
            ];
        }
    }

    static if (command.among(
        DentistCommand.generateDazzlerOptions,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @Option("reference-error")
        @Help(format!"
            estimated error rate in reference (default: %s)
        "(defaultValue!referenceErrorRate))
        @Validate!(value => enforce!CLIException(
            0.0 < value && value <= 0.3,
            "reference error rate must be in (0, 0.3]"
        ))
        double referenceErrorRate = .01;
    }

    static if (command.among(
        DentistCommand.collectPileUps,
    ))
    {
        @Option("existing-gap-bonus")
        @Help(format!q"{
            if a candidate would close an existing gap its size is multipled
            by <double> before conflict resolution
            (see --best-pile-up-margin). (default: %s)
        }"(defaultValue!existingGapBonus))
        @Validate!(value => enforce!CLIException(1.0 <= value, "--existing-gap-bonus must be at least 1.0"))
        double existingGapBonus = 6.0;
    }

    static if (command.among(
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
        TestingCommand.checkResults,
    ))
    {
        @Option("threads", "T")
        @Help("use <uint> threads (defaults to the number of cores)")
        uint numThreads;

        @PostValidate(Priority.high)
        void hookInitThreads()
        {
            if (numThreads > 0)
                defaultPoolThreads = numThreads - 1;

            numThreads = defaultPoolThreads + 1;
        }
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
        DentistCommand.chainLocalAlignments,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @(Validate!(validatePositive!("trace-point-distance", trace_point_t)))
        trace_point_t tracePointDistance;

        @PreValidate(Priority.medium)
        void hookGetTracePointDistance()
        {
            static if (is(typeof(readsDb)) && is(typeof(readsAlignmentFile)))
                auto alignmentFile = readsAlignmentFile;
            else static if (is(typeof(dbAlignmentFile)))
                auto alignmentFile = dbAlignmentFile;

            static if (is(typeof(alignmentFile)))
            {
                try
                {
                    validateLasFile(alignmentFile);
                    tracePointDistance = getTracePointDistance(alignmentFile);
                }
                catch (Exception e)
                {
                    // this will trigger the validation if need be
                    assert(tracePointDistance == 0);
                }
            }
        }
    }

    @Option("quiet", "q")
    @Help("
        reduce output as much as possible reporting only fatal errors. If
        given this option overrides --verbose.
    ")
    OptionFlag quiet;

    @Option("verbose", "v")
    @Help("
        increase output to help identify problems; use up to three times.
        Warning: performance may be drastically reduced if using three times.
    ")
    void increaseVerbosity() pure
    {
        ++verbosity;
    }

    @Option()
    @Validate!(value => enforce!CLIException(
        0 <= value && value <= 3,
        "verbosity must used 0-3 times"
    ))
    size_t verbosity = 0;

    @PostValidate(Priority.high)
    void hookInitLogLevel()
    {
        if (quiet)
            verbosity = 0;

        if (verbosity >= 3)
            logJsonWarn("info", "high level of verbosity may drastically reduce performance");

        switch (verbosity)
        {
        case 3:
            setLogLevel(LogLevel.debug_);
            break;
        case 2:
            setLogLevel(LogLevel.diagnostic);
            break;
        case 1:
            setLogLevel(LogLevel.info);
            break;
        case 0:
        default:
            setLogLevel(LogLevel.error);
            break;
        }
    }

    static if (needWorkdir)
    {
        /**
            Last part of the working directory name. A directory in the temp
            directory as returned by `std.file.tmpDir` with the naming scheme will
            be created to hold all data for the computation.
        */
        enum workdirTemplate = format!"dentist-%s-XXXXXX"(command);

        /// This is a temporary directory to store all working data.
        @Option("workdir", "P")
        @Help("use <string> as a working directory")
        string workdir;

        @PostValidate(Priority.high)
        void hookCreateWorkDir()
        {
            if (workdir is null)
            {
                auto workdirTemplate = buildPath(tempDir(), workdirTemplate);

                workdir = mkdtemp(workdirTemplate);
            }
            else
            {
                try
                {
                    enforce!CLIException(
                        isDir(workdir),
                        "--workdir is not a directory",
                    );

                    logJsonInfo(
                        "info", "using existing workdir",
                        "workdir", workdir,
                    );
                }
                catch (FileException e)
                {
                    // `isDir` failed; just ignore that
                }

                try
                {
                    mkdirRecurse(workdir);
                }
                catch(Exception e)
                {
                    throw new CLIException("could not create workdir: " ~ e.msg, e);
                }

                logJsonInfo(
                    "info", "recursively created workdir",
                    "workdir", workdir,
                );
            }
        }

        @CleanUp(Priority.low)
        void hookCleanWorkDir() const
        {
            if (keepTemp)
                return;

            try
            {
                rmdirRecurse(workdir);
            }
            catch (Exception e)
            {
                log(LogLevel.fatal, "Fatal: " ~ e.msg);
            }
        }
    }

    static if (is(typeof(OptionsFor!command().referenceErrorRate))) {
        @Validate!validateAverageCorrelationRate
        @property auto selfAlignmentOptionsAverageCorrelationRate() const
        {
            return (1 - referenceErrorRate)^^2;
        }

        @property string[] shortVsTrueAssemblyAlignmentOptions() const
        {
            return [
                DamapperOptions.symmetric,
                DamapperOptions.oneDirection,
                format!(DamapperOptions.averageCorrelationRate ~ "%f")(
                    selfAlignmentOptionsAverageCorrelationRate,
                ),
            ];
        }
    }

    static if (
        is(typeof(OptionsFor!command().minAnchorLength)) &&
        is(typeof(OptionsFor!command().referenceErrorRate))
    ) {
        @property string[] selfAlignmentOptions() const
        {
            return [
                DalignerOptions.identity,
                format!(DalignerOptions.minAlignmentLength ~ "%d")(minAnchorLength),
                format!(DalignerOptions.averageCorrelationRate ~ "%f")(
                    selfAlignmentOptionsAverageCorrelationRate,
                ),
            ];
        }

        @property string[] tandemAlignmentOptions() const
        {
            return [
                format!(DatanderOptions.minAlignmentLength ~ "%d")(minAnchorLength),
                format!(DatanderOptions.averageCorrelationRate ~ "%f")(
                    selfAlignmentOptionsAverageCorrelationRate,
                ),
            ];
        }
    }

    static if (
        is(typeof(OptionsFor!command().referenceErrorRate)) &&
        is(typeof(OptionsFor!command().readsErrorRate))
    ) {
        @Validate!validateAverageCorrelationRate
        @property auto refVsReadsAlignmentOptionsAverageCorrelationRate() const
        {
            return (1 - referenceErrorRate) * (1 - readsErrorRate);
        }

        @property string[] refVsReadsAlignmentOptions() const
        {
            return [
                DamapperOptions.symmetric,
                DamapperOptions.oneDirection,
                DamapperOptions.bestMatches ~ ".7",
                format!(DamapperOptions.averageCorrelationRate ~ "%f")(
                    refVsReadsAlignmentOptionsAverageCorrelationRate,
                ),
            ];
        }
    }

    static if (
        is(typeof(OptionsFor!command().readsErrorRate)) &&
        is(typeof(OptionsFor!command().numAuxiliaryThreads))
    ) {
        enum forceLargeTracePointType = 200;

        @Validate!validateAverageCorrelationRate
        @property auto pileUpAlignmentOptionsAverageCorrelationRate() const
        {
            return (1 - readsErrorRate)^^2;
        }

        @property string[] pileUpAlignmentOptions() const
        {
            return [
                DalignerOptions.numThreads ~ numAuxiliaryThreads.to!string,
                DalignerOptions.bridge,
                format!(DalignerOptions.tracePointDistance ~ "%d")(forceLargeTracePointType),
                format!(DalignerOptions.minAlignmentLength ~ "%d")(minAnchorLength),
                DalignerOptions.masks ~ "dust",
                format!(DalignerOptions.averageCorrelationRate ~ "%f")(
                    pileUpAlignmentOptionsAverageCorrelationRate,
                ),
            ];
        }

        @property string[] pileUpDustOptions() const
        {
            return [];
        }
    }

    static if (
        is(typeof(OptionsFor!command().numAuxiliaryThreads)) &&
        is(typeof(OptionsFor!command().tracePointDistance)) &&
        is(typeof(OptionsFor!command().workdir))
    ) {
        enum flankingContigsRepeatMaskName = "rep";
        enum postConsensusAlignmentOptionsAverageCorrelationRate = 0.7;

        @property string[] postConsensusAlignmentOptions() const
        {
            return [
                DalignerOptions.asymmetric,
                DalignerOptions.numThreads ~ numAuxiliaryThreads.to!string,
                DalignerOptions.masks ~ dbdustMaskName,
                DalignerOptions.masks ~ flankingContigsRepeatMaskName,
                format!(DalignerOptions.minAlignmentLength ~ "%d")(tracePointDistance),
                format!(DalignerOptions.averageCorrelationRate ~ "%f")(
                    postConsensusAlignmentOptionsAverageCorrelationRate,
                ),
            ];
        }
    }

    static if (command.among(
        DentistCommand.collectPileUps,
    ))
    {
        static struct AnchorSkippingPileUpsOptions
        {
            string[] damapperOptions;
            string[] dbsplitOptions;
            string workdir;
        }

        @property string[] intermediateContigsAlignmentOptions() const
        {
            return [
                DamapperOptions.symmetric,
                DamapperOptions.oneDirection,
                DamapperOptions.numThreads ~ numAuxiliaryThreads.to!string,
                format!(DamapperOptions.averageCorrelationRate ~ "%f")(
                    refVsReadsAlignmentOptionsAverageCorrelationRate,
                ),
            ];
        }

        @property auto anchorSkippingPileUpsOptions() const
        {
            return const(AnchorSkippingPileUpsOptions)(
                // dalignerOptions
                intermediateContigsAlignmentOptions,
                // dbsplitOptions
                [],
                // workdir
                workdir,
            );
        }
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        static struct ConsensusOptions
        {
            string[] daccordOptions;
            string[] dalignerOptions;
            string[] dbsplitOptions;
            string[] lasFilterAlignmentsOptions;
            string[] dbdustOptions;
            string workdir;
            coord_t properAlignmentAllowance;
        }

        @property auto consensusOptions() const
        {
            return const(ConsensusOptions)(
                // daccordOptions
                [
                    DaccordOptions.produceFullSequences,
                    DaccordOptions.numberOfThreads ~ numAuxiliaryThreads.to!string,
                ],
                // dalignerOptions
                pileUpAlignmentOptions,
                // dbsplitOptions
                [],
                // lasFilterAlignmentsOptions
                [
                    LasFilterAlignmentsOptions.errorThresold ~ (2.0 * readsErrorRate).to!string,
                ],
                // dbdustOptions
                [],
                // workdir
                workdir,
                // properAlignmentAllowance
                properAlignmentAllowance,
            );
        }
    }

    static auto defaultValue(alias property)(in uint precision = 2) pure nothrow
    {
        OptionsFor!command defaultOptions;

        auto value = __traits(getMember, defaultOptions, property.stringof);

        static if (isFloatingPoint!(typeof(value)))
            return value.toString(precision);
        else
            return value;
    }

    static auto numArguments() pure nothrow
    {
        alias ThisOptions = OptionsFor!command;
        size_t lowerBound;
        size_t upperBound;

        foreach (member; __traits(allMembers, ThisOptions))
        {
            alias symbol = Alias!(__traits(getMember, ThisOptions, member));
            alias argUDAs = getUDAs!(symbol, Argument);

            static if (argUDAs.length > 0)
            {
                lowerBound += argUDAs[0].lowerBound;
                upperBound += argUDAs[0].upperBound;
            }
        }

        return tuple!("lowerBound", "upperBound")(lowerBound, upperBound);
    }
}

private bool hasValidOptionNames(Options)() pure nothrow
{
    foreach (member; __traits(allMembers, Options))
    {
        alias symbol = Alias!(__traits(getMember, Options, member));
        alias argUDAs = getUDAs!(symbol, Argument);

        static if (argUDAs.length > 0)
        {
            enum argName = argUDAs[0].name;

            static assert(argName[0] == '<', "argument `" ~ argName ~ "` must start with `<`");
            static assert(argName[$ - 1] == '>', "argument `" ~ argName ~ "` must end with `>`");
        }
    }

    return true;
}

static foreach (command; EnumMembers!DentistCommand)
    static assert(hasValidOptionNames!(OptionsFor!command));

unittest
{
    static foreach (command; EnumMembers!DentistCommand)
    {
        static assert(is(OptionsFor!command));
    }
}

/// A short summary for each command to be output underneath the usage.
template commandSummary(DentistCommand command)
{
    static if (command == DentistCommand.validateConfig)
        enum commandSummary = q"{
            Validate config file. Exit with non-zero status if errors are
            found.
        }".wrap;
    else static if (command == TestingCommand.translocateGaps)
        enum commandSummary = q"{
            Translocate gaps from first assembly to second assembly.
        }".wrap;
    else static if (command == TestingCommand.filterMask)
        enum commandSummary = q"{
            Filter a Dazzler mask.
        }".wrap;
    else static if (command == TestingCommand.buildPartialAssembly)
        enum commandSummary = q"{
            Build a partial assembly from a mask.
        }".wrap;
    else static if (command == TestingCommand.findClosableGaps)
        enum commandSummary = q"{
            Find which gaps are closable, ie. the true alignment of the reads
            provides sufficient spanning reads.
        }".wrap;
    else static if (command == DentistCommand.generateDazzlerOptions)
        enum commandSummary = q"{
            Generate a set of options to pass to `daligner` and `damapper`
            needed for the input alignments.
        }".wrap;
    else static if (command == DentistCommand.maskRepetitiveRegions)
        enum commandSummary = q"{
            Mask regions that have a alignment coverage that is out of bounds.
        }".wrap;
    else static if (command == DentistCommand.showMask)
        enum commandSummary = q"{
            Show a short summary of the mask.
        }".wrap;
    else static if (command == DentistCommand.chainLocalAlignments)
        enum commandSummary = q"{
            Chain local alignment.
        }".wrap;
    else static if (command == DentistCommand.collectPileUps)
        enum commandSummary = q"{
            Build pile ups.
        }".wrap;
    else static if (command == DentistCommand.showPileUps)
        enum commandSummary = q"{
            Show a short summary of the pile ups.
        }".wrap;
    else static if (command == DentistCommand.processPileUps)
        enum commandSummary = q"{
            Process pile ups.
        }".wrap;
    else static if (command == DentistCommand.showInsertions)
        enum commandSummary = q"{
            Show a short summary of the insertions.
        }".wrap;
    else static if (command == DentistCommand.mergeInsertions)
        enum commandSummary = q"{
            Merge multiple insertions files into a single one.
        }".wrap;
    else static if (command == DentistCommand.output)
        enum commandSummary = q"{
            Write output.
        }".wrap;
    else static if (command == DentistCommand.translateCoords)
        enum commandSummary = q"{
            Translate coordinates of result assembly to coordinates of
            input assembly.
        }".wrap;
    else static if (command == TestingCommand.checkResults)
        enum commandSummary = q"{
            Check results of some gap closing procedure.
        }".wrap;
    else
        static assert(0, "missing commandSummary for " ~ command.to!string);
}

unittest
{
    static foreach (command; EnumMembers!DentistCommand)
    {
        static assert(is(typeof(commandSummary!command)));
    }
}

/// This describes the basic, ie. non-command-specific, options of `dentist`.
struct BaseOptions
{
    @Option("dependencies", "d")
    @Help("Print a list of external binaries that must be available on PATH.")
    OptionFlag listDependencies;

    mixin HelpOption;

    @Option("version")
    @Help("Print software version.")
    OptionFlag version_;

    @Argument("<command>")
    @Help(format!q"{
        Execute <command>. Available commands are: %-(%s, %). Use
        `dentist <command> --help` to get help for a specific command.
        <command> may be abbreviated by using a unique prefix of the full
        command string.
    }"([dentistCommands]))
    DentistCommand command;

    @Argument("<options...>", Multiplicity.optional)
    @Help("Command specific options")
    string commandOptions;
}

class CLIException : Exception
{
    ///
    mixin basicExceptionCtors;
}

private
{
    ReturnCode runCommand(DentistCommand command)(in string[] args)
    {
        alias Options = OptionsFor!command;
        enum usage = usageString!Options(executableName ~ " " ~ Options.commandName);

        Options options;

        try
        {
            assertExternalToolsAvailable();
            options = processOptions(parseArgs!Options(args[1 .. $]));
        }
        catch (ArgParseHelp e)
        {
            // Help was requested
            stderr.write(usage);
            stderr.writeln();
            stderr.writeln(commandSummary!command);
            stderr.writeln();
            stderr.write(helpString!Options);

            return ReturnCode.ok;
        }
        catch (UsageRequested e)
        {
            stderr.write(usage);

            return ReturnCode.ok;
        }
        catch (Exception e)
        {
            stderr.writeln("Error: " ~ (true || shouldLog(LogLevel.diagnostic)
                ? e.to!string
                : e.msg));
            stderr.writeln();
            stderr.write(usage);

            return ReturnCode.commandlineError;
        }

        const finalOptions = options;
        logInfo(finalOptions.serializeToJsonString());

        scope (exit) cast(void) cleanUp(finalOptions);

        try
        {
            mixin("import dentist.commands." ~ command.to!string ~ " : execute;");
            execute(finalOptions);

            return ReturnCode.ok;
        }
        catch (Exception e)
        {
            stderr.writeln("Error: " ~ (shouldLog(LogLevel.diagnostic)
                ? e.to!string
                : e.msg));

            return ReturnCode.runtimeError;
        }
    }

    enum getUDA(alias symbol, T) = getUDAs!(symbol, T)[0];

    struct Validate(alias _validate, bool isEnabled = true) {
        static if (isEnabled)
            alias validate = _validate;
        else
            alias validate = __truth;

        static bool __truth(T)(T) { return true; }
    }

    enum Priority
    {
        low,
        medium,
        high,
    }

    struct PreValidate {
        Priority priority;
    }

    struct PostValidate {
        Priority priority;
    }

    struct CleanUp {
        Priority priority;
    }

    template cmpPriority(T)
    {
        enum cmpPriority(alias a, alias b) = getUDA!(a, T).priority > getUDA!(b, T).priority;
    }

    unittest
    {
        struct Tester
        {
            @PostValidate(Priority.low)
            void priorityLow() { }

            @PostValidate(Priority.medium)
            void priorityMedium() { }

            @PostValidate(Priority.high)
            void priorityHigh() { }
        }

        alias compare = cmpPriority!PostValidate;

        static assert(compare!(
            Tester.priorityHigh,
            Tester.priorityLow,
        ));
        static assert(!compare!(
            Tester.priorityLow,
            Tester.priorityHigh,
        ));
        static assert(!compare!(
            Tester.priorityMedium,
            Tester.priorityMedium,
        ));
    }

    Options processOptions(Options)(Options options)
    {
        alias preValidateQueue = staticSort!(
            cmpPriority!PreValidate,
            getSymbolsByUDA!(Options, PreValidate),
        );

        static foreach (alias symbol; preValidateQueue)
        {
            mixin("options." ~ __traits(identifier, symbol) ~ "();");
        }

        static foreach (alias symbol; getSymbolsByUDA!(Options, Validate))
        {{
            alias validateUDAs = getUDAs!(symbol, Validate);

            foreach (validateUDA; validateUDAs)
            {
                alias validate = validateUDA.validate;
                auto value = __traits(getMember, options, __traits(identifier, symbol));
                alias Value = typeof(value);
                alias Validator = typeof(validate);

                try
                {
                    static if (is(typeof(validate(value))))
                        cast(void) validate(value);
                    else static if (is(typeof(validate(value, options))))
                        cast(void) validate(value, options);
                    else
                        static assert(0, format!q"{
                            validator for %s.%s should have a signature of
                            `void (T value);` or `void (T value, Options options);` -
                            maybe the validator does not compile?
                        }"(Options.stringof, symbol.stringof).wrap(size_t.max));
                }
                catch (Exception cause)
                {
                    enum isOption = getUDAs!(symbol, Option).length > 0;
                    enum isArgument = getUDAs!(symbol, Argument).length > 0;

                    static if (isOption)
                    {
                        enum thing = "option";
                        enum name = getUDAs!(symbol, Option)[0].toString();
                    }
                    else static if (isArgument)
                    {
                        enum thing = "argument";
                        enum name = getUDAs!(symbol, Argument)[0].name;
                    }
                    else
                    {
                        enum thing = "property";
                        enum name = __traits(identifier, symbol);
                    }

                    throw new CLIException("invalid " ~ thing ~ " " ~ name ~ ": " ~ cause.msg, cause);
                }
            }
        }}

        alias postValidateQueue = staticSort!(
            cmpPriority!PostValidate,
            getSymbolsByUDA!(Options, PostValidate),
        );

        static foreach (alias symbol; postValidateQueue)
        {
            mixin("options." ~ __traits(identifier, symbol) ~ "();");
        }

        return options;
    }

    unittest
    {
        import std.exception : assertThrown;

        struct Tester
        {
            @Validate!(value => enforce!Exception(value == 1))
            int a = 1;

            @Validate!((value, options) => enforce!Exception(value == 2 * options.a))
            int b = 2;

            string[] calls;

            @PostValidate(Priority.low)
            void priorityLow() {
                calls ~= "priorityLow";
            }

            @PostValidate(Priority.medium)
            void priorityMedium() {
                calls ~= "priorityMedium";
            }

            @PostValidate(Priority.high)
            void priorityHigh() {
                calls ~= "priorityHigh";
            }
        }

        Tester options;

        options = processOptions(options);

        assert(options.calls == [
            "priorityHigh",
            "priorityMedium",
            "priorityLow",
        ]);

        options.a = 2;

        assertThrown!Exception(processOptions(options));
    }

    Options cleanUp(Options)(Options options)
    {
        alias cleanUpQueue = staticSort!(
            cmpPriority!CleanUp,
            getSymbolsByUDA!(Options, CleanUp),
        );

        static foreach (alias symbol; cleanUpQueue)
        {
            mixin("options." ~ __traits(identifier, symbol) ~ "();");
        }

        return options;
    }

    unittest
    {
        import std.exception : assertThrown;

        struct Tester
        {
            string[] calls;

            @CleanUp(Priority.low)
            void priorityLow() {
                calls ~= "priorityLow";
            }

            @CleanUp(Priority.medium)
            void priorityMedium() {
                calls ~= "priorityMedium";
            }

            @CleanUp(Priority.high)
            void priorityHigh() {
                calls ~= "priorityHigh";
            }
        }

        Tester options;

        options = cleanUp(options);

        assert(options.calls == [
            "priorityHigh",
            "priorityMedium",
            "priorityLow",
        ]);
    }

    void parseRange(alias dest, string msg = "ill-formatted range")(in string rangeString) pure
            if (isStaticArray!(typeof(dest)) && dest.length == 2)
    {
        try
        {
            rangeString[].formattedRead!"%d..%d"(dest[0], dest[1]);
        }
        catch (Exception e)
        {
            throw new CLIException(msg);
        }
    }

    DestType parseRange(DestType, string msg = "ill-formatted range")(in string rangeString) pure
            if (isStaticArray!DestType && DestType.init.length == 2)
    {
        try
        {
            DestType dest;

            rangeString[].formattedRead!"%d..%d"(dest[0], dest[1]);

            return dest;
        }
        catch (Exception e)
        {
            throw new CLIException(msg);
        }
    }

    void validatePositive(string option, V)(V value)
    {
        enforce!CLIException(
            0 < value,
            option ~ " must be greater than zero",
        );
    }

    void validateAverageCorrelationRate(V)(V value)
    {
        enforce!CLIException(
            0.7 <= value && value < 1.0,
            "-e option of daligner/damapper must be in [0.7, 1) - some error rate(s) are too high",
        );
    }

    void validateCoverageBounds(DestType, string option)(in string coverageBoundsString)
    {
        auto coverageBounds = parseRange!DestType(coverageBoundsString);
        auto from = coverageBounds[0];
        auto to = coverageBounds[1];

        enforce!CLIException(
            coverageBounds == coverageBounds.init || 0 <= from && from < to,
            "invalid coverage bounds (--" ~ option ~ "); check that 0 <= <from> < <to>"
        );
    }

    void validateFilesExist(string msg = null)(in string[] files)
    {
        foreach (file; files)
        {
            static if (msg is null)
                validateFileExists(file);
            else
                validateFileExists!msg(file);
        }
    }

    void validateFileExists(string msg = "cannot open file `%s`")(in string file)
    {
        enforce!CLIException(file.exists, format!msg(file));
    }

    alias typeOf(alias T) = typeof(T);

    void validateFileExtension(
        string msg = "expected %-(%s or %) but got %s",
        extensions...
    )(in string file)
            if (allSatisfy!(isSomeString, staticMap!(typeOf, extensions)))
    {
        enum defaultMsg = "expected %-(%s or %) but got %s";

        enforce!CLIException(
            file.endsWith(extensions),
            format!(msg !is null ? msg : defaultMsg)([extensions], file)
        );
    }

    void validateDB(string extension = null)(in string dbFile)
        if (extension is null || extension.among(".dam", ".db"))
    {
        static if (extension is null)
            enum extensions = AliasSeq!(".dam", ".db");
        else
            enum extensions = AliasSeq!(extension);

        validateFileExtension!(null, extensions)(dbFile);
        validateFileExists(dbFile);

        foreach (hiddenDbFile; getHiddenDbFiles(dbFile))
        {
            validateFileExists!"cannot open hidden database file `%s`"(hiddenDbFile);
        }
    }

    void validateLasFile(in string lasFile, in string dbA=null, in string dbB=null)
    {
        auto cwd = getcwd.absolutePath;

        enforce!CLIException(
            lasFile.endsWith(".las"),
            format!"expected .las file, got `%s`"(lasFile),
        );
        validateFileExists(lasFile);
        enforce!CLIException(
            !lasEmpty(lasFile),
            format!"empty alignment file `%s`"(lasFile),
        );
    }

    void validateInputMasks(in string dbFile, in string[] maskDestinations)
    {
        foreach (maskDestination; maskDestinations)
            validateInputMask(dbFile, maskDestination);
    }

    void validateInputMask(in string dbFile, in string maskDestination)
    {
        foreach (maskFile; getMaskFiles(dbFile, maskDestination))
        {
            validateFileExists!"cannot open hidden mask file `%s`"(maskFile);
        }
    }

    void validateOutputMask(in string dbFile, in string maskDestination)
    {
        foreach (maskFile; getMaskFiles(dbFile, maskDestination))
        {
            validateFileWritable!"cannot write hidden mask file `%s`: %s"(maskFile);
        }
    }

    void validateFileWritable(string msg = "cannot open file `%s` for writing: %s")(string fileName)
    {
        auto deleteAfterwards = !fileName.exists;

        try
        {
            cast(void) File(fileName, "a");
        }
        catch (ErrnoException e)
        {
            throw new CLIException(format!msg(fileName, e.msg));
        }

        if (deleteAfterwards)
        {
            try
            {
                remove(fileName);
            }
            catch (FileException e)
            {
                logJsonWarn(
                    "info", "failed to delete file after testing",
                    "error", e.toString(),
                    "file", fileName,
                );
            }
        }
    }

    void validateCoordStrings(string[] coordStrings)
    {
        foreach (coordString; coordStrings)
            cast(void) parseCoordString(coordString);
    }

    OutputCoordinate parseCoordString(string coordString)
    {
        enum coordRegex = ctRegex!`^(scaffold/(?P<scaffoldId>\d+)/)?(contig/(?P<contigId>\d+)/)?(?P<coord>\d+)$`;
        OutputCoordinate coord;

        auto matches = coordString.matchFirst(coordRegex);

        enforce!CLIException(cast(bool) matches, "ill-formatted coord-string");

        coord.coord = matches["coord"].to!(typeof(coord.coord));
        enforce!CLIException(coord.coord > 0, "<coord> is 1-based");

        if (matches["contigId"] != "")
        {
            coord.contigId = matches["contigId"].to!(typeof(coord.contigId));
            enforce!CLIException(coord.contigId > 0, "<contig-id> is 1-based");
        }

        if (matches["scaffoldId"] != "")
        {
            coord.scaffoldId = matches["scaffoldId"].to!(typeof(coord.scaffoldId));
            enforce!CLIException(coord.scaffoldId > 0, "<scaffold-id> is 1-based");
        }

        version (unittest) { } else
            enforce!CLIException(
                coord.originType == OutputCoordinate.OriginType.scaffold,
                "not yet implemented; use format `scaffold/<uint:scaffold-id>/<uint:coord>`",
            );

        return coord;
    }

    unittest
    {
        import std.exception : assertThrown;

        {
            auto coord = parseCoordString(`scaffold/1125/42`);
            assert(coord.scaffoldId == 1125);
            assert(coord.contigId == 0);
            assert(coord.coord == 42);
            assert(coord.idx == 41);
            assert(coord.originType == OutputCoordinate.OriginType.scaffold);
        }
        {
            auto coord = parseCoordString(`scaffold/1125/contig/13/42`);
            assert(coord.scaffoldId == 1125);
            assert(coord.contigId == 13);
            assert(coord.coord == 42);
            assert(coord.idx == 41);
            assert(coord.originType == OutputCoordinate.OriginType.scaffoldContig);
        }
        {
            auto coord = parseCoordString(`contig/7/42`);
            assert(coord.scaffoldId == 0);
            assert(coord.contigId == 7);
            assert(coord.coord == 42);
            assert(coord.idx == 41);
            assert(coord.originType == OutputCoordinate.OriginType.contig);
        }
        {
            auto coord = parseCoordString(`42`);
            assert(coord.scaffoldId == 0);
            assert(coord.contigId == 0);
            assert(coord.coord == 42);
            assert(coord.idx == 41);
            assert(coord.originType == OutputCoordinate.OriginType.global);
        }
        {
            assertThrown!CLIException(parseCoordString(`scaffold/0/42`));
            assertThrown!CLIException(parseCoordString(`scaffold/1125/0`));
            assertThrown!CLIException(parseCoordString(`scaffold/0/contig/1/42`));
            assertThrown!CLIException(parseCoordString(`scaffold/1125/contig/0/42`));
            assertThrown!CLIException(parseCoordString(`scaffold/1125/contig/1/0`));
            assertThrown!CLIException(parseCoordString(`contig/0/42`));
            assertThrown!CLIException(parseCoordString(`contig/7/0`));
            assertThrown!CLIException(parseCoordString(`0`));
        }
    }

    void execIf(alias fun)(bool test)
    {
        if (test)
            fun();
    }

    void execUnless(alias fun)(bool test)
    {
        if (!test)
            fun();
    }
}

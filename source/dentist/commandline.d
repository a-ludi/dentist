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
import dentist.common.binio : PileUpDb;
import dentist.common.scaffold : JoinPolicy;
import dentist.dazzler :
    DaccordOptions,
    DalignerOptions,
    DamapperOptions,
    getHiddenDbFiles,
    getMaskFiles,
    getTracePointDistance,
    lasEmpty,
    LasFilterAlignmentsOptions,
    provideDamFileInWorkdir,
    provideLasFileInWorkdir,
    ProvideMethod,
    provideMethods;
import dentist.swinfo :
    copyright,
    executableName,
    description,
    license,
    version_;
import dentist.util.algorithm : staticPredSwitch;
import dentist.util.log;
import dentist.util.tempfile : mkdtemp;
import std.algorithm :
    among,
    each,
    endsWith,
    filter,
    find,
    map,
    startsWith;
import std.conv;
import std.exception : enforce, ErrnoException;
import std.file : exists, FileException, getcwd, isDir, tempDir, remove, rmdirRecurse;
import std.format : format, formattedRead;
import std.math : ceil, floor, log_e = log;
import std.meta : Alias, AliasSeq, staticMap, staticSort;
import std.parallelism : defaultPoolThreads, totalCPUs;
import std.path : absolutePath, buildPath;
import std.range : only, takeOne;
import std.regex : ctRegex, matchFirst;
import std.stdio : File, stderr;
import std.string : join, tr, wrap;
import std.traits :
    arity,
    EnumMembers,
    getSymbolsByUDA,
    getUDAs,
    isCallable,
    isStaticArray,
    Parameters,
    ReturnType;
import std.typecons : BitFlags, tuple;
import transforms : camelCase, snakeCaseCT;
import vibe.data.json : serializeToJsonString;


/// Possible returns codes of the command line execution.
enum ReturnCode
{
    ok,
    commandlineError,
    runtimeError,
}

/// Possible returns codes of the command line execution.
mixin("enum DentistCommand {" ~
    testingOnly!"translocateGaps," ~
    testingOnly!"findClosableGaps," ~
    "generateDazzlerOptions," ~
    "maskRepetitiveRegions," ~
    "showMask," ~
    "collectPileUps," ~
    "showPileUps," ~
    "processPileUps," ~
    "showInsertions," ~
    "mergeInsertions," ~
    "output," ~
    testingOnly!"translateCoords," ~
    testingOnly!"checkResults," ~
"}");

struct TestingCommand
{
    @disable this();

    static DentistCommand opDispatch(string command)() pure nothrow
    {
        static if (isTesting)
            return mixin("DentistCommand." ~ command);
        else
            return cast(DentistCommand) size_t.max;
    }
}

enum dashCase(string camelCase) = camelCase.snakeCaseCT.tr("_", "-");

private enum dentistCommands = staticMap!(
    dashCase,
    __traits(allMembers, DentistCommand),
);

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
    stderr.writeln(format!"%s %s"(executableName, version_));
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
struct OptionsFor(DentistCommand command)
{
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

    static if (command == DentistCommand.maskRepetitiveRegions)
    {
        @ArgumentsParser
        auto parseArguments(const(string)[] leftOver)
        {
            alias referenceSymbol = Alias!(__traits(getMember, this, "refFile"));
            enum referenceUDA = getUDAs!(referenceSymbol, Argument)[0];

            enforce!ArgParseError(leftOver.length >= numArguments.lowerBound, referenceUDA.multiplicityError(0));
            enforce!ArgParseError(leftOver.length <= numArguments.upperBound, "Missing positional arguments.");

            auto hasReadsFile = (leftOver.length == numArguments.upperBound);
            handleArg!"refFile"(this, leftOver[0]);
            leftOver = leftOver[1 .. $];

            if (hasReadsFile)
            {
                handleArg!"readsFile"(this, leftOver[0]);
                leftOver = leftOver[1 .. $];
            }

            handleArg!"dbAlignmentInputFile"(this, leftOver[0]);
            leftOver = leftOver[1 .. $];

            foreach (member; __traits(allMembers, typeof(this)))
            {
                alias symbol = Alias!(__traits(getMember, this, member));
                alias argUDAs = getUDAs!(symbol, Argument);

                static if (
                    argUDAs.length > 0 &&
                    !member.among("refFile", "readsFile", "dbAlignmentInputFile")
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
        TestingCommand.findClosableGaps,
        TestingCommand.checkResults,
    ))
    {
        @Argument("<in:true-assembly>")
        @Help("the 'true' assembly in .dam format")
        @Validate!(validateDB!".dam")
        string trueAssemblyFile;
        @Option()
        string trueAssemblyDb;

        @PostValidate()
        void hookProvideTrueAssemblyFileInWorkDir()
        {
            trueAssemblyDb = provideDamFileInWorkdir(trueAssemblyFile, provideMethod, workdir);
        }
    }

    static if (command.among(
        TestingCommand.translocateGaps,
    ))
    {
        @Argument("<in:short-read-assembly>")
        @Help("short-read assembly in .dam format")
        @Validate!(validateDB!".dam")
        string shortReadAssemblyFile;
        @Option()
        string shortReadAssemblyDb;

        @PostValidate()
        void hookProvideShortReadAssemblyFileInWorkDir()
        {
            shortReadAssemblyDb = provideDamFileInWorkdir(shortReadAssemblyFile, provideMethod, workdir);
        }
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
        DentistCommand.showMask,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
        DentistCommand.output,
    ))
    {
        @Argument("<in:reference>")
        @Help("reference assembly in .dam format")
        @Validate!(validateDB!".dam")
        string refFile;
        @Option()
        string refDb;

        @PostValidate()
        void hookProvideRefFileInWorkDir()
        {
            refDb = provideDamFileInWorkdir(refFile, provideMethod, workdir);
        }
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        static if (command == DentistCommand.maskRepetitiveRegions)
            enum argReadsMultiplicity = Multiplicity.optional;
        else
            enum argReadsMultiplicity = 1;

        @Argument("<in:reads>", argReadsMultiplicity)
        @Help("set of PacBio reads in .db/.dam format")
        @Validate!(validateReadsDb)
        string readsFile;
        @Option()
        string readsDb;

        static void validateReadsDb(string readsFile)
        {
            if (argReadsMultiplicity == 1 || readsFile !is null)
                validateDB(readsFile);
        }

        @property bool hasReadsDb() const pure nothrow
        {
            return readsFile !is null;
        }

        @PostValidate()
        void hookProvideReadsFileInWorkDir()
        {
            if (hasReadsDb)
                readsDb = provideDamFileInWorkdir(readsFile, provideMethod, workdir);
        }
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Argument("<in:result>")
        @Help("result assembly in .dam format")
        @Validate!(validateDB!".dam")
        string resultFile;
        @Option()
        string resultDb;

        @PostValidate()
        void hookProvideResultFileInWorkDir()
        {
            resultDb = provideDamFileInWorkdir(resultFile, provideMethod, workdir);
        }
    }

    static if (command.among(
        TestingCommand.translocateGaps,
    ))
    {
        @Argument("<in:short-vs-true-read-alignment>")
        @Help(q"{
            locals alignments of the short-read assembly against the 'true'
            assembly in form of a .las file as produced by `daligner`
        }")
        @Validate!((value, options) => validateLasFile(value, options.trueAssemblyFile, options.shortReadAssemblyFile))
        string shortReadAssemblyAlignmentInputFile;
        @Option()
        string shortReadAssemblyAlignmentFile;

        @PostValidate()
        void hookProvideShortReadAssemblyAlignmentInWorkDir()
        {
            shortReadAssemblyAlignmentFile = provideLasFileInWorkdir(
                shortReadAssemblyAlignmentInputFile,
                provideMethod,
                workdir,
            );
        }
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Argument("<in:alignment>")
        @Help("self-alignment of the reference assembly or reads vs. reference alignment")
        @Validate!((value, options) => options.readsFile is null
            ? validateLasFile(value, options.refFile)
            : validateLasFile(value, options.refFile, options.readsFile))
        string dbAlignmentInputFile;
        @Option()
        string dbAlignmentFile;

        @PostValidate()
        void hookProvideShortReadAssemblyAlignmentInWorkDir()
        {
            dbAlignmentFile = provideLasFileInWorkdir(
                dbAlignmentInputFile,
                provideMethod,
                workdir,
            );
        }
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
        @Validate!((value, options) => validateLasFile(value, options.refFile, options.readsFile))
        string readsAlignmentInputFile;
        @Option()
        string readsAlignmentFile;

        @PostValidate()
        void hookProvideReadsAlignmentInWorkDir()
        {
            readsAlignmentFile = provideLasFileInWorkdir(
                readsAlignmentInputFile,
                provideMethod,
                workdir,
            );
        }
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Argument("<in:result-vs-true-alignment>")
        @Help(q"{
            alignments chains of the result assembly against the 'true'
            assembly in form of a .las file as produced by `damapper`
        }")
        @Validate!((value, options) => validateLasFile(value, options.trueAssemblyFile, options.resultFile))
        string resultsAlignmentInputFile;
        @Option()
        string resultsAlignmentFile;

        @PostValidate()
        void hookProvideResultsAlignmentInWorkDir()
        {
            resultsAlignmentFile = provideLasFileInWorkdir(
                resultsAlignmentInputFile,
                provideMethod,
                workdir,
            );
        }
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
        DentistCommand.showMask,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @Argument("<in:repeat-mask>")
        @Help("read <repeat-mask> generated by the `mask-repetitive-regions` command")
        @Validate!((value, options) => validateInputMask(options.refFile, value))
        string repeatMask;
    }

    static if (command.among(
        TestingCommand.findClosableGaps,
        TestingCommand.checkResults,
    ))
    {
        @Argument("<in:mapped-regions-mask>")
        @Help(q"{
            read regions that were kept aka. output contigs from the Dazzler
            mask. Given a path-like string without extension: the `dirname`
            designates the directory to write the mask to. The mask comprises
            two hidden files `.[REFERENCE].[MASK].{anno,data}`.
        }")
        @Validate!((value, options) => validateInputMask(options.trueAssemblyFile, value))
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
        TestingCommand.translateCoords,
    ))
    {
        @Argument("<in:debug-graph>")
        @Help(q"{
            read the assembly graph from <debug-graph> generate
            (see `--debug-graph` of the `output` command)
        }")
        @Validate!validateFileExists
        string assemblyGraphFile;
    }

    static if (command.among(
        TestingCommand.translateCoords,
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
        @Help(q"{
            write regions that were kept aka. output contigs into a Dazzler
            mask. Given a path-like string without extension: the `dirname`
            designates the directory to write the mask to. The mask comprises
            two hidden files `.[REFERENCE].[MASK].{anno,data}`.
        }")
        @Validate!((value, options) => validateOutputMask(options.trueAssemblyFile, value))
        string mappedRegionsMask;
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
        @Help(q"{
            write inferred repeat mask into a Dazzler mask. Given a path-like
            string without extension: the `dirname` designates the directory to
            write the mask to. The mask comprises two hidden files
            `.[REFERENCE].[MASK].{anno,data}`.
        }")
        @Validate!((value, options) => validateOutputMask(options.refFile, value))
        string repeatMask;
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
        @Argument("<out:merged-insertions>")
        @Help("write merged insertion information to <merged-insertions>")
        @Validate!validateFileWritable
        string mergedInsertionsFile;

        @Argument("<in:insertions>", Multiplicity.oneOrMore)
        @Help("merge insertion information from <insertions>... generated by the `processPileUps` command")
        @Validate!validateFilesExist
        @Validate!"a.length >= 2"
        string[] insertionsFiles;
    }

    static if (command.among(
        TestingCommand.translocateGaps,
        DentistCommand.output,
    ))
    {
        @Argument("<out:assembly>", Multiplicity.optional)
        @Help("write output assembly to <assembly> (default: stdout)")
        @Validate!(value => (value is null).execUnless!(() => validateFileWritable(value)))
        string assemblyFile;
    }

    mixin HelpOption;

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("batch", "b")
        @MetaVar("<from>..<to>")
        @Help(q"{
            process only a subset of the pile ups in the given range (excluding <to>).
            <from> and <to> are zero-based indices into the pile up DB or <to>
            may be `$` to indicate the end of the pile up DB.
        }")
        void parsePileUpBatch(string batchString) pure
        {
            try
            {
                if (batchString.endsWith("$"))
                {
                    batchString.formattedRead!"%d..$"(pileUpBatch[0]);
                    pileUpBatch[1] = id_t.max;
                }
                else
                {
                    batchString.formattedRead!"%d..%d"(pileUpBatch[0], pileUpBatch[1]);
                }
            }
            catch (Exception e)
            {
                throw new CLIException("ill-formatted batch range");
            }
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
        @Validate!validateBatchRange
        id_t[2] pileUpBatch;

        static void validateBatchRange(id_t[2] pileUpBatch, OptionsFor!command options)
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

        @PostValidate()
        void hookEnsurePresenceOfBatchRange()
        {
            if (pileUpBatch == pileUpBatch.init || pileUpBatch[1] == id_t.max)
            {
                pileUpBatch[1] = pileUpLength;
            }
        }

        @property id_t pileUpBatchSize() const pure nothrow
        {
            return pileUpBatch[1] - pileUpBatch[0];
        }
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("bucket-size", "b")
        @Help(format!q"{
            bucket size of the gap length histogram; use 0 to disable (default: %d)
        }"(defaultValue!bucketSize))
        coord_t bucketSize = 500;
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("daccord-threads")
        @Help("use <uint> threads for `daccord` (defaults to floor(totalCpus / <threads>) )")
        uint numDaccordThreads;

        @PostValidate(Priority.low)
        void hookInitDaccordThreads()
        {
            if (numDaccordThreads == 0)
            {
                numDaccordThreads = totalCPUs / numThreads;
            }
        }
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        @Option("debug-scaffold")
        @MetaVar("<file>")
        @Help("write the assembly scaffold to <file>; use `show-insertions` to inspect the result")
        @Validate!(value => (value is null).execUnless!(() => validateFileWritable(value)))
        string assemblyGraphFile;
    }


    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("only-contig")
        @Help("restrict analysis to contig <uint> (experimental)")
        id_t onlyContigId;
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("debug-alignment")
        @MetaVar("<file>")
        @Help("write the result alignment to a tabular file <file>")
        @Validate!(value => (value is null).execUnless!(() => validateFileWritable(value)))
        string alignmentTabular;
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("debug-gap-details")
        @MetaVar("<file>")
        @Help("write the statistics for every single gap to a tabular file <file>")
        @Validate!(value => (value is null).execUnless!(() => validateFileWritable(value)))
        string gapDetailsTabular;
    }

    static if (command.among(
        DentistCommand.collectPileUps,
    ))
    {
        @Option("best-pile-up-margin")
        @Help(q"{
            given a set of possibly of conflicting pile ups, if the largest
            has <double> times more reads than the second largest it is
            considered unique
        }")
        @Validate!(value => enforce!CLIException(value > 1.0, "--best-pile-up-margin must be greater than 1.0"))
        double bestPileUpMargin = 3.0;
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        @Option("extend-contigs")
        @Help("if given extend contigs even if no spanning reads can be found")
        OptionFlag shouldExtendContigs;
    }

    static if (command.among(
        TestingCommand.translocateGaps,
        DentistCommand.output,
    ))
    {
        @Option("fasta-line-width", "w")
        @Help(format!"line width for ouput FASTA (default: %d)"(defaultValue!fastaLineWidth))
        @Validate!(value => enforce!CLIException(value > 0, "fasta line width must be greater than zero"))
        size_t fastaLineWidth = 50;
    }

    static if (needWorkdir)
    {
        @Option("input-provide-method", "p")
        @MetaVar(format!"{%-(%s,%)}"([provideMethods]))
        @Help(format!q"{
            use the given method to provide the input files in the working
            directory (default: `%s`)
        }"(defaultValue!provideMethod))
        ProvideMethod provideMethod = ProvideMethod.symlink;
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
        TestingCommand.translateCoords,
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
        DentistCommand.showMask,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @Option("mask", "m")
        @MetaVar("<string>...")
        @Help("additional masks (see <in:repeat-mask>)")
        void addMask(string mask) pure
        {
            additionalMasks ~= mask;
        }

        @Option()
        @Validate!((values, options) => validateInputMasks(options.refFile, values))
        string[] additionalMasks;
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
                "must provide either --read-coverage or --acceptable-coverage-reads",
            );
            enforce!CLIException(
                (maxCoverageReads != maxCoverageReads.init) ^
                (readCoverage != readCoverage.init),
                "must not provide both --read-coverage and --acceptable-coverage-reads",
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
                "minimum anchor length should be greater than --trace-point-spacing"
            ),
            is(typeof(OptionsFor!command().tracePointDistance)),
        )
        size_t minAnchorLength = 500;
    }

    static if (command.among(
        TestingCommand.checkResults,
    ))
    {
        @Option("min-insertion-length")
        @Help(format!q"{
            an insertion must have at least this num ber base pairs to be
            considered as (partial) insertion (default: %d)
        }"(defaultValue!minInsertionLength))
        @Validate!(value => enforce!CLIException(value > 0, "minimum insertion length must be greater than zero"))
        size_t minInsertionLength = 50;
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("min-extension-length")
        @Help(format!q"{
            extensions must have at least <ulong> bps of consensus to be inserted (default: %d)
        }"(defaultValue!minExtensionLength))
        @Validate!(value => enforce!CLIException(value > 0, "minimum extension length must be greater than zero"))
        size_t minExtensionLength = 100;
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
            require at least <uint> spanning reads to close a gap (default: %d)
        }"(defaultValue!minSpanningReads))
        size_t minSpanningReads = defaultMinSpanningReads;
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Option("read-coverage", "C")
        @Help(q"{
            this is used to provide good default values for --acceptable-coverage-reads;
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
        @Help("estimated error rate in reads")
        @Validate!(value => enforce!CLIException(
            0.0 < value && value <= 0.3,
            "reads error rate must be in (0, 0.3]"
        ))
        double readsErrorRate = .15;
    }

    static if (command.among(
        DentistCommand.generateDazzlerOptions,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
    ))
    {
        @Option("reference-error")
        @Help("estimated error rate in reference")
        @Validate!(value => enforce!CLIException(
            0.0 < value && value <= 0.3,
            "reference error rate must be in (0, 0.3]"
        ))
        double referenceErrorRate = .01;
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
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
        DentistCommand.output,
        TestingCommand.checkResults,
    ))
    {
        @Option("trace-point-spacing", "s")
        @Help("trace point spacing used for the ref vs. reads alignment")
        trace_point_t tracePointDistance;

        @PostValidate()
        void hookEnsurePresenceOfTracePointDistance()
        {
            if (tracePointDistance > 0)
                return;

            tracePointDistance = getTracePointDistance();
        }
    }

    @Option("verbose", "v")
    @Help("increase output to help identify problems; use up to three times")
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

    static if (needWorkdir)
    {
        /**
            Last part of the working directory name. A directory in the temp
            directory as returned by `std.file.tmpDir` with the naming scheme will
            be created to hold all data for the computation.
        */
        enum workdirTemplate = format!"dentist-%s-XXXXXX"(command);

        /// This is a temporary directory to store all working data.
        @Option("workdir", "w")
        @Help("use <string> as a working directory")
        @Validate!(value => enforce!CLIException(
            value is null || value.isDir,
            format!"workdir is not a directory: %s"(value),
        ))
        string workdir;

        @PostValidate(Priority.high)
        void hookCreateWorkDir()
        {
            if (workdir !is null)
                return;

            auto workdirTemplate = buildPath(tempDir(), workdirTemplate);

            workdir = mkdtemp(workdirTemplate);
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

    @PostValidate()
    void hookInitLogLevel()
    {
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

    static if (
        is(typeof(OptionsFor!command().minAnchorLength)) &&
        is(typeof(OptionsFor!command().referenceErrorRate))
    ) {
        @Validate!validateAverageCorrelationRate
        @property auto selfAlignmentOptionsAverageCorrelationRate() const
        {
            return (1 - referenceErrorRate)^^2;
        }

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
        is(typeof(OptionsFor!command().minAnchorLength)) &&
        is(typeof(OptionsFor!command().readsErrorRate)) &&
        is(typeof(OptionsFor!command().numDaccordThreads))
    ) {
        @Validate!validateAverageCorrelationRate
        @property auto pileUpAlignmentOptionsAverageCorrelationRate() const
        {
            return (1 - readsErrorRate)^^2;
        }

        @property string[] pileUpAlignmentOptions() const
        {
            return [
                DalignerOptions.identity,
                DalignerOptions.numThreads ~ numDaccordThreads.to!string,
                format!(DalignerOptions.minAlignmentLength ~ "%d")(minAnchorLength),
                format!(DalignerOptions.averageCorrelationRate ~ "%f")(
                    pileUpAlignmentOptionsAverageCorrelationRate,
                ),
            ];
        }
    }

    static if (
        is(typeof(OptionsFor!command().minAnchorLength)) &&
        is(typeof(OptionsFor!command().readsErrorRate)) &&
        is(typeof(OptionsFor!command().referenceErrorRate)) &&
        is(typeof(OptionsFor!command().numDaccordThreads))
    ) {
        @property string[] postConsensusAlignmentOptions() const
        {
            return [
                DalignerOptions.asymmetric,
                DalignerOptions.numThreads ~ numDaccordThreads.to!string,
                format!(DalignerOptions.minAlignmentLength ~ "%d")(tracePointDistance),
                format!(DalignerOptions.averageCorrelationRate ~ "%f")((1 - referenceErrorRate^^minReadsPerPileUp) * (1 - readsErrorRate)),
            ];
        }
    }

    static if (isTesting)
    {
        @property string[] trueAssemblyVsResultAlignmentOptions() const
        {
            return [
                DamapperOptions.symmetric,
                DamapperOptions.oneDirection,
                DamapperOptions.averageCorrelationRate ~ ".7",
            ];
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
            string workdir;
        }

        @property auto consensusOptions() const
        {
            return const(ConsensusOptions)(
                // daccordOptions
                [
                    DaccordOptions.produceFullSequences,
                    DaccordOptions.numberOfThreads ~ numDaccordThreads.to!string,
                ],
                // dalignerOptions
                pileUpAlignmentOptions,
                // dbsplitOptions
                [],
                // lasFilterAlignmentsOptions
                [
                    LasFilterAlignmentsOptions.errorThresold ~ (2.0 * readsErrorRate).to!string,
                ],
                // workdir
                workdir,
            );
        }
    }

    static auto defaultValue(alias property)() pure nothrow
    {
        OptionsFor!command defaultOptions;

        return __traits(getMember, defaultOptions, property.stringof);
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
    static if (command == TestingCommand.translocateGaps)
        enum commandSummary = q"{
            Translocate gaps from first assembly to second assembly.
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
    else static if (command == TestingCommand.translateCoords)
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
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

private
{
    ReturnCode runCommand(DentistCommand command)(in string[] args)
    {
        alias Options = OptionsFor!command;
        enum commandName = command.to!string.snakeCaseCT.tr("_", "-");
        enum usage = usageString!Options(executableName ~ " " ~ commandName);

        Options options;

        try
        {
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
            stderr.writeln("Error: " ~ (shouldLog(LogLevel.diagnostic)
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
        static foreach (alias symbol; getSymbolsByUDA!(Options, Validate))
        {{
            alias validate = getUDAs!(symbol, Validate)[0].validate;
            auto value = __traits(getMember, options, __traits(identifier, symbol));
            alias Value = typeof(value);
            alias Validator = typeof(validate);


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

    import std.meta : allSatisfy, staticMap;
    import std.traits : isSomeString;

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

    void validateLasFile(in string lasFile, in string dbA, in string dbB=null)
    {
        auto cwd = getcwd.absolutePath;

        enforce!CLIException(
            lasFile.endsWith(".las"),
            format!"expected .las file, got `%s`"(lasFile),
        );
        validateFileExists(lasFile);
        enforce!CLIException(
            !lasEmpty(lasFile, dbA, dbB, cwd),
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

/**
    Defines the behavior of the `dentist` command line client.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commandline;

import darg :
    ArgParseHelp,
    Argument,
    Help,
    helpString,
    MetaVar,
    Multiplicity,
    Option,
    OptionFlag,
    parseArgs,
    usageString;
import dentist.common.alignments :
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
import dentist.util.log;
import dentist.util.tempfile : mkdtemp;
import std.algorithm : among, each, endsWith, filter, find, map, startsWith;
import std.conv;
import std.exception : enforce, ErrnoException;
import std.file : exists, FileException, getcwd, isDir, tempDir, remove, rmdirRecurse;
import std.format : format, formattedRead;
import std.meta : AliasSeq, staticMap, staticSort;
import std.parallelism : defaultPoolThreads, totalCPUs;
import std.path : absolutePath, buildPath;
import std.range : only, takeOne;
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
import std.typecons : BitFlags;
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
enum DentistCommand
{
    generateDazzlerOptions,
    maskRepetitiveRegions,
    showMask,
    collectPileUps,
    showPileUps,
    processPileUps,
    showInsertions,
    mergeInsertions,
    output,
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
        stderr.writeln("Error: " ~ e.msg);
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
        stderr.writeln(e.to!string);

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
        DentistCommand.maskRepetitiveRegions,
        DentistCommand.showMask,
        DentistCommand.collectPileUps,
        DentistCommand.processPileUps,
        DentistCommand.output,
    );

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
        @Argument("<in:reads>")
        @Help("set of PacBio reads in .dam format")
        @Validate!validateDB
        string readsFile;
        @Option()
        string readsDb;

        @PostValidate()
        void hookProvideReadsFileInWorkDir()
        {
            readsDb = provideDamFileInWorkdir(readsFile, provideMethod, workdir);
        }
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Argument("<in:self-alignment>")
        @Help(q"{
            local alignments of the reference against itself in form of a .las
            file as produced by `daligner`
        }")
        @Validate!((value, options) => validateLasFile(value, options.refFile))
        string selfAlignmentInputFile;
        @Option()
        string selfAlignmentFile;

        @PostValidate()
        void hookProvideSelfAlignmentInWorkDir()
        {
            selfAlignmentFile = provideLasFileInWorkdir(
                selfAlignmentInputFile,
                provideMethod,
                workdir,
            );
        }
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
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
        @Help("read <repeat-mask> generated by the `maskRepetitiveRegions` command")
        @Validate!((value, options) => validateInputMask(options.refFile, value))
        string repeatMask;
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
        DentistCommand.output,
    ))
    {
        @Argument("<out:assembly>", Multiplicity.optional)
        @Help("write output assembly to <assembly> (default: stdout)")
        @Validate!validateFileWritable
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
            process only a subset of the pile ups in the given range (excluding <to>);
            <from> and <to> are zero-based indices into the pile up DB
        }")
        void parsePileUpBatch(string batchString) pure
        {
            try
            {
                batchString.formattedRead!"%d..%d"(pileUpBatch[0], pileUpBatch[1]);
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
                (0 <= from && from < to && to <= options.pileUpLength),
                format!"invalid batch range; check that 0 <= <from> < <to> <= %d"(
                        options.pileUpLength)
            );
        }

        @PostValidate()
        void hookEnsurePresenceOfBatchRange()
        {
            if (pileUpBatch == pileUpBatch.init)
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
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Option("coverage-reads", "C")
        @MetaVar("<from>..<to>")
        @Help(q"{
            this is used to derive a repeat mask from the ref vs. reads alignment;
            if the alignment coverage is out of the interval [<from>, <to>]
            it will be considered repetitive
        }")
        void parseCoverageBoundsReads(string rangeString) pure
        {
            parseRange!coverageBoundsReads(rangeString);
        }

        @Option()
        @Validate!(validateCoverageBounds!"reads")
        id_t[2] coverageBoundsReads;
    }

    static if (command.among(
        DentistCommand.maskRepetitiveRegions,
    ))
    {
        @Option("coverage-self", "c")
        @MetaVar("<from>..<to>")
        @Help(q"{
            this is used to derive a repeat mask from the self alignment;
            if the alignment coverage is out of the interval [<from>, <to>]
            it will be considered repetitive
        }")
        void parseCoverageBoundsSelf(string rangeString) pure
        {
            parseRange!coverageBoundsSelf(rangeString);
        }

        @Option()
        @Validate!(validateCoverageBounds!"self")
        id_t[2] coverageBoundsSelf;
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
        @Option("extend-contigs")
        @Help("if given extend contigs even if no spanning reads can be found")
        OptionFlag shouldExtendContigs;
    }

    static if (command.among(
        DentistCommand.output,
    ))
    {
        @Option("fasta-line-width", "w")
        @Help(format!"line width for ouput FASTA (default: %d)"(defaultValue!fastaLineWidth))
        @Validate!(value => enforce!CLIException(value > 0, "fasta line width must be greater than zero"))
        size_t fastaLineWidth = 50;
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("good-anchor-length")
        @Help(format!q"{
            alignment anchors with at least this length will get no penalty (default: %d)
        }"(defaultValue!goodAnchorLength))
        @Validate!(value => enforce!CLIException(value > 0, "good anchor length must be greater than zero"))
        size_t goodAnchorLength = 1000;
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
        size_t minAnchorLength = 200;
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("min-extension-length")
        @Help(format!q"{
            extensions must have at least <ulong> bps of processPileUps to be inserted (default: %d)
        }"(defaultValue!minExtensionLength))
        @Validate!(value => enforce!CLIException(value > 0, "good extension length must be greater than zero"))
        size_t minExtensionLength = 100;
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("min-reads-per-pile-up")
        @Help(format!q"{
            alignment anchors with at least this length will get no penalty (default: %d)
        }"(defaultValue!minReadsPerPileUp))
        @Validate!(value => enforce!CLIException(value > 0, "min reads per pile up must be greater than zero"))
        size_t minReadsPerPileUp = 5;
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
            0.0 < value && value < 1.0,
            "reads error rate must be in (0, 1)"
        ))
        double readsErrorRate = .15;
    }

    static if (command.among(
        DentistCommand.generateDazzlerOptions,
        DentistCommand.collectPileUps,
    ))
    {
        @Option("reference-error")
        @Help("estimated error rate in reference")
        @Validate!(value => enforce!CLIException(
            0.0 < value && value < 1.0,
            "reference error rate must be in (0, 1)"
        ))
        double referenceErrorRate = .01;
    }

    static if (command.among(
        DentistCommand.processPileUps,
    ))
    {
        @Option("threads", "T")
        @Help("use <uint> threads (defaults to the number of cores)")
        uint numThreads;

        @PostValidate()
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
        @property string[] selfAlignmentOptions() const
        {
            return [
                DalignerOptions.identity,
                format!(DalignerOptions.minAlignmentLength ~ "%d")(minAnchorLength),
                format!(DalignerOptions.averageCorrelationRate ~ "%f")((1 - referenceErrorRate)^^2),
            ];
        }
    }

    static if (
        is(typeof(OptionsFor!command().referenceErrorRate)) &&
        is(typeof(OptionsFor!command().readsErrorRate))
    ) {
        @property string[] refVsReadsAlignmentOptions() const
        {
            return [
                DamapperOptions.symmetric,
                DamapperOptions.bestMatches ~ ".7",
                format!(DamapperOptions.averageCorrelationRate ~ "%f")((1 - referenceErrorRate) * (1 - readsErrorRate)),
            ];
        }
    }

    static if (
        is(typeof(OptionsFor!command().minAnchorLength)) &&
        is(typeof(OptionsFor!command().readsErrorRate))
    ) {
        @property string[] pileUpAlignmentOptions() const
        {
            return [
                DalignerOptions.identity,
                format!(DalignerOptions.minAlignmentLength ~ "%d")(minAnchorLength),
                format!(DalignerOptions.averageCorrelationRate ~ "%f")((1 - readsErrorRate)^^2),
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
            string workdir;
        }

        @property auto consensusOptions() const
        {
            return const(ConsensusOptions)(
                [
                    DaccordOptions.produceFullSequences,
                    DaccordOptions.numberOfThreads ~ numDaccordThreads.to!string,
                ],
                pileUpAlignmentOptions,
                [],
                workdir,
            );
        }
    }

    static auto defaultValue(alias property)() pure nothrow
    {
        OptionsFor!command defaultOptions;

        return __traits(getMember, defaultOptions, property.stringof);
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
    static if (command == DentistCommand.generateDazzlerOptions)
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
            stderr.writeln("Error: " ~ e.msg);
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
            stderr.writeln("Error: " ~ e.to!string);

            return ReturnCode.runtimeError;
        }
    }

    enum getUDA(alias symbol, T) = getUDAs!(symbol, T)[0];

    struct Validate(alias _validate) {
        alias validate = _validate;
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
            auto value = __traits(getMember, options, symbol.stringof);
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

    void validateCoverageBounds(string what)(id_t[2] coverageBounds)
    {
        auto from = coverageBounds[0];
        auto to = coverageBounds[1];

        enforce!CLIException(
            0 <= from && from < to,
            coverageBounds == coverageBounds.init
                ? "missing --coverage-" ~ what
                : "invalid coverage bounds (" ~ what ~ "); check that 0 <= <from> < <to>"
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

    void validateDB(string extension = null)(in string dbFile)
        if (extension is null || extension.among(".dam", ".db"))
    {
        static if (extension is null)
            enum extensions = AliasSeq!(".dam", ".db");
        else
            enum extensions = AliasSeq!(extension);

        enforce!CLIException(
            dbFile.endsWith(extensions),
            format!("expected " ~ [extensions].join(" or ") ~ " file, got `%s`")(dbFile),
        );
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
}

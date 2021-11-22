/**
    Central logging facility for DENTIST.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.log;

import std.array;
import std.datetime;
import std.format;
import std.process;
import std.range;
import std.stdio;
import core.thread;


private
{
    __gshared LogLevel minLevel = LogLevel.info;
}


/// Sets the minimum log level to be printed.
void setLogLevel(LogLevel level) nothrow
{
    synchronized
        minLevel = level;
}


/// Get the minimum log level to be printed. Use `shouldLog` for conditionals.
LogLevel getLogLevel()
{
    return minLevel;
}


/// Check whether message of `logLevel` should be logged.
bool shouldLog(LogLevel logLevel)
{
    return logLevel >= minLevel;
}


/**
    Logs a message in compressed single-line JSON format.

    Produces a JSON object with the key-value pairs given as `args` and
    default fields `"thread"`, `"timestamp"` and `"logLevel"`.


    Example:
    ---
    logJsonInfo(
        "action", "findTheTruth",
        "answer", 42,
        "elapsedSecs", 1337,
    );

    // --> (real output is compressed in a single line)
    // {
    //     "action": "findTheTruth",
    //     "answer": 42,
    //     "elapsedSecs": 1337
    //     "thread": 123467890,
    //     "timestamp", 123467890
    //     "logLevel", "info"
    // }
    ---

    Params:
        args = pairs of `name` (`string`) and `value`
*/
void logJsonDebug(T...)(lazy T args) nothrow
{
    logJson(LogLevel.debug_, args);
}
/// ditto
void logJsonDiagnostic(T...)(lazy T args) nothrow
{
    logJson(LogLevel.diagnostic, args);
}
/// ditto
void logJsonInfo(T...)(lazy T args) nothrow
{
    logJson(LogLevel.info, args);
}
/// ditto
void logJsonWarn(T...)(lazy T args) nothrow
{
    logJson(LogLevel.warn, args);
}
/// ditto
void logJsonError(T...)(lazy T args) nothrow
{
    logJson(LogLevel.error, args);
}

/// ditto
void logJson(T...)(LogLevel level, lazy T args) nothrow
{
    import dentist.util.range : Chunks, chunks;
    import std.conv : to;
    import std.datetime.systime : Clock;
    import std.traits : isSomeString;
    import vibe.data.json : Json;

    enum threadKey = "thread";
    enum timestampKey = "timestamp";
    enum logLevelKey = "logLevel";

    if (level < minLevel)
        return;

    try
    {
        Json json = Json.emptyObject;

        static foreach (KeyValuePair; Chunks!(2, T))
        {
            static assert(isSomeString!(KeyValuePair.chunks[0]), "missing name");
        }

        try
        {
            json[logLevelKey] = level.to!string;
        }
        catch(Exception e)
        {
            json[logLevelKey] = null;
        }

        json[timestampKey] = Clock.currStdTime;
        json[threadKey] = thisThreadID;

        foreach (keyValuePair; args.chunks!2)
        {
            json[keyValuePair[0]] = keyValuePair[1];
        }

        return log(level, json.to!string);
    }
    catch (Exception e)
    {
        // this is bad but what can we do..
        debug assert(false, e.msg);
    }
}

///
unittest
{
    import std.stdio : File, stderr;
    import vibe.data.json : Json, parseJsonString;

    auto origStderr = stderr;
    stderr = File.tmpfile();
    scope (exit)
    {
        stderr.close();
        stderr = origStderr;
    }

    logJsonError(
        "error", "mysterious observation",
        "secret", 42,
    );

    stderr.rewind();
    auto observed = parseJsonString(stderr.readln);

    assert(observed["thread"].type == Json.Type.int_);
    assert(observed["timestamp"].type == Json.Type.int_);
    assert(observed["error"] == "mysterious observation");
    assert(observed["secret"] == 42);
}

/**
    Logs a message.

    Params:
        level = The log level for the logged message
        fmt = See $(LINK http://dlang.org/phobos/std_format.html#format-string)
*/
void logDebug(T...)(string fmt, lazy T args) nothrow
{
    log(LogLevel.debug_, fmt, args);
}
/// ditto
void logDiagnostic(T...)(string fmt, lazy T args) nothrow
{
    log(LogLevel.diagnostic, fmt, args);
}
/// ditto
void logInfo(T...)(string fmt, lazy T args) nothrow
{
    log(LogLevel.info, fmt, args);
}
/// ditto
void logWarn(T...)(string fmt, lazy T args) nothrow
{
    log(LogLevel.warn, fmt, args);
}
/// ditto
void logError(T...)(string fmt, lazy T args) nothrow
{
    log(LogLevel.error, fmt, args);
}

/// ditto
void log(T...)(LogLevel level, string fmt, lazy T args) nothrow
{
    if (level < minLevel)
        return;

    try
    {
        auto txt = appender!string();
        txt.reserve(256);
        static if (args.length > 0)
        {
            formattedWrite(txt, fmt, args);
        }
        else
        {
            txt ~= fmt;
        }

        if (level >= minLevel)
        {
            File output = stderr;

            synchronized if (output.isOpen)
            {
                output.writeln(txt.data);
                output.flush();
            }
        }
    }
    catch (Exception e)
    {
        // this is bad but what can we do..
        debug assert(false, e.msg);
    }
}

/// Specifies the log level for a particular log message.
enum LogLevel
{
    debug_,
    diagnostic,
    info,
    warn,
    error,
    fatal,
    none
}


/// Do not use directly. Use `mixin(traceExecution)` instead.
struct ExecutionTracer(LogLevel logLevel = LogLevel.diagnostic)
{
    import std.datetime.stopwatch : StopWatch;
    import std.typecons : Yes;

    string functionName;
    StopWatch timer;

    this(int dummy, string fnName = __FUNCTION__)
    {
        this.functionName = fnName;

        logJson(
            logLevel,
            `state`, `enter`,
            `function`, this.functionName,
        );

        this.timer = StopWatch(Yes.autoStart);
    }

    ~this()
    {
        timer.stop();

        logJson(
            logLevel,
            `state`, `exit`,
            `function`, functionName,
            `timeElapsed`, timer.peek().total!`hnsecs`,
        );
    }
}


/// Print JSON log entries upon entering and leaving the function reporting
/// the execution time.
///
/// Example:
/// ---
/// void foo()
/// {
///     mixin(traceExecution);
///
///     logJsonInfo("info", "working on foo()")
/// }
///
/// void main()
/// {
///     foo();
///     // --> (real output is compressed in a single line)
///     // {
///     //     "function": "foo",
///     //     "state": "enter",
///     //     "thread": 123467890,
///     //     "timestamp", 123467890
///     //     "logLevel", "info"
///     // }
///     // {
///     //     "info": "working on foo()",
///     //     "thread": 123467890,
///     //     "timestamp", 123467890
///     //     "logLevel", "info"
///     // }
///     // {
///     //     "function": "foo",
///     //     "state": "exit",
///     //     "timeElapsed": 1234567890,
///     //     "thread": 123467890,
///     //     "timestamp", 123467890
///     //     "logLevel", "info"
///     // }
/// }
/// ---
string traceExecution(LogLevel logLevel = LogLevel.diagnostic)()
{
    import std.conv : to;
    import std.string : replace;
    import std.traits : moduleName;

    return q"{
        static import $thisModule;

        scope __executionTracer = $thisModule.ExecutionTracer!($logLevel)(0);
    }"
        .replace("$thisModule", moduleName!LogLevel)
        .replace("$logLevel", "LogLevel." ~ logLevel.to!string);
}

unittest
{
    import std.regex : ctRegex, matchFirst;
    import std.stdio : File, stderr;
    import vibe.data.json : Json, parseJsonString;

    auto origStderr = stderr;
    stderr = File.tmpfile();
    scope (exit)
    {
        stderr.close();
        stderr = origStderr;
    }

    void doSomething()
    {
        mixin(traceExecution!(LogLevel.error));

        import core.thread : Thread;
        import core.time : dur;

        Thread.getThis.sleep(dur!"hnsecs"(50));
    }

    doSomething();
    stderr.rewind();

    enum functionFQN = ctRegex!`dentist\.util\.log\.__unittest_L[0-9]+_C[0-9]+\.doSomething`;
    auto observed1 = parseJsonString(stderr.readln);
    auto observed2 = parseJsonString(stderr.readln);

    assert(observed1["thread"].type == Json.Type.int_);
    assert(observed1["timestamp"].type == Json.Type.int_);
    assert(observed1["state"] == "enter");
    assert(matchFirst(observed1["function"].to!string, functionFQN));

    assert(observed2["thread"].type == Json.Type.int_);
    assert(observed2["timestamp"].type == Json.Type.int_);
    assert(observed2["state"] == "exit");
    assert(matchFirst(observed2["function"].to!string, functionFQN));
}


/// Tracks progress and outputs information regularly.
///
/// `ProgressMeter.Format` can be used to choose a format suitable for
/// terminals (`ProgressMeter.Format.human`) or for log files
/// (`ProgressMeter.Format.json`).
struct ProgressMeter
{
    import std.datetime.stopwatch : StopWatch;
    import std.algorithm : max;
    import std.stdio :
        File,
        stderr;
    import std.typecons :
        Flag,
        No,
        Tuple,
        Yes;

    private alias UnitSpec = Tuple!(size_t, "multiplier", char, "name");


    /// Available display units for the progress meter.
    enum Unit : UnitSpec
    {
        auto_ = UnitSpec(0, '\0'),
        one = UnitSpec(1, ' '),
        kilo = UnitSpec(10^^3, 'k'),
        mega = UnitSpec(10^^6, 'M'),
        giga = UnitSpec(10^^9, 'G'),
        peta = UnitSpec(10^^12, 'P'),
        min = one,
        max = peta,
    }


    /// Display format of the progress meter.
    enum Format : ubyte
    {
        /// Displays a single line that is updated regularly. This is suitable
        /// for terminal output.
        human,

        /// Produces a series of compressed, single-line JSON object
        /// describing the progress. This is suitable for output to a regular
        /// file.
        json,
    }


    /// Wait at least the amount of milliseconds before updating the status.
    size_t printEveryMsecs = 500;

    /// Display ticks in this unit.
    Unit unit;

    /// Use `precision` digits after the decimal point.
    size_t precision = 3;

    /// Specifies the 100% mark if given. No percentage is displayed if this
    /// is zero.
    size_t totalTicks;

    /// Number of ticks until now.
    size_t numTicks;

    /// Suppress output but still track the progress.
    Flag!"silent" silent;

    /// Choose display format.
    Format format;

    private File _output;
    private bool hasOutput;
    private StopWatch timer;
    private StopWatch lastPrint;


    /// Set the output file to write status updates to. Default it to use
    /// `std.stdio.stderr`.
    @property void output(File output)
    {
        hasOutput = true;
        _output = output;
    }


    /// Get a reference to the output file.
    @property auto ref File output()
    {
        if (!hasOutput)
            output = stderr;

        return _output;
    }


    /// Start the timer.
    ///
    /// This implicitly resets the timer and tick count.
    void start()
    {
        numTicks = 0;
        if (!silent)
        {
            lastPrint.reset();
            lastPrint.start();
            printProgressLine(LineLocation.first);
        }
        timer.reset();
        timer.start();
    }


    /// Add a single tick.
    void tick()
    {
        ++numTicks;

        if (!silent && lastPrint.peek.total!"msecs" > printEveryMsecs)
            printProgressLine(LineLocation.middle);
    }


    /// Stop the timer and print a last status update.
    void stop()
    {
        timer.stop();

        if (!silent)
            printProgressLine(LineLocation.last);
    }


    /// Check if `timeUnit` is allowed by `std.datetime.stopwatch.StopWatch.peek.total`.
    static enum isValidTimeUnit(string timeUnit) = is(typeof(timer.peek.total!timeUnit));
    static assert(isValidTimeUnit!"msecs");


    /// Get the number of elapsed `timeUnit`s.
    @property auto elapsed(string timeUnit)() const nothrow @safe if (isValidTimeUnit!timeUnit)
    {
        return timer.peek.total!timeUnit;
    }


    /// Get the average throughput in ticks per `timeUnit`.
    ///
    /// Bugs: this will cause an arithmetic error if no time has elapsed.
    @property auto ticksPer(string timeUnit)() const nothrow @safe if (isValidTimeUnit!timeUnit)
    {
        return cast(double) numTicks / elapsed!timeUnit;
    }


    /// True if the estimated time of arrival (ETA) can be calculated.
    @property auto hasETA() const nothrow @safe
    {
        return totalTicks > 0 && numTicks > 0;
    }

    alias hasEstimatedTimeOfArrival = hasETA;


    /// Calculate the estimated time of arrival (ETA).
    ///
    /// This simply assumes that the current average throughput
    /// (`ticksPer!timeUnit`) will not change.
    @property auto eta(string timeUnit)() const nothrow @safe if (isValidTimeUnit!timeUnit)
    {
        return (totalTicks - numTicks)/ticksPer!timeUnit;
    }

    /// ditto
    alias estimatedTimeOfArrival = eta;


    /// Select the smallest unit such that the number of decimal digits
    /// is up to three.
    static Unit selectUnitFor(size_t number) pure nothrow @safe
    {
        import std.traits : EnumMembers;

        foreach (unit; EnumMembers!Unit)
            if (unit.multiplier > 0 && number / unit.multiplier < 1000)
                return unit;
        return Unit.max;
    }

    ///
    unittest
    {
        assert(selectUnitFor(1) == Unit.one);
        assert(selectUnitFor(25 * 1024) == Unit.kilo);
        assert(selectUnitFor(13 * 1024*1024) == Unit.mega);
    }


private:


    enum LineLocation : ubyte
    {
        first,
        middle,
        last,
    }


    void printProgressLine(LineLocation lineLocation)
    {
        final switch (format)
        {
            case Format.human:
                printHumanProgressLine(lineLocation);
                break;
            case Format.json:
                printJsonProgressLine(lineLocation);
                break;
        }

        lastPrint.reset();
    }


    void printHumanProgressLine(LineLocation lineLocation)
    {
        enum progressFormat = "\rrecords: %04.*f%c  elapsed: %04.*f sec  rate: %04.*f records/sec";
        enum progressFormatWithTotal = "\rrecords: %04.*f/%04.*f%c (%04.2f%%) eta: %04.*f sec  elapsed: %04.*f sec  rate: %04.*f records/sec";
        auto elapsedSecs = timer.peek.total!"msecs" * 1e-3;

        auto unit = this.unit == Unit.auto_
            ? selectUnitFor(max(numTicks, totalTicks))
            : this.unit;

        if (totalTicks == 0)
            output.writef!progressFormat(
                precision,
                cast(double) numTicks / unit.multiplier,
                unit.name,
                precision,
                elapsedSecs,
                precision,
                cast(double) numTicks / elapsedSecs,
            );
        else
            output.writef!progressFormatWithTotal(
                precision,
                cast(double) numTicks / unit.multiplier,
                precision,
                cast(double) totalTicks / unit.multiplier,
                unit.name,
                (100.0 * numTicks / totalTicks),
                precision,
                eta!"seconds",
                precision,
                elapsedSecs,
                precision,
                cast(double) numTicks / elapsedSecs,
            );

        final switch (lineLocation)
        {
            case LineLocation.first:
            case LineLocation.middle:
                output.flush();
                break;
            case LineLocation.last:
                output.writeln();
                break;
        }
    }


    void printJsonProgressLine(LineLocation lineLocation)
    {
        import vibe.data.json;

        if (lineLocation == LineLocation.first)
            return;

        auto elapsedSecs = timer.peek.total!"msecs" * 1e-3;
        auto unit = this.unit == Unit.auto_
            ? selectUnitFor(max(numTicks, totalTicks))
            : this.unit;

        auto status = Json.emptyObject;

        status["ticks"] = numTicks;
        status["elapsedSecs"] = elapsedSecs;
        status["ticksPerSec"] = cast(double) numTicks / elapsedSecs;
        if (totalTicks > 0)
            status["etaSecs"] = eta!"seconds";

        output.writeln(status.toString());
    }
}

/**
    Central logging facility for dentist.

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
    LogLevel minLevel = LogLevel.info;
}

/// Sets the minimum log level to be printed.
void setLogLevel(LogLevel level) nothrow
{
    minLevel = level;
}

LogLevel getLogLevel()
{
    return minLevel;
}

bool shouldLog(LogLevel level)
{
    return level >= minLevel;
}

/**
    Logs a message in JSON format.
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
    import std.datetime.systime : Clock;
    import std.traits : isSomeString;
    import vibe.data.json : Json;

    enum threadKey = "thread";
    enum timestampKey = "timestamp";

    if (level < minLevel)
        return;

    try
    {
        Json json = Json.emptyObject;

        static foreach (KeyValuePair; Chunks!(2, T))
        {
            static assert(isSomeString!(KeyValuePair.chunks[0]), "missing name");
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

    logJsonError("error", "mysterious observation", "secret", 42);

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
        fmt = See http://dlang.org/phobos/std_format.html#format-string
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
    string pref;
    final switch (level)
    {
    case LogLevel.debug_:
        pref = "TRACE";
        break;
    case LogLevel.diagnostic:
        pref = "DEBUG";
        break;
    case LogLevel.info:
        pref = "INFO";
        break;
    case LogLevel.warn:
        pref = "WARN";
        break;
    case LogLevel.error:
        pref = "ERROR";
        break;
    case LogLevel.fatal:
        pref = "FATAL";
        break;
    case LogLevel.none:
        assert(false);
    }
    auto threadid = () @trusted{ return cast(ulong) cast(void*) Thread.getThis(); }();
    threadid ^= threadid >> 32;

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

string traceExecution(LogLevel logLevel = LogLevel.diagnostic)()
{
    import std.conv : to;
    import std.string : replace;
    import std.traits : moduleName;

    return q"{
        import std.datetime.stopwatch : $prefix__StopWatch = StopWatch;

        $prefix__StopWatch $prefix__enter() {
            import $thisModule : LogLevel;
            import std.traits : fullyQualifiedName;
            import std.typecons : Yes;

            logJson(
                $logLevel,
                `state`, `enter`,
                `function`, $function,
            );

            return $prefix__StopWatch(Yes.autoStart);
        }

        void $prefix__exit($prefix__StopWatch timer) {
            import $thisModule : LogLevel;
            import std.traits : fullyQualifiedName;

            timer.stop();
            logJson(
                $logLevel,
                `state`, `exit`,
                `function`, $function,
                `timeElapsed`, timer.peek().total!`hnsecs`,
            );
        }

        auto $prefix__timer = $prefix__enter();

        scope (exit)
            $prefix__exit($prefix__timer);
    }"
        .replace("$thisModule", moduleName!LogLevel)
        .replace("$function", "fullyQualifiedName!(__traits(parent, $prefix__enter))")
        .replace("$logLevel", "LogLevel." ~ logLevel.to!string)
        .replace("$prefix", `__traceExecution`);
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

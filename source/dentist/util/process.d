/**
    Convenience wrappers for executing subprocesses.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.process;


import dentist.util.log;
import std.algorithm :
    endsWith,
    filter;
import std.array : array;
import std.process :
    kill,
    Redirect,
    Config,
    pipeProcess,
    pipeShell,
    ProcessPipes,
    wait;
import std.range.primitives;
import std.traits : isSomeString;
import vibe.data.json : toJson = serializeToJson;


auto pipeLines(Range)(Range command, in string workdir = null)
        if (isInputRange!Range && isSomeString!(ElementType!Range))
{
    auto sanitizedCommand = command.filter!"a != null".array;

    return new LinesPipe!ProcessInfo(ProcessInfo(sanitizedCommand, workdir));
}

auto pipeLines(in string shellCommand, in string workdir = null)
{
    return new LinesPipe!ShellInfo(ShellInfo(shellCommand, workdir));
}

unittest
{
    import std.algorithm : equal;
    import std.range : only, take;

    auto cheers = pipeLines("yes 'Cheers!'");
    assert(cheers.take(5).equal([
        "Cheers!",
        "Cheers!",
        "Cheers!",
        "Cheers!",
        "Cheers!",
    ]));

    auto helloWorld = pipeLines(only("echo", "Hello World!"));
    assert(helloWorld.equal(["Hello World!"]));
}

private struct ProcessInfo
{
    const(string[]) command;
    const(string) workdir;
}

private struct ShellInfo
{
    const(string) command;
    const(string) workdir;
}

static final class LinesPipe(CommandInfo)
{
    static enum lineTerminator = "\n";

    private CommandInfo processInfo;
    private ProcessPipes process;
    private string currentLine;

    this(CommandInfo processInfo)
    {
        this.processInfo = processInfo;
    }

    ~this()
    {
        if (!(process.pid is null))
            releaseProcess();
    }

    void releaseProcess()
    {
        if (!process.stdout.isOpen)
            return;

        process.stdout.close();

        version (Posix)
        {
            import core.sys.posix.signal : SIGKILL;

            process.pid.kill(SIGKILL);
        }
        else
        {
            static assert(0, "Only intended for use on POSIX compliant OS.");
        }
        process.pid.wait();
    }

    private void ensureInitialized()
    {
        if (!(process.pid is null))
            return;

        process = launchProcess();

        if (!empty)
            popFront();
    }

    static if (is(CommandInfo == ProcessInfo))
        ProcessPipes launchProcess()
        {
            logJsonDiagnostic(
                "action", "execute",
                "type", "pipe",
                "command", processInfo.command.toJson,
                "state", "pre",
            );

            return pipeProcess(
                processInfo.command,
                Redirect.stdout,
                null,
                Config.none,
                processInfo.workdir,
            );
        }
    else static if (is(CommandInfo == ShellInfo))
        ProcessPipes launchProcess()
        {
            logJsonDiagnostic(
                "action", "execute",
                "type", "pipe",
                "shell", processInfo.command,
                "state", "pre",
            );

            return pipeShell(
                processInfo.command,
                Redirect.stdout,
                null,
                Config.none,
                processInfo.workdir,
            );
        }

    void popFront()
    {
        ensureInitialized();
        assert(!empty, "Attempting to popFront an empty LinesPipe");
        currentLine = process.stdout.readln();

        if (currentLine.empty)
        {
            currentLine = null;
            releaseProcess();
        }

        if (currentLine.endsWith(lineTerminator))
            currentLine = currentLine[0 .. $ - lineTerminator.length];
    }

    @property string front()
    {
        ensureInitialized();
        assert(!empty, "Attempting to fetch the front of an empty LinesPipe");

        return currentLine;
    }

    @property bool empty()
    {
        ensureInitialized();

        if (!process.stdout.isOpen || process.stdout.eof)
        {
            releaseProcess();

            return true;
        }
        else
        {
            return false;
        }
    }
}

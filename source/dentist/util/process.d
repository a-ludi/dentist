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
import std.typecons : Flag, No, Yes;
import vibe.data.json : toJson = serializeToJson;


auto pipeLines(Flag!"isBuffered" isBuffered = No.isBuffered, Range)(Range command, in string workdir = null)
        if (isInputRange!Range && isSomeString!(ElementType!Range))
{
    auto sanitizedCommand = command.filter!"a != null".array;

    return new LinesPipe!(ProcessInfo, isBuffered)(ProcessInfo(sanitizedCommand, workdir));
}

auto pipeLines(Flag!"isBuffered" isBuffered = No.isBuffered)(in string shellCommand, in string workdir = null)
{
    return new LinesPipe!(ShellInfo, isBuffered)(ShellInfo(shellCommand, workdir));
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

static final class LinesPipe(CommandInfo, Flag!"isBuffered" isBuffered)
{
    static enum lineTerminator = "\n";

    static if (isBuffered)
        alias line_t = char[];
    else
        alias line_t = string;

    private CommandInfo processInfo;
    private ProcessPipes process;
    private line_t currentLine;


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

        logJsonDiagnostic(
            "action", "execute",
            "type", "pipe",
            "command", processInfo.command.toJson,
            "state", "pre",
        );

        process = launchProcess();

        if (!empty)
            popFront();
    }

    static if (is(CommandInfo == ProcessInfo))
        ProcessPipes launchProcess()
        {
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

        static if (isBuffered)
        {
            process.stdout.readln(currentLine);
        }
        else
        {
            currentLine = process.stdout.readln();
        }

        if (currentLine.length == 0)
        {
            currentLine = null;
            releaseProcess();
        }

        if (currentLine.endsWith(lineTerminator))
            currentLine = currentLine[0 .. $ - lineTerminator.length];
    }

    @property line_t front()
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

/**
    Returns true iff `name` can be executed via the process function in
    `std.process`. By default, `PATH` will be searched if `name` does not
    contain directory separators.

    Params:
        name       = Path to file or name of executable
        searchPath = Determines wether or not the path should be searched.
*/
version (Posix) bool isExecutable(scope string name, Flag!"searchPath" searchPath = Yes.searchPath)
{
    // Implementation is analogous to logic in `std.process.spawnProcessImpl`.
    import std.algorithm : any;
    import std.path : isDirSeparator;

    if (!searchPath || any!isDirSeparator(name))
        return isExecutableFile(name);
    else
        return searchPathFor(name) !is null;
}

version (Posix) unittest
{
    assert(isExecutable("/bin/sh", No.searchPath));
    assert(isExecutable("/bin/sh", Yes.searchPath));
    assert(isExecutable("sh"));
    assert(!isExecutable("does-not-exist-anywhere"));
}


version (Posix) private bool isExecutableFile(scope string path) nothrow
{
    // Implementation is analogous to private function `std.process.isExecutable`.
    import core.sys.posix.unistd : access, X_OK;
    import std.string : toStringz;

    return (access(path.toStringz(), X_OK) == 0);
}


version (Posix) private string searchPathFor(scope string executable)
{
    // Implementation is analogous to private function `std.process.searchPathFor`.
    import std.algorithm.iteration : splitter;
    import std.conv : to;
    import std.path : buildPath;
    static import core.stdc.stdlib;

    auto pathz = core.stdc.stdlib.getenv("PATH");
    if (pathz == null)  return null;

    foreach (dir; splitter(to!string(pathz), ':'))
    {
        auto execPath = buildPath(dir, executable);

        if (isExecutableFile(execPath))
            return execPath;
    }

    return null;
}

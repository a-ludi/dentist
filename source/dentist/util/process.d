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
    ProcessPipes,
    wait;
import std.range.primitives;
import std.traits : isSomeString;
import vibe.data.json : toJson = serializeToJson;


auto pipeLines(Range)(Range command, in string workdir = null)
        if (isInputRange!Range && isSomeString!(ElementType!Range))
{
    static final class LinesPipe
    {
        static enum lineTerminator = "\n";

        private const string[] command;
        private const string workdir;
        private ProcessPipes process;
        private string currentLine;

        this(in string[] command, in string workdir)
        {
            this.command = command;
            this.workdir = workdir;
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
                "command", command.toJson,
                "state", "pre",
            );
            process = pipeProcess(command, Redirect.stdout, null, Config.none, workdir);

            if (!empty)
                popFront();
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

    auto sanitizedCommand = command.filter!"a != null".array;

    return new LinesPipe(sanitizedCommand, workdir);
}

unittest
{
    import std.algorithm : equal;
    import std.range : only, take;

    auto cheers = only("yes", "Cheers!").pipeLines(".");
    assert(cheers.take(5).equal([
        "Cheers!",
        "Cheers!",
        "Cheers!",
        "Cheers!",
        "Cheers!",
    ]));

    auto helloWorld = only("echo", "Hello World!").pipeLines(".");
    assert(helloWorld.equal(["Hello World!"]));
}


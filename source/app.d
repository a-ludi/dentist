/**
    Application entry point.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module app;

import std.conv;
import std.stdio;


version (Posix)
{
    /// Start `dentist` with the given set of arguments.
    int main(string[] args)
    {
        import dentist.commandline;

        return run(args);
    }
}
else
{
    static assert(0, "not compatible with non-POSIX systems.");
}

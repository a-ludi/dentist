/**
    Some additional string functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.util.string;

import std.algorithm : joiner, map;
import std.array : array;
import std.conv : to;
import std.range : chain, cycle, take;
import std.string : lineSplitter;
import std.traits : isSomeString;

/**
    Adds one level of indentation for a multi-line string. Adds indentSize
    spaces to each non-empty line.

    Returns: indented string
*/
S indent(S)(S str, in size_t indentSize = 4) if (isSomeString!S)
{
    immutable lineSep = "\n";
    alias indentLine = (line) => chain(" ".cycle.take(line.length == 0 ? 0 : indentSize), line);

    // dfmt off
    return str[]
        .lineSplitter
        .map!indentLine
        .joiner(lineSep)
        .chain(str[$ - lineSep.length .. $] == lineSep ? lineSep : "")
        .array
        .to!S;
    // dfmt on
}

///
unittest
{
    assert("a\nb".indent == "    a\n    b");
    assert("a\nb\n\n".indent == "    a\n    b\n\n");
    assert("a\nb\n".indent(2) == "  a\n  b\n");
}

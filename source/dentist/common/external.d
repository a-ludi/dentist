/**
    This package holds function for easy verification of external tools'
    existence.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.external;

import dentist.modules : modules;
import std.array : array;
import std.algorithm :
    joiner,
    sort,
    uniq;
import std.meta : Filter, staticMap;
import std.traits : getUDAs;


struct ExternalDependency
{
    string executable;
    string package_;
    string url;

    string toString() const pure nothrow
    {
        if (package_ is null && url is null)
            return executable;
        else if (url is null)
            return executable ~ " (part of `" ~ package_ ~ "`)";
        else if (package_ is null)
            return executable ~ " (see " ~ url ~ ")";
        else
            return executable ~ " (part of `" ~ package_ ~ "`; see " ~ url ~ ")";
    }

    static enum isExternalDependency(alias value) = is(typeof(value) == ExternalDependency);

    unittest
    {
        static assert(isExternalDependency!(ExternalDependency("someTool")));
        static assert(!isExternalDependency!"someTool");
    }

    unittest
    {
        @ExternalDependency("someTool")
        void callSomeTool(in string parameter)
        {
            // calls `someTool`
        }

        static assert(fromSymbol!callSomeTool == [ExternalDependency("someTool")]);
    }
}

static enum ExternalDependency[] fromSymbol(alias symbol) = [Filter!(
    ExternalDependency.isExternalDependency,
    getUDAs!(symbol, ExternalDependency),
)];

ExternalDependency[] getExternalDependencies(Modules...)()
{
    import std.array : array;
    import std.algorithm : joiner, multiSort;
    import std.meta : staticMap;
    import std.traits : getSymbolsByUDA;

    alias _getSymbolsByUDA(alias Module) = getSymbolsByUDA!(Module, ExternalDependency);
    alias byExecutableLt = (a, b) => a.executable < b.executable;
    alias byExecutableEq = (a, b) => a.executable == b.executable;
    alias byUrlLt = (a, b) => a.url < b.url;
    alias byPackageLt = (a, b) => a.package_ < b.package_;

    return [staticMap!(
        fromSymbol,
        staticMap!(
            _getSymbolsByUDA,
            Modules,
        ),
    )]
        .joiner
        .array
        .sort!byExecutableLt
        .release
        .uniq!byExecutableEq
        .array
        .multiSort!(byPackageLt, byUrlLt, byExecutableLt)
        .release;
}

unittest
{
    struct Caller
    {
        @ExternalDependency("someTool")
        void callSomeTool(in string parameter)
        {
            // calls `someTool`
        }

        @ExternalDependency("otherTool")
        void callOtherTool(in string parameter)
        {
            // calls `otherTool`
        }
    }

    static assert(getExternalDependencies!Caller == [
        ExternalDependency("otherTool"),
        ExternalDependency("someTool"),
    ]);
}

enum externalDependencies = getExternalDependencies!modules;

/**
    This package holds function for easy verification of external tools'
    existence.

    See_also: `ExternalDependency`, `externalDependencies`
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


/// Used as a decorator to mark external dependencies. External dependencies
/// are executable that are expected to be on the `PATH`.
///
/// This decorator allows for automatic checks concerning these dependencies
/// at start up of the program rather than waiting for an error during
/// execution.
struct ExternalDependency
{
    /// Name of the executable, e.g. `LAsort`.
    string executable;

    /// Package name that is displayed to the user to aid installation of
    /// missing dependencies, e.g. `DALIGNER`.
    string package_;

    /// Package url/homepage that is displayed to the user to aid
    /// installation of missing dependencies, e.g.
    /// `https://github.com/thegenemyers/DALIGNER`.
    string url;


    /// Build a human-readable string that contains all available information.
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


    private static enum isExternalDependency(alias value) = is(typeof(value) == ExternalDependency);

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


/// Extract external dependencies declared for `symbol`.
private static enum ExternalDependency[] fromSymbol(alias symbol) = [Filter!(
    ExternalDependency.isExternalDependency,
    getUDAs!(symbol, ExternalDependency),
)];


/// Extract all external dependencies of all symbols in the listed `Modules`.
private ExternalDependency[] getExternalDependencies(Modules...)()
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


/// All external dependencies in DENTIST. Note that the actual list depends
/// on the build config (`testing` or not). You can get a valid list by
/// calling `dentist -d`.
///
/// See_also: `ExternalDependency`
enum externalDependencies = getExternalDependencies!modules;

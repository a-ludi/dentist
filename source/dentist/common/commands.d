/**
    Subcommands of dentist.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.common.commands;

import std.meta : staticMap;
import dentist.common :
    isTesting,
    testingOnly;
import dentist.util.string : dashCaseCT;


/// Subcommands of dentist.
mixin("enum DentistCommand {" ~
    "validateConfig," ~
    testingOnly!"translocateGaps," ~
    testingOnly!"buildPartialAssembly," ~
    testingOnly!"findClosableGaps," ~
    "generateDazzlerOptions," ~
    "maskRepetitiveRegions," ~
    "propagateMask," ~
    "filterMask," ~
    "mergeMasks," ~
    "showMask," ~
    "bed2mask," ~
    "chainLocalAlignments," ~
    "collectPileUps," ~
    "showPileUps," ~
    "processPileUps," ~
    "showInsertions," ~
    "mergeInsertions," ~
    "output," ~
    "translateCoords," ~
    "validateRegions," ~
    testingOnly!"checkResults," ~
"}");

/// Helper for subcommands that are only available in testing version.
struct TestingCommand
{
    @disable this();

    static DentistCommand opDispatch(string command)() pure nothrow
    {
        static if (isTesting)
            return mixin("DentistCommand." ~ command);
        else
            return cast(DentistCommand) size_t.max;
    }
}

/// Tuple of dentist subcommand names.
enum dentistCommands = staticMap!(
    dashCaseCT,
    __traits(allMembers, DentistCommand),
);

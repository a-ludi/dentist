/**
    This module contains information about this software.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.swinfo;

import dentist.common : testingOnly;
import std.string : wrap;

enum executableName = "dentist";
enum gitVersion = "v4.0.0-60-ga685cf0-dirty";
enum gitCommit = "a685cf02d6e8de4b805e9f00b5ef416d7c30b376+dirty";
enum version_ = gitVersion ~ testingOnly!"+testing";
enum description = "Close assembly gaps using long-reads with focus on correctness.".wrap;
enum copyright = "Copyright © 2018, Arne Ludwig <arne.ludwig@posteo.de>".wrap;
enum license = q"{
    Subject to the terms of the MIT license, as written in the included
    LICENSE file
}".wrap;

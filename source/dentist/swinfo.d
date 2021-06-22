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
enum gitVersion = "v2.0.0";
enum gitCommit = "5785bb3913734394cc1cd1233713ac7361def180";
enum version_ = gitVersion ~ testingOnly!"+testing";
enum description = "Close assembly gaps using long-reads with focus on correctness.".wrap;
enum copyright = "Copyright © 2018, Arne Ludwig <arne.ludwig@posteo.de>".wrap;
enum license = q"{
    Subject to the terms of the MIT license, as written in the included
    LICENSE file
}".wrap;

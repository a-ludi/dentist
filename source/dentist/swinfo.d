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
enum gitCommit = "68d2a8ddad732c0309f6c4376d0d1582487c9a20";
enum version_ = gitVersion ~ testingOnly!"+testing";
enum description = "Close assembly gaps using long-reads with focus on correctness.".wrap;
enum copyright = "Copyright © 2018, Arne Ludwig <arne.ludwig@posteo.de>".wrap;
enum license = q"{
    Subject to the terms of the MIT license, as written in the included
    LICENSE file
}".wrap;

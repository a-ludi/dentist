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
enum gitVersion = "v1.0.1-30-gfe489f1-dirty";
enum gitCommit = "fe489f11425a70ab50d15dba80929d10e7c1b54c+dirty";
enum version_ = gitVersion ~ testingOnly!"+testing";
enum description = "Close assembly gaps using long-reads with focus on correctness.".wrap;
enum copyright = "Copyright © 2018, Arne Ludwig <arne.ludwig@posteo.de>".wrap;
enum license = q"{
    Subject to the terms of the MIT license, as written in the included
    LICENSE file
}".wrap;

/**
    This is the `validateConfig` command of `dentist`.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.validateConfig;

import dentist.commandline : OptionsFor;
import dentist.common.commands : DentistCommand;
import dentist.common.configfile : validateConfigFile;
import std.stdio : writefln;
import vibe.data.json :
    toJson = serializeToJson,
    toJsonString = serializeToPrettyJson;


/// Options for the `translateCoords` command.
alias Options = OptionsFor!(DentistCommand.validateConfig);

/// Execute the `translateCoords` command with `options`.
void execute(in Options options)
{
    validateConfigFile(options.configFile);
    writefln!"Config `%s` valid!"(options.configFile);
}

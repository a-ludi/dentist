/**
    This is the `validate-config` command of DENTIST.

    Command_Summary:

    ---
    Validate config file. Exit with non-zero status and a descriptive error
    message if errors are found.
    ---

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist.commands.validateConfig;

package(dentist) enum summary = "
    Validate config file. Exit with non-zero status and a descriptive error
    message if errors are found.
";

import dentist.commandline : OptionsFor;
import dentist.common.commands : DentistCommand;
import dentist.common.configfile : validateConfigFile;
import std.stdio : writefln;
import vibe.data.json :
    toJson = serializeToJson,
    toJsonString = serializeToPrettyJson;


/// Options for the `validate-config` command.
alias Options = OptionsFor!(DentistCommand.validateConfig);


/// Execute the `validate-config` command with `options`.
void execute(in Options options)
{
    validateConfigFile(options.configFile);
    writefln!"Config `%s` valid!"(options.configFile);
}

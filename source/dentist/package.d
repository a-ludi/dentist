/**
    Exposes the `run` function to execute DENTIST commands. It takes an array
    of CLI arguments to select the command and construct its options.

    See_also: `dentist.commandline.run`, `dentist.commands`
    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module dentist;

public import dentist.commandline : run;

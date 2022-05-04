#!/usr/bin/env -S jq --slurp -f

## Usage: debug-too-many-open-files.jq <json:log>
##
## Generate statistics about spawned processes from a log file.
##
## Positional arguments:
##  <json:log>   Some log file with appropriate log entries (requires version(Issue31_TooManyOpenFiles))


def getCounts(what):
    reduce .[] as $call ({};
        ($call | what) as $key |
        if has($key) then
            .[$key] += 1
        else
            .[$key] = 1
        end
    );

def caller:
    .stackTrace |
    map(split(" dentist.")[1]) |
    until(.[0] | startswith("util.") | not; .[1:]) |
    "dentist." + .[0] |
    split("(")[0];



map(select(has("stackTrace"))) |
{
  "perFunction": getCounts(.function),
  "perExecutable": getCounts(.command[0]),
  "perCaller": getCounts(caller)
}

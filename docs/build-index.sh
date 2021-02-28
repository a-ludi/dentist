#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'
PROG="$(basename "$0")"
SCRIPTDIR="$(dirname "$0")"


function set_defaults()
{
    README="$(realpath "$SCRIPTDIR/../README.md")"
    TARGET="$(realpath "$SCRIPTDIR/index.md")"
}


function error()
{
    if (( $# > 0 ));
    then
        echo "$PROG: error:" "$@"
    else
        echo
    fi
} >&2


function bail_out()
{
    error "$@"
    exit 1
}


function bail_out_usage()
{
    error "$@"
    error
    print_usage
    exit 1
}


function log()
{
    echo "--" "$@"
} >&2


function print_usage()
{
    echo "USAGE:  $PROG [-h] [<README>]"
} >&2


function print_help()
{
    print_usage
    echo
    echo 'Builds index.md from projects README.md.'
    echo
    echo 'Optional arguments:'
    echo ' --dry-run, -n   Just print the result; do not alter target file.'
    echo ' --help, -h      Prints this help.'
    echo ' --usage         Print a short command summary.'
    echo ' --version       Print software version.'
} >&2


function print_version()
{
    echo "$PROG v0"
    echo
    echo "Copyright © 2019, Arne Ludwig <arne.ludwig@posteo.de>"
} >&2


function parse_args()
{
    ARGS=()
    for ARG in "$@";
    do
        if [[ "${ARG:0:1}" == - ]];
        then
            case "$ARG" in
                -n|--dry-run)
                    DRY_RUN=1
                    ;;
                -h|--help)
                    print_help

                    exit
                    ;;
                --usage)
                    print_usage

                    exit
                    ;;
                --version)
                    print_version

                    exit
                    ;;
                *)
                    bail_out_usage "unkown option $ARG"
                    ;;
            esac
        else
            ARGS+=( "$ARG" )
        fi
    done

    (( ${#ARGS[*]} == 1 )) && README="${ARGS[0]}"
    (( ${#ARGS[*]} <= 1 )) || bail_out_usage "too many arguments"
}


function main()
{
    set_defaults
    parse_args "$@"

    if [[ -v DRY_RUN ]]
    then
        TARGET=/dev/stdout
    fi

    awk '
        (!print_all_lines && $0 ~ /docs\/logo\.png|^>/) {
            skip_lines = 2;
        }

        ($0 == "Table of Contents") {
            print_all_lines = 1;
        }

        {
            if (!print_all_lines && skip_lines > 0) {
                --skip_lines;
            } else {
                print
            }
        }
    ' "$README" > "$TARGET"
}


main "$@"

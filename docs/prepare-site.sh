#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'
PROG="$(basename "$0")"
SCRIPTDIR="$(dirname "$0")"


function set_defaults()
{
    README="$(realpath "$SCRIPTDIR/../README.md")"
    INDEX="$(realpath "$SCRIPTDIR/index.md")"
    LICENSE="$(realpath "$SCRIPTDIR/../LICENSE")"
    LICENSE_MD="$(realpath "$SCRIPTDIR/license.md")"
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


function remove_logo()
{
    sed -E '/docs\/logo\.png/ { N; d }'
}


function remove_short_description()
{
    sed -E '
        1, /Table of Contents/ {
            /^>.*DENTIST.*/ { N; d }
        }
    '
}


function remove_table_of_contents()
{
    sed -E '
        /^(#\s*)?Table of Contents$/,/^(#\s*)?Install$/ {
            /^(#\s*)?Install$/ ! d
        }
    '
}


function adjust_relative_paths()
{
    sed -E '
        s,\./docs(/.+)\.md,%DOT%\1.html,g
        s,\./LICENSE,%DOT%/license.html,g
        s,(`[^`])*\.(/[^`]*`),\1%DOT%\2,g
        s,\./,https://github.com/a-ludi/dentist/blob/develop/,g
        s,%DOT%,.,g
    '
}


function add_license_header()
{
    echo 'License'
    echo '======='
    echo
    cat
}


function adjust_email_links()
{
    sed -E 's/<[^>]+>/\&lt;&\&gt;/'
}


function list_sub_pages()
{
    find "$SCRIPTDIR" -maxdepth 1 -type f -name '*.md' -not -name 'index.md'
}

function add_back_link()
{
    sed -i -E '
        /^={3,}/ {
            n;
            /Back to homepage/ ! {
                i \
[&larr; Back to homepage](/)\

                :print-rest
                n
                b print-rest
            }
        }
    ' "$@"
}


function main()
{
    set_defaults
    parse_args "$@"

    if [[ -v DRY_RUN ]]
    then
        INDEX=/dev/stdout
        LICENSE_MD=/dev/stdout
    fi

    cat "$README" \
    | remove_logo \
    | remove_short_description \
    | remove_table_of_contents \
    | adjust_relative_paths \
    > "$INDEX"

    cat "$LICENSE" \
    | add_license_header \
    | adjust_email_links \
    > "$LICENSE_MD"

    [[ -v DRY_RUN ]] || add_back_link $(list_sub_pages)
}


main "$@"

#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'
PROG="$(basename "$0")"


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
    echo "USAGE:  $PROG [-h]"
} >&2


function print_help()
{
    print_usage
    echo
    echo 'Build DENTIST in release mode and create a tarball.'
    echo
    echo 'Optional arguments:'
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

    (( ${#ARGS[*]} == 0 )) || bail_out_usage "too many arguments"
}


function prepare_dist()
{
    [[ ! -e "$DIST_DIR" ]] || bail_out "could not create dist directory: file already exists: $DIST_DIR"

    mkdir "$DIST_DIR"
    cp dentist README.md CHANGELOG.md LICENSE "$DIST_DIR"

    mkdir "$DIST_DIR/snakemake"
    cp snakemake/cluster.yml "$DIST_DIR/snakemake/cluster.example.yml"
    cp snakemake/snakemake.yml "$DIST_DIR/snakemake/snakemake.example.yml"
    cp snakemake/{profile-slurm.yml,Snakefile,workflow_helper.py} "$DIST_DIR/snakemake"
}


function main()
{
    DENTIST=dentist

    parse_args "$@"

    dub build --build=release
    strip -s "$DENTIST"

    DENTIST_VERSION="$("$DENTIST" --version |& head -n1)"
    DENTIST_VERSION="${DENTIST_VERSION#dentist }"
    ARCH="$(uname -m)"
    TARBALL="dentist.$DENTIST_VERSION.$ARCH.tar.gz"
    DIST_DIR="dentist.$DENTIST_VERSION.$ARCH"

    prepare_dist
    tar -czf "$TARBALL" "$DIST_DIR"

    log "created $TARBALL"
}


main "$@"

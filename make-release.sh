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


function make_tarball()
{
    DENTIST_VERSION="$(dentist --version 2>&1 | head -n1)" && \
    DENTIST_VERSION="${DENTIST_VERSION#dentist }" && \
    DENTIST_VERSION="${DENTIST_VERSION% (*)}" && \
    DENTIST_SRC="/opt/dentist" && \
    ARCH="$(uname -m)" && \
    TARBALL="dentist.$DENTIST_VERSION.$ARCH.tar.gz" && \
    DIST_DIR="dentist.$DENTIST_VERSION.$ARCH" && \
    INCLUDE_BINARIES=( Catrack DAM2fasta DAScover DASqv DB2fasta DBa2b DBb2a
        DBdump DBdust DBmv DBrm DBshow DBsplit DBstats DBtrim DBwipe LAa2b
        LAb2a LAcat LAcheck LAdump LAmerge LAshow LAsort LAsplit TANmask
        computeintrinsicqv daccord daligner damapper datander dentist dumpLA
        fasta2DAM fasta2DB lasfilteralignments rangen simulator )
    for (( I = ${#INCLUDE_BINARIES[*]} - 1; I >= 0; --I ))
    do
        INCLUDE_BINARIES[$I]="$BINDIR/${INCLUDE_BINARIES[$I]}"
    done

    cd /tmp && \
    install -Dt "$DIST_DIR" \
        "$DENTIST_SRC/README.md" \
        "$DENTIST_SRC/CHANGELOG.md" \
        "$DENTIST_SRC/LICENSE" && \
    install -Dt "$DIST_DIR/snakemake" \
        "$DENTIST_SRC/snakemake/cluster.yml" \
        "$DENTIST_SRC/snakemake/snakemake.yml" \
        "$DENTIST_SRC/snakemake/profile-slurm."*".yml" \
        "$DENTIST_SRC/snakemake/Snakefile" && \
    install -Dt "$DIST_DIR/bin" "${INCLUDE_BINARIES[@]}" && \
    tar --remove-files -czf "$TARBALL" "$DIST_DIR" && \
    echo "tarball:$TARBALL"
    realpath "$TARBALL"
}


function main()
{
    parse_args "$@"

    log "building conatiner image"
    trap 'rm -f .docker-build-id' exit
    docker build --iidfile .docker-build-id .

    (
        log "gathering release files"
        CONTAINER_ID=$(docker create -v "$PWD:/opt/dentist:ro" "$(< .docker-build-id)" bash -c "$(declare -f make_tarball); make_tarball") && \
        trap 'docker rm $CONTAINER_ID' exit && \

        TARBALL="$(docker start -a "$CONTAINER_ID" | tee /dev/stderr | grep -E '^tarball:')"
        TARBALL="${TARBALL#tarball:}"

        log "copying tarball $TARBALL"
        docker cp "$CONTAINER_ID:/tmp/$TARBALL" ./

        log "created $TARBALL"
    )
}


main "$@"

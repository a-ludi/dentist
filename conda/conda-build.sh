#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'
PROG="$(basename "$0")"


function set_defaults()
{
    ARTIFACTS_DIR=_artifacts
    BUILD_IMAGE=conda-build
    CONDA_BUILD_ARGS=(
        -c a_ludi
        -c bioconda
    )
    CONTAINER_ID=''
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
    echo -ne '\e[1m'
    echo -n "--" "$@"
    echo -e '\e[0m'
} >&2


function print_usage()
{
    echo "USAGE:  $PROG [-h] <recipes>"
} >&2


function print_help()
{
    print_usage
    echo
    echo 'This program shows a helpful message.'
    echo
    echo 'Positional arguments:'
    echo ' <recipe>        Path to conda recipe.'
    echo
    echo 'Optional arguments:'
    echo ' --artifacts=<dir>  Place build artifacts under <dir> '"(default: $ARTIFACTS_DIR)."
    echo ' --debug            If given make a snapshot of the build image for debugging.'
    echo ' --upload, -u       Automatically upload package to Anaconda.'
    echo ' --help, -h         Prints this help.'
    echo ' --usage            Print a short command summary.'
    echo ' --version          Print software version.'
} >&2


function print_version()
{
    echo "$PROG v0"
    echo
    echo "Copyright © 2021, Arne Ludwig <arne.ludwig@posteo.de>"
} >&2


function parse_args()
{
    ARGS=()
    for ARG in "$@";
    do
        if [[ "${ARG:0:1}" == - ]];
        then
            case "$ARG" in
                --artifacts=*)
                    ARTIFACTS_DIR=${ARG#--artifacts=}
                    ;;
                --debug)
                    DEBUG=1
                    ;;
                -u|--upload)
                    DO_UPLOAD=1
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

    (( ${#ARGS[*]} == 1 )) || bail_out_usage "<recipe> is missing"

    RECIPE="${ARGS[0]}"
}


function clean_up()
{
    rm -f .docker-build-id

    if [[ -n "$CONTAINER_ID" ]]
    then
        if [[ -v DEBUG ]]
        then

            log "creating build image snapshot"
            SNAPSHOT_ID="$(docker container commit -c 'ENTRYPOINT ["/bin/bash"]' "$CONTAINER_ID")"
            SNAPSHOT_ID="${SNAPSHOT_ID#sha256:}"
            log "created snapshot image $SNAPSHOT_ID"
        fi

        docker rm "$CONTAINER_ID" > /dev/null
    fi
}


function main()
{
    set_defaults
    parse_args "$@"

    trap 'clean_up' exit

    log "building container image"
    docker build -t "$BUILD_IMAGE" conda-build

    log "gathering release files"
    CONTAINER_ID="$(docker create -v "$(realpath "$RECIPE"):/recipe:ro" "$BUILD_IMAGE" "${CONDA_BUILD_ARGS[@]}")"

    log "building conda package ..."
    # builds conda package and sets TARBALL
    ENV_EXPORT="$(docker start -a "$CONTAINER_ID")"
    eval "$ENV_EXPORT"
    LOCAL_TARBALL="$ARTIFACTS_DIR/${TARBALL#/opt/miniconda3/conda-bld/}"

    log "copying tarball $TARBALL -> $LOCAL_TARBALL"
    mkdir -p "$(dirname "$LOCAL_TARBALL")"
    docker cp "$CONTAINER_ID:$TARBALL" "$LOCAL_TARBALL"

    if [[ -v DO_UPLOAD ]]
    then
        anaconda upload "$LOCAL_TARBALL"
    else
        log "upload package with:"
        echo -e "\tanaconda upload ${LOCAL_TARBALL@Q}"
    fi
}


main "$@"

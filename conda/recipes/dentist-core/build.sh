#!/usr/bin/env bash

# Copyright: © 2021 Arne Ludwig <arne.ludwig@posteo.de>
# License:   Subject to the terms of the MIT license, as written in the
#            included LICENSE file.
# Authors:   Arne Ludwig <arne.ludwig@posteo.de>


# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail -x
IFS="$(printf '\n\t')"

PROG="$0"
export PREFIX="${PREFIX:-/usr/bin}"
export BINDIR="$PREFIX/bin"

DENTIST_BUILD="${DENTIST_BUILD:-release}"
DENTIST_BUILD_CONFIG="${DENTIST_BUILD_CONFIG:-default}"
DMD_VERSION="${DMD_VERSION:-"2.100.2"}"

DLANG_ROOT="$PWD/.dlang"


log()
{
    echo -ne '\e[1m'
    echo -n "--" "$@"
    echo -e '\e[0m'
} >&2


bail_out()
{
    echo "$PROG: error: $*"
    exit 1
} >&2


main()
{
    [[ -v PREFIX ]] || bail_out "PREFIX environment variable is undefined"
    [[ -z "${PATH##*$BINDIR*}" ]] || PATH="$BINDIR${PATH+:}$PATH"
    [[ "$(uname)" == "Linux" || "$(uname)" == "Darwin" ]] || bail_out "only Linux or OSX platform supported"

    export BINDIR
    mkdir -p "$BINDIR"

    install_dmd

    log "building unit tests..."
    dub build --config="$DENTIST_BUILD_CONFIG" --build=unittest
    mv dentist dentist-unittest

    log "building main executable..."
    dub build --config="$DENTIST_BUILD_CONFIG" --build="$DENTIST_BUILD"

    # clean up build artifacts
    rm -rf .dub
    [[ ! -d "$DLANG_ROOT" ]] || rm -rf "$DLANG_ROOT"

    log "installing binaries..."
    install -t "$BINDIR" dentist dentist-unittest
}


install_dmd()
{
    log "installing DMD $DMD_VERSION..."
    dlang_install install "dmd-$DMD_VERSION"

    log "activating DMD $DMD_VERSION..."
    # Setup PATH, LIBRARY_PATH, LD_LIBRARY_PATH, DMD, DC, and PS1
    # shellcheck disable=SC1090
    source "$(dlang_install "dmd-$DMD_VERSION" -a)"

    export DUB
    DUB="$(dlang_install get-path "dmd-$DMD_VERSION" --dub)"
}


gpg_major_version()
{
    local GPG_VERSION
    GPG_VERSION="$(gpg --version | head -n1)"
    GPG_VERSION="${GPG_VERSION##* }"

    echo "${GPG_VERSION%%.*}"
}


dlang_install()
{
    if ! [[ -e "$DLANG_ROOT/install.sh" ]]
    then
        mkdir -p "$DLANG_ROOT"
        curl https://dlang.org/install.sh > "$DLANG_ROOT/install.sh"
        chmod +x "$DLANG_ROOT/install.sh"

        if (( $(gpg_major_version) < 2 ))
        then
            # NOTE: older GnuPG versions do not support latest trends in
            #       signature algos. This, we sacrifice verification for
            #       for a working installation.
            sed -i -E 's/^verify\(\)\s*\{/\0\n    return/' .dlang/install.sh
        fi
    fi

    "$DLANG_ROOT/install.sh" -p "$DLANG_ROOT" "$@"
}


main "$@"

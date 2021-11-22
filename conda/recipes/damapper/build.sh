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

    # Build parts
    make
    # Install
    make DEST_DIR="$BINDIR" install
}

main "$@"

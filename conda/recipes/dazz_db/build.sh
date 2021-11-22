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
BUILD_PREFIX="${BUILD_PREFIX:-$PREFIX}"
export BUILD_LIBDIR="$BUILD_PREFIX/lib"


bail_out()
{
    echo "$PROG: error: $*"
    exit 1
} >&2


main()
{
    [[ -v PREFIX ]] || bail_out "PREFIX environment variable is undefined"
    [[ -z "${PATH##*$BINDIR*}" ]] || PATH="$BINDIR${PATH+:}${PATH:-}"
    [[ "$(uname)" == "Linux" || "$(uname)" == "Darwin" ]] || bail_out "only Linux or OSX platform supported"

    export BINDIR
    mkdir -p "$BINDIR"

    # Use compiler that was supplied by conda
    sed -i -E 's/^\tgcc\b/\t$(CC)/' Makefile

    # Build parts
    make CFLAGS="-O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing -Wl,-L$BUILD_LIBDIR"
    # Install
    make DEST_DIR="$BINDIR" install
}

main "$@"

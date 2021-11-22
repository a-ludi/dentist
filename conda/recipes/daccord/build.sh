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

    build_install_libmaus2
    build_install_daccord
    make -j2 -C libmaus2 uninstall
}


build_install_libmaus2()
(
    cd libmaus2
    autoupdate
    autoreconf -i -f
    ./configure "--prefix=${PREFIX}" --with-gmp
    make -j2
    make install
)


build_install_daccord()
(
    cd daccord
    autoreconf -i -f
    # statically link libamus2
    export MAKEFLAGS='LDFLAGS=-static-libtool-libs'
    # statically link libunwind
    rm -v "$BUILD_PREFIX"/lib/libunwind*so*
    ./configure "--prefix=${PREFIX}"
    make
    make install
)


main "$@"

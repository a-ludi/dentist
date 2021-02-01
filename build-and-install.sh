#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'

function fail()
{
    echo "error:" "$@"

    exit 1
} >&2

(( $# >= 1 )) || fail 'missing argument <software>'
SOFTWARE="$1"

(( $# >= 2 )) || fail 'missing argument <toolchain>'
TOOLCHAIN="$2"

[[ -v REPO ]] || fail 'missing variable REPO'

BRANCH="${BRANCH:-}"
BUILD="${BUILD:-release}"
PREBUILD="${PREBUILD:-}"
INSTALL_CMD="${INSTALL_CMD:-install -t \"\$BINDIR\" \$(find . -maxdepth 1 -executable -not -type d)}"
BINDIR="${BINDIR:-/usr/local/bin}"
BUILDDIR="${BUILDDIR:-/opt}"
NCPUS="${NCPUS:-1}"
(( NCPUS >= 1 )) || fail 'NCPUS must be integer >= 1'
NCPUS="$(( NCPUS ))"


cd "$BUILDDIR"
[[ -d "$SOFTWARE" ]] || git clone "$REPO" "$SOFTWARE"
cd "$SOFTWARE"
[[ -z "$BRANCH" ]] || git checkout "$BRANCH"
if [[ -f .gitmodules ]]
then
    git submodule init
    git submodule update --recursive --jobs=$NCPUS
fi
[[ -z "$PREBUILD" ]] || eval "$PREBUILD"
if [[ "$TOOLCHAIN" == make ]]
then
    make -j$NCPUS
elif [[ "$TOOLCHAIN" == dub ]]
then
    dub build ${BUILD_CONFIG:+"--config=$BUILD_CONFIG"} --build="$BUILD" $( (( NCPUS == 1 )) || echo '--parallel' )
else
    echo "invalid toolchain: $TOOLCHAIN" >&2
    exit 10
fi
eval "$INSTALL_CMD"
cd "$BUILDDIR"
rm -rf "$SOFTWARE"

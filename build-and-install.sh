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

function mkabsdir() {
    mkdir -p "$1"
    realpath "$1"
}

BRANCH="${BRANCH:-}"
BUILD="${BUILD:-release}"
PREBUILD="${PREBUILD:-}"
INSTALL_CMD="${INSTALL_CMD:-install -t \"\$BINDIR\" \$(find . -maxdepth 1 -executable -not -type d)}"
CLEANBUILD="${CLEANBUILD:-1}"
BINDIR="$(mkabsdir "${BINDIR:-/usr/local/bin}")"
BUILDDIR="$(mkabsdir "${BUILDDIR:-/opt}")"
SOURCEDIR="$(realpath "$BUILDDIR/$SOFTWARE")"
NCPUS="${NCPUS:-1}"
(( NCPUS >= 1 )) || fail 'NCPUS must be integer >= 1'
NCPUS="$(( NCPUS ))"


[[ -d "$SOURCEDIR" ]] || git clone "$REPO" "$SOURCEDIR"

pushd "$SOURCEDIR"

if [[ -n "$BRANCH" ]] && ! git diff --quiet HEAD "$BRANCH"
then
    git checkout "$BRANCH"
fi

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

popd

(( CLEANBUILD == 0 )) || rm -rf "$SOURCEDIR"

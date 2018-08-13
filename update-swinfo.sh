#!/bin/bash

echo -n 'Updating `dentist.swinfo` ... ' >&2

trap 'echo failed >&2' ERR
set -e  # exit when any command fails

INFO_JSON="$(dub describe | jq -c '.packages[0]')"

EXECUTABLE_NAME="$(jq '.targetName' <<<"$INFO_JSON")"
VERSION="$(jq '.version' <<<"$INFO_JSON")"
DESCRIPTION="$(jq '.description' <<<"$INFO_JSON")"
COPYRIGHT="$(jq '.copyright' <<<"$INFO_JSON")"

sed -Ei \
    -e 's/(executableName\s*=\s*)[^;]+;/\1'"$EXECUTABLE_NAME"';/' \
    -e 's/(version_\s*=\s*)[^;]+;/\1'"$VERSION"';/' \
    -e 's/(description\s*=\s*)[^;]+;/\1'"$DESCRIPTION"'.wrap;/' \
    -e 's/(copyright\s*=\s*)[^;]+;/\1'"$COPYRIGHT"'.wrap;/' \
    source/dentist/swinfo.d

echo 'done' >&2

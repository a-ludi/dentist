#!/bin/bash

SWINFO_FILE=source/dentist/swinfo.d

function main()
{
    echo -n 'Updating `dentist.swinfo` ... ' >&2

    get_updated_swinfo "$SWINFO_FILE" > "$SWINFO_FILE~"

    if ! cmp --quiet "$SWINFO_FILE" "$SWINFO_FILE~";
    then
        mv "$SWINFO_FILE~" "$SWINFO_FILE"

        echo 'done' >&2
    else
        echo 'skipped' >&2
    fi
}

function is_git_dirty()
{
    git status --porcelain | grep -qv '\bsource/dentist/swinfo.d$'
}

function get_updated_swinfo()
{
    INFO_JSON="$(dub describe | jq -c '.packages[0]')"

    EXECUTABLE_NAME="$(jq '.targetName' <<<"$INFO_JSON")"
    GIT_VERSION="$(git describe --dirty)"
    GIT_BRANCH="$(git branch --show-current)"
    GIT_COMMIT="$(git rev-parse HEAD)"
    DESCRIPTION="$(jq '.description' <<<"$INFO_JSON")"
    COPYRIGHT="$(jq '.copyright' <<<"$INFO_JSON")"

    if [[ -z "${GIT_BRANCH##release/*}" ]]
    then
        GIT_VERSION="$(git config --get gitflow.prefix.versiontag)${GIT_BRANCH#release/}"
    fi

    if is_git_dirty
    then
        GIT_COMMIT="$GIT_COMMIT+dirty"
    fi

    sed -E \
        -e 's/(executableName\s*=\s*)[^;]+;/\1'"$EXECUTABLE_NAME"';/' \
        -e 's/(gitVersion\s*=\s*)[^;]+;/\1"'"$GIT_VERSION"'";/' \
        -e 's/(gitCommit\s*=\s*)[^;]+;/\1"'"$GIT_COMMIT"'";/' \
        -e 's/(description\s*=\s*)[^;]+;/\1'"$DESCRIPTION"'.wrap;/' \
        -e 's/(copyright\s*=\s*)[^;]+;/\1'"$COPYRIGHT"'.wrap;/' \
        "$1"
}

function clean_up()
{
    rm -f "$SWINFO_FILE~"
}

function on_error()
{
    echo failed >&2
}

trap clean_up EXIT
trap on_error ERR
set -e  # exit when any command fails

main

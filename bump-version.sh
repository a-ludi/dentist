#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'
PROG="$(basename "$0")"


function bail_out()
{
    echo "$PROG: error: $@"
    exit 1
} >&2


function cancel_by_user()
{
    local CLEANUP="${1:-true}"

    echo "$PROG: process cancelled by user"
    eval "$CLEANUP"
    exit 1
} >&2


function get_release_version()
{
    VERSION="$(sed -nE '/^##\s/ {s/^##\s+\[([^]]+)\].*/\1/p ; q}' CHANGELOG.md)"

    [[ "$VERSION" == "Unreleased" ]] && bail_out "update the Changelog first"
    grep -qE '^[0-9]+\.[0-9]+\.[0-9]+(-(alpha|beta)\.[0-9]+)?$' <<<"$VERSION" || \
        bail_out "ill-formatted version not allowed: $VERSION"
}


function get_release_changelog()
{
    CHANGELOG="$(sed -nE '
        /^##\s/ {
            h; # put release title in hold buffer

        :collect
            n; # move to next line

            /^##\s/ { # found next release notes
                x; # get current release notes
                s/\s+$//; # strip trailing whitespace
                p; # print out
                q; # stop
            }
            H; # append line to hold buffer
            b collect; # loop
        }
    ' CHANGELOG.md)"
}


function remove_existing_release_tag()
{
    if git rev-parse --quiet --verify "$RELEASE_TAG" > /dev/null
    then
        git tag -d "$RELEASE_TAG"
    fi
}


function create_release_tag()
{
    git tag "$@" --annotate --cleanup=verbatim "$RELEASE_TAG" --message="$CHANGELOG"
}


function show_release_tag()
{
    git show --format=short --no-patch "$RELEASE_TAG"
}


function user_confirm()
{
    local CLEANUP="${2:-true}"
    local REPLY
    read -rp "$1 (y/N)"

    [[ "$REPLY" == "y" ]] || cancel_by_user "$CLEANUP"
}


function validate_swinfo()
{
    grep -F "enum gitVersion = \"$RELEASE_TAG\";" source/dentist/swinfo.d || \
        bail_out "could not find the correct version number in " \
                 "source/dentist/swinfo.d; is the worktree clean? "\
                 "(git status should be empty)"
}


function main()
{
    get_release_version
    RELEASE_TAG="v$VERSION"
    get_release_changelog

    # create preliminary tag
    remove_existing_release_tag
    create_release_tag
    # validate the version and release notes
    show_release_tag
    user_confirm "Are the version number and changelog correct?" \
        'git tag -d "$RELEASE_TAG"'

    # update version info
    ./update-swinfo.sh
    # validate correct version without `+dirty
    validate_swinfo

    # commit swinfo.d
    git commit source/dentist/swinfo.d -m 'Bump version'

    # update release tag
    create_release_tag -f
    show_release_tag
}

main "$@"

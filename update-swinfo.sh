#!/bin/bash

SWINFO_FILE=source/dentist/swinfo.d

function main()
{
    echo -n 'Updating `dentist.swinfo` ... ' >&2

    get_updated_swinfo "$SWINFO_FILE" > "$SWINFO_FILE~"

    if ! cmp -s "$SWINFO_FILE" "$SWINFO_FILE~";
    then
        mv "$SWINFO_FILE~" "$SWINFO_FILE"

        echo 'done' >&2
    else
        echo 'skipped' >&2
    fi
}

function is_git_clean()
{
    git status --porcelain | \
        awk '
            BEGIN {
                ignore["source/dentist/swinfo.d"] = 1;
                ignore["source/dentist/swinfo.d~"] = 1;
                ignore["update-swinfo.sh"] = 1;

                list_dirty_files = ("LIST_DIRTY_FILES" in ENVIRON && ENVIRON["LIST_DIRTY_FILES"] != 0);
            }

            (FILENAME == ".dockerignore" && substr($0, 1, 1) == "!") {
                current = substr($0, 2);
                gsub(/\*/, "[^/]*", current);
                gsub(/\?/, "[^/]?", current);
                gsub(/\./, "\\.", current);
                current = "^" current "(/|$)";

                if (include_re)
                    include_re = include_re "|" current;
                else
                    include_re = current;
            }

            (FILENAME == "-") {
                if ($2 in ignore || !($2 ~ include_re)) {
                    # ignore
                } else {
                    # there is a change to a file that is not ignored
                    if (list_dirty_files)
                        dirty_files[++n] = $0;
                    else
                        exit ++n;
                }
            }

            END {
                if (list_dirty_files) {
                    print "include_re" "=" include_re > "/dev/stderr";
                    for (i = 1; i <= n; ++i)
                        print dirty_files[i] > "/dev/stderr";
                }

                exit n
            }
        ' .dockerignore -
}

function get_updated_swinfo()
{
    INFO_JSON="$(dub describe | jq -c '.packages[0]')"

    EXECUTABLE_NAME="$(jq '.targetName' <<<"$INFO_JSON")"
    GIT_VERSION="$(git describe)"
    GIT_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
    GIT_COMMIT="$(git rev-parse HEAD)"
    DESCRIPTION="$(jq '.description' <<<"$INFO_JSON")"
    COPYRIGHT="$(jq '.copyright' <<<"$INFO_JSON")"

    if [[ -z "${GIT_BRANCH##release/*}" ]]
    then
        GIT_VERSION="$(git config --get gitflow.prefix.versiontag)${GIT_BRANCH#release/}"
    fi

    if ! is_git_clean
    then
        GIT_VERSION="$GIT_VERSION-dirty"
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

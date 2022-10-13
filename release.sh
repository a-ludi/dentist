#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'
PROG="$(basename "$0")"

STATUS_FILE="$PWD/.release-status"
GUI_EDITOR=subl
SUDO='pkexec --keep-cwd'
DENTIST="${DENTIST:+$PWD}"
GH_PAGES_ROOT="$(realpath "$PWD/../dentist-gh-pages")"
GH_PAGES_BRANCH=gh-pages
PUBLIC_GIT=public


function bail_out()
{
    echo -en "\e[31mError: "
    echo -n "$@"
    echo -e "\e[0m"
    exit 1
} >&2


function sudo()
{
    eval "$SUDO $@"
}


function gui_editor()
{
    eval "$GUI_EDITOR $@"
}


function p()
{
    echo "$@"
}


function part()
{
    if [[ -v NUM_PARTS ]]
    then
        (( ++NUM_PARTS ))
    else
        NUM_PARTS=1
    fi
    echo
    echo
    echo -en "\e[1m"
    echo -n "Part $NUM_PARTS:" "$@"
    echo -e "\e[0m"
    echo
}


function joinstr()
{
    OPTIND=1
    while getopts "s:" OPTION; do
        case "$OPTION" in
            s)
                local SEP="$OPTARG"
                ;;
            *)
                echo "$(basename "$0"): joinstr: unkown option $OPTION" >&2
                ;;
        esac
    done
    shift $(( $OPTIND - 1 ))
    local IFS="${SEP:- }"

    echo -n "$*"
}


function git_commit()
{
    git rev-parse HEAD
}


function task_signature()
{
    eval "$(grep -E '^NUM_TASKS=[0-9]+$' "$STATUS_FILE")"
    TASK_ID=$(( NUM_TASKS + 1 ))
    echo "$TASK_ID:" "${FUNCNAME[2]}" "$@" -- "$(git_commit)"

    sed -i -E "s/^NUM_TASKS=$NUM_TASKS\$/NUM_TASKS=$TASK_ID/" "$STATUS_FILE"
}


function find_task_record()
{
    # strip the current git commit
    local NEEDLE="${TASK_SIGNATURE% -- *}"
    # search for the identifying part of the task signature
    grep -F "$NEEDLE -- " "$STATUS_FILE"
}


function has_task_been_run()
{
    find_task_record &> /dev/null
}


function git_history_contains()
{
    set +o pipefail
    git rev-list HEAD | grep -qF "$1"
    set -o pipefail
}


function was_task_completed()
{
    local CACHE_HIT="$(find_task_record)"

    if [[ -n "$CACHE_HIT" ]]
    then
        # check if the recorded commit is in the current git history
        git_history_contains "${CACHE_HIT##* -- }"
    else
        return 1
    fi
}


function get_task_opts()
{
    OPTIND=1
    while getopts "clx" OPTION; do
        case "$OPTION" in
            c)
                USE_TASK_CACHE=1
                ;;
            l)
                LONG_TASK=1
                ;;
            x)
                EXEC_COMMAND=1
                ;;
            *)
                echo "$PROG:${FUNCNAME[1]}: unkown option $OPTION" >&2
                ;;
        esac
    done
}


function task()
{
    get_task_opts "$@"
    shift $(( OPTIND - 1 ))
    TASK_SIGNATURE="$(task_signature "$@")"

    TASK_OUT="$(echo -n '-' "$@" '... ')"
    echo -n "$TASK_OUT"

    if [[ -v USE_TASK_CACHE ]] && was_task_completed
    then
        skip_task
        exit 0
    else
        trap fail_task ERR
        trap interrupt_task SIGINT
        trap end_task EXIT
    fi

    [[ "${DEBUG:-}" != succeed-all ]] || exit 0
}

function report_task_status()
{
    [[ ! -v LONG_TASK ]] || echo -n "$TASK_OUT"
    echo -e "$@"
}

function skip_task()
{
    TASK_FAILED=1

    local REASON=""
    (( $# == 0 )) || REASON=" ($(joinstr "$@"))"

    echo -e "\e[33mskipped$REASON\e[0m"
}


function fail_task()
{
    [[ ! -v TASK_FAILED ]] || return

    TASK_FAILED=1
    report_task_status "\e[41mFAILED\e[0m"
    return 1
}


function interrupt_task()
{
    TASK_FAILED=1
    report_task_status "\e[33mcancelled\e[0m"
    exit 42
}


function end_task()
{
    if ! [[ -v TASK_FAILED ]]
    then
        report_task_status "\e[32mdone\e[0m"
        has_task_been_run || echo "$TASK_SIGNATURE" >> "$STATUS_FILE"
    else
        unset TASK_FAILED
    fi
}


function prepare_status_file()
{
    touch "$STATUS_FILE"
    sort -n < "$STATUS_FILE" \
    | sed -E 's/^NUM_TASKS=.*$/NUM_TASKS=0/' \
    > "$STATUS_FILE~"
    mv "$STATUS_FILE~" "$STATUS_FILE"
}


function capture_output()
{
    CAPTURE="$(mktemp .capture.XXXXXX)"

    "$@" > "$CAPTURE" 2>&1 || {
        local status=$?
        cat "$CAPTURE" >&2
    }

    rm -f "$CAPTURE"
    return "${status:-0}"
}


function test_prerequisits()
{
    local MISSING=()
    for CMD in make git dub jq curl snakemake conda-build singularity gio anaconda
    do
        command -v "$CMD" &> /dev/null || MISSING+=("$CMD")
    done

    (( ${#MISSING[*]} == 0 )) \
    || bail_out "missing executable(s): $(joinstr "${MISSING[@]}")"

    if ! docker info &> /dev/null
    then
        echo 'Docker daemon is not running. Starting...'
        sudo systemctl start docker
    fi
}


function assert_clean_git_tree()
(
    task "assert GIT tree is in a clean state"

    git status --porcelain=v1 | awk '
        BEGIN {
            ignore["source/dentist/swinfo.d"] = 1;
            ignore["release.sh"] = 1;
            ignore["conda/conda-build/entrypoint.sh"] = 1;
        }

        {
            idx = substr($0, 1, 1);
            tree = substr($0, 2, 1);
            file = substr($0, 4);
        }

        (!(file in ignore)) { ++mod; }

        END { exit mod; }
    '
)

function git_checkout()
(
    task "checkout $1 branch"

    git checkout -q "$1"
)

function clean_up_after_build()
{
    git checkout -q config-schema.json
}

function assert_dub_build()
(
    task -c 'code compiles with `dub build`'

    dub build --verror
)

function assert_dub_build_testing()
(
    task -c 'code compiles with `dub build --config=testing`'

    dub build --verror --config=testing
)

function assert_unittests_pass()
(
    task -c 'unit tests pass with with `dub test`'

    capture_output dub test --verror
)

function assert_unittests_pass_testing()
(
    task -c 'unit tests pass with with `dub test --config=testing`'

    capture_output dub test --verror --config=testing
)

function assert_pre_push_hook_succeeds()
(
    task -c 'pre-push hook succeeds'

    capture_output ./test-current-commit.sh
)


function build_conda_packages_locally()
(
    task -c "build Conda packages locally"

    capture_output \
        conda-build -c a_ludi -c bioconda \
            --clobber=conda/recipes/dentist-core/local-build.yaml \
            conda/recipes/dentist-core
)

function build_singularity_image_locally()
(
    task -c "build Singularity image locally"

    capture_output sudo singularity build --force example/dist/dentist_local.sif singularity/dentist_local.def
)


function assert_example_works_locally()
(
    task -c "assert example works locally:" "$@"

    capture_output make -o dist/checksum.md5 -C example dentist_version=local "$@"
)

function clean_up_after_local_example()
{
    git checkout -q example/dist/snakemake.yml
}

function get_last_release()
{
    COMMIT="${1:-HEAD}${1:+^}"
    LAST_RELEASE="$(git describe --match='v*' --long "$COMMIT")"
    LAST_RELEASE="${LAST_RELEASE%%-*}"

    if (( $# == 0 ))
    then
        echo -e "- detected last release as \e[32m$LAST_RELEASE\e[0m"
    else
        echo -e "- adjusted last release to \e[32m$LAST_RELEASE\e[0m"
    fi
}

function get_release_version()
{
    if grep -qE '^DENTIST_VERSION=' "$STATUS_FILE"
    then
        eval "$(grep -E '^DENTIST_VERSION=' "$STATUS_FILE")"
        echo -e "- current release is \e[32m$DENTIST_VERSION\e[0m"
    else
        read -rp 'Please enter the next version: ' DENTIST_VERSION
        echo "DENTIST_VERSION=$DENTIST_VERSION" >> "$STATUS_FILE"
    fi
}

function ask_yes_or_no()
{
    while true
    do
        read -rp "$1 (yes/no) "

        if [[ "${REPLY,,}" == yes ]]
        then
            REPLY=0
            break
        elif [[ "${REPLY,,}" == no ]]
        then
            REPLY=1
            break
        else
            echo "Please answer with either 'yes' or 'no'."
        fi
    done
}

function question_task()
{
    get_task_opts "$@"
    shift $(( OPTIND - 1 ))
    TASK_SIGNATURE="$(task_signature "$@")"

    if [[ -v USE_TASK_CACHE ]] && was_task_completed
    then
        SKIP=1
    fi
    [[ "${DEBUG:-}" != succeed-all ]] || SKIP=1

    trap interrupt_task SIGINT
    local ALL_OK=0
    while [[ "$1" != "--" ]]
    do
        if [[ ! -v SKIP ]]
        then
            local Q="$1"

            if [[ -v EXEC_COMMAND ]]
            then
                shift
                eval "$1"
            fi

            ask_yes_or_no "$Q"
            ALL_OK=$(( ALL_OK = REPLY ))
        fi
        shift
    done

    shift
    echo -n '-' "$@" '... '
    [[ "${DEBUG:-}" != succeed-all ]] || unset SKIP

    if [[ -v SKIP ]]
    then
        skip_task
    elif (( ALL_OK == 0 ))
    then
        end_task
    else
        fail_task
    fi

    exit "$ALL_OK"
}

function ask_changelog_updated()
(
    question_task -c \
        "Has the changelog been updated with all changes since $LAST_RELEASE?" \
        -- \
        'changelog has been updated'
)

function insert_commits_since_last_release_into_changelog()
(
    task 'insert commits since last release into changelog'

    cp CHANGELOG.md{,~}
    {
        sed -nE '/^##\s+/ q; p' < CHANGELOG.md~
        echo '> Update changelog appropriately!'
        echo "> The full GIT log since $LAST_RELEASE is displayed in the terminal."
        echo '> Return to the terminal and close the log when done.'
        echo
        echo
        echo "## [Unreleased] - $(date '+%Y-%m-%d')"
        git log --pretty=oneline "$LAST_RELEASE"... | sed -E 's/^[0-9a-f]{40}/-/'
        echo
        echo
        sed -nE '/^##\s+/,$ p' < CHANGELOG.md~
    } > CHANGELOG.md
    rm CHANGELOG.md~
)

function manually_update_changelog()
(
    task 'manually update changelog'

    gui_editor "CHANGELOG.md:$(sed -nE '/^##\s+/ q; p' < CHANGELOG.md | wc -l)" &
    git log "$LAST_RELEASE"...
    git commit -q -m 'Updated CHANGELOG' -e CHANGELOG.md
)

function ask_readme_updated()
(
    question_task -c \
        "Has the README been updated with all changes since $LAST_RELEASE?" \
        "Have all version numbers in README been updated to $DENTIST_VERSION?" \
        "Have all release links in README been updated?" \
        -- \
        'README has been updated'
)

function manually_update_readme()
(
    task 'manually update README'

    gui_editor -w README.md
)



function ask_readme_synced_with_example()
(
    question_task -c \
        "Do the example sections in README.md and example/dist/README.md agree?" \
        -- \
        'example sections in README.md and example/dist/README.md agree'
)

function manually_sync_readme_with_example()
(
    task "manually sync'ed README with example"

    gui_editor -w README.md
)

function git_is_modified()
(
    git status --porcelain=v1 | grep -qF -f <( echo "$*" )
)

function git_commit_if_modified()
{
    OPTS=()
    while [[ "$1" != -- ]]
    do
        OPTS+=( "$1" )
        shift
    done
    shift

    (( $# > 0 )) || bail_out "No files to git_commit_if_modified; did you" \
                             "forget to separate file names using --?"

    if git_is_modified "$@"
    then
        git add "$@"
        git commit "${OPTS[@]}"
    fi
}


function git_files_changed_since_last_release()
{
    local NUM_COMMITS="$(git log --oneline "$LAST_RELEASE.." -- "$@" | wc -l)"
    (( NUM_COMMITS > 0 ))
}


function assert_dentist_container_in_snakefile_is_uptodate()
(
    task -c "updated default dentist_container in Snakefile"

    sed -i -E \
        's#(library://a-ludi/default/dentist:)[0-9]+\.[0-9]+\.[0-9]+#\1'"${DENTIST_VERSION#v}"'#' \
        snakemake/Snakefile
    git_commit_if_modified -q \
            -m 'Updated default dentist_container in Snakefile' \
            -- \
            snakemake/Snakefile
)

function assert_dentist_env_in_snakefile_is_uptodate()
(
    task -c "updated default dentist_env in Snakefile"

    sed -i -E \
        's#(envs/dentist_)v[0-9+](\.yml)#\1'"${DENTIST_VERSION%%.*}"'\2#' \
        snakemake/Snakefile
    git_commit_if_modified -q \
        -m 'Updated default dentist_env in Snakefile' \
        -- \
        snakemake/Snakefile
)

function create_new_dentist_env_definition()
(
    task -c "created new dentist_env definition"

    local LAST_MAJOR="${LAST_RELEASE%%.*}"
    LAST_MAJOR="${LAST_MAJOR#v}"
    local CUR_MAJOR="${DENTIST_VERSION%%.*}"
    CUR_MAJOR="${CUR_MAJOR#v}"
    sed -E \
        "s/(dentist-core==)$LAST_MAJOR(\\.\\*\\.\\*)/\\1$CUR_MAJOR\\2/" \
        < "snakemake/envs/dentist_v$LAST_MAJOR.yml" \
        > "snakemake/envs/dentist_v$CUR_MAJOR.yml"

    git_commit_if_modified -q \
        -m 'Added new DENTIST conda environment definition' \
        -- \
        "snakemake/envs/dentist_v$CUR_MAJOR.yml"
)

function create_new_singularity_image_definition()
(
    task -c "created new singularity image definition"

    local S_LAST="${LAST_RELEASE#v}"
    local S_VER="${DENTIST_VERSION#v}"
    sed -E \
        "s/(DENTIST_VERSION=)${S_LAST//./\\.}/\\1$S_VER/" \
        < "singularity/dentist_$S_LAST.def" \
        > "singularity/dentist_$S_VER.def"

    git_commit_if_modified -q \
        -m 'Added new singularity image definition' \
        -- \
        "singularity/dentist_$S_VER.def"
)

function update_conda_recipes()
(
    task -c "updated conda recipes"

    CONDA_RECIPES=( conda/recipes/dentist{,-core}/meta.yaml )
    for RECIPE in "${CONDA_RECIPES[@]}"
    do
        sed -i -E \
            "s/\"${LAST_RELEASE//./\\.}\"/\"$DENTIST_VERSION\"/" \
            "$RECIPE"
    done

    git_commit_if_modified -q \
        -m 'Updated conda recipes' \
        -- \
        "${CONDA_RECIPES[@]}"
)

function assert_version_bumped()
(
    task -cl 'bumped version'

    ./bump-version.sh
)

function build_conda_release_package_locally()
(
    task -c "build Conda release package locally"

    GIT_URL="../../.." \
    DENTIST_VERSION="$DENTIST_VERSION" \
    capture_output \
        conda-build -c a_ludi -c bioconda \
            conda/recipes/dentist-core
)

function make_singularity_image()
(
    task -c 'made singularity image'

    capture_output make dentist_version="$DENTIST_VERSION" singularity-image
)

function make_dist_tarball()
(
    task -c 'made dist tarball and docker image'

    capture_output make dentist_version="$DENTIST_VERSION" clean-dist
    capture_output make dentist_version="$DENTIST_VERSION" dist
)

function verify_correct_software_version()
(
    task -c 'verified correct software version'

    "dentist.$DENTIST_VERSION.x86_64/bin/dentist" --version |& grep -q "$DENTIST_VERSION"
)

function update_list_of_cli_options()
(
    task -c 'updated list of CLI options'

    capture_output make dentist_version="$DENTIST_VERSION" docs/list-of-commandline-options.md
    git_commit_if_modified -q \
        -m 'Updated list of CLI options' \
        -- \
        docs/list-of-commandline-options.md
)

function git_merge_develop()
(
    task -c 'merged develop into master branch'

    git merge -q --ff-only develop
)

function generate_api_docs()
(
    task -c 'generated api docs'

    capture_output make -B docs.json
)

function go_to_gh_pages_checkout()
{
    (
        task "checking out branch $GH_PAGES_BRANCH in $GH_PAGES_ROOT"

        if [[ ! -d "$GH_PAGES_ROOT" ]]
        then
            mkdir -p "$GH_PAGES_ROOT/.."
            git clone --quiet --local --branch $GH_PAGES_BRANCH "$PWD" "$GH_PAGES_ROOT"
        fi
        cd "$GH_PAGES_ROOT"
        git switch --quiet --no-guess "$GH_PAGES_BRANCH"
    )

    cd "$GH_PAGES_ROOT"
    trap 'p "- left directory $GH_PAGES_ROOT"' EXIT
}

function copy_latest_api_docs_json()
(
    task -c 'copy latest api docs.json'

    cp "$DENTIST/docs.json" "docs.$DENTIST_VERSION.json"
)

function update_current_api_version()
(
    task -c 'update current API version'

    ln -sf "$DENTIST_VERSION" api/current
)

function build_gh_pages()
(
    task -c 'build GitHub pages'

    capture_output make
)

function commit_changes_to_gh_pages()
(
    task -c 'commit changes to GitHub pages'

    git add \
        "docs.$DENTIST_VERSION.json" \
        "api/$DENTIST_VERSION"
    git commit -q -am "Updated site to DENTIST $DENTIST_VERSION"
)

function git_push()
(
    local remote="$1"
    shift

    task -c "pushing" "$@" "to $remote"

    capture_output git push --quiet "$remote" "$@"
)

function get_release_changelog()
{
    sed -nE '
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
    ' CHANGELOG.md
}

function urlencode()
{
    python -c 'from urllib.parse import quote; import sys; print(quote(sys.stdin.read()))'
}

function create_public_release_draft()
(
    echo 'Create a public release draft and remember to include the release tarball.'
    question_task -xc \
        'Have you created the release draft?' \
            'gio open "https://github.com/a-ludi/dentist/releases/new?tag=$DENTIST_VERSION&title=$DENTIST_VERSION&body=$(get_release_changelog | urlencode)"' \
        'Have you included the release tarball?' \
            'true' \
        -- \
        'created public release draft'
)

function update_dub_package()
(
    echo 'Please update the DUB package manually.'
    question_task -xc \
        "Is the DUB package updated to $DENTIST_VERSION?" \
            'gio open https://code.dlang.org/my_packages/dentist' \
        -- \
        'updated DUB package'
)

function build_and_upload_conda_package()
(
    task -c "build and upload conda package $1"

    if git_files_changed_since_last_release "conda/recipes/$1"
    then
        cd conda
        capture_output ./conda-build.sh --upload "recipes/$1"
    else
        skip_task 'no changes since last release'
    fi
)


function build_singularity_image()
(
    task -c "build Singularity image"

    capture_output sudo singularity build "dentist_${DENTIST_VERSION#v}.sif" "singularity/dentist_${DENTIST_VERSION#v}.def"
    capture_output sudo chown "$USER:" "dentist_${DENTIST_VERSION#v}.sif"
)


function sign_singularity_image()
(
    task -cl "sign Singularity image"

    singularity sign "dentist_${DENTIST_VERSION#v}.sif"
)


function upload_singularity_image_to_singularity_hub()
(
    task -c "upload Singularity image to Singularity hub"

    capture_output singularity push "dentist_${DENTIST_VERSION#v}.sif" "library://a-ludi/default/dentist:${DENTIST_VERSION#v}"
)


function upload_singularity_image_to_github()
(
    question_task -xc \
        "Have you added the Singularity image ($PWD/dentist_${DENTIST_VERSION#v}.sif) to the release draft?" \
            'gio open "https://github.com/a-ludi/dentist/releases?q=$DENTIST_VERSION"' \
        -- \
        "upload Singularity image"
)


function prepare_example_dist()
(
    task -c "prepare example distribution"

    capture_output make -C example dentist_version="$DENTIST_VERSION" dist
    git commit -q -m "Updated DENTIST example to $DENTIST_VERSION" example
)


function validate_example_dist()
(
    task -c "validate example distribution"

    function _validate()
    {
        cd /tmp
        tar -xzf "$DENTIST/example/dentist-example.tar.gz"
        cd dentist-example/
        PATH="$PWD/bin:$PATH" snakemake --configfile=snakemake.yml --cores=all
        md5sum -c checksum.md5
        [[ -f "dentist_${DENTIST_VERSION#v}.sif" ]] || echo "error: SIF file missing"
    }

    capture_output _validate
)


function add_example_dist_to_release_draft()
(
    echo "Add the example dist tarball ($DENTIST/example/dentist-example.tar.gz) to the release draft."
    question_task -xc \
        'Have you included the example tarball?' \
            'gio open "https://github.com/a-ludi/dentist/releases?q=$DENTIST_VERSION"' \
        -- \
        'add example dist tarball to release draft'
)


function verify_release()
(
    echo 'Check the published content.'

    question_task -xc \
        'Is the released example working?' \
            'gio open https://github.com/a-ludi/dentist#example' \
        'Is the GitHub page correctly updated?' \
            'gio open https://a-ludi.github.io/dentist' \
        'Is the DUB package correctly updated?' \
            'gio open https://code.dlang.org/packages/dentist' \
        'Is the Singularity image correctly updated?' \
            'gio open https://cloud.sylabs.io/library/a-ludi/default/dentist' \
        'Is the Conda package correctly updated?' \
            'gio open https://anaconda.org/a_ludi/dentist-core' \
        -- \
        'verified release'
)


function publish_release()
(
    echo 'Publish the release draft on GitHub.'
    question_task -xc \
        'Is the release published?' \
            'gio open "https://github.com/a-ludi/dentist/releases?q=$DENTIST_VERSION"' \
        -- \
        'published release'
)


function main()
{
    prepare_status_file

    p "Welcome to the DENTIST release script!"
    p "======================================"
    test_prerequisits


    part "Quality Control"

    assert_clean_git_tree
    git_checkout develop
    assert_dub_build
    clean_up_after_build
    assert_dub_build_testing
    clean_up_after_build
    assert_unittests_pass
    assert_unittests_pass_testing
    assert_pre_push_hook_succeeds
    build_conda_packages_locally
    build_singularity_image_locally
    assert_example_works_locally prepare-dist recreate-conda-env test-conda && clean_up_after_local_example
    assert_example_works_locally prepare-dist test-binaries && clean_up_after_local_example
    assert_example_works_locally prepare-dist test-singularity && clean_up_after_local_example


    part "Prepare DENTIST release"

    get_last_release
    get_release_version

    [[ "$LAST_RELEASE" != "$DENTIST_VERSION" ]] || get_last_release "$DENTIST_VERSION"

    if ! ask_changelog_updated
    then
        insert_commits_since_last_release_into_changelog
        manually_update_changelog
    fi
    ask_readme_updated || manually_update_readme
    ask_readme_synced_with_example || manually_sync_readme_with_example

    create_new_dentist_env_definition
    assert_dentist_container_in_snakefile_is_uptodate
    assert_dentist_env_in_snakefile_is_uptodate
    create_new_singularity_image_definition
    update_conda_recipes
    assert_pre_push_hook_succeeds
    assert_version_bumped
    build_conda_release_package_locally
    make_dist_tarball
    verify_correct_software_version
    update_list_of_cli_options
    # make_singularity_image
    assert_pre_push_hook_succeeds
    git_checkout master
    git_merge_develop


    part "Update GitHub Pages"

    generate_api_docs

    (
        go_to_gh_pages_checkout
        copy_latest_api_docs_json
        update_current_api_version
        build_gh_pages
        commit_changes_to_gh_pages
    )


    part "Publish release"

    echo -e '\e[41mAttention:\e[0m The release will be published now.'
    echo 'Make sure that all of the above steps succeed before publishing.'
    echo
    ask_yes_or_no 'Ready to publish?'

    git_push "$PUBLIC_GIT" develop master "$DENTIST_VERSION"
    create_public_release_draft
    update_dub_package
    (
        go_to_gh_pages_checkout
        git_push "$PUBLIC_GIT" gh-pages
    )

    build_and_upload_conda_package libunwind-static
    build_and_upload_conda_package daccord
    build_and_upload_conda_package daligner
    build_and_upload_conda_package damapper
    build_and_upload_conda_package damasker
    build_and_upload_conda_package dascrubber
    build_and_upload_conda_package dazz_db
    build_and_upload_conda_package dentist-core
    # DEPREACTED build_and_upload_conda_package dentist
    build_singularity_image
    sign_singularity_image
    upload_singularity_image_to_singularity_hub
    upload_singularity_image_to_github

    prepare_example_dist
    validate_example_dist
    add_example_dist_to_release_draft
    publish_release
    verify_release

    echo
    echo -e "\e[1mRelease \e[1;32m$DENTIST_VERSION\e[0m\e[1m successfully published!\e[0m"
}


main "$@"

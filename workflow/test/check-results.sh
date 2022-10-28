#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'
PROG="$0"

# Add the included binaries in PATH so we can use them later
PATH="$PWD/bin:$PATH"

# This file will include details of the checks
LOG="check-results.log"

# These files must exist after a successful run of the workflow
RESULT_FILES=(
    gap-closed.agp
    gap-closed.closed-gaps.bed
    gap-closed.fasta
    workdir/gap-closed-preliminary.1.reads.las
    workdir/gap-closed-preliminary.agp
    workdir/.gap-closed-preliminary.bps
    workdir/.gap-closed-preliminary.closed-gaps.anno
    workdir/gap-closed-preliminary.closed-gaps.bed
    workdir/.gap-closed-preliminary.closed-gaps.data
    workdir/gap-closed-preliminary.dam
    workdir/.gap-closed-preliminary.dentist-self.anno
    workdir/.gap-closed-preliminary.dentist-self.data
    workdir/.gap-closed-preliminary.dentist-weak-coverage-1B.anno
    workdir/.gap-closed-preliminary.dentist-weak-coverage-1B.data
    workdir/.gap-closed-preliminary.dentist-weak-coverage.anno
    workdir/.gap-closed-preliminary.dentist-weak-coverage.data
    workdir/.gap-closed-preliminary.dust.anno
    workdir/.gap-closed-preliminary.dust.data
    workdir/gap-closed-preliminary.fasta
    workdir/gap-closed-preliminary.gap-closed-preliminary.las
    workdir/.gap-closed-preliminary.hdr
    workdir/.gap-closed-preliminary.idx
    workdir/gap-closed-preliminary.reads.10.las
    workdir/gap-closed-preliminary.reads.11.las
    workdir/gap-closed-preliminary.reads.12.las
    workdir/gap-closed-preliminary.reads.1.las
    workdir/gap-closed-preliminary.reads.2.las
    workdir/gap-closed-preliminary.reads.3.las
    workdir/gap-closed-preliminary.reads.4.las
    workdir/gap-closed-preliminary.reads.5.las
    workdir/gap-closed-preliminary.reads.6.las
    workdir/gap-closed-preliminary.reads.7.las
    workdir/gap-closed-preliminary.reads.8.las
    workdir/gap-closed-preliminary.reads.9.las
    workdir/gap-closed-preliminary.reads.las
    workdir/.gap-closed-preliminary.tan.anno
    workdir/.gap-closed-preliminary.tan.data
    workdir/insertions/batch.0-50.db
    workdir/insertions/batch.100-146.db
    workdir/insertions/batch.50-100.db
    workdir/insertions.db
    workdir/pile-ups.db
    workdir/reads.10.gap-closed-preliminary.las
    workdir/reads.10.reference.las
    workdir/reads.11.gap-closed-preliminary.las
    workdir/reads.11.reference.las
    workdir/reads.12.gap-closed-preliminary.las
    workdir/reads.12.reference.las
    workdir/reads.1.gap-closed-preliminary.las
    workdir/reads.1.reference.las
    workdir/reads.2.gap-closed-preliminary.las
    workdir/reads.2.reference.las
    workdir/reads.3.gap-closed-preliminary.las
    workdir/reads.3.reference.las
    workdir/reads.4.gap-closed-preliminary.las
    workdir/reads.4.reference.las
    workdir/reads.5.gap-closed-preliminary.las
    workdir/reads.5.reference.las
    workdir/reads.6.gap-closed-preliminary.las
    workdir/reads.6.reference.las
    workdir/reads.7.gap-closed-preliminary.las
    workdir/reads.7.reference.las
    workdir/reads.8.gap-closed-preliminary.las
    workdir/reads.8.reference.las
    workdir/reads.9.gap-closed-preliminary.las
    workdir/reads.9.reference.las
    workdir/.reads.bps
    workdir/reads.db
    workdir/.reads.dentist-reads-1-8B.anno
    workdir/.reads.dentist-reads-1-8B.data
    workdir/.reads.dentist-reads-9-12B.anno
    workdir/.reads.dentist-reads-9-12B.data
    workdir/.reads.dentist-self-1-8B.anno
    workdir/.reads.dentist-self-1-8B.data
    workdir/.reads.dentist-self-9-12B.anno
    workdir/.reads.dentist-self-9-12B.data
    workdir/.reads.idx
    workdir/.reads.tan-1-8B.anno
    workdir/.reads.tan-1-8B.data
    workdir/.reads.tan-9-12B.anno
    workdir/.reads.tan-9-12B.data
    workdir/.reference.1.tan.anno
    workdir/.reference.1.tan.data
    workdir/.reference.bps
    workdir/reference.dam
    workdir/.reference.dentist-reads.anno
    workdir/.reference.dentist-reads.data
    workdir/.reference.dentist-reads-H-1-8B.anno
    workdir/.reference.dentist-reads-H-1-8B.data
    workdir/.reference.dentist-reads-H-9-12B.anno
    workdir/.reference.dentist-reads-H-9-12B.data
    workdir/.reference.dentist-reads-H.anno
    workdir/.reference.dentist-reads-H.data
    workdir/.reference.dentist-self.anno
    workdir/.reference.dentist-self.data
    workdir/.reference.dentist-self-H-1-8B.anno
    workdir/.reference.dentist-self-H-1-8B.data
    workdir/.reference.dentist-self-H-9-12B.anno
    workdir/.reference.dentist-self-H-9-12B.data
    workdir/.reference.dentist-self-H.anno
    workdir/.reference.dentist-self-H.data
    workdir/.reference.dust.anno
    workdir/.reference.dust.data
    workdir/.reference.hdr
    workdir/.reference.idx
    workdir/reference.reads.10.las
    workdir/reference.reads.11.las
    workdir/reference.reads.12.las
    workdir/reference.reads.1.las
    workdir/reference.reads.2.las
    workdir/reference.reads.3.las
    workdir/reference.reads.4.las
    workdir/reference.reads.5.las
    workdir/reference.reads.6.las
    workdir/reference.reads.7.las
    workdir/reference.reads.8.las
    workdir/reference.reads.9.las
    workdir/reference.reads.las
    workdir/reference.reference.las
    workdir/.reference.tan.anno
    workdir/.reference.tan.data
    workdir/.reference.tan-H-1-8B.anno
    workdir/.reference.tan-H-1-8B.data
    workdir/.reference.tan-H-9-12B.anno
    workdir/.reference.tan-H-9-12B.data
    workdir/.reference.tan-H.anno
    workdir/.reference.tan-H.data
    workdir/skip-gaps.txt
    workdir/TAN.gap-closed-preliminary.1.las
    workdir/TAN.gap-closed-preliminary.las
    workdir/TAN.reference.1.las
    workdir/TAN.reference.las
    workdir/validation-report.1.json
    workdir/validation-report.json

)

# Checksums for files that cannot be compare directly but need some form of
# pre-processing
declare -A CHECKSUMS
# [[ -f '.checksumsrc' ]] && source '.checksumsrc'
CHECKSUMS=(
    ['workdir/reference.reference.las']=3852a9d97f4ab3adf597a575d4e500da
    ['workdir/reference.reads.las']=a34110d345aa8da6fda869869dde6a8e
    ['workdir/pile-ups.db']=7ca62b79e38399da42a0f8ae7e3dffc6
    ['workdir/insertions.db']=d36c351322d22a7163e8e3b5c65ee8d9
    ['workdir/gap-closed-preliminary.gap-closed-preliminary.las']=3852a9d97f4ab3adf597a575d4e500da
    ['workdir/gap-closed-preliminary.reads.las']=a34110d345aa8da6fda869869dde6a8e
    ['gap-closed.agp']=c5b1d6df1291c83208fd7892fafc692f
)


# Check failures will be stored here
export ERROR_FILE

# Report the test as skipped; do not record an error; continue execution
ERROR_SKIP=1
# Report given error as failure; continue execution
ERROR_LOG=2
# Report given error as failure; stop execution
ERROR_FAIL=3

# Private: call upon assertions failure; captures the name of the assertion,
# records an error message and return with the given error level (see above).
function __fail()
{
    local ERROR_LEVEL="$1"
    shift
    local ASSERTION="${ASSERTION:-"${FUNCNAME[2]}"}"
    local IFS=' '

    echo "[FAILED] $ASSERTION: $*" >> "$ERROR_FILE"

    return "$ERROR_LEVEL"
}


# Record the error alongside the assertion name and continue the tests.
function skip_test()
{
    return "$ERROR_SKIP"
}


# Record the error alongside the assertion name and continue the tests.
function log_error()
{
    __fail "$ERROR_LOG" "$@"
}


# Record the error alongside the assertion name and stop execution.
function fail_immediately()
{
    __fail "$ERROR_FAIL" "$@"
}


function prepare_output_files()
{
    # [optional] instead of checking md5sums compute them and print a bash
    #   array that can be used by this script
    if [[ -v WRITE_CHECKSUMS ]]
    then
        echo "CHECKSUMS=(" > "$WRITE_CHECKSUMS"
    fi

    # clear the log file
    cat /dev/null > "$LOG"
    # create temp file for errors
    ERROR_FILE="$(mktemp --tmpdir check-results.XXXXXX.err)"

    trap finish_output_files exit
}


function finish_output_files()
{
    # see prepare_output_files
    if [[ -v WRITE_CHECKSUMS ]]
    then
        echo ")" >> "$WRITE_CHECKSUMS"
    fi

    rm -f "$ERROR_FILE"
}


function assert_results_present()
{
    for RESULT_FILE in "${RESULT_FILES[@]}"
    do
        [[ -f "$RESULT_FILE" ]] \
        || fail_immediately "missing result file $RESULT_FILE;" \
            "complete the workflow!"
    done
}


# Check only the given files from checksum.md5
function check_md5sums_partial()
{
    local ON_ERROR="${ON_ERROR:-log_error}"
    local IFS=$' '
    local MESSAGE="${MESSAGE:-"corrupted files among: $*"}"

    local IFS=$'\n'
    awk '
        (FILENAME == "-") { selection[$1] = 1 }
        (FILENAME != "-" && $2 in selection) { print }
    ' - <<<"$*" checksum.md5 \
    | md5sum -c - &>> "$LOG" \
    || ASSERTION="${ASSERTION:-"${FUNCNAME[1]}"}" "$ON_ERROR" "$MESSAGE"
}


# Check the md5sum of stdin against the hash stored in CHECKSUMS[$1]
function check_md5sum_stdin()
{
    local FILE="$1"
    local ERROR="${2:-"corrupted file: $FILE"}"

    if [[ ! -v WRITE_CHECKSUMS ]]
    then
        md5sum -c <(echo "${CHECKSUMS["$FILE"]}  -") &>> "$LOG" \
        || ASSERTION="${ASSERTION:-"${FUNCNAME[1]}"}" log_error "$ERROR"
    else
        {
            echo -n "    [${FILE@Q}]="
            md5sum - | tr -d ' -'
        } >> "$WRITE_CHECKSUMS"
        skip_test
    fi
}


function check_inputs()
{
    ON_ERROR=fail_immediately \
    MESSAGE='corrupted input DBs; try removing the entire directory and extracting dentist-example.tar.gz again' \
    ASSERTION="${FUNCNAME[0]}" \
    check_md5sums_partial \
        workdir/.reference.bps \
        workdir/.reference.hdr \
        workdir/.reference.idx \
        workdir/.reads.bps \
        workdir/.reads.idx \
        workdir/reference.dam \
        workdir/reads.db
}


function check_mask()
{
    local DB="$1"
    shift

    for MASK in "$@"
    do
        ASSERTION="${FUNCNAME[0]} $DB $MASK" \
        check_md5sums_partial "workdir/.$DB.$MASK."{anno,data}
    done
}


function check_alignment()
{
    local IFS=' '
    local DB_A="$1"
    local DB_B="${2:-"$DB_A"}"
    local BLOCK_A="${3:+.}${3:-}"
    local BLOCK_B="${4:+.}${4:-}"
    local LAS="workdir/$DB_A$BLOCK_A.$DB_B$BLOCK_B.las"

    LAdump -cdtl "workdir/$DB_A" "workdir/$DB_B" "$LAS" \
    | ASSERTION="${FUNCNAME[0]} $*" check_md5sum_stdin "$LAS"
}


function check_pile_ups()
{
    dentist show-pile-ups --json -v -v -v workdir/pile-ups.db \
    |& jq -cS 'select(has("readAlignments") or has("numPileUps"))' \
    | check_md5sum_stdin workdir/pile-ups.db
}


function check_insertions()
{
    dentist show-insertions --json -v -v -v workdir/insertions.db \
    |& jq -cS 'select(has("payload") or has("numInsertions"))' \
    | check_md5sum_stdin workdir/insertions.db
}


function check_preliminary_gap_closed()
{
    ASSERTION="${FUNCNAME[0]}" \
    check_md5sums_partial \
        workdir/.gap-closed-preliminary.bps \
        workdir/.gap-closed-preliminary.hdr \
        workdir/.gap-closed-preliminary.idx \
        workdir/gap-closed-preliminary.dam
}


function check_validation_result()
{
    ASSERTION="${FUNCNAME[0]}" \
    check_md5sums_partial \
        workdir/validation-report.json \
        workdir/skip-gaps.txt
    check_mask gap-closed-preliminary dentist-weak-coverage
}


function check_gap_closed_assembly()
{
    grep -vE '^#\s*TOOL:' gap-closed.agp \
    | check_md5sum_stdin gap-closed.agp
    ASSERTION="${FUNCNAME[0]}" \
    check_md5sums_partial \
        gap-closed.closed-gaps.bed \
        gap-closed.fasta
}


function all_checks()
{
    CHECKS=(
        assert_results_present
        check_inputs
        'check_mask reference dust'
        'check_mask reference tan'
        'check_alignment reference'
        'check_mask reference dentist-self'
        'check_alignment reference reads'
        'check_mask reference dentist-reads'
        'check_mask reads tan-1-8B'
        'check_mask reads tan-9-12B'
        'check_mask reads dentist-self-1-8B'
        'check_mask reads dentist-self-9-12B'
        'check_mask reads dentist-reads-1-8B'
        'check_mask reads dentist-reads-9-12B'
        'check_mask reference tan-H'
        'check_mask reference dentist-self-H'
        'check_mask reference dentist-reads-H'
        check_pile_ups
        check_insertions
        check_preliminary_gap_closed
        'check_mask gap-closed-preliminary dust'
        'check_mask gap-closed-preliminary tan'
        'check_alignment gap-closed-preliminary'
        'check_mask gap-closed-preliminary dentist-self'
        'check_alignment gap-closed-preliminary reads'
        'check_mask gap-closed-preliminary closed-gaps'
        check_validation_result
        check_gap_closed_assembly
    )
    NUM_PASSED=0
    NUM_FAILED=0
    NUM_SKIPPED=0

    prepare_output_files
    echo -n "Running ${#CHECKS[*]} checks: "
    for CHECK in "${CHECKS[@]}"
    do
        if eval "$CHECK"
        then
            (( ++NUM_PASSED ))
            echo -n '.'
        else
            local ERROR_LEVEL="$?"

            if (( ERROR_LEVEL == ERROR_SKIP ))
            then
                (( ++NUM_SKIPPED ))
                echo -n 's'
            else
                (( ++NUM_FAILED ))
                echo -n 'f'

                (( ERROR_LEVEL < ERROR_FAIL )) || break
            fi
        fi
    done
    echo

    if (( NUM_FAILED > 0 ))
    then
        cat "$ERROR_FILE"
        echo
        echo "Details can be found in $LOG"
    fi

    echo "Check results: $NUM_PASSED passed, $NUM_SKIPPED skipped, $NUM_FAILED failed"

    (( NUM_FAILED == 0 ))
}

all_checks "$@"

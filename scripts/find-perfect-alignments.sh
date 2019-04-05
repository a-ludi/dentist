#!/bin/bash

function main()
{
    NUM_THREADS=1

    while getopts "kt:T:" OPTION; do
        case "$OPTION" in
            k)
                KEEP_TEMP=1
                ;;
            t)
                TEMP_DIR="$OPTARG"
                ;;
            T)
                NUM_THREADS="$OPTARG"
                ;;
            *)
                echo "$(basename "$0"): unkown option $OPTION" >&2
                echo >&2
                usage
                ;;
        esac
    done
    shift $(($OPTIND - 1))

    REFERENCE="$(realpath "$1")"
    shift
    TARGETS=($*)
    for (( i = 0; i < ${#TARGETS[*]}; i++ )); do
        TARGETS[i]="$(realpath "${TARGETS[i]}")"
    done

    WORKDIR="$(mktemp -d find-perfect-alignments.XXXXXX)"
    cd "$WORKDIR"

    REFERENCE_CONTIGS="$(assert_contigs_fasta "$REFERENCE" 'ref')"
    for TARGET in ${TARGETS[*]}; do
        TARGET_CONTIGS="$(assert_contigs_fasta "$TARGET" 'result')"
        PREFIX="$(db_name "$TARGET")"

        nucmer \
            --threads="$NUM_THREADS" \
            --forward \
            --prefix=$PREFIX \
            "$TARGET_CONTIGS" \
            "$REFERENCE_CONTIGS"
        delta-filter -i 100 -q "$PREFIX.delta" > "$PREFIX.F.delta"

        if [[ $REFERENCE == $TARGET ]]; then
            FILTER='{ if ($3 < $4 && ($5 == $9 || $6 == $9) && $10 != $11) print $0 }'
        else
            FILTER='{ if ($3 < $4 && $5 == $6 && $6 == $9) print $0 }'
        fi

        show-coords -l -T -q "$PREFIX.F.delta" | \
            tail -n+5 | \
            awk -F'\t' "$FILTER"
    done
}

function assert_contigs_fasta()
{
    DB="$1"
    CONTIGS_FASTA="$(db_name "$1").contigs.fasta"

    contigs_fasta "$DB" "$2" > "$CONTIGS_FASTA"

    echo "$CONTIGS_FASTA"
}

function contigs_fasta()
{
    contig_dump "$1" | awk -F'\t' '{ printf ">'"$2"'-%09d %s\n%s\n", $1, substr($2, 2), $3 }'
}

function contig_dump()
{
    DBdump -rhs "$1" | sed -nE '/^R/ {s/^R\s+//;h;n; s/^H[^>]+//;H;n;n; s/^S[^actg]+//;H;x;s/\n/\t/g; p}'
}

function db_name()
{
    basename "$1" ".dam"
}

function mktemp()
{
    if ! [[ -v __TEMP_FILES ]]; then
        declare -a __TEMP_FILES
    fi

    declare -a OPTS=()
    if [[ -v TEMP_DIR ]]; then
        OPTS+=("--tmpdir=$TEMP_DIR")
    else
        OPTS+=("--tmpdir")
    fi

    local TMPFILE="$(/usr/bin/env mktemp "${OPTS[@]}" "$@")"

    local IFS=$'\n'
    __TEMP_FILES+=($TMPFILE)

    echo -n "$TMPFILE"
}

function __cleanup()
{
    if [[ ! -v KEEP_TEMP && -v __TEMP_FILES ]]; then
        local IFS=$'\n'

        for TMPFILE in ${__TEMP_FILES[*]}; do
            rm -rf "$TMPFILE"
        done

        unset __TEMP_FILES
    fi
}

trap __cleanup exit err

function alternate_ifs()
{
    OLD_IFS="$IFS"
    IFS="$1"
}

function reset_ifs()
{
    IFS="$OLD_IFS"
    unset OLD_IFS
}

if ! [[ $- =~ i ]]; then
    # Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
    set -euo pipefail
    IFS=$'\n'

    main "$@"
fi

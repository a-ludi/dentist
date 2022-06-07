#!/bin/bash

SCHEMA_FILE=config-schema.json


#!/bin/bash


function main()
{
    echo -n 'Updating `'"$SCHEMA_FILE"'` ... ' >&2

    ./dentist --config-schema > "$SCHEMA_FILE~"

    if ! cmp -s "$SCHEMA_FILE" "$SCHEMA_FILE~";
    then
        mv "$SCHEMA_FILE~" "$SCHEMA_FILE"

        echo 'done' >&2
    else
        echo 'skipped' >&2
    fi
}

function clean_up()
{
    rm -f "$SCHEMA_FILE~"
}

function on_error()
{
    echo failed >&2
}

trap clean_up EXIT
trap on_error ERR
set -e  # exit when any command fails

main

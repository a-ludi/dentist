#!/bin/bash

SCHEMA_FILE=config-schema.json


#!/bin/bash


function main()
{
    echo -n 'Updating `'"$SCHEMA_FILE"'` ... '

    if [[ "$DUB_BUILD_TYPE" == unittest || "$DUB_BUILD_TYPE" == docs-json ]]
    then
        echo 'skipped'
        return
    fi

    "$DUB_BUILD_PATH/dentist" --config-schema > "$SCHEMA_FILE~"

    if ! cmp -s "$SCHEMA_FILE" "$SCHEMA_FILE~";
    then
        mv "$SCHEMA_FILE~" "$SCHEMA_FILE"

        echo 'done'
    else
        echo 'skipped'
    fi
}

function clean_up()
{
    rm -f "$SCHEMA_FILE~"
}

function on_error()
{
    echo failed
}

trap clean_up EXIT
trap on_error ERR
set -e  # exit when any command fails

main

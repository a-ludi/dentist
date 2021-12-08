#!/bin/bash

# shellcheck source=include.rc, disable=SC1091
source "$(dirname "$0")/include.rc"

cat "$@" \
| remove_logo \
| remove_short_description \
| remove_table_of_contents \
| adjust_relative_paths

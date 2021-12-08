#!/bin/bash

# shellcheck source=include.rc, disable=SC1091
source "$(dirname "$0")/include.rc"

cat "$@" \
| add_back_link

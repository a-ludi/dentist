#!/bin/bash

# shellcheck source=include.rc, disable=SC1091
source "$(dirname "$0")/include.rc"

cat "$@" \
| add_md_header "License" \
| add_back_link \
| adjust_email_links

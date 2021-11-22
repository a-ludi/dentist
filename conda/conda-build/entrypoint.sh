#!/usr/bin/env dash

set -eu

# shellcheck disable=SC1091
. /app/activate-conda.sh

export CALLER_PID=$$
trap 'exit 1' USR1
check_pipe()
{
    eval "$@" || kill -SIGUSR1 $CALLER_PID
}

set -x

check_pipe conda-build "$@" /recipe 2>&1 \
| tee /dev/stderr \
| sed -nE '/^anaconda\s+upload\s+\\/ { n; s/^\s+/TARBALL=/p }'
